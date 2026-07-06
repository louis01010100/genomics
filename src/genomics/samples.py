import logging
import shutil
from pathlib import Path

from pathos.multiprocessing import ProcessPool

from .vcf import Vcf, concat


def export_samples(vcf_files, samples_file, output_dir, n_threads=1):
    """Produce one GT-only, single-sample, sorted, bgzipped, indexed VCF per
    target sample under `output_dir`. Each per-sample VCF concatenates that
    sample's records across every source VCF that contains it. See the
    genomics `samples` spec."""
    output_dir = Path(output_dir)
    samples = _load_samples(Path(samples_file))
    vcf_files = [Path(vcf_file) for vcf_file in vcf_files]

    tmp_dir = output_dir / 'tmp'
    staging = tmp_dir / 'staging'
    tmp_dir.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True, exist_ok=True)

    # Stage 1: up-front presence check -- union of the sources' header samples.
    # Absent-from-all samples warn and are skipped; if none are present, error.
    source_samples = {}
    union = set()
    for vcf_file in vcf_files:
        source_samples[vcf_file] = Vcf(vcf_file, tmp_dir, new_tmp=False).samples
        union |= source_samples[vcf_file]

    present_samples = [sample for sample in samples if sample in union]
    absent = [sample for sample in samples if sample not in union]

    if not present_samples:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise ValueError(
            'all requested samples are absent from every source: '
            f'{sorted(absent)}')
    if absent:
        logging.warning(
            '%d sample(s) absent from all sources, skipped: %s',
            len(absent), ', '.join(absent))

    # Stage 2: per source, strip -> select -> split (one piped pass), run
    # concurrently across sources -- each source is independent. Restricted to
    # the present samples that source actually contains.
    tasks = [(vcf_file, set(present_samples) & source_samples[vcf_file])
             for vcf_file in vcf_files]

    def _process_source(task):
        vcf_file, present_here = task
        if not present_here:
            return {}
        # new_tmp=False keeps the per-source working files under a readable
        # tmp/<source>-samples/ dir instead of an opaque tmp/<md5>/ dir.
        return Vcf(vcf_file, tmp_dir, n_threads,
                   new_tmp=False).strip_select_split(present_here)

    if n_threads > 1 and len(tasks) > 1:
        with ProcessPool(n_threads) as pool:
            per_source = pool.map(_process_source, tasks)
    else:
        per_source = [_process_source(task) for task in tasks]

    # Collect pieces per sample in source order -- deterministic and independent
    # of source completion order.
    sample2pieces = {sample: [] for sample in present_samples}
    for pieces in per_source:
        for sample, piece in pieces.items():
            sample2pieces[sample].append(Path(piece))

    # Stage 3: per sample, concat (+ sort + bgzip + index) into staging.
    if n_threads > 1 and len(present_samples) > 1:
        with ProcessPool(n_threads) as pool:
            pool.map(
                lambda sample: _build_sample(sample, sample2pieces[sample],
                                             staging, tmp_dir),
                present_samples,
            )
    else:
        for sample in present_samples:
            _build_sample(sample, sample2pieces[sample], staging, tmp_dir)

    # Stage 4: publish atomically into a flat samples/ subfolder, then clean up.
    samples_dir = output_dir / 'samples'
    samples_dir.mkdir(parents=True, exist_ok=True)
    for sample in present_samples:
        for suffix in ('.vcf.bgz', '.vcf.bgz.csi', '.vcf.bgz.tbi'):
            piece = staging / f'{sample}{suffix}'
            if piece.exists():
                shutil.move(str(piece), str(samples_dir / piece.name))

    shutil.rmtree(tmp_dir, ignore_errors=True)


def _build_sample(sample, pieces, staging, tmp_dir):
    concat(
        vcf_files=[str(piece) for piece in pieces],
        output_file=staging / f'{sample}.vcf.bgz',
        tmp_dir=tmp_dir / f'concat-{sample}',
        n_threads=1,
    )
    return sample


def _load_samples(samples_file: Path) -> list:
    """De-duplicated target sample names from a one-column file whose first line
    is the header `sample_name`, preserving first-seen order."""
    seen = set()
    ordered = []
    with open(samples_file, 'rt') as fh:
        header = fh.readline().strip()
        if header != 'sample_name':
            raise ValueError(
                "samples-file must be a one-column file with a 'sample_name' "
                f"header (got {header!r})")
        for line in fh:
            name = line.strip()
            if not name or name in seen:
                continue
            seen.add(name)
            ordered.append(name)
    return ordered
