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
    sample_set = set(samples)

    tmp_dir = output_dir / 'tmp'
    staging = tmp_dir / 'staging'
    tmp_dir.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True, exist_ok=True)

    # Stage 2: per source, strip -> select -> split as a single piped pass,
    # restricted to the target samples actually present in that source.
    sample2pieces = {sample: [] for sample in samples}
    for vcf_file in vcf_files:
        source = Vcf(Path(vcf_file), tmp_dir, n_threads)
        present = sample_set & source.samples
        if not present:
            continue
        for sample, piece in source.strip_select_split(present).items():
            sample2pieces[sample].append(Path(piece))

    # Stage 3: fail fast -- a requested sample absent from all sources.
    missing = [sample for sample in samples if not sample2pieces[sample]]
    if missing:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise ValueError(f'samples absent from all sources: {sorted(missing)}')

    # Stage 4-5: per sample, concat (+ sort + bgzip + index) into staging.
    if n_threads > 1 and len(samples) > 1:
        with ProcessPool(n_threads) as pool:
            pool.map(
                lambda sample: _build_sample(sample, sample2pieces[sample],
                                             staging, tmp_dir),
                samples,
            )
    else:
        for sample in samples:
            _build_sample(sample, sample2pieces[sample], staging, tmp_dir)

    # Stage 6: publish atomically (nothing reaches output_dir until every
    # sample built), then remove temporaries.
    for sample in samples:
        for suffix in ('.vcf.bgz', '.vcf.bgz.csi', '.vcf.bgz.tbi'):
            piece = staging / f'{sample}{suffix}'
            if piece.exists():
                shutil.move(str(piece), str(output_dir / piece.name))

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
    """De-duplicated target sample names, preserving first-seen order."""
    seen = set()
    ordered = []
    with open(samples_file, 'rt') as fh:
        for line in fh:
            name = line.strip()
            if not name or name in seen:
                continue
            seen.add(name)
            ordered.append(name)
    return ordered
