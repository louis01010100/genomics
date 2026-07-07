"""Unit tests for the reworked `snv-truth` pure helpers in genomics.truth.

These target the backbone-driven exact-match rework: reading the backbone VCF as
coordinates, grouping into SNP families by ID, loading the two `cram-depth` tables,
the contig-class/gender-aware depth selection, the exact-match key, and the fill
homref/nocall decision.
"""
import gzip
import math
import shutil
import subprocess
from pathlib import Path

import pytest

requires_tabix = pytest.mark.skipif(
    shutil.which('bgzip') is None or shutil.which('tabix') is None,
    reason='bgzip/tabix not on PATH',
)


def _index_depths(tsv_file: Path, text: str) -> Path:
    """Write a depth TSV then bgzip + tabix-index it exactly like the producer,
    returning the .bgz path the loaders random-search."""
    tsv_file.write_text(text)
    bgz_file = Path(f'{tsv_file}.bgz')
    subprocess.run(f'bgzip -f -c {tsv_file} > {bgz_file}', shell=True, check=True)
    subprocess.run(
        ['tabix', '-f', '-s', '1', '-b', '2', '-e', '2', '-S', '1', str(bgz_file)],
        check=True,
    )
    return bgz_file

from genomics.truth import (
    load_coordinates,
    load_snv_family,
    group_families,
    DepthProvider,
    depth_for,
    call_families,
    match_key,
    is_homref,
    fixploidy_lines,
    load_ploidy,
    _require_single_sample,
    prepare_sample,
)

BACKBONE = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr1>\n'
    '##contig=<ID=chrX>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    # a 3-ALT family: full multiallelic + its splits, same ID string repeated
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tC,G,T\t.\t.\t.\n'
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tc\t.\t.\t.\n'
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tG\t.\t.\t.\n'
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tT\t.\t.\t.\n'
    # a lone single-ALT family
    'chr1\t11080559\tAX-9\tAAAAAAC\tA\t.\t.\t.\n'
)


def _write(path: Path, text: str, gz: bool):
    if gz:
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize('gz,name', [(True, 'bb.vcf.bgz'), (False, 'bb.vcf')])
def test_load_coordinates_reads_vcf(tmp_path, gz, name):
    f = _write(tmp_path / name, BACKBONE, gz)
    rows = load_coordinates(f)
    assert len(rows) == 5
    r = rows[0]
    assert r['chrom'] == 'chr1'
    assert r['pos'] == 930165 and isinstance(r['pos'], int)
    assert r['id'] == 'AX-1,AX-2,AX-3'
    assert r['ref'] == 'A'
    assert r['alt'] == 'C,G,T'
    # lowercase alt is uppercased
    assert rows[1]['alt'] == 'C'


def test_group_families_by_id_first_seen_order(tmp_path):
    rows = load_coordinates(_write(tmp_path / 'bb.vcf', BACKBONE, False))
    fams = group_families(rows)
    assert list(fams.keys()) == ['AX-1,AX-2,AX-3', 'AX-9']
    assert len(fams['AX-1,AX-2,AX-3']) == 4
    assert len(fams['AX-9']) == 1


# --- snv-family TSV loader (id = fmid) ---

SNV_FAMILY_TSV = (
    'chrom\tpos\tref\talt\tfmid\tmsid\tasid\tpsid\n'
    'chr1\t930165\tA\tC,G,T\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'
    'chr1\t930165\tA\tc\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'   # lowercase alt
    'chr1\t930165\tA\tG\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'
    'chr1\t930165\tA\tT\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'
    'chr1\t11080559\tAAAAAAC\tA\tFM-2\tNA\tAffx-9\tAX-9\n'
)


def test_load_snv_family_reads_gzip_tsv_id_from_fmid(tmp_path):
    f = _write(tmp_path / 'sf.tsv.gz', SNV_FAMILY_TSV, True)
    rows = load_snv_family(f)
    assert len(rows) == 5
    r = rows[0]
    assert r['chrom'] == 'chr1'
    assert r['pos'] == 930165 and isinstance(r['pos'], int)
    assert r['id'] == 'FM-1'            # id sourced from the fmid column
    assert r['ref'] == 'A' and r['alt'] == 'C,G,T'
    assert rows[1]['alt'] == 'C'        # lowercase alt uppercased
    # msid/asid/psid are read past, not part of the coordinate record
    assert set(r) == {'chrom', 'pos', 'id', 'ref', 'alt'}


def test_load_snv_family_selects_columns_by_name(tmp_path):
    # reordered columns still parse (select by NAME, not position)
    reordered = 'psid\tfmid\tchrom\tpos\tref\talt\nAX-9\tFM-7\tchr1\t500\tA\tG\n'
    rows = load_snv_family(_write(tmp_path / 'sf.tsv', reordered, False))
    assert rows == [{'chrom': 'chr1', 'pos': 500, 'id': 'FM-7', 'ref': 'A', 'alt': 'G'}]


def test_load_snv_family_rejects_missing_column(tmp_path):
    # a file lacking the fmid column (e.g. a legacy backbone VCF) fails fast
    with pytest.raises(Exception):
        load_snv_family(_write(tmp_path / 'bad.tsv', 'chrom\tpos\tref\talt\nchr1\t1\tA\tG\n', False))


def test_group_families_by_fmid(tmp_path):
    rows = load_snv_family(_write(tmp_path / 'sf.tsv.gz', SNV_FAMILY_TSV, True))
    fams = group_families(rows)
    assert list(fams.keys()) == ['FM-1', 'FM-2']
    assert len(fams['FM-1']) == 4 and len(fams['FM-2']) == 1


class _CountingProvider:
    """Duck-typed depth provider for pure call_families logic tests: serves depths
    from in-memory dicts and counts how many times it is consulted."""

    def __init__(self, autosomes=None, sex=None):
        self._autosomes = autosomes or {}
        self._sex = sex or {}
        self.n_queries = 0

    def autosomes(self, chrom, pos):
        self.n_queries += 1
        return self._autosomes.get((chrom, pos))

    def sex(self, chrom, pos):
        self.n_queries += 1
        return self._sex.get((chrom, pos))


@pytest.fixture
def depth_provider(tmp_path):
    a = _index_depths(
        tmp_path / 'autosomes-depth.tsv',
        'chrom\tpos\tdepth_mean\tn_samples\n'
        'chr1\t100\t12.5\t4\n'
        'chrM\t5\t900.0\t4\n'
    )
    s = _index_depths(
        tmp_path / 'sex-depth.tsv',
        'chrom\tpos\tmean_male\tn_male\tmean_female\tn_female\n'
        'chrX\t100\t5.0\t2\t11.0\t2\n'
        'chrY\t200\t8.0\t2\t\t0\n'
    )
    return DepthProvider(a, s)


@requires_tabix
def test_depth_provider_autosomes_memoized(depth_provider):
    assert depth_provider.autosomes('chr1', 100) == pytest.approx(12.5)
    assert depth_provider.autosomes('chrM', 5) == pytest.approx(900.0)
    assert depth_provider.autosomes('chr1', 999) is None   # missing -> None
    # one tabix query per distinct position; repeats served from cache
    assert depth_provider.n_queries == 3
    depth_provider.autosomes('chr1', 100)
    depth_provider.autosomes('chr1', 999)
    assert depth_provider.n_queries == 3


@requires_tabix
def test_depth_provider_sex(depth_provider):
    rec = depth_provider.sex('chrX', 100)
    assert rec['male'] == pytest.approx(5.0) and rec['female'] == pytest.approx(11.0)
    # blank female column -> None
    rec = depth_provider.sex('chrY', 200)
    assert rec['male'] == pytest.approx(8.0) and rec['female'] is None
    # missing position -> None record
    assert depth_provider.sex('chrX', 999) is None


@requires_tabix
def test_depth_for_dispatch_by_contig_and_gender(depth_provider):
    assert depth_for('chr1', 100, 'male', depth_provider) == pytest.approx(12.5)
    assert depth_for('chrM', 5, 'female', depth_provider) == pytest.approx(900.0)
    assert depth_for('chrX', 100, 'male', depth_provider) == pytest.approx(5.0)
    assert depth_for('chrX', 100, 'female', depth_provider) == pytest.approx(11.0)
    assert depth_for('chrY', 200, 'female', depth_provider) is None   # blank female
    assert depth_for('chr1', 999, 'male', depth_provider) is None     # missing
    assert depth_for('chrX', 999, 'female', depth_provider) is None


def test_call_families_no_depth_query_when_matched_or_observed():
    """A matched family and an observed 0/0 fill must not consult depth — depth is
    fetched only when a site truly needs it to infer the genotype."""
    families = {
        'M': [{'chrom': 'chr1', 'pos': 100, 'id': 'M', 'ref': 'A', 'alt': 'C'}],
        'O': [{'chrom': 'chr1', 'pos': 200, 'id': 'O', 'ref': 'A', 'alt': 'G'}],
    }
    lookup = {match_key('chr1', 100, 'A', 'C'): '0/1'}
    nonvariant = {('chr1', 200): 'homref'}
    provider = _CountingProvider()
    matched, fills = call_families(families, lookup, nonvariant, 'male', provider, min_depth=2)
    assert provider.n_queries == 0
    assert len(matched) == 1
    gt = {line.split('\t')[1]: line.split('\t')[9] for line in fills}
    assert gt['200'] == '0/0'


def test_match_key_altset_order_and_case_invariant():
    assert match_key('chr1', 930165, 'A', 'C,G') == match_key('chr1', 930165, 'A', 'g,c')
    assert match_key('chr1', '930165', 'A', 'C') == match_key('chr1', 930165, 'a', 'c')
    assert match_key('chr1', 930165, 'A', 'C,G') != match_key('chr1', 930165, 'A', 'C')
    assert match_key('chr1', 930165, 'A', 'C') != match_key('chr1', 930165, 'A', 'T')


def test_is_homref_threshold():
    assert is_homref(4.0, 2) is True
    assert is_homref(2.0, 2) is True
    assert is_homref(1.9, 2) is False
    assert is_homref(None, 2) is False
    assert is_homref(math.nan, 2) is False


def test_fill_depth_is_per_member_position():
    """Members of one family can sit at different positions (axiom backbone
    re-normalizes/left-aligns expanded records), so the homref/nocall decision must
    be looked up per member, not once for the whole family."""
    families = {
        'AX-1,AX-2': [
            {'chrom': 'chr1', 'pos': 100, 'id': 'AX-1,AX-2', 'ref': 'A', 'alt': 'C'},
            {'chrom': 'chr1', 'pos': 500, 'id': 'AX-1,AX-2', 'ref': 'A', 'alt': 'G'},
        ]
    }
    # pos 100 homref, pos 500 nocall
    provider = _CountingProvider(autosomes={('chr1', 100): 10.0, ('chr1', 500): 1.0})
    matched, fills = call_families(families, {}, {}, 'male', provider, min_depth=2)
    assert matched == []
    assert provider.n_queries == 2   # one depth consult per unobserved member
    gt = {line.split('\t')[1]: line.split('\t')[9] for line in fills}
    assert gt['100'] == '0/0'
    assert gt['500'] == './.'


def test_build_lookup_splits_variant_and_nonvariant(tmp_path):
    """build_lookup returns the exact-match variant lookup plus a non-variant lookup
    (chrom,pos) -> homref|nocall built from the retained ALT='.' records."""
    from genomics.truth import build_lookup
    vcf = tmp_path / 'prepared.vcf.bgz'
    with gzip.open(vcf, 'wt') as fh:
        fh.write(
            '##fileformat=VCFv4.2\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
            'chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\n'    # variant
            'chr1\t200\t.\tA\t.\t.\t.\t.\tGT\t0/0\n'    # observed homref
            'chr1\t300\t.\tA\t.\t.\t.\t.\tGT\t./.\n'    # observed nocall
        )
    lookup, nonvariant = build_lookup(vcf)
    assert lookup[match_key('chr1', 100, 'A', 'C')] == '0/1'
    # non-variant records do not pollute the variant lookup
    assert match_key('chr1', 200, 'A', '.') not in lookup
    assert nonvariant[('chr1', 200)] == 'homref'
    assert nonvariant[('chr1', 300)] == 'nocall'


def test_call_families_honors_observed_calls():
    """An explicitly observed 0/0 / ./. is honored over the depth-based fill; only a
    site with no sample record falls back to depth."""
    families = {
        'AX-1': [{'chrom': 'chr1', 'pos': 100, 'id': 'AX-1', 'ref': 'A', 'alt': 'C'}],
        'AX-2': [{'chrom': 'chr1', 'pos': 200, 'id': 'AX-2', 'ref': 'A', 'alt': 'G'}],
        'AX-3': [{'chrom': 'chr1', 'pos': 300, 'id': 'AX-3', 'ref': 'A', 'alt': 'T'}],
    }
    nonvariant = {('chr1', 100): 'homref', ('chr1', 200): 'nocall'}
    # depth deliberately DISAGREES with the observed calls; pos 300 absent
    provider = _CountingProvider(autosomes={('chr1', 100): 1.0, ('chr1', 200): 10.0})
    matched, fills = call_families(families, {}, nonvariant, 'male', provider, min_depth=2)
    assert matched == []
    # only pos 300 (unobserved) falls back to depth; observed sites never consult it
    assert provider.n_queries == 1
    gt = {line.split('\t')[1]: line.split('\t')[9] for line in fills}
    assert gt['100'] == '0/0'   # observed homref honored despite low depth
    assert gt['200'] == './.'   # observed nocall honored despite high depth
    assert gt['300'] == './.'   # no observation, no depth -> nocall fill


HEADER = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr1,length=250000000>\n'
    '##contig=<ID=chrX,length=156040895>\n'
    '##contig=<ID=chrY,length=57227415>\n'
    '##contig=<ID=chrM,length=16569>\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
)


@pytest.mark.skipif(shutil.which('bcftools') is None, reason='bcftools not on PATH')
def test_fixploidy_grch38_gender_aware(tmp_path):
    """Fill records get GRCh38 gender-aware ploidy at real non-PAR coordinates via
    the hg38 ploidy resource: male sex-chrom/chrM haploid, female chrY nocall."""
    header_file = tmp_path / 'header.vcf'
    header_file.write_text(HEADER)
    ploidy_file = tmp_path / 'grch38.ploidy'
    ploidy_file.write_text(load_ploidy('hg38'))

    fills = [
        'chr1\t100\ta\tA\tC\t.\t.\t.\tGT\t0/0',       # autosome -> diploid
        'chrX\t3000000\tb\tA\tC\t.\t.\t.\tGT\t0/0',   # non-PAR chrX
        'chrY\t3000000\tc\tA\tC\t.\t.\t.\tGT\t0/0',   # MSY chrY
        'chrM\t50\td\tA\tC\t.\t.\t.\tGT\t0/0',        # mito
        'chrY\t3000001\te\tA\tC\t.\t.\t.\tGT\t./.',   # nocall on MSY
    ]

    (tmp_path / 'm').mkdir()
    (tmp_path / 'f').mkdir()
    male = fixploidy_lines(fills, header_file, 'S1', 'male', ploidy_file, tmp_path / 'm')
    female = fixploidy_lines(fills, header_file, 'S1', 'female', ploidy_file, tmp_path / 'f')

    def gt(lines, i):
        return lines[i].split('\t')[9]

    # male: autosome diploid; non-PAR chrX/chrY and chrM haploid; MSY nocall -> '.'
    assert [gt(male, i) for i in range(5)] == ['0/0', '0', '0', '0', '.']
    # female: chrX diploid (default 2), chrY absent (ploidy 0 -> nocall), chrM haploid
    assert [gt(female, i) for i in range(5)] == ['0/0', '0/0', '.', '0', '.']


def test_load_ploidy_assemblies_and_validation():
    for assembly in ('hg38', 'hg19'):
        text = load_ploidy(assembly)
        assert 'chrX' in text and 'chrY' in text and 'chrM' in text
    # hg38 vs hg19 use different non-PAR chrX coordinates
    assert '2781480' in load_ploidy('hg38')
    assert '2699521' in load_ploidy('hg19')
    with pytest.raises(ValueError):
        load_ploidy('grch38')


@pytest.mark.skipif(shutil.which('bcftools') is None, reason='bcftools not on PATH')
def test_fixploidy_hg19_par_vs_nonpar(tmp_path):
    """The hg19 ploidy resource carves out hg19 PAR: a male PAR1 site stays diploid
    while a non-PAR chrX site collapses to haploid (validates the hg19 coordinates)."""
    header_file = tmp_path / 'header.vcf'
    header_file.write_text(HEADER)
    ploidy_file = tmp_path / 'hg19.ploidy'
    ploidy_file.write_text(load_ploidy('hg19'))

    fills = [
        'chrX\t1000000\ta\tA\tC\t.\t.\t.\tGT\t0/0',   # hg19 PAR1 (60001-2699520) -> diploid
        'chrX\t3000000\tb\tA\tC\t.\t.\t.\tGT\t0/0',   # hg19 non-PAR -> haploid
    ]
    (tmp_path / 'm').mkdir()
    male = fixploidy_lines(fills, header_file, 'S1', 'male', ploidy_file, tmp_path / 'm')
    assert [line.split('\t')[9] for line in male] == ['0/0', '0']


@requires_tabix
@pytest.mark.skipif(shutil.which('bcftools') is None, reason='bcftools not on PATH')
def test_require_single_sample(tmp_path):
    """Each input VCF must contain exactly one sample; multi-sample -> error."""
    def build(name, samples):
        header = (
            '##fileformat=VCFv4.2\n'
            '##contig=<ID=chr1,length=1000000>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            + '\t'.join(samples) + '\n')
        gts = '\t'.join('0/1' for _ in samples)
        plain = tmp_path / f'{name}.vcf'
        plain.write_text(header + f'chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t{gts}\n')
        bgz = tmp_path / f'{name}.vcf.bgz'
        subprocess.run(f'bgzip -c {plain} > {bgz}', shell=True, check=True)
        subprocess.run(['bcftools', 'index', '-f', str(bgz)], check=True)
        return bgz

    assert _require_single_sample(build('one', ['S1'])) == 'S1'
    with pytest.raises(ValueError):
        _require_single_sample(build('two', ['S1', 'S2']))


@requires_tabix
def test_prepare_sample_is_lean_single_pass(tmp_path):
    # The combined input is already single-sample, GT-only, INFO-stripped (from lean
    # extraction), so prepare_sample must not run subset/keep-format/drop-info passes,
    # and trim+normalize must be one streamed pass (no intermediate trimmed VCF).
    if shutil.which('samtools') is None:
        pytest.skip('samtools not on PATH')

    genome = tmp_path / 'genome.fa'
    with genome.open('w') as fh:
        fh.write('>chr1\n' + 'A' * 300 + '\n')
    subprocess.run(['samtools', 'faidx', str(genome)], check=True)

    # single-sample, GT-only, INFO='.'; ALT has an unused allele (C, trimmed) and a
    # lowercase allele (g, kept as-is now that prepare is trim+normalize only);
    # GT uses only allele 1.
    combined = tmp_path / 'combined.vcf'
    combined.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t100\t.\tA\tg,C\t.\t.\t.\tGT\t0/1\n'
    )
    subprocess.run(f'bgzip -f {combined}', shell=True, check=True)
    combined_bgz = tmp_path / 'combined.vcf.bgz'
    (tmp_path / 'combined.vcf.gz').rename(combined_bgz)

    tmp_base = tmp_path / 'work'
    tmp_base.mkdir()
    prepared = prepare_sample(combined_bgz, genome, tmp_base)

    assert not Path(f'{prepared}.csi').exists(), 'prepared VCF should not be indexed'

    # correctness: single sample, ALT trimmed to the used allele (case preserved)
    assert subprocess.run(['bcftools', 'query', '-l', str(prepared)],
                          capture_output=True, text=True, check=True).stdout.split() == ['S1']
    rows = subprocess.run(['bcftools', 'view', '-H', str(prepared)],
                          capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(rows) == 1
    cols = rows[0].split('\t')
    assert cols[3] == 'A'                 # REF
    assert cols[4] == 'g'                 # ALT: C trimmed, g case preserved (no uppercase)
    assert cols[9].split(':')[0] == '0/1'  # GT preserved

    # lean: none of the removed/redundant passes ran (no intermediates written)
    for marker in ('-samples', '-format', '-info', '-trim_alt', '-uppercase'):
        stray = list(tmp_base.rglob(f'*{marker}*.vcf.bgz'))
        assert not stray, f'unexpected intermediate ({marker}): {stray}'
