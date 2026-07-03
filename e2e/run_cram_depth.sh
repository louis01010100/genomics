#!/usr/bin/env bash
# E2E test: run `genomics cram-depth` against the persisted e2e fixture (no rebuild)
# and assert every edge-case combination + thread-invariance. Exits non-zero on failure.
#
# Usage:   ./run_cram_depth.sh [output-dir]
# Env:     GENOMICS_SRC=<path>   PYTHON=<python>   N_THREADS=<n>
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FIXTURE="$HERE/cram_depth_e2e_fixture"

# Walk up to the repo container (holds .bare and wt/) to locate the code under test.
container="$HERE"
while [ "$container" != "/" ] && ! { [ -d "$container/.bare" ] && [ -d "$container/wt" ]; }; do
    container="$(dirname "$container")"
done

SRC="${GENOMICS_SRC:-}"
if [ -z "$SRC" ]; then
    for name in feat-cram-depth main; do
        if [ -d "$container/wt/$name/src/genomics" ]; then SRC="$container/wt/$name/src"; break; fi
    done
fi
[ -n "$SRC" ] || { echo "could not locate genomics source (set GENOMICS_SRC)" >&2; exit 1; }

PY="${PYTHON:-$container/wt/main/.venv/bin/python}"
command -v samtools >/dev/null 2>&1 || export PATH="/home/louis/applications/samtools:$PATH"
command -v samtools >/dev/null 2>&1 || { echo "samtools not found" >&2; exit 1; }

OUT="${1:-$HERE/workspace}"
rm -rf "$OUT"

work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
crams="$work/crams.tsv"
printf 'sample\tcram\n' > "$crams"
for c in "$FIXTURE"/*.cram; do
    printf '%s\t%s\n' "$(basename "$c" .cram)" "$c" >> "$crams"
done

echo "genomics source : $SRC"
echo "fixture         : $FIXTURE"
echo "output dir      : $OUT"
echo

run_cli () { # out_dir n_threads
    PYTHONPATH="$SRC" "$PY" -m genomics cram-depth \
        --crams-file "$crams" \
        --genders-file "$FIXTURE/genders.tsv" \
        --genome-file "$FIXTURE/miniref.fa" \
        --output-dir "$1" \
        --n-threads "$2"
}

threads="${N_THREADS:-4}"
run_cli "$OUT" "$threads"
run_cli "$work/out1" 1          # second run for the thread-invariance check

echo
echo "=== $OUT ==="; ls -1 "$OUT"
echo "=== autosomes-depth.tsv (head) ==="; head -5 "$OUT/autosomes-depth.tsv"
echo "=== sex-depth.tsv (head) ==="; head -5 "$OUT/sex-depth.tsv"
echo

A="$OUT/autosomes-depth.tsv"; S="$OUT/sex-depth.tsv"
fail=0
assert () { if eval "$2"; then echo "PASS: $1"; else echo "FAIL: $1"; fail=1; fi; }

# files / cleanup
assert "exactly autosomes-depth.tsv, depth.log, sex-depth.tsv produced" \
  '[ "$(ls -1 "$OUT" | sort | tr "\n" ,)" = "autosomes-depth.tsv,depth.log,sex-depth.tsv," ]'
assert "no per-sample per-position depth artifacts left" \
  '[ ! -d "$OUT/depths" ] && ! ls "$OUT"/*.tsv 2>/dev/null | grep -qvE "autosomes-depth.tsv|sex-depth.tsv"'
# schemas
assert "autosomes header exact" 'head -1 "$A" | grep -qx "chrom	pos	depth_mean	n_samples"'
assert "sex header exact"       'head -1 "$S" | grep -qx "chrom	pos	mean_male	n_male	mean_female	n_female"'
# classification: autosome + mito -> autosomes file
assert "autosomes file contains exactly chr1 and chrM" \
  '[ "$(awk -F"\t" "NR>1{print \$1}" "$A" | sort -u | tr "\n" ,)" = "chr1,chrM," ]'
assert "chr1 has 2000 rows (whole window via -aa)" '[ "$(awk -F"\t" "\$1==\"chr1\"" "$A" | wc -l)" -eq 2000 ]'
assert "chrM has 2000 rows"                        '[ "$(awk -F"\t" "\$1==\"chrM\"" "$A" | wc -l)" -eq 2000 ]'
assert "autosomes n_samples==4 on every row"       '! awk -F"\t" "NR>1 && \$4!=4{print}" "$A" | grep -q .'
# mito treated like an autosome, never a sex chromosome
assert "chrM present in autosomes file" 'grep -qP "^chrM\t" "$A"'
assert "chrM absent from sex file"      '! grep -qP "^chrM\t" "$S"'
# sex chromosomes only in sex file, per-gender counts
assert "sex file contains exactly chrX and chrY" \
  '[ "$(awk -F"\t" "NR>1{print \$1}" "$S" | sort -u | tr "\n" ,)" = "chrX,chrY," ]'
assert "no chrX/chrY leak into autosomes file" '! grep -qP "^(chrX|chrY)\t" "$A"'
assert "sex n_male==2 and n_female==2 on every row" \
  '! awk -F"\t" "NR>1 && (\$4!=2 || \$6!=2){print}" "$S" | grep -q .'
# gender x chrX: female (diploid) > male (haploid)
assert "chrX: avg female depth > 1.5x avg male depth" \
  'awk -F"\t" "\$1==\"chrX\"{m+=\$3; f+=\$5} END{exit !(f > 1.5*m && m>0)}" "$S"'
# gender x chrY: male present, female ~0, zero-depth rows present
assert "chrY: avg male depth > 5x avg female depth" \
  'awk -F"\t" "\$1==\"chrY\"{m+=\$3; f+=\$5} END{exit !(m > 5*f && m>5)}" "$S"'
assert "chrY: some female positions are zero-depth (-aa emits zeros)" \
  'awk -F"\t" "\$1==\"chrY\" && \$5==0{c++} END{exit !(c>0)}" "$S"'
# exclusion
assert "decoy contig excluded from both files" '! grep -q "decoy" "$A" "$S"'
# parallelism determinism
assert "output byte-identical for --n-threads 1 vs $threads (autosomes)" 'cmp -s "$A" "$work/out1/autosomes-depth.tsv"'
assert "output byte-identical for --n-threads 1 vs $threads (sex)"       'cmp -s "$S" "$work/out1/sex-depth.tsv"'

echo "================ RESULT ================"
[ $fail -eq 0 ] && echo "E2E: ALL PASS" || { echo "E2E: FAILURES"; exit 1; }
