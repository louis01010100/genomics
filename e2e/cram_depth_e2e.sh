#!/usr/bin/env bash
# Self-contained end-to-end test for `genomics cram-depth`.
#
# Builds the SMALLEST fixture that still exercises every edge-case combination,
# by slicing tiny windows out of real 1000G GRCh38 CRAMs + reference and
# coordinate-shifting the real reads onto a mini-reference, then runs the CLI
# and asserts the outputs.
#
# Edge cases covered (2 male + 2 female samples, ~2 kb per contig):
#   chr1  autosome         -> single mean over all samples      (autosomes-depth.tsv)
#   chrM  mitochondria     -> treated like an autosome          (autosomes-depth.tsv)
#   chrX  non-PAR          -> per-gender; female (diploid) > male (haploid)
#   chrY  MSY              -> per-gender; male present, female ~0 + zero-depth rows
#   decoy (chrUn_...)      -> excluded from both outputs
#   + n_samples/n_male/n_female counts, -aa zero-depth rows,
#     and thread-invariance (identical output for --n-threads 1 vs 4).
#
# Source data is machine-local (~15 GB CRAMs); this is a manual/local e2e.
set -euo pipefail

# ---- source data (real) ----
REF=/affx/louis/tfs/genomics/human/variants/kgp/grch38/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
GENDERS=/affx/louis/tfs/genomics/human/variants/kgp/samples/genders.tsv
CRAMDIR=/home/louis/tfs/genomics/human/variants/kgp/grch38/align/raw

# ---- code under test ----
REPO=/affx/louis/tfs/projects/personal/github/genomics
SRC=${GENOMICS_SRC:-$REPO/wt/feat-cram-depth/src}
PY=${PYTHON:-$REPO/wt/main/.venv/bin/python}
command -v samtools >/dev/null 2>&1 || export PATH="/home/louis/applications/samtools:$PATH"
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }

WORK=$(mktemp -d /tmp/cram_depth_e2e.XXXXXX)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# sample gender cram   (2 male + 2 female with locally-available CRAMs)
SAMPLES=(
  "HG00403 male   $CRAMDIR/HG00403.final.cram"
  "HG00406 male   $CRAMDIR/HG00406.final.cram"
  "HG00404 female $CRAMDIR/HG00404.final.cram"
  "HG00419 female $CRAMDIR/HG00419.final.cram"
)
# newname contig start end subsample(1=all)
WINDOWS=(
  "chr1                       chr1                       1000001 1002000 1"
  "chrX                       chrX                       5000001 5002000 1"
  "chrY                       chrY                       6000001 6002000 1"
  "chrM                       chrM                       1       2000    0.003"
  "chrUn_JTFH01000960v1_decoy chrUn_JTFH01000960v1_decoy 1       1000    1"
)

echo "### building minimal fixture from real data ..."
: > miniref.fa
for w in "${WINDOWS[@]}"; do
  read -r name contig start end sub <<<"$w"
  echo ">$name" >> miniref.fa
  samtools faidx "$REF" "$contig:$start-$end" | tail -n +2 >> miniref.fa
done
samtools faidx miniref.fa

: > crams.tsv; printf 'sample\tcram\n' >> crams.tsv
for row in "${SAMPLES[@]}"; do
  read -r sample gender cram <<<"$row"
  printf '%s\t%s/%s.cram\n' "$sample" "$WORK" "$sample" >> crams.tsv
  sam="$sample.sam"
  { printf '@HD\tVN:1.6\tSO:unsorted\n'
    awk '{print "@SQ\tSN:"$1"\tLN:"$2}' miniref.fa.fai; } > "$sam"
  for w in "${WINDOWS[@]}"; do
    read -r name contig start end sub <<<"$w"
    off=$((start - 1)); sflag=(); [ "$sub" != "1" ] && sflag=(-s "$sub")
    samtools view "${sflag[@]}" -T "$REF" "$cram" "$contig:$start-$end" \
      | awk -v OFS='\t' -v start="$start" -v endw="$end" -v off="$off" -v name="$name" '
          { pos=$4; cig=$6; rlen=0; num="";
            for (i=1;i<=length(cig);i++){ c=substr(cig,i,1);
              if (c ~ /[0-9]/) num=num c;
              else { if (c=="M"||c=="D"||c=="N"||c=="="||c=="X") rlen+=num+0; num="" } }
            rend=pos+rlen-1;
            if (pos < start || rend > endw) next;      # keep reads fully inside window
            $3=name; $4=pos-off; $7="*"; $8=0; $9=0; print }' >> "$sam"
  done
  samtools sort -O cram --reference miniref.fa -o "$sample.cram" "$sam" 2>/dev/null
  samtools index "$sample.cram"
done

run_cli () { # out_dir n_threads
  PYTHONPATH="$SRC" "$PY" -m genomics cram-depth \
    --crams-file crams.tsv --genders-file "$GENDERS" \
    --genome-file miniref.fa --output-dir "$1" --n-threads "$2" 2>"$1.err" \
    || { echo "CLI FAILED (n-threads=$2)"; cat "$1.err"; exit 1; }
}

echo "### running cram-depth (n-threads=4, and n-threads=1 for invariance) ..."
run_cli "$WORK/out"  4
run_cli "$WORK/out1" 1
A="$WORK/out/autosomes-depth.tsv"; S="$WORK/out/sex-depth.tsv"

echo "================ autosomes-depth.tsv (head) ================"; head -3 "$A"
echo "================ sex-depth.tsv (head) ================";       head -3 "$S"

fail=0
assert () { if eval "$2"; then echo "PASS: $1"; else echo "FAIL: $1"; fail=1; fi; }

# --- files / cleanup ---
assert "exactly autosomes-depth.tsv, depth.log, sex-depth.tsv produced" \
  '[ "$(ls -1 "$WORK/out" | sort | tr "\n" ,)" = "autosomes-depth.tsv,depth.log,sex-depth.tsv," ]'
assert "no per-sample per-position depth artifacts left" \
  '[ ! -d "$WORK/out/depths" ] && ! ls "$WORK/out"/*.tsv 2>/dev/null | grep -qvE "autosomes-depth.tsv|sex-depth.tsv"'
# --- schemas ---
assert "autosomes header exact" 'head -1 "$A" | grep -qx "chrom	pos	depth_mean	n_samples"'
assert "sex header exact"       'head -1 "$S" | grep -qx "chrom	pos	mean_male	n_male	mean_female	n_female"'
# --- classification: autosome + mito -> autosomes file ---
assert "autosomes file contains exactly chr1 and chrM" \
  '[ "$(awk -F"\t" "NR>1{print \$1}" "$A" | sort -u | tr "\n" ,)" = "chr1,chrM," ]'
assert "chr1 has 2000 rows (whole window via -aa)" '[ "$(awk -F"\t" "\$1==\"chr1\"" "$A" | wc -l)" -eq 2000 ]'
assert "chrM has 2000 rows"                        '[ "$(awk -F"\t" "\$1==\"chrM\"" "$A" | wc -l)" -eq 2000 ]'
assert "autosomes n_samples==4 on every row"       '! awk -F"\t" "NR>1 && \$4!=4{print}" "$A" | grep -q .'
# --- mito is treated like an autosome, never a sex chromosome ---
assert "chrM present in autosomes file" 'grep -qP "^chrM\t" "$A"'
assert "chrM absent from sex file"      '! grep -qP "^chrM\t" "$S"'
# --- sex chromosomes only in sex file, per-gender counts ---
assert "sex file contains exactly chrX and chrY" \
  '[ "$(awk -F"\t" "NR>1{print \$1}" "$S" | sort -u | tr "\n" ,)" = "chrX,chrY," ]'
assert "no chrX/chrY leak into autosomes file" '! grep -qP "^(chrX|chrY)\t" "$A"'
assert "sex n_male==2 and n_female==2 on every row" \
  '! awk -F"\t" "NR>1 && (\$4!=2 || \$6!=2){print}" "$S" | grep -q .'
# --- gender x chrX: female diploid > male haploid ---
assert "chrX: avg female depth > 1.5x avg male depth (diploid vs haploid)" \
  'awk -F"\t" "\$1==\"chrX\"{m+=\$3; f+=\$5; n++} END{exit !(f > 1.5*m && m>0)}" "$S"'
# --- gender x chrY: male present, female ~0, zero-depth rows present ---
assert "chrY: avg male depth > 5x avg female depth (female lacks chrY)" \
  'awk -F"\t" "\$1==\"chrY\"{m+=\$3; f+=\$5; n++} END{exit !(m > 5*f && m>5)}" "$S"'
assert "chrY: some female positions are zero-depth (-aa emits zeros)" \
  'awk -F"\t" "\$1==\"chrY\" && \$5==0{c++} END{exit !(c>0)}" "$S"'
# --- exclusion ---
assert "decoy contig excluded from both files" '! grep -q "decoy" "$A" "$S"'
# --- parallelism determinism ---
assert "output byte-identical for --n-threads 1 vs 4 (autosomes)" 'cmp -s "$A" "$WORK/out1/autosomes-depth.tsv"'
assert "output byte-identical for --n-threads 1 vs 4 (sex)"       'cmp -s "$S" "$WORK/out1/sex-depth.tsv"'

echo "================ RESULT ================"
[ $fail -eq 0 ] && echo "E2E: ALL PASS" || { echo "E2E: FAILURES"; exit 1; }
