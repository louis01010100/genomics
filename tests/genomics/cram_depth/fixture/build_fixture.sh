#!/usr/bin/env bash
# Regenerate the small real-derived e2e fixture for `genomics cram-depth`.
#
# Slices tiny windows out of real 1000G GRCh38 CRAMs + reference, coordinate-shifts
# the real reads onto a mini-reference, and writes small CRAMs. The windows are
# chosen so the fixture exercises every classification / gender edge case:
#   chr1  (autosome)      -> single mean over all samples, autosomes-depth.tsv
#   chrM  (mitochondria)  -> treated like autosome, autosomes-depth.tsv
#   chrX  (non-PAR)       -> per-gender; female (diploid) > male (haploid)
#   chrY  (MSY)           -> per-gender; male present, female ~0 (+ zero-depth rows)
#   chrUn_...decoy        -> excluded from both outputs
#
# Source data is machine-local (~15GB CRAMs); this is a regenerate tool, not run in CI.
# The committed fixture (miniref.fa, *.cram) is what the e2e test consumes.
set -euo pipefail
export PATH="/home/louis/applications/samtools:$PATH"

REF=/affx/louis/tfs/genomics/human/variants/kgp/grch38/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
CRAMDIR=/home/louis/tfs/genomics/human/variants/kgp/grch38/align/raw
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

# sample gender cram
SAMPLES=(
  "HG00403 male   $CRAMDIR/HG00403.final.cram"
  "HG00406 male   $CRAMDIR/HG00406.final.cram"
  "HG00404 female $CRAMDIR/HG00404.final.cram"
  "HG00419 female $CRAMDIR/HG00419.final.cram"
)

# newname contig start end subsample
WINDOWS=(
  "chr1                       chr1                       1000001 1002000 1"
  "chrX                       chrX                       5000001 5002000 1"
  "chrY                       chrY                       6000001 6002000 1"
  "chrM                       chrM                       1       2000    0.003"
  "chrUn_JTFH01000960v1_decoy chrUn_JTFH01000960v1_decoy 1       1000    1"
)

# ---- mini-reference (rename each slice to its target contig) ----
: > miniref.fa
for w in "${WINDOWS[@]}"; do
  read -r name contig start end sub <<<"$w"
  echo ">$name" >> miniref.fa
  samtools faidx "$REF" "$contig:$start-$end" | tail -n +2 >> miniref.fa
done
samtools faidx miniref.fa

# ---- per-sample mini CRAMs (real reads, shifted onto mini-reference) ----
: > genders.tsv
printf 'sample\tgender\n' >> genders.tsv
for row in "${SAMPLES[@]}"; do
  read -r sample gender cram <<<"$row"
  printf '%s\t%s\n' "$sample" "$gender" >> genders.tsv

  sam="$sample.sam"
  { printf '@HD\tVN:1.6\tSO:unsorted\n'
    awk '{print "@SQ\tSN:"$1"\tLN:"$2}' miniref.fa.fai; } > "$sam"

  for w in "${WINDOWS[@]}"; do
    read -r name contig start end sub <<<"$w"
    off=$((start - 1))
    sflag=(); [ "$sub" != "1" ] && sflag=(-s "$sub")
    samtools view "${sflag[@]}" -T "$REF" "$cram" "$contig:$start-$end" \
      | awk -v OFS='\t' -v start="$start" -v endw="$end" -v off="$off" -v name="$name" '
          {
            pos=$4; cig=$6; rlen=0; num="";
            for (i=1;i<=length(cig);i++) {
              c=substr(cig,i,1);
              if (c ~ /[0-9]/) { num=num c }
              else { if (c=="M"||c=="D"||c=="N"||c=="="||c=="X") rlen+=num+0; num="" }
            }
            rend=pos+rlen-1;
            if (pos < start || rend > endw) next;   # keep reads fully inside window
            $3=name; $4=pos-off; $7="*"; $8=0; $9=0;
            print
          }' >> "$sam"
  done

  samtools sort -O cram --reference miniref.fa -o "$sample.cram" "$sam"
  samtools index "$sample.cram"
  rm -f "$sam"
done

echo "=== fixture built in $HERE ==="
ls -lh miniref.fa *.cram
echo "=== mini header (contigs) ==="
samtools view -H HG00403.cram | awk '$1=="@SQ"{print $2, $3}'
