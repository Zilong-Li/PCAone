#!/usr/bin/env bash
# Reproduce the projection benchmark in docs/biobank/projection.md: PCAone
# --project 1/2 against plink2 --score on a simulated biobank-sized target.
#
#   bash projection_vs_plink2.sh [workdir] [path/to/PCAone] [path/to/plink2]
#
# 205,000 samples and 100,000 variants with 1% missing calls (plink2 --dummy):
# the PCs of the first 5,000 (variants with MAF >= 1% in them, 20 PCs) and the
# other 200,000 projected onto them, from a BED of 4.6 GB and its PGEN. Each
# tool scores the target on its own reference PCA. The files are read once
# before the runs, so the times are those of the page cache.
set -euo pipefail
export LC_ALL=C   # decimal points from time and awk

WORK="${1:-$PWD/projection-benchmark}"
PCAONE="${2:-$(command -v PCAone || echo "$PWD/PCAone")}"
PLINK2="${3:-$(command -v plink2)}"
THREADS="${THREADS:-8}"
K=20
mkdir -p "$WORK" && cd "$WORK"

if [ ! -f tgt.pgen ]; then
  "$PLINK2" --dummy 205000 100000 0.01 acgt --seed 7 --threads "$THREADS" --make-bed --out all > /dev/null
  head -5000 all.fam | awk '{print $1, $2}' > ref.keep
  "$PLINK2" --bfile all --keep ref.keep --maf 0.01 --make-bed --out ref --threads "$THREADS" > /dev/null
  "$PLINK2" --bfile all --remove ref.keep --extract ref.bim --make-bed --out tgt --threads "$THREADS" > /dev/null
  "$PLINK2" --bfile tgt --make-pgen --out tgt --threads "$THREADS" > /dev/null
  rm -f all.*
fi
"$PCAONE" -b ref -k $K -d 3 -n "$THREADS" -o pcaone_ref > /dev/null
"$PLINK2" --bfile ref --freq counts --pca $K allele-wts --threads "$THREADS" --out plink2_ref > /dev/null

# wall seconds and peak memory (GB) of a command
measure() {
  local name=$1; shift
  if [ "$(uname)" = Darwin ]; then
    /usr/bin/time -l "$@" > /dev/null 2> time.log
    awk -v n="$name" '/ real /{r=$1} /maximum resident/{m=$1/2^30} END{printf "%-24s %8.2f s %6.2f GB\n", n, r, m}' time.log
  else
    /usr/bin/time -f "%e %M" "$@" > /dev/null 2> time.log
    tail -1 time.log | awk -v n="$name" '{printf "%-24s %8.2f s %6.2f GB\n", n, $1, $2/2^20}'
  fi
}
cat tgt.bed tgt.pgen > /dev/null
NC=$(head -1 plink2_ref.eigenvec.allele | awk '{print NF}')   # #CHROM ID REF ALT [PROVISIONAL_REF?] A1 PC1..
A1=$((NC - K))
for f in bfile pfile; do
  measure "plink2 --score ($f)" "$PLINK2" --$f tgt --read-freq plink2_ref.acount \
    --score plink2_ref.eigenvec.allele 2 $A1 header-read no-mean-imputation variance-standardize \
    --score-col-nums $((A1 + 1))-$NC --threads "$THREADS" --out plink2_$f
done
for m in 1 2; do
  measure "PCAone --project $m (bed)" "$PCAONE" -b tgt -P pcaone_ref --project $m -n "$THREADS" -o pcaone_p${m}_bed
  measure "PCAone --project $m (pgen)" "$PCAONE" -p tgt -P pcaone_ref --project $m -n "$THREADS" -o pcaone_p${m}_pgen
done
