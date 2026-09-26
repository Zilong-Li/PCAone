#!/usr/bin/env bash
# Reproduce the --evaladmix benchmark in docs/evaladmix.md.
#
#   bash run_all.sh [workdir] [path/to/PCAone]
#
# Fetches the data and the comparison software, runs every method, and leaves
# the outputs in $WORK/out for benchmark.R to read. Nothing is simulated: the
# dataset ships with relateAdmix, and step 6 only masks some of its calls.
set -euo pipefail

WORK="${1:-$PWD/evaladmix-benchmark}"
PCAONE="${2:-$(command -v PCAone || echo "$PWD/PCAone")}"
THREADS="${THREADS:-8}"
HERE="$(cd "$(dirname "$0")" && pwd)"   # before the cd below: $0 may be relative

[ -x "$PCAONE" ] || { echo "PCAone not found at $PCAONE -- pass it as argument 2" >&2; exit 1; }
# both may be relative, as in the docs (./PCAone); the script cds below
PCAONE="$(cd "$(dirname "$PCAONE")" && pwd)/$(basename "$PCAONE")"
mkdir -p "$WORK"/{src,out}
WORK="$(cd "$WORK" && pwd)"
cd "$WORK/src"
echo "==> PCAone: $PCAONE"

# ---------------------------------------------------------------- 1. data ---
# 126 individuals, 104,290 autosomal SNPs, no missing genotypes, simulated from
# two source populations with a known pedigree. smallPlink.2.{P,Q} are ADMIXTURE
# output for K=2 shipped with the data, so ADMIXTURE itself is never run and no
# method is advantaged by a better admixture fit.
[ -d relateAdmix ] || git clone --quiet https://github.com/aalbrechtsen/relateAdmix.git
DATA="$WORK/src/relateAdmix/data"
echo "==> data: $DATA/smallPlink.{bed,bim,fam} + smallPlink.2.{P,Q}"

# ------------------------------------------------------- 2. RelateAdmix ML ---
# its estimates for unrelated pairs depend on -P: THREADS=2 gives an RMSE of
# 0.00046 on them instead of 0.00054 (related pairs are unaffected)
if [ ! -x relateAdmix/src/relateAdmix ]; then
  ( cd relateAdmix/src && { [ -f Makefile ] || cp CPP_Makefile Makefile; } && make -s -j"$THREADS" )
fi
( cd "$DATA" && ../src/relateAdmix -plink smallPlink -f smallPlink.2.P \
      -q smallPlink.2.Q -P "$THREADS" >/dev/null 2>&1 )
cp "$DATA/output.k" "$WORK/out/relateadmix.k"
echo "==> RelateAdmix done"

# ----------------------------------------------------------- 3. evalAdmix ---
# v1.0 or later, for both -method evalAdmix (EM) and -method corrected
[ -d evalAdmix ] || git clone --quiet https://github.com/GenisGE/evalAdmix.git
[ -x evalAdmix/evalAdmix ] || ( cd evalAdmix && make -s -j"$THREADS" )
EA="$WORK/src/evalAdmix/evalAdmix"
"$EA" -plink "$DATA/smallPlink" -fname "$DATA/smallPlink.2.P" -qname "$DATA/smallPlink.2.Q" \
      -P "$THREADS" -o "$WORK/out/evaladmix_em.corres" >/dev/null 2>&1
"$EA" -plink "$DATA/smallPlink" -fname "$DATA/smallPlink.2.P" -qname "$DATA/smallPlink.2.Q" \
      -method corrected -P "$THREADS" -o "$WORK/out/evaladmix_proj.corres" >/dev/null 2>&1
echo "==> evalAdmix (EM + projection) done"

# ------------------------------------------- 4. evalPopStructure R reference ---
[ -d evalPopStructure ] || git clone --quiet https://github.com/popgenDK/evalPopStructure.git
echo "==> evalPopStructure R reference fetched"

# ------------------------------------------------- 5. PCAone --evaladmix -----
# K=2 -> K-1 = 1 PC. -d 0 (IRAM) for a deterministic run.
"$PCAONE" -b "$DATA/smallPlink" -k 1 -d 0 --evaladmix --maf 0.05 \
      -n "$THREADS" -o "$WORK/out/pcaone" >/dev/null 2>&1
# multi-column .eigvecs projected on its first column: must equal the run above
"$PCAONE" -b "$DATA/smallPlink" -k 4 -d 0 --evaladmix --evaladmix-k 1 --maf 0.05 \
      -n "$THREADS" -o "$WORK/out/pcaone_k4" >/dev/null 2>&1
# out-of-core (PCAone rejects --maf out-of-core, so this one runs on all sites)
"$PCAONE" -b "$DATA/smallPlink" -k 1 -d 0 --evaladmix -m 0.002 \
      -n "$THREADS" -o "$WORK/out/pcaone_ooc" >/dev/null 2>&1
"$PCAONE" -b "$DATA/smallPlink" -k 1 -d 0 --evaladmix \
      -n "$THREADS" -o "$WORK/out/pcaone_ic"  >/dev/null 2>&1
echo "==> PCAone --evaladmix done"

# ------------------------------------------------------ 6. missing genotypes ---
# The same data with calls set to missing: at random (5, 10, 20%), at a rate
# that varies by sample (0-40%), and by batch (two random halves of the samples,
# each missing its own 25% of the sites). evalAdmix skips a site for a pair
# unless both are genotyped; PCAone imputes and rescales each pair by the sites
# both are genotyped at. Same .P/.Q and K-1 = 1 PC as above.
MISS="$WORK/out/missing"; mkdir -p "$MISS"
seed=1
for design in mcar0.05 mcar0.1 mcar0.2 varying batch; do
  seed=$((seed + 1))
  [ -f "$MISS/$design.bed" ] || Rscript "$HERE/make_missing.R" "$DATA/smallPlink" "$MISS/$design" "$design" "$seed"
  "$PCAONE" -b "$MISS/$design" -k 1 -d 0 --evaladmix --maf 0.05 -n "$THREADS" -o "$MISS/pcaone_$design" >/dev/null 2>&1
  [ -f "$MISS/evaladmix_em_$design.corres" ] ||
    "$EA" -plink "$MISS/$design" -fname "$DATA/smallPlink.2.P" -qname "$DATA/smallPlink.2.Q" \
          -P "$THREADS" -o "$MISS/evaladmix_em_$design.corres" >/dev/null 2>&1
done
echo "==> missing-genotype runs done"

echo
echo "outputs in $WORK/out. now run:"
echo "  Rscript $HERE/benchmark.R $WORK $DATA"
