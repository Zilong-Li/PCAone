#!/bin/sh
# Invalid option values must stop PCAone with a message naming the option,
# before any input is read. No data is needed: x.bed etc. never get opened.
set -u

cd "$(dirname "$0")/.."
PCAONE=${PCAONE:-./PCAone}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0

# bad <text expected in the error> <args...>
bad() {
  msg=$1
  shift
  if out=$("$PCAONE" "$@" -o "$tmp/out" 2>&1); then
    echo "FAIL (accepted): $*"
    fail=1
  elif ! printf '%s\n' "$out" | grep -qF -- "$msg"; then
    echo "FAIL (wrong error): $*"
    printf '%s\n' "$out" | sed 's/^/    /'
    fail=1
  fi
}

# input
bad "no input file" -k 3
bad "only one of -b/--bfile" -b x -c y
bad "unexpected argument: 3" -b x 3
# unsigned options refuse a sign instead of wrapping to 4294967295
bad "--maxp: '-1'" -b x --maxp -1
bad "-n: '-2'" -b x -n -2
# ranges
bad "--maxp" -b x --maxp 0
bad "-n/--threads" -b x -n 0
bad "-v/--verbose" -b x -v 4
bad "-m/--memory" -b x -m -1
bad "-k/--pc" -b x -k 0
bad "-d/--svd" -b x -d 7
bad "-C/--scale" -b x -C 5
bad "-w/--batches" -b x -w 12
bad "-w/--batches" -b x -w 2
bad "--scale-factor" -c x --scale-factor 0
bad "--buffer" -b x --buffer 0
bad "--imaxiter" -b x --imaxiter 0
bad "--itol" -b x --itol 0
bad "--ncv" -b x -k 10 --ncv 5
bad "--rand" -b x --rand 2
bad "--tol-rsvd" -b x --tol-rsvd 0
bad "--tol-em" -b x --tol-em -1
bad "--tol-maf" -b x --tol-maf 0
bad "--emu cannot be used" -b x --emu --pcangsd
bad "--emu cannot be used" -G x --emu
bad "--maf" -b x --maf -0.1
bad "--maf" -b x --maf 0.5
bad "--project" -b x --project -1
bad "--project" -b x --project 4
bad "--project-bootstrap needs" -b x -P p --project 2 --project-bootstrap 1
bad "--project-bootstrap-save requires" -b x -P p --project 2 --project-bootstrap-save
bad "--inbreed" -b x --inbreed 2
bad "--selection" -b x --selection 3
bad "--evaladmix-k must be" -b x --evaladmix --evaladmix-k -1
bad "--evaladmix-k requires" -b x --evaladmix-k 2
bad "--evaladmix-k cannot" -b x -k 3 --evaladmix --evaladmix-k 4
bad "--ld-r2" -b x --ld-r2 1.5
bad "--ld-r2" -b x --ld-r2 -0.2
bad "--ld-bp" -b x --ld-bp 0
bad "--ld-stats" -b x --ld-stats 2
bad "--clump-p1 cannot" -b x --clump-p1 0.05 --clump-p2 0.01
bad "--clump-p1" -b x --clump-p1 0
bad "--clump-p2" -b x --clump-p2 1.5
bad "--clump-r2" -b x --clump-r2 0
bad "--clump-bp" -b x --clump-bp 0
bad "--clump-names" -b x --clump-names CHR,BP
bad "--clump-names" -b x --clump-names CHR,,P
# input types a mode cannot read (these crashed or ran something else)
bad "--pcangsd supports only" -p x --pcangsd
bad "--pcangsd supports only" -g x --pcangsd
bad "--pcangsd supports only" -c x --pcangsd
bad "--emu supports only" -c x --emu
bad "--emu with --bgen" -g x --emu -m 1
bad "--project supports only" -g x -P p --project 1
bad "--project supports only" -c x -P p --project 1
bad "--project 3 requires" -b x -P p --project 3
bad "--evaladmix supports only" -g x --evaladmix
bad "--evaladmix supports only" -c x --evaladmix
bad "--inbreed supports only" -g x -P p --inbreed 1

if [ "$fail" -eq 0 ]; then
  echo "SUCCESS: invalid options are rejected."
fi
exit "$fail"
