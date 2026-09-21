# Two robustness bugs in the LD path: NaN R² and a segfault on a missing `.mbim`

Both were found while validating the [`read_usv` fix](read-usv-fix.md) and are
independent of it.

## 1. Zero-variance variants produce `NaN` R², which spreads

`ld_r2_big()`, `ld_prune_big()` and the clumping path all begin with

```cpp
Arr1D sds = 1.0 / calc_sds(G);
```

A variant with **no residual variance** — monomorphic, or with genotypes lying
entirely in the span of the PCs being adjusted for — has `sd = 0`, so this is
`inf`. The correlation is then

```cpp
double r = G.col(i).dot(G.col(k)) * (sds(i) * sds(k) * df);
```

and since the dot product with a zero column is `0`, the result is `0 * inf`,
i.e. `NaN`. It is written to `.ld.gz` verbatim as `-nan`, for **every pair**
involving that variant.

`calc_sds()` already carried the comment `// return 1e-9 when sd is 0`, so the
case was anticipated — the guard was just never written.

On the 126-sample test data, 4 of 104,290 variants are monomorphic (AF = 0, all
126 samples the same homozygote). Each falls in a window with ~61 neighbours, so
those 4 variants produced **244 `-nan` rows**. PCAone does warn `sites with
MAF=0 found! remove them first!` during reading, but then proceeds and emits the
NaNs anyway.

`-nan` in a numeric column breaks naive downstream parsing, and R's `as.numeric`
turns it into `NA` that silently propagates through any summary statistic.

### The fix

A new `calc_inv_sds()` returns **0** for zero-variance columns instead of `inf`,
so `r` comes out as exactly `0`:

```cpp
Arr1D calc_inv_sds(const Mat2D& X, int& nzero) {
  const Arr1D sds = calc_sds(X);
  ...
  if (sds(i) > 1e-9) inv(i) = 1.0 / sds(i);
  else { inv(i) = 0.0; ++nzero; }
}
```

Reporting 0 states what is actually true: a variant with no variance carries no
LD information. All three call sites use it and emit one shared warning naming
the count.

### Effect

```
4 of 104290 variants have no variance after ancestry adjustment (monomorphic, or
fully explained by the PCs). their R2 is reported as 0 rather than NaN.
consider --maf to remove them.
```

| | before | after |
|---|---|---|
| `-nan` rows in `.ld.gz` | 244 | **0** |
| all other R² values | — | **bit-identical** (4,611,432 pairs, max diff 0) |

## 2. A missing `.mbim` segfaults instead of reporting the problem

`-P/--USV <prefix>` sets `filebim = <prefix>.mbim`, and `Data::prepare()` calls
`read_frq()` on it. `read_frq()` opened the stream without checking:

```cpp
std::ifstream fin(path);
std::string line;
while (getline(fin, line)) { ... }
return Eigen::Map<Mat1D>(V.data(), V.size());
```

A nonexistent file simply reads as empty, `F` stays size 0, and the first
`F(snp_idx)` downstream indexes out of bounds — **segfault, exit 139, no
message**.

This is easy to hit, because **only a run with `-D`/`--ld` writes a `.mbim`**. A
prefix from a plain PCA run looks perfectly valid and has the `.eigvecs`,
`.sigvals` and `.loadings` the flag advertises, but no `.mbim`:

```bash
PCAone -b plink -k 2 -o pcs                             # no .mbim written
PCAone -b plink -P pcs -R --ld-bp 1000000 -o ld         # segfault, before the fix
```

### The fix

Check the stream, and check that something was parsed:

```cpp
if (!fin.is_open()) cao.error("can not open the allele frequency file\n => " + path);
...
if (V.empty()) cao.error("no allele frequencies found in\n => " + path);
```

The same command now names the missing file and exits cleanly. Use `-D` when
generating the prefix you intend to pass to `-P`.

## Regression

All other behaviour is unchanged: `--evaladmix` output is bit-identical before
and after, plain PCA and `-D` runs are unaffected, and every non-NaN R² value in
the 4.6M-pair comparison is unchanged to the last bit.
