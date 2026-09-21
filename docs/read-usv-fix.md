# A bug in `read_usv()`: matrices with more than one column are transposed

`Utils::read_usv()` parses a whitespace-separated numeric matrix from a text
file. It is used to read PCAone's own `.eigvecs` output back in. For any file
with **more than one column it returns the wrong matrix**, silently.

## The cause

The parser walks the file line by line and pushes each value onto a flat
`std::vector<double>`, so the buffer ends up in **row-major** order. The last
line then hands that buffer to Eigen:

```cpp
return Eigen::Map<Mat2D>(V.data(), j, k);
```

`Mat2D` is `Eigen::MatrixXd`, which is **column-major**. Mapping a row-major
buffer as column-major reinterprets the memory: element `(r, c)` of the result
is taken from buffer position `r + c*j`, where the value actually stored there
is the one for row `r + c*j / k` — so the contents are scrambled, equivalent to
reading the transpose of a `k x j` matrix.

It is correct only when `k == 1`, where the two layouts coincide. That is why
the bug has gone unnoticed: the common case is a single-column file.

A concrete consequence in this branch: `--evaladmix-k 1` reads a 4-column
`.eigvecs` and takes its first column as PC1. Before the fix, the resulting
kinship matrix correlated with the correct answer at only **r = 0.64**, instead
of agreeing exactly. PC1 itself was identical between the two runs
(`cor = 1.000`), so the error was entirely in the reader.

## The fix

Map the buffer with the layout it was written in:

```cpp
using MatRowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
return Eigen::Map<MatRowMajor>(V.data(), j, k);
```

Eigen converts to the column-major `Mat2D` on return, so no caller has to
change.

Two guards were added alongside it. The column count `k` is only assigned
inside a branch that requires at least two lines to have been read, so a
**single-row file returned `k = 0`** and an empty matrix; `k` is now taken from
the row length in that case. And the total element count is checked against
`j * k`, so a ragged file is reported rather than mapped out of bounds.

## Who was affected

- **`src/LD.cpp:493`**, in the ancestry-adjusted LD path. It does

  ```cpp
  Mat2D U = read_usv(params.fileU);
  data->G = (Mat2D::Identity(U.rows(), U.rows()) - U * U.transpose()) * data->G;
  ```

  so with `-k > 1` the projector `U * U'` was built from a scrambled `U`, and the
  residuals — and therefore every downstream LD R², pruning and clumping result —
  were computed against the wrong subspace. With `-k 1` it was correct. This is
  the practically important case, since ancestry-adjusted LD with a single PC is
  unusual.

- `src/EvalAdmix.cpp` (this branch), via `--evaladmix-k`.

## Measured effect on ancestry-adjusted LD

Built both versions and ran the affected path on the 126-sample test data with
**K = 2**:

```bash
PCAone -b smallPlink -k 2 -D -o pcs2                       # 2-column .eigvecs + .mbim
PCAone -b smallPlink -P pcs2 -R --ld-bp 1000000 -o ld      # the path that calls read_usv
```

Both runs emit the same 4,611,677 SNP pairs with identical keys, so only the R²
column moves — but it moves almost everywhere:

| | |
|---|---|
| finite pairs compared | 4,611,432 |
| R² values that differ | 4,609,439 (**99.96%**) |
| correlation between old and new R² | **0.694** |
| mean R² | 0.01240 (old) → 0.01014 (new) |
| mean absolute difference | 0.0072 |
| largest single change | 0.2226 (rs2772239–rs16854621: **0.2556 → 0.0329**) |

The old code **systematically inflates LD**, because the scrambled `U` spans the
wrong subspace and so removes the wrong component of the genotypes. The effect
on pruning is not subtle:

| pairs with R² > 0.2 | count |
|---|---|
| old | 370 |
| new | 58 |
| flagged by old but not new | 322 |
| flagged by new but not old | 10 |

So at an `--ld-r2 0.2` threshold the pre-fix code would have pruned roughly six
times as many pairs as it should.

**Which one is right.** R² was recomputed independently in R for the six
most-affected pairs, replicating what `LD.cpp` does — centre `g/2` by the `.mbim`
allele frequency, residualise on `U U'`, correlate — using the same `.eigvecs`
read correctly:

| SNP A | SNP B | R² old | R² new | independent | matches |
|---|---|---|---|---|---|
| rs2772239 | rs16854621 | 0.255570 | 0.032934 | 0.032934 | **new** |
| rs3011777 | rs10822288 | 0.234191 | 0.020818 | 0.020818 | **new** |
| rs4713994 | rs9369005 | 0.218378 | 0.005801 | 0.005801 | **new** |
| rs6443784 | rs675924 | 0.221030 | 0.009264 | 0.009264 | **new** |
| rs739856 | rs2732439 | 0.242211 | 0.035463 | 0.035463 | **new** |
| rs10991574 | rs4742988 | 0.351003 | 0.144807 | 0.144807 | **new** |

Exact agreement with the fixed version in every case.

Anyone who has run `-D`/`-R`/`--ld-r2`/`--clump` with `-P` and `-k > 1` should
regenerate those results.

## Regression test

`--evaladmix-k` reads a multi-column `.eigvecs` and projects on a subset of it,
so it exercises the reader directly:

```bash
PCAone -b plink -k 1 -d 0 --evaladmix --maf 0.05 -o one
PCAone -b plink -k 4 -d 0 --evaladmix --evaladmix-k 1 --maf 0.05 -o four
# one.kinship and four.kinship must be identical
```

Both runs use the same PC1 and must therefore produce the same matrix. Before
the fix the maximum absolute difference was 8.8e-02; after it, 0.
