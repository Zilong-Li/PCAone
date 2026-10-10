# Relatedness at biobank scale

The dense `--evaladmix` output does not scale: at N = 245,000 (All of Us
srWGS) the Gram matrix alone is 480 GB of RAM, and each N x N text file is
about 540 GB. `--evaladmix-kin <cutoff>` computes the same statistic within
`-m` and writes only the pairs whose kinship reaches the cutoff, plus a maximal
set of unrelated samples:

```bash
PCAone -b cohort -k 16 -m 64 -o pcs
PCAone -b cohort -P pcs -k <K-1> --evaladmix --evaladmix-kin 0.0442 -m 64 -o rel
```

The statistic, its accuracy and the choice of `K-1` PCs are described in
[Relatedness](../small-n/relatedness.md); everything there applies here. The
PCs must come from the same samples, in the same order. `--maf` is not
available out-of-core, so filter rare variants beforehand, e.g. with
`plink2 --maf 0.01 --make-bed`.

## Output

- `rel.kin0`, one line per pair with kinship >= the cutoff, ordered by the
  first sample and then the second, both in file order:
  `#FID1 IID1 FID2 IID2 NSNP KINSHIP` (`#IID1 IID2 NSNP KINSHIP` for a `.psam`
  without FID). `KINSHIP` is the entry `.kinship` would have, to the printed
  digit. `NSNP` is the number of sites where both samples are genotyped. A
  missing call at a site whose frequency is exactly 0, 0.5 or 1 cannot be
  detected (see [Missing genotypes](../small-n/relatedness.md#missing-genotypes-how-they-are-handled)),
  so it counts as genotyped;
- `rel.unrelated`, a maximal set of samples without a pair at or above
  `--evaladmix-unrelated` (default: the cutoff), as `FID IID` lines for
  `--keep`. It uses the greedy rule of Hail's `maximal_independent_set` (which
  All of Us used) and plink2's `--king-cutoff`: drop the sample with the most
  relatives left until none has any, breaking ties by more missing calls, then
  take back every dropped sample whose relatives were all dropped. A sample
  genotyped at none of the sites has no estimate, so it is left out, with a
  warning.

With `--evaladmix-ibd`, `rel.kin0` gains the probabilities of sharing 0, 1 and
2 alleles IBD as three more columns, `... NSNP KINSHIP K0 K1 K2`, and the log
counts how many first-degree pairs have `k0 < 0.125` (parent–offspring) and
how many have more (full sibs). See
[IBD sharing](../small-n/relatedness.md#ibd-sharing---evaladmix-ibd).

## Choosing the cutoff

A cutoff of 0.0442 keeps 3rd degree and closer; All of Us used 0.1 for its
unrelated set. `.kin0` can only report pairs at or above the cutoff it was run
with, and a rerun repeats the whole computation. So choose the cutoff with the
rerun in mind: for example, write the pairs from 0.0221 and build the unrelated
set at 0.1 with `--evaladmix-unrelated 0.1`.

## Checks in the log

The log reports each stripe with an estimate of the time left, and counts the
pairs per KING degree bin. It also compares the scatter of unrelated pairs with
chance. Kinship of unrelated pairs has an sd of about `0.5 / sqrt(Meff)` when
the PCs fit, where `Meff = (sum v)^2 / sum v^2` with `v = f(1-f)` per site
(divided by the mean share of the sites a sample is genotyped at).
The pairs below 0 contain no relatives, so their root mean square measures the
scatter. On the simulations below the two agree to within 0.3% (0.00761
against 0.00761 at N = 20,000). That gives two warnings:

- the scatter is more than 1.3 times chance: the PCs leave structure, as with
  one PC too few (1.47 times), or fine-scale or founder groups — see
  [Choosing the number of PCs](../small-n/relatedness.md#choosing-the-number-of-pcs);
- chance alone should put more than a tenth as many unrelated pairs above the
  cutoff as were written: the cutoff is within the noise of the sites. At N =
  8,000 with 5,000 sites, a cutoff of 0.0221 is 2.9 sd, and chance predicts
  60,285 of the 62,505 pairs written. Raise the cutoff, or use more sites.

## Method details

### How the pairs are computed

In `corres_ij = sb_i sb_j S_ij - (L R')_ij` only `S_ij` is a pair
quantity. Every other term is per sample, including `T = S Q = G(G'Q) - M gbar (gbar'Q)`,
so one pass over the genotypes collects them in `O(NMr)`. The pairs are then
computed in stripes of samples, each holding `S(j, i)` for the stripe's samples
`i` and every `j > i`. Each stripe is one matrix product of the stripe's rows
of `G` with the rows from the stripe down, and one more pass over the genotypes:
from RAM without `-m`, from the file with `-m`. The stripes are sized to fit in
`-m`, so the total work is still the dense path's single `O(N^2 M)` Gram
product, split into pieces. With missing genotypes, the sites a pair has in
common, `n_ij = n - m_i - m_j + mm_ij`, need `mm_ij`, the sites both samples
miss. That count comes from the lists of missing calls, not from a second
`N x N` product. At biobank call rates nearly every site has a missing call,
but only a few samples miss each one. A site missing in a large share of the
stripe falls back to a product of missing indicators. The pair counts are held
in float, so `--evaladmix-kin` refuses more than 2^24 (16.8 million) sites with
missing calls.

With `--evaladmix-ibd`, each stripe holds a second matrix of pairs, for the
dominance residuals, so the same `-m` makes more stripes.

### Memory

`-m` bounds the stripes, the genotype block and the working buffers. Peak RSS
was at most 20 MB above `-m` in every run below. Without `-m` the genotypes are
held in RAM, as in the dense path, and the stripes take up to 2 GB. A small
`-m` only costs extra passes over the genotypes.

### Checks

On every test set, the `.kin0` of a cutoff of -0.5 (all pairs) matches
`.kinship` character for character. That held in-core, with `-m` (1 to 9
stripes), with and without missing genotypes (including sites missing 30% of
calls), and for BED and PGEN input: 32 million pairs at N = 8,000, plus the
smaller sets in `tests/test_evaladmix_pairs.py`.

### Performance

Measured with 2 threads on x86-64. The data are simulated 3-way admixture with
5,000 unlinked sites, 1.4% missing calls (one sample in five at three times the
rate) and planted relatives, and `-k 2`. Time and peak RSS are for the
`--evaladmix` run alone:

| samples | dense `--evaladmix` | `--evaladmix-kin 0.0442` |
|---|---|---|
| 8,000 | 12.4 s, 1.45 GB, plus two 0.6 GB files | 6.4 s at `-m 1` (1.01 GB, 2 stripes); 11.7 s at `-m 0.3` (0.32 GB, 9 stripes) |
| 20,000 | needs 6 GB of RAM for its two N x N matrices, and writes two 3.6 GB files | 39.5 s at `-m 1` (1.02 GB, 5 stripes); 32.4 s at `-m 2` (2.00 GB, 3 stripes) |

At N = 20,000 every planted pair is found: 40/40 duplicates (mean 0.4996),
2,400/2,400 parent–offspring and full-sib pairs (0.2499), and 100/100 half sibs
(0.1244). One of the 2 x 10^8 unrelated pairs reaches 0.0442 (0.0446), which is
expected with 5,000 sites. The `.unrelated` set drops 1,041 of the 2,182 samples
with relatives. On the planted families that is the minimum possible: one per
duplicate, and the children of each family.

With `--evaladmix-ibd`, on N = 3,000 samples and 20,000 sites with 1% missing
and 2 threads, `--evaladmix-kin 0.05 -m 0.2` took 17.8 s instead of 11.4 s, in
the same 0.22 GB, with 4 stripes instead of 3.

### Extrapolation to All of Us (not measured)

The cost is `N^2 M` flops: 9 x 10^15 for N = 245,000 and M = 150,000
LD-pruned sites. The stripes above ran at about 28 GFLOP/s per core, counting
reading, imputation and output. At that rate, 64 cores would take about 1.4
hours, and `-m 64` makes 5 to 10 passes over the 9 GB `.bed`. Linking MKL or
OpenBLAS (see [Install](../installation.md)) usually speeds up the products
further.
