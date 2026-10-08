# Projection

Project new samples onto existing PCs is supported with `--project` option.
First, we run PCAone on a set of reference samples, which writes the loadings
and `.mbim` by default:

```shell
PCAone -b example/ref -k 10 -o ref
```

Pass the reference prefix with `-P/--USV` to read its loadings, singular
values and allele frequencies. PCAone matches overlapping markers against
`.mbim` and corrects flipped alleles. Both datasets must use the same genome
build and compatible alleles.

Below is a command using example PLINK data to project new target samples onto
the reference coordinates.

```shell
## --USV: prefix to .eigvecs, .sigvals, .loadings, .mbim
## --project: check the manual on projection methods
PCAone -b example/new \
       --USV ref \
       --project 2 \
       -o new
```

The projection uses every PC in the reference; add `-k` to use only the leading
ones (e.g. `-k 3`). This holds for every analysis that reads a reference with
`-P/--USV`: `--project`, `--selection`, `--inbreed`, `--evaladmix` and the
[LD analyses](ld.md), so one reference run with a generous `-k` serves them all.

To visualize uncertainty in projected coordinates, use `--project-bootstrap` to
resample matched SNPs with replacement after the baseline `--project 2`
projection. PCAone writes per-sample/per-PC summaries to `.proj.bootstrap.tsv`
and PC-pair covariance matrices to `.proj.bootstrap.cov.tsv`, which can be used
to draw uncertainty ellipses around the baseline estimates. Add
`--project-bootstrap-save` to write raw bootstrap replicate coordinates to
`.proj.bootstrap.eigvecs`.

```shell
PCAone -b example/new \
       --USV ref \
       --project 2 \
       --project-bootstrap 100 \
       -o new
```

For BEAGLE genotype likelihood input, use `--project 3` to perform
genotype-likelihood aware projection with an EM procedure:

```shell
PCAone -G example/target.beagle.gz \
       --USV ref \
       --project 3 \
       -o target
```

The target is scaled as the reference PCA was (recorded in its `.sigvals`), so
`--scale` is not needed here.

**NB:** The two BEAGLE alleles must be the two alleles of the reference `.mbim`,
in either order. BEAGLE likelihoods count `allele2`, which PCAone matches to the
reference's counted allele `A1` (5th column of `.mbim`), and sites where they
are swapped are flipped automatically; sites with other alleles are dropped.
Use angsd `-doMajorMinor 3` with a `-sites` file taken from the reference
`.mbim` to fix the alleles.
