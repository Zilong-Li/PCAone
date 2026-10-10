# Small-N analysis

These guides are for studies of up to a few thousand samples: a reference panel
such as 1000 Genomes, a population-genetics sample of a non-model species, or
low-depth sequencing analysed from genotype likelihoods. At this size an
`N x N` matrix fits in memory, so PCAone can compute the exact PCA and the full
matrices of pairwise statistics.

Every analysis after the PCA reads the PCs back with `-P/--USV`, so one PCA run
serves them all. Compute it once with a generous `-k`, then select the leading
PCs each analysis needs with `-k`:

```shell
./PCAone -b example/plink -d 3 -k 10 -o ref               # the PCs
./PCAone -b example/plink -P ref -k 3 --ld-r2 0.2 -o adj  # LD pruning on 3 of them
```

| Guide | What it does |
|---|---|
| [Exact PCA and EM-PCA](pca.md) | `--svd 3`; `--emu` for missing calls; `--pcangsd` for genotype likelihoods |
| [Projection](projection.md) | place new samples on the PCs of a reference (`--project`) |
| [Selection scan](selection.md) | per-variant statistics for differentiation along the PCs (`--selection`) |
| [HWE and inbreeding](hwe.md) | per-site HWE tests and per-sample F under population structure (`--inbreed`) |
| [Relatedness](relatedness.md) | kinship and IBD sharing from the residuals of the PCs (`--evaladmix`) |
| [Ancestry-adjusted LD](ld.md) | LD on the residuals of the PCs: R2, pruning and clumping |

Each guide ends with **Method details**: the model, how PCAone computes it, and
how it was checked. The commands use the [example data](../guide/input-output.md#example-data).
For tens of thousands of samples or more, see [Biobank analysis](../biobank/index.md).
