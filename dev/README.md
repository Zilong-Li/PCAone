# Developer notes

Architecture and bug-fix write-ups for people working on the PCAone source.
They are kept in the repository but are not part of the documentation website,
whose sources are in [`docs/`](../docs).

- [Architecture](architecture.md): class diagram of the readers and PCA methods.
- [HWE test fixes](hwe-lrt-fix.md): two bugs in `--inbreed 1`, fixed in v0.8.0.
- [LD robustness fixes](ld-robustness-fixes.md): NaN R2 and a segfault on a missing `.mbim`.
- [`read_usv()` fix](read-usv-fix.md): multi-column `.eigvecs` were transposed.
- [v0.8.0 change notes](changes-v0.8.0.md): the full notes behind the v0.8.0 change log entry.
