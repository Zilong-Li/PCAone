#!/usr/bin/env python3
"""Out-of-core --project 1 and 2 (-m, and the default) must give the in-core coordinates (-m 0).

With -m, the projection streams the target in blocks of sites and keeps sums
per sample (run_projection_blockwise), instead of reading the N x M genotypes.
It is checked against the in-core run on a target that has

- missing calls at rates from 0 to 99.9%, so some samples have fewer called
  sites than PCs (NA), and some miss most of a block,
- reference sites it lacks, sites the reference lacks, and swapped alleles,

for a standardized and an unscaled reference, with blocks of a few dozen sites,
the --project-bootstrap with its replicates in several passes over the
samples, and --no-prefetch. The scores agree to the 6 digits written.

Both also write .eigvecs2, the scores with a #FID IID PC1 .. header and the
FID and IID of the .fam, and the bootstrap files have FID and IID after the
sample index.

Only the standard library is needed: the PLINK files are written here.

    python3 tests/test_project_out_of_core.py
"""

from __future__ import annotations

import math
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PCAONE = str(ROOT / "PCAone")
failures: list[str] = []

NREF, NTGT, M, K = 300, 700, 2000, 4
CODE = {2: 0b00, -1: 0b01, 1: 0b10, 0: 0b11}  # copies of A1 -> BED code


def write_plink(prefix: Path, geno: list[list[int]], bim: list[tuple[str, int, str, str]]) -> None:
    """geno[j][i] = copies of A1 (0/1/2) or -1 for missing; SNP-major"""
    n = len(geno[0])
    with open(f"{prefix}.bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        for row in geno:
            out = bytearray((n + 3) // 4)
            for i, g in enumerate(row):
                out[i // 4] |= CODE[g] << (2 * (i % 4))
            f.write(out)
    with open(f"{prefix}.bim", "w") as f:
        for chrom, bp, a1, a2 in bim:
            f.write(f"{chrom}\t{chrom}_{bp}\t0\t{bp}\t{a1}\t{a2}\n")
    with open(f"{prefix}.fam", "w") as f:
        for i in range(n):
            f.write(f"fam{i // 3} ind{i} 0 0 0 -9\n")


def simulate(tmp: Path) -> None:
    rng = random.Random(2026)
    pops = 3
    sites = [("1" if j < M // 2 else "2", 1000 + 37 * j) for j in range(M)]
    freq = []
    for _ in sites:
        p = rng.uniform(0.05, 0.95)
        freq.append([min(max(rng.gauss(p, 0.12), 0.01), 0.99) for _ in range(pops)])

    def draw(f: float) -> int:
        return (rng.random() < f) + (rng.random() < f)

    # the reference: complete genotypes
    ref = [[draw(freq[j][i % pops]) for i in range(NREF)] for j in range(M)]
    write_plink(tmp / "ref", ref, [(c, bp, "A", "G") for c, bp in sites])

    # the target: admixed samples with missing rates from 0 to 99.9%
    rates = [0.0, 0.01, 0.1, 0.6, 0.97, 0.999]
    miss = [rates[i % len(rates)] for i in range(NTGT)]
    w = [[rng.random() for _ in range(pops)] for _ in range(NTGT)]
    w = [[x / sum(row) for x in row] for row in w]
    keep = [j for j in range(M) if rng.random() > 0.05]  # sites the target lacks
    bim, geno = [], []
    for j in keep:
        row = []
        for i in range(NTGT):
            f = sum(w[i][p] * freq[j][p] for p in range(pops))
            row.append(-1 if rng.random() < miss[i] else draw(f))
        c, bp = sites[j]
        if rng.random() < 0.2:  # the alleles the other way round: count G
            bim.append((c, bp, "G", "A"))
            row = [g if g < 0 else 2 - g for g in row]
        else:
            bim.append((c, bp, "A", "G"))
        geno.append(row)
        if rng.random() < 0.03:  # a site the reference lacks
            bim.append((c, bp + 1, "C", "T"))
            geno.append([draw(0.3) for _ in range(NTGT)])
    write_plink(tmp / "tgt", geno, bim)


def pcaone(args: list[str]) -> str:
    p = subprocess.run([PCAONE, *args, "-n", "2"], capture_output=True, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stdout + p.stderr)
        raise SystemExit(f"PCAone failed: {' '.join(args)}")
    return p.stdout + p.stderr


def read_rows(path: Path, skip_header: bool = False) -> list[list[str]]:
    lines = path.read_text().splitlines()
    return [line.split("\t") for line in (lines[1:] if skip_header else lines)]


def same_numbers(name: str, a: Path, b: Path, header: bool = False, rtol: float = 2e-5) -> None:
    ra, rb = read_rows(a, header), read_rows(b, header)
    if len(ra) != len(rb):
        failures.append(f"{name}: {len(ra)} rows vs {len(rb)}")
        return
    worst = 0.0
    for x, y in zip(ra, rb):
        if len(x) != len(y):
            failures.append(f"{name}: a row has {len(x)} fields vs {len(y)}")
            return
        for u, v in zip(x, y):
            if u == v:
                continue
            try:
                fu, fv = float(u), float(v)
            except ValueError:
                failures.append(f"{name}: {u!r} vs {v!r}")
                return
            if math.isnan(fu) or math.isnan(fv) or math.isinf(fu) or math.isinf(fv):
                failures.append(f"{name}: {u!r} vs {v!r}")
                return
            # %g with 6 digits: a relative error, or an absolute one near 0
            worst = max(worst, abs(fu - fv) / max(abs(fu), abs(fv), 1e-3))
    if worst > rtol:
        failures.append(f"{name}: largest relative difference {worst:.3g}")
    else:
        print(f"ok  {name} ({len(ra)} rows, max rel diff {worst:.2g})")


def check_ids(name: str, out: Path, fam: Path, bootstrap: bool = False) -> None:
    """.eigvecs2: the .eigvecs with the FID and IID of the .fam; the bootstrap files: FID IID after sample"""
    ids = [line.split()[:2] for line in fam.read_text().splitlines()]
    plain = read_rows(Path(f"{out}.eigvecs"))
    lines = Path(f"{out}.eigvecs2").read_text().splitlines()
    want = ["#FID", "IID"] + [f"PC{k + 1}" for k in range(len(plain[0]))]
    rows = [line.split("\t") for line in lines[1:]]
    if lines[0].split("\t") != want:
        failures.append(f"{name}: .eigvecs2 header {lines[0]!r}")
    elif [r[:2] for r in rows] != ids or [r[2:] for r in rows] != plain:
        failures.append(f"{name}: .eigvecs2 is not the .eigvecs with the IDs of the .fam")
    else:
        print(f"ok  {name} .eigvecs2")
    if not bootstrap:
        return
    for ext, first in ((".proj.bootstrap.tsv", 0), (".proj.bootstrap.cov.tsv", 0), (".proj.bootstrap.eigvecs", 1)):
        rows = read_rows(Path(f"{out}{ext}"))
        if rows[0][first : first + 3] != ["sample", "FID", "IID"]:
            failures.append(f"{name} {ext}: header {rows[0]}")
        elif any(r[first + 1 : first + 3] != ids[int(r[first]) - 1] for r in rows[1:]):
            failures.append(f"{name} {ext}: FID IID do not match the sample index")
        else:
            print(f"ok  {name} {ext} IDs")


def main() -> int:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        simulate(tmp)
        for ref, scale in (("std", []), ("raw", ["-C", "0"])):
            pcaone(["-b", str(tmp / "ref"), "-k", str(K), *scale, "-o", str(tmp / ref)])
            for mode in ("1", "2"):
                base = ["-b", str(tmp / "tgt"), "-P", str(tmp / ref), "--project", mode]
                incore, ooc = tmp / f"{ref}.p{mode}.in", tmp / f"{ref}.p{mode}.m"
                pcaone([*base, "-m", "0", "-o", str(incore)])
                log = pcaone([*base, "-m", "0.001", "-o", str(ooc)])
                if "out-of-core projection" not in log:
                    failures.append(f"{ref} --project {mode}: -m did not stream the target")
                same_numbers(f"{ref} --project {mode}", Path(f"{incore}.eigvecs"), Path(f"{ooc}.eigvecs"))
                # without -m, the target is streamed from its 2-bit calls as well
                dflt = tmp / f"{ref}.p{mode}.default"
                log = pcaone([*base, "-o", str(dflt)])
                if "from the 2-bit calls" not in log:
                    failures.append(f"{ref} --project {mode}: without -m, the target was not streamed")
                same_numbers(f"{ref} --project {mode} default", Path(f"{incore}.eigvecs"), Path(f"{dflt}.eigvecs"))
                for tag, out in (("in-core", incore), ("-m", ooc)):
                    check_ids(f"{ref} --project {mode} {tag}", out, tmp / "tgt.fam")
                if mode == "2":
                    na = sum(r[0] == "NA" for r in read_rows(Path(f"{ooc}.eigvecs")))
                    if na == 0:
                        failures.append(f"{ref} --project 2: no sample with fewer called sites than PCs")
                    # --no-prefetch reads the BED with the subset of sites as well
                    nop = tmp / f"{ref}.p2.nop"
                    pcaone([*base, "-m", "0.001", "--no-prefetch", "-o", str(nop)])
                    same_numbers(f"{ref} --project 2 --no-prefetch", Path(f"{incore}.eigvecs"), Path(f"{nop}.eigvecs"))

        # the bootstrap, its samples in groups: one pass over the genotypes each
        base = ["-b", str(tmp / "tgt"), "-P", str(tmp / "std"), "--project", "2", "--project-bootstrap", "20",
                "--project-bootstrap-save", "--seed", "7"]
        incore, ooc = tmp / "boot.in", tmp / "boot.m"
        pcaone([*base, "-m", "0", "-o", str(incore)])
        log = pcaone([*base, "-m", "0.001", "-o", str(ooc)])
        if "passes = 1\n" in log or "passes =" not in log:
            failures.append("--project-bootstrap: -m 0.001 should take the samples in several passes")
        same_numbers("bootstrap .eigvecs", Path(f"{incore}.eigvecs"), Path(f"{ooc}.eigvecs"))
        for tag, out in (("in-core", incore), ("-m", ooc)):
            check_ids(f"bootstrap {tag}", out, tmp / "tgt.fam", bootstrap=True)
        for ext in (".proj.bootstrap.tsv", ".proj.bootstrap.cov.tsv", ".proj.bootstrap.eigvecs"):
            # the variance of nearly flat replicates is a difference of sums
            same_numbers(f"bootstrap {ext}", Path(f"{incore}{ext}"), Path(f"{ooc}{ext}"), header=True, rtol=1e-3)

    if failures:
        print("\n".join(["FAILED:"] + failures))
        return 1
    print("all out-of-core projection checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
