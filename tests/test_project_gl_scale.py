#!/usr/bin/env python3
"""Check that --project 3 honours the scale of the reference PCA.

Uses only the Python standard library. One reference is written three ways:
centred 0..1 (scale=0, gscale=1), centred dosages as pcangsd decomposes them
(scale=0, gscale=2: same V, twice the S), and standardised (scale=-9). With
near-certain GLs the projection must match the oracle x V / S of the true
genotypes, x being the matrix the reference decomposed; with noisy GLs the
output must be a fixed point of the EM that maps the reconstruction back as
pi = f + recon / a_j (PCAngsd's Algorithm 1 is a_j = 2). Writing the BEAGLE
alleles in the other order (GL columns reversed) must not change the result.
"""

import gzip
import math
import random
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PCAONE = ROOT / "PCAone"
N = 30
M = 300
S01 = 5.0


def site_scale(kind, f):
    if kind == "r01":
        return 1.0
    if kind == "r02":
        return 2.0
    return math.sqrt(2.0) / math.sqrt(f * (1.0 - f))


def write_reference(tmp, kind, f, delta):
    a = [site_scale(kind, fj) for fj in f]
    v = [aj * dj for aj, dj in zip(a, delta)]
    norm = math.sqrt(sum(x * x for x in v))
    v = [x / norm for x in v]
    s = 2.0 * S01 if kind == "r02" else S01
    scale, gscale = {"r01": (0, 1), "r02": (0, 2), "rstd": (-9, 1)}[kind]
    prefix = tmp / kind
    Path(str(prefix) + ".loadings").write_text("".join(f"{x:.12g}\n" for x in v))
    Path(str(prefix) + ".sigvals").write_text(
        f"#{N},{M},scale={scale},ploidy=2,gscale={gscale}\n{s:.12g}\n"
    )
    Path(str(prefix) + ".eigvals").write_text(f"{s * s / M:.12g}\n")
    Path(str(prefix) + ".eigvecs").write_text("0\n" * N)
    # PCAone's .mbim counts A1, so A1 = C and F = freq(C): the allele the GLs count
    Path(str(prefix) + ".mbim").write_text(
        "".join(f"1\tsnp{j}\t0\t{j + 1}\tC\tA\t{fj:.12g}\n" for j, fj in enumerate(f))
    )
    return prefix, a, v, s


def write_beagle(path, gls, swapped=False):
    """GLs are for AA, AC, CC (C counted); swapped writes allele1=C, allele2=A"""
    with gzip.open(path, "wt") as out:
        out.write("marker\tallele1\tallele2\t" + "\t".join(f"s{i}\ts{i}\ts{i}" for i in range(N)) + "\n")
        alleles = "1\t0" if swapped else "0\t1"
        for j in range(M):
            row = "\t".join(f"{x:.12g}" for i in range(N) for x in (gls[i][j][::-1] if swapped else gls[i][j]))
            out.write(f"1_{j + 1}\t{alleles}\t{row}\n")


def project(tmp, beagle, prefix, tag):
    out = tmp / f"{prefix.name}.{tag}"
    subprocess.run(
        [str(PCAONE), "-G", str(beagle), "-P", str(prefix), "--project", "3", "-k", "1",
         "--maxiter", "500", "--tol-em", "1e-9", "-n", "1", "-o", str(out)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
    )
    return [float(line.split()[0]) for line in Path(str(out) + ".eigvecs").read_text().split("\n") if line]


def em_step(u, gl, f, a, v, s):
    """one E+M step of the GL projection for K=1, on the reference's scale"""
    num = 0.0
    for j in range(M):
        pt = min(max(f[j] + u * s * v[j] / a[j], 1e-4), 1.0 - 1e-4)
        p0 = gl[j][0] * (1.0 - pt) ** 2
        p1 = gl[j][1] * 2.0 * pt * (1.0 - pt)
        p2 = gl[j][2] * pt * pt
        g = a[j] * ((p1 + 2.0 * p2) / (2.0 * (p0 + p1 + p2)) - f[j])
        num += g * v[j]
    return num / s  # sum(v^2) == 1


def main():
    rng = random.Random(2024)
    base = [rng.uniform(0.15, 0.85) for _ in range(M)]
    delta = [rng.uniform(-0.12, 0.12) for _ in range(M)]
    pops = [[b + d for b, d in zip(base, delta)], [b - d for b, d in zip(base, delta)]]
    f = base  # the reference is an even mix of both populations
    geno = [[int(rng.random() < p) + int(rng.random() < p) for p in pops[i % 2]] for i in range(N)]

    sharp = [[[1.0 - 2e-6 if g == k else 1e-6 for k in range(3)] for g in row] for row in geno]
    noisy = []
    for row in geno:
        out = []
        for g in row:
            depth = rng.randint(0, 3)
            alt = sum(rng.random() < {0: 0.01, 1: 0.5, 2: 0.99}[g] for _ in range(depth))
            lik = [0.01 ** alt * 0.99 ** (depth - alt), 0.5 ** depth, 0.99 ** alt * 0.01 ** (depth - alt)]
            out.append([x / sum(lik) for x in lik])
        noisy.append(out)

    with tempfile.TemporaryDirectory(prefix="pcaone-project-gl-") as directory:
        tmp = Path(directory)
        write_beagle(tmp / "sharp.beagle.gz", sharp)
        write_beagle(tmp / "noisy.beagle.gz", noisy)
        write_beagle(tmp / "swapped.beagle.gz", noisy, swapped=True)
        refs = {kind: write_reference(tmp, kind, f, delta) for kind in ("r01", "r02", "rstd")}

        results = {}
        for kind, (prefix, a, v, s) in refs.items():
            u = project(tmp, tmp / "sharp.beagle.gz", prefix, "sharp")
            oracle = [sum(a[j] * (row[j] / 2.0 - f[j]) * v[j] for j in range(M)) / s for row in geno]
            span = max(abs(x) for x in oracle)
            for x, y in zip(u, oracle):
                assert abs(x - y) <= 0.03 * span, (kind, "sharp GLs vs true genotypes", x, y)
            results[kind] = project(tmp, tmp / "noisy.beagle.gz", prefix, "noisy")
            span = max(abs(x) for x in results[kind])
            for i, x in enumerate(results[kind]):
                y = em_step(x, noisy[i], f, a, v, s)
                assert abs(x - y) <= 1e-3 * span, (kind, "not a fixed point of the EM", i, x, y)
            swapped = project(tmp, tmp / "swapped.beagle.gz", prefix, "swapped")
            for x, y in zip(results[kind], swapped):
                assert abs(x - y) <= 1e-5 * max(1.0, abs(x)), (kind, "BEAGLE allele order changes the result", x, y)
            print("PASS:", kind)

        for x, y in zip(results["r01"], results["r02"]):
            assert abs(x - y) <= 1e-5 * max(1.0, abs(x)), ("0..1 and 0..2 references disagree", x, y)
        print("PASS: 0..1 and 0..2 references agree")


if __name__ == "__main__":
    main()
