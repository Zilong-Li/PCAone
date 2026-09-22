#!/usr/bin/env python3
"""Check that in-core EM results are invariant to winSVD's SNP permutation.

Uses only the Python standard library and small generated BED/BEAGLE inputs.
The full sample-space sketch avoids confusing approximation error with a
change in the EM model. Compare reconstructed matrices to allow SVD sign flips.
"""

import gzip
import math
import random
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PCAONE = ROOT / "PCAone"
N = 12
M = 48


def matrix(path):
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()
            if line and not line.startswith("#")]


def close(left, right, label, tolerance=2e-4):
    assert len(left) == len(right), label
    assert all(len(a) == len(b) for a, b in zip(left, right)), label
    for a, b in zip(left, right):
        for x, y in zip(a, b):
            assert math.isfinite(x) and math.isfinite(y), label
            assert abs(x - y) <= tolerance * max(1.0, abs(x), abs(y)), (label, x, y)


def reconstructed(prefix):
    u = matrix(Path(str(prefix) + ".eigvecs"))
    s = [row[0] for row in matrix(Path(str(prefix) + ".sigvals"))]
    v = matrix(Path(str(prefix) + ".loadings"))
    return [[sum(ui[k] * s[k] * vj[k] for k in range(len(s))) for vj in v] for ui in u]


def fixture(tmp):
    rng = random.Random(491)
    bed = bytearray(b"\x6c\x1b\x01")
    bim = []
    gl = ["marker allele1 allele2 " + " ".join(f"s{i} s{i} s{i}" for i in range(N))]
    for j in range(M):
        calls = []
        likelihoods = []
        for i in range(N):
            p = min(0.94, 0.03 + (j % 8) * 0.085 + (0.19 if i < N // 2 else 0))
            g = int(rng.random() < p) + int(rng.random() < p)
            missing = i == (j * 5) % N
            calls.append(1 if missing else {0: 3, 1: 2, 2: 0}[g])
            probs = [1 / 3] * 3 if missing else [0.03, 0.03, 0.03]
            if not missing:
                probs[g] = 0.94
            likelihoods.extend(f"{x:.12g}" for x in probs)
        for i in range(0, N, 4):
            bed.append(sum(calls[i + k] << (2 * k) for k in range(4)))
        bim.append(f"1\tsnp{j}\t0\t{j + 1}\tA\tC\n")
        gl.append(f"1_{j + 1} C A " + " ".join(likelihoods))
    prefix = tmp / "input"
    prefix.with_suffix(".bed").write_bytes(bed)
    prefix.with_suffix(".bim").write_text("".join(bim))
    prefix.with_suffix(".fam").write_text("".join(f"s{i} s{i} 0 0 0 -9\n" for i in range(N)))
    beagle = tmp / "input.beagle.gz"
    with gzip.open(beagle, "wt") as out:
        out.write("\n".join(gl) + "\n")
    Path(str(beagle) + ".bim").write_text("".join(bim))
    return prefix, beagle


def run(tmp, name, input_args, extra, solver):
    prefix = tmp / f"{name}_{solver}"
    cmd = [str(PCAONE), *input_args, "-k", "2", "--oversamples", "10", "-w", "4",
           "--maxp", "20", "--tol-rsvd", "1e-12", "--tol-em", "0", "--maxiter", "4",
           "-V", "-n", "1", "-v", "0", "-o", str(prefix), *extra]
    cmd += ["-d", "1"] if solver == "single" else ["-d", "2"]
    if solver == "ordered":
        cmd.append("-S")
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, (cmd, result.stdout, result.stderr)
    if solver == "shuffled":
        # The same permutation must survive every EM iteration.
        assert Path(str(prefix) + ".log").read_text().count("permuting data matrix") == 1
    return prefix


def main():
    with tempfile.TemporaryDirectory(prefix="pcaone-em-order-") as directory:
        tmp = Path(directory)
        bed, beagle = fixture(tmp)
        cases = [
            ("pca", ["-b", str(bed)], []),
            ("emu_final_scale", ["-b", str(bed)], ["--emu", "--maxiter", "0"]),
            ("emu", ["-b", str(bed)], ["--emu"]),
            ("emu_filtered", ["-b", str(bed)], ["--emu", "--maf", "0.15"]),
            ("emu_centered", ["-b", str(bed)], ["--emu", "-C", "0"]),
            ("pcangsd", ["-G", str(beagle)], []),
            ("pcangsd_filtered", ["-G", str(beagle)], ["--maf", "0.15"]),
        ]
        for name, inputs, extra in cases:
            outputs = [run(tmp, name, inputs, extra, solver)
                       for solver in ("single", "ordered", "shuffled")]
            reference = outputs[0]
            ref_mbim = Path(str(reference) + ".mbim").read_text()
            if "filtered" in name:
                assert 16 <= len(ref_mbim.splitlines()) < M, "fixture must exercise filtering"
            for output in outputs[1:]:
                close(matrix(Path(str(reference) + ".sigvals")),
                      matrix(Path(str(output) + ".sigvals")), name + " singular values")
                close(reconstructed(reference), reconstructed(output), name + " USV")
                assert ref_mbim == Path(str(output) + ".mbim").read_text(), name + " SNP frequencies"
                if name.startswith("pcangsd"):
                    close(matrix(Path(str(reference) + ".cov")),
                          matrix(Path(str(output) + ".cov")), name + " covariance")
            print("PASS:", name)


if __name__ == "__main__":
    main()
