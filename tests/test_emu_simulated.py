#!/usr/bin/env python3
"""EMU must recover the simulated population structure under heavy missingness.

The panel comes from tests/simulate_emu_data.py, a scaled-down copy of the
simulation study in Meisner et al. (2021), Bioinformatics 37:1868 -- three
low-Fst populations plus admixed individuals, pseudo-haploid single-read
calls, and the paper's scenario 1 missingness (5-50% per individual, ~29% of
all calls dropped).  Mean imputation alone smears that much missingness into
the leading PCs, so the checks below only pass if the EM step is doing its
job.

What is asserted, in the top-2 PC space of `PCAone --emu`:
  * every unadmixed individual is closest to its own population centroid;
  * each pair of populations is separated by more than twice the pooled
    within-population spread, including the hardest (Han/Dai-like) pair;
  * the admixed individuals land between the source populations, where their
    simulated ancestry proportions say they should.

Standard library only.
"""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulate_emu_data import POP_NAMES, build  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PCAONE = ROOT / "PCAone"

NSNPS = 20000
PER_POP = 84
NADMIXED = 48
SEED = 2021
K = 2

# Calibrated on this fixture over five seeds: `--emu` gives accuracy 1.00,
# separation 4.3-4.8 and correlations 0.95/0.82-0.94, while plain mean-imputed
# PCA on the same panel only reaches separation 1.7-2.0.  The separation bound
# is therefore what distinguishes a working EM step from a broken one.
MIN_ASSIGNMENT_ACCURACY = 0.98
MIN_SEPARATION = 3.0  # centroid distance over pooled within-population spread
MIN_ADMIXTURE_CORRELATION = 0.75


def read_matrix(path: Path) -> list[list[float]]:
    rows = [[float(x) for x in line.split()]
            for line in path.read_text().splitlines() if line.strip()]
    for row in rows:
        for x in row:
            assert math.isfinite(x), f"non-finite value in {path.name}"
    return rows


def standardize(pcs: list[list[float]]) -> list[list[float]]:
    """Put every PC on a common scale so distances weight them equally."""
    n = len(pcs)
    out = [row[:] for row in pcs]
    for j in range(len(pcs[0])):
        mean = sum(row[j] for row in pcs) / n
        sd = math.sqrt(sum((row[j] - mean) ** 2 for row in pcs) / n)
        assert sd > 0, f"PC{j + 1} is constant"
        for row in out:
            row[j] = (row[j] - mean) / sd
    return out


def centroid(points: list[list[float]]) -> list[float]:
    return [sum(p[j] for p in points) / len(points) for j in range(len(points[0]))]


def distance(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def spread(points: list[list[float]], center: list[float]) -> float:
    return math.sqrt(sum(distance(p, center) ** 2 for p in points) / len(points))


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    assert dx > 0 and dy > 0
    return num / (dx * dy)


def run_pcaone(prefix: Path, out: Path) -> Path:
    cmd = [str(PCAONE), "-b", str(prefix), "--emu", "-k", str(K),
           "-n", "2", "-v", "0", "-o", str(out)]
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True,
                            timeout=900)
    assert result.returncode == 0, (cmd, result.stdout, result.stderr)
    return out


def check(prefix: Path, out: Path) -> None:
    labels = [line.split("\t")[1] for line in
              Path(str(prefix) + ".pop").read_text().splitlines() if line.strip()]
    proportions = read_matrix(Path(str(prefix) + ".q"))
    pcs = read_matrix(Path(str(out) + ".eigvecs"))
    assert len(pcs) == len(labels), "one row of eigvecs per individual"
    assert len(pcs[0]) == K, f"expected {K} PCs"
    pcs = standardize(pcs)

    groups = {name: [pcs[i] for i, lab in enumerate(labels) if lab == name]
              for name in POP_NAMES}
    centroids = {name: centroid(points) for name, points in groups.items()}
    spreads = {name: spread(points, centroids[name]) for name, points in groups.items()}

    # 1. every unadmixed individual is nearest to its own population centroid
    correct = sum(1 for name in POP_NAMES for point in groups[name]
                  if min(POP_NAMES, key=lambda o: distance(point, centroids[o])) == name)
    total = sum(len(points) for points in groups.values())
    accuracy = correct / total
    assert accuracy >= MIN_ASSIGNMENT_ACCURACY, f"assignment accuracy {accuracy:.3f}"

    # 2. the populations are separated by more than their internal spread
    worst = min(
        (distance(centroids[a], centroids[b]) / (spreads[a] + spreads[b]), a, b)
        for i, a in enumerate(POP_NAMES) for b in POP_NAMES[i + 1:])
    assert worst[0] >= MIN_SEPARATION, f"{worst[1]}/{worst[2]} separation {worst[0]:.2f}"

    # 3. admixed individuals sit where their true ancestry proportions predict
    admixed = [i for i, lab in enumerate(labels) if lab == "ADMIX"]
    assert admixed, "fixture must contain admixed individuals"
    for j in range(K):
        expected = [sum(proportions[i][k] * centroids[name][j]
                        for k, name in enumerate(POP_NAMES)) for i in admixed]
        observed = [pcs[i][j] for i in admixed]
        r = pearson(expected, observed)
        assert r >= MIN_ADMIXTURE_CORRELATION, f"PC{j + 1} admixture r = {r:.3f}"
        print(f"PC{j + 1} admixture correlation: {r:.3f}")

    print(f"assignment accuracy: {accuracy:.3f}")
    print(f"closest pair {worst[1]}/{worst[2]}: separation {worst[0]:.2f}")


def main() -> None:
    assert PCAONE.is_file(), f"build PCAone first: {PCAONE} not found"
    with tempfile.TemporaryDirectory(prefix="pcaone-emu-sim-") as directory:
        tmp = Path(directory)
        prefix = build(tmp / "sim", nsnps=NSNPS, per_pop=PER_POP,
                       nadmixed=NADMIXED, seed=SEED)
        check(prefix, run_pcaone(prefix, tmp / "emu"))
    print("PASS: emu recovers the simulated structure")


if __name__ == "__main__":
    main()
