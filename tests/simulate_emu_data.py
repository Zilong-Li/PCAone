#!/usr/bin/env python3
"""Simulate a small EMU-paper-style dataset in PLINK format.

Scaled-down version of the simulation study in Meisner, Liu, Huang and
Albrechtsen (2021), "Large-scale inference of population structure in presence
of missingness using PCA", Bioinformatics 37:1868-1875.

How the paper's design maps onto this script
--------------------------------------------
* Three source populations with low pairwise Fst.  The paper used allele
  frequencies estimated from the HGDP Han, Dai and Uygur samples; here the
  same shape is produced with the Balding-Nichols model on a small tree:
  an East Asian branch splits into two closely related populations (Han- and
  Dai-like), and the third population (Uygur-like) is a 50/50 mixture of that
  branch with a diverged West Eurasian branch.  PC1 then separates the
  Uygur-like population from the two East Asian ones and PC2 splits the
  low-Fst pair, which is the configuration the paper's figures show.
* Admixed individuals on top of the source populations (the paper simulates
  150 admixed individuals out of 900); ancestry proportions are Dirichlet.
* Genotypes are Binomial(2, pi_ij) with pi_ij = sum_k q_ik * f_kj, then turned
  into pseudo-haploid single-read calls exactly as the paper does.
* Sites with minor allele frequency below 5% are discarded, as in the paper.
* Missingness scenario 1 of the paper: each individual draws a missingness
  rate uniformly from 5-50% and every one of its calls is dropped
  independently with that probability.

Single-read shortcut: the paper keeps a homozygote as-is and calls a
heterozygote 0 or 1 with equal probability, so
    P(read = 1) = pi^2 + 2*pi*(1 - pi) * 0.5 = pi .
The pseudo-haploid call is therefore just Bernoulli(pi_ij) and no diploid
genotype has to be drawn first.  Pass --diploid to keep true diploid
genotypes instead (the EM model in PCAone is the same either way).

Outputs <prefix>.bed/.bim/.fam plus three truth files used by the tests:
<prefix>.pop (population label per individual), <prefix>.q (true ancestry
proportions) and <prefix>.miss (target and realised missingness rate).

Standard library only.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

# PLINK .bed codes, keyed by the dosage of the second (alt) allele.
BED_CODE = {0: 3, 1: 2, 2: 0}
BED_MISSING = 1
BED_MAGIC = b"\x6c\x1b\x01"

# Drift parameters of the population tree, in Balding-Nichols Fst units.
# Chosen so the pairwise Fst values land near the HGDP trio used in the paper
# (Han-Dai ~ 0.006, Han-Uygur ~ 0.019, Dai-Uygur ~ 0.024).
DRIFT_EAST_ASIA = 0.020   # global ancestral -> East Asian branch
DRIFT_WEST_EURASIA = 0.080  # global ancestral -> West Eurasian branch
# The tip drift is larger than the paper's closest pair (Fst ~ 0.006 between
# Han and Dai).  This panel carries ~20K sites against the paper's ~350K, so at
# the paper's value the Han/Dai axis sits at the edge of detectability and the
# test would be flaky; 0.010 puts that pair at Fst ~ 0.02, still well inside the
# low-Fst regime.  Pass --drift 0.4 to recover the paper's distances.
DRIFT_TIP = 0.010         # branch -> population tip
UYGUR_WEST_FRACTION = 0.5  # West Eurasian ancestry of the third population

POP_NAMES = ("HAN", "DAI", "UYG")


def balding_nichols(rng: random.Random, p: float, fst: float) -> float:
    """Draw a descendant frequency from ancestral frequency p with drift fst."""
    if fst <= 0:
        return p
    scale = (1.0 - fst) / fst
    alpha = max(p * scale, 1e-6)
    beta = max((1.0 - p) * scale, 1e-6)
    return rng.betavariate(alpha, beta)


def site_frequencies(rng: random.Random, drift: float) -> list[float]:
    """Frequencies of the three source populations at one independent site."""
    ancestral = rng.uniform(0.05, 0.95)
    east = balding_nichols(rng, ancestral, DRIFT_EAST_ASIA * drift)
    west = balding_nichols(rng, ancestral, DRIFT_WEST_EURASIA * drift)
    han = balding_nichols(rng, east, DRIFT_TIP * drift)
    dai = balding_nichols(rng, east, DRIFT_TIP * drift)
    uygur = ((1.0 - UYGUR_WEST_FRACTION) * balding_nichols(rng, east, DRIFT_TIP * drift)
             + UYGUR_WEST_FRACTION * balding_nichols(rng, west, DRIFT_TIP * drift))
    return [han, dai, uygur]


def ancestry_proportions(rng: random.Random, per_pop: int, nadmixed: int,
                         alpha: float) -> tuple[list[list[float]], list[str]]:
    """Per-individual ancestry proportions and population labels."""
    proportions: list[list[float]] = []
    labels: list[str] = []
    for k, name in enumerate(POP_NAMES):
        for _ in range(per_pop):
            q = [0.0, 0.0, 0.0]
            q[k] = 1.0
            proportions.append(q)
            labels.append(name)
    for _ in range(nadmixed):
        gammas = [rng.gammavariate(alpha, 1.0) for _ in POP_NAMES]
        total = sum(gammas)
        proportions.append([g / total for g in gammas])
        labels.append("ADMIX")
    return proportions, labels


def simulate(nsnps: int, per_pop: int, nadmixed: int, seed: int, maf: float,
             alpha: float, drift: float, diploid: bool,
             miss_low: float, miss_high: float):
    """Return (rows, labels, proportions, rates, realised) for the whole panel.

    Each element of `rows` is the list of PLINK codes for one kept site.
    """
    rng = random.Random(seed)
    proportions, labels = ancestry_proportions(rng, per_pop, nadmixed, alpha)
    nsamples = len(labels)
    rates = [rng.uniform(miss_low, miss_high) for _ in range(nsamples)]

    rows: list[list[int]] = []
    missing = [0] * nsamples
    attempts = 0
    while len(rows) < nsnps:
        attempts += 1
        freqs = site_frequencies(rng, drift)
        dosages = []
        for q in proportions:
            pi = q[0] * freqs[0] + q[1] * freqs[1] + q[2] * freqs[2]
            if diploid:
                dosages.append((rng.random() < pi) + (rng.random() < pi))
            else:
                # single-read sampling collapses to one Bernoulli draw
                dosages.append(2 * (rng.random() < pi))
        # the paper filters on the complete data, before missingness is added
        freq = sum(dosages) / (2.0 * nsamples)
        if min(freq, 1.0 - freq) < maf:
            continue
        codes = []
        for i, dosage in enumerate(dosages):
            if rng.random() < rates[i]:
                codes.append(BED_MISSING)
                missing[i] += 1
            else:
                codes.append(BED_CODE[dosage])
        rows.append(codes)

    realised = [m / float(nsnps) for m in missing]
    if attempts > 50 * nsnps:  # pragma: no cover - guards a pathological --maf
        raise SystemExit("MAF filter rejected too many sites; lower --maf")
    return rows, labels, proportions, rates, realised


def write_plink(prefix: Path, rows, labels) -> None:
    nsamples = len(labels)
    padded = (4 - nsamples % 4) % 4
    bed = bytearray(BED_MAGIC)
    bim = []
    for j, codes in enumerate(rows):
        packed = codes + [0] * padded
        for i in range(0, len(packed), 4):
            bed.append(sum(packed[i + b] << (2 * b) for b in range(4)))
        bim.append(f"1\tsnp{j + 1}\t0\t{j + 1}\tA\tC\n")
    # not Path.with_suffix: a prefix may legitimately contain a dot
    Path(str(prefix) + ".bed").write_bytes(bytes(bed))
    Path(str(prefix) + ".bim").write_text("".join(bim))
    Path(str(prefix) + ".fam").write_text(
        "".join(f"{labels[i]}_{i + 1} {labels[i]}_{i + 1} 0 0 0 -9\n"
                for i in range(nsamples)))


def write_truth(prefix: Path, labels, proportions, rates, realised) -> None:
    Path(str(prefix) + ".pop").write_text(
        "".join(f"{labels[i]}_{i + 1}\t{labels[i]}\n" for i in range(len(labels))))
    Path(str(prefix) + ".q").write_text(
        "".join("\t".join(f"{x:.6f}" for x in q) + "\n" for q in proportions))
    Path(str(prefix) + ".miss").write_text(
        "#target\trealised\n"
        + "".join(f"{rates[i]:.6f}\t{realised[i]:.6f}\n" for i in range(len(rates))))


def build(prefix: Path, nsnps: int = 20000, per_pop: int = 84, nadmixed: int = 48,
          seed: int = 2021, maf: float = 0.05, alpha: float = 1.0,
          drift: float = 1.0, diploid: bool = False,
          miss_low: float = 0.05, miss_high: float = 0.50) -> Path:
    """Simulate and write the panel; returns the PLINK prefix."""
    rows, labels, proportions, rates, realised = simulate(
        nsnps, per_pop, nadmixed, seed, maf, alpha, drift, diploid,
        miss_low, miss_high)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_plink(prefix, rows, labels)
    write_truth(prefix, labels, proportions, rates, realised)
    return prefix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("prefix", type=Path, help="output PLINK prefix")
    parser.add_argument("--nsnps", type=int, default=20000,
                        help="sites kept after the MAF filter (default 20000)")
    parser.add_argument("--per-pop", type=int, default=84,
                        help="unadmixed individuals per population (default 84)")
    parser.add_argument("--nadmixed", type=int, default=48,
                        help="admixed individuals (default 48)")
    parser.add_argument("--seed", type=int, default=2021)
    parser.add_argument("--maf", type=float, default=0.05,
                        help="minor allele frequency threshold (default 0.05)")
    parser.add_argument("--alpha", type=float, default=1.0,
                        help="Dirichlet concentration for admixed individuals")
    parser.add_argument("--drift", type=float, default=1.0,
                        help="multiplier on all Fst values; >1 makes the "
                             "populations easier to separate")
    parser.add_argument("--diploid", action="store_true",
                        help="keep diploid genotypes instead of the paper's "
                             "pseudo-haploid single-read calls")
    parser.add_argument("--miss-low", type=float, default=0.05)
    parser.add_argument("--miss-high", type=float, default=0.50,
                        help="per-individual missingness is uniform on "
                             "[--miss-low, --miss-high]; the paper's scenario 1 "
                             "is 0.05-0.50 and its scenario 2 is 0.90-0.99")
    args = parser.parse_args()

    build(args.prefix, nsnps=args.nsnps, per_pop=args.per_pop,
          nadmixed=args.nadmixed, seed=args.seed, maf=args.maf,
          alpha=args.alpha, drift=args.drift, diploid=args.diploid,
          miss_low=args.miss_low, miss_high=args.miss_high)
    nsamples = 3 * args.per_pop + args.nadmixed
    print(f"wrote {args.prefix}.bed/.bim/.fam: {nsamples} individuals x "
          f"{args.nsnps} sites")


if __name__ == "__main__":
    main()
