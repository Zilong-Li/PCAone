#!/usr/bin/env python3
"""--em-k sets the PCs that model the individual allele frequencies in the
EM-PCA iterations (--emu, --pcangsd); -k sets the PCs written from the final
matrix.

The panel is that of tests/simulate_emu_data.py (three populations plus
admixed individuals, ~29% missing calls), as PLINK calls for --emu and
--pcangsd and as genotype likelihoods (BEAGLE) for PCAngsd. With every solver,
in-core and with -m where the input allows it:

  * no --em-k, and --em-k equal to -k, give the same bytes as before the
    option existed, i.e. the EM models with -k PCs;
  * --em-k e -k K writes K PCs (.eigvecs, .eigvals, .sigvals, and .eigvecs2
    from PLINK input) whatever e is, and its leading PCs are those of
    --em-k e -k e: the EM with e PCs is the same run, only the decomposition
    of the final matrix keeps more or fewer PCs. That holds both for K > e
    (model with few PCs, write many) and K < e. For PCAngsd it holds to the
    EM's tolerance, as -k e writes the EM's last decomposition itself;
  * BEAGLE input: the .cov depends on the model only, so it is the same bytes
    for any -k, and its K eigenvectors in .eigvecs2 lead with those of -k e;
  * the model is what --em-k says: -k 5 --em-k 2 and -k 5 differ.

An exact numpy EMU and PCAngsd (EM with e PCs, then the K PCs of the final
matrix and, for BEAGLE, the covariance) were checked by hand against every
path at a matched --tol-em of 1e-9: IRAM agreed to 0.05 deg, the RSVDs to
0.45 deg, and the .cov to 1e-5.

Standard library only.

    python3 tests/test_em_k.py
"""

from __future__ import annotations

import gzip
import math
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulate_emu_data import build, simulate  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PCAONE = ROOT / "PCAone"

OUT_OF_CORE = ["-m", "0.002"]  # several blocks of the 160 x 6000 panel
SOLVERS = (
    ("iram", ["-d", "0"]),
    ("iram-ooc", ["-d", "0", *OUT_OF_CORE]),
    ("winsvd", ["-d", "2"]),
    ("winsvd-ooc", ["-d", "2", *OUT_OF_CORE]),
    ("ssvd", ["-d", "1"]),
    ("ssvd-ooc", ["-d", "1", *OUT_OF_CORE]),
)
# leading PCs of two final decompositions of the same EM fit. IRAM solves each
# to --itol ("tight"); the RSVDs stop at --tol-rsvd, which leaves them ~0.3 deg
# apart ("loose"). PCAngsd's -k e run writes the EM's own last decomposition,
# one iteration before the fit that the -k K run decomposes, so the two agree
# to the EM's tolerance: "loose" for every solver, with --tol-em tightened.
MIN_COR = {"tight": 0.999999, "loose": 0.9999}
MAX_SIGVAL_REL = {"tight": 1e-5, "loose": 1e-4}
PCANGSD_TOL = ["--tol-em", "1e-8"]

failures: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(("ok   " if ok else "FAIL ") + name + (f": {detail}" if detail else ""))
    if not ok:
        failures.append(name)


def run(out: Path, args: list[str]) -> Path:
    cmd = [str(PCAONE), "-n", "2", "-v", "0", "-o", str(out), *args]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, (cmd, r.stdout, r.stderr)
    return out


def read(path: Path, skip: int = 0) -> list[list[float]]:
    return [[float(x) for x in line.split()[skip:]] for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def columns(rows: list[list[float]]) -> list[list[float]]:
    return [list(c) for c in zip(*rows)]


def abs_cor(x: list[float], y: list[float]) -> float:
    """|cos| of two eigenvectors"""
    num = sum(a * b for a, b in zip(x, y))
    return abs(num) / math.sqrt(sum(a * a for a in x) * sum(b * b for b in y))


def same_bytes(a: Path, b: Path, suffixes: tuple[str, ...]) -> bool:
    return all(Path(f"{a}{s}").read_bytes() == Path(f"{b}{s}").read_bytes() for s in suffixes)


def leading_agree(label: str, kind: str, a: Path, b: Path, n: int) -> None:
    """the leading n PCs of runs a and b are the same"""
    ua, ub = columns(read(Path(f"{a}.eigvecs"))), columns(read(Path(f"{b}.eigvecs")))
    sa, sb = [r[0] for r in read(Path(f"{a}.sigvals"))], [r[0] for r in read(Path(f"{b}.sigvals"))]
    cor = min(abs_cor(ua[j], ub[j]) for j in range(n))
    rel = max(abs(sa[j] - sb[j]) / sb[j] for j in range(n))
    check(f"{label}: leading {n} PCs", cor >= MIN_COR[kind] and rel <= MAX_SIGVAL_REL[kind],
          f"min |cor| {cor:.7f}, max sigval rel. diff {rel:.1e}")


def leading_cov_vecs(label: str, a: Path, b: Path, n: int) -> None:
    """BEAGLE: the same .cov, and .eigvecs2 of a leads with that of b"""
    va = [line.split()[: 2 + n] for line in Path(f"{a}.eigvecs2").read_text().splitlines()]
    vb = [line.split()[: 2 + n] for line in Path(f"{b}.eigvecs2").read_text().splitlines()]
    check(f"{label}: same .cov, leading {n} .eigvecs2", same_bytes(a, b, (".cov",)) and va == vb)


def shape(label: str, out: Path, nsamples: int, k: int, eigvecs2: bool) -> None:
    vecs = read(Path(f"{out}.eigvecs"))
    ok = (len(vecs) == nsamples and all(len(r) == k for r in vecs) and len(read(Path(f"{out}.eigvals"))) == k
          and len(read(Path(f"{out}.sigvals"))) == k)
    if eigvecs2:
        header = Path(f"{out}.eigvecs2").read_text().splitlines()[0].split()
        ok = ok and header[2:] == [f"PC{j + 1}" for j in range(k)]
    check(f"{label}: writes {k} PCs", ok)


def write_beagle(path: Path, nsnps: int, per_pop: int, nadmixed: int, seed: int) -> int:
    """the diploid panel as genotype likelihoods: 0.94 on the call, flat if missing"""
    rows, labels, *_ = simulate(nsnps, per_pop, nadmixed, seed, 0.05, 1.0, 1.0, True, 0.05, 0.50)
    copies = {3: 0, 2: 1, 0: 2}  # PLINK code -> copies of allele 2
    lines = ["marker allele1 allele2 " + " ".join(f"s{i} s{i} s{i}" for i in range(len(labels)))]
    for j, codes in enumerate(rows):
        gl = []
        for c in codes:
            if c == 1:
                gl += ["0.333333"] * 3
            else:
                p = ["0.03"] * 3
                p[copies[c]] = "0.94"
                gl += p
        lines.append(f"1_{j + 1} A C " + " ".join(gl))
    with gzip.open(path, "wt") as out:
        out.write("\n".join(lines) + "\n")
    return len(labels)


def run_cases(label: str, kind: str, args: list[str], nsamples: int, tmp: Path, beagle: bool) -> None:
    o = lambda name: tmp / f"{label.replace(' ', '-')}-{name}"  # noqa: E731
    sufs = (".eigvecs", ".eigvals", ".sigvals") + ((".cov", ".eigvecs2") if beagle else ())
    k2 = run(o("k2"), ["-k", "2", *args])
    k2e2 = run(o("k2e2"), ["-k", "2", "--em-k", "2", *args])
    check(f"{label}: --em-k equal to -k changes nothing", same_bytes(k2, k2e2, sufs))

    # model with 2 PCs, write 5
    k5e2 = run(o("k5e2"), ["-k", "5", "--em-k", "2", *args])
    shape(f"{label} -k 5 --em-k 2", k5e2, nsamples, 5, not beagle)
    leading_agree(f"{label} -k 5 --em-k 2 vs -k 2", kind, k5e2, k2, 2)
    if beagle:
        leading_cov_vecs(f"{label} -k 5 --em-k 2 vs -k 2", k5e2, k2, 2)

    # model with 4 PCs, write 2
    k4 = run(o("k4"), ["-k", "4", *args])
    k2e4 = run(o("k2e4"), ["-k", "2", "--em-k", "4", *args])
    shape(f"{label} -k 2 --em-k 4", k2e4, nsamples, 2, not beagle)
    leading_agree(f"{label} -k 2 --em-k 4 vs -k 4", kind, k2e4, k4, 2)
    if beagle:
        leading_cov_vecs(f"{label} -k 2 --em-k 4 vs -k 4", k2e4, k4, 2)

    # and the model is what --em-k says: 5 PCs model other frequencies than 2
    k5 = run(o("k5"), ["-k", "5", *args])
    s5 = [r[0] for r in read(Path(f"{k5}.sigvals"))]
    s5e2 = [r[0] for r in read(Path(f"{k5e2}.sigvals"))]
    rel = max(abs(a - b) / b for a, b in zip(s5[:2], s5e2[:2]))
    check(f"{label}: -k 5 --em-k 2 is not -k 5", rel > 10 * MAX_SIGVAL_REL[kind], f"leading sigvals differ by {rel:.1e}")


def main() -> int:
    assert PCAONE.is_file(), f"build PCAone first: {PCAONE} not found"
    with tempfile.TemporaryDirectory(prefix="pcaone-em-k-") as directory:
        tmp = Path(directory)
        prefix = build(tmp / "sim", nsnps=6000, per_pop=40, nadmixed=40, seed=2021)
        nsamples = len(Path(f"{prefix}.fam").read_text().splitlines())
        beagle = tmp / "sim.beagle.gz"
        nbeagle = write_beagle(beagle, nsnps=3000, per_pop=40, nadmixed=40, seed=2022)
        for solver, extra in SOLVERS:
            kind = "tight" if solver.startswith("iram") else "loose"
            run_cases(f"emu {solver}", kind, ["-b", str(prefix), "--emu", *extra], nsamples, tmp, False)
            run_cases(f"pcangsd plink {solver}", "loose", ["-b", str(prefix), "--pcangsd", *PCANGSD_TOL, *extra],
                      nsamples, tmp, False)
            if not solver.endswith("-ooc"):  # BEAGLE runs in-core only
                run_cases(f"pcangsd beagle {solver}", "loose", ["-G", str(beagle), *PCANGSD_TOL, *extra], nbeagle,
                          tmp, True)
    if failures:
        print(f"FAILED: {len(failures)} check(s)")
        return 1
    print("SUCCESS: --em-k models with its PCs and -k sets the PCs written.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
