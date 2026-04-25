#!/usr/bin/env python3

from __future__ import annotations

import gzip
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PCAONE = ROOT / "PCAone"


def run(cmd: list[str], cwd: Path = ROOT, timeout: int = 180) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise AssertionError(f"command failed with exit code {result.returncode}: {' '.join(cmd)}")
    return result


def maybe_make_pgen_from_plink(tmp: Path) -> Path | None:
    plink2 = shutil.which("plink2")
    if plink2 is None:
        print("SKIP: plink2 is not available; cannot synthesize matched PGEN input")
        return None

    out = tmp / "plink_as_pgen"
    cmd = [
        plink2,
        "--bfile",
        str(ROOT / "example/plink"),
        "--make-pgen",
        "--out",
        str(out),
    ]
    try:
        result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=120)
    except subprocess.TimeoutExpired:
        print("SKIP: plink2 timed out while converting PLINK to PGEN")
        return None
    if result.returncode != 0:
        print("SKIP: plink2 failed while converting PLINK to PGEN")
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        return None
    for suffix in (".pgen", ".pvar", ".psam"):
        if not Path(f"{out}{suffix}").exists():
            print(f"SKIP: plink2 did not create {out}{suffix}")
            return None
    return out


def read_matrix(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                rows.append([float(x) for x in line.split()])
    return rows


def read_vector(path: Path) -> list[float]:
    return [row[0] for row in read_matrix(path)]


def assert_close_scalar(left: float, right: float, tol: float, label: str) -> None:
    if not math.isfinite(left) or not math.isfinite(right) or abs(left - right) > tol:
        raise AssertionError(f"{label}: {left} != {right} within {tol}")


def assert_close_vector(left: Path, right: Path, tol: float) -> None:
    xs = read_vector(left)
    ys = read_vector(right)
    if len(xs) != len(ys):
        raise AssertionError(f"length mismatch: {left}={len(xs)}, {right}={len(ys)}")
    for i, (x, y) in enumerate(zip(xs, ys), start=1):
        assert_close_scalar(x, y, tol, f"{left.name}:{i}")


def assert_close_matrix(left: Path, right: Path, tol: float) -> None:
    xs = read_matrix(left)
    ys = read_matrix(right)
    if len(xs) != len(ys):
        raise AssertionError(f"row count mismatch: {left}={len(xs)}, {right}={len(ys)}")
    for i, (xr, yr) in enumerate(zip(xs, ys), start=1):
        if len(xr) != len(yr):
            raise AssertionError(f"column count mismatch at row {i}: {left}={len(xr)}, {right}={len(yr)}")
        for j, (x, y) in enumerate(zip(xr, yr), start=1):
            assert_close_scalar(x, y, tol, f"{left.name}:{i}:{j}")


def column_signs(reference: Path, observed: Path) -> list[float]:
    ref = read_matrix(reference)
    obs = read_matrix(observed)
    if not ref or len(ref) != len(obs) or len(ref[0]) != len(obs[0]):
        raise AssertionError(f"cannot align signs for {reference} and {observed}")
    signs: list[float] = []
    for col in range(len(ref[0])):
        dot = sum(ref[row][col] * obs[row][col] for row in range(len(ref)))
        signs.append(-1.0 if dot < 0.0 else 1.0)
    return signs


def assert_close_matrix_with_signs(left: Path, right: Path, signs: list[float], tol: float) -> None:
    xs = read_matrix(left)
    ys = read_matrix(right)
    if len(xs) != len(ys):
        raise AssertionError(f"row count mismatch: {left}={len(xs)}, {right}={len(ys)}")
    for i, (xr, yr) in enumerate(zip(xs, ys), start=1):
        if len(xr) != len(yr) or len(xr) != len(signs):
            raise AssertionError(f"column count mismatch at row {i}: {left}, {right}")
        for j, (x, y) in enumerate(zip(xr, yr), start=1):
            assert_close_scalar(x, signs[j - 1] * y, tol, f"{left.name}:{i}:{j}")


def assert_mbim_equal(left: Path, right: Path, tol: float = 1e-8) -> None:
    with left.open() as lf, right.open() as rf:
        for line_no, (lline, rline) in enumerate(zip(lf, rf), start=1):
            lfields = lline.strip().split()
            rfields = rline.strip().split()
            if lfields[:6] != rfields[:6]:
                raise AssertionError(f"mbim metadata mismatch at line {line_no}: {lfields[:6]} != {rfields[:6]}")
            assert_close_scalar(float(lfields[6]), float(rfields[6]), tol, f"mbim AF line {line_no}")


def assert_hwe_equal(left: Path, right: Path, tol: float = 1e-6) -> None:
    with left.open() as lf, right.open() as rf:
        lheader = next(lf)
        rheader = next(rf)
        if lheader != rheader:
            raise AssertionError("HWE header mismatch")
        for line_no, (lline, rline) in enumerate(zip(lf, rf), start=2):
            lfields = lline.strip().split()
            rfields = rline.strip().split()
            if lfields[0] != rfields[0]:
                raise AssertionError(f"HWE ID mismatch at line {line_no}: {lfields[0]} != {rfields[0]}")
            for col in range(1, 4):
                assert_close_scalar(float(lfields[col]), float(rfields[col]), tol, f"HWE line {line_no} col {col}")


def assert_gzip_text_equal(left: Path, right: Path) -> None:
    with gzip.open(left, "rt") as lf, gzip.open(right, "rt") as rf:
        ltext = lf.read()
        rtext = rf.read()
    if ltext != rtext:
        raise AssertionError(f"gzip text mismatch: {left} != {right}")


def out(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def run_pca_pair(tmp: Path, pgen_prefix: Path) -> Path:
    plink = tmp / "plink_pca"
    pgen = tmp / "pgen_pca"
    run([str(PCAONE), "-b", "example/plink", "-k", "3", "-d", "3", "-V", "-n", "4", "-o", str(plink), "-v", "0"])
    run([str(PCAONE), "-p", str(pgen_prefix), "--hardcall", "-k", "3", "-d", "3", "-V", "-n", "4", "-o", str(pgen), "-v", "0"])
    assert_close_vector(out(plink, ".eigvals"), out(pgen, ".eigvals"), 1e-6)
    assert_close_vector(out(plink, ".sigvals"), out(pgen, ".sigvals"), 1e-6)
    assert_mbim_equal(out(plink, ".mbim"), out(pgen, ".mbim"))
    signs = column_signs(out(plink, ".eigvecs"), out(pgen, ".eigvecs"))
    assert_close_matrix_with_signs(out(plink, ".eigvecs"), out(pgen, ".eigvecs"), signs, 2e-5)
    assert_close_matrix_with_signs(out(plink, ".loadings"), out(pgen, ".loadings"), signs, 2e-5)
    return plink


def run_projection_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    plink = tmp / "plink_project"
    pgen = tmp / "pgen_project"
    run([str(PCAONE), "-b", "example/plink", "--USV", str(ref), "--project", "2", "-k", "3", "-o", str(plink), "-v", "0"])
    run([str(PCAONE), "-p", str(pgen_prefix), "--hardcall", "--USV", str(ref), "--project", "2", "-k", "3", "-o", str(pgen), "-v", "0"])
    assert_close_matrix(out(plink, ".eigvecs"), out(pgen, ".eigvecs"), 2e-5)


def run_selection_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    for method, suffixes in (
        ("1", [".galinsky", ".galinsky.pval"]),
        ("2", [".zscore", ".pcadapt", ".pcadapt.chi2", ".pcadapt.pval", ".pcadapt.gif"]),
    ):
        plink = tmp / f"plink_selection_{method}"
        pgen = tmp / f"pgen_selection_{method}"
        run([str(PCAONE), "-b", "example/plink", "--USV", str(ref), "--selection", method, "-k", "3", "-o", str(plink), "-v", "0"])
        run([str(PCAONE), "-p", str(pgen_prefix), "--hardcall", "--USV", str(ref), "--selection", method, "-k", "3", "-o", str(pgen), "-v", "0"])
        for suffix in suffixes:
            assert_close_matrix(out(plink, suffix), out(pgen, suffix), 2e-5)


def run_inbreeding_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    plink = tmp / "plink_inbreed"
    pgen = tmp / "pgen_inbreed"
    run([str(PCAONE), "-b", "example/plink", "--USV", str(ref), "--inbreed", "1", "-k", "3", "-o", str(plink), "-v", "0"])
    run([str(PCAONE), "-p", str(pgen_prefix), "--hardcall", "--USV", str(ref), "--inbreed", "1", "-k", "3", "-o", str(pgen), "-v", "0"])
    assert_hwe_equal(out(plink, ".hwe"), out(pgen, ".hwe"), 2e-5)


def run_ld_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    plink_r2 = tmp / "plink_r2"
    pgen_r2 = tmp / "pgen_r2"
    run([str(PCAONE), "-b", "example/plink", "--USV", str(ref), "--print-r2", "--ld-bp", "1000", "-o", str(plink_r2), "-v", "0"])
    run([str(PCAONE), "-p", str(pgen_prefix), "--hardcall", "--USV", str(ref), "--print-r2", "--ld-bp", "1000", "-o", str(pgen_r2), "-v", "0"])
    assert_gzip_text_equal(out(plink_r2, ".ld.gz"), out(pgen_r2, ".ld.gz"))

    plink_prune = tmp / "plink_prune"
    pgen_prune = tmp / "pgen_prune"
    run([str(PCAONE), "-b", "example/plink", "--USV", str(ref), "--ld-r2", "0.8", "--ld-bp", "1000000", "-o", str(plink_prune), "-v", "0"])
    run([str(PCAONE), "-p", str(pgen_prefix), "--hardcall", "--USV", str(ref), "--ld-r2", "0.8", "--ld-bp", "1000000", "-o", str(pgen_prune), "-v", "0"])
    for suffix in (".ld.prune.in", ".ld.prune.out"):
        if out(plink_prune, suffix).read_text() != out(pgen_prune, suffix).read_text():
            raise AssertionError(f"LD prune mismatch for {suffix}")


def main() -> int:
    if not PCAONE.exists():
        raise SystemExit(f"Missing binary: {PCAONE}")
    with tempfile.TemporaryDirectory(prefix="pcaone-pgen-plink-") as tmpdir:
        tmp = Path(tmpdir)
        pgen_prefix = maybe_make_pgen_from_plink(tmp)
        if pgen_prefix is None:
            return 0

        ref = run_pca_pair(tmp, pgen_prefix)
        run_projection_pair(tmp, pgen_prefix, ref)
        run_selection_pair(tmp, pgen_prefix, ref)
        run_inbreeding_pair(tmp, pgen_prefix, ref)
        run_ld_pair(tmp, pgen_prefix, ref)

    print("PGEN/PLINK equivalence checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
