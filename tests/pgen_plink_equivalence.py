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


def run(
    cmd: list[str], cwd: Path = ROOT, timeout: int = 180
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        cmd, cwd=cwd, text=True, capture_output=True, timeout=timeout
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise AssertionError(
            f"command failed with exit code {result.returncode}: {' '.join(cmd)}"
        )
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
        result = subprocess.run(
            cmd, cwd=ROOT, text=True, capture_output=True, timeout=120
        )
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
            raise AssertionError(
                f"column count mismatch at row {i}: {left}={len(xr)}, {right}={len(yr)}"
            )
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


def assert_close_matrix_with_signs(
    left: Path, right: Path, signs: list[float], tol: float
) -> None:
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
                raise AssertionError(
                    f"mbim metadata mismatch at line {line_no}: {lfields[:6]} != {rfields[:6]}"
                )
            assert_close_scalar(
                float(lfields[6]), float(rfields[6]), tol, f"mbim AF line {line_no}"
            )


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
                raise AssertionError(
                    f"HWE ID mismatch at line {line_no}: {lfields[0]} != {rfields[0]}"
                )
            for col in range(1, 4):
                assert_close_scalar(
                    float(lfields[col]),
                    float(rfields[col]),
                    tol,
                    f"HWE line {line_no} col {col}",
                )


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
    run(
        [
            str(PCAONE),
            "-b",
            "example/plink",
            "-k",
            "3",
            "-d",
            "3",
            "-V",
            "-n",
            "4",
            "-o",
            str(plink),
            "-v",
            "0",
        ]
    )
    run(
        [
            str(PCAONE),
            "-p",
            str(pgen_prefix),
            "--hardcall",
            "-k",
            "3",
            "-d",
            "3",
            "-V",
            "-n",
            "4",
            "-o",
            str(pgen),
            "-v",
            "0",
        ]
    )
    assert_close_vector(out(plink, ".eigvals"), out(pgen, ".eigvals"), 1e-6)
    assert_close_vector(out(plink, ".sigvals"), out(pgen, ".sigvals"), 1e-6)
    assert_mbim_equal(out(plink, ".mbim"), out(pgen, ".mbim"))
    signs = column_signs(out(plink, ".eigvecs"), out(pgen, ".eigvecs"))
    assert_close_matrix_with_signs(
        out(plink, ".eigvecs"), out(pgen, ".eigvecs"), signs, 2e-5
    )
    assert_close_matrix_with_signs(
        out(plink, ".loadings"), out(pgen, ".loadings"), signs, 2e-5
    )
    return plink


def run_projection_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    plink = tmp / "plink_project"
    pgen = tmp / "pgen_project"
    run(
        [
            str(PCAONE),
            "-b",
            "example/plink",
            "--USV",
            str(ref),
            "--project",
            "2",
            "-k",
            "3",
            "-o",
            str(plink),
            "-v",
            "0",
        ]
    )
    run(
        [
            str(PCAONE),
            "-p",
            str(pgen_prefix),
            "--hardcall",
            "--USV",
            str(ref),
            "--project",
            "2",
            "-k",
            "3",
            "-o",
            str(pgen),
            "-v",
            "0",
        ]
    )
    assert_close_matrix(out(plink, ".eigvecs"), out(pgen, ".eigvecs"), 2e-5)


def run_selection_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    for method, suffixes in (
        ("1", [".galinsky", ".galinsky.pval"]),
        (
            "2",
            [".zscore", ".pcadapt", ".pcadapt.chi2", ".pcadapt.pval", ".pcadapt.gif"],
        ),
    ):
        plink = tmp / f"plink_selection_{method}"
        pgen = tmp / f"pgen_selection_{method}"
        run(
            [
                str(PCAONE),
                "-b",
                "example/plink",
                "--USV",
                str(ref),
                "--selection",
                method,
                "-k",
                "3",
                "-o",
                str(plink),
                "-v",
                "0",
            ]
        )
        run(
            [
                str(PCAONE),
                "-p",
                str(pgen_prefix),
                "--hardcall",
                "--USV",
                str(ref),
                "--selection",
                method,
                "-k",
                "3",
                "-o",
                str(pgen),
                "-v",
                "0",
            ]
        )
        for suffix in suffixes:
            assert_close_matrix(out(plink, suffix), out(pgen, suffix), 2e-5)


def write_toy_bed(prefix: Path, n: int, m: int, n_mono: int, seed: int = 7) -> None:
    """A small PLINK trio whose first `n_mono` sites are monomorphic."""
    import random

    rng = random.Random(seed)
    # two populations, so the PCs have something to find
    freqs = [(rng.uniform(0.15, 0.85), rng.uniform(0.15, 0.85)) for _ in range(m)]
    with open(f"{prefix}.bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        for j in range(m):
            codes = []
            for i in range(n):
                if j < n_mono:
                    dose = 0  # every sample homozygous for the same allele
                else:
                    p = freqs[j][0 if i < n // 2 else 1]
                    dose = sum(1 for _ in range(2) if rng.random() < p)
                codes.append({2: 0, 1: 2, 0: 3}[dose])  # PLINK bed coding
            for i in range(0, n, 4):
                byte = 0
                for b, c in enumerate(codes[i : i + 4]):
                    byte |= c << (2 * b)
                f.write(bytes([byte]))
    with open(f"{prefix}.bim", "w") as f:
        for j in range(m):
            f.write(f"1\trs{j}\t0\t{j + 1}\tA\tG\n")
    with open(f"{prefix}.fam", "w") as f:
        for i in range(n):
            f.write(f"F{i} I{i} 0 0 0 -9\n")


def run_selection_zero_variance(tmp: Path) -> None:
    """Sites with no variance must be reported as NA, not as p = 1.

    They have an undefined z-score, so leaving them in put them into the robust
    covariance, the genomic-inflation median and the output as if they were
    ordinary sites. pcadapt drops them (`zscores[pass, ]`) and reports NA.
    """
    n, m, n_mono = 300, 1000, 10
    toy = tmp / "zerovar"
    write_toy_bed(toy, n, m, n_mono)

    ref = tmp / "zerovar_ref"
    # --svd 0: winSVD needs more sites than this toy set has
    run([str(PCAONE), "-b", str(toy), "-k", "2", "--svd", "0", "-o", str(ref), "-v", "0"])

    for method, suffixes in (
        ("1", [".galinsky", ".galinsky.pval"]),
        ("2", [".zscore", ".pcadapt", ".pcadapt.chi2", ".pcadapt.pval"]),
    ):
        target = tmp / f"zerovar_selection_{method}"
        run(
            [
                str(PCAONE), "-b", str(toy), "--USV", str(ref), "--selection", method,
                "-k", "2", "-o", str(target), "-v", "0",
            ]
        )
        for suffix in suffixes:
            rows = [
                line.split()
                for line in out(target, suffix).read_text().splitlines()
                if line and not line.startswith("#")
            ]
            if len(rows) != m:
                raise AssertionError(f"{suffix}: expected {m} rows, got {len(rows)}")
            for i, row in enumerate(rows):
                is_na = all(v == "NA" for v in row)
                any_na = any(v == "NA" for v in row)
                if i < n_mono and not is_na:
                    raise AssertionError(f"{suffix}: site {i} is monomorphic but reads {row}")
                if i >= n_mono and any_na:
                    raise AssertionError(f"{suffix}: site {i} has variance but reads NA")


def run_expect_error(cmd: list[str], needle: str) -> None:
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=180)
    if result.returncode == 0:
        raise AssertionError(f"expected a failure, but the command succeeded: {' '.join(cmd)}")
    blob = result.stdout + result.stderr
    if needle not in blob:
        raise AssertionError(f"expected {needle!r} in the error output, got:\n{blob[-1500:]}")


def run_selection_bad_eigvecs(tmp: Path) -> None:
    """A mismatched .eigvecs must be refused, not read past the end of the matrix.

    read_eigvecs() used to fill M(j, k1) for every line in the file and only
    compare the row count with n afterwards, so a .eigvecs from a bigger cohort
    wrote past the end of the matrix -- a 200x overflow segfaults. Its inner loop
    stopped only on k1 == k, so a line with fewer than K fields ran off the end
    of the string and the run finished with rc=0 and garbage.
    """
    n, m = 200, 800
    toy = tmp / "badusv"
    write_toy_bed(toy, n, m, 0)
    ref = tmp / "badusv_ref"
    run([str(PCAONE), "-b", str(toy), "-k", "3", "--svd", "0", "-o", str(ref), "-v", "0"])

    lines = out(ref, ".eigvecs").read_text().split("\n")
    lines = [ln for ln in lines if ln.strip()]
    cases = {
        "toomany": lines * 2,
        "toofew": lines[: n // 2],
        "narrow": ["\t".join(ln.split()[:2]) for ln in lines],
    }
    for name, rows in cases.items():
        bad = tmp / f"badusv_{name}"
        out(bad, ".eigvecs").write_text("\n".join(rows) + "\n")
        for suffix in (".eigvals", ".sigvals"):
            out(bad, suffix).write_text(out(ref, suffix).read_text())
        for method in ("1", "2"):
            run_expect_error(
                [
                    str(PCAONE), "-b", str(toy), "--USV", str(bad), "--selection", method,
                    "-k", "3", "-o", str(tmp / f"badusv_out_{name}_{method}"), "-v", "0",
                ],
                f"badusv_{name}.eigvecs",
            )


def run_selection_ooc(tmp: Path) -> None:
    """--selection must give the same answer in-core and out-of-core.

    Needs no plink2: it compares a PGEN run against itself. Out-of-core PGEN
    selection used to segfault, because params.perm is on by default while the
    permutation it announces is only built on the PCA path, which --selection
    skips -- FilePgen then indexed an empty permutation. PLINK was unaffected,
    so a PLINK-only check would not have caught it.
    """
    pgen = ROOT / "example" / "plink2.chr1"
    if not Path(f"{pgen}.pgen").exists():
        print("SKIP: example/plink2.chr1.pgen is not available")
        return

    ref = tmp / "ooc_ref"
    run([str(PCAONE), "-p", str(pgen), "-k", "3", "-o", str(ref), "-v", "0"])

    for method, suffixes in (
        ("1", [".galinsky", ".galinsky.pval"]),
        (
            "2",
            [".zscore", ".pcadapt", ".pcadapt.chi2", ".pcadapt.pval", ".pcadapt.gif"],
        ),
    ):
        incore = tmp / f"ooc_selection_{method}_incore"
        blocked = tmp / f"ooc_selection_{method}_blocked"
        for target, extra in ((incore, []), (blocked, ["-m", "0.002"])):
            run(
                [
                    str(PCAONE),
                    "-p",
                    str(pgen),
                    "--USV",
                    str(ref),
                    "--selection",
                    method,
                    "-k",
                    "3",
                    "-o",
                    str(target),
                    "-v",
                    "0",
                ]
                + extra
            )
        for suffix in suffixes:
            assert_close_matrix(out(incore, suffix), out(blocked, suffix), 0.0)


def run_inbreeding_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    plink = tmp / "plink_inbreed"
    pgen = tmp / "pgen_inbreed"
    run(
        [
            str(PCAONE),
            "-b",
            "example/plink",
            "--USV",
            str(ref),
            "--inbreed",
            "1",
            "-k",
            "3",
            "-o",
            str(plink),
            "-v",
            "0",
        ]
    )
    run(
        [
            str(PCAONE),
            "-p",
            str(pgen_prefix),
            "--hardcall",
            "--USV",
            str(ref),
            "--inbreed",
            "1",
            "-k",
            "3",
            "-o",
            str(pgen),
            "-v",
            "0",
        ]
    )
    assert_hwe_equal(out(plink, ".hwe"), out(pgen, ".hwe"), 2e-5)


def run_ld_pair(tmp: Path, pgen_prefix: Path, ref: Path) -> None:
    plink_r2 = tmp / "plink_r2"
    pgen_r2 = tmp / "pgen_r2"
    run(
        [
            str(PCAONE),
            "-b",
            "example/plink",
            "--USV",
            str(ref),
            "--print-r2",
            "--ld-bp",
            "1000",
            "-o",
            str(plink_r2),
            "-v",
            "0",
        ]
    )
    run(
        [
            str(PCAONE),
            "-p",
            str(pgen_prefix),
            "--hardcall",
            "--USV",
            str(ref),
            "--print-r2",
            "--ld-bp",
            "1000",
            "-o",
            str(pgen_r2),
            "-v",
            "0",
        ]
    )
    assert_gzip_text_equal(out(plink_r2, ".ld.gz"), out(pgen_r2, ".ld.gz"))

    plink_prune = tmp / "plink_prune"
    pgen_prune = tmp / "pgen_prune"
    run(
        [
            str(PCAONE),
            "-b",
            "example/plink",
            "--USV",
            str(ref),
            "--ld-r2",
            "0.8",
            "--ld-bp",
            "1000000",
            "-o",
            str(plink_prune),
            "-v",
            "0",
        ]
    )
    run(
        [
            str(PCAONE),
            "-p",
            str(pgen_prefix),
            "--hardcall",
            "--USV",
            str(ref),
            "--ld-r2",
            "0.8",
            "--ld-bp",
            "1000000",
            "-o",
            str(pgen_prune),
            "-v",
            "0",
        ]
    )
    for suffix in (".ld.prune.in", ".ld.prune.out"):
        if out(plink_prune, suffix).read_text() != out(pgen_prune, suffix).read_text():
            raise AssertionError(f"LD prune mismatch for {suffix}")


def main() -> int:
    if not PCAONE.exists():
        raise SystemExit(f"Missing binary: {PCAONE}")
    with tempfile.TemporaryDirectory(prefix="pcaone-pgen-plink-") as tmpdir:
        tmp = Path(tmpdir)
        run_selection_ooc(tmp)
        run_selection_zero_variance(tmp)
        run_selection_bad_eigvecs(tmp)

        pgen_prefix = maybe_make_pgen_from_plink(tmp)
        if pgen_prefix is None:
            print("PGEN out-of-core, zero-variance and bad-.eigvecs selection checks passed")
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
