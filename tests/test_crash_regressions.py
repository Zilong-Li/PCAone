#!/usr/bin/env python3
"""Regression tests for inputs that used to crash PCAone (segfault, abort,
SIGFPE, std::bad_alloc or an uncaught exception).

Each case must now either succeed or stop with a clean "Error:" message and
exit code 1. Only the standard library is needed: the PLINK files are written
here. The CSV cases need the `zstd` command or the `zstandard` module, and the
PGEN case needs `plink2`; each is skipped when missing.

    python3 tests/test_crash_regressions.py
"""

from __future__ import annotations

import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PCAONE = str(ROOT / "PCAone")
CRASH_TEXT = ("terminate called", "Segmentation", "bad_alloc", "Aborted", "double free", "corrupted")
failures: list[str] = []


def run(name: str, args: list[str], expect_ok: bool, needle: str | None = None) -> str:
    """Run PCAone; fail on any crash, or on a result other than the expected one."""
    p = subprocess.run([PCAONE, *args, "-n", "2"], capture_output=True, text=True)
    out = p.stdout + p.stderr
    crashed = p.returncode < 0 or p.returncode >= 128 or any(t in out for t in CRASH_TEXT)
    if crashed:
        failures.append(f"{name}: crashed (rc={p.returncode})")
    elif expect_ok and p.returncode != 0:
        failures.append(f"{name}: failed (rc={p.returncode}): {out.strip().splitlines()[-1:]}")
    elif not expect_ok and p.returncode != 1:
        failures.append(f"{name}: expected a clean error (rc=1), got rc={p.returncode}")
    elif needle is not None and needle not in out:
        failures.append(f"{name}: message lacks {needle!r}")
    else:
        print(f"ok  {name}")
    return out


def write_bed(prefix: Path, geno: list[list[int]], bim: list[tuple[str, int]], n: int) -> None:
    """geno[j][i] = copies of A1 (0/1/2) or -1 for missing; PLINK SNP-major."""
    code = {2: 0b00, -1: 0b01, 1: 0b10, 0: 0b11}
    with open(f"{prefix}.bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        for row in geno:
            out = bytearray((n + 3) // 4)
            for i, g in enumerate(row):
                out[i // 4] |= code[g] << (2 * (i % 4))
            f.write(out)
    with open(f"{prefix}.bim", "w") as f:
        for j, (chrom, bp) in enumerate(bim):
            f.write(f"{chrom}\tsnp{j}\t0\t{bp}\tA\tG\n")
    with open(f"{prefix}.fam", "w") as f:
        for i in range(n):
            f.write(f"ind{i} ind{i} 0 0 0 -9\n")


def make_plink(tmp: Path, n: int = 60, per_chr: int = 300) -> tuple[Path, list[tuple[str, int]]]:
    rng = random.Random(11)
    bim = [("1", 1000 + 100 * j) for j in range(per_chr)] + [("2", 500 + 100 * j) for j in range(per_chr)]
    pops = [i % 3 for i in range(n)]
    geno = []
    for _ in bim:
        f = [min(max(rng.gauss(0.4, 0.15), 0.05), 0.95) for _ in range(3)]
        row = []
        for i in range(n):
            if rng.random() < 0.02:
                row.append(-1)
            else:
                p = f[pops[i]]
                row.append((rng.random() < p) + (rng.random() < p))
        geno.append(row)
    prefix = tmp / "t"
    write_bed(prefix, geno, bim, n)
    # the first 300 sites only, to pair with a reference made from all 600
    write_bed(tmp / "t_small", geno[:per_chr], bim[:per_chr], n)
    return prefix, bim


def zstd_compress(src: Path, dst: Path) -> bool:
    if shutil.which("zstd"):
        subprocess.run(["zstd", "-q", "-f", str(src), "-o", str(dst)], check=True)
        return True
    try:
        import zstandard  # type: ignore
    except ImportError:
        return False
    dst.write_bytes(zstandard.ZstdCompressor().compress(src.read_bytes()))
    return True


def main() -> int:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        t, bim = make_plink(tmp)
        o = str(tmp / "o")

        # winSVD with fewer than bands^2 sites: the trailing blocks are empty,
        # and their size underflowed to ~2^64 (std::bad_alloc)
        run("winsvd_small_m", ["-b", str(t), "-k", "2", "-o", o], True)
        run("winsvd_w1024", ["-b", str(t), "-k", "2", "-w", "1024", "-o", o], True)
        run("reference", ["-b", str(t), "-k", "2", "-V", "-o", str(tmp / "ref")], True)

        # clumping over two chromosomes. The first row of chr2 is the top hit;
        # it was dropped in-core, and out-of-core the extra column per
        # chromosome wrote past the matrix (segfault).
        rows = ["CHR\tBP\tP"]
        for j, (chrom, bp) in enumerate(bim):
            p = 1e-12 if j == 300 else (1e-8 if j == 10 else 0.5)
            rows.append(f"{chrom}\t{bp}\t{p:g}")
        (tmp / "assoc.txt").write_text("\n".join(rows) + "\n")
        rows_na = rows.copy()
        rows_na[5] = rows_na[5].rsplit("\t", 1)[0] + "\tNA"
        (tmp / "assoc_na.txt").write_text("\n".join(rows_na) + "\n")
        clump = ["-b", str(t), "--ld-stats", "1", "--clump-r2", "0.9"]
        run("clump_incore", clump + ["--clump", str(tmp / "assoc.txt"), "-o", str(tmp / "c0")], True)
        run("clump_ooc", clump + ["--clump", str(tmp / "assoc.txt"), "-m", "0.0001", "-o", str(tmp / "c1")], True)
        c0 = (tmp / "c0.p0.clump").read_text() if (tmp / "c0.p0.clump").exists() else ""
        c1 = (tmp / "c1.p0.clump").read_text() if (tmp / "c1.p0.clump").exists() else None
        if "\t500\t1e-12" not in c0:
            failures.append("clump_incore: the first row of chr2 (the top hit) is missing")
        if c1 is not None and c1 != c0:
            failures.append("clump_ooc: differs from the in-core result")
        run("clump_na_pvalue", clump + ["--clump", str(tmp / "assoc_na.txt"), "-o", o], True)

        # duplicated position at the end of the last chromosome
        dup = tmp / "dup"
        for suf in (".bed", ".fam"):
            shutil.copy(f"{t}{suf}", f"{dup}{suf}")
        lines = Path(f"{t}.bim").read_text().splitlines()
        last = lines[-1].split("\t")
        last[3] = lines[-2].split("\t")[3]
        lines[-1] = "\t".join(last)
        Path(f"{dup}.bim").write_text("\n".join(lines) + "\n")
        run("ld_dup_position", ["-b", str(dup), "--ld-stats", "1", "-R", "--ld-bp", "1000", "-o", o], True)

        # --inbreed on a target with fewer sites than the reference (segfault)
        run("inbreed_wrong_ref", ["-b", str(tmp / "t_small"), "-P", str(tmp / "ref"), "--inbreed", "1", "-k", "2",
                                  "-o", o], False, "sites")

        # a missing input used to abort through an uncaught exception
        run("missing_input", ["-b", str(tmp / "nope"), "-o", o], False, "Error:")

        # input types a mode cannot read are refused up front
        run("pcangsd_pgen", ["-p", str(t), "--pcangsd", "-o", o], False, "--pcangsd supports only")
        run("pcangsd_bgen", ["-g", str(t), "--pcangsd", "-o", o], False, "--pcangsd supports only")
        run("emu_bgen_ooc", ["-g", str(t), "--emu", "-m", "1", "-o", o], False, "--emu with --bgen")
        run("project_bgen", ["-g", str(t), "-P", str(tmp / "ref"), "--project", "1", "-o", o], False,
            "--project supports only")
        run("evaladmix_bgen", ["-g", str(t), "--evaladmix", "-o", o], False, "--evaladmix supports only")

        # PGEN: out-of-core --inbreed indexed an empty permutation (segfault)
        plink2 = shutil.which("plink2")
        if plink2:
            subprocess.run([plink2, "--bfile", str(t), "--make-pgen", "--out", str(tmp / "tp"), "--silent"], check=True)
            run("pgen_reference", ["-p", str(tmp / "tp"), "-k", "2", "-V", "-o", str(tmp / "pref")], True)
            run("pgen_inbreed_ooc", ["-p", str(tmp / "tp"), "-P", str(tmp / "pref"), "--inbreed", "1", "-k", "2",
                                     "-m", "0.0001", "-o", o], True)
        else:
            print("skip pgen cases: plink2 not found")

        # CSV
        rng = random.Random(5)
        mat = [[rng.randint(0, 9) for _ in range(20)] for _ in range(300)]
        def csv(name: str, lines: list[str]) -> Path | None:
            raw = tmp / f"{name}.csv"
            raw.write_text("".join(line + "\n" for line in lines))
            dst = tmp / f"{name}.csv.zst"
            return dst if zstd_compress(raw, dst) else None
        good = csv("good", [",".join(map(str, r)) for r in mat])
        if good is None:
            print("skip csv cases: neither zstd nor zstandard found")
        else:
            run("csv_iram", ["-c", str(good), "-k", "2", "-d", "0", "-o", o], True)
            run("csv_buffer4", ["-c", str(good), "-k", "2", "-m", "0.0001", "--buffer", "4", "-o", o], True)
            run("csv_M_smaller", ["-c", str(good), "-k", "2", "--M", "200", "--N", "20", "-o", o], False, "lines")
            run("csv_N_smaller", ["-c", str(good), "-k", "2", "--M", "300", "--N", "10", "-o", o], False, "columns")
            run("csv_emu", ["-c", str(good), "--emu", "-o", o], False, "--emu supports only")
            wide = [",".join(map(str, r)) for r in mat]
            wide[1] += ",1,2"  # the second line was never checked
            run("csv_wide_line2", ["-c", str(csv("wide", wide)), "-k", "2", "-o", o], False, "columns")
            text = [",".join(map(str, r)) for r in mat]
            text[7] = "NA," + text[7].split(",", 1)[1]
            run("csv_non_numeric", ["-c", str(csv("text", text)), "-k", "2", "-o", o], False, "not a number")
            run("csv_non_numeric_ooc", ["-c", str(tmp / "text.csv.zst"), "-k", "2", "-m", "0.0001", "-o", o], False,
                "not a number")
            run("csv_empty", ["-c", str(csv("empty", [])), "-k", "2", "-o", o], False, "empty")

    if failures:
        print("FAILED:")
        for f in failures:
            print("  " + f)
        return 1
    print("SUCCESS: no crashes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
