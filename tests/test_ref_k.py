#!/usr/bin/env python3
"""-k/--pc picks the leading PCs of the -P/--USV reference in every two-stage
analysis: LD (-R, --ld-r2), --evaladmix, --project, --selection and --inbreed.

For each one, -k 3 on a 5-PC reference must give the same output as a copy of
that reference cut to its first 3 PCs, no -k must equal -k 5 (all PCs), and
-k 6 must be refused.

    python3 tests/test_ref_k.py
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_crash_regressions import PCAONE, make_plink  # noqa: E402

failures: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(("ok  " if ok else "FAIL ") + name + (f": {detail}" if detail and not ok else ""))
    if not ok:
        failures.append(name)


def run(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run([PCAONE, *[str(a) for a in args], "-n", "2"], capture_output=True, text=True)


def outputs(prefix: Path) -> dict[str, bytes]:
    """every output of a run but its log; a .gz is compared decompressed"""
    res = {}
    for p in prefix.parent.glob(prefix.name + ".*"):
        suffix = p.name[len(prefix.name):]
        if suffix == ".log":
            continue
        res[suffix] = gzip.decompress(p.read_bytes()) if suffix.endswith(".gz") else p.read_bytes()
    return res


def cut_reference(src: Path, dst: Path, k: int) -> None:
    """a copy of the reference PCA src with only its first k PCs"""
    for suf in (".eigvecs", ".loadings"):
        rows = Path(f"{src}{suf}").read_text().splitlines()
        Path(f"{dst}{suf}").write_text("".join("\t".join(r.split("\t")[:k]) + "\n" for r in rows))
    Path(f"{dst}.eigvals").write_text("".join(line + "\n" for line in Path(f"{src}.eigvals").read_text().splitlines()[:k]))
    sig = Path(f"{src}.sigvals").read_text().splitlines()  # a #N,M,... header, then one value per PC
    Path(f"{dst}.sigvals").write_text("".join(line + "\n" for line in sig[: k + 1]))
    shutil.copy(f"{src}.mbim", f"{dst}.mbim")


def main() -> int:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        t, _ = make_plink(tmp)
        ref, cut = tmp / "ref5", tmp / "cut3"
        p = run(["-b", t, "-k", "5", "-o", ref])
        check("reference_pca", p.returncode == 0, p.stderr[-300:])
        cut_reference(ref, cut, 3)

        analyses = {
            "ld_r2": ["-R", "--ld-bp", "1000"],
            "ld_prune": ["--ld-r2", "0.2", "--ld-bp", "1000"],
            "evaladmix": ["--evaladmix"],
            "project": ["--project", "2"],
            "selection1": ["--selection", "1"],
            "selection2": ["--selection", "2"],
            "inbreed": ["--inbreed", "1"],
        }
        for name, extra in analyses.items():
            o = tmp / name
            runs = {
                "k3": run(["-b", t, *extra, "-P", ref, "-k", "3", "-o", f"{o}_k3"]),
                "cut3": run(["-b", t, *extra, "-P", cut, "-o", f"{o}_cut3"]),
                "k5": run(["-b", t, *extra, "-P", ref, "-k", "5", "-o", f"{o}_k5"]),
                "all": run(["-b", t, *extra, "-P", ref, "-o", f"{o}_all"]),
            }
            bad = {k: (p.stdout + p.stderr)[-300:] for k, p in runs.items() if p.returncode != 0}
            check(f"{name}_runs", not bad, str(bad))
            if bad:
                continue
            k3, cut3 = outputs(Path(f"{o}_k3")), outputs(Path(f"{o}_cut3"))
            k5, all_ = outputs(Path(f"{o}_k5")), outputs(Path(f"{o}_all"))
            check(f"{name}_leading_pcs", k3 and k3 == cut3,
                  f"-k 3 differs from the reference cut to 3 PCs in {sorted(s for s in k3 if k3[s] != cut3.get(s))}")
            check(f"{name}_default_all", all_ and all_ == k5, "no -k differs from -k 5")
            check(f"{name}_k_matters", k3 != all_, "-k 3 gives the output of all 5 PCs")
            p = run(["-b", t, *extra, "-P", ref, "-k", "6", "-o", f"{o}_k6"])
            check(f"{name}_too_many", p.returncode != 0 and "larger than the 5 PCs" in p.stdout + p.stderr,
                  (p.stdout + p.stderr)[-300:])

    if failures:
        print("FAILED:", ", ".join(failures))
        return 1
    print("SUCCESS: -k selects the leading reference PCs in every two-stage analysis.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
