#!/usr/bin/env python3
"""BGEN input: every encoding of the same genotypes gives the same PCs.

One panel of hard-called genotypes (three populations, 2% missing calls) is
written as BGEN layout 2 with no, zlib and (when a zstd module is available)
zstd compression, at 8, 12 and 16 bits, phased and unphased, and as layout 1
with no and zlib compression. 0 and 1 are exact at any bit depth, so every
file holds the same dosages, and each must give the bytes of .eigvals,
.eigvecs and .loadings of the reference file (layout 2, zlib, 8 bits), run the
same way: in-core, with -m (the variants of a block read in the shuffled
order) and with -m -S. The out-of-core PCs agree with the in-core ones.

A phased file is read as the sum of the haplotypes' allele probabilities; the
bgen library took the second haplotype for a heterozygote and gave a dosage of
1 to 0|1 but of 2 to 1|0.

Standard library only.

    python3 tests/test_bgen_formats.py
"""

from __future__ import annotations

import math
import random
import struct
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PCAONE = ROOT / "PCAone"

try:
    from compression import zstd as _zstd  # Python 3.14

    def zstd_compress(b: bytes) -> bytes:
        return _zstd.compress(b)
except ImportError:
    try:
        import zstandard

        def zstd_compress(b: bytes) -> bytes:
            return zstandard.ZstdCompressor().compress(b)
    except ImportError:
        zstd_compress = None

N, M, K = 150, 3000, 3


def simulate(seed: int = 7):
    """haplotypes (copies of the second allele) and missingness, per variant"""
    rng = random.Random(seed)
    pop = [i % K for i in range(N)]
    variants = []
    for _ in range(M):
        p = rng.uniform(0.05, 0.95)
        f = 0.1
        pk = [min(max(rng.betavariate(p * (1 - f) / f, (1 - p) * (1 - f) / f), 1e-3), 1 - 1e-3) for _ in range(K)]
        h = [(int(rng.random() < pk[pop[i]]), int(rng.random() < pk[pop[i]])) for i in range(N)]
        miss = [rng.random() < 0.02 for _ in range(N)]
        variants.append((h, miss))
    return variants


def pack(vals, bits):
    acc = nbits = 0
    out = bytearray()
    for v in vals:
        acc |= v << nbits
        nbits += bits
        while nbits >= 8:
            out.append(acc & 0xFF)
            acc >>= 8
            nbits -= 8
    if nbits:
        out.append(acc & 0xFF)
    return bytes(out)


def write_bgen(path: Path, variants, layout=2, compression="zlib", bits=8, phased=False):
    comp = {"none": 0, "zlib": 1, "zstd": 2}[compression]
    squeeze = {"none": lambda b: b, "zlib": lambda b: zlib.compress(b, 6), "zstd": zstd_compress}[compression]
    with open(path, "wb") as f:
        f.write(struct.pack("<I", 20))
        f.write(struct.pack("<III", 20, len(variants), N) + b"bgen" + struct.pack("<I", comp | (layout << 2)))
        for j, (h, miss) in enumerate(variants):
            vid = f"snp{j + 1}".encode()
            rec = bytearray()
            if layout == 1:
                rec += struct.pack("<I", N)
            for s in (vid, vid, b"1"):
                rec += struct.pack("<H", len(s)) + s
            rec += struct.pack("<I", j + 1)
            if layout == 2:
                rec += struct.pack("<H", 2)
            for a in (b"A", b"G"):
                rec += struct.pack("<I", len(a)) + a
            g = [a + b for a, b in h]
            if layout == 1:
                probs = bytearray()
                for i in range(N):
                    v = [0, 0, 0] if miss[i] else [32768 * (g[i] == c) for c in range(3)]
                    probs += struct.pack("<HHH", *v)
                if compression == "none":
                    rec += probs
                else:
                    c = squeeze(bytes(probs))
                    rec += struct.pack("<I", len(c)) + c
            else:
                mx = (1 << bits) - 1
                vals = []
                for i in range(N):
                    if phased:  # P(first allele) of each haplotype
                        vals += [0, 0] if miss[i] else [mx * (1 - h[i][0]), mx * (1 - h[i][1])]
                    else:  # P(0 and 1 copies of the second allele)
                        vals += [0, 0] if miss[i] else [mx * (g[i] == 0), mx * (g[i] == 1)]
                ploidy = bytes(2 | (0x80 if m else 0) for m in miss)
                raw = struct.pack("<IHBB", N, 2, 2, 2) + ploidy + struct.pack("<BB", int(phased), bits) + pack(vals, bits)
                if compression == "none":
                    rec += struct.pack("<I", len(raw)) + raw
                else:
                    c = squeeze(raw)
                    rec += struct.pack("<II", len(c) + 4, len(raw)) + c
            f.write(rec)


def run(bgen: Path, out: Path, *args) -> None:
    cmd = [str(PCAONE), "-g", str(bgen), "-k", "3", "-n", "2", "-o", str(out), *map(str, args)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"FAIL: {' '.join(cmd)}\n{r.stdout}\n{r.stderr}")


def read(path: Path):
    return [[float(x) for x in line.split()] for line in path.read_text().split("\n") if line.strip()]


def same_pcs(a: Path, b: Path, tol: float) -> str | None:
    ea, eb = read(Path(f"{a}.eigvals")), read(Path(f"{b}.eigvals"))
    for x, y in zip(ea, eb):
        if abs(x[0] - y[0]) > tol * abs(y[0]):
            return f"eigvals {x[0]} vs {y[0]}"
    va, vb = read(Path(f"{a}.eigvecs")), read(Path(f"{b}.eigvecs"))
    for k in range(3):
        x = [r[k] for r in va]
        y = [r[k] for r in vb]
        dot = sum(p * q for p, q in zip(x, y))
        cos = abs(dot) / math.sqrt(sum(p * p for p in x) * sum(q * q for q in y))
        if cos < 1 - tol:
            return f"PC{k + 1} |cos| = {cos}"
    return None


def main() -> int:
    if not PCAONE.exists():
        print(f"build {PCAONE} first")
        return 1
    variants = simulate()
    encodings = [
        ("l2-zlib-8", dict()),
        ("l2-none-8", dict(compression="none")),
        ("l2-zlib-12", dict(bits=12)),
        ("l2-zlib-16", dict(bits=16)),
        ("l2-zlib-8-phased", dict(phased=True)),
        ("l2-none-16-phased", dict(compression="none", bits=16, phased=True)),
        ("l1-zlib", dict(layout=1)),
        ("l1-none", dict(layout=1, compression="none")),
    ]
    if zstd_compress:
        encodings += [("l2-zstd-8", dict(compression="zstd")), ("l2-zstd-16-phased", dict(compression="zstd", bits=16, phased=True))]
    else:
        print("skip zstd: no zstd module")
    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for name, kw in encodings:
            write_bgen(tmp / f"{name}.bgen", variants, **kw)
        modes = (("ic", []), ("ooc", ["-m", "0.001"]), ("ooc-S", ["-m", "0.001", "-S"]))
        for mode, args in modes:
            ref = tmp / f"l2-zlib-8.{mode}"
            run(tmp / "l2-zlib-8.bgen", ref, *args)
            for name, _ in encodings[1:]:
                out = tmp / f"{name}.{mode}"
                run(tmp / f"{name}.bgen", out, *args)
                # the same dosages, the same run: the same bytes
                differ = [suf for suf in (".eigvals", ".eigvecs", ".loadings")
                          if Path(f"{out}{suf}").read_bytes() != Path(f"{ref}{suf}").read_bytes()]
                err = f"{' '.join(differ)} differ from the reference" if differ else None
                print(f"{'FAIL' if err else 'ok  '} {name:20s} {mode:6s} {err or ''}")
                failures += err is not None
        # the out-of-core runs read the blocks (in a shuffled order with -m)
        # and stop at their own tolerance: the same PCs, not the same bytes
        for mode in ("ooc", "ooc-S"):
            err = same_pcs(tmp / f"l2-zlib-8.{mode}", tmp / "l2-zlib-8.ic", 1e-2)
            print(f"{'FAIL' if err else 'ok  '} {'in-core vs ' + mode:27s} {err or ''}")
            failures += err is not None
    print("all passed" if not failures else f"{failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
