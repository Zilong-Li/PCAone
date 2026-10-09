#!/usr/bin/env python3
"""Out-of-core reads with the next block read in the background (default) give
the same bytes as without (--no-prefetch), so the same PCs, for every reader and
method; and the row-panel H = X * G and H += X * G give the PCs of one thread
for any number of threads, with enough samples for the panels to be used."""
import random
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# reuse the BED fixture of the permutation test
source_text = (ROOT / 'tests' / 'test_bed_permutation.py').read_text()
helpers = source_text[:source_text.index('with tempfile.TemporaryDirectory')]
ns = {'__file__': str(ROOT / 'tests' / 'test_bed_permutation.py')}
exec(compile(helpers, 'test_bed_permutation.py', 'exec'), ns)
fixture = ns['fixture']


def pcaone(*args):
    result = subprocess.run([str(ROOT / 'PCAone'), *map(str, args)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def same_outputs(a, b, suffixes=('.eigvals', '.eigvecs', '.loadings')):
    for s in suffixes:
        pa, pb = Path(str(a) + s), Path(str(b) + s)
        if pa.exists() or pb.exists():
            assert pa.read_bytes() == pb.read_bytes(), (a, b, s)


def write_bed(prefix, n, m, seed):
    """n samples in two populations, m SNPs, 2 % missing"""
    rng = random.Random(seed)
    width = (n + 3) // 4
    bed = bytearray(b'\x6c\x1b\x01')
    for j in range(m):
        codes = []
        for i in range(n):
            p = 0.05 + (j % 9) * 0.05 + (0.3 if i < n // 3 else 0)
            g = int(rng.random() < p) + int(rng.random() < p)
            codes.append(1 if rng.random() < 0.02 else {0: 3, 1: 2, 2: 0}[g])
        bed += bytes(sum(codes[i + k] << (2 * k) for k in range(min(4, n - i))) for i in range(0, n, 4))
    assert len(bed) == 3 + m * width
    Path(str(prefix) + '.bed').write_bytes(bytes(bed))
    Path(str(prefix) + '.bim').write_text(''.join(f'1\trs{j}\t0\t{j + 1}\tA\tC\n' for j in range(m)))
    Path(str(prefix) + '.fam').write_text(''.join(f'F{i} I{i} 0 0 0 -9\n' for i in range(n)))
    return prefix


checked = []
with tempfile.TemporaryDirectory(prefix='pcaone-prefetch-') as name:
    directory = Path(name)
    source, *_ = fixture(directory)

    def both(label, *args):
        """the same run with and without the background reads"""
        on, off = directory / (label + '-on'), directory / (label + '-off')
        pcaone(*args, '-o', on)
        pcaone(*args, '--no-prefetch', '-o', off)
        same_outputs(on, off)
        checked.append(label)
        return on

    bed = ('--bfile', source, '-k', 3, '-m', '0.00001')
    both('winsvd', *bed, '-n', 2)
    both('winsvd-noshuffle', *bed, '-n', 2, '-S')
    both('emu', *bed, '-n', 2, '--emu', '--maxiter', 3)
    both('ssvd', *bed, '-n', 2, '--svd', 1)
    both('iram', *bed, '-n', 2, '--svd', 0)
    both('exact', *bed, '-n', 2, '--svd', 3)
    # The row panels of H = X * G and H += X * G (mul_X_Y) start at 192
    # samples; the fixture has 41. With 300 samples and 3 threads, X is cut into
    # two panels (0-143, 144-299: the last 12 rows join the second panel), each
    # multiplied on one thread, and the PCs must be those of one thread: the
    # in-core sSVD assigns H, the out-of-core winSVD accumulates block by block.
    big = write_bed(directory / 'big', 300, 1500, seed=5)
    for label, args in [('threads-insvd', ('--svd', 1)),
                        ('threads-winsvd', ('-m', '0.0005')),
                        ('threads-winsvd-noshuffle', ('-m', '0.0005', '-S'))]:
        one, three = directory / (label + '-n1'), directory / (label + '-n3')
        pcaone('--bfile', big, '-k', 4, *args, '-n', 1, '-o', one)
        pcaone('--bfile', big, '-k', 4, *args, '-n', 3, '-o', three)
        same_outputs(one, three)
        checked.append(label)

    # PGEN (pgenlib's readers in place of PgenReader's), when the example is there
    pgen = ROOT / 'example' / 'plink2.chr1'
    if Path(str(pgen) + '.pgen').exists():
        pg = ('--pgen', pgen, '-k', 3, '-m', '0.001', '-n', 2)
        both('pgen-dosage', *pg)
        both('pgen-hardcall', *pg, '--hardcall')
        both('pgen-emu', *pg, '--hardcall', '--emu', '--maxiter', 2)
        both('pgen-ssvd', *pg, '--svd', 1)
    else:
        print('skip PGEN: example/plink2.chr1.pgen is not available')

    # CSV, read through its binary copy (next block requested into the page cache)
    rng = random.Random(3)
    raw = directory / 'counts.csv'
    raw.write_text(''.join(','.join(str(rng.randint(0, 9)) for _ in range(30)) + '\n' for _ in range(500)))
    zst = Path(str(raw) + '.zst')
    if shutil.which('zstd'):
        subprocess.run(['zstd', '-q', '-f', str(raw), '-o', str(zst)], check=True)
    else:
        try:
            import zstandard
            zst.write_bytes(zstandard.ZstdCompressor().compress(raw.read_bytes()))
        except ImportError:
            zst = None
    if zst:
        both('csv', '-c', zst, '-k', 3, '-m', '0.00002', '-n', 2)
    else:
        print('skip CSV: neither zstd nor zstandard found')

    # BGEN (read ahead into the page cache), when the bgen module can write one
    try:
        import numpy as np
        from bgen import BgenWriter
    except ImportError:
        print('skip BGEN: the python bgen module is not available')
    else:
        path = directory / 'small.bgen'
        rng = np.random.default_rng(4)
        with BgenWriter(str(path), 40, samples=[f'S{i}' for i in range(40)], compression='zstd', layout=2) as bw:
            for j in range(600):
                p = rng.dirichlet([1, 1, 1], size=40)
                bw.add_variant(f'v{j}', f'v{j}', '1', j + 1, ['A', 'C'], p, ploidy=2, phased=False, bit_depth=8)
        both('bgen', '--bgen', path, '-k', 3, '-m', '0.00002', '-n', 2)

print('prefetch: identical outputs with and without background reads:', ', '.join(checked))
