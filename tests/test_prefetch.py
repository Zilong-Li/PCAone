#!/usr/bin/env python3
"""Out-of-core reads with the next block read in the background (default) give
the same bytes as without (--no-prefetch), so the same PCs, for every reader and
method; and the hand-threaded H += X * G of winSVD gives the same PCs for any
number of threads."""
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
    # the row panels of H += X * G: the PCs do not depend on the number of threads
    one, three = directory / 'n1', directory / 'n3'
    pcaone(*bed, '-n', 1, '-o', one)
    pcaone(*bed, '-n', 3, '-o', three)
    same_outputs(one, three, ('.eigvals', '.eigvecs'))
    checked.append('threads')

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

    # CSV, read through its binary copy
    import random
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
