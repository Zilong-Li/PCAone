#!/usr/bin/env python3
"""Out-of-core BED permutation: bytes and metadata, --seed, the invariance of
winSVD/EMU to the SNP order within a -w band, and the output guards."""
import math
import os
import random
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
N, M = 41, 1003  # partial final blocks, missing values, and population structure
WIDTH = (N + 3) // 4


def plink_file(prefix, suffix):
    return Path(str(prefix) + suffix)  # with_suffix would replace the '.perm' of 'x.perm'


def write_plink(prefix, bed, bim, fam):
    plink_file(prefix, '.bed').write_bytes(bed)
    plink_file(prefix, '.bim').write_text(''.join(line + '\n' for line in bim))
    plink_file(prefix, '.fam').write_bytes(fam)
    return prefix


def fixture(directory):
    rng = random.Random(192)
    bed = bytearray(b'\x6c\x1b\x01')
    for j in range(M):
        calls = []
        for i in range(N):
            p = 0.1 + (j % 5) * 0.1 + (0.25 if i < N // 2 else 0)
            g = int(rng.random() < p) + int(rng.random() < p)
            calls.append(1 if rng.random() < 0.025 else {0: 3, 1: 2, 2: 0}[g])
        for i in range(0, N, 4):
            bed.append(sum(calls[i + k] << (2 * k) for k in range(min(4, N - i))))
    bim = [f'1\trs{j}\t0\t{j + 1}\tA\tC' for j in range(M)]
    fam = ''.join(f'F{i} I{i} 0 0 0 -9\n' for i in range(N)).encode()
    return write_plink(directory / 'source', bytes(bed), bim, fam), bytes(bed), bim, fam


def reordered(directory, label, rows, bed, bim, fam):
    """A copy of the source with its SNPs in the order `rows` (source indices)."""
    body = b''.join(bed[3 + i * WIDTH:3 + (i + 1) * WIDTH] for i in rows)
    return write_plink(directory / label, bed[:3] + body, [bim[i] for i in rows], fam)


def matrix(path):
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()
            if line and not line.startswith('#')]


def compare(left, right, right_rows):
    """PCs of `left` (a permuted run, loadings in source order) and `right`
    (an unpermuted run whose row r of the loadings is source SNP right_rows[r])."""
    for suffix in ['eigvals', 'sigvals']:
        a, b = matrix(Path(str(left) + '.' + suffix)), matrix(Path(str(right) + '.' + suffix))
        for ra, rb in zip(a, b):
            for x, y in zip(ra, rb):
                assert math.isclose(x, y, rel_tol=2e-5, abs_tol=2e-7), (suffix, x, y)
    # printed with limited precision; the sample Gram is sign and rotation invariant
    u, v = matrix(Path(str(left) + '.eigvecs')), matrix(Path(str(right) + '.eigvecs'))
    for i in range(N):
        for j in range(N):
            a = sum(x * y for x, y in zip(u[i], u[j]))
            b = sum(x * y for x, y in zip(v[i], v[j]))
            assert abs(a - b) < 2e-5, ('PC subspace', a, b)
    a, rows = matrix(Path(str(left) + '.loadings')), matrix(Path(str(right) + '.loadings'))
    b = [None] * M
    for r, source in enumerate(right_rows):
        b[source] = rows[r]
    for pc in range(len(u[0])):
        sign = 1 if sum(x[pc] * y[pc] for x, y in zip(u, v)) >= 0 else -1
        assert max(abs(x[pc] - sign * y[pc]) for x, y in zip(a, b)) < 2e-5, ('loadings', pc)


with tempfile.TemporaryDirectory(prefix='pcaone-bed-') as name:
    directory = Path(name)
    source, bed, bim, fam = fixture(directory)

    def pcaone(*args, cwd=None):
        return subprocess.run([str(ROOT / 'PCAone'), *map(str, args)], capture_output=True, text=True, cwd=cwd)

    def run(label, bfile=source, seed=42, memory='0.000003', extra=(), verbose=3):
        prefix = directory / label
        result = pcaone('--bfile', bfile, '-m', memory, '-k', 2, '--oversamples', 2, '-n', 2, '-w', 4,
                        '--maxp', 5, '-v', verbose, '-V', '--seed', seed, '-o', prefix, *extra)
        assert result.returncode == 0, result.stdout + result.stderr
        log = Path(str(prefix) + '.log').read_text()
        if verbose < 3 or '-S' in extra:
            for suffix in ['bed', 'bim', 'fam']:
                assert not Path(str(prefix) + '.perm.' + suffix).exists()
            return prefix, [], log
        perm_bim = Path(str(prefix) + '.perm.bim').read_text().splitlines()
        order = [int(line.split()[1][2:]) for line in perm_bim]
        assert sorted(order) == list(range(M))
        assert perm_bim == [bim[i] for i in order]
        expected = bed[:3] + b''.join(bed[3 + i * WIDTH:3 + (i + 1) * WIDTH] for i in order)
        assert Path(str(prefix) + '.perm.bed').read_bytes() == expected
        assert Path(str(prefix) + '.perm.fam').read_bytes() == fam
        mbim = Path(str(prefix) + '.mbim').read_text().splitlines()
        assert [line.split()[:6] for line in mbim] == [line.split()[:6] for line in bim]
        return prefix, order, log

    # --seed decides the order; --buffer does not
    _, order, _ = run('seed')
    _, again, _ = run('again')
    _, changed, _ = run('changed', seed=43)
    _, buffered, _ = run('buffered', extra=('--buffer', 1))
    assert order == again == buffered and order != changed
    # PortableRng: std::shuffle gave another order for the same --seed on libc++ (macOS)
    assert order[:8] == [5, 9, 11, 14, 19, 23, 24, 25]
    assert sum(i * source for i, source in enumerate(order)) == 272509205

    rng = random.Random(7)
    for label, memory, extra in [('many', '0.000003', ()), ('few', '0.0002', ()),
                                 ('emu', '0.000003', ('--emu', '--maxiter', 2))]:
        permuted, order, log = run(label, memory=memory, extra=extra)
        band = int(re.search(r'SNPs per band:\s*(\d+)', log)[1])
        blocksize, nblocks, factor = map(int, re.search(
            r'after adjustment by PCAone: .*blocksize = (\d+) , nblocks = (\d+) , factor = (\d+)', log).groups())
        assert band == blocksize * factor and M % band != 0
        if label == 'few':
            assert factor == 1  # a band is a read block
        else:
            # bands of 26 read blocks, and fewer blocks than 4 bands of them
            assert factor > 1 and nblocks < 4 * factor
        bands = [order[first:first + band] for first in range(0, M, band)]
        assert all(chunk == sorted(chunk) for chunk in bands)
        # Any order within each band, across its read blocks, gives the same PCs
        # as the permuted run: winSVD changes Omg only between bands.
        within = []
        for chunk in bands:
            chunk = chunk[:]
            rng.shuffle(chunk)
            within += chunk
        copy = reordered(directory, label + '-within', within, bed, bim, fam)
        result, _, _ = run(label + '-within-pca', copy, memory=memory, extra=('-S', *extra))
        compare(permuted, result, within)
        # Moving SNPs between bands changes them, so the comparison can tell.
        across = order[:]
        rng.shuffle(across)
        copy = reordered(directory, label + '-across', across, bed, bim, fam)
        result, _, _ = run(label + '-across-pca', copy, memory=memory, extra=('-S', *extra))
        try:
            compare(permuted, result, across)
        except AssertionError:
            pass
        else:
            raise AssertionError('a different band membership gave the same PCs: ' + label)

    run('cleanup', verbose=1)
    run('disabled', extra=('-S',))
    result = pcaone('-b', source, '-m', '0.001', '--bed-shuffle', 'block', '-o', directory / 'removed')
    assert result.returncode != 0

    for label, bad_bed, bad_bim in [('header', b'\x6c', bim), ('body', bed[:-1], bim), ('empty', bed[:3], [])]:
        broken = write_plink(directory / ('bad-' + label), bad_bed, bad_bim, fam)
        out = directory / ('failed-' + label)
        result = pcaone('-b', broken, '-m', '0.000003', '-k', 2, '-n', 2, '-o', out)
        assert result.returncode != 0
        assert not Path(str(out) + '.perm.bed').exists()

    # The output must never overwrite the input, whatever the spelling or link.
    alias = write_plink(directory / 'alias.perm', bed, bim, fam)
    result = pcaone('-b', './alias.perm', '-m', '0.000003', '-k', 2, '-n', 2, '-o', 'alias', cwd=directory)
    assert result.returncode != 0 and 'overwrite the input' in result.stdout + result.stderr
    assert plink_file(alias, '.bed').read_bytes() == bed and plink_file(alias, '.fam').read_bytes() == fam
    os.symlink(plink_file(source, '.bed'), directory / 'link.perm.bed')
    result = pcaone('-b', source, '-m', '0.000003', '-k', 2, '-n', 2, '-o', directory / 'link')
    assert result.returncode != 0 and 'overwrite the input' in result.stdout + result.stderr
    assert (directory / 'link.perm.bed').is_symlink()

    # A failure part way removes what the permutation wrote, and only that.
    (directory / 'blocked.perm.bim').mkdir()
    result = pcaone('-b', source, '-m', '0.000003', '-k', 2, '-n', 2, '-o', directory / 'blocked')
    assert result.returncode != 0 and 'blocked.perm.bim' in result.stdout + result.stderr
    assert not (directory / 'blocked.perm.bed').exists() and (directory / 'blocked.perm.bim').is_dir()
    assert plink_file(source, '.bed').read_bytes() == bed
print('BED permutation: bytes, metadata, seed, band invariance (winSVD/EMU), and guards passed')
