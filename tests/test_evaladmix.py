#!/usr/bin/env python3
"""Two-stage evalAdmix: reference reuse, validation, streaming and PGEN parity."""
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / 'PCAone'


def run(*args, error=None, binary=BIN):
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
    if error:
        assert result.returncode != 0 and error in result.stdout + result.stderr, result
    else:
        assert result.returncode == 0, result.stdout + result.stderr
    return result


def values(prefix, suffix):
    return [float(v) for line in Path(str(prefix) + suffix).read_text().splitlines()[1:]
            if line and not line.startswith('#') for v in line.split()]


def agree(a, b):
    for suffix in ['.corres', '.kinship']:
        x, y = values(a, suffix), values(b, suffix)
        assert len(x) == len(y)
        assert all(math.isclose(v, w, abs_tol=2e-6) for v, w in zip(x, y)), suffix


with tempfile.TemporaryDirectory(prefix='pcaone-evaladmix-') as tmp:
    t = Path(tmp)
    source = t / 'data'
    rng = random.Random(21)
    n, m = 32, 401
    bed = bytearray(b'\x6c\x1b\x01')
    for j in range(m):
        calls = []
        for i in range(n):
            p = .15 + .35 * (i < n // 2)
            g = int(rng.random() < p) + int(rng.random() < p)
            calls.append(1 if rng.random() < .04 else {0: 3, 1: 2, 2: 0}[g])
        for i in range(0, n, 4):
            bed.append(sum(calls[i + k] << (2 * k) for k in range(4)))
    source.with_suffix('.bed').write_bytes(bed)
    source.with_suffix('.bim').write_text(''.join(f'1 rs{j} 0 {j+1} A C\n' for j in range(m)))
    source.with_suffix('.fam').write_text(''.join(f'F I{i} 0 0 0 -9\n' for i in range(n)))
    ref = t / 'pcs'
    run('-b', source, '-k', 3, '-d', 0, '-n', 1, '-o', ref)
    reference = ref.with_suffix('.eigvecs').read_bytes()
    common = ['-b', source, '--evaladmix', '-n', 1]
    out = t / 'analysis'
    run(*common, '-P', ref, '-o', out)  # all 3 PCs without -k
    assert sorted(p.suffix for p in t.glob('analysis.*')) == ['.corres', '.kinship', '.log']
    assert reference == ref.with_suffix('.eigvecs').read_bytes()
    scores_only = t / 'scores_only'
    scores_only.with_suffix('.eigvecs').write_bytes(reference)
    run(*common, '-P', scores_only, '-o', scores_only)
    agree(out, scores_only)
    assert scores_only.with_suffix('.eigvecs').read_bytes() == reference
    for label, extra in [('direct', ['--read-U', ref.with_suffix('.eigvecs')]),
                         ('stream', ['-P', ref, '-m', .00001]),
                         ('full_stream', ['-P', ref, '-m', .00001, '-d', 3])]:
        run(*common, *extra, '-o', t / label)
        agree(out, t / label)
    # A leading subset must equal a reference containing only those columns.
    subset = t / 'subset.eigvecs'
    subset.write_text('\n'.join(line.split()[0] for line in reference.decode().splitlines()) + '\n')
    run(*common, '-P', ref, '-k', 1, '-o', t / 'subset_a')
    run(*common, '--read-U', subset, '-o', t / 'subset_b')
    agree(t / 'subset_a', t / 'subset_b')
    run(*common, '-o', out, error='please use -P/--USV')
    run(*common, '-P', ref, '-k', 4, '-o', out, error='larger than the 3 PCs')
    subset.write_text('\n'.join(subset.read_text().splitlines()[:-1]) + '\n')
    run(*common, '--read-U', subset, '-o', out, error='31 rows')
    run(*common, '-P', t / 'absent', '-o', out, error='cannot read PC scores')
    # Optional numerical regression against the former combined PCA + analysis.
    if os.environ.get('PCAONE_BASELINE'):
        for label, extra in [('baseline', []), ('baseline_maf', ['--maf', .05])]:
            run('-b', source, '-k', 3, '-d', 0, '--evaladmix', '-n', 1,
                '-o', t / label, *extra, binary=os.environ['PCAONE_BASELINE'])
            run(*common, '-P', t / label, '-o', t / (label + '_new'), *extra)
            agree(t / label, t / (label + '_new'))
    plink2 = shutil.which('plink2')
    if plink2:
        pg = t / 'pgen'
        run('--bfile', source, '--make-pgen', '--out', pg, '--silent', binary=plink2)
        for label, extra in [('pgen_ic', []), ('pgen_ooc', ['-m', .00001])]:
            run('-p', pg, '--evaladmix', '-P', ref, '-n', 1, '-o', t / label, *extra)
            agree(out, t / label)
    else:
        print('SKIP PGEN parity: plink2 unavailable')
print('PASS two-stage evalAdmix')
