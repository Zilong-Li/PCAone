#!/usr/bin/env python3
"""--evaladmix-kin: the stripes agree with the dense matrix, related pairs and
the unrelated set, BED and PGEN sample IDs, and the option guards."""
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / 'PCAone'


def run(*args, error=None):
    result = subprocess.run([str(BIN), *map(str, args)], capture_output=True, text=True)
    if error:
        assert result.returncode != 0 and error in result.stdout + result.stderr, result
    else:
        assert result.returncode == 0, result.stdout + result.stderr
    return result


def dense(prefix):
    lines = Path(str(prefix) + '.kinship').read_text().splitlines()
    return lines[0].split('\t'), [row.split('\t') for row in lines[1:]]


def pairs(prefix):
    lines = Path(str(prefix) + '.kin0').read_text().splitlines()
    head = lines[0].lstrip('#').split('\t')
    rows = [dict(zip(head, line.split('\t'))) for line in lines[1:]]
    return lines[0], rows


with tempfile.TemporaryDirectory(prefix='pcaone-evaladmix-pairs-') as tmp:
    t = Path(tmp)
    rng = random.Random(7)
    n, m = 150, 1500
    # two populations; the last 4 samples duplicate samples 0..3 and samples
    # 140..143 are children of 10..17. A few sites miss 30% of the calls, so both
    # ways of counting the sites a pair misses together are exercised.
    freq = [(rng.uniform(.1, .9), rng.uniform(.1, .9)) for _ in range(m)]
    bed = bytearray(b'\x6c\x1b\x01')
    for j in range(m):
        hap = [[int(rng.random() < freq[j][i < n // 2]) for _ in range(2)] for i in range(n)]
        for c in range(4):
            hap[140 + c] = [rng.choice(hap[10 + 2 * c]), rng.choice(hap[11 + 2 * c])]
        for c in range(4):
            hap[n - 1 - c] = list(hap[c])
        rate = .3 if j % 50 == 0 else .02
        calls = [1 if rng.random() < rate else {0: 3, 1: 2, 2: 0}[sum(h)] for h in hap]
        for i in range(0, n, 4):
            bed.append(sum(calls[i + k] << (2 * k) for k in range(4) if i + k < n))
    src = t / 'data'
    src.with_suffix('.bed').write_bytes(bed)
    src.with_suffix('.bim').write_text(''.join(f'1 rs{j} 0 {j+1} A C\n' for j in range(m)))
    src.with_suffix('.fam').write_text(''.join(f'F{i} I{i} 0 0 0 -9\n' for i in range(n)))
    ref = t / 'pcs'
    run('-b', src, '-k', 1, '-n', 2, '-o', ref)
    common = ['-b', src, '-P', ref, '--evaladmix', '-n', 2]
    run(*common, '-o', t / 'dense')
    ids, K = dense(t / 'dense')
    # every pair, in-core (one stripe) and from the file in stripes of 64 samples
    for label, extra in [('incore', []), ('stream', ['-m', .0001])]:
        out = run(*common, '--evaladmix-kin', -.5, *extra, '-o', t / label)
        if label == 'stream':
            assert '3 stripe(s)' in out.stdout, out.stdout
        head, rows = pairs(t / label)
        assert head == '#FID1\tIID1\tFID2\tIID2\tNSNP\tKINSHIP', head
        assert len(rows) == n * (n - 1) // 2
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                r = rows[k]
                k += 1
                assert (r['IID1'], r['IID2'], r['FID1']) == (ids[i], ids[j], 'F' + ids[i][1:]), r
                assert abs(float(r['KINSHIP']) - float(K[i][j])) <= 2e-6, (label, i, j, r, K[i][j])
                assert 0 < int(r['NSNP']) <= m
    # related pairs and the unrelated set
    run(*common, '--evaladmix-kin', .1, '-m', .0001, '-o', t / 'rel')
    _, rows = pairs(t / 'rel')
    got = {(r['IID1'], r['IID2']) for r in rows}
    planted = {(f'I{c}', f'I{n - 1 - c}') for c in range(4)}
    planted |= {(f'I{p}', f'I{140 + c}') for c in range(4) for p in (10 + 2 * c, 11 + 2 * c)}
    assert got == planted, got ^ planted
    keep = [line.split('\t') for line in (t / 'rel.unrelated').read_text().splitlines()]
    kept = {iid for _, iid in keep}
    assert all(fid == 'F' + iid[1:] for fid, iid in keep)
    assert not any(a in kept and b in kept for a, b in got)  # no related pair is kept
    dropped = set(ids) - kept
    assert len(dropped) == 8, dropped  # one per duplicate, the child of each trio
    assert all(any((d, o) in got or (o, d) in got for o in kept) for d in dropped)  # maximal
    # a stricter cutoff for the unrelated set keeps the trios
    run(*common, '--evaladmix-kin', .1, '--evaladmix-unrelated', .4, '-o', t / 'dup')
    assert len((t / 'dup.unrelated').read_text().splitlines()) == n - 4
    # PGEN (a PLINK 1 .bed is a valid .pgen), with and without an FID column
    for label, psam in [('pfid', '#FID\tIID\tSEX\n' + ''.join(f'F{i}\tI{i}\tNA\n' for i in range(n))),
                        ('piid', '#IID\tSEX\n' + ''.join(f'I{i}\tNA\n' for i in range(n)))]:
        pg = t / label
        shutil.copy(src.with_suffix('.bed'), pg.with_suffix('.pgen'))
        pg.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' +
                                           ''.join(f'1\t{j+1}\trs{j}\tC\tA\n' for j in range(m)))
        pg.with_suffix('.psam').write_text(psam)
        run('-p', pg, '-P', ref, '--evaladmix', '-n', 2, '--evaladmix-kin', .1, '-m', .0001, '-o', t / (label + '_o'))
    assert (t / 'pfid_o.kin0').read_text() == (t / 'rel.kin0').read_text()
    assert (t / 'pfid_o.unrelated').read_text() == (t / 'rel.unrelated').read_text()
    piid = (t / 'piid_o.kin0').read_text().splitlines()
    assert piid[0] == '#IID1\tIID2\tNSNP\tKINSHIP'
    assert piid[1:] == ['\t'.join(f[1::2][:2] + f[4:]) for f in
                        (line.split('\t') for line in (t / 'rel.kin0').read_text().splitlines()[1:])]
    assert (t / 'piid_o.unrelated').read_text() == ''.join(line.split('\t')[1] + '\n' for line in
                                                          (t / 'rel.unrelated').read_text().splitlines())
    # the log compares the scatter of unrelated pairs with chance; scores that
    # miss the two populations scatter them far more
    log = (t / 'rel.log').read_text()
    assert 'chance alone gives about' in log and 'times as much as chance' not in log, log
    noise = t / 'noise.eigvecs'
    noise.write_text(''.join(f'{rng.gauss(0, 1):.6f}\n' for _ in range(n)))
    out = run('-b', src, '--read-U', noise, '--evaladmix', '--evaladmix-kin', .1, '-o', t / 'noise')
    assert 'times as much as chance alone' in out.stdout + out.stderr, out.stdout
    # the scores must be of the same samples in the same order: .eigvecs2 has their IDs
    rev = t / 'rev'
    rev.with_suffix('.eigvecs').write_bytes(ref.with_suffix('.eigvecs').read_bytes())
    lines = ref.with_suffix('.eigvecs2').read_text().splitlines()
    (t / 'rev.eigvecs2').write_text('\n'.join(lines[:1] + lines[:0:-1]) + '\n')
    run('-b', src, '-P', rev, '--evaladmix', '--evaladmix-kin', .1, '-o', t / 'x', error='same samples in the same order')
    run('-b', src, '-P', rev, '--evaladmix', '-o', t / 'x', error='same samples in the same order')
    (t / 'rev.eigvecs2').unlink()  # without it, only the count is checked
    assert 'only their number is checked' in run('-b', src, '-P', rev, '--evaladmix', '-o', t / 'x').stdout
    # guards
    run('-b', src, '-P', ref, '--evaladmix-kin', .1, '-o', t / 'x', error='requires --evaladmix')
    run(*common, '--evaladmix-unrelated', .1, '-o', t / 'x', error='requires --evaladmix-kin')
    run(*common, '--evaladmix-kin', .1, '--evaladmix-unrelated', .05, '-o', t / 'x', error='[--evaladmix-kin, 0.5]')
    run(*common, '--evaladmix-kin', .7, '-o', t / 'x', error='in [-0.5, 0.5]')
print('PASS evalAdmix pairs (--evaladmix-kin)')
