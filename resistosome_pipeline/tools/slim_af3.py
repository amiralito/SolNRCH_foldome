#!/usr/bin/env python3
"""
Slim AF3 tar.gz outputs for the SNI pipeline (runs with numpy only; orjson optional).
For each <name>.tar.gz writes into an output tar:
  <name>/<name>_model.cif.gz
  <name>/<name>_summary_confidences.json
  <name>/<name>_paemask.npz   (token_chain_ids, per protein-chain token ranges, and for every
                               protein chain pair the packed boolean mask min(PAE_ij, PAE_ji^T) < cutoff)
Usage: slim_af3.py OUT_TAR PAE_CUTOFF file1.tar.gz file2.tar.gz ...
"""
import sys, os, io, tarfile, gzip, time
import numpy as np
try:
    import orjson as _j; loads = _j.loads
except ImportError:
    import json as _j; loads = _j.loads

out_tar, cutoff, files = sys.argv[1], float(sys.argv[2]), sys.argv[3:]
mode = 'a' if os.path.exists(out_tar) else 'w'
with tarfile.open(out_tar, mode) as tout:
    done = set(m.split('/')[0] for m in tout.getnames()) if mode == 'a' else set()
    for f in files:
        name = os.path.basename(f)[:-7]
        if name in done:
            continue
        t0 = time.time()
        want = {f'{name}/{name}_model.cif': 'cif', f'{name}/{name}_summary_confidences.json': 'sum',
                f'{name}/{name}_confidences.json': 'conf'}
        got = {}
        with tarfile.open(f, 'r|gz') as tin:
            for m in tin:
                k = want.get(m.name)
                if k:
                    got[k] = tin.extractfile(m).read()
                    if len(got) == 3:
                        break
        conf = loads(got.pop('conf'))
        tci = conf['token_chain_ids']
        pae = np.asarray(conf['pae'], dtype=np.float32); del conf
        ranges, start = {}, 0
        for i in range(1, len(tci) + 1):
            if i == len(tci) or tci[i] != tci[i - 1]:
                ranges[tci[i - 1]] = (start, i); start = i
        prot = [c for c, (a, b) in ranges.items() if b - a >= 30]
        arrays = {'token_chain_ids': np.array(tci), 'prot_chains': np.array(prot),
                  'prot_ranges': np.array([ranges[c] for c in prot]), 'cutoff': np.array(cutoff)}
        for i, a in enumerate(prot):
            for b in prot[i + 1:]:
                ra, rb = ranges[a], ranges[b]
                blk = np.minimum(pae[ra[0]:ra[1], rb[0]:rb[1]], pae[rb[0]:rb[1], ra[0]:ra[1]].T) < cutoff
                arrays[f'mask_{a}_{b}'] = np.packbits(blk, axis=None)
                arrays[f'shape_{a}_{b}'] = np.array(blk.shape)
        bio = io.BytesIO(); np.savez_compressed(bio, **arrays)
        payload = {f'{name}/{name}_model.cif.gz': gzip.compress(got['cif'], 6),
                   f'{name}/{name}_summary_confidences.json': got['sum'],
                   f'{name}/{name}_paemask.npz': bio.getvalue()}
        for p, data in payload.items():
            ti = tarfile.TarInfo(p); ti.size = len(data); tout.addfile(ti, io.BytesIO(data))
        print(f'{name}\t{time.time()-t0:.1f}s', flush=True)
