#!/usr/bin/env python3
"""
Run resistosome_pipeline_v2 on 'slim' AF3 inputs (model.cif + summary_confidences +
packed PAE<cutoff masks for every protein-chain pair) produced by slim_af3.py.
Equivalent to the full run for the PAE filter at the cutoff used when slimming.
Usage: run_pipeline_slim.py <pipeline_dir> <run_pipeline.py args...>
"""
import sys, os, json, runpy
import numpy as np
pipe = sys.argv[1]; sys.argv = [os.path.join(pipe, 'run_pipeline.py')] + sys.argv[2:]
sys.path.insert(0, pipe)
from resistosome import af3, pipeline

_orig = af3.load_full_confidences
def load_full_confidences(path):
    with open(path) as fh:
        d = json.load(fh)
    if '_slim_npz' not in d:
        return d
    z = np.load(os.path.join(os.path.dirname(path), d['_slim_npz']))
    tci = [str(x) for x in z['token_chain_ids']]
    n = len(tci)
    pae = np.full((n, n), 99.0, dtype=np.float32)
    rng = {c: tuple(r) for c, r in zip([str(x) for x in z['prot_chains']], z['prot_ranges'])}
    prot = list(rng)
    for i, a in enumerate(prot):
        for b in prot[i + 1:]:
            shp = tuple(z[f'shape_{a}_{b}'])
            m = np.unpackbits(z[f'mask_{a}_{b}'])[: shp[0] * shp[1]].reshape(shp).astype(bool)
            (a0, a1), (b0, b1) = rng[a], rng[b]
            blk = pae[a0:a1, b0:b1]; blk[m] = 0.0
            blkT = pae[b0:b1, a0:a1]; blkT[m.T] = 0.0
    return {'pae': pae, 'token_chain_ids': tci, '_slim_cutoff': float(z['cutoff'])}

af3.load_full_confidences = load_full_confidences
pipeline.load_full_confidences = load_full_confidences
# get_pae_matrix: accept ndarray directly
_gpm = af3.get_pae_matrix
def get_pae_matrix(fc):
    if fc is not None and isinstance(fc.get('pae'), np.ndarray):
        return fc['pae']
    return _gpm(fc)
af3.get_pae_matrix = get_pae_matrix
pipeline.get_pae_matrix = get_pae_matrix
runpy.run_path(sys.argv[0], run_name='__main__')
