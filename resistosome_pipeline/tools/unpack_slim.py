#!/usr/bin/env python3
"""Unpack *.slim.tar (or bundles of them) into AF3-like directories for run_pipeline_slim.py."""
import sys, os, tarfile, gzip, io, json, glob
src, dst = sys.argv[1], sys.argv[2]
os.makedirs(dst, exist_ok=True)
def unpack_one(tf):
    for m in tf.getmembers():
        if not m.isfile(): continue
        data = tf.extractfile(m).read()
        name, fn = m.name.split('/')[-2], m.name.split('/')[-1]
        d = os.path.join(dst, name); os.makedirs(d, exist_ok=True)
        if fn.endswith('_model.cif.gz'):
            open(os.path.join(d, fn[:-3]), 'wb').write(gzip.decompress(data))
        else:
            open(os.path.join(d, fn), 'wb').write(data)
            if fn.endswith('_paemask.npz'):
                json.dump({'_slim_npz': fn}, open(os.path.join(d, fn.replace('_paemask.npz', '_confidences.json')), 'w'))
n = 0
for f in sorted(glob.glob(os.path.join(src, '*.tar'))):
    with tarfile.open(f) as t:
        inner = [m for m in t.getmembers() if m.name.endswith('.slim.tar')]
        if inner:
            for m in inner:
                with tarfile.open(fileobj=io.BytesIO(t.extractfile(m).read())) as t2: unpack_one(t2); n += 1
        else:
            unpack_one(t); n += 1
print('unpacked', n)
