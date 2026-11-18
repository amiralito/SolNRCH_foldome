#!/usr/bin/env python3
"""Aggregate v2 per-structure JSONs (possibly from several partial runs) into one
per-NLR summary, using the v2 pipeline's own aggregate_replicates/build_summary_table.
Usage: aggregate_v2_json.py --pipeline_dir P --out OUT json_dir1 [json_dir2 ...]"""
import argparse, glob, json, os, sys, shutil
ap = argparse.ArgumentParser(); ap.add_argument('--pipeline_dir', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--k', type=float, default=1.96)
ap.add_argument('dirs', nargs='+'); a = ap.parse_args()
sys.path.insert(0, a.pipeline_dir)
from resistosome.pipeline import aggregate_replicates, build_summary_table, provenance
jd = os.path.join(a.out, 'per_structure_json'); os.makedirs(jd, exist_ok=True)
groups = {}
for d in a.dirs:
    for f in glob.glob(os.path.join(d, '*_results.json')):
        r = json.load(open(f)); groups.setdefault(r['base_id'], {})[r.get('seed')] = r
        shutil.copy(f, jd)
aggs = []
for bid in sorted(groups):
    reps = [groups[bid][s] for s in sorted(groups[bid], key=lambda x: x or 0)]
    ag = aggregate_replicates(reps, k=a.k); ag.update({'pipeline_version': reps[0].get('pipeline_version'),
                                                       'git_commit': reps[0].get('git_commit')})
    json.dump(ag, open(os.path.join(jd, f'{bid}_aggregated.json'), 'w'), indent=2, default=str)
    aggs.append(ag)
print(sum(len(v) for v in groups.values()), 'structures,', len(groups), 'NLRs ->', build_summary_table(aggs, a.out))
