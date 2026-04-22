#!/usr/bin/env python3
"""
Test the resistosome pipeline with a synthetic 4-protomer ring structure.
"""
import os
import sys
import math
import json
import tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def create_synthetic_pdb(output_path: str, n_chains: int = 4, n_residues: int = 80):
    """
    Create a synthetic PDB file with n_chains arranged in a ring.
    Each chain has helical segments to test DSSP/helix detection.
    """
    chain_ids = [chr(ord('A') + i) for i in range(n_chains)]
    ring_radius = 30.0  # Angstroms

    # Standard alpha-helix parameters
    helix_rise = 1.5  # Å per residue
    helix_radius = 2.3  # Å
    residues_per_turn = 3.6
    ang_per_res = 2 * math.pi / residues_per_turn

    # Simple amino acid sequence with helical propensity
    aa_seq = 'AEQLKRIA' * (n_residues // 8 + 1)
    aa_seq = aa_seq[:n_residues]

    three_letter = {
        'A': 'ALA', 'E': 'GLU', 'Q': 'GLN', 'L': 'LEU',
        'K': 'LYS', 'R': 'ARG', 'I': 'ILE', 'D': 'ASP',
    }

    lines = []
    atom_num = 1

    for ci, chain_id in enumerate(chain_ids):
        # Position chain around ring
        theta = 2 * math.pi * ci / n_chains
        cx = ring_radius * math.cos(theta)
        cy = ring_radius * math.sin(theta)
        cz = 0.0

        # Helix axis along z, with slight outward tilt
        for ri in range(n_residues):
            aa1 = aa_seq[ri]
            aa3 = three_letter.get(aa1, 'ALA')
            resnum = ri + 1

            # Alpha-helix Cα position
            helix_angle = ri * ang_per_res
            x = cx + helix_radius * math.cos(helix_angle)
            y = cy + helix_radius * math.sin(helix_angle)
            z = cz + ri * helix_rise - n_residues * helix_rise / 2

            # Write CA atom
            lines.append(
                f"ATOM  {atom_num:5d}  CA  {aa3} {chain_id}{resnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 30.00           C  "
            )
            atom_num += 1

            # Write N atom (slightly offset)
            nx = x + 0.5
            ny = y - 0.3
            nz = z - 0.7
            lines.append(
                f"ATOM  {atom_num:5d}  N   {aa3} {chain_id}{resnum:4d}    "
                f"{nx:8.3f}{ny:8.3f}{nz:8.3f}  1.00 30.00           N  "
            )
            atom_num += 1

            # Write C atom
            ccx = x - 0.5
            ccy = y + 0.3
            ccz = z + 0.7
            lines.append(
                f"ATOM  {atom_num:5d}  C   {aa3} {chain_id}{resnum:4d}    "
                f"{ccx:8.3f}{ccy:8.3f}{ccz:8.3f}  1.00 30.00           C  "
            )
            atom_num += 1

            # Write O atom
            ox = x - 1.2
            oy = y + 0.8
            oz = z + 0.5
            lines.append(
                f"ATOM  {atom_num:5d}  O   {aa3} {chain_id}{resnum:4d}    "
                f"{ox:8.3f}{oy:8.3f}{oz:8.3f}  1.00 30.00           O  "
            )
            atom_num += 1

        lines.append(f"TER   {atom_num:5d}      {aa3} {chain_id}{n_residues:4d}")
        atom_num += 1

    lines.append("END")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def test_pipeline():
    """Run the full pipeline test."""
    work = tempfile.mkdtemp(prefix='resisto_test_')
    input_dir = os.path.join(work, 'input')
    output_dir = os.path.join(work, 'output')
    os.makedirs(input_dir)
    os.makedirs(output_dir)

    print("=" * 60)
    print("Testing Resistosome Analysis Pipeline")
    print("=" * 60)

    # Create synthetic PDB files (2 replicates)
    for seed in [1, 2]:
        pdb_path = os.path.join(input_dir, f'test_protein_4protomers_seed{seed}.pdb')
        create_synthetic_pdb(pdb_path, n_chains=4, n_residues=80)
        print(f"Created synthetic PDB: {pdb_path}")

    # Create dummy annotation CSV
    import pandas as pd
    mhd_csv = os.path.join(work, 'mhd.csv')
    ploop_csv = os.path.join(work, 'ploop.csv')
    pd.DataFrame({
        'ID_normalized': ['test_protein'],
        'start': [60],
        'end': [62],
    }).to_csv(mhd_csv, index=False)
    pd.DataFrame({
        'ID_normalized': ['test_protein'],
        'start': [25],
        'end': [33],
    }).to_csv(ploop_csv, index=False)

    # Run pipeline
    from resistosome.pipeline import run_pipeline
    import logging
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s')

    summary_path = run_pipeline(
        input_dir=input_dir,
        output_dir=output_dir,
        mhd_csv=mhd_csv,
        ploop_csv=ploop_csv,
        work_dir=os.path.join(work, 'tmp'),
    )

    # Check outputs
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)

    if summary_path and os.path.exists(summary_path):
        df = pd.read_excel(summary_path)
        print("\nSummary table:")
        print(df.to_string())
    else:
        print("WARNING: No summary table generated")

    # Check per-structure JSONs
    json_dir = os.path.join(output_dir, 'per_structure_json')
    if os.path.exists(json_dir):
        for fname in sorted(os.listdir(json_dir)):
            fpath = os.path.join(json_dir, fname)
            with open(fpath) as f:
                data = json.load(f)
            if 'aggregated' in fname:
                print(f"\nAggregated results ({fname}):")
                for k, v in data.items():
                    if k not in ('base_id', 'seeds'):
                        print(f"  {k}: {v}")
            else:
                print(f"\nPer-structure ({fname}):")
                for k in ['iptm', 'contacts_lcb_within', 'bsa_lcb_within',
                           'sd_theta_rot', 'symmetry_rmsd', 'd_apex_mean',
                           'theta_apex_mean', 'l_apex_mean', 'h_abs_mean',
                           'mu_h_mean', 'd_mhd_p_mean']:
                    print(f"  {k}: {data.get(k)}")

    # Cleanup
    import shutil
    shutil.rmtree(work, ignore_errors=True)
    print("\n✓ Test complete")


if __name__ == '__main__':
    test_pipeline()
