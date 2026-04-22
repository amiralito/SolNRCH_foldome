#!/usr/bin/env python3
"""
Determine the spatial ring order of protomer chains in a resistosome structure.

For each protein chain, computes the Cα centroid, then sorts chains by their
angular position around the ring centre. This gives the true adjacency order
regardless of chain naming conventions.

Usage:
    python get_ring_order.py structure.cif [--format cif|pdb]

Output:
    Prints the chain IDs in ring order (comma-separated) and a mapping file.
"""

import argparse
import numpy as np
import gemmi
import sys


def get_protein_chain_centroids(structure, use_auth_id=True):
    """
    Extract Cα centroids for each protein chain.
    Returns dict of chain_id -> centroid (np.array of shape (3,)).
    """
    model = structure[0]  # first model
    centroids = {}

    for chain in model:
        chain_id = chain.name  # this is auth_asym_id by default in gemmi
        ca_coords = []

        for residue in chain:
            # Skip non-protein residues (ligands, water, etc.)
            if not residue.find_atom("CA", "\0"):
                continue
            # Check it's a standard amino acid
            ri = gemmi.find_tabulated_residue(residue.name)
            if not ri.is_amino_acid():
                continue

            ca = residue.find_atom("CA", "\0")
            if ca:
                pos = ca.pos
                ca_coords.append([pos.x, pos.y, pos.z])

        if len(ca_coords) > 20:  # must have substantial protein content
            centroids[chain_id] = np.mean(ca_coords, axis=0)

    return centroids


def determine_ring_order(centroids):
    """
    Given chain centroids, determine their angular order around the ring centre.
    Returns list of chain IDs in ring-adjacent order.
    """
    chain_ids = list(centroids.keys())
    coords = np.array([centroids[c] for c in chain_ids])

    # Ring centre = mean of all centroids
    centre = np.mean(coords, axis=0)

    # Fit plane through centroids using PCA to get the ring normal
    centred = coords - centre
    _, _, Vt = np.linalg.svd(centred)
    normal = Vt[-1]  # smallest component = plane normal

    # Project centroids onto the plane
    # Use the two principal axes as x, y in the plane
    axis_x = Vt[0]
    axis_y = Vt[1]

    angles = []
    for i, chain_id in enumerate(chain_ids):
        vec = centred[i]
        x_proj = np.dot(vec, axis_x)
        y_proj = np.dot(vec, axis_y)
        angle = np.arctan2(y_proj, x_proj)
        angles.append((angle, chain_id))

    # Sort by angle to get ring order
    angles.sort(key=lambda x: x[0])
    ring_order = [chain_id for _, chain_id in angles]

    return ring_order, centre, normal


def print_adjacency_pairs(ring_order):
    """Print adjacent pairs in the ring."""
    n = len(ring_order)
    pairs = []
    for i in range(n):
        c1 = ring_order[i]
        c2 = ring_order[(i + 1) % n]
        pairs.append(f"{c1}--{c2}")
    return pairs


def main():
    parser = argparse.ArgumentParser(
        description="Determine spatial ring order of protomer chains"
    )
    parser.add_argument("structure", help="Path to CIF or PDB file")
    parser.add_argument(
        "--format", choices=["cif", "pdb"], default=None,
        help="File format (auto-detected from extension if not specified)"
    )
    args = parser.parse_args()

    # Auto-detect format
    fmt = args.format
    if fmt is None:
        if args.structure.endswith(".cif") or args.structure.endswith(".mmcif"):
            fmt = "cif"
        else:
            fmt = "pdb"

    # Read structure
    if fmt == "cif":
        doc = gemmi.cif.read(args.structure)
        structure = gemmi.make_structure_from_block(doc.sole_block())
    else:
        structure = gemmi.read_pdb(args.structure)

    structure.setup_entities()

    # Get centroids
    centroids = get_protein_chain_centroids(structure)

    if len(centroids) < 3:
        print(f"ERROR: Only {len(centroids)} protein chains found. Need at least 3 for ring detection.", file=sys.stderr)
        sys.exit(1)

    print(f"Structure: {args.structure}")
    print(f"Protein chains found: {list(centroids.keys())}")
    print(f"N protomers: {len(centroids)}")
    print()

    # Determine ring order
    ring_order, centre, normal = determine_ring_order(centroids)

    # Compute inter-centroid distances to validate adjacency
    print(f"Ring order: {','.join(ring_order)}")
    print(f"Adjacent pairs: {' | '.join(print_adjacency_pairs(ring_order))}")
    print()

    # Print distances between adjacent pairs for validation
    print("Adjacent pair distances (Å):")
    n = len(ring_order)
    for i in range(n):
        c1 = ring_order[i]
        c2 = ring_order[(i + 1) % n]
        dist = np.linalg.norm(centroids[c1] - centroids[c2])
        print(f"  {c1}--{c2}: {dist:.1f}")

    # Also show what alphabetical order would give
    alpha_order = sorted(centroids.keys())
    print(f"\nAlphabetical order: {','.join(alpha_order)}")
    if alpha_order == ring_order:
        print("✓ Alphabetical order matches ring order")
    else:
        print("✗ Alphabetical order does NOT match ring order — use ring order for adjacency!")


if __name__ == "__main__":
    main()
