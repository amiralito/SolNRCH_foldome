"""
Geometric measurements: interface angles, symmetry RMSD, N-terminal distances.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
from Bio.PDB.Chain import Chain

from .structure import chain_centroid, get_nterm_coord, STANDARD_AA

logger = logging.getLogger(__name__)


def _get_ca_at_resnum(chain: Chain, resnum: int) -> Optional[np.ndarray]:
    """Get Cα coordinate at a specific residue number."""
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        if residue.id[1] == resnum and 'CA' in residue:
            return residue['CA'].get_vector().get_array()
    return None


def compute_ring_centroid(chains: List[Chain]) -> np.ndarray:
    """Compute the overall centroid of the ring from chain centroids."""
    centroids = np.array([chain_centroid(c) for c in chains])
    return centroids.mean(axis=0)


def compute_interface_angles(chains: List[Chain]) -> List[float]:
    """
    Compute the angle subtended at the ring centroid by each adjacent pair.
    For a perfect n-mer ring, each angle should be 360/n degrees.
    Returns angles for pairs: (0,1), (1,2), ..., (n-1,0).
    """
    n = len(chains)
    centroids = np.array([chain_centroid(c) for c in chains])
    ring_center = centroids.mean(axis=0)

    angles = []
    for i in range(n):
        j = (i + 1) % n
        v1 = centroids[i] - ring_center
        v2 = centroids[j] - ring_center

        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_deg = np.degrees(np.arccos(cos_angle))
        angles.append(round(angle_deg, 2))

    return angles


def compute_symmetry_rmsd(chains: List[Chain]) -> float:
    """
    Compute symmetry RMSD: RMSD of chain-centroid distances from their mean.
    Low values indicate high rotational symmetry.
    """
    centroids = np.array([chain_centroid(c) for c in chains])
    ring_center = centroids.mean(axis=0)

    distances = np.linalg.norm(centroids - ring_center, axis=1)
    mean_dist = distances.mean()

    rmsd = np.sqrt(np.mean((distances - mean_dist) ** 2))
    return round(rmsd, 3)


def compute_nterm_distances_allpairs(chains: List[Chain],
                                      nterm_resnum: Optional[int] = None
                                      ) -> Tuple[Dict[str, float], float]:
    """
    Compute all pairwise N-terminal distances between chains.

    Parameters
    ----------
    chains : list of Bio.PDB.Chain
    nterm_resnum : optional residue number to use as the N-terminal reference
                   (e.g. from MADA HMM hit). If None, uses the first residue.

    Returns (distance_dict, mean_distance).
    """
    n = len(chains)
    coords = {}
    for chain in chains:
        if nterm_resnum is not None:
            coord = _get_ca_at_resnum(chain, nterm_resnum)
        else:
            coord = get_nterm_coord(chain)
        if coord is not None:
            coords[chain.id] = coord

    distances = {}
    dist_values = []
    for i in range(n):
        for j in range(i + 1, n):
            ci = chains[i].id
            cj = chains[j].id
            if ci in coords and cj in coords:
                d = float(np.linalg.norm(coords[ci] - coords[cj]))
                distances[f"{ci}-{cj}"] = round(d, 2)
                dist_values.append(d)

    mean_dist = float(np.mean(dist_values)) if dist_values else 0.0
    return distances, round(mean_dist, 2)


def compute_adjacent_nterm_distances(chains: List[Chain]) -> List[float]:
    """
    Compute N-terminal distances for adjacent chain pairs only (circular).
    Returns list of distances for pairs (0,1), (1,2), ..., (n-1,0).
    """
    n = len(chains)
    distances = []
    for i in range(n):
        j = (i + 1) % n
        c1 = get_nterm_coord(chains[i])
        c2 = get_nterm_coord(chains[j])
        if c1 is not None and c2 is not None:
            distances.append(float(np.linalg.norm(c1 - c2)))
        else:
            distances.append(np.nan)
    return distances


def compute_pore_axis(chains: List[Chain]) -> Optional[np.ndarray]:
    """
    Compute the central pore/symmetry axis of the resistosome ring.

    Uses PCA on the protomer centroids: the pore axis is the normal to the
    best-fit plane of the centroids (i.e. the smallest principal component).
    For a well-formed ring, the centroids lie roughly in a plane and the
    axis perpendicular to that plane is the pore axis.

    Returns a unit vector, or None if <3 chains.
    """
    centroids = np.array([chain_centroid(c) for c in chains])
    if len(centroids) < 3:
        return None

    centered = centroids - centroids.mean(axis=0)
    _, _, Vt = np.linalg.svd(centered)
    # Last component = normal to best-fit plane = pore axis
    pore_axis = Vt[-1]
    return pore_axis / np.linalg.norm(pore_axis)


def determine_ring_order(chains: List[Chain]) -> List[Chain]:
    """
    Determine the spatial ring order of protomer chains by angular position.

    Computes Cα centroids, fits a plane via PCA, then sorts chains by their
    angle around the ring centre. This gives the true adjacency order
    regardless of chain naming conventions.

    For AlphaFold 3 outputs (chains A-F), this matches alphabetical order.
    For PDB/RCSB structures with arbitrary chain letters, this correctly
    identifies which chains are actually adjacent in the ring.

    Falls back to input order if < 3 chains.
    """
    if len(chains) < 3:
        return chains

    # Compute Cα centroid for each chain
    centroids = np.array([chain_centroid(c) for c in chains])

    # Ring centre
    centre = centroids.mean(axis=0)
    centred = centroids - centre

    # PCA to get the plane of the ring
    _, _, Vt = np.linalg.svd(centred)
    axis_x = Vt[0]  # first principal component (in-plane)
    axis_y = Vt[1]  # second principal component (in-plane)

    # Compute angular position of each chain around the ring
    angles = []
    for i, chain in enumerate(chains):
        vec = centred[i]
        x_proj = np.dot(vec, axis_x)
        y_proj = np.dot(vec, axis_y)
        angle = np.arctan2(y_proj, x_proj)
        angles.append((angle, chain))

    # Sort by angle → ring-adjacent order
    angles.sort(key=lambda x: x[0])
    ordered = [chain for _, chain in angles]

    # Log if order differs from input
    input_ids = [c.id for c in chains]
    output_ids = [c.id for c in ordered]
    if input_ids != output_ids:
        logger.info(f"  Ring order differs from input: {input_ids} → {output_ids}")
    else:
        logger.debug(f"  Ring order matches input: {output_ids}")

    return ordered


def validate_ring_adjacency(chains: List[Chain]) -> List[float]:
    """
    Compute centroid-to-centroid distances for adjacent pairs in current order.
    Useful for QC: adjacent distances should all be similar.
    """
    centroids = np.array([chain_centroid(c) for c in chains])
    n = len(chains)
    distances = []
    for i in range(n):
        j = (i + 1) % n
        d = float(np.linalg.norm(centroids[i] - centroids[j]))
        distances.append(round(d, 1))
    return distances
