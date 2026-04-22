"""
AlphaFold 3 output parsing: ipTM, PAE matrix, confidence data.
"""

import json
import logging
import numpy as np
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


def load_summary_confidences(filepath: str) -> Optional[Dict[str, Any]]:
    """Load AF3 summary_confidences.json."""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"Could not load summary confidences from {filepath}: {e}")
        return None


def get_iptm(summary_conf: Optional[Dict]) -> Optional[float]:
    """Extract ipTM from summary confidences."""
    if summary_conf is None:
        return None
    # AF3 uses 'iptm'; some versions may use 'interface_ptm'
    for key in ('iptm', 'interface_ptm', 'ipTM'):
        iptm = summary_conf.get(key)
        if iptm is not None:
            return float(iptm)
    return None


def load_full_confidences(filepath: str) -> Optional[Dict[str, Any]]:
    """
    Load AF3 full confidences JSON.
    This contains the PAE matrix among other confidence metrics.
    """
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"Could not load full confidences from {filepath}: {e}")
        return None


def get_pae_matrix(full_conf: Optional[Dict]) -> Optional[np.ndarray]:
    """
    Extract the PAE matrix from full confidences.
    Returns a 2D numpy array of shape (N_residues, N_residues).
    """
    if full_conf is None:
        return None
    pae = full_conf.get('pae')
    if pae is None:
        return None
    try:
        return np.array(pae, dtype=np.float32)
    except Exception as e:
        logger.warning(f"Could not parse PAE matrix: {e}")
        return None


def get_chain_residue_ranges(full_conf: Optional[Dict], n_chains: int,
                              chain_lengths: list) -> Optional[Dict[str, tuple]]:
    """
    Compute residue index ranges for each chain in the PAE matrix.
    If token_chain_ids is available, use it for accurate mapping (handles ligands).
    Otherwise, assume concatenated protein chains in order.
    Returns {chain_id: (start_idx, end_idx)} for indexing into the PAE matrix.
    """
    if full_conf is None:
        return None

    # Try to use token_chain_ids for accurate mapping
    token_chain_ids = full_conf.get('token_chain_ids')
    if token_chain_ids is not None:
        ranges = {}
        current_chain = None
        start_idx = 0
        for i, cid in enumerate(token_chain_ids):
            if cid != current_chain:
                if current_chain is not None:
                    ranges[current_chain] = (start_idx, i)
                current_chain = cid
                start_idx = i
        if current_chain is not None:
            ranges[current_chain] = (start_idx, len(token_chain_ids))
        return ranges

    # Fallback: assume chains are concatenated in alphabetical order
    ranges = {}
    chain_ids = [chr(ord('A') + i) for i in range(n_chains)]
    offset = 0
    for i, cid in enumerate(chain_ids):
        length = chain_lengths[i]
        ranges[cid] = (offset, offset + length)
        offset += length
    return ranges


def get_pae_for_chain_pair(pae_matrix: np.ndarray,
                            chain_ranges: Dict[str, tuple],
                            chain1_id: str, chain2_id: str) -> Optional[np.ndarray]:
    """
    Extract the PAE submatrix for a pair of chains.
    Returns the inter-chain PAE block (residues of chain1 vs chain2).
    """
    if pae_matrix is None or chain_ranges is None:
        return None
    if chain1_id not in chain_ranges or chain2_id not in chain_ranges:
        return None
    r1 = chain_ranges[chain1_id]
    r2 = chain_ranges[chain2_id]
    # PAE[i,j] = predicted error in j given i's position
    # We want the block for chain1 rows, chain2 columns AND vice versa
    block_12 = pae_matrix[r1[0]:r1[1], r2[0]:r2[1]]
    block_21 = pae_matrix[r2[0]:r2[1], r1[0]:r1[1]]
    # Return the element-wise minimum (most confident direction)
    return np.minimum(block_12, block_21.T)
