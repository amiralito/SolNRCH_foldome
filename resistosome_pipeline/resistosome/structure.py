"""
Structure parsing: load PDB/CIF, extract protein chains and coordinates.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
from Bio.PDB import PDBParser, MMCIFParser, Selection
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

logger = logging.getLogger(__name__)

# Standard amino acid 3-letter codes
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEC',  # selenomethionine, selenocysteine
}

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'SEC': 'C',
}


def load_structure(filepath: str, fmt: str = 'cif') -> Optional[Structure]:
    """Load a structure from PDB or CIF file."""
    try:
        if fmt == 'cif':
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        structure = parser.get_structure('model', filepath)
        return structure
    except Exception as e:
        logger.error(f"Failed to load structure from {filepath}: {e}")
        return None


def get_protein_chains(structure: Structure) -> List[Chain]:
    """
    Extract protein chains from the first model.
    Filters to chains that have at least some standard amino acid residues.
    """
    model = structure[0]
    protein_chains = []
    for chain in model.get_chains():
        aa_residues = [r for r in chain.get_residues()
                       if r.get_resname() in STANDARD_AA and r.id[0] == ' ']
        if len(aa_residues) >= 10:  # at least 10 residues to be a protein chain
            protein_chains.append(chain)
    return protein_chains


def get_ca_coords(chain: Chain) -> Tuple[np.ndarray, List[int], List[str]]:
    """
    Get Cα coordinates, residue numbers, and residue names for a chain.
    Returns (coords [N,3], residue_numbers, residue_names_1letter).
    """
    coords = []
    res_nums = []
    res_names = []
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        if 'CA' in residue:
            coords.append(residue['CA'].get_vector().get_array())
            res_nums.append(residue.id[1])
            res_names.append(THREE_TO_ONE.get(residue.get_resname(), 'X'))
    return np.array(coords), res_nums, res_names


def get_all_atom_coords(chain: Chain) -> List[Tuple[str, int, str, np.ndarray]]:
    """
    Get all heavy-atom coordinates for a chain.
    Returns list of (chain_id, res_num, atom_name, coord).
    """
    atoms = []
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        for atom in residue.get_atoms():
            if atom.element == 'H':
                continue
            atoms.append((chain.id, residue.id[1], atom.get_name(),
                          atom.get_vector().get_array()))
    return atoms


def chain_centroid(chain: Chain) -> np.ndarray:
    """Compute the Cα centroid of a chain."""
    coords, _, _ = get_ca_coords(chain)
    return coords.mean(axis=0) if len(coords) > 0 else np.zeros(3)


def get_nterm_coord(chain: Chain) -> Optional[np.ndarray]:
    """Get the coordinate of the N-terminal residue (first Cα)."""
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        if 'CA' in residue:
            return residue['CA'].get_vector().get_array()
    return None


def get_sequence(chain: Chain) -> Tuple[str, List[int]]:
    """Get the 1-letter sequence and residue numbers for a chain."""
    seq = []
    nums = []
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        seq.append(THREE_TO_ONE.get(residue.get_resname(), 'X'))
        nums.append(residue.id[1])
    return ''.join(seq), nums
