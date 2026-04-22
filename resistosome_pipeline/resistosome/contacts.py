"""
Inter-chain contact counting and BSA (buried surface area) calculations.
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import freesasa

from .config import CONTACT_CUTOFF_A, PAE_CUTOFF_A, SASA_PROBE_RADIUS
from .structure import STANDARD_AA

logger = logging.getLogger(__name__)


def count_contacts(chain1: Chain, chain2: Chain,
                   cutoff: float = CONTACT_CUTOFF_A,
                   pae_block: Optional[np.ndarray] = None,
                   pae_cutoff: float = PAE_CUTOFF_A,
                   chain1_resnums: Optional[List[int]] = None,
                   chain2_resnums: Optional[List[int]] = None) -> int:
    """
    Count inter-chain atomic contacts within cutoff distance.
    If pae_block is provided, only count contacts where PAE < pae_cutoff.
    Uses KDTree for efficient spatial search.
    """
    from scipy.spatial import cKDTree

    # Build residue-to-index mappings for PAE lookup
    if pae_block is not None and chain1_resnums is not None and chain2_resnums is not None:
        c1_res_to_idx = {rn: i for i, rn in enumerate(chain1_resnums)}
        c2_res_to_idx = {rn: i for i, rn in enumerate(chain2_resnums)}
        use_pae = True
    else:
        use_pae = False

    # Collect heavy atoms for each chain with residue number tracking
    def get_atoms_with_resnum(chain):
        coords = []
        resnums = []
        for residue in chain.get_residues():
            if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
                continue
            rnum = residue.id[1]
            for atom in residue.get_atoms():
                if atom.element != 'H':
                    coords.append(atom.get_vector().get_array())
                    resnums.append(rnum)
        return np.array(coords) if coords else np.empty((0, 3)), resnums

    coords1, resnums1 = get_atoms_with_resnum(chain1)
    coords2, resnums2 = get_atoms_with_resnum(chain2)

    if len(coords1) == 0 or len(coords2) == 0:
        return 0

    # KDTree for fast neighbor search
    tree2 = cKDTree(coords2)
    pairs = tree2.query_ball_point(coords1, cutoff)

    contact_count = 0
    if not use_pae:
        # Simply count all contacts
        contact_count = sum(len(p) for p in pairs)
    else:
        # Filter by PAE
        for i, neighbors in enumerate(pairs):
            rn1 = resnums1[i]
            idx1 = c1_res_to_idx.get(rn1)
            if idx1 is None or idx1 >= pae_block.shape[0]:
                contact_count += len(neighbors)  # no PAE info → count anyway
                continue
            for j in neighbors:
                rn2 = resnums2[j]
                idx2 = c2_res_to_idx.get(rn2)
                if idx2 is None or idx2 >= pae_block.shape[1]:
                    contact_count += 1
                elif pae_block[idx1, idx2] < pae_cutoff:
                    contact_count += 1

    return contact_count


def extract_contact_details(chain1: Chain, chain2: Chain,
                             cutoff: float = CONTACT_CUTOFF_A,
                             pae_block: Optional[np.ndarray] = None,
                             pae_cutoff: float = PAE_CUTOFF_A,
                             chain1_resnums: Optional[List[int]] = None,
                             chain2_resnums: Optional[List[int]] = None) -> List[Dict]:
    """
    Extract detailed contact information between two chains with
    interaction type classification.

    Returns a list of dicts, each describing one atomic contact:
      {chain1, resnum1, resname1, atom1, chain2, resnum2, resname2, atom2,
       distance, contact_type}

    Contact types:
      - salt_bridge:  charged pair (Arg/Lys/His + ↔ Asp/Glu −), < 4.0 Å
      - hydrogen_bond: N/O donor ↔ N/O/S acceptor, < 3.5 Å (heavy-atom distance)
      - disulfide:    Cys SG ↔ Cys SG, < 2.5 Å
      - hydrophobic:  C ↔ C between apolar atoms/residues, < 4.0 Å
      - vdw:          any other contact within cutoff
    """
    from scipy.spatial import cKDTree
    from .structure import THREE_TO_ONE

    if pae_block is not None and chain1_resnums is not None and chain2_resnums is not None:
        c1_res_to_idx = {rn: i for i, rn in enumerate(chain1_resnums)}
        c2_res_to_idx = {rn: i for i, rn in enumerate(chain2_resnums)}
        use_pae = True
    else:
        use_pae = False

    # Collect atoms with full annotation (including element and 3-letter code)
    def get_annotated_atoms(chain):
        entries = []
        for residue in chain.get_residues():
            rname3 = residue.get_resname()
            if rname3 not in STANDARD_AA or residue.id[0] != ' ':
                continue
            rn = residue.id[1]
            rname1 = THREE_TO_ONE.get(rname3, 'X')
            for atom in residue.get_atoms():
                if atom.element != 'H':
                    entries.append({
                        'coord': atom.get_vector().get_array(),
                        'resnum': rn,
                        'resname': rname1,
                        'resname3': rname3,
                        'atom': atom.get_name(),
                        'element': atom.element.upper(),
                    })
        return entries

    atoms1 = get_annotated_atoms(chain1)
    atoms2 = get_annotated_atoms(chain2)

    if not atoms1 or not atoms2:
        return []

    coords1 = np.array([a['coord'] for a in atoms1])
    coords2 = np.array([a['coord'] for a in atoms2])

    # Use max of cutoff and 4.0 Å to catch salt bridges and hydrophobic contacts
    search_cutoff = max(cutoff, 4.0)
    tree2 = cKDTree(coords2)
    neighbor_lists = tree2.query_ball_point(coords1, search_cutoff)

    contacts = []
    for i, neighbors in enumerate(neighbor_lists):
        a1 = atoms1[i]
        for j in neighbors:
            a2 = atoms2[j]

            # PAE filter
            if use_pae:
                idx1 = c1_res_to_idx.get(a1['resnum'])
                idx2 = c2_res_to_idx.get(a2['resnum'])
                if (idx1 is not None and idx2 is not None and
                        idx1 < pae_block.shape[0] and idx2 < pae_block.shape[1]):
                    if pae_block[idx1, idx2] >= pae_cutoff:
                        continue

            dist = float(np.linalg.norm(coords1[i] - coords2[j]))

            # Classify the contact type
            ctype = _classify_contact(a1, a2, dist)
            if ctype is None:
                continue  # outside all type-specific cutoffs

            contacts.append({
                'chain1': chain1.id,
                'resnum1': a1['resnum'],
                'resname1': a1['resname'],
                'atom1': a1['atom'],
                'chain2': chain2.id,
                'resnum2': a2['resnum'],
                'resname2': a2['resname'],
                'atom2': a2['atom'],
                'distance': round(dist, 2),
                'contact_type': ctype,
            })

    # Sort by chain1 resnum, then chain2 resnum
    contacts.sort(key=lambda c: (c['resnum1'], c['resnum2'], c['atom1'], c['atom2']))
    return contacts


# ── Contact type classification helpers ──

# Atoms that carry formal positive charge
_POSITIVE_ATOMS = {
    'ARG': {'NH1', 'NH2', 'NE'},
    'LYS': {'NZ'},
    'HIS': {'ND1', 'NE2'},
}

# Atoms that carry formal negative charge
_NEGATIVE_ATOMS = {
    'ASP': {'OD1', 'OD2'},
    'GLU': {'OE1', 'OE2'},
}

# H-bond donor atoms (N and O that can donate, backbone + sidechain)
_DONOR_ELEMENTS = {'N', 'O'}

# H-bond acceptor elements
_ACCEPTOR_ELEMENTS = {'N', 'O', 'S'}

# Hydrophobic residues
_HYDROPHOBIC_RES = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'}

# Aromatic residues (for pi-stacking — classified as hydrophobic here)
_AROMATIC_RES = {'PHE', 'TYR', 'TRP', 'HIS'}

# Cutoffs
_SALT_BRIDGE_CUTOFF = 4.0
_HBOND_CUTOFF = 3.5
_DISULFIDE_CUTOFF = 2.5
_HYDROPHOBIC_CUTOFF = 4.0


def _classify_contact(a1: Dict, a2: Dict, dist: float) -> Optional[str]:
    """
    Classify a single atomic contact.
    Returns contact type string, or None if outside all cutoffs.

    Priority: disulfide > salt_bridge > hydrogen_bond > hydrophobic > vdw
    """
    r1, r2 = a1['resname3'], a2['resname3']
    at1, at2 = a1['atom'], a2['atom']
    e1, e2 = a1['element'], a2['element']

    # Disulfide bond: Cys SG — Cys SG
    if (r1 == 'CYS' and r2 == 'CYS' and at1 == 'SG' and at2 == 'SG'
            and dist < _DISULFIDE_CUTOFF):
        return 'disulfide'

    # Salt bridge: positive ↔ negative charged atoms
    if dist < _SALT_BRIDGE_CUTOFF:
        pos1 = at1 in _POSITIVE_ATOMS.get(r1, set())
        neg1 = at1 in _NEGATIVE_ATOMS.get(r1, set())
        pos2 = at2 in _POSITIVE_ATOMS.get(r2, set())
        neg2 = at2 in _NEGATIVE_ATOMS.get(r2, set())
        if (pos1 and neg2) or (neg1 and pos2):
            return 'salt_bridge'

    # Hydrogen bond: N/O donor to N/O/S acceptor (heavy-atom distance proxy)
    if dist < _HBOND_CUTOFF:
        if (e1 in _DONOR_ELEMENTS and e2 in _ACCEPTOR_ELEMENTS) or \
           (e2 in _DONOR_ELEMENTS and e1 in _ACCEPTOR_ELEMENTS):
            # Exclude backbone C=O to C=O (both are O but no H to donate)
            # In practice, N-O and O-N pairs are the dominant H-bonds
            return 'hydrogen_bond'

    # Hydrophobic: C-C contacts between apolar residues/atoms
    if dist < _HYDROPHOBIC_CUTOFF:
        if e1 == 'C' and e2 == 'C':
            # Both residues are hydrophobic/aromatic, or both atoms are
            # sidechain carbons (not backbone C, CA, or carbonyl C for polar res)
            if (r1 in _HYDROPHOBIC_RES or r1 in _AROMATIC_RES or
                    r2 in _HYDROPHOBIC_RES or r2 in _AROMATIC_RES):
                return 'hydrophobic'
            # C-C between non-hydrophobic residues — still apolar interaction
            # if both are sidechain carbons (CB, CG, CD, CE, CZ)
            backbone_c = {'C', 'CA'}
            if at1 not in backbone_c and at2 not in backbone_c:
                return 'hydrophobic'

    # General van der Waals — catch-all within the original cutoff
    if dist < _SALT_BRIDGE_CUTOFF:  # use widest cutoff for vdw
        return 'vdw'

    return None  # outside all cutoffs


def compute_bsa_freesasa(model, chain1_id: str, chain2_id: str,
                          probe_radius: float = SASA_PROBE_RADIUS) -> Optional[float]:
    """
    Compute BSA between two chains using FreeSASA.
    BSA = SASA(chain1) + SASA(chain2) - SASA(chain1+chain2)
    """
    try:
        # We need to write temporary selections for FreeSASA
        # Get the structure as a freesasa-compatible object
        atoms_1, atoms_2, atoms_12 = [], [], []
        coords_1, coords_2, coords_12 = [], [], []
        radii_1, radii_2, radii_12 = [], [], []

        # Default VdW radii (simplified)
        VDW_RADII = {
            'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'SE': 1.9,
            'H': 1.2, 'P': 1.8,
        }

        def collect_atoms(chain, atom_list, coord_list, radius_list):
            for residue in chain.get_residues():
                if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
                    continue
                for atom in residue.get_atoms():
                    if atom.element == 'H':
                        continue
                    atom_list.append(atom)
                    coord_list.append(atom.get_vector().get_array())
                    radius_list.append(VDW_RADII.get(atom.element.upper(), 1.7))

        chain1 = model[chain1_id]
        chain2 = model[chain2_id]

        collect_atoms(chain1, atoms_1, coords_1, radii_1)
        collect_atoms(chain2, atoms_2, coords_2, radii_2)

        if not coords_1 or not coords_2:
            return None

        # SASA of chain1 alone
        coords_arr_1 = np.array(coords_1).flatten().tolist()
        result_1 = freesasa.calcCoord(coords_arr_1, radii_1,
                                       freesasa.Parameters({'probe-radius': probe_radius}))
        sasa_1 = result_1.totalArea()

        # SASA of chain2 alone
        coords_arr_2 = np.array(coords_2).flatten().tolist()
        result_2 = freesasa.calcCoord(coords_arr_2, radii_2,
                                       freesasa.Parameters({'probe-radius': probe_radius}))
        sasa_2 = result_2.totalArea()

        # SASA of complex
        coords_combined = coords_arr_1 + coords_arr_2
        radii_combined = radii_1 + radii_2
        result_12 = freesasa.calcCoord(coords_combined, radii_combined,
                                        freesasa.Parameters({'probe-radius': probe_radius}))
        sasa_12 = result_12.totalArea()

        bsa = (sasa_1 + sasa_2 - sasa_12) / 2.0  # Divide by 2 for interface area
        return round(bsa, 1)

    except Exception as e:
        logger.error(f"BSA calculation failed for {chain1_id}-{chain2_id}: {e}")
        return None
