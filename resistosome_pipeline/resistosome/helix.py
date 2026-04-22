"""
Secondary structure analysis: DSSP, helix identification, helix geometry.

Handles AlphaFold 3 CIF files (which use non-standard mmCIF dictionaries)
by converting to PDB format before running DSSP.  Falls back to phi/psi-
based assignment if DSSP still fails.
"""

import logging
import tempfile
import subprocess
import os
import math
import numpy as np
from typing import Dict, List, Tuple, Optional
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select as PDBIOSelect

from .config import HELIX_CODES
from .structure import STANDARD_AA, get_ca_coords, get_sequence

logger = logging.getLogger(__name__)


class HelixSegment:
    """Represents a single helix segment in a chain."""

    def __init__(self, start_resnum: int, end_resnum: int,
                 residue_indices: List[int], sequence: str):
        self.start_resnum = start_resnum
        self.end_resnum = end_resnum
        self.residue_indices = residue_indices  # indices into the chain's residue list
        self.sequence = sequence
        self.length = len(residue_indices)

    def __repr__(self):
        return f"Helix({self.start_resnum}-{self.end_resnum}, len={self.length})"


def _structure_to_pdb(structure: Structure, output_path: str) -> str:
    """
    Convert a BioPython Structure to PDB format.
    Handles AF3 CIF files that DSSP cannot read directly.
    Only exports protein chains with single-character chain IDs
    (AF3 ligand chains with multi-character IDs like 'AA' break PDB format).
    """
    from .structure import STANDARD_AA

    class ProteinChainSelect(PDBIOSelect):
        """Only select protein chains with single-character IDs."""
        def accept_chain(self, chain):
            if len(chain.id) > 1:
                return False
            # Check if chain has protein residues
            for res in chain.get_residues():
                if res.get_resname() in STANDARD_AA and res.id[0] == ' ':
                    return True
            return False

        def accept_residue(self, residue):
            # Accept standard AA residues + HETATM that are part of protein
            return (residue.get_resname() in STANDARD_AA and
                    residue.id[0] == ' ')

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, select=ProteinChainSelect())
    # Ensure HEADER line exists
    with open(output_path, 'r') as f:
        content = f.read()
    if not content.startswith('HEADER'):
        with open(output_path, 'w') as f:
            f.write('HEADER    RESISTOSOME PREDICTION\n')
            f.write(content)
    return output_path


def run_dssp(structure: Structure, model_path: str, fmt: str = 'cif') -> Optional[Dict]:
    """
    Run DSSP on a structure and return the DSSP dictionary.
    For CIF files (especially AF3), converts to PDB first since DSSP/mkdssp
    often cannot parse AF3 mmCIF files (missing mmcif_ma.dic).

    Falls back to Cα-distance-based helix assignment if DSSP fails entirely.

    Returns {(chain_id, res_num): {ss, aa, acc, phi, psi}} or None on failure.
    """
    import tempfile
    tmp_files = []

    def _make_tmp_pdb(suffix='_dssp.pdb'):
        """Create a temp PDB file in the system temp directory."""
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        tmp_files.append(path)
        return path

    try:
        # ── Strategy 1: For CIF, convert to PDB then run DSSP ──
        if fmt == 'cif':
            try:
                tmp_pdb = _make_tmp_pdb('_cif2pdb.pdb')
                _structure_to_pdb(structure, tmp_pdb)
                result = _run_dssp_on_pdb(structure, tmp_pdb)
                if result:
                    logger.info("  DSSP succeeded (CIF→PDB conversion)")
                    return result
            except Exception as e:
                logger.debug(f"DSSP via CIF→PDB conversion failed: {e}")

            # Also try mkdssp subprocess on the converted PDB
            try:
                tmp_pdb2 = _make_tmp_pdb('_cif2pdb_sub.pdb')
                _structure_to_pdb(structure, tmp_pdb2)
                result = _run_dssp_subprocess(structure, tmp_pdb2)
                if result:
                    logger.info("  DSSP succeeded (CIF→PDB subprocess)")
                    return result
            except Exception as e:
                logger.debug(f"Subprocess DSSP on converted CIF failed: {e}")

        # ── Strategy 2: For PDB, run DSSP directly (with HEADER fix) ──
        if fmt == 'pdb':
            actual_path = model_path
            try:
                with open(model_path, 'r') as f:
                    first_line = f.readline()
                if not first_line.startswith('HEADER'):
                    tmp_pdb = _make_tmp_pdb('_header.pdb')
                    with open(model_path, 'r') as f_in, open(tmp_pdb, 'w') as f_out:
                        f_out.write('HEADER    RESISTOSOME PREDICTION\n')
                        f_out.write(f_in.read())
                    actual_path = tmp_pdb
            except Exception:
                pass

            try:
                result = _run_dssp_on_pdb(structure, actual_path)
                if result:
                    logger.info("  DSSP succeeded (direct PDB)")
                    return result
            except Exception as e:
                logger.debug(f"Direct DSSP on PDB failed: {e}")

        # ── Strategy 3: Fallback to Cα-distance-based assignment ──
        logger.info("  Falling back to Cα-distance-based secondary structure assignment")
        result = _assign_ss_from_phipsi(structure)
        if result:
            return result

        logger.warning("All secondary structure assignment methods failed")
        return None

    finally:
        # Always clean up temp files
        for f in tmp_files:
            try:
                if os.path.exists(f):
                    os.remove(f)
            except OSError:
                pass


def _run_dssp_on_pdb(structure: Structure, pdb_path: str) -> Optional[Dict]:
    """
    Run DSSP via BioPython on a PDB file.
    Re-parses the PDB to ensure model matches the file (important when
    PDB was converted from CIF with chain filtering).
    Returns dict or None.
    """
    import warnings
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        pdb_struct = parser.get_structure('dssp_input', pdb_path)
        model = pdb_struct[0]

        # Suppress noisy mkdssp/libcifpp warnings about mmCIF parsing
        # (mkdssp tries CIF first, fails, falls back to PDB — all harmless)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UserWarning,
                                     module=r'Bio\.PDB\.DSSP')
            dssp = DSSP(model, pdb_path, dssp='mkdssp')
        result = {}
        for key in dssp.keys():
            chain_id, res_id = key
            dssp_data = dssp[key]
            result[(chain_id, res_id[1])] = {
                'ss': dssp_data[2],
                'aa': dssp_data[1],
                'acc': dssp_data[3],
                'phi': dssp_data[4],
                'psi': dssp_data[5],
            }
        if result:
            return result
    except Exception as e:
        logger.debug(f"_run_dssp_on_pdb failed: {e}")
    return None


def _run_dssp_subprocess(structure: Structure, pdb_path: str) -> Optional[Dict]:
    """
    Run mkdssp as a subprocess and parse its output directly.
    More robust than BioPython's DSSP wrapper for edge cases.
    """
    try:
        # Run mkdssp
        proc = subprocess.run(
            ['mkdssp', '--output-format', 'dssp', pdb_path],
            capture_output=True, text=True, timeout=60
        )
        if proc.returncode != 0:
            # Try older mkdssp syntax
            proc = subprocess.run(
                ['mkdssp', '-i', pdb_path],
                capture_output=True, text=True, timeout=60
            )
        if proc.returncode != 0:
            return None

        # Parse DSSP output
        result = {}
        in_data = False
        for line in proc.stdout.split('\n'):
            if line.startswith('  #  RESIDUE'):
                in_data = True
                continue
            if not in_data or len(line) < 17:
                continue

            # DSSP fixed-width format
            try:
                resnum_str = line[5:10].strip()
                if not resnum_str:
                    continue
                chain_id = line[11]
                aa = line[13]
                ss = line[16]
                if ss == ' ':
                    ss = '-'

                resnum = int(resnum_str)
                result[(chain_id, resnum)] = {
                    'ss': ss,
                    'aa': aa,
                    'acc': 0,
                    'phi': 0.0,
                    'psi': 0.0,
                }
            except (ValueError, IndexError):
                continue

        return result if result else None
    except Exception as e:
        logger.debug(f"_run_dssp_subprocess failed: {e}")
        return None


def _assign_ss_from_phipsi(structure: Structure) -> Optional[Dict]:
    """
    Fallback secondary structure assignment based on backbone dihedral angles.

    Alpha-helix: phi ≈ -57°, psi ≈ -47° (tolerance ±30°)
    Requires at least 3 consecutive helical residues.

    Uses Cα distances as a proxy when dihedrals can't be computed:
    Alpha-helix: Cα(i) to Cα(i+3) ≈ 5.0-5.5 Å
    """
    from .structure import THREE_TO_ONE

    result = {}
    model = structure[0]

    for chain in model.get_chains():
        residues = []
        for res in chain.get_residues():
            if res.get_resname() not in STANDARD_AA or res.id[0] != ' ':
                continue
            if 'CA' in res:
                residues.append(res)

        if len(residues) < 4:
            continue

        # Get CA coordinates
        ca_coords = np.array([r['CA'].get_vector().get_array() for r in residues])

        # Method: Use i to i+3 Cα distance as helix indicator
        # In an alpha-helix, d(i, i+3) ≈ 5.0-5.5 Å
        # In extended/sheet, d(i, i+3) ≈ 9.5-10.5 Å
        # In coil, variable but typically > 6 Å

        is_helical = [False] * len(residues)

        for i in range(len(residues) - 3):
            d_i3 = float(np.linalg.norm(ca_coords[i] - ca_coords[i + 3]))
            # Also check i to i+2 distance for tighter constraint
            d_i2 = float(np.linalg.norm(ca_coords[i] - ca_coords[i + 2]))

            # Alpha-helix: d(i,i+3) ~ 5.0-5.5, d(i,i+2) ~ 5.2-5.8
            if 4.5 <= d_i3 <= 6.2 and 4.5 <= d_i2 <= 6.5:
                is_helical[i] = True
                is_helical[i + 1] = True
                is_helical[i + 2] = True
                is_helical[i + 3] = True

        # Require at least 3 consecutive helical residues
        # Smooth: remove isolated helical residues
        smoothed = list(is_helical)
        for i in range(len(smoothed)):
            if smoothed[i]:
                # Check if at least 2 neighbors are also helical
                neighbors = 0
                if i > 0 and is_helical[i - 1]:
                    neighbors += 1
                if i < len(is_helical) - 1 and is_helical[i + 1]:
                    neighbors += 1
                if neighbors == 0:
                    smoothed[i] = False

        for i, res in enumerate(residues):
            chain_id = chain.id
            resnum = res.id[1]
            aa = THREE_TO_ONE.get(res.get_resname(), 'X')
            ss = 'H' if smoothed[i] else '-'
            result[(chain_id, resnum)] = {
                'ss': ss,
                'aa': aa,
                'acc': 0,
                'phi': 0.0,
                'psi': 0.0,
            }

    if result:
        # Log summary
        n_helix = sum(1 for v in result.values() if v['ss'] in HELIX_CODES)
        n_total = len(result)
        logger.info(f"  Phi/psi fallback assigned {n_helix}/{n_total} "
                     f"residues as helical ({100*n_helix/n_total:.1f}%)")

    return result if result else None


def identify_helices(chain: Chain, dssp_data: Dict,
                     min_helix_len: int = 3) -> List[HelixSegment]:
    """
    Identify helix segments in a chain using DSSP assignments.
    No merging is applied — DSSP defines helix boundaries directly.

    Parameters
    ----------
    chain : Bio.PDB.Chain
    dssp_data : dict from run_dssp
    min_helix_len : minimum helix length to keep
    """
    chain_id = chain.id
    residues = []
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        rnum = residue.id[1]
        from .structure import THREE_TO_ONE
        aa = THREE_TO_ONE.get(residue.get_resname(), 'X')
        ss = dssp_data.get((chain_id, rnum), {}).get('ss', '-')
        is_helix = ss in HELIX_CODES
        residues.append((rnum, aa, is_helix))

    if not residues:
        return []

    # Find helix segments (contiguous runs of helix residues)
    raw_segments = []
    current_start = None
    for i, (rnum, aa, is_helix) in enumerate(residues):
        if is_helix:
            if current_start is None:
                current_start = i
        else:
            if current_start is not None:
                raw_segments.append((current_start, i - 1))
                current_start = None
    if current_start is not None:
        raw_segments.append((current_start, len(residues) - 1))

    if not raw_segments:
        return []

    # Build HelixSegment objects, filter by min length
    helices = []
    for start_idx, end_idx in raw_segments:
        indices = list(range(start_idx, end_idx + 1))
        resnums = [residues[i][0] for i in indices]
        seq = ''.join(residues[i][1] for i in indices)
        if len(indices) >= min_helix_len:
            helices.append(HelixSegment(
                start_resnum=resnums[0],
                end_resnum=resnums[-1],
                residue_indices=indices,
                sequence=seq,
            ))

    return helices


def get_cc_helices(helices: List[HelixSegment],
                   nbarc_start: Optional[int] = None,
                   n_helices: int = 4) -> List[HelixSegment]:
    """
    Get helices belonging to the CC domain.

    If nbarc_start is provided, returns all helices that START before
    nbarc_start (i.e. within the CC domain).

    If nbarc_start is not available, falls back to taking the first
    n_helices helices.
    """
    if nbarc_start is not None:
        return [h for h in helices if h.start_resnum < nbarc_start]
    return helices[:n_helices]


def split_helix_at_kinks(chain: Chain, helix: HelixSegment,
                          angle_threshold: float = 30.0,
                          window: int = 5,
                          min_segment_len: int = 5) -> List[HelixSegment]:
    """
    Split a helix at kink points where the backbone direction changes sharply.

    This handles the common DSSP artefact where α1 and α2 are fused into
    one continuous helix when the linker between them is very short (2-3 residues).

    Algorithm:
      1. Compute local direction vectors in sliding windows of `window` Cα atoms
      2. At each position, measure the angle between consecutive direction vectors
      3. Split at positions where the angle exceeds `angle_threshold`

    Parameters
    ----------
    chain : Bio.PDB.Chain
    helix : HelixSegment to potentially split
    angle_threshold : degrees; direction change above this triggers a split
    window : number of Cα atoms for local direction estimation
    min_segment_len : minimum residues in a resulting segment

    Returns
    -------
    list of HelixSegment : the original helix if no kink found, or sub-helices
    """
    ca_coords, res_nums, _ = get_ca_coords(chain)
    if len(ca_coords) == 0:
        return [helix]

    resnum_to_idx = {rn: i for i, rn in enumerate(res_nums)}

    # Gather Cα coordinates for this helix
    helix_resnums = []
    helix_coords = []
    helix_chain_indices = []  # indices into chain's residue list
    for rn in range(helix.start_resnum, helix.end_resnum + 1):
        if rn in resnum_to_idx:
            helix_resnums.append(rn)
            helix_coords.append(ca_coords[resnum_to_idx[rn]])
            helix_chain_indices.append(resnum_to_idx[rn])

    n = len(helix_coords)
    if n < 2 * window:
        return [helix]  # Too short to detect kinks

    coords = np.array(helix_coords)

    # Compute local direction vectors using sliding windows
    directions = []
    for i in range(n - window + 1):
        segment = coords[i:i + window]
        # Direction = vector from first to last point in window
        d = segment[-1] - segment[0]
        norm = np.linalg.norm(d)
        if norm > 0:
            directions.append(d / norm)
        else:
            directions.append(np.zeros(3))

    # Compute angle between consecutive direction vectors
    angles = []
    for i in range(len(directions) - 1):
        dot = np.clip(np.dot(directions[i], directions[i + 1]), -1.0, 1.0)
        angle_deg = np.degrees(np.arccos(dot))
        angles.append(angle_deg)

    # Find the position of maximum angle above threshold
    if not angles:
        return [helix]

    max_angle = max(angles)
    if max_angle < angle_threshold:
        return [helix]  # No significant kink

    # Split at the position of maximum kink
    # The angle at index i corresponds to the transition between
    # direction[i] (centered around residue i+window//2) and
    # direction[i+1] (centered around residue i+1+window//2)
    # Split point is at residue index i + window//2 + 1
    max_idx = angles.index(max_angle)
    split_pos = max_idx + window // 2 + 1  # index into helix_resnums

    # Build sub-segments
    from .structure import THREE_TO_ONE as T2O
    segments = []
    ranges = [(0, split_pos), (split_pos, n)]

    for start, end in ranges:
        if end - start < min_segment_len:
            continue
        seg_resnums = helix_resnums[start:end]
        seg_indices = list(range(helix.residue_indices[start],
                                  helix.residue_indices[min(end, len(helix.residue_indices)) - 1] + 1))

        # Get sequence for this segment
        seg_seq = []
        for residue in chain.get_residues():
            if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
                continue
            rn = residue.id[1]
            if rn in set(seg_resnums):
                aa = T2O.get(residue.get_resname(), 'X')
                seg_seq.append(aa)

        segments.append(HelixSegment(
            start_resnum=seg_resnums[0],
            end_resnum=seg_resnums[-1],
            residue_indices=seg_indices,
            sequence=''.join(seg_seq),
        ))

    if not segments:
        return [helix]

    logger.debug(f"Split helix {helix.start_resnum}-{helix.end_resnum} "
                  f"(len={helix.length}) at kink angle {max_angle:.1f}° → "
                  f"{[f'{s.start_resnum}-{s.end_resnum}' for s in segments]}")

    return segments


def helix_axis_vector(chain: Chain, helix: HelixSegment) -> Optional[np.ndarray]:
    """
    Compute the helix axis as a unit vector using PCA on Cα positions.
    The first principal component approximates the helix axis.
    """
    ca_coords, res_nums, _ = get_ca_coords(chain)
    if len(ca_coords) == 0:
        return None

    # Map residue numbers to indices in the ca_coords array
    resnum_to_idx = {rn: i for i, rn in enumerate(res_nums)}

    helix_coords = []
    for rnum in range(helix.start_resnum, helix.end_resnum + 1):
        if rnum in resnum_to_idx:
            helix_coords.append(ca_coords[resnum_to_idx[rnum]])

    if len(helix_coords) < 3:
        return None

    coords = np.array(helix_coords)
    coords_centered = coords - coords.mean(axis=0)

    # PCA: first principal component
    _, _, Vt = np.linalg.svd(coords_centered)
    axis = Vt[0]  # first principal component
    return axis / np.linalg.norm(axis)


def helix_linearity(chain: Chain, helix: HelixSegment) -> Optional[float]:
    """
    Compute the linearity of Cα atoms in a helix segment as the fraction
    of variance explained by the first principal component (R²).

    A perfectly straight helix → R² ≈ 1.0
    A kinked/disordered region → R² << 1.0

    Returns R² between 0.0 and 1.0, or None if too few residues.
    """
    ca_coords, res_nums, _ = get_ca_coords(chain)
    if len(ca_coords) == 0:
        return None

    resnum_to_idx = {rn: i for i, rn in enumerate(res_nums)}

    helix_coords = []
    for rnum in range(helix.start_resnum, helix.end_resnum + 1):
        if rnum in resnum_to_idx:
            helix_coords.append(ca_coords[resnum_to_idx[rnum]])

    if len(helix_coords) < 3:
        return None

    coords = np.array(helix_coords)
    coords_centered = coords - coords.mean(axis=0)

    _, S, _ = np.linalg.svd(coords_centered)
    # Variance explained = S² (eigenvalues are singular values squared)
    var_explained = S ** 2
    total_var = var_explained.sum()
    if total_var == 0:
        return None

    r_squared = float(var_explained[0] / total_var)
    return round(r_squared, 4)


def compute_helix_angle(chain: Chain, helix1: HelixSegment,
                         helix2: HelixSegment) -> Optional[float]:
    """
    Compute the inclination angle between two helix axes (α1 and α2).
    Returns angle in degrees (0° = parallel, 90° = perpendicular).
    """
    v1 = helix_axis_vector(chain, helix1)
    v2 = helix_axis_vector(chain, helix2)
    if v1 is None or v2 is None:
        return None

    cos_angle = abs(np.dot(v1, v2))  # abs because helix direction is arbitrary
    cos_angle = min(cos_angle, 1.0)
    angle_deg = np.degrees(np.arccos(cos_angle))
    return round(angle_deg, 2)


def compute_helix_pore_angle(chain: Chain, helix: HelixSegment,
                              pore_axis: np.ndarray) -> Optional[float]:
    """
    Compute the tilt angle of a helix relative to the resistosome pore axis.

    This measures how much α1 deviates from the central symmetry axis:
      0° = helix runs parallel to the pore axis (tight funnel)
     90° = helix is perpendicular to the pore axis (flat, splayed out)

    Parameters
    ----------
    chain : protomer chain containing the helix
    helix : HelixSegment (typically α1)
    pore_axis : unit vector of the resistosome central symmetry axis

    Returns angle in degrees.
    """
    h_axis = helix_axis_vector(chain, helix)
    if h_axis is None:
        return None

    # abs(dot) because helix direction is arbitrary (N→C or C→N)
    cos_angle = abs(np.dot(h_axis, pore_axis))
    cos_angle = min(cos_angle, 1.0)
    angle_deg = np.degrees(np.arccos(cos_angle))
    return round(angle_deg, 2)


def compute_helix_length_angstrom(chain: Chain, helix: HelixSegment) -> Optional[float]:
    """
    Compute the length of a helix in Angstroms (distance between first and last Cα).
    """
    ca_coords, res_nums, _ = get_ca_coords(chain)
    if len(ca_coords) == 0:
        return None

    resnum_to_idx = {rn: i for i, rn in enumerate(res_nums)}

    first_idx = resnum_to_idx.get(helix.start_resnum)
    last_idx = resnum_to_idx.get(helix.end_resnum)

    if first_idx is None or last_idx is None:
        return None

    dist = float(np.linalg.norm(ca_coords[first_idx] - ca_coords[last_idx]))
    return round(dist, 2)


def get_helix_residue_sequence(chain: Chain, helix: HelixSegment) -> str:
    """Get the 1-letter amino acid sequence of a helix segment."""
    return helix.sequence
