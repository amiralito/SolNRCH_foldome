"""
HMM-based N-terminal motif detection for α1-helix boundary definition.

Extracts sequences from structures, runs hmmsearch, and parses
domain table output to define α1-helix boundaries per protein.
"""

import os
import re
import subprocess
import logging
from typing import Dict, List, Optional, Tuple

from Bio.PDB.Chain import Chain

from .structure import STANDARD_AA, THREE_TO_ONE

logger = logging.getLogger(__name__)


# ── Fast sequence extraction (no full structure parsing) ──

def _extract_seq_from_cif_fast(filepath: str, chain_id: str = 'A') -> Optional[str]:
    """
    Extract amino acid sequence for a chain directly from a CIF file
    by scanning for Cα ATOM records. Avoids full structure parsing.
    """
    three_to_one = THREE_TO_ONE
    standard = STANDARD_AA
    residues = {}  # resnum -> resname

    # First pass: find column order from _atom_site loop header
    col_names = []
    with open(filepath) as f:
        in_header = False
        for line in f:
            ls = line.strip()
            if ls.startswith('_atom_site.'):
                col_names.append(ls.split('.')[1])
                in_header = True
            elif in_header:
                break  # End of header, data starts

    if not col_names:
        return None

    col_map = {c: i for i, c in enumerate(col_names)}
    group_idx = col_map.get('group_PDB')
    atom_idx = col_map.get('label_atom_id', col_map.get('auth_atom_id'))
    resname_idx = col_map.get('label_comp_id', col_map.get('auth_comp_id'))
    chain_idx = col_map.get('auth_asym_id', col_map.get('label_asym_id'))
    resnum_idx = col_map.get('auth_seq_id', col_map.get('label_seq_id'))
    ncols = len(col_names)

    if any(x is None for x in [group_idx, atom_idx, resname_idx, chain_idx, resnum_idx]):
        return None

    # Second pass: scan only ATOM lines for chain A CA atoms
    with open(filepath) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            fields = line.split()
            if len(fields) < ncols:
                continue
            if fields[chain_idx] != chain_id:
                continue
            atom = fields[atom_idx].strip('"')
            if atom != 'CA':
                continue
            resname = fields[resname_idx]
            if resname not in standard:
                continue
            try:
                resnum = int(fields[resnum_idx])
            except ValueError:
                continue
            residues[resnum] = resname

    if not residues:
        return None

    return ''.join(three_to_one.get(residues[rn], 'X') for rn in sorted(residues))


def _extract_seq_from_pdb_fast(filepath: str, chain_id: str = 'A') -> Optional[str]:
    """
    Extract amino acid sequence for a chain directly from a PDB file.
    """
    three_to_one = THREE_TO_ONE
    standard = STANDARD_AA
    residues = {}

    with open(filepath) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            chain = line[21]
            if chain != chain_id:
                continue
            resname = line[17:20].strip()
            if resname not in standard:
                continue
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            residues[resnum] = resname

    if not residues:
        return None

    seq = []
    for rn in sorted(residues):
        aa = three_to_one.get(residues[rn], 'X')
        seq.append(aa)
    return ''.join(seq)


def extract_sequence_fast(filepath: str, fmt: str = 'cif',
                           chain_id: str = 'A') -> Optional[str]:
    """
    Fast sequence extraction without BioPython structure parsing.
    Reads only Cα ATOM records for the specified chain.
    """
    if fmt == 'cif':
        return _extract_seq_from_cif_fast(filepath, chain_id)
    else:
        return _extract_seq_from_pdb_fast(filepath, chain_id)


def extract_sequence_from_chain(chain: Chain) -> str:
    """
    Extract the amino acid sequence from a chain (standard residues only).
    Returns 1-letter code sequence.
    """
    seq = []
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        aa = THREE_TO_ONE.get(residue.get_resname(), 'X')
        seq.append(aa)
    return ''.join(seq)


def extract_sequences_from_structures(structure_data: Dict[str, Chain],
                                       output_fasta: str) -> Dict[str, str]:
    """
    Extract sequences from chain A of each structure and write a FASTA file.

    Parameters
    ----------
    structure_data : dict mapping base_id -> Chain (chain A of the structure)
    output_fasta : path to write the FASTA file

    Returns
    -------
    dict mapping base_id -> sequence string
    """
    sequences = {}
    with open(output_fasta, 'w') as f:
        for base_id, chain in sorted(structure_data.items()):
            seq = extract_sequence_from_chain(chain)
            if seq:
                sequences[base_id] = seq
                f.write(f">{base_id}\n{seq}\n")

    logger.info(f"Extracted {len(sequences)} sequences to {output_fasta}")
    return sequences


def run_hmmsearch(hmm_path: str, fasta_path: str,
                  output_dir: str,
                  evalue_cutoff: float = 0.01) -> Optional[str]:
    """
    Run hmmsearch and return path to the domain table output.

    Parameters
    ----------
    hmm_path : path to HMM model file
    fasta_path : path to FASTA file with sequences
    output_dir : directory for output files
    evalue_cutoff : E-value threshold for reporting

    Returns
    -------
    str : path to domain table file, or None if hmmsearch failed
    """
    domtbl_path = os.path.join(output_dir, 'hmmsearch_domtbl.out')
    out_path = os.path.join(output_dir, 'hmmsearch.out')

    cmd = [
        'hmmsearch',
        '--domtblout', domtbl_path,
        '-E', str(evalue_cutoff),
        '--noali',
        hmm_path,
        fasta_path,
    ]

    logger.info(f"Running hmmsearch: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600
        )
        if result.returncode != 0:
            logger.error(f"hmmsearch failed (exit {result.returncode}): {result.stderr}")
            return None

        # Save full output
        with open(out_path, 'w') as f:
            f.write(result.stdout)

        if os.path.exists(domtbl_path):
            n_hits = sum(1 for line in open(domtbl_path) if not line.startswith('#'))
            logger.info(f"hmmsearch completed: {n_hits} domain hits in {domtbl_path}")
            return domtbl_path
        else:
            logger.error("hmmsearch produced no domain table output")
            return None

    except FileNotFoundError:
        logger.error("hmmsearch not found. Install HMMER3: "
                      "conda install -c bioconda hmmer  or  brew install hmmer")
        return None
    except subprocess.TimeoutExpired:
        logger.error("hmmsearch timed out after 600 seconds")
        return None
    except Exception as e:
        logger.error(f"hmmsearch error: {e}")
        return None


def parse_domtbl(domtbl_path: str,
                 evalue_threshold: float = 0.01) -> Dict[str, int]:
    """
    Parse hmmsearch domain table to extract α1-helix end coordinates.

    For each protein, takes the best N-terminal domain hit (domain closest
    to residue 1 with E-value below threshold) and returns the alignment
    end coordinate (ali_to).

    α1-helix is then defined as residue 1 to ali_to.

    Parameters
    ----------
    domtbl_path : path to --domtblout file from hmmsearch
    evalue_threshold : maximum per-domain i-Evalue to consider

    Returns
    -------
    dict mapping target_name (lowercase) -> ali_to (int)
        The end residue of the α1-helix for each protein.
    """
    # domtbl columns (0-indexed, space-delimited):
    #  0: target name
    #  1: target accession
    #  2: tlen
    #  3: query name
    #  4: query accession
    #  5: qlen
    #  6: E-value (full seq)
    #  7: score (full seq)
    #  8: bias (full seq)
    #  9: domain # (this domain)
    # 10: of (total domains)
    # 11: c-Evalue
    # 12: i-Evalue
    # 13: score (domain)
    # 14: bias (domain)
    # 15: hmm from
    # 16: hmm to
    # 17: ali from
    # 18: ali to
    # 19: env from
    # 20: env to
    # 21: acc
    # 22+: description

    hits = {}  # target_name -> list of (ali_from, ali_to, i_evalue)

    with open(domtbl_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue

            target = fields[0]
            try:
                i_evalue = float(fields[12])
                ali_from = int(fields[17])
                ali_to = int(fields[18])
            except (ValueError, IndexError):
                continue

            if i_evalue > evalue_threshold:
                continue

            target_lower = target.lower()
            if target_lower not in hits:
                hits[target_lower] = []
            hits[target_lower].append((ali_from, ali_to, i_evalue))

    # For each protein, pick the hit closest to residue 1
    alpha1_bounds = {}
    for target, domain_hits in hits.items():
        # Sort by ali_from (closest to N-terminus first), then by E-value
        domain_hits.sort(key=lambda x: (x[0], x[2]))
        best = domain_hits[0]
        alpha1_bounds[target] = (best[0], best[1])  # (ali_from, ali_to)

    logger.info(f"Parsed {len(alpha1_bounds)} α1-helix boundaries from {domtbl_path}")

    # Log some stats
    if alpha1_bounds:
        starts = [v[0] for v in alpha1_bounds.values()]
        ends = [v[1] for v in alpha1_bounds.values()]
        logger.info(f"  α1 start range: {min(starts)}-{max(starts)}, "
                     f"end range: {min(ends)}-{max(ends)}")

    return alpha1_bounds


def get_alpha1_bounds_for_structure(alpha1_map: Dict[str, Tuple[int, int]],
                                    base_id: str,
                                    protein_names: Optional[List[str]] = None
                                    ) -> Optional[Tuple[int, int]]:
    """
    Look up the α1-helix HMM alignment bounds for a structure.

    Tries multiple matching strategies:
      1. Exact match on base_id
      2. base_id is substring of a domtbl key
      3. domtbl key is substring of base_id
      4. If protein_names provided, try matching through those

    Parameters
    ----------
    alpha1_map : dict from parse_domtbl (target_name_lower -> (ali_from, ali_to))
    base_id : normalized base ID of the structure
    protein_names : optional list of known protein names for matching

    Returns
    -------
    (int, int) or None : (ali_from, ali_to) of α1-helix, or None if not found
    """
    bid = base_id.lower()

    # 1. Exact match
    if bid in alpha1_map:
        return alpha1_map[bid]

    # 2. base_id is substring of a domtbl key
    for key, val in alpha1_map.items():
        if bid in key:
            return val

    # 3. domtbl key is subset of base_id tokens (word-boundary matching)
    #    e.g. "nrc0" should NOT match "slnrc0_sa" because "nrc0" ∉ {"slnrc0","sa"}
    bid_tokens = set(bid.replace('-', '_').split('_'))
    best_match = None
    best_len = 0
    for key, val in alpha1_map.items():
        key_tokens = set(key.replace('-', '_').split('_'))
        if key_tokens and key_tokens.issubset(bid_tokens) and len(key) > best_len:
            best_match = val
            best_len = len(key)
    if best_match is not None:
        return best_match

    # 4. Try via protein_names (token-based to avoid partial matches)
    if protein_names:
        bid_tokens = set(bid.replace('-', '_').split('_'))
        for name in protein_names:
            nl = name.lower()
            name_tokens = set(nl.replace('-', '_').split('_'))
            if name_tokens and name_tokens.issubset(bid_tokens):
                # Find this protein name in the domtbl
                for key, val in alpha1_map.items():
                    key_tokens = set(key.replace('-', '_').split('_'))
                    if name_tokens.issubset(key_tokens):
                        return val

    return None


def build_alpha1_segment(chain: Chain, end_resnum: int,
                          start_resnum: Optional[int] = None):
    """
    Build a HelixSegment for α1 from start_resnum to end_resnum
    using the chain's actual residue data.

    Parameters
    ----------
    chain : Bio.PDB.Chain
    end_resnum : last residue of α1 (from HMM ali_to)
    start_resnum : first residue of α1 (e.g. true N-terminal Met).
                   If None, uses the chain's first residue.

    Returns a HelixSegment, or None if residues can't be found.
    """
    from .helix import HelixSegment

    residues = []
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        rn = residue.id[1]
        aa = THREE_TO_ONE.get(residue.get_resname(), 'X')
        residues.append((rn, aa))

    if not residues:
        return None

    # Find residues from start_resnum to end_resnum
    indices = []
    resnums = []
    seq_chars = []
    for i, (rn, aa) in enumerate(residues):
        if start_resnum is not None and rn < start_resnum:
            continue
        if rn <= end_resnum:
            indices.append(i)
            resnums.append(rn)
            seq_chars.append(aa)

    if not indices:
        return None

    return HelixSegment(
        start_resnum=resnums[0],
        end_resnum=resnums[-1],
        residue_indices=indices,
        sequence=''.join(seq_chars),
    )


def find_true_nterm_resnum(chain: Chain, mada_ali_from: int) -> int:
    """
    Find the "true" N-terminal residue number based on MADA HMM hit.

    Gene models sometimes have spurious N-terminal extensions before the
    real MADA motif. This function finds the correct start:
      1. If residue at ali_from is Met → use it
      2. Otherwise, find the closest Met upstream of ali_from
      3. If no upstream Met, use ali_from itself

    Parameters
    ----------
    chain : Bio.PDB.Chain
    mada_ali_from : residue number where the MADA HMM alignment starts

    Returns
    -------
    int : residue number of the true N-terminal residue
    """
    residues = []
    for residue in chain.get_residues():
        if residue.get_resname() not in STANDARD_AA or residue.id[0] != ' ':
            continue
        rn = residue.id[1]
        aa = THREE_TO_ONE.get(residue.get_resname(), 'X')
        residues.append((rn, aa))

    if not residues:
        return mada_ali_from

    # Check if the residue at ali_from is Met
    for rn, aa in residues:
        if rn == mada_ali_from:
            if aa == 'M':
                return mada_ali_from
            break

    # Find the closest Met upstream of ali_from
    best_met = None
    for rn, aa in residues:
        if rn >= mada_ali_from:
            break
        if aa == 'M':
            best_met = rn  # Keep updating; last one is closest to ali_from

    if best_met is not None:
        return best_met

    # No upstream Met found — use ali_from itself
    return mada_ali_from
