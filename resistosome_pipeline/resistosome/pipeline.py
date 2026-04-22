"""
Main pipeline: per-structure analysis and cross-replicate aggregation.
"""

import json
import logging
import os
import math
import numpy as np
import pandas as pd
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple

from .config import CONTACT_CUTOFF_A, PAE_CUTOFF_A, K_CONFIDENCE
from .io_utils import (StructureInput, discover_inputs,
                       normalize_id_from_filename, load_protein_names)
from .structure import (load_structure, get_protein_chains, get_ca_coords,
                        get_sequence, chain_centroid, get_nterm_coord)
from .af3 import (load_summary_confidences, get_iptm, load_full_confidences,
                   get_pae_matrix, get_chain_residue_ranges, get_pae_for_chain_pair)
from .contacts import count_contacts, extract_contact_details, compute_bsa_freesasa
from .geometry import (compute_interface_angles, compute_symmetry_rmsd,
                       compute_nterm_distances_allpairs, compute_pore_axis,
                       determine_ring_order, validate_ring_adjacency)
from .helix import (run_dssp, identify_helices, compute_helix_angle,
                     compute_helix_pore_angle, helix_linearity,
                     compute_helix_length_angstrom, get_cc_helices,
                     split_helix_at_kinks)
from .hmm import (extract_sequence_from_chain, extract_sequences_from_structures,
                  run_hmmsearch, parse_domtbl, get_alpha1_bounds_for_structure,
                  build_alpha1_segment, find_true_nterm_resnum)
from .biophysics import (cc_domain_hydrophobicity, hydrophobic_moment_eisenberg)
from .stats import lcb, ucb, mean_sd

logger = logging.getLogger(__name__)


def _mean_plddt_for_segment(chain, start_resnum: int, end_resnum: int) -> Optional[float]:
    """Compute mean pLDDT (from Cα B-factors) for a residue range."""
    plddt_vals = []
    for residue in chain.get_residues():
        if residue.id[0] != ' ':
            continue
        rn = residue.id[1]
        if rn < start_resnum or rn > end_resnum:
            continue
        for atom in residue.get_atoms():
            if atom.name == 'CA':
                plddt_vals.append(atom.bfactor)
                break
    if not plddt_vals:
        return None
    return round(float(np.mean(plddt_vals)), 2)


def load_annotations(mhd_csv: Optional[str], ploop_csv: Optional[str],
                     nbarc_csv: Optional[str] = None) -> Tuple[Dict, Dict, Dict]:
    """
    Load MHD, P-loop, and NB-ARC annotation CSVs.
    Returns three dicts mapping normalized_id -> (start, end).
    """
    mhd_annot = {}
    ploop_annot = {}
    nbarc_annot = {}

    if mhd_csv and os.path.exists(mhd_csv):
        df = pd.read_csv(mhd_csv)
        for _, row in df.iterrows():
            nid = str(row['ID_normalized']).strip().lower()
            mhd_annot[nid] = (int(row['start']), int(row['end']))

    if ploop_csv and os.path.exists(ploop_csv):
        df = pd.read_csv(ploop_csv)
        for _, row in df.iterrows():
            nid = str(row['ID_normalized']).strip().lower()
            ploop_annot[nid] = (int(row['start']), int(row['end']))

    if nbarc_csv and os.path.exists(nbarc_csv):
        df = pd.read_csv(nbarc_csv)
        for _, row in df.iterrows():
            nid = str(row['ID_normalized']).strip().lower()
            nbarc_annot[nid] = (int(row['start']), int(row['end']))

    logger.info(f"Loaded annotations: {len(mhd_annot)} MHD, "
                f"{len(ploop_annot)} P-loop, {len(nbarc_annot)} NB-ARC entries")
    return mhd_annot, ploop_annot, nbarc_annot


def compute_motif_distance(chain, mhd_range: Tuple[int, int],
                            ploop_range: Tuple[int, int]) -> Optional[float]:
    """
    Compute distance between MHD motif centroid and P-loop motif centroid (Cα atoms).
    """
    ca_coords, res_nums, _ = get_ca_coords(chain)
    if len(ca_coords) == 0:
        return None

    rn_to_idx = {rn: i for i, rn in enumerate(res_nums)}

    # MHD centroid
    mhd_coords = []
    for rn in range(mhd_range[0], mhd_range[1] + 1):
        if rn in rn_to_idx:
            mhd_coords.append(ca_coords[rn_to_idx[rn]])
    if not mhd_coords:
        return None

    # P-loop centroid
    ploop_coords = []
    for rn in range(ploop_range[0], ploop_range[1] + 1):
        if rn in rn_to_idx:
            ploop_coords.append(ca_coords[rn_to_idx[rn]])
    if not ploop_coords:
        return None

    mhd_centroid = np.mean(mhd_coords, axis=0)
    ploop_centroid = np.mean(ploop_coords, axis=0)

    dist = float(np.linalg.norm(mhd_centroid - ploop_centroid))
    return round(dist, 2)


def analyze_single_structure(si: StructureInput, work_dir: str,
                              mhd_annot: Dict, ploop_annot: Dict,
                              nbarc_annot: Dict,
                              alpha1_map: Optional[Dict] = None,
                              protein_names: Optional[List[str]] = None,
                              hmm_fallback_dssp: bool = False,
                              hmm_requested: bool = False) -> Optional[Dict]:
    """
    Analyze a single structure (one replicate) and return raw per-structure metrics.
    """
    logger.info(f"Analyzing: {os.path.basename(si.source_path)} (seed={si.seed})")

    # Prepare files (decompress / extract)
    si.prepare(work_dir)

    if si.model_path is None:
        logger.error(f"No model file found for {si.source_path}")
        return None

    # Load structure
    structure = load_structure(si.model_path, si.model_format)
    if structure is None:
        return None

    model = structure[0]
    chains = get_protein_chains(structure)
    n_chains = len(chains)

    if n_chains < 2:
        logger.warning(f"Only {n_chains} protein chain(s) found, skipping")
        return None

    chain_ids = [c.id for c in chains]
    logger.info(f"  Found {n_chains} protein chains: {chain_ids}")

    # Determine spatial ring order (handles non-alphabetical PDB chain naming)
    chains = determine_ring_order(chains)
    chain_ids = [c.id for c in chains]
    logger.info(f"  Ring-ordered chains: {chain_ids}")

    # Validate adjacency distances (QC)
    adj_dists = validate_ring_adjacency(chains)
    logger.info(f"  Adjacent centroid distances: {adj_dists}")
    dist_cv = float(np.std(adj_dists) / np.mean(adj_dists)) if np.mean(adj_dists) > 0 else 0.0
    if dist_cv > 0.15:
        logger.warning(f"  High adjacency distance CV ({dist_cv:.2f}) — ring order may be unreliable")

    # ── ipTM (AF3 only) ──
    iptm = None
    if si.summary_conf_path:
        summary_conf = load_summary_confidences(si.summary_conf_path)
        iptm = get_iptm(summary_conf)

    # ── PAE matrix (AF3 only) ──
    pae_matrix = None
    chain_ranges = None
    if si.confidences_path:
        full_conf = load_full_confidences(si.confidences_path)
        pae_matrix = get_pae_matrix(full_conf)
        if pae_matrix is not None:
            chain_lengths = []
            for c in chains:
                _, rn, _ = get_ca_coords(c)
                chain_lengths.append(len(rn))
            chain_ranges = get_chain_residue_ranges(full_conf, n_chains, chain_lengths)

    # ── Per-interface metrics (adjacent pairs, circular) ──
    contact_counts = []
    contact_details = {}  # per-interface contact details
    bsa_values = []
    adjacent_pairs = []

    for i in range(n_chains):
        j = (i + 1) % n_chains
        c1, c2 = chains[i], chains[j]
        pair_label = f"{c1.id}-{c2.id}"
        adjacent_pairs.append(pair_label)

        # Contacts
        pae_block = None
        c1_resnums = None
        c2_resnums = None
        if pae_matrix is not None and chain_ranges is not None:
            pae_block = get_pae_for_chain_pair(pae_matrix, chain_ranges, c1.id, c2.id)
            _, c1_resnums, _ = get_ca_coords(c1)
            _, c2_resnums, _ = get_ca_coords(c2)

        n_contacts = count_contacts(c1, c2, cutoff=CONTACT_CUTOFF_A,
                                     pae_block=pae_block, pae_cutoff=PAE_CUTOFF_A,
                                     chain1_resnums=c1_resnums,
                                     chain2_resnums=c2_resnums)
        contact_counts.append(n_contacts)

        # Detailed contacts
        details = extract_contact_details(c1, c2, cutoff=CONTACT_CUTOFF_A,
                                           pae_block=pae_block, pae_cutoff=PAE_CUTOFF_A,
                                           chain1_resnums=c1_resnums,
                                           chain2_resnums=c2_resnums)
        contact_details[pair_label] = details

        # Update total contact count from classified contacts
        contact_counts[-1] = len(details)

        # BSA
        bsa = compute_bsa_freesasa(model, c1.id, c2.id)
        bsa_values.append(bsa)

    # ── Per-type contact counts across interfaces ──
    CONTACT_TYPES = ['hydrogen_bond', 'salt_bridge', 'hydrophobic', 'disulfide', 'vdw']
    per_type_counts = {ct: [] for ct in CONTACT_TYPES}  # per-interface counts
    for pair in adjacent_pairs:
        details = contact_details[pair]
        type_counts = {}
        for c in details:
            ct = c['contact_type']
            type_counts[ct] = type_counts.get(ct, 0) + 1
        for ct in CONTACT_TYPES:
            per_type_counts[ct].append(type_counts.get(ct, 0))

    # ── Interface angles ──
    interface_angles = compute_interface_angles(chains)
    sd_theta = float(np.std(interface_angles, ddof=0)) if interface_angles else None

    # ── Symmetry RMSD ──
    sym_rmsd = compute_symmetry_rmsd(chains)

    # ── HMM-based α1 boundary (needed before N-terminal distances) ──
    hmm_alpha1_end = None
    hmm_alpha1_start = None
    true_nterm_resnum = None  # MADA-based N-terminal override for D_APEX
    if alpha1_map:
        bounds = get_alpha1_bounds_for_structure(
            alpha1_map, si.base_id, protein_names=protein_names)
        if bounds is not None:
            hmm_alpha1_start, hmm_alpha1_end = bounds
            # Find true N-terminal: Met at MADA start or closest upstream Met
            # Use chain A as reference (all protomers are identical in sequence)
            true_nterm_resnum = find_true_nterm_resnum(chains[0], hmm_alpha1_start)
            if true_nterm_resnum != hmm_alpha1_start:
                logger.info(f"  MADA HMM hit: residues {hmm_alpha1_start}-{hmm_alpha1_end}, "
                             f"true N-terminal Met at residue {true_nterm_resnum}")
            logger.info(f"  HMM-defined α1: residues {true_nterm_resnum}-{hmm_alpha1_end}")
        else:
            if hmm_fallback_dssp:
                logger.warning(f"  No HMM hit for {si.base_id}, falling back to DSSP for α1")
            else:
                logger.warning(f"  No HMM hit for {si.base_id}, α1 metrics will be NA")

    # ── N-terminal distances (all pairwise) ──
    # Use MADA-based N-terminal if available, otherwise first residue
    nterm_dist_dict, nterm_mean = compute_nterm_distances_allpairs(
        chains, nterm_resnum=true_nterm_resnum)

    # ── DSSP and helix analysis ──
    dssp_data = run_dssp(structure, si.model_path, si.model_format)

    # Compute pore axis for THETA_APEX (α1 tilt relative to central symmetry axis)
    pore_axis = compute_pore_axis(chains)

    # Look up NB-ARC boundary for this protein
    nbarc_range = nbarc_annot.get(si.base_id)
    nbarc_start = nbarc_range[0] if nbarc_range else None
    if nbarc_start is not None:
        logger.info(f"  NB-ARC starts at residue {nbarc_start}, CC domain: 1-{nbarc_start - 1}")
    else:
        logger.warning(f"  No NB-ARC annotation found for {si.base_id}, "
                        "falling back to first 4 helices for CC domain")

    theta_apex_values = []  # per-protomer α1 tilt vs pore axis
    l_apex_values = []      # per-protomer α1 length (residue count)
    h_abs_values = []       # per-protomer CC-domain hydrophobicity
    mu_h_values = []        # per-protomer α1 amphipathic moment
    alpha1_linearity_values = []  # per-protomer R² of α1 Cα linearity
    alpha1_plddt_values = []      # per-protomer mean pLDDT of α1 region
    d_cc_nb_linker_values = []  # per-protomer CC-to-NB-ARC linker distance (Å)
    l_cc_nb_linker_values = []  # per-protomer CC-to-NB-ARC linker length (residues)
    cc_helix_info = {}      # per-chain CC-domain helix metadata

    # Determine α1 source
    if hmm_alpha1_end is not None:
        alpha1_source = 'hmm'
    elif hmm_requested and not hmm_fallback_dssp:
        # HMM was requested but no hit for this protein, no fallback → all NA
        alpha1_source = 'none'
    else:
        alpha1_source = 'dssp'

    # Linearity threshold: below this R², THETA_APEX is unreliable
    LINEARITY_THRESHOLD = 0.9

    for chain in chains:
        # ── α1 definition ──
        alpha1 = None

        if alpha1_source == 'hmm':
            alpha1 = build_alpha1_segment(chain, hmm_alpha1_end,
                                          start_resnum=true_nterm_resnum)
        elif alpha1_source == 'dssp':
            # Fall back to DSSP first helix, with kink splitting
            if dssp_data is not None:
                all_helices_tmp = identify_helices(chain, dssp_data)
                cc_helices_tmp = get_cc_helices(all_helices_tmp, nbarc_start=nbarc_start)
                if len(cc_helices_tmp) >= 1:
                    # Split first helix at kinks to separate α1 from α2
                    sub_helices = split_helix_at_kinks(chain, cc_helices_tmp[0])
                    alpha1 = sub_helices[0]
        # alpha1_source == 'none': alpha1 stays None

        # ── DSSP-based CC-domain helices (for H_ABS, cc_helix_info, linker) ──
        # Apply kink splitting to all CC helices for more accurate helix counting
        cc_helices = []
        if dssp_data is not None:
            all_helices = identify_helices(chain, dssp_data)
            raw_cc = get_cc_helices(all_helices, nbarc_start=nbarc_start)
            for h in raw_cc:
                cc_helices.extend(split_helix_at_kinks(chain, h))

        # ── α1-based metrics ──
        if alpha1 is not None:
            # Linearity: R² of PCA on α1 Cα atoms
            lin = helix_linearity(chain, alpha1)
            alpha1_linearity_values.append(lin)

            # Mean pLDDT of α1 region (from Cα B-factors)
            plddt = _mean_plddt_for_segment(chain, alpha1.start_resnum, alpha1.end_resnum)
            alpha1_plddt_values.append(plddt)

            # THETA_APEX: only if linearity is above threshold
            if pore_axis is not None and lin is not None and lin >= LINEARITY_THRESHOLD:
                angle = compute_helix_pore_angle(chain, alpha1, pore_axis)
            else:
                angle = None
                if lin is not None and lin < LINEARITY_THRESHOLD:
                    logger.debug(f"  Chain {chain.id}: α1 linearity {lin:.3f} < {LINEARITY_THRESHOLD}, "
                                  "THETA_APEX set to NA")
            theta_apex_values.append(angle)

            # L_APEX: residue count of α1
            l_apex_values.append(alpha1.length)

            # MU_H: amphipathic moment of α1
            mu = hydrophobic_moment_eisenberg(alpha1.sequence)
            mu_h_values.append(mu)
        else:
            theta_apex_values.append(None)
            l_apex_values.append(None)
            mu_h_values.append(None)
            alpha1_linearity_values.append(None)
            alpha1_plddt_values.append(None)

        # ── CC-domain metrics (always DSSP-based) ──
        if len(cc_helices) >= 1:
            cc_seqs = [h.sequence for h in cc_helices]
            h_abs = cc_domain_hydrophobicity(cc_seqs, n_helices=len(cc_helices))
            h_abs_values.append(h_abs)

            # CC-to-NB-ARC linker metrics
            if nbarc_start is not None:
                last_cc_helix = cc_helices[-1]
                ca_coords_arr, res_nums, _ = get_ca_coords(chain)
                rn_to_idx = {rn: i for i, rn in enumerate(res_nums)}

                last_helix_end_idx = rn_to_idx.get(last_cc_helix.end_resnum)
                nbarc_start_idx = rn_to_idx.get(nbarc_start)

                if last_helix_end_idx is not None and nbarc_start_idx is not None:
                    d_linker = float(np.linalg.norm(
                        ca_coords_arr[last_helix_end_idx] - ca_coords_arr[nbarc_start_idx]))
                    d_cc_nb_linker_values.append(round(d_linker, 2))
                else:
                    d_cc_nb_linker_values.append(None)

                l_linker = nbarc_start - last_cc_helix.end_resnum - 1
                l_cc_nb_linker_values.append(l_linker)
            else:
                d_cc_nb_linker_values.append(None)
                l_cc_nb_linker_values.append(None)
        else:
            h_abs_values.append(None)
            d_cc_nb_linker_values.append(None)
            l_cc_nb_linker_values.append(None)

        # CC-domain helix info (metadata only, always DSSP)
        chain_helix_data = []
        for hi, helix in enumerate(cc_helices):
            chain_helix_data.append({
                'helix_index': hi + 1,
                'start_resnum': helix.start_resnum,
                'end_resnum': helix.end_resnum,
                'length': helix.length,
                'sequence': helix.sequence,
            })
        # Prepend HMM α1 info if available
        if hmm_alpha1_end is not None and alpha1 is not None:
            a1_lin = alpha1_linearity_values[-1] if alpha1_linearity_values else None
            a1_plddt = alpha1_plddt_values[-1] if alpha1_plddt_values else None
            chain_helix_data.insert(0, {
                'helix_index': 0,
                'label': 'α1_hmm',
                'start_resnum': alpha1.start_resnum,
                'end_resnum': alpha1.end_resnum,
                'length': alpha1.length,
                'sequence': alpha1.sequence,
                'linearity_r2': a1_lin,
                'mean_plddt': a1_plddt,
            })
        cc_helix_info[chain.id] = chain_helix_data

    if not theta_apex_values:
        theta_apex_values = [None] * n_chains
        l_apex_values = [None] * n_chains
        h_abs_values = [None] * n_chains
        mu_h_values = [None] * n_chains
        alpha1_linearity_values = [None] * n_chains
        alpha1_plddt_values = [None] * n_chains
        d_cc_nb_linker_values = [None] * n_chains
        l_cc_nb_linker_values = [None] * n_chains

    # ── MHD-to-P-loop distance ──
    d_mhd_p_values = []
    mhd_range = mhd_annot.get(si.base_id)
    ploop_range = ploop_annot.get(si.base_id)

    if mhd_range and ploop_range:
        for chain in chains:
            d = compute_motif_distance(chain, mhd_range, ploop_range)
            d_mhd_p_values.append(d)
    else:
        d_mhd_p_values = [None] * n_chains

    # ── Assemble per-structure results ──
    # Within-structure summaries
    contacts_lcb_within = lcb([float(c) for c in contact_counts])
    bsa_lcb_within = lcb([b for b in bsa_values if b is not None])

    # Per-type contact LCBs within structure
    per_type_lcb = {}
    per_type_raw = {}
    for ct in CONTACT_TYPES:
        counts = [float(c) for c in per_type_counts[ct]]
        per_type_lcb[ct] = lcb(counts)
        per_type_raw[ct] = dict(zip(adjacent_pairs, per_type_counts[ct]))

    def safe_mean(vals):
        clean = [v for v in vals if v is not None and not math.isnan(v)]
        return float(np.mean(clean)) if clean else None

    def safe_sd(vals):
        clean = [v for v in vals if v is not None and not math.isnan(v)]
        return round(float(np.std(clean, ddof=0)), 4) if len(clean) >= 2 else None

    results = {
        'model_name': os.path.basename(si.source_path),
        'base_id': si.base_id,
        'seed': si.seed,
        'n_chains': n_chains,
        'chain_ids': chain_ids,  # in ring-adjacent order
        'ring_adjacency_distances': adj_dists,
        'ring_adjacency_cv': round(dist_cv, 4),
        'timestamp': datetime.now().isoformat(),
        'alpha1_source': alpha1_source,  # 'hmm', 'dssp', or 'none'
        'true_nterm_resnum': true_nterm_resnum,  # MADA-based N-terminal, or None

        # Raw per-interface values
        'raw_contact_counts': dict(zip(adjacent_pairs, contact_counts)),
        'raw_bsa_values': dict(zip(adjacent_pairs,
                                    [b if b is not None else None for b in bsa_values])),
        'raw_interface_angles': dict(zip(adjacent_pairs, interface_angles)),

        # Within-structure penalized/summarized values
        'iptm': iptm,
        'contacts_lcb_within': contacts_lcb_within,
        'bsa_lcb_within': bsa_lcb_within,
        'sd_theta_rot': round(sd_theta, 4) if sd_theta is not None else None,
        'symmetry_rmsd': sym_rmsd,
        'd_apex_mean': nterm_mean,

        # Per-type contact LCBs within structure
        'contacts_hbond_lcb': per_type_lcb.get('hydrogen_bond'),
        'contacts_salt_bridge_lcb': per_type_lcb.get('salt_bridge'),
        'contacts_hydrophobic_lcb': per_type_lcb.get('hydrophobic'),
        'contacts_disulfide_lcb': per_type_lcb.get('disulfide'),
        'contacts_vdw_lcb': per_type_lcb.get('vdw'),

        # Raw per-type counts per interface
        'raw_contact_counts_by_type': per_type_raw,

        # Per-protomer values (for cross-rep aggregation)
        'theta_apex_per_protomer': theta_apex_values,
        'l_apex_per_protomer': l_apex_values,
        'h_abs_per_protomer': h_abs_values,
        'mu_h_per_protomer': mu_h_values,
        'alpha1_linearity_per_protomer': alpha1_linearity_values,
        'alpha1_plddt_per_protomer': alpha1_plddt_values,
        'd_mhd_p_per_protomer': d_mhd_p_values,
        'd_cc_nb_linker_per_protomer': d_cc_nb_linker_values,
        'l_cc_nb_linker_per_protomer': l_cc_nb_linker_values,

        # Within-structure summaries for per-protomer metrics
        'sd_theta_apex': safe_sd(theta_apex_values),
        'theta_apex_mean': safe_mean(theta_apex_values),
        'l_apex_mean': safe_mean(l_apex_values),
        'h_abs_mean': safe_mean(h_abs_values),
        'mu_h_mean': safe_mean(mu_h_values),
        'alpha1_linearity_mean': safe_mean(alpha1_linearity_values),
        'alpha1_plddt_mean': safe_mean(alpha1_plddt_values),
        'd_mhd_p_mean': safe_mean(d_mhd_p_values),
        'd_cc_nb_linker_mean': safe_mean(d_cc_nb_linker_values),
        'd_cc_nb_linker_sd': safe_sd(d_cc_nb_linker_values),
        'l_cc_nb_linker_mean': safe_mean(l_cc_nb_linker_values),
        'l_cc_nb_linker_sd': safe_sd(l_cc_nb_linker_values),

        # Detailed N-terminal distances
        'nterm_distances': nterm_dist_dict,

        # Detailed contact lists per interface
        'contact_details': contact_details,

        # Contact type summary per interface
        'contact_type_summary': {
            pair: {ctype: sum(1 for c in contacts if c['contact_type'] == ctype)
                   for ctype in ['salt_bridge', 'hydrogen_bond', 'disulfide',
                                 'hydrophobic', 'vdw']
                   if any(c['contact_type'] == ctype for c in contacts)}
            for pair, contacts in contact_details.items()
        },

        # CC-domain helix info (first 4 helices per protomer)
        'cc_helix_info': cc_helix_info,
    }

    return results


def aggregate_replicates(rep_results: List[Dict]) -> Dict:
    """
    Aggregate metrics across replicates with LCB/UCB penalties.
    """
    base_id = rep_results[0]['base_id']
    n_reps = len(rep_results)
    seeds = [r['seed'] for r in rep_results]

    out = {
        'base_id': base_id,
        'n_replicates': n_reps,
        'seeds': seeds,
        'n_chains': rep_results[0]['n_chains'],
    }

    # ── ipTM_LCB: LCB across replicates ──
    iptm_vals = [r['iptm'] for r in rep_results if r['iptm'] is not None]
    out['ipTM_LCB'] = lcb(iptm_vals) if iptm_vals else None
    m, s = mean_sd(iptm_vals)
    out['ipTM_mean'] = m
    out['ipTM_sd'] = s

    # ── Sum_CONTACTS: LCB of within-structure LCBs ──
    contacts_vals = [r['contacts_lcb_within'] for r in rep_results
                     if r['contacts_lcb_within'] is not None]
    out['Sum_CONTACTS'] = lcb(contacts_vals) if contacts_vals else None
    m, s = mean_sd(contacts_vals)
    out['Sum_CONTACTS_mean'] = m
    out['Sum_CONTACTS_sd'] = s

    # ── Per-type contacts: LCB of within-structure LCBs ──
    _TYPE_KEYS = [
        ('contacts_hbond_lcb', 'CONTACTS_HBOND'),
        ('contacts_salt_bridge_lcb', 'CONTACTS_SALT_BRIDGE'),
        ('contacts_hydrophobic_lcb', 'CONTACTS_HYDROPHOBIC'),
        ('contacts_disulfide_lcb', 'CONTACTS_DISULFIDE'),
        ('contacts_vdw_lcb', 'CONTACTS_VDW'),
    ]
    for src_key, out_prefix in _TYPE_KEYS:
        vals = [r[src_key] for r in rep_results if r.get(src_key) is not None]
        out[out_prefix] = lcb(vals) if vals else None
        m, s = mean_sd(vals)
        out[f'{out_prefix}_mean'] = m
        out[f'{out_prefix}_sd'] = s

    # ── BSA_INT_PROTO: LCB of within-structure LCBs ──
    bsa_vals = [r['bsa_lcb_within'] for r in rep_results
                if r['bsa_lcb_within'] is not None]
    out['BSA_INT_PROTO'] = lcb(bsa_vals) if bsa_vals else None
    m, s = mean_sd(bsa_vals)
    out['BSA_INT_PROTO_mean'] = m
    out['BSA_INT_PROTO_sd'] = s

    # ── SD_THETA_ROT: UCB across replicates ──
    theta_sd_vals = [r['sd_theta_rot'] for r in rep_results
                     if r['sd_theta_rot'] is not None]
    out['SD_THETA_ROT'] = ucb(theta_sd_vals) if theta_sd_vals else None
    m, s = mean_sd(theta_sd_vals)
    out['SD_THETA_ROT_mean'] = m
    out['SD_THETA_ROT_sd'] = s

    # ── S_PROTO: UCB across replicates ──
    sym_vals = [r['symmetry_rmsd'] for r in rep_results
                if r['symmetry_rmsd'] is not None]
    out['S_PROTO'] = ucb(sym_vals) if sym_vals else None
    m, s = mean_sd(sym_vals)
    out['S_PROTO_mean'] = m
    out['S_PROTO_sd'] = s

    # ── D_APEX: UCB across replicates ──
    dapex_vals = [r['d_apex_mean'] for r in rep_results
                  if r['d_apex_mean'] is not None]
    out['D_APEX'] = ucb(dapex_vals) if dapex_vals else None
    m, s = mean_sd(dapex_vals)
    out['D_APEX_mean'] = m
    out['D_APEX_sd'] = s

    # ── SD_THETA_APEX: UCB across replicates (variability penalty) ──
    sd_theta_apex_vals = [r['sd_theta_apex'] for r in rep_results
                          if r['sd_theta_apex'] is not None]
    out['SD_THETA_APEX'] = ucb(sd_theta_apex_vals) if sd_theta_apex_vals else None
    m, s = mean_sd(sd_theta_apex_vals)
    out['SD_THETA_APEX_mean'] = m
    out['SD_THETA_APEX_sd'] = s

    # ── THETA_APEX mean (supplementary, across replicates) ──
    theta_apex_mean_vals = [r['theta_apex_mean'] for r in rep_results
                            if r['theta_apex_mean'] is not None]
    m, s = mean_sd(theta_apex_mean_vals)
    out['THETA_APEX_mean'] = m
    out['THETA_APEX_sd'] = s

    # ── L_APEX: mean ± SD (no penalty) ──
    l_apex_vals = [r['l_apex_mean'] for r in rep_results
                   if r['l_apex_mean'] is not None]
    m, s = mean_sd(l_apex_vals)
    out['L_APEX_mean'] = m
    out['L_APEX_sd'] = s

    # ── H_ABS: mean ± SD (no penalty) ──
    h_abs_vals = [r['h_abs_mean'] for r in rep_results
                  if r['h_abs_mean'] is not None]
    m, s = mean_sd(h_abs_vals)
    out['H_ABS_mean'] = m
    out['H_ABS_sd'] = s

    # ── MU_H: mean ± SD (no penalty) ──
    mu_h_vals = [r['mu_h_mean'] for r in rep_results
                 if r['mu_h_mean'] is not None]
    m, s = mean_sd(mu_h_vals)
    out['MU_H_mean'] = m
    out['MU_H_sd'] = s

    # ── ALPHA1_LINEARITY: mean ± SD (quality indicator) ──
    lin_vals = [r['alpha1_linearity_mean'] for r in rep_results
                if r['alpha1_linearity_mean'] is not None]
    m, s = mean_sd(lin_vals)
    out['ALPHA1_LINEARITY_mean'] = m
    out['ALPHA1_LINEARITY_sd'] = s

    # ── ALPHA1_PLDDT: mean ± SD (quality indicator) ──
    plddt_vals = [r['alpha1_plddt_mean'] for r in rep_results
                  if r['alpha1_plddt_mean'] is not None]
    m, s = mean_sd(plddt_vals)
    out['ALPHA1_PLDDT_mean'] = m
    out['ALPHA1_PLDDT_sd'] = s

    # ── D_MHD_P: mean ± SD (no penalty) ──
    d_mhd_p_vals = [r['d_mhd_p_mean'] for r in rep_results
                    if r['d_mhd_p_mean'] is not None]
    m, s = mean_sd(d_mhd_p_vals)
    out['D_MHD_P_mean'] = m
    out['D_MHD_P_sd'] = s

    # ── D_CC_NB_LINKER: mean ± SD (no penalty) ──
    d_linker_vals = [r['d_cc_nb_linker_mean'] for r in rep_results
                     if r['d_cc_nb_linker_mean'] is not None]
    m, s = mean_sd(d_linker_vals)
    out['D_CC_NB_LINKER_mean'] = m
    out['D_CC_NB_LINKER_sd'] = s

    # ── L_CC_NB_LINKER: mean ± SD (no penalty) ──
    l_linker_vals = [r['l_cc_nb_linker_mean'] for r in rep_results
                     if r['l_cc_nb_linker_mean'] is not None]
    m, s = mean_sd(l_linker_vals)
    out['L_CC_NB_LINKER_mean'] = m
    out['L_CC_NB_LINKER_sd'] = s

    return out


def _analyze_worker(si: StructureInput, work_dir: str,
                     kwargs: Dict) -> Optional[Dict]:
    """
    Worker function for parallel structure analysis.
    analyze_single_structure handles prepare() internally.
    We just ensure cleanup() runs after.
    """
    try:
        result = analyze_single_structure(si, work_dir, **kwargs)
        return result
    finally:
        si.cleanup()


def run_pipeline(input_dir: str, output_dir: str,
                 mhd_csv: Optional[str] = None,
                 ploop_csv: Optional[str] = None,
                 nbarc_csv: Optional[str] = None,
                 protein_names_file: Optional[str] = None,
                 hmm_path: Optional[str] = None,
                 hmm_domtbl: Optional[str] = None,
                 hmm_fallback_dssp: bool = False,
                 workers: Optional[int] = None,
                 work_dir: Optional[str] = None) -> str:
    """
    Run the full resistosome analysis pipeline.

    Parameters
    ----------
    input_dir : str, directory containing structure files
    output_dir : str, directory for output files
    mhd_csv : str, path to MHD annotation CSV
    ploop_csv : str, path to P-loop annotation CSV
    nbarc_csv : str, path to NB-ARC domain annotation CSV
    protein_names_file : str, path to text file with protein names (one per line)
    hmm_path : str, path to HMM model for α1-helix detection (runs hmmsearch)
    hmm_domtbl : str, path to pre-computed hmmsearch --domtblout file
    work_dir : str, temporary working directory (default: output_dir/tmp)

    Returns
    -------
    str : path to the summary Excel file
    """
    os.makedirs(output_dir, exist_ok=True)
    if work_dir is None:
        work_dir = os.path.join(output_dir, 'tmp')
    os.makedirs(work_dir, exist_ok=True)
    json_dir = os.path.join(output_dir, 'per_structure_json')
    os.makedirs(json_dir, exist_ok=True)

    # Load protein names if provided
    protein_names = None
    if protein_names_file and os.path.exists(protein_names_file):
        protein_names = load_protein_names(protein_names_file)

    # Load annotations
    mhd_annot, ploop_annot, nbarc_annot = load_annotations(mhd_csv, ploop_csv, nbarc_csv)

    # Discover inputs
    groups = discover_inputs(input_dir, protein_names=protein_names)

    if not groups:
        logger.error(f"No structure files found in {input_dir}")
        return ""

    # ── Determine worker count ──
    import multiprocessing
    if workers is None:
        workers = min(multiprocessing.cpu_count(), 8)
    workers = max(1, workers)

    # ── HMM-based α1-helix detection ──
    alpha1_map = None
    if hmm_domtbl and os.path.exists(hmm_domtbl):
        # Use pre-computed domain table
        logger.info(f"Loading pre-computed HMM domain table: {hmm_domtbl}")
        alpha1_map = parse_domtbl(hmm_domtbl)
    elif hmm_path and os.path.exists(hmm_path):
        # Extract sequences from first replicate of each group, run hmmsearch
        logger.info(f"Running HMM-based α1 detection with: {hmm_path}")
        hmm_dir = os.path.join(output_dir, 'hmm')
        os.makedirs(hmm_dir, exist_ok=True)

        # Fast parallel sequence extraction from first replicate of each group
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from .hmm import extract_sequence_fast

        tasks = []  # (base_id, si)
        for base_id, inputs in groups.items():
            tasks.append((base_id, inputs[0]))

        logger.info(f"  Extracting sequences from {len(tasks)} structures...")

        sequences = {}  # base_id -> sequence

        # Use ThreadPoolExecutor for I/O-bound sequence extraction (no pickle issues)
        from concurrent.futures import ThreadPoolExecutor, as_completed as as_completed_t

        def _do_extract(args):
            base_id, si = args
            try:
                si.prepare(work_dir)
                if si.model_path is None:
                    return base_id, None
                seq = extract_sequence_fast(
                    si.model_path, fmt=si.model_format, chain_id='A')
                return base_id, seq
            except Exception as e:
                return base_id, None
            finally:
                si.cleanup()

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(_do_extract, t) for t in tasks]
            done = 0
            for future in as_completed_t(futures):
                done += 1
                base_id, seq = future.result()
                if seq:
                    sequences[base_id] = seq
                if done % 50 == 0:
                    logger.info(f"    {done}/{len(tasks)} sequences extracted")

        logger.info(f"  Extracted {len(sequences)} sequences")

        if sequences:
            fasta_path = os.path.join(hmm_dir, 'sequences.fasta')
            with open(fasta_path, 'w') as f:
                for base_id, seq in sorted(sequences.items()):
                    f.write(f">{base_id}\n{seq}\n")

            domtbl_path = run_hmmsearch(hmm_path, fasta_path, hmm_dir)
            if domtbl_path:
                alpha1_map = parse_domtbl(domtbl_path)
            else:
                logger.warning("hmmsearch failed, falling back to DSSP for α1")
        else:
            logger.warning("No sequences extracted, falling back to DSSP for α1")

    if alpha1_map:
        logger.info(f"α1-helix boundaries available for {len(alpha1_map)} proteins (HMM-based)")
    elif hmm_path or hmm_domtbl:
        logger.warning("HMM was provided but no α1 boundaries found — "
                        "check sequence names match between structures and HMM input")
    else:
        logger.info("Using DSSP for α1-helix detection (no --hmm or --hmm_domtbl provided)")

    # Track whether HMM mode was requested (even if no hits found)
    hmm_requested = bool(hmm_path or hmm_domtbl)

    # Flatten all structures into a task list
    all_tasks = []
    for base_id, inputs in groups.items():
        for si in inputs:
            all_tasks.append((base_id, si))

    total_structures = len(all_tasks)
    logger.info(f"Processing {total_structures} structures across "
                f"{len(groups)} groups with {workers} worker(s)")

    # ── Analyze structures (parallel or serial) ──
    from concurrent.futures import ProcessPoolExecutor, as_completed

    # Shared kwargs for each analysis call
    analysis_kwargs = dict(
        mhd_annot=mhd_annot,
        ploop_annot=ploop_annot,
        nbarc_annot=nbarc_annot,
        alpha1_map=alpha1_map,
        protein_names=protein_names,
        hmm_fallback_dssp=hmm_fallback_dssp,
        hmm_requested=hmm_requested,
    )

    # Collect results keyed by base_id
    results_by_group: Dict[str, List[Tuple[int, Dict]]] = {}
    n_done = 0

    if workers == 1:
        # Serial mode (simpler logging, easier debugging)
        for base_id, si in all_tasks:
            try:
                result = _analyze_worker(si, work_dir, analysis_kwargs)
                if result is not None:
                    results_by_group.setdefault(base_id, []).append(
                        (si.seed, result))
            except Exception as e:
                logger.error(f"Failed to analyze {si.source_path}: {e}",
                              exc_info=True)
            n_done += 1
            if n_done % 10 == 0 or n_done == total_structures:
                logger.info(f"  Progress: {n_done}/{total_structures}")
    else:
        # Parallel mode
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_info = {}
            for base_id, si in all_tasks:
                future = executor.submit(
                    _analyze_worker, si, work_dir, analysis_kwargs)
                future_to_info[future] = (base_id, si)

            for future in as_completed(future_to_info):
                base_id, si = future_to_info[future]
                n_done += 1
                try:
                    result = future.result()
                    if result is not None:
                        results_by_group.setdefault(base_id, []).append(
                            (si.seed, result))
                except Exception as e:
                    logger.error(f"Failed to analyze {si.source_path}: {e}",
                                  exc_info=True)
                if n_done % 10 == 0 or n_done == total_structures:
                    logger.info(f"  Progress: {n_done}/{total_structures}")

    # ── Save per-structure JSONs and aggregate ──
    all_aggregated = []

    for base_id in groups:
        seed_results = results_by_group.get(base_id, [])
        if not seed_results:
            continue

        # Sort by seed for deterministic output
        seed_results.sort(key=lambda x: x[0] if x[0] is not None else 0)
        rep_results = [r for _, r in seed_results]

        # Save per-structure JSONs
        for seed, result in seed_results:
            seed_label = f"_seed{seed}" if seed is not None else ""
            json_fname = f"{base_id}{seed_label}_results.json"
            json_path = os.path.join(json_dir, json_fname)
            with open(json_path, 'w') as f:
                json.dump(result, f, indent=2, default=_json_default)

        # Aggregate replicates
        aggregated = aggregate_replicates(rep_results)
        all_aggregated.append(aggregated)

        # Save aggregated JSON
        agg_json = os.path.join(json_dir, f"{base_id}_aggregated.json")
        with open(agg_json, 'w') as f:
            json.dump(aggregated, f, indent=2, default=_json_default)

    # Build summary table
    if all_aggregated:
        summary_path = build_summary_table(all_aggregated, output_dir)
        logger.info(f"\nSummary table saved to: {summary_path}")
    else:
        summary_path = ""
        logger.warning("No results to summarize")

    # Cleanup temp dir
    import shutil
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir, ignore_errors=True)

    return summary_path


def build_summary_table(aggregated_results: List[Dict], output_dir: str) -> str:
    """Build the summary Excel table from aggregated results."""
    rows = []
    for agg in aggregated_results:
        row = {
            'ID': agg['base_id'],
            'N_replicates': agg['n_replicates'],
            'N_chains': agg['n_chains'],
            'Seeds': str(agg['seeds']),
            # Penalized scores
            'ipTM_LCB': agg.get('ipTM_LCB'),
            'Sum_CONTACTS': agg.get('Sum_CONTACTS'),
            'CONTACTS_HBOND': agg.get('CONTACTS_HBOND'),
            'CONTACTS_SALT_BRIDGE': agg.get('CONTACTS_SALT_BRIDGE'),
            'CONTACTS_HYDROPHOBIC': agg.get('CONTACTS_HYDROPHOBIC'),
            'CONTACTS_DISULFIDE': agg.get('CONTACTS_DISULFIDE'),
            'CONTACTS_VDW': agg.get('CONTACTS_VDW'),
            'BSA_INT_PROTO': agg.get('BSA_INT_PROTO'),
            'SD_THETA_ROT': agg.get('SD_THETA_ROT'),
            'S_PROTO': agg.get('S_PROTO'),
            'D_APEX': agg.get('D_APEX'),
            'SD_THETA_APEX': agg.get('SD_THETA_APEX'),
            # Mean ± SD scores
            'THETA_APEX_mean': agg.get('THETA_APEX_mean'),
            'THETA_APEX_sd': agg.get('THETA_APEX_sd'),
            'L_APEX_mean': agg.get('L_APEX_mean'),
            'L_APEX_sd': agg.get('L_APEX_sd'),
            'H_ABS_mean': agg.get('H_ABS_mean'),
            'H_ABS_sd': agg.get('H_ABS_sd'),
            'MU_H_mean': agg.get('MU_H_mean'),
            'MU_H_sd': agg.get('MU_H_sd'),
            'ALPHA1_LINEARITY_mean': agg.get('ALPHA1_LINEARITY_mean'),
            'ALPHA1_LINEARITY_sd': agg.get('ALPHA1_LINEARITY_sd'),
            'ALPHA1_PLDDT_mean': agg.get('ALPHA1_PLDDT_mean'),
            'ALPHA1_PLDDT_sd': agg.get('ALPHA1_PLDDT_sd'),
            'D_MHD_P_mean': agg.get('D_MHD_P_mean'),
            'D_MHD_P_sd': agg.get('D_MHD_P_sd'),
            'D_CC_NB_LINKER_mean': agg.get('D_CC_NB_LINKER_mean'),
            'D_CC_NB_LINKER_sd': agg.get('D_CC_NB_LINKER_sd'),
            'L_CC_NB_LINKER_mean': agg.get('L_CC_NB_LINKER_mean'),
            'L_CC_NB_LINKER_sd': agg.get('L_CC_NB_LINKER_sd'),
            # Raw means and SDs for penalized metrics
            'ipTM_mean': agg.get('ipTM_mean'),
            'ipTM_sd': agg.get('ipTM_sd'),
            'Sum_CONTACTS_mean': agg.get('Sum_CONTACTS_mean'),
            'Sum_CONTACTS_sd': agg.get('Sum_CONTACTS_sd'),
            'CONTACTS_HBOND_mean': agg.get('CONTACTS_HBOND_mean'),
            'CONTACTS_HBOND_sd': agg.get('CONTACTS_HBOND_sd'),
            'CONTACTS_SALT_BRIDGE_mean': agg.get('CONTACTS_SALT_BRIDGE_mean'),
            'CONTACTS_SALT_BRIDGE_sd': agg.get('CONTACTS_SALT_BRIDGE_sd'),
            'CONTACTS_HYDROPHOBIC_mean': agg.get('CONTACTS_HYDROPHOBIC_mean'),
            'CONTACTS_HYDROPHOBIC_sd': agg.get('CONTACTS_HYDROPHOBIC_sd'),
            'CONTACTS_DISULFIDE_mean': agg.get('CONTACTS_DISULFIDE_mean'),
            'CONTACTS_DISULFIDE_sd': agg.get('CONTACTS_DISULFIDE_sd'),
            'CONTACTS_VDW_mean': agg.get('CONTACTS_VDW_mean'),
            'CONTACTS_VDW_sd': agg.get('CONTACTS_VDW_sd'),
            'BSA_INT_PROTO_mean': agg.get('BSA_INT_PROTO_mean'),
            'BSA_INT_PROTO_sd': agg.get('BSA_INT_PROTO_sd'),
            'SD_THETA_ROT_mean': agg.get('SD_THETA_ROT_mean'),
            'SD_THETA_ROT_sd': agg.get('SD_THETA_ROT_sd'),
            'SD_THETA_APEX_mean': agg.get('SD_THETA_APEX_mean'),
            'SD_THETA_APEX_sd': agg.get('SD_THETA_APEX_sd'),
            'S_PROTO_mean': agg.get('S_PROTO_mean'),
            'S_PROTO_sd': agg.get('S_PROTO_sd'),
            'D_APEX_mean': agg.get('D_APEX_mean'),
            'D_APEX_sd': agg.get('D_APEX_sd'),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    xlsx_path = os.path.join(output_dir, 'resistosome_analysis_summary.xlsx')
    df.to_excel(xlsx_path, index=False, sheet_name='Summary')
    csv_path = os.path.join(output_dir, 'resistosome_analysis_summary.csv')
    df.to_csv(csv_path, index=False)
    return xlsx_path


def _json_default(obj):
    """JSON serializer for numpy types."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, float) and math.isnan(obj):
        return None
    return str(obj)
