"""
Configuration constants for the resistosome analysis pipeline.
"""

# ── Contact & distance thresholds ──
CONTACT_CUTOFF_A = 3.5       # Angstroms for inter-chain contacts
PAE_CUTOFF_A = 10.0          # Angstroms, PAE filter for AF3 contacts

# ── Penalty parameters ──
K_CONFIDENCE = 1.96          # z-value for 95% CI in LCB/UCB

# ── Secondary structure ──
MIN_LINKER_LENGTH = 2        # Minimum non-helix residues between helices to keep them separate

# ── DSSP helix codes ──
HELIX_CODES = {'H', 'G', 'I'}  # alpha, 3-10, pi helices

# ── Hydrophobicity: Kyte-Doolittle (for H_ABS) ──
KYTE_DOOLITTLE = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}

# ── Hydrophobicity: Eisenberg consensus (for MU_H amphipathic moment) ──
EISENBERG = {
    'A':  0.620, 'R': -2.530, 'N': -0.780, 'D': -0.900, 'C':  0.290,
    'Q': -0.850, 'E': -0.740, 'G':  0.480, 'H': -0.400, 'I':  1.380,
    'L':  1.060, 'K': -1.500, 'M':  0.640, 'F':  1.190, 'P':  0.120,
    'S': -0.180, 'T': -0.050, 'W':  0.810, 'Y':  0.260, 'V':  1.080,
}

# ── Helical wheel angle for amphipathic moment (alpha helix) ──
HELIX_ANGULAR_FREQ = 100.0  # degrees per residue

# ── SASA probe radius ──
SASA_PROBE_RADIUS = 1.4  # Angstroms (water probe)
