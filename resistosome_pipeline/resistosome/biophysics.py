"""
Biophysical property calculations: hydrophobicity and amphipathic moment.
"""

import math
import numpy as np
from typing import List, Optional

from .config import KYTE_DOOLITTLE, EISENBERG, HELIX_ANGULAR_FREQ


def mean_hydrophobicity_kd(sequence: str) -> Optional[float]:
    """
    Compute mean Kyte-Doolittle hydrophobicity for a sequence.
    Used for H_ABS (CC-domain core hydrophobicity).
    """
    values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in sequence.upper() if aa in KYTE_DOOLITTLE]
    if not values:
        return None
    return float(np.mean(values))


def hydrophobic_moment_eisenberg(sequence: str,
                                  angle: float = HELIX_ANGULAR_FREQ) -> Optional[float]:
    """
    Compute the hydrophobic moment (μH) of a sequence assuming alpha-helical structure.
    Uses the Eisenberg consensus hydrophobicity scale.

    μH = (1/N) * sqrt( (Σ Hi·sin(i·δ))² + (Σ Hi·cos(i·δ))² )

    where δ = angular frequency (100° for alpha helix), Hi = hydrophobicity of residue i.

    Parameters
    ----------
    sequence : str, amino acid sequence
    angle : float, angular frequency in degrees per residue (default 100° for α-helix)
    """
    delta_rad = math.radians(angle)
    h_values = [EISENBERG.get(aa, 0.0) for aa in sequence.upper()]
    n = len(h_values)
    if n == 0:
        return None

    sum_sin = 0.0
    sum_cos = 0.0
    for i, h in enumerate(h_values):
        theta = i * delta_rad
        sum_sin += h * math.sin(theta)
        sum_cos += h * math.cos(theta)

    mu_h = math.sqrt(sum_sin ** 2 + sum_cos ** 2) / n
    return round(mu_h, 4)


def cc_domain_hydrophobicity(helix_sequences: List[str], n_helices: int = 4) -> Optional[float]:
    """
    Compute mean hydrophobicity of the CC domain (first n_helices helices).
    Concatenates the sequences of the first n helices and computes mean KD.
    """
    if not helix_sequences:
        return None
    cc_seq = ''.join(helix_sequences[:n_helices])
    if not cc_seq:
        return None
    return mean_hydrophobicity_kd(cc_seq)
