"""
Statistical penalty functions: Lower/Upper Confidence Bound.
"""

import math
import numpy as np
from typing import List, Optional, Tuple

from .config import K_CONFIDENCE


def lcb(values: List[float], k: float = K_CONFIDENCE) -> Optional[float]:
    """
    Lower Confidence Bound: LCB = mean - k * SE
    where SE = SD / sqrt(n).
    Used for metrics where higher is better (penalizes noisy high values).
    """
    clean = [v for v in values if v is not None and not math.isnan(v)]
    if not clean:
        return None
    n = len(clean)
    mean = float(np.mean(clean))
    if n == 1:
        return mean
    sd = float(np.std(clean, ddof=1))
    se = sd / math.sqrt(n)
    return round(mean - k * se, 6)


def ucb(values: List[float], k: float = K_CONFIDENCE) -> Optional[float]:
    """
    Upper Confidence Bound: UCB = mean + k * SE
    Used for metrics where lower is better (penalizes noisy low values).
    """
    clean = [v for v in values if v is not None and not math.isnan(v)]
    if not clean:
        return None
    n = len(clean)
    mean = float(np.mean(clean))
    if n == 1:
        return mean
    sd = float(np.std(clean, ddof=1))
    se = sd / math.sqrt(n)
    return round(mean + k * se, 6)


def mean_sd(values: List[float]) -> Tuple[Optional[float], Optional[float]]:
    """
    Simple mean ± SD (no penalty).
    Returns (mean, sd).
    """
    clean = [v for v in values if v is not None and not math.isnan(v)]
    if not clean:
        return None, None
    mean = float(np.mean(clean))
    sd = float(np.std(clean, ddof=1)) if len(clean) > 1 else 0.0
    return round(mean, 6), round(sd, 6)
