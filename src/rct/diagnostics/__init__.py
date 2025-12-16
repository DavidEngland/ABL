"""Diagnostics for bias assessment, bootstrap confidence intervals, and regime classification.

Modules:
  bias: Bias ratio B = Ri_g(z_g) / Ri_b and error quantification.
  stability: Regime classification (unstable, near-neutral, stable, supercritical).
  bootstrap: Bootstrap resampling for confidence interval estimation.
"""

from .bias import bias_ratio_diagnostic, geometric_mean_height
from .stability import classify_regime, regime_thresholds
from .bootstrap import bootstrap_bias

__all__ = [
    "bias_ratio_diagnostic",
    "geometric_mean_height",
    "classify_regime",
    "regime_thresholds",
    "bootstrap_bias",
]
