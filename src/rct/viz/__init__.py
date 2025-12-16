"""Visualization and reference table generation.

Modules:
  plots: Matplotlib-based plots (bias phase, Ri profiles, curvature).
  tables: Reference table generation over parameter grids.
"""

from .plots import plot_bias_phase, plot_ri_profile, plot_curvature
from .tables import make_reference_table

__all__ = [
    "plot_bias_phase",
    "plot_ri_profile",
    "plot_curvature",
    "make_reference_table",
]
