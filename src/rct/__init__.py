"""Richardson Curvature Toolkit (RCT) — Curvature-aware Richardson number corrections.

A reproducible, pedagogical package for diagnosing and correcting systematic bias
in gradient and bulk Richardson number calculations arising from vertical discretization,
curvature, and stability transitions in the atmospheric boundary layer.

Version: 0.0.1-alpha
License: MIT (pending)
"""

__version__ = "0.0.1-alpha"
__author__ = "David England, Dick McNider (UAH), Arastoo Pour-Biazar (UAH)"
__description__ = (
    "Curvature-aware Richardson number corrections for stable and complex boundary layers."
)

from . import core, diagnostics, data, viz, utils

__all__ = [
    "core",
    "diagnostics",
    "data",
    "viz",
    "utils",
    "__version__",
]
