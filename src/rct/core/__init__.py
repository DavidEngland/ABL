"""Core Richardson number estimation and curvature diagnostics.

Modules:
  derivatives: Unified first/second derivative and curvature calculations.
  ri_estimators: Gradient & bulk Richardson number computations.
  profiles: MOST stability function families and Newton inversion (ζ from Ri_g).
  curvature: High-level curvature proxy and phase analysis.
  correction_ode: ODE solver for grid-aware correction factors C(z).
"""

from .derivatives import central_with_curvature, second_derivative
from .ri_estimators import ri_gradient, ri_bulk, bias_ratio
from .profiles import make_profile, zeta_from_ri_gradient, profile_catalog
from .curvature import curvature_proxy, compute_neutral_curvature

__all__ = [
    "central_with_curvature",
    "second_derivative",
    "ri_gradient",
    "ri_bulk",
    "bias_ratio",
    "make_profile",
    "zeta_from_ri_gradient",
    "profile_catalog",
    "curvature_proxy",
    "compute_neutral_curvature",
]
