from __future__ import annotations
import numpy as np
from typing import Tuple

def diffs(y: np.ndarray, z: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Centered derivatives dy/dz and d²y/dz² on irregular z using np.gradient.
    """
    y1 = np.gradient(y, z, edge_order=2)
    y2 = np.gradient(y1, z, edge_order=2)
    return y1, y2

def dzeta_terms(z: np.ndarray, L: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute dζ/dz and d²ζ/dz² for ζ=z/L(z) given L(z) on levels z.
    """
    L1, L2 = diffs(L, z)
    dzeta_dz = (L - z * L1) / (L * L)
    d2zeta_dz2 = (-2.0 * L1) / (L * L) - (z * L2) / (L * L) + 2.0 * z * (L1 * L1) / (L * L * L)
    return dzeta_dz, d2zeta_dz2

def map_curvature_z(
    z: np.ndarray,
    L: np.ndarray,
    dRi_dzeta: np.ndarray,
    d2Ri_dzeta2: np.ndarray,
) -> np.ndarray:
    """
    Map ζ-space derivatives to height-space curvature:
      ∂²Ri/∂z² = (dζ/dz)² ∂²Ri/∂ζ² + (d²ζ/dz²) ∂Ri/∂ζ
    """
    dzeta_dz, d2zeta_dz2 = dzeta_terms(z, L)
    return (dzeta_dz * dzeta_dz) * d2Ri_dzeta2 + d2zeta_dz2 * dRi_dzeta

def omit_error_metric(
    dzeta_dz: np.ndarray,
    d2zeta_dz2: np.ndarray,
    dRi_dzeta: np.ndarray,
    d2Ri_dzeta2: np.ndarray,
) -> np.ndarray:
    """
    E_omit = |(d²ζ/dz² * dRi/dζ) / ((dζ/dz)² * d²Ri/dζ²)|
    If E_omit < eps, constant-L shortcut is acceptable at that level.
    """
    denom = (dzeta_dz * dzeta_dz) * d2Ri_dzeta2
    # Avoid division by zero in flat-curvature zones
    safe = np.where(np.abs(denom) > 0, denom, np.nan)
    return np.abs(d2zeta_dz2 * dRi_dzeta) / np.abs(safe)

def map_curvature_auto(
    z: np.ndarray,
    L: np.ndarray,
    dRi_dzeta: np.ndarray,
    d2Ri_dzeta2: np.ndarray,
    eps: float = 0.05,
) -> np.ndarray:
    """
    Auto-select constant-L or full mapping per level using E_omit threshold.
    """
    dzeta_dz, d2zeta_dz2 = dzeta_terms(z, L)
    E = omit_error_metric(dzeta_dz, d2zeta_dz2, dRi_dzeta, d2Ri_dzeta2)
    # constant-L curvature = (1/L²) d²Ri/dζ²
    curv_const = d2Ri_dzeta2 / (L * L)
    curv_full = (dzeta_dz * dzeta_dz) * d2Ri_dzeta2 + d2zeta_dz2 * dRi_dzeta
    use_const = (E < eps) | ~np.isfinite(E)
    return np.where(use_const, curv_const, curv_full)

def map_curvature_point(
    z: float,
    L: float,
    dL: float,
    d2L: float,
    dRi_dzeta: float,
    d2Ri_dzeta2: float,
) -> float:
    """
    Pointwise mapping when L, dL, d2L are already known at z.
    """
    dzeta_dz = (L - z * dL) / (L * L)
    d2zeta_dz2 = (-2.0 * dL) / (L * L) - (z * d2L) / (L * L) + 2.0 * z * (dL * dL) / (L * L * L)
    return (dzeta_dz * dzeta_dz) * d2Ri_dzeta2 + d2zeta_dz2 * dRi_dzeta
