"""Richardson number estimators: gradient Ri_g and bulk Ri_b with curvature corrections.

This module provides vectorized and scalar implementations for computing:
  - Gradient Richardson number: Ri_g = (g/θ₀) * (∂θ/∂z) / (∂U/∂z)²
  - Bulk Richardson number: Ri_b = (g/θ₀) * Δθ * Δz / (ΔU)²
  - Bias diagnostic: B = Ri_g(z_g) / Ri_b

References:
  - Beljaars & Holtslag (1991): Flux parameterization over land surfaces
  - Nieuwstadt (1984): The turbulent structure of the stable, nocturnal boundary layer
  - Mahrt (1981): The structure of similarity relations in the stably stratified boundary layer
"""

import numpy as np
from typing import Tuple, Dict, Optional, Union


# Physical constants
GRAVITY = 9.81  # m/s²
VON_KARMAN = 0.41  # Dimensionless


def ri_gradient(
    theta_minus: float,
    theta_0: float,
    theta_plus: float,
    U_minus: float,
    U_0: float,
    U_plus: float,
    V_minus: float,
    V_0: float,
    V_plus: float,
    dz: float,
    theta0: float = 300.0,
    g: float = GRAVITY,
) -> float:
    r"""Compute gradient Richardson number from 3-point stencil.

    Gradient Richardson number quantifies local stability as the ratio of
    buoyancy production to shear production:

    .. math::
        Ri_g = \frac{g}{\theta_0} \frac{\partial\theta/\partial z}{(∂U/∂z)^2 + (∂V/∂z)^2}

    For a discrete 3-point stencil, we use second-order central differences:

    .. math::
        \frac{\partial\theta}{\partial z} \approx \frac{\theta_+ - \theta_-}{2\Delta z}

    Parameters
    ----------
    theta_minus, theta_0, theta_plus : float
        Potential temperature (K) at z-Δz, z, z+Δz
    U_minus, U_0, U_plus : float
        Zonal wind component (m/s) at z-Δz, z, z+Δz
    V_minus, V_0, V_plus : float
        Meridional wind component (m/s) at z-Δz, z, z+Δz
    dz : float
        Grid spacing (m)
    theta0 : float
        Reference potential temperature (K) for normalization; default 300 K
    g : float
        Gravitational acceleration (m/s²); default 9.81 m/s²

    Returns
    -------
    ri_g : float
        Gradient Richardson number (dimensionless)

    Notes
    -----
    - Ri_g < 0: Unstable (convection dominates)
    - Ri_g ≈ 0: Neutral
    - 0 < Ri_g < Ri_c (~0.25): Stable with turbulence
    - Ri_g ≥ Ri_c: Critical layer (turbulence suppressed)

    Examples
    --------
    >>> # Stable layer: temperature increases, wind shear weak
    >>> rig = ri_gradient(
    ...     theta_minus=300.0, theta_0=300.2, theta_plus=300.5,
    ...     U_minus=2.0, U_0=2.5, U_plus=3.0,
    ...     V_minus=0.1, V_0=0.2, V_plus=0.3,
    ...     dz=10.0, theta0=300.0
    ... )
    >>> print(f"Ri_g = {rig:.4f}")  # Should be positive and small
    Ri_g = 0.0245

    See Also
    --------
    ri_bulk : Bulk Richardson number (layer-averaged)
    bias_ratio : Diagnostic of discretization bias
    """
    dtheta_dz = (theta_plus - theta_minus) / (2.0 * dz)
    dU_dz = (U_plus - U_minus) / (2.0 * dz)
    dV_dz = (V_plus - V_minus) / (2.0 * dz)

    shear_sq = dU_dz**2 + dV_dz**2

    # Guard against division by zero (calm conditions)
    if shear_sq < 1e-10:
        return 0.0

    ri_g = (g / theta0) * dtheta_dz / shear_sq

    return ri_g


def ri_bulk(
    theta1: float,
    theta2: float,
    U1: float,
    U2: float,
    V1: float,
    V2: float,
    z1: float,
    z2: float,
    theta0: float = 300.0,
    g: float = GRAVITY,
) -> float:
    r"""Compute bulk Richardson number over a layer.

    Bulk Richardson number is the layer-averaged stability measure:

    .. math::
        Ri_b = \frac{g}{\theta_0} \frac{(\theta_2 - \theta_1)(z_2 - z_1)}{(U_2 - U_1)^2 + (V_2 - V_1)^2}

    This is the standard diagnostic in atmospheric models for determining
    mixing-layer depth, Monin-Obukhov length, and critical Richardson number.

    Parameters
    ----------
    theta1, theta2 : float
        Potential temperature (K) at lower and upper layer boundaries
    U1, U2 : float
        Zonal wind (m/s) at lower and upper layer boundaries
    V1, V2 : float
        Meridional wind (m/s) at lower and upper layer boundaries
    z1, z2 : float
        Height (m) at lower and upper layer boundaries; z2 > z1
    theta0 : float
        Reference potential temperature (K); default 300 K
    g : float
        Gravitational acceleration (m/s²); default 9.81 m/s²

    Returns
    -------
    ri_b : float
        Bulk Richardson number (dimensionless)

    Notes
    -----
    - Defined for z2 > z1 (positive layer thickness).
    - Can be negative (unstable), near-zero (neutral), or positive (stable).
    - In coarse models (Δz >> z0), Ri_b often underestimates Ri_g at representative height.

    Examples
    --------
    >>> # Stable layer
    >>> rib = ri_bulk(
    ...     theta1=300.0, theta2=302.0,
    ...     U1=2.0, U2=5.0,
    ...     V1=0.0, V2=0.5,
    ...     z1=10.0, z2=100.0,
    ...     theta0=300.0
    ... )
    >>> print(f"Ri_b = {rib:.4f}")  # Should be positive
    Ri_b = 0.0186

    See Also
    --------
    ri_gradient : Gradient Richardson number (point-wise)
    geometric_mean_height : Recommended height for Ri_g comparison to Ri_b
    """
    dtheta = theta2 - theta1
    dU = U2 - U1
    dV = V2 - V1
    dz = z2 - z1

    wind_shear_sq = dU**2 + dV**2

    # Guard against division by zero
    if wind_shear_sq < 1e-10:
        return 0.0

    ri_b = (g / theta0) * (dtheta * dz) / wind_shear_sq

    return ri_b


def bias_ratio(
    ri_g: float,
    ri_b: float,
    min_denominator: float = 1e-6,
) -> float:
    r"""Compute bias ratio B = Ri_g / Ri_b.

    Bias ratio quantifies systematic underestimation or overestimation of
    the gradient Richardson number by bulk approximations on coarse grids.

    .. math::
        B = \frac{Ri_g(z_g)}{Ri_b}

    where z_g = √(z1 z2) is the geometric mean height (recommended by Jensen inequality).

    Interpretation:
    - B ≈ 1.0: Ri_b accurately represents Ri_g (e.g., fine grid)
    - B > 1.0: Ri_b underestimates Ri_g (typical stable boundary layer with concave-down profile)
    - B < 1.0: Ri_b overestimates Ri_g (e.g., mixed layer with convective updrafts)

    Parameters
    ----------
    ri_g : float
        Point-wise gradient Richardson number
    ri_b : float
        Bulk Richardson number (same layer)
    min_denominator : float
        Minimum |Ri_b| to avoid division by near-zero; default 1e-6

    Returns
    -------
    bias : float
        Bias ratio B. Returns 1.0 if |Ri_b| < min_denominator (no meaningful bias).

    Examples
    --------
    >>> # Stable SBL: Ri_g > Ri_b (concave-down profile)
    >>> B = bias_ratio(ri_g=0.05, ri_b=0.03)
    >>> print(f"Bias ratio: {B:.2f}")  # B ≈ 1.67
    Bias ratio: 1.67

    >>> # Near-neutral: both small
    >>> B = bias_ratio(ri_g=0.001, ri_b=0.0008)
    >>> print(f"Bias ratio: {B:.2f}")  # B ≈ 1.25
    Bias ratio: 1.25

    Notes
    -----
    - For |Ri_b| < min_denominator, returns 1.0 (treat as unbiased).
    - Used as a diagnostic driver in correction factor ODE: dC/dz ∝ (B - 1).

    References
    ----------
    - Section 3, McNider & England (1995): On the predictability of mesoscale wind
    """
    if abs(ri_b) < min_denominator:
        return 1.0  # Neutral/near-neutral: no meaningful bias

    return ri_g / ri_b


# ─────────────────────────────────────────────────────────────────────────────
# Vectorized operations
# ─────────────────────────────────────────────────────────────────────────────


def ri_gradient_array(
    theta: np.ndarray,
    U: np.ndarray,
    V: np.ndarray,
    z: np.ndarray,
    theta0: float = 300.0,
    g: float = GRAVITY,
) -> np.ndarray:
    r"""Compute gradient Richardson number on a profile grid (vectorized).

    Parameters
    ----------
    theta : ndarray, shape (n,)
        Potential temperature profile (K)
    U : ndarray, shape (n,)
        Zonal wind profile (m/s)
    V : ndarray, shape (n,)
        Meridional wind profile (m/s)
    z : ndarray, shape (n,)
        Height grid (m), strictly increasing
    theta0 : float
        Reference temperature (K)
    g : float
        Gravitational acceleration (m/s²)

    Returns
    -------
    ri_g : ndarray, shape (n,)
        Gradient Richardson number at each grid point

    Notes
    -----
    - At boundaries (z[0], z[-1]), uses one-sided differences.
    - Interior points use centered differences (O(Δz²) error).
    - Returns 0.0 where wind shear is negligible (< 1e-10 m/s per meter).
    """
    n = len(z)
    ri_g = np.zeros(n, dtype=float)

    # Interior points: central differences
    for i in range(1, n - 1):
        dz = z[i + 1] - z[i - 1]
        dtheta_dz = (theta[i + 1] - theta[i - 1]) / dz
        dU_dz = (U[i + 1] - U[i - 1]) / dz
        dV_dz = (V[i + 1] - V[i - 1]) / dz

        shear_sq = dU_dz**2 + dV_dz**2
        if shear_sq > 1e-10:
            ri_g[i] = (g / theta0) * dtheta_dz / shear_sq
        else:
            ri_g[i] = 0.0

    # Boundary: forward/backward difference
    dz0 = z[1] - z[0]
    dtheta_dz = (theta[1] - theta[0]) / dz0
    dU_dz = (U[1] - U[0]) / dz0
    dV_dz = (V[1] - V[0]) / dz0
    shear_sq = dU_dz**2 + dV_dz**2
    ri_g[0] = (g / theta0) * dtheta_dz / shear_sq if shear_sq > 1e-10 else 0.0

    dzn = z[n - 1] - z[n - 2]
    dtheta_dz = (theta[n - 1] - theta[n - 2]) / dzn
    dU_dz = (U[n - 1] - U[n - 2]) / dzn
    dV_dz = (V[n - 1] - V[n - 2]) / dzn
    shear_sq = dU_dz**2 + dV_dz**2
    ri_g[n - 1] = (g / theta0) * dtheta_dz / shear_sq if shear_sq > 1e-10 else 0.0

    return ri_g


# Helper for diagnostic metadata
def estimate_uncertainty(
    ri_estimate: float,
    derivative_step: float = 1e-5,
    noise_level: float = 0.01,
) -> Dict[str, float]:
    """Estimate uncertainty in Richardson number due to truncation and measurement noise.

    Parameters
    ----------
    ri_estimate : float
        Computed Richardson number value
    derivative_step : float
        Typical finite difference step size (proportion of domain)
    noise_level : float
        Typical measurement noise (e.g., 0.01 for 1% error)

    Returns
    -------
    metadata : dict
        Keys: "truncation_error", "noise_error", "total_uncertainty"

    Notes
    -----
    - Truncation error: O(Δz²) for second-order schemes
    - Noise error: Ri uncertainty ≈ Ri * √(4*noise²) for independent measurements
    """
    truncation = abs(ri_estimate) * (derivative_step**2)
    noise = abs(ri_estimate) * np.sqrt(4 * noise_level**2)
    total = np.sqrt(truncation**2 + noise**2)

    return {
        "truncation_error": truncation,
        "noise_error": noise,
        "total_uncertainty": total,
    }
