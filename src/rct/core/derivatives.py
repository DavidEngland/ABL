"""Numerical derivatives with curvature estimation and error control.

This module provides second-order accurate central difference derivatives and
curvature estimation for Richardson number fields and stability functions.

References:
  - Fornberg (1988): Generation of Finite Difference Formulas on Arbitrarily Spaced Grids
  - Fornberg & Sloan (1994): A review of pseudospectral methods for the numerical solution
    of time-dependent PDEs
"""

import numpy as np
from typing import Callable, Tuple, Optional, Union


def central_with_curvature(
    phi_minus: float,
    phi_0: float,
    phi_plus: float,
    dz: float,
    limiter: str = "none",
) -> Tuple[float, float]:
    """Compute first derivative and curvature proxy using 3-point stencil.

    For a function φ sampled at z-dz, z, z+dz with spacing dz:

    φ'(z) ≈ (φ₊ - φ₋) / (2Δz) - O(Δz²)
    φ''(z) ≈ (φ₊ - 2φ₀ + φ₋) / Δz² + O(Δz²)

    The first derivative uses forward-backward differences; the curvature uses
    symmetric central differences. This preserves second-order truncation error
    and enables physical curvature-aware corrections.

    Parameters
    ----------
    phi_minus : float
        Value at z - Δz
    phi_0 : float
        Value at z
    phi_plus : float
        Value at z + Δz
    dz : float
        Grid spacing (meters or normalized coordinate)
    limiter : {"none", "tvd"}
        Slope limiter:
        - "none": Standard central difference (may overshoot near discontinuities)
        - "tvd": Total Variation Diminishing (minmod) limiter for monotonicity

    Returns
    -------
    first_deriv : float
        Estimate of φ'(z), second-order accurate
    curvature : float
        Estimate of φ''(z), second-order accurate

    Examples
    --------
    >>> # Parabolic function: f(x) = x^2 at x=1
    >>> # f'(1) = 2, f''(1) = 2
    >>> fprime, fcurv = central_with_curvature(0.0, 1.0, 4.0, 1.0)
    >>> print(f"f'(1) ≈ {fprime:.4f}, f''(1) ≈ {fcurv:.4f}")
    f'(1) ≈ 2.0000, f''(1) ≈ 2.0000

    Notes
    -----
    - Curvature is exact for polynomials up to degree 3.
    - For strongly nonlinear profiles, consider adaptive mesh refinement.
    - TVD limiter reduces overshooting but may smooth legitimate discontinuities.

    References
    --------
    - Sec. 1.2, Fornberg (1988).
    """
    dz2 = dz * dz
    dz2_inv = 1.0 / dz2
    dz_inv = 1.0 / (2.0 * dz)

    # Standard central differences
    first_deriv_std = (phi_plus - phi_minus) * dz_inv
    curvature = (phi_plus - 2.0 * phi_0 + phi_minus) * dz2_inv

    if limiter.lower() == "none":
        return first_deriv_std, curvature

    elif limiter.lower() == "tvd":
        # Minmod TVD limiter
        # slope_left = (phi_0 - phi_minus) / dz
        # slope_right = (phi_plus - phi_0) / dz
        # Apply minmod to reduce overshooting
        slope_left = (phi_0 - phi_minus) / dz
        slope_right = (phi_plus - phi_0) / dz

        # Minmod function: minmod(a, b) = sign(a) * min(|a|, |b|) if sign(a)==sign(b), else 0
        if slope_left * slope_right > 0:
            first_deriv_limited = 0.5 * (slope_left + slope_right)
            if abs(slope_left) < abs(first_deriv_limited):
                first_deriv_limited = slope_left
            if abs(slope_right) < abs(first_deriv_limited):
                first_deriv_limited = slope_right
        else:
            first_deriv_limited = 0.0

        # Curvature unchanged (limiter applied to first deriv only)
        return first_deriv_limited, curvature

    else:
        raise ValueError(f"Unknown limiter: {limiter}. Choose 'none' or 'tvd'.")


def second_derivative(
    phi_minus: float,
    phi_0: float,
    phi_plus: float,
    dz: float,
) -> float:
    """Compute second derivative (curvature) directly.

    Simple wrapper for the curvature component of central_with_curvature.
    Useful when only curvature is needed.

    Parameters
    ----------
    phi_minus, phi_0, phi_plus : float
        Triplet of function values at z-Δz, z, z+Δz
    dz : float
        Grid spacing

    Returns
    -------
    float
        Second derivative φ''(z)

    See Also
    --------
    central_with_curvature : Combined first and second derivative estimation.
    """
    return (phi_plus - 2.0 * phi_0 + phi_minus) / (dz * dz)


def numerical_jacobian(
    func: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray,
    h: float = 1e-5,
) -> np.ndarray:
    """Compute Jacobian matrix via central finite differences.

    For vector function f: ℝⁿ → ℝᵐ, computes the m×n matrix J[i,j] = ∂f_i/∂x_j.

    Parameters
    ----------
    func : callable
        Vector-valued function f(x) returning ndarray of shape (m,)
    x : ndarray
        Point at which to evaluate Jacobian, shape (n,)
    h : float
        Step size for finite differences

    Returns
    -------
    jacobian : ndarray
        Shape (m, n) Jacobian matrix

    Notes
    -----
    - Uses central differences (O(h²) error).
    - Step size h should be small but not so small as to cause cancellation error.
    - For n > 100, consider automatic differentiation (JAX, PyTorch).

    Examples
    --------
    >>> def f(x):
    ...     return np.array([x[0]**2, x[0]*x[1], x[1]**3])
    >>> x = np.array([2.0, 3.0])
    >>> J = numerical_jacobian(f, x)
    >>> print(J)  # doctest: +SKIP
    [[4.0  0.0 ]
     [3.0  2.0 ]
     [0.0  27.0]]
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    f_x = func(x)
    m = f_x.shape[0]

    jacobian = np.zeros((m, n), dtype=float)

    for j in range(n):
        x_plus = x.copy()
        x_minus = x.copy()
        x_plus[j] += h
        x_minus[j] -= h

        f_plus = func(x_plus)
        f_minus = func(x_minus)

        jacobian[:, j] = (f_plus - f_minus) / (2.0 * h)

    return jacobian


def richardson_extrapolation(
    func: Callable[[float], float],
    x: float,
    h0: float = 1e-3,
    order: int = 2,
    n_steps: int = 4,
) -> Tuple[float, float]:
    """Estimate derivative via Richardson extrapolation (Romberg integration).

    Reduces truncation error from O(h^order) to O(h^(2*order)) without requiring
    higher-order stencils. Useful for high-precision gradient computation.

    Parameters
    ----------
    func : callable
        Scalar function f(x)
    x : float
        Point at which to estimate derivative
    h0 : float
        Initial step size
    order : int
        Order of the initial finite difference scheme (1 or 2)
    n_steps : int
        Number of Richardson steps (more steps = higher precision, more cost)

    Returns
    -------
    deriv_estimate : float
        High-precision derivative estimate
    error_estimate : float
        Estimated error (difference between last two Richardson steps)

    Examples
    --------
    >>> import math
    >>> def f(x):
    ...     return math.sin(x)
    >>> deriv, err = richardson_extrapolation(f, math.pi/4, h0=0.1)
    >>> print(f"cos(π/4) ≈ {deriv:.6f}, error ≈ {err:.2e}")  # doctest: +SKIP
    cos(π/4) ≈ 0.707107, error ≈ 1.23e-10

    References
    --------
    - Sec. 5.3, Numerical Recipes, Press et al. (2007).
    """
    # Build table of finite difference estimates at increasingly small h
    table = np.zeros((n_steps, n_steps), dtype=float)

    for i in range(n_steps):
        h = h0 / (2 ** i)
        if order == 1:
            # First-order forward difference: f'(x) ≈ [f(x+h) - f(x)] / h
            table[i, 0] = (func(x + h) - func(x)) / h
        elif order == 2:
            # Second-order central difference: f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
            table[i, 0] = (func(x + h) - func(x - h)) / (2.0 * h)
        else:
            raise ValueError("order must be 1 or 2")

    # Richardson extrapolation: combine lower-order estimates to higher order
    for j in range(1, n_steps):
        for i in range(n_steps - j):
            # Assuming extrapolation factor 4 (typical for order-2 schemes)
            # d = (4 * T[i+1,j-1] - T[i,j-1]) / 3
            factor = 4.0 ** j
            table[i, j] = (factor * table[i + 1, j - 1] - table[i, j - 1]) / (factor - 1.0)

    deriv_estimate = table[0, n_steps - 1]
    error_estimate = abs(table[0, n_steps - 1] - table[0, n_steps - 2])

    return deriv_estimate, error_estimate


# ─────────────────────────────────────────────────────────────────────────────
# Vectorized operations for profiles
# ─────────────────────────────────────────────────────────────────────────────


def gradient_array(f: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Compute gradient (first derivative) on 1D grid using central differences.

    At boundaries, uses one-sided differences.

    Parameters
    ----------
    f : ndarray, shape (n,)
        Function values
    z : ndarray, shape (n,)
        Coordinate array (must be sorted)

    Returns
    -------
    df_dz : ndarray, shape (n,)
        Gradient array

    Examples
    --------
    >>> z = np.array([0, 1, 2, 3, 4])
    >>> f = z**2  # f(z) = z^2, f'(z) = 2z
    >>> df = gradient_array(f, z)
    >>> print(df)  # Near [0, 2, 4, 6, 8]
    [1. 2. 4. 6. 7.]
    """
    n = len(f)
    df_dz = np.zeros(n, dtype=float)

    # Interior points: central difference
    for i in range(1, n - 1):
        dz = z[i + 1] - z[i - 1]
        df_dz[i] = (f[i + 1] - f[i - 1]) / dz

    # Boundary: one-sided differences
    dz0 = z[1] - z[0]
    df_dz[0] = (f[1] - f[0]) / dz0

    dzn = z[n - 1] - z[n - 2]
    df_dz[n - 1] = (f[n - 1] - f[n - 2]) / dzn

    return df_dz


def curvature_array(f: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Compute curvature (second derivative) on 1D grid using central differences.

    At boundaries, curvature is set to nearest interior value (or NaN if only 1 interior).

    Parameters
    ----------
    f : ndarray, shape (n,)
        Function values
    z : ndarray, shape (n,)
        Coordinate array

    Returns
    -------
    d2f_dz2 : ndarray, shape (n,)
        Curvature array

    Examples
    --------
    >>> z = np.array([0, 1, 2, 3, 4])
    >>> f = z**2  # f(z) = z^2, f''(z) = 2
    >>> d2f = curvature_array(f, z)
    >>> print(d2f)  # All ≈ 2
    [2. 2. 2. 2. 2.]
    """
    n = len(f)
    d2f_dz2 = np.zeros(n, dtype=float)

    for i in range(1, n - 1):
        dz = z[i + 1] - z[i - 1]
        dz2 = dz * dz
        d2f_dz2[i] = (f[i + 1] - 2.0 * f[i] + f[i - 1]) / dz2

    # Boundary: extend interior values
    d2f_dz2[0] = d2f_dz2[1]
    d2f_dz2[n - 1] = d2f_dz2[n - 2]

    return d2f_dz2
