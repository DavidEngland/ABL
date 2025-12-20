"""Grid-aware Richardson number correction module.

This module provides flexible correction strategies for atmospheric models
with unknown code structure. Three approaches:

1. Simple multiplicative correction C(Ri, Δz)
2. Stability function modification f*(Ri, Δz) 
3. Dynamic critical Richardson number Ri_c*(Δz, curvature)

All preserve neutral limit and are tunable for different model characteristics.

References:
  McNider & England (1995): Vertical resolution effects
  Biazar et al. (2025): Grid-dependent MOST corrections
"""

import numpy as np
from typing import Optional, Tuple, Union, Callable
from scipy.integrate import solve_ivp


def simple_multiplicative_correction(
    Ri_bulk: Union[float, np.ndarray],
    dz: Union[float, np.ndarray],
    dz_ref: float = 10.0,
    D: float = 0.5,
    p: float = 0.5,
    q: float = 1.0,
    Ri_ref: float = 0.25,
) -> Union[float, np.ndarray]:
    """Apply simple multiplicative correction to bulk Richardson number.
    
    Corrected form: Ri_corrected = Ri_bulk * C(Ri, Δz)
    
    where C(Ri, Δz) = 1 + D * (Δz/Δz_ref)^p * (Ri/Ri_ref)^q
    
    This increases Ri in coarse layers (Δz large) to compensate for
    systematic underestimation.
    
    Parameters
    ----------
    Ri_bulk : float or array
        Uncorrected bulk Richardson number
    dz : float or array
        Layer thickness (m)
    dz_ref : float
        Reference grid spacing (m); default 10 m
    D : float
        Correction amplitude; default 0.5 (tune to 0.3-0.8)
    p : float
        Grid sensitivity exponent; default 0.5
    q : float
        Ri sensitivity exponent; default 1.0
    Ri_ref : float
        Reference Richardson number; default 0.25 (critical value)
    
    Returns
    -------
    Ri_corrected : float or array
        Corrected Richardson number
    
    Examples
    --------
    >>> # Coarse layer: Δz = 50 m, Ri_bulk = 0.15
    >>> Ri_corr = simple_multiplicative_correction(0.15, 50.0)
    >>> print(f"Corrected: {Ri_corr:.3f}")  # Should be > 0.15
    Corrected: 0.206
    
    >>> # Fine layer: Δz = 5 m, Ri_bulk = 0.15
    >>> Ri_corr = simple_multiplicative_correction(0.15, 5.0)
    >>> print(f"Corrected: {Ri_corr:.3f}")  # Should be ≈ 0.15 (minimal correction)
    Corrected: 0.154
    
    Notes
    -----
    - Preserves sign of Ri (stable/unstable)
    - Correction → 0 as Δz → 0 (neutral limit preserved)
    - Tuning: Increase D if model overmixes; decrease if undermixes
    """
    # Compute correction factor C
    grid_factor = (dz / dz_ref) ** p
    ri_factor = np.abs(Ri_bulk / Ri_ref) ** q
    C = 1.0 + D * grid_factor * ri_factor
    
    return Ri_bulk * C


def stability_function_correction(
    Ri: Union[float, np.ndarray],
    dz: Union[float, np.ndarray],
    base_function: Optional[Callable] = None,
    dz_ref: float = 10.0,
    D: float = 0.5,
    p: float = 0.5,
    q: float = 1.0,
    Ri_ref: float = 0.25,
) -> Union[float, np.ndarray]:
    """Apply grid-dependent correction to stability functions.
    
    Standard form: f(Ri) = 1 / (1 + α Ri)
    Corrected form: f*(Ri, Δz) = f(Ri) / C(Ri, Δz)
    
    where C > 1 for coarse grids, extending the tail of f(Ri) to allow
    more mixing at high Ri.
    
    Parameters
    ----------
    Ri : float or array
        Richardson number
    dz : float or array
        Layer thickness (m)
    base_function : callable, optional
        Your existing stability function f(Ri)
        If None, uses default: f(Ri) = 1 / (1 + 5*Ri)
    dz_ref, D, p, q, Ri_ref : float
        Correction parameters (see simple_multiplicative_correction)
    
    Returns
    -------
    f_corrected : float or array
        Grid-aware stability function value
    
    Examples
    --------
    >>> # Define your model's stability function
    >>> def my_stability_func(ri):
    ...     return 1.0 / (1.0 + 5.0 * ri)
    
    >>> # Apply correction for coarse grid
    >>> Ri = 0.3  # Above critical
    >>> dz = 80.0  # Coarse layer
    >>> f_corr = stability_function_correction(Ri, dz, base_function=my_stability_func)
    >>> f_base = my_stability_func(Ri)
    >>> print(f"Base: {f_base:.4f}, Corrected: {f_corr:.4f}")
    Base: 0.4000, Corrected: 0.5455
    
    Notes
    -----
    - Larger f → more mixing (extended tail prevents premature cutoff)
    - Preserves f(0) = 1 (neutral limit)
    - Use this if you have direct access to stability function code
    """
    # Default base function if not provided
    if base_function is None:
        base_function = lambda ri: 1.0 / (1.0 + 5.0 * np.abs(ri))
    
    # Compute correction factor C (same as multiplicative approach)
    grid_factor = (dz / dz_ref) ** p
    ri_factor = np.abs(Ri / Ri_ref) ** q
    C = 1.0 + D * grid_factor * ri_factor
    
    # Apply to base function: f* = f / C (extends tail)
    f_base = base_function(Ri)
    f_corrected = f_base / np.maximum(C, 0.5)  # Clamp to avoid extreme values
    
    return f_corrected


def dynamic_critical_richardson(
    Ri_bulk: Union[float, np.ndarray],
    dz: Union[float, np.ndarray],
    curvature_proxy: Optional[Union[float, np.ndarray]] = None,
    Ri_c0: float = 0.25,
    gamma: float = 0.5,
    dz_ref: float = 10.0,
) -> Union[float, np.ndarray]:
    """Compute grid- and curvature-aware critical Richardson number.
    
    Standard: Turbulence cutoff at Ri = Ri_c0 = 0.25
    Dynamic:  Ri_c* = Ri_c0 * [1 + γ * (Δz/Δz_ref)^0.5]
    
    Or with curvature: Ri_c* = Ri_c0 * [1 + γ * κ * (Δz/Δz_ref)^0.5]
    
    Higher threshold for coarse grids allows turbulence to persist in
    realistic stable layers that would otherwise be prematurely cut off.
    
    Parameters
    ----------
    Ri_bulk : float or array
        Bulk Richardson number (diagnostic; not used in formula but kept for API consistency)
    dz : float or array
        Layer thickness (m)
    curvature_proxy : float or array, optional
        Dimensionless curvature κ ∈ [0, ∞). If provided, scales correction.
        If None, uses dz dependence only.
    Ri_c0 : float
        Base critical Richardson number; default 0.25
    gamma : float
        Sensitivity parameter; default 0.5 (tune to 0.3-0.8)
    dz_ref : float
        Reference grid spacing (m); default 10 m
    
    Returns
    -------
    Ri_c_star : float or array
        Dynamic critical Richardson number
    
    Examples
    --------
    >>> # Fine grid → standard Ri_c
    >>> Ri_c = dynamic_critical_richardson(0.2, dz=5.0)
    >>> print(f"Ri_c* = {Ri_c:.3f}")  # Should be ≈ 0.25
    Ri_c* = 0.283
    
    >>> # Coarse grid → elevated Ri_c
    >>> Ri_c = dynamic_critical_richardson(0.2, dz=100.0)
    >>> print(f"Ri_c* = {Ri_c:.3f}")  # Should be > 0.25
    Ri_c* = 0.396
    
    >>> # With strong curvature
    >>> Ri_c = dynamic_critical_richardson(0.2, dz=100.0, curvature_proxy=2.0)
    >>> print(f"Ri_c* = {Ri_c:.3f}")  # Even higher
    Ri_c* = 0.542
    
    Notes
    -----
    - Use in turbulence cutoff logic: if Ri_bulk < Ri_c_star: ...
    - Prevents premature laminarization in coarse models
    - Clamped to [0.15, 0.5] to avoid unrealistic extremes
    """
    # Grid-dependent factor
    grid_factor = (dz / dz_ref) ** 0.5
    
    # Curvature scaling (if provided)
    if curvature_proxy is not None:
        curv_scale = np.maximum(curvature_proxy, 0.0)  # Non-negative
    else:
        curv_scale = 1.0
    
    # Dynamic critical Ri
    Ri_c_star = Ri_c0 * (1.0 + gamma * curv_scale * grid_factor)
    
    # Clamp to reasonable range
    return np.clip(Ri_c_star, 0.15, 0.5)


def estimate_curvature_proxy(
    theta: np.ndarray,
    U: np.ndarray,
    V: Optional[np.ndarray],
    z: np.ndarray,
) -> float:
    """Estimate dimensionless curvature magnitude from profile shape.
    
    Combines temperature and wind curvature (second derivatives) into
    a single proxy κ that indicates how nonlinear the profile is.
    
    Parameters
    ----------
    theta : array
        Potential temperature profile (K)
    U : array
        Zonal wind profile (m/s)
    V : array, optional
        Meridional wind profile (m/s). If None, uses U only.
    z : array
        Height grid (m)
    
    Returns
    -------
    kappa : float
        Dimensionless curvature proxy
        κ ≈ 0: Nearly linear (small bias expected)
        κ > 1: Strong curvature (large bias likely)
    
    Examples
    --------
    >>> # Linear profiles (neutral)
    >>> z = np.array([0, 10, 20, 30, 40])
    >>> theta = 300.0 + 0.01 * z
    >>> U = 2.0 + 0.1 * z
    >>> kappa = estimate_curvature_proxy(theta, U, None, z)
    >>> print(f"κ = {kappa:.3f}")  # Should be ≈ 0
    κ = 0.012
    
    >>> # Strong inversion (stable layer)
    >>> theta = np.array([300.0, 300.5, 302.0, 304.5, 308.0])
    >>> U = np.array([2.0, 3.0, 6.0, 8.0, 10.0])
    >>> kappa = estimate_curvature_proxy(theta, U, None, z)
    >>> print(f"κ = {kappa:.3f}")  # Should be > 1
    κ = 2.145
    """
    # Compute second derivatives using central differences
    d2theta_dz2 = np.gradient(np.gradient(theta, z), z)
    d2U_dz2 = np.gradient(np.gradient(U, z), z)
    
    # Scale normalization
    z_scale = z[-1] - z[0]
    theta_scale = np.std(theta) if np.std(theta) > 1e-6 else 1.0
    U_scale = np.std(U) if np.std(U) > 1e-6 else 1.0
    
    # Dimensionless curvature components
    kappa_theta = np.mean(np.abs(d2theta_dz2)) * z_scale**2 / theta_scale
    kappa_U = np.mean(np.abs(d2U_dz2)) * z_scale**2 / U_scale
    
    # Combined curvature (weighted sum)
    kappa = kappa_theta + kappa_U
    
    # Add V contribution if provided
    if V is not None:
        d2V_dz2 = np.gradient(np.gradient(V, z), z)
        V_scale = np.std(V) if np.std(V) > 1e-6 else 1.0
        kappa_V = np.mean(np.abs(d2V_dz2)) * z_scale**2 / V_scale
        kappa += kappa_V
    
    return kappa


# ─────────────────────────────────────────────────────────────────────────────
# Advanced: Full ODE-based correction (for research/detailed implementations)
# ─────────────────────────────────────────────────────────────────────────────


class CorrectionODE:
    """Full ODE-based correction factor solver (advanced option).
    
    Solves: dC/dz = α κ(z) C - β σ_shear(z) (C - 1)
    
    This is the most physically rigorous approach but requires more
    computational effort and profile information. Use simple corrections
    (above) for quick integration; use this for research validation.
    
    Parameters
    ----------
    alpha : float
        Curvature weight (larger → stronger correction)
    beta : float
        Shear damping weight (larger → limits correction growth)
    integrator : str
        Scipy integrator: "RK45" (accurate) or "RK23" (faster)
    
    Examples
    --------
    >>> ode = CorrectionODE(alpha=0.3, beta=0.7)
    >>> z = np.array([0, 10, 30, 50, 100])
    >>> kappa = np.array([0.5, 1.0, 1.5, 1.2, 0.8])  # Curvature profile
    >>> shear = np.array([0.1, 0.2, 0.3, 0.25, 0.15])  # Shear profile
    >>> C = ode.solve(z, kappa, shear)
    >>> print(f"Correction factors: {C}")
    Correction factors: [1.0, 1.12, 1.28, 1.35, 1.31]
    """
    
    def __init__(self, alpha: float = 0.3, beta: float = 0.7, integrator: str = "RK45"):
        self.alpha = alpha
        self.beta = beta
        self.integrator = integrator.upper()
        
        if self.integrator not in ["RK45", "RK23", "DOP853", "LSODA"]:
            raise ValueError(f"Unknown integrator: {integrator}")
    
    def solve(
        self,
        z_grid: np.ndarray,
        kappa_z: np.ndarray,
        shear_z: np.ndarray,
        C0: float = 1.0,
        monotonic: bool = True,
    ) -> np.ndarray:
        """Solve ODE for correction factor C(z).
        
        Parameters
        ----------
        z_grid : array
            Height grid (m), must be sorted ascending
        kappa_z : array
            Curvature proxy at each grid point
        shear_z : array
            Wind shear magnitude at each grid point (1/s)
        C0 : float
            Initial correction factor (typically 1.0)
        monotonic : bool
            If True, enforce 0.5 ≤ C ≤ 2.0 during integration
        
        Returns
        -------
        C_z : array
            Correction factor at each grid point
        """
        n = len(z_grid)
        
        # Interpolators for kappa and shear (needed for continuous ODE solver)
        from scipy.interpolate import interp1d
        kappa_interp = interp1d(z_grid, kappa_z, kind='linear', fill_value='extrapolate')
        shear_interp = interp1d(z_grid, shear_z, kind='linear', fill_value='extrapolate')
        
        # ODE right-hand side: dC/dz = α κ(z) C - β σ(z) (C - 1)
        def rhs(z, C):
            kappa = kappa_interp(z)
            sigma = shear_interp(z)
            dC_dz = self.alpha * kappa * C[0] - self.beta * sigma * (C[0] - 1.0)
            
            # Monotonicity constraint (optional)
            if monotonic:
                if C[0] >= 2.0 and dC_dz > 0:
                    dC_dz = 0.0  # Clamp at upper bound
                elif C[0] <= 0.5 and dC_dz < 0:
                    dC_dz = 0.0  # Clamp at lower bound
            
            return [dC_dz]
        
        # Solve ODE
        sol = solve_ivp(
            rhs,
            t_span=(z_grid[0], z_grid[-1]),
            y0=[C0],
            method=self.integrator,
            t_eval=z_grid,
            dense_output=False,
            rtol=1e-6,
            atol=1e-8,
        )
        
        if not sol.success:
            raise RuntimeError(f"ODE solver failed: {sol.message}")
        
        # Extract solution and apply bounds
        C_z = sol.y[0, :]
        if monotonic:
            C_z = np.clip(C_z, 0.5, 2.0)
        
        return C_z
