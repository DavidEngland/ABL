"""Profile generation and Ri↔ζ transformation utilities for MOST stability functions.

This module provides:
- A registry (make_profile) of commonly used atmospheric surface-layer (ABL / MOST) stability
  function families (momentum φ_m and heat φ_h) with parameter-driven callables.
- Conversions between gradient Richardson number Ri_g and nondimensional height ζ = z/L via
  near-neutral series + Newton refinement.
- Generic curvature helper (in companion files) relying on logarithmic derivatives.

Profile families (refs):
  BD_PL      : Power-law Businger–Dyer form near neutral (Businger et al. 1971; Dyer 1974).
  BD_CLASSIC : Piecewise (unstable power-law, stable linear) classical composite.
  HOG88      : Linear stable functions (Högström 1988).
  QSBL       : Quadratic truncated stable polynomial (cf. Beljaars & Holtslag 1991 style surrogate).
  CB         : Cheng–Brutsaert monotone (Cheng & Brutsaert 2005) |γ ζ| form.
  RPL        : Regularized power law (pole removed by rational mapping).
  VEXP       : Variable exponent (curvature tuning) extension.
  DTP        : Dynamic turbulent Prandtl (φ_h = Pr_t(Ri)*φ_m) (e.g., Grachev & Fairall 1997 concept).
  URC        : Direct Ri-based closure f_m(Ri) without explicit ζ (analytic inversion family).

Key equations:
  MOST power-law baseline: φ(ζ) = (1 - β ζ)^(-α) with domain ζ < 1/β (stable branch).
  Gradient Richardson number: Ri_g(ζ) = ζ φ_h / φ_m^2.
  Near-neutral series (ζ→0): Ri_g = ζ + Δ ζ^2 + 0.5(Δ^2 + c1) ζ^3 + …,
      Δ = α_h β_h - 2 α_m β_m,  c1 = α_h β_h^2 - 2 α_m β_m^2.

Inversion strategy:
  1. Series seed ζ₀ = Ri - Δ Ri^2 + (1.5Δ^2 - 0.5 c1) Ri^3 (valid for small Ri_g).
  2. Newton iteration on f(ζ)=ζ F(ζ) - Ri_target with derivative f'(ζ)=F + ζ F' = F(1+ζ V_log).
     Stop when |Δζ| < tol.

Numerical notes:
  - Central differences: O(h^2) truncation; using h=1e-6 balances roundoff vs. truncation for moderate ζ.
  - Instability near poles (ζ→1/β): denominator (1 - β ζ) amplifies derivative noise; guard sequences
    upstream (outside this module) should mask such ζ.
  - For URC profiles (Ri-based), φ_m, φ_h returned as functions of Ri directly; ζ-transformation omitted.
  - For DTP profiles, dynamic Pr_t applied in ζ→Ri wrappers, not in base φ generation.

References (short list):
  Businger et al. 1971; Dyer 1974; Högström 1988; Beljaars & Holtslag 1991; Cheng & Brutsaert 2005;
  Grachev & Fairall 1997; Monin & Obukhov 1954 (original similarity theory).
"""

import math
from typing import Callable, Dict, Tuple, Optional, Union

# --- Numerical derivative helpers -------------------------------------------------
def _central_diff(f: Callable[[float], float], x: float, h: float = 1e-6) -> float:
    """Second-order central difference approximation of f'(x).

    Error: O(h^2). Choose h small enough to limit truncation yet avoid cancellation.
    Not vectorized; intended for scalar ζ evaluations inside Newton refinement.
    """
    return (f(x + h) - f(x - h)) / (2 * h)

def _second_diff(f: Callable[[float], float], x: float, h: float = 1e-6) -> float:
    """Second-order central difference approximation of f''(x).

    Used only for auxiliary logging (W_log) in inversion diagnostics; not required for Newton step.
    Error: O(h^2).
    """
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h)

# --- Profile registry -------------------------------------------------------------
def make_profile(tag: str, pars: Dict[str, Union[float, int]]) -> Tuple[Callable[[float], float], Callable[[float], float]]:
    """Return (phi_m(ζ), phi_h(ζ)) callables for given profile family tag and parameter dict.

    Parameters
    ----------
    tag : str
        Registry key (case-insensitive). See module docstring for families.
    pars : dict
        Parameter set required by each family:
          BD_PL:     am, bm, ah, bh (α_m, β_m, α_h, β_h)
          BD_CLASSIC:a (unstable β), cm, ch, pm_exp, ph_exp
          HOG88:     cm, ch, c0h
          QSBL:      am, bm, ah, bh (linear/quadratic coefficients)
          CB:        gm, pm, gh, ph (γ_m, p_m, γ_h, p_h)
          RPL:       alpha_m, beta_m, delta_m, alpha_h, beta_h, delta_h
          VEXP:      alpha_m, beta_m, eta_m, alpha_h, beta_h, eta_h
          DTP:       base_tag, base_pars, a1, a2 (base φ before dynamic Pr_t)
          URC:       b_m, Ri_c, e_m [, b_h, e_h, Ri_c_h]

    Returns
    -------
    (phi_m, phi_h) : tuple of callables
        Each callable expects ζ (except URC where Ri is expected; see ri_to_phi_wrappers).

    Notes
    -----
    - Domains: Power-law families require 1 - β ζ > 0 for stability branch.
    - For URC, φ_m (and optionally φ_h) provided as Ri-based closure functions directly.
    - No masking or validation performed here; caller should enforce ζ domain limits.
    """
    tag = tag.upper()

    if tag == 'BD_PL':  # Businger–Dyer near-neutral power-law
        am, bm, ah, bh = pars['am'], pars['bm'], pars['ah'], pars['bh']
        return (lambda z: (1 - bm * z) ** (-am),
                lambda z: (1 - bh * z) ** (-ah))

    if tag == 'BD_CLASSIC':  # Piecewise classical form
        a = pars.get('a', 16.0)
        cm = pars.get('cm', 5.0)
        ch = pars.get('ch', 7.0)
        pm_e = pars.get('pm_exp', -0.25)
        ph_e = pars.get('ph_exp', -0.5)

        def phi_m(z: float) -> float:
            return (1 - a * z) ** (pm_e) if z < 0.0 else (1 + cm * z)

        def phi_h(z: float) -> float:
            return (1 - a * z) ** (ph_e) if z < 0.0 else (1 + ch * z)

        return phi_m, phi_h

    if tag == 'HOG88':  # Högström linear stable
        cm = pars.get('cm', 5.0)
        ch = pars.get('ch', 7.8)
        c0h = pars.get('c0h', 0.95)
        return (lambda z: 1.0 + cm * z,
                lambda z: c0h + ch * z)

    if tag == 'QSBL':  # Quadratic truncation stable surrogate
        am, bm, ah, bh = pars['am'], pars['bm'], pars['ah'], pars['bh']
        return (lambda z: 1.0 + am * z + bm * z * z,
                lambda z: 1.0 + ah * z + bh * z * z)

    if tag == 'CB':  # Cheng–Brutsaert monotone
        gm, pm, gh, ph = pars['gm'], pars['pm'], pars['gh'], pars['ph']
        return (lambda z: (1 + gm * abs(z)) ** pm,
                lambda z: (1 + gh * abs(z)) ** ph)

    if tag == 'RPL':  # Regularized power law
        am, bm, dm = pars['alpha_m'], pars['beta_m'], pars['delta_m']
        ah, bh, dh = pars['alpha_h'], pars['beta_h'], pars['delta_h']

        def g(b: float, d: float, z: float) -> float:
            return (b * z) / (1 + d * b * z)

        return (lambda z: (1 + g(bm, dm, z)) ** am,
                lambda z: (1 + g(bh, dh, z)) ** ah)

    if tag == 'VEXP':  # Variable exponent
        am, bm, em = pars['alpha_m'], pars['beta_m'], pars['eta_m']
        ah, bh, eh = pars['alpha_h'], pars['beta_h'], pars['eta_h']
        return (lambda z: (1 - bm * z) ** (-am * (1 + em * z)),
                lambda z: (1 - bh * z) ** (-ah * (1 + eh * z)))

    if tag == 'DTP':  # Return base profile; dynamic Pr_t applied later
        base_tag = pars['base_tag']
        base_pars = pars['base_pars']
        return make_profile(base_tag, base_pars)

    if tag == 'URC':  # Ri-based closure (Ri argument expected later)
        b_m, Ri_c, e_m = pars['b_m'], pars['Ri_c'], pars['e_m']
        fm = lambda Ri: (1 + b_m * Ri / Ri_c) ** (-e_m)
        if 'b_h' in pars and 'e_h' in pars:
            b_h, e_h = pars['b_h'], pars['e_h']
            Ri_c_h = pars.get('Ri_c_h', Ri_c)
            fh = lambda Ri: (1 + b_h * Ri / Ri_c_h) ** (-e_h)
        else:
            fh = lambda Ri: float('nan')  # heat unspecified
        return fm, fh

    raise ValueError(f"unknown profile tag '{tag}'")

# --- Ri / ζ utilities -------------------------------------------------------------
def F_from(phi_m: Callable[[float], float],
           phi_h: Callable[[float], float]) -> Callable[[float], float]:
    """Return F(ζ) = φ_h(ζ) / φ_m(ζ)^2 used in Ri_g(ζ) = ζ F(ζ)."""
    return lambda z: phi_h(z) / (phi_m(z) ** 2)

def ri_from_zeta(zeta: float,
                 phi_m: Callable[[float], float],
                 phi_h: Callable[[float], float]) -> float:
    """Compute Ri_g(ζ) given φ_m, φ_h.

    Ri_g(ζ) = ζ * φ_h / φ_m^2
    """
    F = F_from(phi_m, phi_h)
    return zeta * F(zeta)

def zeta_from_ri_series(Ri: float, Delta: float, c1: float) -> float:
    """Near-neutral inversion series ζ(Ri) truncated at O(Ri^3).

    ζ ≈ Ri - Δ Ri^2 + (1.5 Δ^2 - 0.5 c1) Ri^3

    Validity: Ri small (|Ri| ≪ 1/β). Used as initial guess for Newton refinement.

    Parameters
    ----------
    Ri    : gradient Richardson number (target)
    Delta : α_h β_h - 2 α_m β_m
    c1    : α_h β_h^2 - 2 α_m β_m^2
    """
    return Ri - Delta * Ri * Ri + (1.5 * Delta * Delta - 0.5 * c1) * (Ri ** 3)

def zeta_from_ri_newton(Ri_target: float,
                        phi_m: Callable[[float], float],
                        phi_h: Callable[[float], float],
                        z0: float,
                        tol: float = 1e-10,
                        maxit: int = 20) -> float:
    """Refine ζ solving f(ζ)=ζ F(ζ) - Ri_target = 0 by Newton's method.

    Derivative: f'(ζ) = F(ζ) + ζ F'(ζ) = F(ζ) (1 + ζ V_log),
      with V_log = (φ_h'/φ_h) - 2(φ_m'/φ_m).

    Numerical derivative for V_log uses central diff; F'(ζ) not formed explicitly.

    Parameters
    ----------
    Ri_target : desired gradient Richardson number
    phi_m, phi_h : profile callables
    z0        : initial guess (series seed)
    tol       : absolute update tolerance
    maxit     : iteration cap

    Returns
    -------
    ζ estimate. May be unreliable if near singular pole (1 - β ζ ≈ 0).

    Failure modes
    -------------
    - fp ≈ 0: derivative vanishes (extreme curvature); returns last iterate.
    - Divergence: does not check bracket; caller should validate domain.
    """
    F = F_from(phi_m, phi_h)
    z = z0
    for _ in range(maxit):
        # V_log via log(F) derivative; W_log is not needed for Newton step.
        Vlog = _central_diff(lambda zz: math.log(F(zz)), z)
        f = z * F(z) - Ri_target
        fp = F(z) + z * F(z) * Vlog  # F(1 + z V_log)
        if fp == 0.0:
            break
        dz = f / fp
        z -= dz
        if abs(dz) < tol:
            break
    return z

def ri_to_phi_wrappers(tag: str,
                       pars: Dict[str, Union[float, int]],
                       Delta: Optional[float] = None,
                       c1: Optional[float] = None
                       ) -> Tuple[Callable[[float], float], Callable[[float], float]]:
    """Return (f_m(Ri), f_h(Ri)) closures from ζ-space profile.

    If tag == 'URC' already Ri-based, simply returns its (fm, fh).

    For other profiles:
      1. Build φ_m(ζ), φ_h(ζ).
      2. Invert Ri → ζ using series + Newton (zeta_from_ri_series + zeta_from_ri_newton).
      3. Compose φ_m(ζ(Ri)), φ_h(ζ(Ri)).

    Dynamic turbulent Prandtl (DTP):
      - Base profile constructed; heat closure uses Pr_t(Ri) = 1 + a1 Ri + a2 Ri^2.

    Parameters
    ----------
    tag   : profile family
    pars  : parameter dict (see make_profile)
    Delta : neutral curvature coefficient (optional; improves series seed)
    c1    : second neutral coefficient (optional)

    Returns
    -------
    fm, fh : callables of Ri. fh may be NaN if heat form undefined (URC without heat).

    Notes
    -----
    - Series seed quality depends on supplied Delta, c1; if absent, fallback ζ₀ = Ri.
    - For large Ri (outside near-neutral), Newton may require damping (not implemented here).
    """
    if tag.upper() == 'URC':
        return make_profile(tag, pars)  # Already Ri-based

    phi_m, phi_h = make_profile(tag, pars)

    def zeta_of_Ri(Ri: float) -> float:
        z0 = zeta_from_ri_series(Ri, Delta, c1) if (Delta is not None and c1 is not None) else Ri
        return zeta_from_ri_newton(Ri, phi_m, phi_h, z0)

    if tag.upper() == 'DTP':
        base_tag = pars['base_tag']
        base_pars = pars['base_pars']
        a1 = pars.get('a1', 0.0)
        a2 = pars.get('a2', 0.0)
        base_m, _ = make_profile(base_tag, base_pars)

        def fm(Ri: float) -> float:
            z = zeta_of_Ri(Ri)
            return base_m(z)

        def fh(Ri: float) -> float:
            z = zeta_of_Ri(Ri)
            Pr_t = 1 + a1 * Ri + a2 * Ri * Ri
            return Pr_t * base_m(z)

        return fm, fh

    def fm(Ri: float) -> float:
        z = zeta_of_Ri(Ri)
        return phi_m(z)

    def fh(Ri: float) -> float:
        z = zeta_of_Ri(Ri)
        return phi_h(z)

    return fm, fh

def compute_vlog_wlog(zeta: float,
                      phi_m: Callable[[float], float],
                      phi_h: Callable[[float], float],
                      h: float = 1e-6) -> tuple[float, float, float, float]:
    """Compute per-component logs v_m, v_h and combined V_log, W_log at ζ.
    v_m = (φ_m'/φ_m); v_h = (φ_h'/φ_h); V_log = v_h - 2 v_m; W_log = V_log'.
    Central differences used if analytic derivatives unavailable.
    """
    pm = phi_m(zeta); ph = phi_h(zeta)
    v_m = _central_diff(phi_m, zeta, h) / max(pm, 1e-300)
    v_h = _central_diff(phi_h, zeta, h) / max(ph, 1e-300)
    # V_log derivative: differentiate V_log directly
    def V_of(zz: float) -> float:
        pmt = phi_m(zz); pht = phi_h(zz)
        vm = _central_diff(phi_m, zz, h) / max(pmt, 1e-300)
        vh = _central_diff(phi_h, zz, h) / max(pht, 1e-300)
        return vh - 2*vm
    V_log = v_h - 2*v_m
    W_log = _central_diff(V_of, zeta, h)
    return v_m, v_h, V_log, W_log

def rig_derivatives_zeta(
    zeta: float,
    phi_m: callable,
    phi_h: callable,
) -> tuple[float, float, float, float]:
    """
    Return (dRi_g/dζ, d²Ri_g/dζ², F, V_log) at ζ for given φ_m, φ_h.
    Uses unified logarithmic derivatives:
      F = φ_h / φ_m²
      V_log = v_h - 2 v_m
      W_log = dV_log/dζ
      dRi/dζ = F (1 + ζ V_log)
      d²Ri/dζ² = F [ 2 V_log + ζ (V_log² - W_log) ]
    """
    F_fun = F_from(phi_m, phi_h)
    F_val = F_fun(zeta)
    # compute combined logs
    _, _, V_log, W_log = compute_vlog_wlog(zeta, phi_m, phi_h)
    dRi_dzeta = F_val * (1.0 + zeta * V_log)
    d2Ri_dzeta2 = F_val * (2.0 * V_log + zeta * (V_log * V_log - W_log))
    return dRi_dzeta, d2Ri_dzeta2, F_val, V_log