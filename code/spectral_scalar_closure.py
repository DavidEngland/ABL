"""Spectral scalar closures for ABL similarity functions.

This module turns the Gegenbauer/Legendre analysis into an executable model.

Core idea
---------
Represent the momentum similarity function in a bounded spectral coordinate,

    xi(zeta) = tanh(alpha * zeta),

and obtain heat or other scalar similarity functions by applying a mode filter
to the same momentum coefficients. This yields a shared dynamical backbone for
Prandtl number, Schmidt-number analogs, and Richardson mappings.

Summary of the Gegenbauer / spherical analysis used here
--------------------------------------------------------
1. Unstable MOST power laws are generating functions of ultraspherical families.
   - phi_m ~ (1 - b_m zeta)^(-1/4) corresponds to Gegenbauer parameter lambda=1/4.
   - phi_h ~ (1 - b_h zeta)^(-1/2) corresponds to the Legendre case lambda=1/2.
2. In the matched canonical unstable case (a_h = 1, b_m = b_h = 16),
   phi_h = phi_m^2 exactly, so Ri_g = zeta and Pr_t = phi_m.
3. Away from that anchor case, scalar transport can be modeled as filtered
   momentum structure in spectral space rather than as an unrelated curve fit.

What is identifiable from theory vs data
----------------------------------------
- Usually fixed or weakly constrained by theory:
  * the basis family (lambda), often near 1/2 or 1/4 depending on the target;
  * the compactification map xi(zeta) and its width alpha;
  * the matched-case identities used as anchor constraints.
- Usually requires data:
  * scalar neutral ratios S0;
  * filter scales n_c for heat, humidity, CO2, or other tracers;
  * site/regime dependence under strong stability or intermittency.

Machine-learning guidance
-------------------------
ML can help infer residual structure or predict n_c from metadata, but the
physics-constrained spectral model should remain the backbone. Use ML to learn
parameter variation, not to replace core monotonicity and anchor identities.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np
from scipy.optimize import least_squares
from scipy.special import eval_gegenbauer


ArrayLike = Iterable[float] | np.ndarray


@dataclass(frozen=True)
class SpectralScalarClosure:
    """Shared momentum-basis closure with scalar-specific modal filters.

    Parameters
    ----------
    coeffs:
        Momentum coefficients a_n in the Gegenbauer basis.
    lambda_:
        Gegenbauer parameter. Use 1/2 for Legendre, 1/4 for the canonical
        momentum generating-function family.
    alpha:
        Compactification width in xi(zeta) = tanh(alpha * zeta).
    """

    coeffs: ArrayLike
    lambda_: float = 0.5
    alpha: float = 0.5

    def __post_init__(self) -> None:
        coeffs = np.asarray(self.coeffs, dtype=float)
        if coeffs.ndim != 1 or coeffs.size == 0:
            raise ValueError("coeffs must be a non-empty 1D array")
        if self.lambda_ <= 0.0:
            raise ValueError("lambda_ must be positive")
        if self.alpha <= 0.0:
            raise ValueError("alpha must be positive")
        object.__setattr__(self, "coeffs", coeffs)

    def xi(self, zeta: ArrayLike) -> np.ndarray:
        zeta = np.asarray(zeta, dtype=float)
        return np.tanh(self.alpha * zeta)

    def _basis_matrix(self, zeta: ArrayLike) -> np.ndarray:
        xi = self.xi(zeta)
        basis = [eval_gegenbauer(n, self.lambda_, xi) for n in range(self.coeffs.size)]
        return np.vstack(basis)

    def phi_m(self, zeta: ArrayLike) -> np.ndarray:
        basis = self._basis_matrix(zeta)
        return np.tensordot(self.coeffs, basis, axes=(0, 0))

    def filter_weights(self, n_c: float) -> np.ndarray:
        if n_c <= 0.0:
            raise ValueError("n_c must be positive")
        n = np.arange(self.coeffs.size, dtype=float)
        return np.exp(-n / n_c)

    def phi_scalar(self, zeta: ArrayLike, s0: float, n_c: float) -> np.ndarray:
        weights = self.filter_weights(n_c)
        basis = self._basis_matrix(zeta)
        return s0 * np.tensordot(self.coeffs * weights, basis, axes=(0, 0))

    def scalar_ratio(self, zeta: ArrayLike, s0: float, n_c: float) -> np.ndarray:
        phi_m = self.phi_m(zeta)
        phi_s = self.phi_scalar(zeta, s0=s0, n_c=n_c)
        return phi_s / phi_m

    def prandtl(self, zeta: ArrayLike, pr0: float = 0.85, n_c: float = 1.2) -> np.ndarray:
        return self.scalar_ratio(zeta, s0=pr0, n_c=n_c)

    def ri_g(self, zeta: ArrayLike, pr0: float = 0.85, n_c: float = 1.2) -> np.ndarray:
        zeta = np.asarray(zeta, dtype=float)
        phi_m = self.phi_m(zeta)
        phi_h = self.phi_scalar(zeta, s0=pr0, n_c=n_c)
        return zeta * phi_h / (phi_m ** 2)

    def fit_scalar_filter(
        self,
        zeta_obs: ArrayLike,
        ratio_obs: ArrayLike,
        *,
        s0_init: float = 0.85,
        n_c_init: float = 1.2,
        fit_s0: bool = True,
        bounds_s0: Tuple[float, float] = (0.2, 2.5),
        bounds_n_c: Tuple[float, float] = (0.1, 25.0),
    ) -> dict:
        """Fit scalar neutral ratio and/or filter scale from observed ratios.

        Parameters
        ----------
        zeta_obs, ratio_obs:
            Observed stability points and corresponding scalar ratios, e.g.
            Pr_t, Schmidt-number analogs, or phi_s / phi_m.
        fit_s0:
            If False, hold s0 fixed and solve only for n_c.
        """

        zeta_obs = np.asarray(zeta_obs, dtype=float)
        ratio_obs = np.asarray(ratio_obs, dtype=float)
        if zeta_obs.shape != ratio_obs.shape:
            raise ValueError("zeta_obs and ratio_obs must have the same shape")

        if fit_s0:
            x0 = np.array([s0_init, n_c_init], dtype=float)
            lower = np.array([bounds_s0[0], bounds_n_c[0]], dtype=float)
            upper = np.array([bounds_s0[1], bounds_n_c[1]], dtype=float)

            def residuals(params: np.ndarray) -> np.ndarray:
                s0, n_c = params
                return self.scalar_ratio(zeta_obs, s0=s0, n_c=n_c) - ratio_obs

        else:
            x0 = np.array([n_c_init], dtype=float)
            lower = np.array([bounds_n_c[0]], dtype=float)
            upper = np.array([bounds_n_c[1]], dtype=float)

            def residuals(params: np.ndarray) -> np.ndarray:
                (n_c,) = params
                return self.scalar_ratio(zeta_obs, s0=s0_init, n_c=n_c) - ratio_obs

        result = least_squares(residuals, x0=x0, bounds=(lower, upper))
        fitted = result.x
        if fit_s0:
            s0_fit, n_c_fit = fitted
        else:
            s0_fit, n_c_fit = s0_init, fitted[0]

        return {
            "s0": float(s0_fit),
            "n_c": float(n_c_fit),
            "cost": float(result.cost),
            "success": bool(result.success),
            "message": result.message,
        }


def exact_ubl_anchor(zeta: ArrayLike, b: float = 16.0) -> dict:
    """Return the exact matched unstable anchor solution.

    This is the analytically solvable case with

        phi_m = (1 - b zeta)^(-1/4),
        phi_h = phi_m^2,
        Ri_g = zeta,
        Pr_t = phi_m.
    """

    zeta = np.asarray(zeta, dtype=float)
    phi_m = (1.0 - b * zeta) ** (-0.25)
    phi_h = phi_m ** 2
    return {
        "phi_m": phi_m,
        "phi_h": phi_h,
        "pr_t": phi_m,
        "ri_g": zeta.copy(),
    }


def recommend_calibration_strategy() -> dict:
    """Summarize what to infer from theory, data, and ML.

    The return value is intended for notebooks and scripts that want a compact,
    machine-readable checklist.
    """

    return {
        "known_from_analysis": [
            "Use a shared momentum-basis expansion in Gegenbauer or Legendre space",
            "Enforce the matched unstable anchor phi_h = phi_m^2, Ri_g = zeta when applicable",
            "Constrain xi(zeta)=tanh(alpha zeta) and choose lambda to control resolution placement",
        ],
        "fit_from_data": [
            "Neutral scalar ratio s0",
            "Filter scale n_c for heat or other scalars",
            "Any site- or regime-dependent departures under strong stability",
        ],
        "ml_role": [
            "Predict n_c from regime metadata, roughness, season, cloud state, or intermittency markers",
            "Learn residual corrections after physics-constrained fitting",
            "Do not replace the monotone constrained backbone with an unconstrained black box",
        ],
    }
