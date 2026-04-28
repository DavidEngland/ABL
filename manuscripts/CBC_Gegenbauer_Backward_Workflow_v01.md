# CBC-Gegenbauer Backward Workflow v01

Status: working note for manuscript integration  
Date: 2026-04-26

## Purpose

Provide a compact, publication-ready workflow that starts from implemented Fortran
stability functions and works backward to the polynomial/spectral structure:

1. MOST unstable forms in zeta = z/L
2. CBC and Gegenbauer expansions in eta = -zeta
3. Legendre link via P_n(cos(theta))
4. Ri_g, Ri, Ri_c, and Pr implications by regime

## Implementation Anchors

- examples/module_cbc_legendre_most.F90
- examples/module_most_profile_utils.F90
- examples/wrf_integration_example.F90
- examples/driver_cbc_gegenbauer_errors.F90
- hw/CBC.md
- hw/Legendre.md

## Core Mathematical Chain

For heat (unstable branch):
$$
\phi_h(\zeta) = (1-b_h\zeta)^{-1/2}, \quad \eta=-\zeta>0
$$
$$
\phi_h(-\eta)=\sum_{n=0}^{\infty}\binom{2n}{n}\left(\frac{b_h\eta}{4}\right)^n
=\sum_{n=0}^{\infty}P_{2n}(0)(b_h\eta)^n
$$
with
$$
P_{2n}(0)=(-1)^n\frac{\binom{2n}{n}}{4^n}.
$$

For momentum:
$$
\phi_m(\zeta)=(1-b_m\zeta)^{-1/4}
$$
$$
\phi_m(-\eta)=\sum_{n=0}^{\infty} \left|C_{2n}^{(1/4)}(0)\right| (b_m\eta)^n
$$
(Gegenbauer equatorial representation).

Degenerate identity when b_m=b_h:
$$
\phi_h=\phi_m^2 \implies Ri_g=\zeta.
$$

Non-degenerate case b_m!=b_h:

- exact inversion requires iterative or approximation path,
- curvature coefficient Delta != 0,
- Ri-space closures diverge between heat and momentum branches.

## Regime-Aware Critical Limits

Parameter-aware limits to use in documentation and code comments:

- UBL momentum limit:

$$
Ri_{c,UBL,m}=-\frac{1}{b_m}
$$

- UBL heat limit:

$$
Ri_{c,UBL,h}=-\frac{1}{b_h}
$$

- SBL linear limit:

$$
Ri_{c,SBL}=\frac{1}{\beta}
$$
(or component-wise 1/beta_m, 1/beta_h).

Suggested conservative single UBL value for mixed closure:
$$
Ri_{c,UBL}=-\frac{1}{\max(b_m,b_h)}.
$$

## Driver Program Deliverable

The Fortran driver examples/driver_cbc_gegenbauer_errors.F90 tabulates:

- exact phi_h vs CBC approximation errors,
- exact phi_m vs Gegenbauer-equator approximation errors,
- recurrence consistency check for P_{2n}(0),
- parameter-derived Ri_c limits.

Recommended appendix table fields:

- zeta
- N (truncation order)
- phi_exact
- phi_approx
- relative error

Recommended figure:

- semilog error vs N for selected zeta values
- separate panels for heat and momentum

## Manuscript Integration Checklist

1. Theory section: add CBC/Gegenbauer/Legendre subsection and explicitly define eta=-zeta.
2. Methods section: state hybrid numerical evaluation rule (series near neutral, direct radicals away from neutral).
3. Results appendix: include driver tabulation and convergence plots.
4. Discussion: interpret spectral meaning of b_m vs b_h mismatch and impact on Ri inversion.
5. Conclusion: include parameter-aware Ri_c limits as operational guidance.
