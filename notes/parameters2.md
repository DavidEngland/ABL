# Unified MOST Parameter Framework and Flux-Profile Structure

Draft note for internal review by McNider and Biazar.
Status: polished Markdown pre-draft for later LaTeX conversion.

## 1. Scope and Positioning

This note consolidates a unified parameter framework for Monin-Obukhov similarity theory (MOST) closures, with explicit links among:

- profile functions $\phi_m(\zeta)$ and $\phi_h(\zeta)$,
- flux-profile integrals $\psi_m(\zeta)$ and $\psi_h(\zeta)$,
- Richardson mappings $Ri_g(\zeta)$,
- regime-aware critical thresholds ($Ri_c^{UBL}$, $Ri_c^{SBL}$),
- and inversion strategy for Arctic/nocturnal stable boundary layers (SBL).

The target is a manuscript-ready theory/method section with parameter tables and implementation hooks for WRF and related models [@Businger_1971; @Dyer1974; @BeljaarsHoltslag1991JAMC; @Holtslag_2013; @Grachev_2012b].

## 2. Canonical Definitions

Define

$$
\zeta = \frac{z}{L}, \qquad
\phi_m(\zeta) = \frac{\kappa z}{u_*}\frac{\partial U}{\partial z}, \qquad
\phi_h(\zeta) = \frac{\kappa z}{\theta_*}\frac{\partial \Theta}{\partial z}.
$$

The flux-profile integrals are

$$
\psi_m(\zeta) = \int_0^\zeta \left(1 - \frac{1}{\phi_m(\zeta')}\right)\frac{d\zeta'}{\zeta'},
\qquad
\psi_h(\zeta) = \int_0^\zeta \left(1 - \frac{1}{\phi_h(\zeta')}\right)\frac{d\zeta'}{\zeta'}.
$$

The gradient Richardson number relation is

$$
Ri_g(\zeta) = \zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}.
$$

## 3. Unified Parameterization Families

### 3.1 Canonical Businger-Dyer Type (Kansas lineage)

Unstable branch ($\zeta<0$):

$$
\phi_m = (1-b_m\zeta)^{-1/4},
\qquad
\phi_h = a_h^{-1}(1-b_h\zeta)^{-1/2}.
$$

Stable branch ($\zeta>0$):

$$
\phi_m = 1+\beta_m\zeta,
\qquad
\phi_h = a_h^{-1}+\beta_h\zeta.
$$

Parameter table:

| Family / source | Branch | $\kappa$ | $a_h^{-1}$ | $b_m$ | $b_h$ | $\beta_m$ | $\beta_h$ | Nominal range | Notes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| Businger et al. (1971) [@Businger_1971] | Unstable | 0.35 | 1.35 | 15.0 | 9.0 | - | - | $-2 \lesssim \zeta < 0$ | Scalar-momentum mismatch |
| Businger et al. (1971) [@Businger_1971] | Stable | 0.35 | 1.35 | - | - | 4.7 | 6.35 | $0 < \zeta \lesssim 1$ | Asymmetric stable slopes |
| Dyer (1974) / Brutsaert (1982) [@Dyer1974] | Unstable | 0.40 | 1.00 | 16.0 | 16.0 | - | - | $-2 \lesssim \zeta < 0$ | Degenerate case: $\phi_h=\phi_m^2$ |
| Dyer (1974) / Brutsaert (1982) [@Dyer1974] | Stable | 0.40 | 1.00 | - | - | 5.0 | 5.0 | $0 < \zeta \lesssim 1$ | Common WRF/LES default |
| Hogstrom (1996) | Unstable | 0.40 | 1.00 | 19.0 | 11.6 | - | - | $-2 \lesssim \zeta < 0$ | Stronger momentum curvature |
| Hogstrom (1996) | Stable | 0.40 | 1.00 | - | - | 5.3 | 5.3 | $0 < \zeta \lesssim 1$ | Slightly steeper stable branch |

### 3.2 Regime map

| Regime | $\zeta$ range | Typical $Ri$ range | Physical character |
| --- | --- | --- | --- |
| Near-neutral unstable | $-1 < \zeta < 0$ | $Ri \lesssim 0$ | Weak buoyant enhancement |
| Strongly unstable | $\zeta < -1$ | $Ri \ll 0$ | Convective plume dominance |
| Moderately stable | $0 < \zeta \le 1$ | $0.02 < Ri < 0.12$ | Continuous but suppressed turbulence |
| Very/extremely stable | $\zeta > 1$ | $0.12 < Ri < 0.7$ | Intermittency, z-less behavior |
| Critical transition | - | $Ri_c \approx 0.21\text{--}0.25$ | Onset of laminarization risk |

### 3.3 Specialized stable forms

These families were introduced to avoid stable-branch over-suppression and runaway decoupling [@BeljaarsHoltslag1991JAMC; @Holtslag_2013].

| Scheme | Variable | Core constants | Nominal range | Practical role |
| --- | --- | --- | --- | --- |
| Beljaars-Holtslag (1991) [@BeljaarsHoltslag1991JAMC] | $m,h$ | $a=1, b=0.667, c=5, d=0.35$ | $0<\zeta\lesssim10$ | Smooth very-stable transition |
| Cheng-Brutsaert (2005) | $m$ | $a=6.1, b=2.5$ | $0<\zeta\lesssim20$ | Stronger curvature |
| Cheng-Brutsaert (2005) | $h$ | $c=5.3, d=1.1$ | $0<\zeta\lesssim20$ | Strong scalar damping |
| SHEBA (2007) | $m$ | $a_m=5, b_m=0.3$ | $0<\zeta\lesssim100$ | Arctic z-less behavior |
| SHEBA (2007) | $h$ | $a_h=5, b_h=0.4$ | $0<\zeta\lesssim100$ | Polar inversion scalar suppression |

## 4. Flux-Profile Integrals and Structural Limits

### 4.1 Unstable Kansas integrals

With

$$
x=(1-b_m\zeta)^{1/4}, \qquad y=(1-b_h\zeta)^{1/2}, \qquad \zeta<0,
$$

use

$$
\psi_m = 2\ln\frac{1+x}{2} + \ln\frac{1+x^2}{2} - 2\arctan x + \frac{\pi}{2},
$$

$$
\psi_h = 2\ln\frac{1+y}{2}.
$$

These are the standard Kansas forms used in many operational and research codes [@Businger_1971; @Dyer1974].

### 4.2 Stable linear branch integrals

For

$$
\phi_m=1+\beta_m\zeta, \qquad \phi_h=a_h^{-1}+\beta_h\zeta,
$$

one obtains

$$
\psi_m = -\beta_m\zeta,
\qquad
\psi_h = -a_h\beta_h\zeta
$$

(assuming neutral matching constants set conventionally).

### 4.3 CBC/Gegenbauer structure and dynamic limits

Using $\eta=-\zeta>0$ on the unstable branch,

$$
\phi_h(-\eta) = (1+b_h\eta)^{-1/2}
= \sum_{n=0}^{\infty} \binom{2n}{n}\left(\frac{b_h\eta}{4}\right)^n.
$$

For momentum,

$$
\phi_m(-\eta) = (1+b_m\eta)^{-1/4}
$$

has the corresponding Gegenbauer-equatorial expansion (parameter $\lambda=1/4$), while heat corresponds to Legendre ($\lambda=1/2$) via $P_n(\cos\theta)$.

Natural parameter-aware Richardson bounds are:

$$
Ri_{c,UBL,m}=-\frac{1}{b_m},
\qquad
Ri_{c,UBL,h}=-\frac{1}{b_h},
\qquad
Ri_{c,SBL}=\frac{1}{\beta}
$$

(or component-wise $1/\beta_m$, $1/\beta_h$).

Degenerate identity case (Dyer/Brutsaert unstable branch, $b_m=b_h=b$, $a_h^{-1}=1$):

$$
\phi_h=\phi_m^2 \implies Ri_g=\zeta.
$$

This limit is analytically valuable for validating inversion solvers and closure consistency.

## 5. Three Peer Reviews (Simulated) and Critique

### Review 1: Boundary-layer physics reviewer

Main critique:

- Strength: clear unification of Kansas and very-stable families.
- Concern: some regime boundaries are presented as universal, but should be framed as dataset/model dependent.
- Concern: explicit treatment of intermittency and wave-turbulence coupling is still too brief for Arctic applications.

Required revisions:

- Add wording that $Ri$ thresholds are "reference ranges" not hard cutoffs.
- Add one paragraph on intermittency and non-local transport caveats in very stable SBL.
- Clarify where MOST assumptions are formally violated.

### Review 2: Numerical methods / model implementation reviewer

Main critique:

- Strength: equations are implementation-oriented.
- Concern: no explicit algorithm statement for switching among linear, specialized stable, and MOST-inversion forms.
- Concern: no mention of safeguarding Newton inversion and clipping/bounds strategy.

Required revisions:

- Add pseudo-algorithm for branch selection and fallback.
- Define recommended bounds for $\zeta$, $Ri$, and iteration tolerances.
- Add one short reproducibility table: defaults, limits, switches.

### Review 3: Applied-math / spectral-theory reviewer

Main critique:

- Strength: CBC to Legendre/Gegenbauer link is novel and potentially high impact.
- Concern: relationship to $P_n(\cos\theta)$ is stated but not explicitly derived in this note.
- Concern: notation switches among $\gamma$, $b$, and $\beta$ are still too easy to misread.

Required revisions:

- Add an appendix derivation block from generating function to CBC coefficients.
- Use one canonical notation throughout body: $b_m,b_h,\beta_m,\beta_h$.
- Keep legacy symbols ($\gamma_m,\gamma_h$) only in a single mapping note/table cell.

## 6. Actioned Improvements in This Markdown Draft

Already implemented in this version:

- All math now uses strict `$...$` and `$$...$$` delimiters.
- Unified notation uses $b_m,b_h,\beta_m,\beta_h$ as primary symbols.
- Added parameter-aware bounds $Ri_{c,UBL,m}$, $Ri_{c,UBL,h}$, $Ri_{c,SBL}$.
- Recast tables into publication-oriented Markdown format.
- Converted narrative to BibTeX-ready citation style (`[@key]`).

Remaining before LaTeX conversion:

- Add explicit short derivation of Legendre/Gegenbauer generating-function mapping.
- Add implementation pseudo-code and solver safeguards subsection.
- Cross-check citation keys against `manuscripts/references.bib` for exact naming.

## 7. Suggested BibTeX Key Set (to normalize)

Use a single canonical key set across manuscript files:

- `@Businger_1971`
- `@Dyer1974`
- `@BeljaarsHoltslag1991JAMC`
- `@Holtslag_2013`
- `@Cuxart2006BLM`
- `@Bosveld2014BLM`
- `@Grachev_2012b`
- `@EnglandMcNider1995BLM`

If key names differ in source `.bib` files, normalize during merge into `manuscripts/references.bib`.

## 8. Editorial Recommendation to Coauthors

This draft is now suitable for scientific-content review by McNider and Biazar before LaTeX typesetting. The next pass should focus on:

- physical framing in very-stable/intermittent regimes,
- model-implementation details and safeguards,
- and citation-key normalization.
