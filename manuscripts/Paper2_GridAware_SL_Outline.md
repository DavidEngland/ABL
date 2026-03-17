# Paper 2 Outline: A Grid-Aware Surface-Layer Parameterization Unifying Log-Profile Asymptotics, Displacement Height, and Regime Transitions

**Authors:** David E. England¹, Richard T. McNider¹, Arastoo P. Biazar¹

**Affiliations:**
¹ Department of Atmospheric & Earth Science, University of Alabama in Huntsville
  (England: Research Engineer III)

**Corresponding Author:** David E. England (david.england@uah.edu)

**Target journal (candidate):** *Journal of the Atmospheric Sciences* or
*Boundary-Layer Meteorology*

**Status:** Outline / pre-draft — March 2026
**Companion design document:** `param/SCAFFOLDING.md`
**Prior paper (Paper 1):** `manuscripts/Grid_Curvature_SBL_v01.md`

---

## Abstract (placeholder)

Standard Monin–Obukhov Similarity Theory (MOST) surface-layer parameterizations
contain three interacting sources of systematic error that become large in practical
atmospheric models: (1) finite vertical grid spacing ($\Delta z$) that violates the
implicit small-parameter assumption in log-profile shear estimation; (2) omission
or inconsistent application of the canopy displacement height $d$; and (3) evaluation
of stability functions outside their valid asymptotic domain, particularly near the
critical Richardson number $Ri_c$ and in the free-convection limit.  We show that
all three errors are mathematically unified as failures of the first-order Taylor
approximation $\ln(1+\varepsilon) \approx \varepsilon$ when $\varepsilon$ is not
small.  Building on the grid-curvature correction of England et al. (Paper 1), we
develop a grid-aware surface-layer parameterization that corrects shear estimation,
enforces displacement-height consistency, and provides monotone, differentiable
regime transitions across all four stability states (unstable, near-neutral, stable,
laminar).  The scheme is validated against GABLS1 LES, 60-m tower observations,
and nocturnal stable boundary layer cases and demonstrates [X]% reduction in
near-surface flux bias at operational grid spacings without additional computational
cost.

---

## 1. Introduction

### 1.1 Motivation: Three connected failures

- Brief review of MOST and the neutral log profile as a baseline
- The well-known SBL warm bias in NWP (Holtslag et al. 2013; McNider et al. 1995)
- Paper 1 identified grid-curvature bias as one key mechanism
- Two additional, related mechanisms not addressed in Paper 1:
  - Displacement height inconsistency
  - Log-profile asymptotic domain violations across regimes
- The mathematical thread connecting all three: $\ln(1+\varepsilon)$ for non-small $\varepsilon$

### 1.2 Scope and contribution

- Extends Paper 1 from the SBL curvature correction to a full four-regime
  surface-layer scheme
- Provides the first unified mathematical treatment of grid, displacement, and
  regime-transition errors within the MOST framework
- Delivers a practical, numerically stable parameterization with documented
  interface contract for NWP integration (targeting WRF / MPAS / 1DBLM)

### 1.3 Organization

- Section 2: Mathematical unification — the $\ln(1+\varepsilon)$ theorem
- Section 3: Displacement height theory and error analysis
- Section 4: Regime catalog — four stability states, boundaries, and transitions
- Section 5: The grid-aware parameterization design
- Section 6: Validation
- Section 7: Sensitivity and discussion
- Section 8: Conclusions

---

## 2. Mathematical Unification: The $\ln(1+\varepsilon)$ Theorem

### 2.1 The core identity

$$
\ln(z + H) = \ln z + \frac{H}{z} - \frac{H^2}{2z^2} + \cdots
\quad \text{valid for } H/z < 1
$$

Three distinct physical realizations of failure when $\varepsilon \equiv H/z \gtrsim 1$:

| Context | $\varepsilon$ | Consequence |
|---------|--------------|-------------|
| Grid shear estimation | $\Delta z / z_k$ | Shear underestimate; Ri_g overestimate |
| Displacement height | $d / z$ | Shear underestimate; Ri_B overestimate (stable bias) |
| Profile offset | $(z+H)/z$ or $(z-d)/z_0$ | Profile breakdown near canopy top |

### 2.2 The three height regimes

Tripartite structure of $\ln(z+H)$:
- $z \gg H$: far-field log, first-order correction $+H/z$
- $z \sim H$: transition — neither asymptotic valid, use exact form
- $z \ll H$: near-surface linear, $\ln(z+H) \approx \ln H + z/H$

Physical implication: the log law is an asymptotic result valid for $z \gg$ all
length scales ($\Delta z$, $d$, $z_0$).  Current models routinely violate this at
their lowest 2–3 levels.

### 2.3 Unified error bounds

Derive leading-order error in $\partial U/\partial z$, $Ri_g$, and $H_s$ as a
function of $\varepsilon$ for each context.  Show that errors in $Ri_g$ are
quadratic in $\varepsilon$ (shear appears squared) and errors in $K_M$ are
amplified by the convexity of $f_m(Ri)$ near $Ri_c$.

---

## 3. Displacement Height: Systematic Analysis

### 3.1 Physical meaning and empirical estimates

- $d \approx 2h_c/3$ for vegetation; $d \approx 0.5$–$0.8 h_c$ for urban
- Three-scale ordering: $z_0 < d < h_c$
- Valid domain of the displaced log profile: $z > d + 2z_0$ at minimum

### 3.2 Log expansion with displacement

$$
\ln\frac{z-d}{z_0} = \ln\frac{z}{z_0} - \frac{d}{z} - \frac{d^2}{2z^2} - \cdots
\quad \text{(all terms negative)}
$$

Contrast with the shifted-origin case (all terms positive).

### 3.3 Shear error table

Reproduce the table from `drafts/log-expansions.md` §4.3 (forest with $d = 10$ m):
shear error ranges from 500% at $z = 12$ m to 11% at $z = 100$ m.

### 3.4 Bulk Richardson number bias

Show:
$$
Ri_B^{\text{naive}} / Ri_B^{\text{correct}} = \frac{\ln(z_2/z_1)}{\ln((z_2-d)/(z_1-d))} > 1
$$

The naive estimate overestimates $Ri_B$ — spuriously stabilizing.  Provide
concrete numerical examples at three land-use categories.

### 3.5 Implications for NWP

Most operational models use $d = 0$ globally or set $d$ as a static lookup table
uncoupled from the flux scheme.  Quantify the resulting bias in $C_D$, $C_H$, and
$H_s$ for urban, forest, and cropland surfaces.

---

## 4. Stability Regimes: A Four-State Framework

### 4.1 Unstable (daytime convective)

- Businger–Dyer gradient functions, Paulson (1970) $\psi_m$, $\psi_h$
- Key result: $Ri_g = \zeta$ exactly for B-D forms — not an approximation
- Free-convection limit ($|L| \to 0$): breakdown of wind-based scaling, $w_*$ regime
- Pitfall: $\psi_m$ is as large as $\ln(z/z_0)$ under strong instability — non-perturbative

### 4.2 Near-neutral transition

- Linearized expressions for $\psi_m$ and $\psi_h$
- Why the asymmetry between stable and unstable corrections matters over rough terrain
- Obukhov length range for near-neutral: $|L| > z/0.05$
- NWP sign-of-stability oscillation: measurement noise, timestep sensitivity

### 4.3 Stable (weakly to strongly stable)

- England–McNider (1995) $f_m(Ri) = (1 - Ri/Ri_c)^2$
- Critical Richardson number as asymptote, not threshold: $Ri_c = 1/\beta$
- Iteration for $\zeta \leftrightarrow Ri$ inversion; near-singularity behavior
- Compound failure: grid-curvature shear underestimate + MOST amplification near $Ri_c$

### 4.4 Laminar

- Turbulence collapse; $f_m \to 0$
- Molecular floor: $K_M \geq \nu$, $K_H \geq \alpha$ prevents runaway cooling
- Intermittency: ensemble mean flux ≈ 20–40% of turbulent value
- Flip-flop bi-stability (McNider et al.): two equilibria, path-dependence,
  implications for NWP forecast determinism

---

## 5. Grid-Aware Parameterization Design

### 5.1 Interface specification

Inputs, outputs, and diagnostic fields (reproduce §3 of `param/SCAFFOLDING.md`).

### 5.2 Regime classifier with grid-curvature correction

The curvature-corrected $Ri_g$:

$$
Ri_g^{\text{corr}} = Ri_g^{\text{FD}} \cdot \left(1 - \frac{\Delta z}{2z_{\text{eff}}}\right)^2
$$

where $z_{\text{eff}} = z_k - d$.  Show this removes the leading-order bias before
regime classification — prevents spurious laminar transitions.

### 5.3 Grid-aware shear: exact log ratio

$$
\left\langle \frac{\partial U}{\partial z} \right\rangle_{\Delta z}
= \frac{u_*}{\kappa \Delta z} \ln\frac{z_{\text{eff},k+1}}{z_{\text{eff},k}}
$$

versus the Taylor-approximated (standard) form.  Table: exact vs. approximate
at five grid spacings from 2 m to 100 m.

### 5.4 Displacement consistency contract

Formal statement of the five-point contract (reproduce from `param/SCAFFOLDING.md`
§8).  This is the operational heart of the scheme — its role in eliminating
silent code bugs.

### 5.5 Stability function implementation

- Branch selection (pseudocode from §6 of SCAFFOLDING.md)
- Monotonicity guarantee across regime boundaries
- Molecular floor implementation

### 5.6 Iterative flux inversion

- Algorithm (pseudocode §7 of SCAFFOLDING.md)
- Near-$Ri_c$ damping
- Free-convection guard
- Convergence criterion and non-convergence handling

---

## 6. Validation

### 6.1 Analytic test cases

- Neutral: recover exact $u_*$ from synthetic two-level profiles
- Stable: recover $\theta_*$ and $L$ for synthetic SBL; compare corrected vs.
  uncorrected $Ri_B$
- Laminar: confirm $K_M \geq \nu$ and no runaway cooling over 6-hour integration

### 6.2 GABLS1 LES benchmark

- Target: GEWEX ABL Study 1 idealized stable BL (Cuxart et al. 2006)
- Metrics: surface temperature, $u_*$, $H_s$ over 9-hour integration
- Compare: standard scheme (no corrections) vs. grid-aware scheme at
  $\Delta z$ = 10, 50, 100, 200 m

### 6.3 60-m tower observations

- Use tower data from `implementations/Met Tower (60 m).md` and `data/`
- Compare: inferred fluxes (corrected vs. uncorrected) against eddy-covariance reference
- Expected: largest improvement in stably stratified nocturnal cases with $z_0 > 0.1$ m

### 6.4 Sensitivity to $d$ and $Ri_c$

- Vary $d$ from 0 to $2h_c/3$: show $Ri_B$ and $H_s$ sensitivity
- Vary $Ri_c$ from 0.167 ($\beta = 6$) to 0.25 ($\beta = 4$): show effect on
  turbulence collapse timing and nocturnal minimum temperature

---

## 7. Discussion

### 7.1 When do the corrections matter?

Rank by regime and surface type:

| Surface | Key correction | Magnitude |
|---------|---------------|-----------|
| Urban canyon | Displacement $d$: 40–80% of $h_c$; shear error 50–200% within 2$h_c$ | Large |
| Forest / tall crops | Displacement + grid curvature: compound effect | Large |
| Short grass | Grid curvature dominant; $d$ small | Moderate |
| Open ocean | $d = 0$; grid curvature only; $z_0$ very small | Small–moderate |
| Arctic tundra | No displacement; SBL very shallow; near-$Ri_c$ behavior dominant | Moderate–large |

### 7.2 Implications for the SBL warm bias

- Paper 1 showed grid-curvature correction reduces SBL bias by 40–55% (GABLS1)
- Paper 2 adds displacement consistency and regime-transition corrections
- Expected additional bias reduction: quantify from validation
- Discuss residual bias: what is NOT fixed by this scheme?
  (Intermittency, gravity waves, RSL effects, heterogeneity)

### 7.3 Implications for polar and high-latitude modeling

- Arctic surface layers are shallow; $\Delta z/z$ is large nearly everywhere
- Low vegetation; $d$ small but $z_0$ variable
- $Ri > Ri_c$ common; laminar branch activated frequently
- Connection to Arctic Amplification (`Arctic Amplification.md`): correct
  nocturnal cooling is critical for sea-ice–atmosphere feedbacks

### 7.4 Computational cost

- Exact $\ln$ calls replace Taylor approximations: negligible cost ($< 1$% overhead)
- Iterative inversion: typically 3–5 iterations; convergence faster than standard
  because corrected $Ri_g$ is a better initial guess for $\zeta$

### 7.5 Extension to dynamic $Ri_c$

- Brief description of the `Dynamic_Ric_Hybrid_MOST_Ri_Draft.md` framework
- How the interface contract (§5.4) facilitates plug-in of dynamic $Ri_c$
- This is the natural Paper 3 topic

---

## 8. Conclusions

1. Log-profile errors from grid spacing, displacement height, and asymptotic
   domain violations are mathematically unified by the $\ln(1+\varepsilon)$
   small-parameter breakdown.

2. The displacement height bias biases $Ri_B$ toward stability (overestimate)
   even for modest $d = 1$ m; at forest scales the effect is large at all
   operational model levels.

3. A formal displacement-consistency contract eliminates silent code bugs that
   mix $z$ and $z-d$ across profile and flux routines — a class of error not
   addressed in prior parameterization literature.

4. The grid-aware scheme recovers the exact finite-difference log-ratio shear,
   corrects $Ri_g$ before regime classification, applies the appropriate
   four-state stability function, and enforces a molecular diffusivity floor.

5. Validation shows [X]% reduction in SBL flux bias at operational grid spacing
   with no additional computational cost.

---

## Appendix A — Log-Expansion Proofs

Full derivations of:
- $\ln(z+H)$ series and convergence domain
- $\ln(z-d)$ series with sign analysis
- Grid shear bias $\Delta z/(2z_k)$ and quadratic propagation to $Ri_g$
- Bulk $Ri_B$ ratio with/without $d$

## Appendix B — Stability Function Derivations

- Paulson (1970) $\psi_m$, $\psi_h$ integrals (full derivation)
- England–McNider $f_m = (1-Ri/Ri_c)^2$ derivation from $\phi_m = 1+\beta\zeta$
- Molecular floor limit and its physical justification

## Appendix C — Pseudocode Listing

Full pseudocode for: `classify_regime`, `compute_shear`, `psi_m_unstable`,
`psi_h_unstable`, `f_m_stable`, `K_floor`, `invert_fluxes` — cross-referenced
to `param/SCAFFOLDING.md`.

---

## Figures (planned)

1. **Fig. 1** — Schematic of the three height regimes and the four stability
   regimes; the two-axis structure governing this paper.
2. **Fig. 2** — Shear error table (log scale) as a function of $z/d$ and $\Delta z/z$;
   shows the compound error surface.
3. **Fig. 3** — $Ri_B$ overestimate as a function of layer depth and $d$, for three
   surface types.
4. **Fig. 4** — Four-panel stability function summary: $\phi_m$, $\psi_m$, $f_m$,
   $K_M/(\kappa u_* z)$ vs. $\zeta$ (or $Ri$) across all regimes.
5. **Fig. 5** — GABLS1 validation: surface temperature and $u_*$ time series,
   corrected vs. uncorrected, at three grid spacings.
6. **Fig. 6** — Tower comparison: scatter plot of corrected vs. observed fluxes
   stratified by regime.
7. **Fig. 7** — Sensitivity: nocturnal minimum temperature vs. $d$ and $Ri_c$;
   shows operational importance of correct values.

---

## Action Items for McNider & Biazar Meeting

1. **McNider:** Confirm target journal — BLM vs. JAS vs. JAMES?
   - BLM: right scope, 6–8 month review cycle
   - JAS: higher impact; expects broader atmospheric science framing
   - JAMES: open-access; good for parameterization papers

2. **Biazar:** Review displacement height formulation (§3) — is the three-scale
   ordering ($z_0 < d < h_c$) standard in UAH's surface datasets?  Any existing
   $d$ estimates for the tower sites?

3. **McNider:** The flip-flop bi-stability discussion (§4.4) — should we cite
   the McNider et al. original result more prominently in the intro, or let the
   math speak for itself?

4. **All:** Decide on target NWP model for demonstration:
   - 1DBLM (McNider's column model) — fastest to implement
   - WRF — broadest community impact
   - MPAS — relevant for global applications

5. **England:** Stand up GABLS1 test harness using existing `src/rct/` and
   `implementations/` code before next meeting; get baseline uncorrected numbers.

6. **All:** Review `param/SCAFFOLDING.md` and resolve the five open design
   questions (§11) before moving to implementation.

---

*England, McNider, Biazar — UAH Atmospheric & Earth Science — March 2026*
