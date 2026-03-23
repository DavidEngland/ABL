# Grid-Aware Surface-Layer Parameterization: Conceptual Scaffolding
*England, McNider, Biazar — UAH Atmospheric & Earth Science*
*Research Engineer III: David E. England*
*March 2026 — Design document, no production code yet*

---

## 0. Purpose and Provenance

This document captures the **conceptual architecture** for a new surface-layer
parameterization that corrects the systematic failures identified in:

- `drafts/log-expansions.md` — log-expansion breakdown, displacement height errors, 17 pitfalls
- `manuscripts/Grid_Curvature_SBL_v01.md` — grid-curvature bias and Ri_g curvature correction (Paper 1)
- `hw/neutral-log-profiles.md` — pedagogical derivations across all four regimes

The central unifying insight is:

> **MOST is not wrong — it is being used outside its asymptotic domain.  The
> resulting errors in shear, Richardson number, stability functions, and surface
> fluxes are mathematically unified as failures of the small-parameter
> approximation ln(1 + ε) ≈ ε when ε = Δz/z (grid), ε = d/z (displacement),
> or ε = H/z (profile offset) is not small.**

The new parameterization operationalizes this insight as a single, modular,
grid-aware surface-layer scheme.

---

## 1. Scope and Design Goals

### What this parameterization does

1. **Corrects the log-profile for finite grid spacing** — uses the exact
   $\ln(1 + \Delta z/z_k)$ in shear estimation, not the Taylor-linearized
   $\Delta z/z_k$.

2. **Incorporates displacement height $d$ consistently** — every profile
   evaluation, tendency, and flux uses $z - d$, never $z$.

3. **Spans all four stability regimes deterministically** — neutral (log),
   unstable (Businger–Dyer $\psi$-corrected), near-neutral (log-linear
   transition), stable (England–McNider $f_m$), and laminar (molecular floor).

4. **Is grid-aware** — the scheme knows its own $\Delta z$ and adjusts the
   curvature correction accordingly, rather than assuming the grid is fine.

5. **Has a clean regime-transition architecture** — no hard cliffs; transitions
   are monotone and differentiable.

6. **Is self-consistent** — the same $z - d$, $z_0$, $\zeta$, $L$ values are
   used for profiles, fluxes, and Richardson number.  No split bookkeeping.

### What this parameterization does NOT do

- Replace a full PBL scheme above the surface layer.
- Resolve the roughness sublayer (that requires a separate canopy submodel).
- Make multi-column or 3D corrections.
- Handle heterogeneous surfaces (blending-height theory is a separate extension).

---

## 2. Module Decomposition

```
param/
├── SCAFFOLDING.md              ← this file
├── DESIGN_DECISIONS.md         ← rationale log (to be written as design matures)
│
├── core/
│   ├── displacement.md         ← displacement height: definitions, valid domain
│   ├── regimes.md              ← regime detection: stability classifier logic
│   ├── profiles.md             ← profile functions: exact, corrected, composite
│   ├── gradients.md            ← shear and temperature gradient: exact vs. FD
│   ├── fluxes.md               ← u*, θ*, τ, Hs, C_D, C_H derivations
│   ├── richardson.md           ← Ri_g, Ri_B, curvature correction, near-Ric behavior
│   └── mixing_length.md        ← l_m hierarchy: Prandtl, Blackadar, stability-modified
│
├── corrections/
│   ├── grid_curvature.md       ← Δz/(2z_k) correction to shear
│   ├── displacement_bias.md    ← Ri_B overestimate from missing d
│   └── asymptotic_composite.md ← matched expansion across height regimes
│
├── stability_functions/
│   ├── neutral.md              ← φ_m = φ_h = 1 baseline
│   ├── unstable_BD.md          ← Businger–Dyer; ψ_m Paulson form
│   ├── stable_linear.md        ← England–McNider; f_m = (1 - Ri/Ri_c)²
│   └── laminar_floor.md        ← K_M = max(f_m·κu*z, ν)
│
├── iteration/
│   ├── solver.md               ← iterative ζ ↔ Ri inversion strategy
│   └── near_neutral_guard.md   ← sign-of-stability guard rails
│
└── validation/
    ├── test_cases.md           ← analytic test cases (known exact solutions)
    └── gabls1_comparison.md    ← target: GABLS1 LES benchmark
```

The `.md` files in each subfolder are **design specifications**, not code.
When implementation begins, each `.md` becomes a companion to a `.py` (or `.f90`)
module of the same name in `src/rct/` or `implementations/`.

---

## 3. Data Flow and Interface Contract

### 3.1 Inputs (per layer interface call)

```
INPUTS (at every model timestep, per column):
  z_k         : height of level k above ground (m)         -- geometric height
  z_k1        : height of level k+1                        -- always > z_k
  dz          : layer thickness = z_k1 - z_k               -- Δz, the grid spacing
  d           : displacement height (m)                    -- surface property
  z0m         : aerodynamic roughness length (m)           -- momentum
  z0h         : scalar roughness length (m)                -- heat/moisture
  U_k, U_k1  : wind speed at k, k+1 (m/s)
  Th_k, Th_k1 : potential temperature at k, k+1 (K)
  Th_bar      : mean layer potential temperature (K)
  g           : gravitational acceleration (m/s²)
  kappa       : von Kármán constant (≈ 0.40)
  nu          : kinematic molecular viscosity (m²/s)       -- laminar floor
```

Note: `z_k` and `z_k1` are already verified to satisfy `z_k > d` before the call.
The caller is responsible for this guard — the scheme is undefined below `z = d`.

### 3.2 Primary outputs

```
OUTPUTS:
  u_star      : friction velocity (m/s)
  theta_star  : temperature scale (K)
  tau         : kinematic momentum flux = u_star² (m²/s²)
  Hs          : kinematic sensible heat flux = -u_star·theta_star (K·m/s)
  K_M         : eddy diffusivity for momentum (m²/s)
  K_H         : eddy diffusivity for heat (m²/s)
  Ri_g_corr   : grid-corrected gradient Richardson number
  regime      : integer flag  {UNSTABLE=-1, NEAR_NEUTRAL=0, STABLE=1, LAMINAR=2}
  zeta        : stability parameter z/L at level k
  L           : Obukhov length (m)
```

### 3.3 Derived diagnostics (optional, for analysis)

```
DIAGNOSTICS (optional):
  shear_exact         : u_*/[κ(z_k - d)]
  shear_FD            : (U_k1 - U_k) / dz     [numerically computed]
  shear_bias          : shear_FD / shear_exact - 1   [grid-curvature error]
  Ri_B_naive          : Ri_B ignoring d
  Ri_B_corrected      : Ri_B with d included in log ratio
  psi_m, psi_h        : stability correction integrals
  f_m, f_h            : stability reduction functions
  C_D, C_H            : bulk transfer coefficients
```

---

## 4. Regime Classifier (pseudocode)

The regime classification gates which stability function branch is used.  It must
be monotone and must not oscillate on consecutive timesteps.

```
FUNCTION classify_regime(Ri_g, zeta, dz, z_k, d):

    -- Step 1: compute effective stability parameter
    z_eff = z_k - d                         -- always use displaced coordinate

    -- Step 2: grid-curvature correction to Ri_g
    --   True shear at z_eff ≈ u*/(κ·z_eff)
    --   FD shear has multiplicative bias (1 - dz/(2·z_eff))
    --   So FD Ri_g is biased HIGH by 1/(1 - dz/(2z_eff))^2 - 1
    bias_factor = dz / (2 * z_eff)
    IF bias_factor >= 1.0 THEN
        WARN "grid too coarse for this level: dz >= 2*(z_k - d)"
        -- clamp and flag; do not attempt MOST
        RETURN regime = UNDEFINED
    END IF
    Ri_g_corr = Ri_g * (1 - bias_factor)^2   -- remove FD overestimate

    -- Step 3: assign regime
    IF Ri_g_corr < -0.05 THEN
        regime = UNSTABLE
    ELSE IF Ri_g_corr < 0.05 THEN
        regime = NEAR_NEUTRAL
    ELSE IF Ri_g_corr < Ri_c THEN
        regime = STABLE
    ELSE
        regime = LAMINAR
    END IF

    RETURN regime, Ri_g_corr
END FUNCTION
```

The thresholds ±0.05 are a design parameter; they should be validated against
the $|\zeta| < 0.05$ near-neutral criterion from `hw/neutral-log-profiles.md`
Problem 7B.

---

## 5. Profile and Shear Evaluation (pseudocode)

### 5.1 Wind shear — exact vs. grid-aware finite difference

```
FUNCTION compute_shear(U_k, U_k1, z_k, z_k1, d, u_star, kappa):

    z_eff_k  = z_k  - d
    z_eff_k1 = z_k1 - d

    -- Analytic (MOST) shear at z_eff_k   [exact for log profile]
    shear_analytic = u_star / (kappa * z_eff_k)

    -- Finite-difference shear  [available from model state]
    dz = z_k1 - z_k
    shear_FD = (U_k1 - U_k) / dz

    -- Exact finite-difference from log profile (no Taylor approximation):
    --   = (u*/κdz) · ln(z_eff_k1 / z_eff_k)
    shear_FD_exact = (u_star / (kappa * dz)) * ln(z_eff_k1 / z_eff_k)

    -- Grid-curvature bias from first-order Taylor error
    --   shear_FD ≈ shear_analytic · (1 - dz/(2·z_eff_k) + ...)
    curvature_bias = dz / (2 * z_eff_k)

    RETURN shear_analytic, shear_FD, shear_FD_exact, curvature_bias
END FUNCTION
```

**Design note:** `shear_FD_exact` uses the *exact* log ratio rather than its
Taylor approximation — this is the central correction.  The ratio
$\ln(z_{\mathrm{eff},k+1}/z_{\mathrm{eff},k}) / (z_{k+1} - z_k)$ is the
corrected mean gradient over the layer.

### 5.2 Log-profile with $\psi$ correction (general)

```
FUNCTION log_profile_wind(z, u_star, z0m, d, kappa, psi_m):
    z_eff = z - d
    ASSERT z_eff > z0m    -- below this, profile undefined
    U = (u_star / kappa) * (ln(z_eff / z0m) - psi_m)
    RETURN U
END FUNCTION
```

The `psi_m` is computed by the appropriate stability-function module and passed
in — the profile function itself is regime-agnostic.

---

## 6. Mixing Length: Theory, Hierarchy, and Grid-Aware Role

### 6.1 Why mixing length matters here

The eddy diffusivity in K-theory can be written in two equivalent ways:

$$
K_M = \kappa\, u_*\, z_{\text{eff}} \cdot f_m(Ri)
\qquad \text{(MOST / surface-layer form)}
$$

$$
K_M = l_m^2 \left|\frac{\partial U}{\partial z}\right| \cdot f_m(Ri)
\qquad \text{(Prandtl mixing-length form)}
$$

These are identical in the log layer because $\partial U/\partial z = u_*/(\kappa z_{\text{eff}})$
implies $l_m = \kappa z_{\text{eff}}$.  But above the surface layer the two
forms diverge: MOST assumes a constant flux layer that no longer holds, while the
Prandtl form can be extended to the full PBL if the mixing length is specified
properly across height.  The mixing length thus bridges the surface-layer scheme
to the overlying PBL closure and determines whether the scheme is
self-consistent when stacked in a column model.

### 6.2 The Prandtl (surface-layer) mixing length

In the log layer under neutral conditions:

$$
l_m = \kappa\,(z - d)
$$

This grows linearly without bound, which is unphysical above the surface layer.
The shear from the log profile gives:

$$
K_M = l_m^2 \frac{\partial U}{\partial z}
    = [\kappa(z-d)]^2 \cdot \frac{u_*}{\kappa(z-d)}
    = \kappa\, u_*\,(z-d)
    \quad \checkmark
$$

So the MOST form and the mixing-length form give identical $K_M$ in the surface
layer — they are the same closure, written differently.

### 6.3 The Blackadar (1962) composite form

To prevent $l_m = \kappa z$ growing without limit, Blackadar proposed blending
the near-surface scaling with an **asymptotic mixing length** $\lambda$:

$$
\frac{1}{l_m} = \frac{1}{\kappa(z - d)} + \frac{1}{\lambda}
$$

Rearranging:

$$
l_m = \frac{\kappa(z-d)}{1 + \kappa(z-d)/\lambda}
$$

Limiting behaviour:
- $z \to d$ (near surface): $l_m \to \kappa(z-d)$ — recovers Prandtl / log law
- $z \to \infty$: $l_m \to \lambda$ — asymptotic value
- Transition height: $z - d \sim \lambda/\kappa$, typically $\sim 50$–$200\,\text{m}$

Common choices for $\lambda$:

| Choice | Expression | Typical value | Context |
|--------|-----------|---------------|---------|
| Blackadar (1962) | $\lambda \propto \int_0^h |U|\,dz / |U|_{\max}$ | 30–80 m | Neutral BL |
| Fixed fraction of BL height | $\lambda = c_\lambda h$ | $c_\lambda \approx 0.1$, $h \sim 1\,\text{km}ight \Rightarrow \lambda \sim 100\,\text{m}$ | Daytime ML |
| Obukhov-length limited (stable) | $\lambda = c_s |L|$ | $c_s \approx 0.5$ | SBL; Nieuwstadt 1984 |
| TKE-based | $\lambda = c_\mu e^{1/2}/N$ | depends on TKE $e$, $N$ | Second-order closures |

For a surface-layer-only scheme, the most tractable choice is the **Obukhov-length
limited** form in the stable branch and the **fixed fraction of BL height** in the
unstable branch.  The BL height $h$ must then be an additional input from the
model's PBL scheme.

### 6.4 Stability modification

With stability, the physically correct mixing length is reduced.  Two formulations
are in common use:

**Option A — Multiplicative:** apply $f_m$ after computing the geometric mixing length:
$$
l_m^{(\text{eff})} = l_m^{(\text{Blackadar})} \cdot \sqrt{f_m(Ri)}
$$
then $K_M = {l_m^{(\text{eff})}}^2 |\partial U/\partial z|$.  This is consistent
with the MOST surface-layer form because $f_m$ enters $K_M$ linearly and $l_m^2$
enters quadratically — the $\sqrt{} $ keeps the product correct.

**Option B — Length-scale limited (Nieuwstadt 1984):** blend with a buoyancy
mixing length $l_b = c_s L$ for stable conditions:

$$
\frac{1}{l_m} = \frac{1}{\kappa(z-d)} + \frac{1}{\lambda} + \frac{1}{c_s L}
\qquad (L > 0,\; \text{stable})
$$

The third term activates only in the stable regime and limits $l_m$ to the
Obukhov scale as $L$ shrinks.  This is the **three-term reciprocal** form the
user notes ($1/l_m = 1/(\kappa z) + a_1/l_2 + ...$).

The three-term form unifies the Prandtl near-surface limit, the Blackadar PBL
scale-height limit, and the Obukhov stability limit into one smooth expression:

$$
\boxed{
\frac{1}{l_m} = \frac{1}{\kappa(z-d)} + \frac{1}{\lambda} + \frac{c_s^{-1}}{L}
}
\qquad (\text{stable}; \; L > 0)
$$

For unstable and neutral conditions, only the first two terms are used
($L < 0$ would add a negative third term, increasing $l_m$ — physically correct
but usually absorbed into the $\psi_m$ correction instead).

### 6.5 Connection to the grid-curvature correction

The mixing length sets the **intrinsic turbulent length scale** at each level.
When the **grid spacing $\Delta z > l_m$**, the model cannot resolve the
turbulence structure that the mixing length describes — the column sees only a
layer-average, not the local gradient.

Near the surface: $l_m \approx \kappa(z-d)$.  Take the first model level at
$z = 5\,\text{m}$, $d = 0$, $\kappa = 0.4$: $l_m \approx 2\,\text{m}$.  A model
with $\Delta z = 20\,\text{m}$ at that level has $\Delta z / l_m \approx 10$
— the grid is an order of magnitude coarser than the mixing length.  This is
precisely the regime where the $-\Delta z/(2z_k)$ shear bias is large.

The grid-curvature correction can be reframed as: **use the mixing-length-consistent
shear estimate** (exact log ratio) rather than the Taylor-approximated FD shear.
In pseudocode:

```
FUNCTION K_M_mixing_length(z_k, z_k1, d, kappa, lambda, L, u_star, Ri_g, Ri_c):

    z_eff = z_k - d

    -- Blackadar + Obukhov-limited mixing length (three-term reciprocal)
    inv_lm = 1.0 / (kappa * z_eff)
    IF lambda > 0.0 THEN
        inv_lm = inv_lm + 1.0 / lambda
    IF L > 0.0 THEN         -- stable
        inv_lm = inv_lm + 1.0 / (c_s * L)
    l_m = 1.0 / inv_lm

    -- Stability-reduced effective mixing length
    f_m_val  = f_m_stable(Ri_g, Ri_c)     -- uses §7.3 formula
    l_m_eff  = l_m * sqrt(f_m_val)

    -- Grid-aware shear (exact log ratio, no Taylor error)
    dz = z_k1 - z_k
    z_eff_k1 = z_k1 - d
    shear = (u_star / (kappa * dz)) * ln(z_eff_k1 / z_eff)

    -- Eddy diffusivity from mixing-length form
    K_M = l_m_eff^2 * shear

    -- Apply molecular floor
    K_M = max(K_M, nu)

    RETURN K_M, l_m, l_m_eff, shear
END FUNCTION
```

This pseudocode deliberately uses the exact log-ratio shear (§5.1), not the
native FD estimate, so the mixing-length form and the MOST form give the same
$K_M$ in the surface layer even on a coarse grid.

### 6.6 Mixing length above the surface layer

The surface-layer scheme hands off to the PBL scheme at height $h_{SL}$
(typically the top of the constant-flux layer, $\approx 0.1 h$).  The mixing
length should match continuously:

- At $z = h_{SL}$ the MOST form gives $l_m = \kappa(h_{SL} - d)/(1 + \kappa h_{SL}/\lambda)$.
- The PBL scheme should use the same $\lambda$ so there is no discontinuity in $K_M$.

In the McNider 1DBLM, the column mixing-length profile is computed from ground
to PBL top before any flux calculation.  The new surface-layer scheme should
expose $l_m(z_k)$ as a diagnostic output so the 1DBLM can use it directly
rather than recomputing it independently.

### 6.7 Summary: what $l_m$ adds to the scheme

| Without explicit $l_m$ | With explicit $l_m$ |
|---|---|
| $K_M = \kappa u_* z_{\text{eff}} f_m$ (MOST; valid only in surface layer) | $K_M = l_m^2 |S| f_m$ (extends smoothly above surface layer) |
| No awareness of BL height | $\lambda$ embeds $h$ information; $K_M$ grows then saturates |
| Stable: abrupt collapse at $Ri_c$ | Stable: additional $l_m$ suppression via $1/(c_s L)$ term |
| Grid shear error propagates to $K_M$ | Exact log-ratio shear makes $K_M$ grid-independent |
| No pathway to TKE coupling | $l_m$ is the natural link: $K_M = c_\mu l_m \sqrt{e}$ in TKE schemes |

---

## 7. Stability Function Branches (pseudocode)

### 6.1 Unstable — Businger–Dyer / Paulson

```
FUNCTION psi_m_unstable(zeta):
    -- zeta < 0
    x = (1 - 16*zeta)^(1/4)
    psi = 2*ln((1+x)/2) + ln((1+x^2)/2) - 2*arctan(x) + PI/2
    RETURN psi

FUNCTION psi_h_unstable(zeta):
    x = (1 - 16*zeta)^(1/4)    -- same x
    psi = 2*ln((1 + x^2)/2)    -- simpler form
    RETURN psi
```

Guard: if `|zeta| > 5` (strongly convective), flag as free-convection limit;
return a capped value and set diagnostic `FREE_CONVECTION_FLAG = TRUE`.
Wind-based flux estimates unreliable; model should switch to convective scaling.

### 6.2 Near-neutral — linear transition

```
FUNCTION psi_m_near_neutral(zeta):
    -- Linear interpolation between unstable and stable expressions
    -- at |zeta| < 0.05, both collapse to ≈ -β·zeta (stable) or +4·zeta (unstable)
    -- Use the appropriate sign:
    IF zeta >= 0 THEN
        RETURN -beta * zeta        -- stable linear: psi is negative (reduces wind)
    ELSE
        RETURN -4 * zeta           -- unstable linearized: psi is positive (increases wind)
    END IF
    -- Note: this is the SAME as the leading-order limit of the full expressions above
    -- No discontinuity at zeta = 0
```

### 6.3 Stable — England–McNider $f_m$ in $Ri$-space

```
FUNCTION f_m_stable(Ri_g, Ri_c):
    -- Ri_g in [0, Ri_c)
    IF Ri_g < 0 THEN
        RETURN 1.0                    -- should not be called in stable branch
    IF Ri_g >= Ri_c THEN
        RETURN 0.0                    -- laminar; handled below
    RETURN (1 - Ri_g / Ri_c)^2

FUNCTION f_h_stable(Ri_g, Ri_c, gamma, beta):
    -- gamma = stable φ_h slope, beta = stable φ_m slope
    -- For Pr_t = 1 (gamma = beta): f_h = f_m
    -- For Pr_t ≠ 1:
    Ri_c_h = 1 / gamma              -- effective Ri_c for heat
    IF Ri_g >= Ri_c_h THEN RETURN 0.0
    -- general form requires separate derivation; see hw/neutral-log-profiles.md Prob 8B(3)
    RETURN (1 - Ri_g / Ri_c_h) * (1 - Ri_g / Ri_c)
```

**Design question (open):** Should $Ri_c$ be static (= $1/\beta = 0.2$) or dynamic
(as in `Dynamic_Ric_Hybrid_MOST_Ri_Draft.md`)?  The dynamic $Ri_c$ framework is
a natural extension but should be a **separate module** that can be swapped in.
See Section 9.

### 6.4 Laminar floor

```
FUNCTION K_floor(z_eff, u_star=None):
    -- Enforce physical minimum: molecular values
    K_M_mol = nu                    -- kinematic viscosity ≈ 1.5e-5 m²/s
    K_H_mol = alpha                 -- thermal diffusivity ≈ 2.1e-5 m²/s
    RETURN K_M_mol, K_H_mol
```

The eddy diffusivities at the end of the scheme are:

```
    K_M = max(kappa * u_star * z_eff * f_m(Ri), K_M_mol)
    K_H = max(kappa * u_star * z_eff * f_h(Ri), K_H_mol)
```

This ensures `K_M > 0` always and prevents runaway cooling (Pitfall 15).

---

## 8. Iterative Flux Inversion (pseudocode)

Given observed (or model-state) bulk profiles $(U_1, U_2)$, $(\Theta_1, \Theta_2)$
at $(z_1, z_2)$, recover $u_*$, $\theta_*$, $L$ iteratively.

```
FUNCTION invert_fluxes(U1, U2, Th1, Th2, z1, z2, d, z0m, z0h, kappa, max_iter=50):

    -- Displaced effective heights
    ze1 = z1 - d;  ze2 = z2 - d
    dU  = U2 - U1;  dTh = Th2 - Th1

    -- Iteration state
    psi_m = 0;  psi_h = 0          -- neutral first guess
    converged = FALSE

    FOR i = 1 TO max_iter:

        -- Step A: estimate u_star, theta_star from corrected log ratio
        log_ratio_m = ln(ze2/ze1) - (psi_m(ze2/L) - psi_m(ze1/L))
        log_ratio_h = ln(ze2/ze1) - (psi_h(ze2/L) - psi_h(ze1/L))

        u_star     = kappa * dU / log_ratio_m
        theta_star = kappa * dTh / log_ratio_h

        -- Step B: guard against u_star → 0 (free-convection singularity)
        IF u_star < u_star_min THEN
            SET free_convection_flag = TRUE
            BREAK    -- wind-based inversion failed; use convective scaling
        END IF

        -- Step C: update Obukhov length
        L_new = - (u_star^3 * Th_bar) / (kappa * g * (-u_star * theta_star))
              = u_star^2 * Th_bar / (kappa * g * theta_star)

        -- Step D: update stability parameter and psi functions
        zeta1 = (z1 - d) / L_new
        zeta2 = (z2 - d) / L_new
        regime = classify_regime(...)
        psi_m = psi_m_for_regime(zeta1, regime)
        psi_h = psi_h_for_regime(zeta1, regime)

        -- Step E: convergence check
        IF |L_new - L_old| / max(|L_new|, 1.0) < tol THEN
            converged = TRUE;  BREAK
        END IF
        L_old = L_new

    END FOR

    IF NOT converged THEN
        WARN "flux inversion did not converge after max_iter steps"
        -- return last iterate; flag in diagnostics

    RETURN u_star, theta_star, L_new, regime, converged
END FUNCTION
```

**Near-$Ri_c$ guard:** before the iteration enters the stable branch, check
whether `Ri_g_corr > 0.9 * Ri_c`.  If so, the $\zeta = Ri/(1-\beta Ri)$
formula is near-singular (Pitfall 14) and the iteration step should be damped.

---

## 9. Displacement Height Consistency Contract

This is the most operationally important convention.  It is encoded here as a
formal contract:

```
CONTRACT: displacement_consistency

1. EVERY height argument passed to any profile, gradient, or flux function
   in this scheme is the GEOMETRIC height z (above ground).

2. The displaced height ze = z - d is computed INSIDE each core function,
   never passed in pre-computed.  This eliminates the possibility of mixing
   raw and displaced heights between callers.

3. The roughness lengths z0m and z0h are defined relative to z = d (i.e., the
   log argument is (z-d)/z0m, not z/z0m).

4. The validity domain check  ze > z0m  is performed at the top of every
   function that calls ln(ze/z0m).  Violation raises a runtime error (or
   returns a MASKED/MISSING value), not a silently wrong result.

5. Ri_g and Ri_B are computed using the displaced shear u*/(κ·ze), never u*/(κ·z).
   The constraint that Ri_g is height-independent in the log layer is only
   satisfied when ze is used consistently — see log-expansions.md §4.4.
```

---

## 10. Extension Points

These are design stubs for future work — **not in scope for Paper 2**, but the
architecture should not preclude them.

### 10.1 Dynamic $Ri_c$

The static $Ri_c = 1/\beta$ is a consequence of the linear MOST forms.  The
`Dynamic Ri_c` framework (`Dynamic_Ric_Hybrid_MOST_Ri_Draft.md`) allows $Ri_c$
to evolve with local shear, stratification, and wave activity.  This can be
supported by making `Ri_c` an **input parameter** to all stability functions
rather than a hardcoded constant.

### 10.2 Roughness Sublayer (RSL) Correction

For $z < 5h_c$, $\phi_m \neq 1$ even at neutral.  An RSL correction multiplies
the MOST-derived fluxes by a factor $\Phi_{\mathrm{RSL}}(z/h_c)$ that tapers to
1 above the RSL.  This is an additive module that wraps the surface-layer call.

### 10.3 Blending Height and Heterogeneous Surfaces

Over patchy terrain, an effective roughness length at the blending height
$z_b \approx 100z_0$ can replace $z_0$ for the grid-scale calculation.  The
parameterization interface should accept an `effective_z0` override.

### 10.4 Matched Asymptotic Profile (research)

From `log-expansions.md` §8, question 4: can we build a composite $\ln(z+H)$
expansion valid across all three height regimes (near-surface linear, transition,
far-field log)?  This would replace the piecewise log with a single smooth
analytic expression and eliminate the need for the regime-matching guard in the
profile function.  This is a research item — see `asymptotics.txt` for summary.

---

## 11. Relationship to Existing Code and Manuscripts

| Existing artifact | Role in this scheme |
|---|---|
| `src/rct/core/derivatives.py` | Supplies finite-difference shear; will be replaced/amended by §5.1 grid-aware shear |
| `src/rct/core/ri_estimators.py` | Existing Ri_g, Ri_B; will be wrapped with curvature correction from §4 |
| `implementations/McNider_1DBLM_fc_module.f90` | Target host model for integration; interface contract §3 is designed to match 1DBLM column call |
| `manuscripts/Grid_Curvature_SBL_v01.md` | Paper 1: curvature bias and neutral-preserving correction — this scheme implements and extends that theory |
| `drafts/log-expansions.md` | Provides the mathematical unification argument and the 17-pitfall table that motivates Sections 4, 5, and 8 of this scheme |
| `hw/neutral-log-profiles.md` | Analytic ground truth for test cases (§validation/) |
| `dynamic_Ric_strategy.md` | Design notes for Extension §9.1 |

---

## 12. Open Design Questions

1. **One function or branching?** Should the scheme be one monolithic function
   that branches internally on regime, or separate public functions per regime
   that the caller selects?  Monolithic is easier for NWP integration; branched
   is easier to test independently.

2. **What is $d$ for urban canopy in NWP?** Current operational models (WRF,
   MPAS) often set $d = 0$ globally.  Should the scheme accept $d = 0$ as a
   valid input, or require $d \geq 0$ with a warning when $d = 0$ and $z_0 > 0.1$?

3. **How to handle $z_k < d$?** The deepest model level may be below the
   canopy displacement height (e.g., a 2 m first level with a 3 m displacement).
   Option A: return missing/undefined and let the caller handle it.
   Option B: apply a canopy-submodel.  Currently Option A is preferred.

4. **Ri_c: 0.20, 0.25, or dynamic?** The static value depends on $\beta$.  The
   literature has $\beta \in [4.7, 6]$, giving $Ri_c \in [0.167, 0.213]$.
   Paper 2 should benchmark sensitivity to this choice explicitly.

5. **Coupling to TKE schemes:** The scheme as described is a first-order
   closure (K-theory).  In higher-order TKE-based schemes (MYNN, etc.), $K_M$
   is diagnosed from TKE rather than from $u_*$.  The $f_m(Ri)$ correction
   could be applied as a **modifier** to TKE-derived $K_M$, but the interaction
   requires careful thought to avoid double-counting.

6. **What asymptotic mixing length $\lambda$?** Three choices are in play:
   (a) fixed fraction of BL height $\lambda = c_\lambda h$ — requires $h$ from
   the PBL scheme; (b) Blackadar integral formula — self-contained but
   expensive; (c) Obukhov-length limited $\lambda = c_s L$ — naturally collapses
   under stable conditions but undefined at neutral.  A composite option:
   use (c) in the stable branch, (a) in the neutral/unstable branch, matched
   at $|\zeta| = 0.05$.  Which is most consistent with the McNider 1DBLM's
   existing mixing-length profile, and does the 1DBLM expose $h$ to the surface
   scheme?

---

*Next step: once design questions 1–3 and 6 are resolved, begin Paper 2 outline
(see `manuscripts/Paper2_GridAware_SL_Outline.md`).*
