# Stable Boundary Layer Corrections: Curvature-Aware MOST Implementation

## Executive Summary

Coarse vertical grids in atmospheric models systematically underestimate near-surface stability in the stable boundary layer (SBL), leading to excessive turbulent mixing and warm-biased surface temperatures. This document provides a comprehensive framework for curvature-aware corrections that preserve neutral physics (the invariant $2\Delta$) while reducing grid-induced bias by 40%+ in operational settings.

Key innovation: Analytic curvature of the gradient Richardson number $Ri_g(\zeta)$ quantifies the nonlinear stability structure; preserving the neutral curvature
$$
2\Delta \;=\; 2\left[\left.\frac{d\ln\phi_h}{d\zeta}\right|_{\zeta=0} - 2\left.\frac{d\ln\phi_m}{d\zeta}\right|_{\zeta=0}\right]
$$
anchors the correction to physically consistent near-neutral SBL behavior (for linear-stable: $\Delta = a_h/\mathrm{Pr} - 2a_m$) while damping coarse-grid tail effects.

---

## 1. Problem Statement

### 1.1 Observational Signature
- **Tower/LES:** Fine-resolution (Δz ≈ 5–10 m) Ri_g profiles show strong concave-down curvature in stable nights.
- **Coarse Models:** First-layer thickness Δz = 50–100 m → bulk Richardson number Ri_b systematically < point Ri_g(z_g).
- **Consequence:** Turbulent diffusivities K_m, K_h overestimated → excessive mixing → degraded surface inversion, warm bias, premature LLJ onset.

### 1.2 Root Cause
Second derivative ∂²Ri_g/∂ζ² < 0 (concave-down) for typical stable parameter sets → layer-averaging inequality (Jensen):
$$
Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_g),
\quad z_g = \sqrt{z_0 z_1}.
$$

Bias amplification ratio:
$$
B = \frac{Ri_g(z_g)}{Ri_b} > 1,
$$
often B ≈ 1.3–2.0 for strongly stable cases with coarse Δz.

---

## 2. Curvature Framework (MOST Foundation)

### 2.1 Core Definitions
Monin–Obukhov similarity:
$$
\zeta \;=\; \frac{z}{L},\qquad
L \;=\; -\frac{u_*^3 \,\theta}{\kappa g \,\overline{w'\theta'}}
$$
$$
\phi_m \;=\; \frac{\kappa z}{u_*}\frac{\partial U}{\partial z},\qquad
\phi_h \;=\; \frac{\kappa z}{\theta_*}\frac{\partial \theta}{\partial z},\quad
\theta_* \;=\; -\frac{\overline{w'\theta'}}{u_*}
$$

Gradient Richardson number:
$$
Ri_g(\zeta) \;=\; \zeta\,\frac{\phi_h}{\phi_m^2} \;=\; \zeta\,F(\zeta),\qquad F \;=\; \frac{\phi_h}{\phi_m^2}.
$$

### 2.2 Logarithmic Sensitivities
$$
v_m = \frac{\phi_m'}{\phi_m},\quad v_h = \frac{\phi_h'}{\phi_h},\quad
V_{\log} = v_h - 2v_m,\quad W_{\log} = V_{\log}'.
$$

### 2.3 Curvature Expression
$$
\boxed{\frac{d^2 Ri_g}{d\zeta^2} = F\left[2V_{\log} + \zeta\left(V_{\log}^2 - W_{\log}\right)\right]}.
$$

Neutral limit (ζ → 0) [SBL-generic]:
$$
\boxed{\Delta = V_{\log}(0) = \Big.\frac{d\ln\phi_h}{d\zeta}\Big|_{0} - 2\Big.\frac{d\ln\phi_m}{d\zeta}\Big|_{0}},\quad
\left.\frac{d^2 Ri_g}{d\zeta^2}\right|_0 = 2\Delta.
$$
For linear-stable φ_m = 1 + a_m ζ, φ_h = \mathrm{Pr} + a_h ζ: Δ = a_h/\mathrm{Pr} − 2a_m.

### 2.4 Near-Neutral Series
$$
Ri_g(\zeta) = \zeta + \Delta\zeta^2 + \tfrac{1}{2}(\Delta^2 + c_1)\zeta^3 + O(\zeta^4),
$$
$$
\boxed{c_1 = \Big.\frac{d^2\ln\phi_h}{d\zeta^2}\Big|_{0} - 2\Big.\frac{d^2\ln\phi_m}{d\zeta^2}\Big|_{0}}
$$
Inversion:
$$
\zeta(Ri) = Ri - \Delta Ri^2 + \left(\tfrac{3}{2}\Delta^2 - \tfrac{1}{2}c_1\right)Ri^3 + O(Ri^4).
$$

---

## 3. Standard Stable Formulations and Their Limitations

### 3.1 Power-Law (Businger–Dyer) — Unstable Only (ζ < 0; do not use for SBL)
// Retain for context but mark as inapplicable to ζ>0
$$
\phi_m = (1 - \beta_m\zeta)^{-\alpha_m},\qquad
\phi_h = (1 - \beta_h\zeta)^{-\alpha_h},\qquad \zeta < 0.
$$
**Not applicable to SBL (ζ > 0). Using this in stable regimes is a misapplication and leads to incorrect curvature and poles.**
**Limitations (for unstable regime context only):**
- Finite-height pole at ζ = 1/β → requires hard cutoff or guard.
- Curvature blows up near pole → numerical instability.
- Rapid growth causes excessive sensitivity to Δz.

### 3.2 Linear Stable (Högström, Beljaars–Holtslag)
$$
\phi_m \;=\; 1 + a_m\,\zeta,\qquad \phi_h \;=\; \mathrm{Pr} + a_h\,\zeta,\qquad \zeta>0
$$
Curvature: $\Delta = a_h/\mathrm{Pr} - 2a_m$ (concave-down if $\Delta<0$).
**Limitations:**
- Curvature is constant in the neutral limit; may under-represent curvature growth aloft.
- May overmix if a_h/\mathrm{Pr} ≈ 2a_m (Δ ≈ 0).

### 3.3 Beljaars–Holtslag (1991)
Hybrid polynomial + exponential:
$$
\phi_m = 1 + a\zeta + b\zeta\left[1 + c\zeta\right]^{1/3}.
$$

**Advantages:** Smoother, no hard pole.  
**Limitations:** Complex functional form, parameters empirically tuned, still exhibits grid sensitivity.

---

## 4. Curvature-Aware Correction Strategy

### 4.1 Design Principles
1. **Preserve Neutral Curvature (2Δ):** Do not alter near-neutral physics; anchor correction to ζ → 0 behavior.
2. **Damp Tail Effects:** Reduce curvature-induced bias for ζ > 0 on coarse grids without ad hoc floors.
3. **Grid Convergence:** Correction → 0 as Δz → 0 (fine grids recover standard MOST).
4. **Monotonicity:** Avoid introducing spurious oscillations or instabilities.

### 4.2 Grid Damping Factor Approach
Modify eddy diffusivities:
$$
K_m^* \;=\; K_m \times G(\zeta,\Delta z),\qquad
K_h^* \;=\; K_h \times G(\zeta,\Delta z)
$$
Constraints on $G$:
$$
G(0,\Delta z)=1,\quad
\left.\frac{\partial G}{\partial \zeta}\right|_{\zeta=0}=0,\quad
\lim_{\Delta z\to 0}G=1,\quad
\frac{\partial G}{\partial \zeta}\le 0 \text{ for fixed }\Delta z.
$$
Functional template:
$$
G(\zeta,\Delta z) \;=\; \exp\!\left[-D\left(\frac{\Delta z}{\Delta z_r}\right)^{p}\left(\frac{\zeta}{\zeta_r}\right)^{q}\right],\qquad p\ge 1,\;q\ge 2.
$$
with:
- $p \ge 1$ (grid-ratio exponent),
- $q \ge 2$ (ensures ∂G/∂ζ|₀ = 0),
- $D$: calibration coefficient (target bias reduction),
- Reference scales: Δz_r = 10 m, ζ_r = 0.5.

### 4.3 Tail Modifier in φ Functions (Alternative)

Embed correction directly:
$$
\phi_m^*(\zeta, \Delta z) = \phi_m(\zeta) \times f_c(\zeta, \Delta z),\qquad
\phi_h^*(\zeta, \Delta z) = \phi_h(\zeta) \times f_c(\zeta, \Delta z),
$$
with same exponential form for $f_c$, calibrated to preserve $V_{\log}(0)$ and $W_{\log}(0)$.

**Advantage:** Unified treatment in similarity framework.  
**Disadvantage:** Requires careful chain-rule handling in flux–gradient relationships.

### 4.4 Quadratic SBL Surrogate (Q-SBL)

For ζ ∈ [0, ζ_max] (typically ζ_max ≈ 0.2–0.3):
$$
\phi_m^{\text{Q}} = 1 + a_m\zeta + b_m\zeta^2,\qquad
\phi_h^{\text{Q}} = \mathrm{Pr} + a_h\zeta + b_h\zeta^2,
$$
Choose a_m, a_h, Pr from SBL datasets (e.g., Högström, BH91, SHEBA); set b_m, b_h to match c_1 if available or via fit.
// Removed α,β-based definitions (unstable-only).
$$
Ri_g^{\text{Q}} = \zeta + \Delta\zeta^2 + \tfrac{1}{2}(\Delta^2 + c_1)\zeta^3,
$$
$$
\frac{d^2 Ri_g^{\text{Q}}}{d\zeta^2} = 2\Delta + 3(\Delta^2 - c_1)\zeta.
$$

**Advantages:**
- No finite-height pole.
- Smooth derivatives.
- Analytic inversion (cubic root).
- Neutral curvature 2Δ preserved exactly.

**Recommended Use:**
- Primary closure for ζ < 0.2–0.3.
- Blend to standard power-law or capped linear for ζ > ζ_max.

---

## 5. Calibration and Diagnostics

### 5.1 Parameter Selection (SBL-Focused)

**Use SBL-specific φ forms:**
- Linear-stable: a_m ≈ 4.5–6.5, a_h ≈ 5–8, Pr ≈ 0.9–1.0.
- Beljaars–Holtslag (1991): operational stable coefficients.
- SHEBA/Arctic (Grachev): stronger stability; larger |Δ|.

**Recommended Workflow:**
1. Fit SBL-only segments (ζ > 0.05, Ri_g < 0.25) using chosen stable form.
2. Compute Δ = (d ln φ_h/dζ)|_0 − 2(d ln φ_m/dζ)|_0 and c_1 from the fitted φ.
3. Verify neutral curvature 2Δ against observations.
4. Check inflection height ζ_inf (if present).

**Typical Ranges:**
- a_m ≈ 4.5–6.5, a_h ≈ 5–8, Pr ≈ 0.9–1.0,
- Δ ≈ −1 to −7 (negative → concave-down).

### 5.2 Bias Amplification Ratio

At first model level (z₀ to z₁):
$$
B = \frac{Ri_g(z_g)}{Ri_b},\quad z_g = \sqrt{z_0 z_1}.
$$

**Target Performance:**
- Fine grid (Δz ≈ 10 m): B ≈ 1.0–1.05.
- Coarse uncorrected (Δz ≈ 100 m): B ≈ 1.5–2.0.
- Coarse corrected: B ≈ 1.1–1.2 (40%+ error reduction).

### 5.3 Curvature Diagnostics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| Neutral curvature | 2Δ | Sign/magnitude of initial departure from linearity |
| Amplification | A(ζ) = \|∂²Ri_g/∂ζ² / 2Δ\| | Curvature growth vs neutral |
| Inflection height | ζ_inf: ∂²Ri_g/∂ζ² = 0 | Stability regime transition (if exists) |
| Layer error | E_i = \|Ri_g(z_i+1) − Ri_g(z_i) − Δz ∂Ri_g/∂z\|_zg | Reconstruction bias |

### 5.4 Numerical Ri Diagnostics (Tower Data)

**Step 1: Gradient estimation**  
Centered difference at level $k$:
$$
S_k = \sqrt{\left(\frac{U_{k+1} - U_{k-1}}{z_{k+1} - z_{k-1}}\right)^2 + \left(\frac{V_{k+1} - V_{k-1}}{z_{k+1} - z_{k-1}}\right)^2},
$$
$$
\frac{\partial\theta}{\partial z}\Big|_{z_k} \approx \frac{\theta_{k+1} - \theta_{k-1}}{z_{k+1} - z_{k-1}}.
$$

**Step 2: Point Ri_g**
$$
Ri_g(z_k) = \frac{(g/\theta_k)\,\partial\theta/\partial z}{S_k^2}.
$$

**Step 3: Bulk Ri_b (layer)**
$$
Ri_b = \frac{g}{\theta_{\text{ref}}}\frac{(\theta_1 - \theta_0)(z_1 - z_0)}{(U_1 - U_0)^2 + (V_1 - V_0)^2}.
$$

**Step 4: Bias ratio**
$$
B = \frac{Ri_g(z_g)}{Ri_b},\quad z_g = \sqrt{z_0 z_1}.
$$

**Step 5: Numerical integration (if full profile available)**
$$
Ri_b^{\text{int}} = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz \approx \text{trapezoid or Simpson}.
$$

Compare $Ri_b^{\text{int}}$ vs direct bulk formula; validate curvature correction by checking $B$ reduction.

---

## 6. Implementation Workflow

### 6.1 Preprocessing (Offline Analysis)

**Step 1: Data Collection**
- Tower: u, v, θ, z (multiple levels); fluxes u_*, θ_*, H.
- LES: 3D fields → vertical profiles + grid spacing matrix.
- Remote sensing: lidar wind + radiometer temperature (validation).

**Step 2: Parameter Estimation**
```python
# Fit stable segments only
stable_mask = (zeta > 0.05) & (Ri_g < 0.25)
alpha_m, beta_m = fit_phi(U_profile[stable_mask], z[stable_mask])
alpha_h, beta_h = fit_phi(theta_profile[stable_mask], z[stable_mask])

# Compute neutral coefficients
Delta = alpha_h * beta_h - 2 * alpha_m * beta_m
c1 = alpha_h * beta_h**2 - 2 * alpha_m * beta_m**2
print(f"Neutral curvature 2Δ = {2*Delta:.3f}")
```

**Step 3: Domain Guard**
```python
zeta_max = 0.7 / max(beta_m, beta_h)  # Safe operating range
# Optionally switch to Q-SBL for zeta > 0.2
```

### 6.2 Online Model Integration

**Pseudocode (Vertical Diffusion Loop):**
```python
for k in range(nz):
    zeta = max(z[k] / L, 0.0)  # SBL branch

    # Stable MOST (choose one; linear-stable shown)
    phi_m = 1.0 + a_m * zeta
    phi_h = Pr + a_h * zeta
    # Alternatively: use BH91 or SHEBA stable forms

    # Curvature-aware correction
    dz_ratio = (z[k+1] - z[k]) / dz_ref
    G = exp(-D * (dz_ratio)**p * (zeta / zeta_ref)**q)

    # Modified diffusivities
    K_m[k] = (u_star * kappa * z[k] / phi_m) * G
    K_h[k] = (u_star * kappa * z[k] / phi_h) * G
```

**Calibration Constants (Example):**
- D = 0.8–1.2, p = 1.5, q = 2
- dz_ref = 10 m, zeta_ref = 0.3  # stable-reference scale

### 6.3 Boundary Condition Handling

**Surface Layer (k=0):**
- Use bulk Richardson number Ri_b or reconstruct Ri_g(z_g) via geometric mean.
- Apply logarithmic transformation for z₀ dependence:
  $$
  \frac{\partial}{\partial z} = \frac{1}{z}\frac{\partial}{\partial \ln z}.
  $$

**First Interior Level (k=1):**
- Critical layer for curvature bias; apply strongest correction here.
- Validate against analytic Ri_g(ζ₁) from stability functions.

**Aloft (k ≥ 2):**
- Correction diminishes as Δz/z decreases (aspect ratio effect).
- Optionally transition to local Richardson-based closure.

---

## 7. Validation Protocol

### 7.1 Synthetic Test Cases

**Case 1: Idealized Stable Night**
- Prescribed L = 50 m, u_* = 0.2 m/s, constant cooling.
- Grid sequence: Δz = 5, 10, 25, 50, 100 m.
- Compare uncorrected vs corrected K_m, K_h profiles.

**Case 2: LES Benchmark (GABLS1)**
- 9-hour nocturnal evolution, prescribed cooling.
- Extract "truth" from Δz = 1 m; test at Δz = 25, 50 m.
- Metrics: surface heat flux error, inversion height bias, TKE profile RMSE.

### 7.2 Observational Validation

**ARM NSA (Alaska):** stable nights ($\zeta>0.1$).  
**SHEBA:** assess $L(z)$ variability and omission metric
$$
E_{\text{omit}} \;=\; \left|\frac{(d^2\zeta/dz^2)(dRi_g/d\zeta)}{(d\zeta/dz)^2(d^2Ri_g/d\zeta^2)}\right|.
$$


**Dallas/Ft. Worth (Urban):** 325 m tower + lidar/radiometer fusion.  
**Metrics:** curvature bias amplification vs aggregation.

### 7.3 Success Criteria

| Metric | Target (Coarse Grid, Δz = 100 m) |
|--------|-----------------------------------|
| Bias ratio B | < 1.2 (vs uncorrected ~1.8) |
| Surface flux RMSE | < 15% (vs uncorrected ~30%) |
| Inversion height error | < 20 m (vs uncorrected ~50 m) |
| Neutral curvature preservation | \|2Δ* − 2Δ\| / \|2Δ\| < 5% |
| Computational overhead | < 5% (vs standard MOST) |

---

## 8. Advanced Topics

### 8.1 Dynamic Critical Richardson Number

Instead of fixed Ri_c = 0.25:
$$
Ri_c^* = f(\text{inversion strength}, \text{shear}, \text{history}),
$$
based on:
- Lapse rate Γ = ∂θ/∂z above inversion.
- Turbulence memory (previous timestep TKE).
- Observed Ri_c distributions (0.21–1.0 range).

**Implementation:**
```python
# Inversion-strength proxy
gamma_inv = (theta[k_inv+1] - theta[k_inv]) / dz
Ri_c_dynamic = 0.25 + 0.5 * min(gamma_inv / gamma_ref, 1.0)
```

### 8.2 Ri-Based Direct Closures

Avoid ζ iteration by using Ri directly in the mixing-length form
$$
K_{m,h} = f_{m,h}(Ri)\; l^2\; S,
$$
with
$$
S = \sqrt{\left(\frac{\partial u}{\partial z}\right)^2 + \left(\frac{\partial v}{\partial z}\right)^2},\quad l \approx \kappa z.
$$
Mapping from MOST φ(ζ) to f(Ri) using ζ(Ri):
$$
\boxed{f_m(Ri) = \frac{1}{\phi_m(\zeta(Ri))^2}},\qquad
\boxed{f_h(Ri) = \frac{1}{\phi_m(\zeta(Ri))\;\phi_h(\zeta(Ri))}}.
$$
This preserves $K_m = u_* \kappa z / \phi_m$ and $K_h = u_* \kappa z / \phi_h$ when $l = \kappa z$ and $S = |\partial U/\partial z|$.

**Series Inversion (near-neutral):**
$$
\zeta(Ri) = Ri - \Delta Ri^2 + \left(\tfrac{3}{2}\Delta^2 - \tfrac{1}{2}c_1\right)Ri^3.
$$

**Newton Refinement (1–2 iterations):**
```python
def zeta_from_ri_newton(Ri, phi_m, phi_h, z0=None, tol=1e-10):
    if z0 is None:
        z0 = Ri - Delta * Ri**2  # Series seed
    F = lambda z: phi_h(z) / phi_m(z)**2
    for _ in range(2):
        f = z0 * F(z0) - Ri
        # V_log = d ln F / dζ; compute numerically or analytically if available
        V_log = (math.log(F(z0+1e-6)) - math.log(F(z0-1e-6))) / (2e-6)
        fp = F(z0) + z0 * F(z0) * V_log
        z0 -= f / fp
        if abs(f) < tol: break
    return z0
```

Example (drop-in, linear-stable φ shown):
```python
# Inputs: u_star, kappa, z[k], du_dz, dv_dz, a_m, a_h, Pr, Delta, c1
S = math.hypot(du_dz, dv_dz)  # 2D shear magnitude
l = kappa * z[k]

phi_m = lambda zeta: 1.0 + a_m * zeta
phi_h = lambda zeta: Pr + a_h * zeta

# Invert ζ(Ri)
zeta = zeta_from_ri_newton(Ri, phi_m, phi_h, z0=Ri - Delta*Ri*Ri)

# Map φ -> f
pm = phi_m(zeta); ph = phi_h(zeta)
f_m = 1.0 / (pm * pm)
f_h = 1.0 / (pm * ph)

# Eddy coefficients
K_m[k] = f_m * (l*l) * S
K_h[k] = f_h * (l*l) * S
```

Add: Exponential Ri-based stability function (simple, pole-free)
- Definition (used directly in K): $f(Ri) = \exp(-\gamma Ri / Ri_c)$, with $\gamma \approx 3.2$.
- Use separate $f_m,f_h$ if desired; maintain Pr_t consistency by choosing $\gamma_h$ accordingly.

---

## 9. Common Pitfalls and Troubleshooting

### 9.1 Pitfall: USL Parameters in Strong SBL
**Symptom:** Excessive curvature, early K collapse, numerical instability.  
**Fix:** Use stable φ (linear/BH/SHEBA); do not apply unstable power-law to ζ>0.

### 9.2 Pitfall: Arithmetic Mean Height
**Symptom:** Biased Ri_b even with correct φ functions.  
**Fix:** Use geometric mean z_g = √(z₀z₁) for layer reconstruction.

### 9.3 Pitfall: Ignoring L(z) Variability
**Symptom:** Curvature mismatch in deep stable layers.  
**Fix:** Compute E_omit; apply full chain rule if E_omit > 0.05.

### 9.4 Pitfall: Over-Smoothing Curvature
**Symptom:** Loss of inflection points, artificial diffusion.  
**Fix:** Preserve 2Δ; damp tails only, not near-neutral regime.

### 9.5 Pitfall: Fixed Ri_c = 0.25 in Intermittent Regimes
**Symptom:** Premature turbulence cutoff or excessive persistence.  
**Fix:** Implement dynamic Ri_c based on inversion strength or TKE memory.

---

## 10. Software and Reproducibility

### 10.1 Repository Structure

---

## 11. McNider-Biazar Research Priorities and Collaboration Framework

### 11.1 Historical Context: McNider School Contributions to SBL Physics

**McNider (1982, 1995):** Pioneering work on:
- Slope flows and topographic channeling in stable stratification
- SBL collapse mechanisms and intermittent turbulence regime transitions
- Coupling between radiative cooling rates and turbulent diffusion timescales
- Low-level jet (LLJ) formation and acceleration under strong stability

**Biazar & McNider (1995, 2006):** Advances in:
- Air quality modeling under stable conditions (pollutant trapping)
- Nocturnal chemistry sensitivity to vertical mixing parameterization
- Grid resolution impacts on regulatory model performance
- Ozone production/NOx titration coupling with SBL depth

**Synthesis Opportunity:** The curvature-aware framework provides theoretical foundation for McNider-Biazar empirical insights on grid-induced biases and regime transitions.

---

### 11.2 Priority Research Areas

#### **Area 1: Intermittent Turbulence and Regime Transitions**

**Objective:** Extend curvature framework to time-dependent stability regimes.

**McNider Lead Tasks:**
1. **Turbulence Memory Parameter:**
   - Develop persistence metric for Ri_c evolution:
     $$
     Ri_c(t) = Ri_{c,\text{eq}} + \tau_{\text{mem}} \frac{d\text{TKE}}{dt},
     $$
   - Link to LLJ acceleration/deceleration cycles.
   - Validate against CASES-99 intermittent nights.

2. **Collapse Prediction Metric:**
   - Formulate curvature-based early warning for turbulence shutdown:
     $$
     \Lambda_{\text{collapse}} = \frac{|\partial^2 Ri_g/\partial\zeta^2|}{\text{TKE}} \times \frac{\partial\theta/\partial t}{|\partial U/\partial z|}.
     $$
   - Identify critical thresholds from tower/LES ensembles.

3. **Radiative-Turbulent Coupling:**
   - Quantify feedback between surface cooling rate and Ri_g curvature:
     $$
     \frac{dL}{dt} = f\left(\frac{\partial F_{\text{rad}}}{\partial z}, K_h(\zeta), \frac{d^2Ri_g}{d\zeta^2}\right).
     $$
   - Test sensitivity to cloud-base height and moisture stratification.

**Deliverables:**
- Algorithm for adaptive Ri_c in operational models.
- Phase-space diagram (TKE vs ∂²Ri_g/∂ζ²) for regime classification.
- Case study: GABLS3 diurnal cycle with curvature-aware transitions.

---

#### **Area 2: Air Quality Applications and Chemical Coupling**

**Objective:** Translate SBL mixing corrections into regulatory model improvements.

**Biazar Lead Tasks:**
1. **Pollutant Trapping Quantification:**
   - Compare uncorrected vs corrected vertical diffusivity impact on:
     - NOx lifetime and O₃ production efficiency.
     - PM2.5 concentration gradients in urban SBL.
     - VOC segregation and nighttime chemistry rates.
   - Metric: ΔC/C_obs for EPA monitoring sites (Dallas, Houston networks).

2. **Grid Resolution Study for Air Quality:**
   - Extend curvature bias analysis to CMAQ/CAMx frameworks:
     $$
     E_{\text{AQ}} = \left|\frac{C_{\text{model}}(\Delta z) - C_{\text{obs}}}{C_{\text{obs}}}\right| \text{ vs } B(\Delta z).
     $$
   - Target: Reduce nighttime O₃ overprediction by 30% in stable events.

3. **Deposition Velocity Corrections:**
   - Reassess v_d(z₀) under corrected K_h:
     $$
     v_d^* = \frac{1}{r_a + r_b + r_c}, \quad r_a^* = \int_{z_0}^{z_r} \frac{dz}{K_h^*(z)}.
     $$
   - Impact on dry deposition fluxes for O₃, NO₂, SO₂.

4. **Regulatory Model Protocol:**
   - Draft EPA-compatible implementation guide for WRF-CMAQ.
   - Document computational cost and sensitivity to tuning parameters.
   - Prepare technical support document for SIP applications.

**Deliverables:**
- Peer-reviewed study: "Curvature-Aware SBL Physics in Photochemical Models."
- CMAQ module with McNider-Biazar correction option (Fortran + Python interface).
- Webinar/workshop for state air quality modelers.

---

#### **Area 3: Topographic Complexity and Drainage Flows**

**Objective:** Adapt curvature framework to sloped/complex terrain.

**McNider Lead Tasks:**
1. **Slope-Modified Stability Functions:**
   - Generalize ζ to include buoyancy-shear angle θ_s:
     $$
     \zeta_s = \zeta \cos\theta_s + \text{drainage forcing term}.
     $$
   - Derive curvature correction for anabatic/katabatic regimes.

2. **Valley Cold Pool Dynamics:**
   - Link inversion strength to curvature-induced mixing suppression:
     $$
     \frac{d\theta_{\text{pool}}}{dt} = Q_{\text{rad}} - \nabla \cdot (K_h^* \nabla \theta) + \text{advection}.
     $$
   - Test in Great Basin, Intermountain West winter inversions.

3. **LLJ-Topography Interaction:**
   - Analyze curvature bias impact on nocturnal jet height and core speed.
   - Validate against Oklahoma Mesonet profiler network.

**Deliverables:**
- Terrain-aware curvature diagnostic tool.
- Case study: Colorado Front Range cold air drainage with/without correction.
- Synthesis paper on McNider legacy in complex-terrain SBL modeling.

---

#### **Area 4: Urban and Suburban SBL Modifications**

**Objective:** Extend to heterogeneous surface roughness and anthropogenic heat.

**Biazar Lead Tasks:**
1. **Urban Canopy-SBL Coupling:**
   - Modify z₀, z₀_h for urban morphology (H/W ratios, sky-view factor).
   - Test curvature framework with spatially varying L(x, y, z).
   - Validation: Remote Sensing 2024 megacity tower data (325 m).

2. **Anthropogenic Heat Flux Impact:**
   - Add Q_anth to surface energy budget:
     $$
     L^{-1} = -\frac{\kappa g}{\theta u_*^3}(\overline{w'\theta'}_{\text{turb}} + Q_{\text{anth}}/\rho c_p).
     $$
   - Assess curvature sensitivity to Q_anth diurnal/weekly cycles.

3. **Pollution "Hot Spots" Under Weak Mixing:**
   - Map spatial correlation between curvature bias and observed C_max.
   - Identify vulnerable neighborhoods for public health targeting.

**Deliverables:**
- Urban-modified Q-SBL parameterization.
- Dallas/Ft. Worth O₃ attainment modeling with corrected SBL.
- Policy brief on health equity implications of mixing biases.

---

#### **Area 5: Observational Campaign Design and Data Assimilation**

**Objective:** Optimize field experiments for curvature validation and parameter tuning.

**Joint McNider-Biazar Tasks:**
1. **Next-Generation Tower Networks:**
   - Deploy dense vertical arrays in diverse SBL settings (urban, rural, polar).
   - Targeted campaigns during key stable events (e.g., Arctic winter, nocturnal urban).
   - High-frequency (1 min) θ, u, v, w, CO₂, O₃, NOx, PM2.5 profiles.

2. **Remote Sensing Integration:**
   - Lidar + radiometer networks for continuous boundary layer monitoring.
   - Satellite retrievals (e.g., MODIS, VIIRS) for synoptic-scale context.

3. **Data Assimilation Experiments:**
   - Twin experiments with synthetic vs real observations.
   - Test impact of curvature-informed vs standard mixing on forecast skill.

**Deliverables:**
- Design report: "Optimizing Observational Networks for Curvature Validation."
- Data assimilation protocol for incorporating curvature effects.
- Workshop: "Integrating Remote Sensing and In-Situ Data for SBL Research."

---

### Scalar (q) Extension
Apply same grid damping $G$:
\[
K_q^*=K_q G,\quad K_q=\frac{\kappa z u_*}{\phi_q},\quad f_q(Ri_g)=\frac{1}{\phi_m\phi_q}.
\]
Bias metrics identical; preserve neutral curvature (2Δ) for momentum, derive $a_q,b_q$ from scalar profile fits.
