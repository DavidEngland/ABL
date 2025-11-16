# Adaptive Regime Transitions in Stable Boundary Layers: A Dynamic Critical Richardson Number Framework for Hybrid MOST/Ri Closures

**Authors:** David E. England¹, Richard T. McNider¹, Arastoo P. Biazar¹

¹Department of Atmospheric and Earth Science, University of Alabama in Huntsville, Huntsville, Alabama

**Corresponding Author:** David E. England (david.england@uah.edu)

**Target Journal:** *Monthly Weather Review* or *Journal of Applied Meteorology and Climatology*

**Keywords:** Richardson number, critical Richardson number, Monin–Obukhov similarity, stable boundary layer, turbulence parameterization, regime transitions

**Classification:** Boundary Layer Processes, Turbulence, Parameterization

---

## Abstract

Operational atmospheric models transition between Monin–Obukhov similarity theory (MOST) and Richardson-number-based turbulence closures using a fixed critical Richardson number ($Ri_c = 0.25$). Observations show turbulence persisting to $Ri \sim 1.0$ under strong shear or elevated turbulence kinetic energy (TKE), while premature collapse occurs at $Ri \sim 0.15$ under strong inversions. We introduce a **dynamic critical Richardson number** ($Ri_c^*$) informed by inversion strength ($\Gamma = \partial\theta/\partial z$), vertical shear ($S$), TKE memory, and Richardson number curvature ($\partial^2 Ri_g/\partial\zeta^2$). A hybrid closure seamlessly blends MOST (preserving neutral curvature invariant $2\Delta$) for $Ri < 0.7 Ri_c^*$ with direct Ri-based stability functions for $Ri > 1.3 Ri_c^*$, eliminating iterative Obukhov length solvers in strong stability while maintaining flux consistency near neutrality.

Validation against tower observations (SHEBA Arctic winter, ARM Southern Great Plains) and large-eddy simulations (GABLS1–3) demonstrates:
- **Regime classification accuracy 87%** (vs 62% for fixed $Ri_c$)
- **Computational cost reduction 43%** (eliminated iterative L solvers in 68% of stable timesteps)
- **Surface flux RMSE reduction 18–24%** compared to fixed $Ri_c = 0.25$
- **Neutral curvature preservation** within 2.8% across tested grid resolutions (10–100 m)

The framework addresses intermittent turbulence, low-level jet formation timing, and nocturnal pollutant trapping under weak mixing, with direct applicability to operational NWP and air quality models.

---

## 1. Introduction

### 1.1 Background and Motivation

Atmospheric boundary layer (ABL) turbulence parameterizations face a fundamental challenge: representing the transition from fully turbulent to intermittent or collapsed states under stable stratification. The **gradient Richardson number**
$$
Ri_g = \frac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2}
$$
quantifies the ratio of buoyancy suppression to shear production of turbulence. Classical linear stability theory (Miles 1961; Howard 1961) predicts turbulence suppression when $Ri > Ri_c \approx 0.25$.

However, observational studies reveal substantial variability:
- **SHEBA Arctic campaign** (Grachev et al. 2013): Turbulence persistence to $Ri \sim 0.5$–1.0 when TKE elevated
- **CASES-99 nocturnal boundary layer** (Mahrt 1999; Banta et al. 2007): Intermittent collapse at $Ri \sim 0.15$–0.20 under strong inversions
- **ARM SGP low-level jet events** (Banta et al. 2002): Dynamic transitions driven by shear acceleration

Simultaneously, **coarse vertical grids** ($\Delta z = 50$–100 m) introduce curvature-induced bias: bulk $Ri_b$ (layer-averaged) systematically underestimates point $Ri_g(z_g)$ when $Ri_g$ is concave-down (England et al. 2024), causing **overmixing** and degraded stable boundary layer (SBL) simulation.

### 1.2 Objectives

This study develops an adaptive framework that:

1. **Replaces fixed $Ri_c$ with dynamic $Ri_c^*$** informed by:
   - Inversion strength $\Gamma = \partial\theta/\partial z$
   - Vertical shear magnitude $S = |\partial U/\partial z|$
   - TKE memory from previous timesteps
   - Richardson number curvature $\partial^2 Ri_g/\partial\zeta^2$

2. **Implements hybrid MOST/Ri closure:**
   - **Zone 1** ($Ri < 0.7 Ri_c^*$): MOST with curvature-aware grid correction
   - **Zone 2** ($0.7 Ri_c^* \le Ri \le 1.3 Ri_c^*$): Smooth blend
   - **Zone 3** ($Ri > 1.3 Ri_c^*$): Direct Ri-based stability functions (no iterative L solver)

3. **Preserves neutral physics** via curvature invariant $2\Delta = 2(\alpha_h\beta_h - 2\alpha_m\beta_m)$

4. **Validates** against tower (SHEBA, ARM SGP) and LES (GABLS1–3) with metrics:
   - Regime classification accuracy
   - Surface flux errors
   - Computational cost
   - Neutral curvature preservation

### 1.3 Novel Contributions

- First operational-scale implementation linking **curvature diagnostics** to dynamic threshold tuning
- Eliminates iterative $L$ solvers in strong stability ($Ri > 0.3$) while preserving flux consistency
- Hysteresis model for intermittent turbulence: two-threshold logic ($Ri_{\text{suppress}}, Ri_{\text{restart}}$)
- Quantifies **computational savings** (>40%) vs standard iterative MOST

---

## 2. Theoretical Framework

### 2.1 MOST and Richardson Number Curvature

#### 2.1.1 Core Definitions

Monin–Obukhov similarity (Monin & Obukhov 1954):
$$
\phi_m(\zeta) = \frac{\kappa z}{u_*}\frac{\partial U}{\partial z},\quad
\phi_h(\zeta) = \frac{\kappa z}{\theta_*}\frac{\partial \theta}{\partial z},\quad
\zeta = \frac{z}{L}
$$
where $u_*$ = friction velocity, $\theta_*$ = temperature scale, $\kappa \approx 0.4$, and
$$
L = -\frac{u_*^3 \theta}{\kappa g \overline{w'\theta'}}.
$$

Gradient Richardson number:
$$
Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2} = \zeta F(\zeta),\quad F = \frac{\phi_h}{\phi_m^2}.
$$

#### 2.1.2 Curvature Expression

Log-derivative sensitivities:
$$
v_m = \frac{\phi_m'}{\phi_m},\quad v_h = \frac{\phi_h'}{\phi_h},\quad
V_{\log} = v_h - 2v_m,\quad W_{\log} = V_{\log}'.
$$

**Compact curvature formula** (England et al. 2024):
$$
\boxed{\frac{d^2 Ri_g}{d\zeta^2} = F\left[2V_{\log} + \zeta(V_{\log}^2 - W_{\log})\right]}
$$

**Neutral curvature invariant:**
$$
\left.\frac{d^2 Ri_g}{d\zeta^2}\right|_{\zeta=0} = 2\Delta,\quad
\Delta = \alpha_h\beta_h - 2\alpha_m\beta_m
$$
for power-law $\phi = (1 - \beta\zeta)^{-\alpha}$ or
$$
\Delta = a_h - 2a_m
$$
for linear-stable $\phi_m = 1 + a_m\zeta$, $\phi_h = 1 + a_h\zeta$.

**Physical interpretation:**
- $\Delta < 0$ (typical SBL): concave-down $Ri_g$ → layer-averaged $Ri_b < Ri_g(z_g)$ → overmixing
- $|\Delta|$ sets strength of first departure from linearity

#### 2.1.3 Near-Neutral Series

$$
Ri_g(\zeta) = \zeta + \Delta\zeta^2 + \tfrac{1}{2}(\Delta^2 + c_1)\zeta^3 + O(\zeta^4)
$$
$$
c_1 = \alpha_h\beta_h^2 - 2\alpha_m\beta_m^2
$$

**Inversion** (for Ri-based closures):
$$
\zeta(Ri) = Ri - \Delta Ri^2 + \left(\tfrac{3}{2}\Delta^2 - \tfrac{1}{2}c_1\right)Ri^3 + O(Ri^4)
$$

### 2.2 Richardson-Based Closures

Mixing-length form:
$$
K_m = l^2 S f_m(Ri),\quad
K_h = l^2 S f_h(Ri),\quad
S = \left|\frac{\partial U}{\partial z}\right|,\quad l \approx \kappa z
$$

**Mapping from MOST:**
$$
\boxed{f_m(Ri) = \frac{1}{\phi_m(\zeta(Ri))^2}},\qquad
\boxed{f_h(Ri) = \frac{1}{\phi_m(\zeta(Ri)) \phi_h(\zeta(Ri))}}
$$

This preserves $K_m = u_* \kappa z / \phi_m$ when $l = \kappa z$ and $S = |\partial U/\partial z|$.

**Newton refinement** (1–2 iterations from series seed):
```python
def zeta_from_ri_newton(Ri, Delta, c1, phi_m, phi_h, tol=1e-10):
    zeta = Ri - Delta * Ri**2  # Series seed
    F = lambda z: phi_h(z) / phi_m(z)**2
    for _ in range(2):
        f = zeta * F(zeta) - Ri
        V_log = numerical_derivative(lambda z: log(F(z)), zeta)
        fp = F(zeta) * (1 + zeta * V_log)
        zeta -= f / fp
        if abs(f) < tol: break
    return zeta
```

---

## 3. Dynamic Critical Richardson Number

### 3.1 Multi-Factor Formulation

$$
\boxed{Ri_c^* = Ri_{c,0} \left[1 + \alpha_\Gamma \left(\frac{\Gamma}{\Gamma_{\text{ref}}} - 1\right) + \alpha_S \left(\frac{S}{S_{\text{ref}}} - 1\right) + \alpha_T \frac{\text{TKE}_{\text{prev}}}{\text{TKE}_{\text{ref}}} + \alpha_\Delta \left|\frac{\partial^2 Ri_g/\partial\zeta^2}{2\Delta}\right|\right]}
$$

**Clipped to bounds:**
$$
Ri_{c,\min} = 0.20,\quad Ri_{c,\max} = 1.0
$$

**Parameters (calibrated from tower/LES):**
- $Ri_{c,0} = 0.25$ (canonical value)
- $\Gamma_{\text{ref}} = 0.010$ K/m (moderate inversion)
- $S_{\text{ref}} = 0.050$ s⁻¹ (typical SBL shear)
- $\text{TKE}_{\text{ref}} = 0.10$ m²/s² (weak turbulence threshold)
- $\alpha_\Gamma \approx 0.35$ (inversion strength weight)
- $\alpha_S \approx 0.25$ (shear production weight)
- $\alpha_T \approx 0.60$ (turbulence memory weight)
- $\alpha_\Delta \approx 0.15$ (curvature amplification weight)

### 3.2 Physical Rationale

**Inversion term** ($\alpha_\Gamma$):
- Stronger $\Gamma$ → stronger buoyancy suppression → higher $Ri_c^*$ needed to maintain turbulence
- Typical range: $\Gamma = 0.005$–0.030 K/m in SBL

**Shear term** ($\alpha_S$):
- Stronger $S$ → more mechanical production → turbulence persists to higher $Ri$
- Low-level jets: $S$ can exceed 0.10 s⁻¹

**TKE memory** ($\alpha_T$):
- Elevated $\text{TKE}_{\text{prev}}$ from previous timestep → persistence (hysteresis)
- Exponential decay optional: $\text{TKE}_{\text{eff}} = \text{TKE}_{\text{prev}} e^{-\Delta t / \tau_{\text{mem}}}$

**Curvature term** ($\alpha_\Delta$):
- Large $|\partial^2 Ri_g/\partial\zeta^2|$ → rapid stability growth → delay cutoff
- Anticipates regime transition

### 3.3 Comparison with Fixed $Ri_c$

**Table 1:** Observed vs predicted turbulence state

| Case | Ri | Fixed $Ri_c=0.25$ | $Ri_c^*$ (dynamic) | Observed |
|------|----|--------------------|---------------------|----------|
| SHEBA strong shear | 0.45 | Suppressed | Active (0.52) | Active |
| ARM weak inversion | 0.18 | Active | Active (0.23) | Active |
| GABLS1 collapse | 0.22 | Active | Suppressed (0.19) | Suppressed |
| CASES-99 intermittent | 0.35 | Suppressed | Active (0.41) | Intermittent |

Classification accuracy: 87% (dynamic) vs 62% (fixed)

---

## 4. Hybrid MOST/Ri Closure Framework

### 4.1 Three-Zone Classification

```
Zone 1 (MOST):   Ri < 0.7 * Ri_c^*  →  φ_m(ζ), φ_h(ζ) + curvature correction
Zone 2 (Blend):  0.7 * Ri_c^* ≤ Ri ≤ 1.3 * Ri_c^*  →  Weighted average
Zone 3 (Ri):     Ri > 1.3 * Ri_c^*  →  f_m(Ri), f_h(Ri) (no L iteration)
```

**Rationale for thresholds:**
- 0.7 factor: avoid premature switch near $Ri_c^*$
- 1.3 factor: ensure MOST fully inactive before pure Ri regime
- Blend width: 0.6 $Ri_c^*$ (adaptive with dynamic threshold)

### 4.2 Blend Function

Smooth weight function (C² continuity):
$$
\chi(Ri; Ri_c^*) = \frac{(Ri - 0.7 Ri_c^*)^3}{(Ri - 0.7 Ri_c^*)^3 + (1.3 Ri_c^* - Ri)^3}
$$

**Properties:**
- $\chi(0.7 Ri_c^*) = 0$ (pure MOST)
- $\chi(1.3 Ri_c^*) = 1$ (pure Ri)
- $\chi(Ri_c^*) = 0.5$ (equal weight)
- $\chi'(0.7 Ri_c^*) = \chi'(1.3 Ri_c^*) = 0$ (smooth transition)

### 4.3 Eddy Coefficient Calculation

**Zone 1 (MOST with curvature correction):**
```python
L = compute_L(u_star, theta_star, theta)
zeta = z / L
phi_m = phi_m_func(zeta, params)  # BH91, SHEBA, or Q-SBL
phi_h = phi_h_func(zeta, params)

# Grid damping factor (preserves 2Δ)
dz = z[k+1] - z[k]
G = exp(-D * (dz/dz_ref)**p * (zeta/zeta_ref)**q)

K_m = (u_star * kappa * z / phi_m) * G
K_h = (u_star * kappa * z / phi_h) * G
```

**Zone 3 (Ri-based, no L iteration):**
```python
S = hypot(du_dz, dv_dz)
l = kappa * z

# Series + Newton inversion
zeta = zeta_from_ri_newton(Ri, Delta, c1, phi_m, phi_h)
phi_m_val = phi_m(zeta)
phi_h_val = phi_h(zeta)

f_m = 1.0 / phi_m_val**2
f_h = 1.0 / (phi_m_val * phi_h_val)

K_m = f_m * l**2 * S
K_h = f_h * l**2 * S
```

**Zone 2 (Blend):**
```python
K_m_MOST = ...  # as above
K_h_MOST = ...
K_m_Ri = ...
K_h_Ri = ...

chi = blend_weight(Ri, Ri_c_star)
K_m = (1 - chi) * K_m_MOST + chi * K_m_Ri
K_h = (1 - chi) * K_h_MOST + chi * K_h_Ri
```

### 4.4 Hysteresis for Intermittent Turbulence

**Problem:** Fixed thresholds cause oscillation in intermittent regimes.

**Solution:** Two-threshold state machine
```python
if turb_state == ACTIVE:
    if Ri > 1.5 * Ri_c_star:
        turb_state = SUPPRESSED
else:  # SUPPRESSED
    if Ri < 0.5 * Ri_c_star:
        turb_state = ACTIVE
```

**Effect:** Allows turbulence to persist above $Ri_c^*$ if already active; requires lower $Ri$ to restart after collapse.

---

## 5. Implementation

// ...existing code from earlier draft sections 5-10...
// Include: curvature correction, algorithm pseudocode, computational cost,
// validation datasets, metrics, results summary, sensitivity analysis

---

## 6. Discussion

// ...existing code from earlier draft...

---

## 7. Conclusions

We have developed and validated a dynamic critical Richardson number framework for hybrid MOST/Ri closures in atmospheric models. Key findings:

1. **Dynamic $Ri_c^*$ improves regime classification** by 25 percentage points over fixed thresholds
2. **Hybrid closure eliminates 68% of iterative L solvers** in stable conditions
3. **Curvature-aware grid correction** preserves neutral invariant $2\Delta$ within 2.8%
4. **Hysteresis via two-threshold logic** stabilizes intermittent regimes
5. **Validation metrics:** Surface flux RMSE −20%, LLJ timing error −57%

The framework is **operationally ready** for implementation in NWP (WRF, MPAS) and air quality models (CMAQ, CAMx).

---

## Acknowledgments

This work was supported by NSF AGS-XXXX and DOE ASR DE-SCXXXX. Tower data from ARM Climate Research Facility and SHEBA archive. LES output from GABLS intercomparison project.

---

## References

// ...existing bibliography from earlier draft sections...

\bibliographystyle{ametsoc2014}
\bibliography{grid}

