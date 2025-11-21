# Richardson-Number Curvature Corrections: Practical Guide for McNider et al.

## Executive Summary

Coarse vertical grids systematically underestimate local stability in the stable boundary layer (SBL) because the gradient Richardson number Ri_g(z) is **concave-down** (d²Ri_g/dz² < 0). This causes bulk layer estimates Ri_b to fall below point values Ri_g(z_g), producing excessive turbulent mixing. We present a **minimal, neutral-preserving multiplicative correction** fc that reduces this bias while maintaining correct near-neutral physics.

**Core principle:** Apply a damping factor fc < 1 to eddy diffusivities K_m, K_h (or stability function f_m) when the bias ratio B = Ri_g(z_g)/Ri_b exceeds a threshold (~1.05), with fc→1 as layer thickness Δz→0 and as ζ→0 (neutral conditions).

---

## 1. Physical Origin of the Bias

### 1.1 Curvature in MOST
Monin–Obukhov similarity gives:
$$
Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}, \quad \zeta = z/L.
$$

Near-neutral Taylor expansion (for smooth φ):
$$
Ri_g(\zeta) = \zeta + \Delta \zeta^2 + O(\zeta^3), \quad \Delta = a_h - 2a_m,
$$
where $a_{m,h}$ are the linear near-neutral slopes of $\phi_{m,h}$. The **neutral curvature** is:
$$
\left.\frac{d^2 Ri_g}{d\zeta^2}\right|_0 = 2\Delta.
$$

**Typical SBL:** Businger–Dyer (BD) coefficients give $a_m \approx 4.7$, $a_h \approx 7.8$ → $\Delta \approx -1.6$ → **concave-down** profile (d²Ri_g/dz² < 0).

### 1.2 Jensen's Inequality and Grid Bias
For a concave-down function over layer [z₀, z₁]:
$$
Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_{\text{rep}}),
$$
where $z_{\text{rep}}$ is any representative height in the layer. The **geometric mean** $z_g = \sqrt{z_0 z_1}$ is the natural choice because:
- It is the **midpoint in log-space**: $\ln z_g = \frac{1}{2}(\ln z_0 + \ln z_1)$.
- Profiles scaling with ln z (log-wind law, MOST) are best represented at z_g.
- Minimizes mean-square error for layer-averaged log-linear quantities.

**Bias ratio:**
$$
B = \frac{Ri_g(z_g)}{Ri_b}.
$$
- If B > 1: bulk estimate underestimates local stability → **over-mixing**.
- If B ≈ 1: negligible curvature bias.
- Typical coarse grids (Δz = 60–100 m): B = 1.2–2.0 in strong SBL.

### 1.3 Shear Denominator Bias
Bulk Richardson number:
$$
Ri_b = \frac{g}{\theta_{\text{ref}}} \frac{\Delta\theta \cdot \Delta z}{(\Delta U)^2}.
$$
If shear is estimated using **arithmetic mean** height instead of **log-mean** $z_L = (z_1 - z_0)/\ln(z_1/z_0)$, the denominator (ΔU)² is **overestimated** (shear underestimated), further **lowering** Ri_b and **amplifying** bias.

**Recommendation:** Use z_L for shear reconstruction when wind follows log-law; use z_g for evaluating point Ri_g(z_g).

---

## 2. Correction Strategy

### 2.1 Design Requirements
A physically sound correction must:
1. **Preserve neutral curvature** 2Δ: fc(ζ=0) = 1 and ∂fc/∂ζ|₀ = 0.
2. **Reduce coarse-grid bias**: fc < 1 when B > threshold.
3. **Vanish on fine grids**: fc → 1 as Δz → 0.
4. **Maintain stability**: no spurious oscillations or negative K.

### 2.2 Multiplicative Damping
Apply correction to diffusivities:
$$
K_m^* = K_m \cdot f_c(B, \Delta z, \zeta), \quad K_h^* = K_h \cdot f_c.
$$
Alternatively, apply to stability function:
$$
f_m^* = f_m \cdot f_c \quad \text{(before computing K)}.
$$

**Why multiplicative?**
- Preserves neutral limit (fc=1 at ζ=0).
- Physically interpretable: reduces turbulent transport proportionally.
- Avoids introducing new length/velocity scales.
- Easy to localize (apply only where bias detected).

---

## 3. Functional Forms for fc

### 3.1 Exponential (Recommended)
$$
f_c = \exp\left[ -\alpha \cdot (B-1) \cdot \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p \cdot \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q \right].
$$

**Parameters:**
- α: strength (start with α = 1.0; tune 0.5–1.5).
- p: Δz scaling exponent (typically p = 1).
- q: ζ scaling exponent (use q ≥ 2 to enforce ∂fc/∂ζ|₀ = 0).
- Δz_ref: reference layer thickness (10 m).
- ζ_ref: reference stability (0.5).

**Properties:**
- fc → 1 as Δz → 0 (fine-grid limit).
- fc → 1 as ζ → 0 (neutral preservation).
- fc → exp(−αΔz/Δz_ref) · exp(−constant) at large ζ (smooth asymptote).
- Always positive; smooth derivatives.

### 3.2 Rational (Alternative)
$$
f_c = \frac{1}{1 + \alpha \cdot \max(0, B-1) \cdot \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p \cdot \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q}.
$$

**Properties:**
- Bounded: fc ∈ (0, 1].
- Milder tail behavior than exponential.
- Easier to interpret coefficient α (fc ≈ 1/(1+α·...) for small corrections).

**Default starting values:** α = 0.8, p = 1, q = 2.

### 3.3 Safety Limits
Apply hard floor to prevent over-damping:
$$
fc_{\text{final}} = \max(f_c, fc_{\text{min}}), \quad fc_{\text{min}} \ge 0.2.
$$

---

## 4. Computational Procedure

### Step 1: Compute Representative Heights
For layer [z₀, z₁]:
$$
z_g = \sqrt{z_0 z_1} \quad \text{(geometric mean)},
$$
$$
z_L = \frac{z_1 - z_0}{\ln(z_1 / z_0)} \quad \text{(log-mean, for shear)}.
$$

### Step 2: Evaluate Point Ri_g(z_g)
- **If φ-functions available:** Compute ζ_g = z_g/L, evaluate:
  $$
  Ri_g(z_g) = \zeta_g \frac{\phi_h(\zeta_g)}{\phi_m(\zeta_g)^2}.
  $$
- **If φ unknown:** Use centered finite differences at z_g:
  $$
  Ri_g(z_g) = \frac{g}{\theta(z_g)} \frac{(\partial\theta/\partial z)|_{z_g}}{(\partial U/\partial z)^2|_{z_g}}.
  $$

### Step 3: Compute Bulk Ri_b
**Method A (Direct):**
$$
Ri_b = \frac{g}{\theta_{\text{ref}}} \frac{(\theta_1 - \theta_0)(z_1 - z_0)}{(U_1 - U_0)^2}.
$$
- Use θ_ref = 0.5(θ₀ + θ₁).
- For log-wind consistency, reconstruct ΔU using z_L (if applicable).

**Method B (Simpson integration, if fine Ri_g profile available):**
$$
Ri_b \approx \frac{1}{6}\left[Ri_g(z_0) + 4 Ri_g(z_g) + Ri_g(z_1)\right].
$$

### Step 4: Compute Bias Ratio
$$
B = \frac{Ri_g(z_g)}{Ri_b}.
$$
Guard against division by zero: if Ri_b ≤ ε (small threshold, e.g., 10⁻⁶), set B = 1 (neutral default).

### Step 5: Decision and Correction
```
if B ≤ B_threshold:  # typically 1.05
    f_c = 1.0  # no correction
else:
    # Compute ζ = z_g / L
    ζ = z_g / L
    # Choose template (exponential recommended)
    fc = exp(-α * (B - 1) * (Δz / Δz_ref)^p * (ζ / ζ_ref)^q)
    fc = max(fc, fc_min)  # apply floor
end
```

### Step 6: Apply Correction
$$
K_m^* = K_m \cdot f_c, \quad K_h^* = K_h \cdot f_c.
$$
Use K*_m, K*_h in vertical diffusion update.

---

## 5. Where to Apply (Operational Notes)

### 5.1 K-Based Models (Interior Flux from K)
- Apply fc to K_m, K_h at the **first interior level** (where K is actually used to compute flux).
- Do NOT alter surface-prescribed fluxes directly (those come from surface scheme).

### 5.2 Surface-Flux Models (Prescribed τ_sfc, H_sfc)
- Surface stress and heat flux are **control quantities** (not computed from K at first level).
- Apply fc to **interior K only** (second level and above) to limit vertical intrusion of surface fluxes into elevated stable layers.

### 5.3 Mixing-Length Alternative
Instead of modifying K directly, reduce effective mixing length:
$$
l^* = \frac{l}{1 + \beta (B-1) (\Delta z / \Delta z_{\text{ref}})^p}.
$$
Then K* = (l*)² S f_m(Ri). Choose β ≈ 0.5–1.0; calibrate for target bias reduction.

---

## 6. Tuning and Validation

### 6.1 Neutral-Preservation Tests
**Unit test 1:** fc(ζ=0, any Δz) = 1 ± ε (ε < 10⁻⁶).

**Unit test 2:** ∂fc/∂ζ|₀ ≈ 0 (verify numerically with small ζ perturbation).

### 6.2 Coarse-Grid Validation
**Reference case:** Fine-grid LES or tower data (Δz ≈ 5–10 m) → B_ref ≈ 1.0–1.05.

**Coarse baseline:** Same case with Δz = 60–100 m (no correction) → B_before (measure).

**Corrected:** Apply fc → B_after (should reduce toward 1 by ~30–50%).

**Target metric:**
$$
\text{Bias reduction} = \frac{B_{\text{before}} - B_{\text{after}}}{B_{\text{before}} - 1} \times 100\%.
$$
Aim for 40–60% reduction without degrading neutral or convective performance.

### 6.3 Sensitivity and Diagnostics
**Log per timestep (or hourly):**
- B (bias ratio).
- fc (correction factor).
- K_old / K_new (fractional change).
- Surface fluxes (τ, H).
- Inversion height (if diagnosed).

**Monitor:**
- No spurious oscillations in K or fluxes.
- Surface energy balance closure (if coupled).
- Comparison with uncorrected run and observations.

### 6.4 Parameter Space Exploration
Start with defaults (α=1.0, p=1, q=2, Δz_ref=10 m, ζ_ref=0.5, fc_min=0.2). Vary:
- α: 0.5, 0.8, 1.0, 1.2, 1.5 (strength).
- q: 1.5, 2.0, 2.5 (neutral-preservation sharpness).
- fc_min: 0.15, 0.2, 0.25 (over-damping guard).

Select combination minimizing bias while maintaining stable integration and neutral fidelity.

---

## 7. Summary: Quick Reference

| **Quantity** | **Purpose** | **Formula / Notes** |
|--------------|-------------|---------------------|
| z_g | Representative height (point Ri) | √(z₀ z₁) |
| z_L | Log-mean (shear reconstruction) | (z₁−z₀)/ln(z₁/z₀) |
| Ri_g(z_g) | Point gradient Richardson | ζ_g φ_h / φ_m² or centered differences |
| Ri_b | Bulk Richardson (layer average) | g/θ_ref · Δθ Δz / (ΔU)² or Simpson |
| B | Bias ratio | Ri_g(z_g) / Ri_b |
| fc | Correction factor | exp[−α(B−1)(Δz/Δz_ref)^p(ζ/ζ_ref)^q] |
| K* | Corrected diffusivity | K · fc |

**Decision threshold:** Apply correction if B > 1.05.

**Neutral preservation:** fc(ζ=0) = 1, ∂fc/∂ζ|₀ = 0 (use q ≥ 2).

**Convergence:** fc → 1 as Δz → 0.

---

## 8. Physical Justification (Expanded)

### 8.1 Why Concave-Down Implies Overmixing
In the SBL, turbulence is suppressed by buoyancy. A concave-down Ri_g(z) means:
- **Near surface (small z):** Ri_g grows slowly (turbulence still active).
- **Aloft (larger z):** Ri_g accelerates upward (rapid stratification, turbulence collapse).

On a coarse grid, the single bulk value Ri_b **averages over this rapid transition**, yielding a value below the local Ri_g at the layer's dynamically important height (z_g). The model then uses this **too-low Ri_b** to compute K, resulting in **too-large K** and **excessive mixing**.

### 8.2 Why fc Must Preserve 2Δ
The neutral curvature 2Δ = 2(a_h − 2a_m) is a **fundamental MOST invariant** tied to:
- The **turbulent Prandtl number** near-neutral slope: dPr_t/dζ|₀ = a_h − a_m.
- The **initial nonlinear response** of stability to stratification.

Any correction that alters 2Δ **breaks the correct neutral-limit physics**, potentially:
- Changing the neutral shear/heat flux ratio.
- Introducing spurious vertical momentum/heat transport at ζ=0.
- Violating established empirical MOST coefficients (Businger–Dyer, etc.).

By requiring fc(ζ=0)=1 and ∂fc/∂ζ|₀=0, we ensure the correction is **inactive near neutral** and only engages where curvature-driven bias is significant (stable conditions, coarse grids).

### 8.3 Why Multiplicative vs. Additive
**Additive correction** (e.g., K* = K − ΔK):
- Risks negative K (unstable).
- Introduces new dimensional scales (ΔK units).
- Hard to make neutral-preserving without case-by-case tuning.

**Multiplicative correction** (K* = K · fc):
- Always positive if fc ∈ (0,1].
- Dimensionless (fc is a pure number).
- Scales naturally with base K (large K → large reduction; small K → small change).
- Easy neutral preservation: fc=1 at ζ=0.

### 8.4 Why Geometric Mean z_g
For quantities scaling with ln z (log-wind, MOST):
$$
\ln z_g = \frac{1}{2}(\ln z_0 + \ln z_1) = \frac{1}{\Delta z}\int_{z_0}^{z_1} \ln z \, dz.
$$
Thus z_g is the **log-space average**, minimizing representation error for log-linear profiles. The arithmetic mean z_a = (z₀+z₁)/2 biases **high** (towards top of layer) for large z₁/z₀, systematically overestimating representative height and underestimating local stability gradients.

**Practical impact:** Using z_a instead of z_g can inflate B by 10–20% in coarse grids, further amplifying overmixing.

---

## 9. Connection to McNider ODE (Optional Advanced Note)

McNider proposed a differential constraint for f_c:
$$
\frac{d}{d(\Delta z)}(f_s(Ri) \cdot f_s(Ri,\Delta z)) = 0,
$$
where f_s is a reference stability function. This enforces that the **product** (fc · layer-thickness · stability) remains invariant under grid refinement—a mathematical expression of the requirement that the **integrated turbulent transport** should be grid-independent.

Taking logarithmic derivative:
$$
\frac{d \ln f_c}{d \ln \Delta z} = -\alpha (B-1) \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q,
$$
which integrates to:
$$
f_c(\Delta z, \zeta) = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}.
$$

For small exponents, this approximates:
$$
f_c \approx \exp\left[-\alpha (B-1) \ln\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right) \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q\right],
$$
which is equivalent to the exponential form with p=1 and logarithmic Δz scaling. The power-law form:
$$
\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p
$$
is a **practical generalization** allowing non-logarithmic Δz dependence (p≠1) for empirical tuning.

**Takeaway:** The exponential fc is not ad hoc—it is a **closed-form solution** to a physically motivated ODE expressing grid-invariance of integrated turbulent flux.

---

## 10. References and Contacts

**Derivation and diagnostics:**
- Ri.md (neutral curvature, Jensen inequality, z_g justification).
- Curvature_Demo.ipynb (interactive sliders, synthetic NetCDF generator).

**Implementation examples:**
- notebooks/Curvature_Demo.ipynb (Python, fully worked).

**Contact:** David E. England
Provide case-specific data (z₀, z₁, L, observed Ri_b, K) for tuned fc suggestions.

---

## Appendix: Pseudocode (Complete Example)

```python
# Inputs: z0, z1 (layer bounds), U0, U1, theta0, theta1, L (Obukhov length)
# Optional: phi_m, phi_h callables if MOST available

import math

# Parameters (defaults)
alpha = 1.0
p = 1.0
q = 2.0
Dz_ref = 10.0  # m
zeta_ref = 0.5
fc_min = 0.2
B_threshold = 1.05
g = 9.81       # m/s²

# Step 1: Representative heights
z_g = math.sqrt(z0 * z1)
z_L = (z1 - z0) / math.log(z1 / z0)
Dz = z1 - z0

# Step 2: Point Ri_g(z_g)
# Option A: MOST (if phi available)
zeta_g = z_g / L
phi_m_g = phi_m(zeta_g)  # user-supplied function
phi_h_g = phi_h(zeta_g)
Ri_g_zg = zeta_g * phi_h_g / (phi_m_g**2)

# Option B: Finite difference (if phi unavailable)
# dU_dz_g = (U1 - U0) / (z1 - z0)  # crude; better: centered at z_g
# dtheta_dz_g = (theta1 - theta0) / (z1 - z0)
# theta_g = 0.5 * (theta0 + theta1)
# Ri_g_zg = (g / theta_g) * (dtheta_dz_g / (dU_dz_g**2))

# Step 3: Bulk Ri_b
theta_ref = 0.5 * (theta0 + theta1)
Delta_theta = theta1 - theta0
Delta_U = U1 - U0
Ri_b = (g / theta_ref) * (Delta_theta * Dz) / (Delta_U**2) if Delta_U != 0 else 1e-10

# Step 4: Bias ratio
B = Ri_g_zg / Ri_b if Ri_b > 1e-10 else 1.0

# Step 5: Correction factor
if B <= B_threshold:
    fc = 1.0
else:
    zeta = z_g / L
    exponent = -alpha * (B - 1) * (Dz / Dz_ref)**p * (zeta / zeta_ref)**q
    fc = math.exp(exponent)
    fc = max(fc, fc_min)

# Step 6: Apply (assume K_old computed from uncorrected closure)
K_m_new = K_m_old * fc
K_h_new = K_h_old * fc

# Diagnostics
print(f"z_g={z_g:.2f}, Ri_g(z_g)={Ri_g_zg:.4f}, Ri_b={Ri_b:.4f}, B={B:.3f}, fc={fc:.3f}")
```

---

**End of Document**
