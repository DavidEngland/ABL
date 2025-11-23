Let’s carefully find the first derivative with respect to \(z\) of
\[
f(z) = k\,z\left(1 + 0.47\,\frac{z}{L}\right)^{-1}.
\]

---

### Step 1. Write it more cleanly
\[
f(z) = k\,z \left(1 + \frac{0.47 z}{L}\right)^{-1}.
\]

---

### Step 2. Use the product rule
If \(f(z) = k\,[z \cdot g(z)]\) with \(g(z) = \left(1 + 0.47 z / L\right)^{-1}\), then
\[
\frac{df}{dz} = k \left[g(z) + z\,g'(z)\right].
\]

---

### Step 3. Find \(g'(z)\)
\[
g(z) = \left(1 + \frac{0.47 z}{L}\right)^{-1}.
\]
Using the chain rule,
\[
g'(z) = -\left(1 + \frac{0.47 z}{L}\right)^{-2}\cdot \frac{0.47}{L}
= -\frac{0.47/L}{\left(1 + 0.47 z / L\right)^{2}}.
\]

---

### Step 4. Substitute back
\[
\frac{df}{dz}
= k \left[
\left(1 + \frac{0.47 z}{L}\right)^{-1}
- z\,\frac{0.47/L}{\left(1 + 0.47 z / L\right)^{2}}
\right].
\]

---

### Step 5. Simplify
Put everything over the common denominator \(\left(1 + 0.47 z / L\right)^2\):
\[
\frac{df}{dz}
= k \,\frac{
\left(1 + \frac{0.47 z}{L}\right) - \frac{0.47 z}{L}
}{
\left(1 + \frac{0.47 z}{L}\right)^{2}
}
= k \,\frac{1}{\left(1 + 0.47 z / L\right)^{2}}.
\]

---

✅ Final derivative:
\[
\boxed{
\frac{df}{dz}
= \frac{k}{\left(1 + 0.47\,\dfrac{z}{L}\right)^{2}}
}.
\]

# Pros and Cons — McNider Similarity-Based Ri Curvature (Quick Brief)

Context

- Idea: use MOST stability functions near the surface (φ_m, φ_h; e.g., linear or exponential approximations) to compute Ri_g(ζ) = ζ φ_h/φ_m² and its curvature d²Ri_g/dζ² analytically; map to height via 1/L² (constant L) or full chain rule for L(z).

Pros

- Physics-grounded: directly ties curvature to φ_m, φ_h; preserves neutral invariant 2Δ if φ near-neutral slopes are respected.
- Closed-form and cheap: no ζ-iteration; easy to evaluate at representative height (prefer z_g).
- Drop-in: compatible with existing MOST codes; can be used to set curvature-aware lower-boundary conditions.
- Transparent diagnostics: 2Δ, c1, V_log, W_log available for QC and tuning.

Cons

- Sensitivity to φ choice: linear vs exponential vs BH/SHEBA fits change Δ and c1; misfit propagates to curvature.
- Domain/guards: power-law poles (ζ<1/β) or exponential tails need guards; variable L(z) can bias height mapping if ignored (use E_omit).
- Near-surface grid bias remains if evaluated at arithmetic means; must use geometric/log means for layers.
- Limited aloft: similarity-only curvature may miss secondary inflections or nonlocal effects above the surface layer.
- Pr_t dependence: assuming Pr_t≈1 simplifies algebra but can bias heat curvature in strong SBL unless calibrated.

Recommendations

- Use similarity-based curvature at z_g for the first layer; preserve 2Δ; apply constant-L map only if E_omit<0.05.
- Prefer Q‑SBL (quadratic) surrogate for ζ≤0.2–0.3 in very stable cases to avoid pole artifacts.
- If residual bias at coarse Δz, add neutral-preserving grid damping G(ζ,Δz); keep G(0)=1 and ∂G/∂ζ|0=0.

Decision asks

- Confirm φ baseline (linear/Q‑SBL/BH) for experiments.
- Approve z_g usage and E_omit threshold (0.05) for height mapping choice.

# Pros and Cons of Similarity-Based Near-Surface Ri Curvature (McNider Approach)

Scope

- Compute Ri_g(ζ)=ζ φ_h/φ_m² and d²Ri_g/dζ² analytically from chosen MOST φ_m, φ_h (linear, BH91-like, exponential surrogate). Map to height via 1/L² (constant L) or full chain rule for L(z). Use at a representative height (z_g).

Pros

- Neutral fidelity: preserves neutral curvature 2Δ exactly when φ near-neutral slopes are matched.
- Analytical clarity: closed-form V_log, W_log enable fast, stable evaluation and clear diagnostics (Δ, c1, ζ_inf).
- Computationally light: no ζ-iteration; suitable for operational first-layer BCs.
- Integrates with grid-aware fixes: pairs naturally with geometric/log means and optional neutral-preserving damping G(ζ,Δz).

Cons

- Parameter dependence: curvature quality hinges on φ calibration in the SBL (Pr_t, αβ or linear slopes); unstable-derived fits can overstate |Δ|.
- Domain issues: power-law poles and over-rapid exponential tails require guards/blends (ζ≤0.2–0.3 recommended for direct use).
- Mapping error if L varies: constant-L shortcut biases height curvature when E_omit is not small; must switch to full mapping.
- Nonlocal/aloft limitations: similarity-only curvature misses nonlocal transport or elevated inflections; needs blending or separate aloft treatment.

Practical guidance

- First layer: evaluate at z_g=√(z₀z₁); report 2Δ and use constant-L map only if E_omit<0.05; else apply full chain rule.
- Very stable: prefer Q‑SBL (quadratic) φ up to ζ≈0.2–0.3; blend to capped form aloft.
- Coarse grids: if bias persists, apply G(ζ,Δz)=exp[−D(Δz/Δz_r)^p(ζ/ζ_r)^q] with G(0)=1 and ∂G/∂ζ|₀=0 to preserve neutrality.

Notes for Arctic use

- Expect stronger |Δ| and higher ζ; prioritize SBL-calibrated φ (SHEBA/Arctic fits) and variable-L mapping checks (E_omit diagnostic).

## Special case: Pr = 1 and φ(ζ) = 1 + β ζ — closed forms and operational recipes

Assume φ_m(ζ)=φ_h(ζ)=1+βζ and Pr=1. Then
- Stability ratio:
  \[
  F(\zeta)=\frac{\phi_h}{\phi_m^2}=\frac{1}{1+\beta\zeta},
  \]
  and the gradient Richardson number is
  \[
  Ri_g(\zeta)=\frac{\zeta}{1+\beta\zeta}.
  \]

- Near‑neutral series (ζ small):
  \[
  Ri_g(\zeta)=\zeta-\beta\zeta^2+\beta^2\zeta^3-\beta^3\zeta^4+\dots
  \]

- Derivatives (closed form):
  \[
  \frac{dRi_g}{d\zeta}=\frac{1}{(1+\beta\zeta)^2},\qquad
  \frac{d^2Ri_g}{d\zeta^2}=-\frac{2\beta}{(1+\beta\zeta)^3}.
  \]
  Note: for β>0 the curvature is negative (concave‑down) as typical in the SBL.

- Bulk (layer average) over ζ∈[0,ζ] (use ζ as nondimensional layer thickness; map back via ζ = z/L if needed):
  \[
  Ri_b(\zeta)=\frac{1}{\zeta}\int_0^\zeta\frac{\xi}{1+\beta\xi}\,d\xi
  =\frac{1}{\zeta\beta^2}\Big(\beta\zeta-\ln(1+\beta\zeta)\Big).
  \]

- Bias ratio (point vs bulk) for the layer:
  \[
  B(\zeta)=\frac{Ri_g(\zeta)}{Ri_b(\zeta)}
  =\frac{\zeta^2\beta^2}{(1+\beta\zeta)\big(\beta\zeta-\ln(1+\beta\zeta)\big)}.
  \]
  Evaluate B at the geometric‑mean height by using ζ_g = z_g/L when converting heights.

- Neutral‑preserving grid correction template (exponential form; q≥2 enforces zero slope at ζ→0):
  \[
  f_c(\Delta z,\zeta)=\exp\!\Big[-D\,(B(\zeta)-1)\Big(\frac{\Delta z}{\Delta z_{\rm ref}}\Big)^p\Big(\frac{\zeta}{\zeta_{\rm ref}}\Big)^q\Big],
  \]
  with suggested defaults Δz_ref=10 m, ζ_ref=0.5, p=1, q=2, D tuned to target bias reduction (e.g., D≈0.3–1.0).

- Dynamic critical Richardson number (simple operational form):
  \[
  Ri_c^*(\zeta)=Ri_{c,0}\Big[1+\gamma\big(B(\zeta)-1\big)\Big],
  \]
  with Ri_{c,0}≈0.25 and γ∈[0,1] (γ>0 raises Ri_c where curvature bias is large). Use smoothed B (time or vertical average) to avoid noisy fluctuations.

- Quick pseudocode (φ‑agnostic implementation using the special-case formulas):
```python
# inputs: z0,z1,L,beta,Delta_z,Delta_z_ref=10.0,zeta_ref=0.5
z_g = sqrt(z0*z1)
zeta_g = z_g / L
beta_z = beta * zeta_g
Ri_g = zeta_g / (1.0 + beta_z)
Ri_b = (1.0 / (zeta_layer * beta**2)) * (beta*zeta_layer - math.log(1.0 + beta*zeta_layer))
B = (Ri_g) / (Ri_b if Ri_b>1e-12 else Ri_g)
# correction factor (defaults)
D = 0.6; p = 1.0; q = 2.0
fc = math.exp(-D*(B-1.0)*(Delta_z/Delta_z_ref)**p*(zeta_g/zeta_ref)**q)
Ri_c_star = Ri_c0*(1.0 + gamma*(B-1.0))
# apply fc to K (or to mixing length)
K_new = K_old * fc
```

Operational notes:
- Use ζ mapping (ζ=z/L) with local or bulk L consistently. If L varies strongly, compute ζ_g using local L(z_g) and guard with the omission metric E_omit (see main text).
- Smooth B in time (running mean) and vertically (small stencil) before applying fc or Ri_c^* to avoid noisy activation.
- The special-case rational form gives closed analytic expressions for Ri_g, Ri_b and B useful for rapid testing and tuning; when φ deviates from identical linear form, replace Ri_g/Ri_b with numerically evaluated point/integral forms but retain the same correction template.
