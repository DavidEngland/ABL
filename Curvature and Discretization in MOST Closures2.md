# Technical Note: Curvature and Vertical Discretization in MOST-Based Closures

## 0. Abstract
This note refines the interpretation of the gradient Richardson number curvature
$$
Ri_g(\zeta)=\zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2},\qquad \zeta=\frac{z}{L},
$$
and its implications for stability-function formulation and vertical discretization in Monin–Obukhov Similarity Theory (MOST). We formalize near-neutral curvature control via the coefficient \(\Delta\), provide quantitative criteria for when a constant-\(L\) approximation is acceptable, and justify the geometric-mean height for layer evaluation. Implementation diagnostics and a ready figure caption are included.

## 1. Core Curvature Expression
Let
$$
F(\zeta)=\frac{\phi_h}{\phi_m^2},\quad
V_{\log}=\frac{\phi_h'}{\phi_h}-2\frac{\phi_m'}{\phi_m},\quad
W_{\log}=\frac{dV_{\log}}{d\zeta}.
$$
Then
$$
\frac{d^2Ri_g}{d\zeta^2}=F(\zeta)\Big[2V_{\log}+\zeta\big(V_{\log}^2-W_{\log}\big)\Big].
$$

Near neutrality (\(\zeta\to0\)):
$$
\Delta=V_{\log}(0)=\alpha_h\beta_h-2\alpha_m\beta_m,\quad
\left.\frac{d^2Ri_g}{d\zeta^2}\right|_{0}=2\Delta.
$$

Interpretation:
- \(\Delta>0\): initial concave-up (heat corrections dominate).
- \(\Delta<0\): initial concave-down (momentum corrections dominate).
- \(|\Delta|\) sets the strength of first departure from linear \(Ri_g(\zeta)\).

## 1A. Unified Log Terms
Component logs:
$$
v_m=\frac{\phi_m'}{\phi_m},\ v_h=\frac{\phi_h'}{\phi_h},\ V_{\log}=v_h-2v_m,\ W_{\log}=V_{\log}'.
$$
Curvature:
$$
\partial_\zeta^2 Ri_g = F[2V_{\log}+\zeta(V_{\log}^2-W_{\log})],\quad F=\phi_h/\phi_m^2.
$$
This replaces earlier generic G,G′ (G=V_log, G′=W_log).

## 2. Elevated Regime Behavior
Away from neutral (\(\zeta\gesssim 0.1\) or approaching stable growth), the \(\zeta(V_{\log}^2-W_{\log})\) term drives rapid curvature escalation or decay. Parameter sets that keep \(V_{\log}^2\) and \(W_{\log}\) in partial balance delay curvature amplification, improving numerical stability in coarse vertical grids.

## 3. Discretization and Representative Height Selection
Given two model levels $z_1<z_2$:
- Arithmetic mean $z_a=(z_1+z_2)/2$ biases MOST because profiles are log-like.
- Geometric mean $z_g=\sqrt{z_1 z_2}$ satisfies
$$
\ln z_g \;=\; \frac{1}{z_2-z_1}\int_{z_1}^{z_2} \ln z\,dz,
$$
minimizing midpoint error for log-linear structure.

Approximate curvature preservation error when using $z_a$ instead of $z_g$:
$$
\epsilon_{\text{geom}}\approx \frac{1}{2}\Big|\frac{d^2Ri_g}{dz^2}\Big|_{z_g}\left(z_a-z_g\right)^2
$$
while for $z_g$ the first-order logarithmic bias term cancels.

Recommendation: use $z_g$ for evaluating \(\phi_{m,h}\), \(Ri_g\), and curvature when \(\Delta z / z_1>0.2\).

## 3.5 Numerical Integration of Ri_g for Ri_b

If full profile $Ri_g(z)$ is available (analytic or interpolated):
$$
Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz.
$$

**Trapezoid rule (2-point):**
$$
Ri_b \approx \frac{1}{2}[Ri_g(z_0) + Ri_g(z_1)].
$$

**Simpson's rule (3-point):**
$$
Ri_b \approx \frac{1}{6}[Ri_g(z_0) + 4Ri_g(z_g) + Ri_g(z_1)].
$$

**Adaptive quadrature:**  
Use SciPy `quad` or Gauss–Legendre nodes for smooth $Ri_g(z)$.

**Comparison:**  
Validate bulk formula $Ri_b = g\Delta\theta\Delta z / (\Delta U)^2$ against numerical integral to assess curvature-induced bias.

## 4. Height Mapping with Variable $L(z)$
General chain rule:
$$
\frac{\partial^2Ri_g}{\partial z^2} \;=\; \left(\frac{d\zeta}{dz}\right)^2\frac{d^2Ri_g}{d\zeta^2} + \frac{d^2\zeta}{dz^2}\frac{dRi_g}{d\zeta},\qquad \zeta=\frac{z}{L(z)}.
$$
Variation metrics:
$$
\varepsilon_1=\frac{z|L'|}{L},\qquad
\chi=\left|\frac{(d^2\zeta/dz^2)(dRi_g/d\zeta)}{(d\zeta/dz)^2(d^2Ri_g/d\zeta^2)}\right|.
$$
Use constant-\(L\) shortcut (\(\partial_z^2Ri_g\approx (1/L^2)\partial_\zeta^2Ri_g\)) if \(\varepsilon_1<0.05\) and \(\chi<0.05\); otherwise include full mapping.

## 5. Discretization Error Diagnostic
Layer reconstruction for a vertical interval \([z_i,z_{i+1}]\):
$$
E_i = \Big|\frac{Ri_g(z_{i+1}) - Ri_g(z_i)}{z_{i+1}-z_i} - \left.\frac{dRi_g}{dz}\right|_{z_g}\Big|
$$
Estimate using \(\frac{dRi_g}{dz} \approx (d\zeta/dz)dRi_g/d\zeta\). Flag layers with \(E_i / |dRi_g/dz|_{z_g} > \eta\) (e.g. \(\eta=0.2\)) for adaptive refinement or alternative φ-form (e.g. quadratic surrogate).

## 6. Practical Implementation Steps
1. Compute φ_m, φ_h at each level; evaluate \(F, V_{\log}, W_{\log}\).
2. Determine \(\Delta\) and neutral curvature; store for run metadata.
3. For each layer, compute \(z_g\); evaluate curvature at \(z_g\).
4. Assess \(L(z)\) variability; choose mapping logic (Section 4).
5. Compute \(E_i\) (Section 5); apply smoothing only if curvature spikes are noise (avoid suppressing true inflection).
6. Diagnostics output: \(\Delta, \partial_\zeta^2Ri_g|_{0},\) max \(|\partial_z^2Ri_g|\), first \(\zeta_{\text{inf}}\) if real and inside domain, fraction of layers with \(\chi>0.05\).

## 7. Blended / Piecewise Stability Functions
To enforce curvature continuity from near-neutral to elevated stable layers:
- Match neutral curvature: ensure surrogate (quadratic or regularized) has same \(2\Delta\).
- Continuity target: \(\left|\Delta_{\text{curv}}\right| = \left|(\partial_\zeta^2Ri_g^{\text{base}} - \partial_\zeta^2Ri_g^{\text{surrogate}})/\partial_\zeta^2Ri_g^{\text{base}}\right| < 0.05\) at blend height \(\zeta_b\).
- Use sigmoid blend in ζ (e.g. \(\sigma(\zeta)=1/(1+e^{-s(\zeta-\zeta_b)})\)) if a sharp switch creates curvature artifacts.

Add: Exponential Ri-based K-closure (RI_EXP)
- Form: f(Ri) = exp(−γ Ri / Ri_c) for use in K = u_* l f(Ri).
- Small-Ri expansion: f ≈ 1 − (γ/Ri_c) Ri + 0.5 (γ/Ri_c)^2 Ri^2.
  - If you target near-neutral series coefficients (f ≈ 1 + s1 Ri + s2 Ri^2), match:
    - s1 = −γ / Ri_c
    - s2 = 0.5 (γ / Ri_c)^2
  - Choose (γ, Ri_c) to fit observed/desired decay of K with Ri in the weakly stable regime.
- Practical defaults: γ ≈ 3.2; Ri_c ≈ 0.25–0.5; allow γ_h ≠ γ_m and Ri_c_h ≠ Ri_c_m.

Blending suggestion
- Use RI_EXP when Ri is prognosed/diagnosed reliably; otherwise fallback to ζ-based φ with curvature-aware correction.
- Ensure smooth switch by a sigmoid in Ri around Ri ≈ 0.1–0.2 to avoid kinks in K.

## 7A. Mapping φ(ζ) to f(Ri) for Mixing-Length Closures

When using the mixing-length form
$$
K_{m,h} = f_{m,h}(Ri)\; l^2\; S,
$$
with shear magnitude
$$
S = \sqrt{\left(\frac{\partial u}{\partial z}\right)^2 + \left(\frac{\partial v}{\partial z}\right)^2},\quad l \approx \kappa z,
$$
the MOST-to-Ri mapping is
$$
f_m(Ri) = \frac{1}{\phi_m(\zeta(Ri))^2},\qquad
f_h(Ri) = \frac{1}{\phi_m(\zeta(Ri))\,\phi_h(\zeta(Ri))}.
$$
Derivation uses MOST identities:
\(\partial U/\partial z = (u_*/\kappa z)\phi_m\),
\(K_m=u_* \kappa z/\phi_m\), \(K_h=u_* \kappa z/\phi_h\),
and S = |\partial U/\partial z|.

Generalized lengths l_m, l_h yield
$$
f_m = \frac{u_* \kappa z}{\phi_m\,l_m^2\,S},\qquad
f_h = \frac{u_* \kappa z}{\phi_h\,l_h^2\,S}.
$$

## 8. Recommended Diagnostics Summary
| Metric | Purpose | Threshold |
|--------|---------|-----------|
| \(2\Delta\) | Neutral curvature sign/magnitude | Report |
| \(\zeta_{\text{inf}}\) | Inflection location | Only if real & < domain limit |
| \(\max\lvert \partial_z^2 Ri_g \rvert\) | Vertical stability stress | Watch for spikes |
| \(\chi\) fraction | Variable-L impact | &lt;10% desirable |
| \(E_i\) distribution | Layer reconstruction error | &gt;0.2 flagged |
| Blend mismatch | Curvature continuity | &lt;5% |

## 9. Figure Caption (Curvature vs. \(\zeta\))
Caption candidate:
“Gradient Richardson number curvature \(\partial_\zeta^2Ri_g\) as a function of nondimensional height \(\zeta=z/L\). Neutral curvature \(2\Delta\) sets the initial concavity; rapid growth governed by \(\zeta(V_{\log}^2-W_{\log})\). Shaded band shows layers where \(\chi>0.05\), indicating significant \(L(z)\) variation. Dashed curve: quadratic stable surrogate matched on neutral curvature.”

## 10. Common Pitfalls
- Using arithmetic mean height with large \(\Delta z\): biases curvature and suppresses inflection detection.
- Over-smoothing derivatives: can erase legitimate curvature sign changes.
- Ignoring \(L(z)\) variability when \(\chi\) is large: misrepresents physical curvature scaling.
- Applying power-law φ beyond pole proximity without guard: can inflate curvature and collapse flux iterates.

## 11. Inflection Flag
Compute ζ_inf from curvature root; if 0<ζ_inf<ζ_top:
- Split layer; report B1 (concave-down) and B2 (concave-up).
- Apply damping only where curvature negative.
Avoid single B metric in mixed-concavity layers; document both contributions for diagnostics.

## 12. Summary
Curvature structure is controlled near neutrality by a single coefficient \(\Delta\); representative-height choice (geometric mean) reduces discretization bias for log-structured profiles. Quantitative metrics (\(\varepsilon_1,\chi,E_i\)) guide when constant-\(L\) assumptions and layer compression are valid. Implementing curvature-aware blending stabilizes SBL behavior and improves diagnostic fidelity.

## 13. (Optional) Quick Code Stub
```python
def layer_curvature(z1,z2,L1,L2,phi_m,phi_h):
    import math
    z_g = math.sqrt(z1*z2)
    # user supplies phi_m/h callables at ζ=z/L(z); approximate L(z_g) by geometric mean
    L_g = math.sqrt(L1*L2)
    zeta_g = z_g/L_g
    # central differences for derivatives (small h)
    h = 1e-6
    def F(zeta): return phi_h(zeta)/phi_m(zeta)**2
    def logF(zeta): return math.log(F(zeta))
    V = (logF(zeta_g+h)-logF(zeta_g-h))/(2*h)
    W = (logF(zeta_g+h)-2*logF(zeta_g)+logF(zeta_g-h))/(h*h)
    curv_zeta = F(zeta_g)*(2*V + zeta_g*(V*V - W))
    return curv_zeta
```

## 14. Critique of Original Draft (For Record)
- Mixed plain-text and LaTeX markers (e.g., `\partial^2 Ri_g/\partial \zeta^2`) without consistent math environments.
- Missing explicit definitions of \(V_{\log}, W_{\log}, \Delta\).
- No quantitative criteria for when constant-\(L\) approximation is valid.
- Geometric mean justification stated qualitatively only; now formalized via integral of \(\ln z\).
- Lacked discretization error metric; now \(E_i\) defined.
- No guidance for figure captions or diagnostics table.

All enhancements preserve original intent while adding formal rigor and actionable criteria.

## 21. Summary Table (Key Neutral Coefficients)
| Symbol | Description | Value |
|--------|-------------|-------|
| \(\alpha_m, \beta_m\) | Momentum closure coefficients | Typically \(\approx 0.45\)–\(0.55\), \(14\)–\(16\) |
| \(\alpha_h, \beta_h\) | Heat closure coefficients | Typically \(\approx 0.45\)–\(0.55\), \(14\)–\(16\) |
| \(a_m, b_m\) | Quadratic SBL coefficients for momentum | \(6.3\) to \(8.8\), \(64\) to \(110\) |
| \(a_h, b_h\) | Quadratic SBL coefficients for heat | \(6.3\) to \(8.8\), \(64\) to \(110\) |
| \(\Delta\) | Near-neutral curvature coefficient | \(\alpha_h \beta_h - 2 \alpha_m \beta_m\) |
| \(c_1\) | Cubic SBL coefficient | \(\alpha_h \beta_h^2 - 2 \alpha_m \beta_m^2\) |

Note (precedence). The explicit surface-layer composition \(Ri_g(\zeta)=\zeta\,\phi_h(\zeta)/\phi_m(\zeta)^2\) appears explicitly in England & McNider (1995, BLM). Earlier MOST works define the ingredients but may not show the compact equality in print; see the curvature note for verification steps.

## 22. Parameter Selection and Curvature Scaling
// renamed from 21A to maintain monotonic section numbering
Recommended near-neutral parameter ranges for MOST:
- \(\alpha_m, \alpha_h\): \(0.45\) to \(0.55\)
- \(\beta_m, \beta_h\): \(14\) to \(16\)

Curvature scaling:
- Compute \(\Delta\) using selected \(\alpha, \beta\).
- For initial \(Ri_g\) estimates, use \(2\Delta\) from above.

## 23. Enhanced MOST / Ri Formulations (Model-Oriented)
For model compatibility and improved stability:
- Implement quadratic or rational function forms for φ_m, φ_h based on \(Ri\) or \(Ri_g\).
- Ensure continuity and smoothness at stability function transition points.

## 24. Appendix—Planetary Scaling and Polar Use
Planetary boundary layer (PBL) adaptation:
- Retain curvature formulation but adjust \(L\) scaling to match planetary conditions.
- Near-neutral coefficients \(\alpha, \beta\) may differ; calibrate using local PBL data.

Polar application notes:
- Be cautious of stability function applicability near poles; curvature effects may differ.
- Consider blending with alternative formulations in extreme latitudes.

## 25. Quadratic SBL Truncation (Q‑SBL)
// (keep single Q‑SBL section; moved earlier duplicate to this consolidated position)
### 25A. Typical Q‑SBL Coefficients (Stable SBL, ζ ≥ 0)
Given α,β pairs in the common range α_{m,h} ≈ 0.45–0.55, β_{m,h} ≈ 14–16, the Q‑SBL coefficients
a_m = α_m β_m, b_m = 0.5 α_m(α_m + 1) β_m^2, and similarly for heat, are typically:
- a_m, a_h ≈ 6.3 to 8.8
- b_m, b_h ≈ 64 to 110

Ready‑to‑use examples
- Symmetric “BD near‑neutral” set: α_m = α_h = 0.50, β_m = β_h = 16
  - a_m = a_h = 8.00, b_m = b_h = 96.0
  - Δ = α_h β_h − 2 α_m β_m = 8 − 16 = −8
  - c1 = α_h β_h^2 − 2 α_m β_m^2 = 128 − 256 = −128
- Cross (momentum smaller, heat larger): α_m = 0.45, β_m = 15; α_h = 0.55, β_h = 15
  - a_m = 6.75, b_m ≈ 73.1; a_h = 8.25, b_h ≈ 95.9
  - Δ = 8.25 − 2·6.75 = −5.25
  - c1 = 0.55·15² − 2·0.45·15² = −78.75

Operational guardrails
- Use Q‑SBL for 0 ≤ ζ ≤ ζ_max with ζ_max ≈ 0.2–0.3 (site dependent).
- Optional linear cap c_m ≈ c_h ≈ 5 if extreme stability segments require bounded growth.
- For diagnostics, report neutral curvature 2Δ and whether ζ_inf (if any) lies below ζ_max.

### 25B. Pure Richardson-Number Closures f_m(Ri), f_h(Ri)
Goal: express φ_m, φ_h as functions of Ri only, avoiding explicit ζ. Two practical approaches:

1) Near‑neutral series (Ri ≲ 0.1)
- From Section 25, φ_m(ζ) ≈ 1 + a_m ζ + b_m ζ²; φ_h(ζ) ≈ 1 + a_h ζ + b_h ζ².
- Invert ζ(Ri) to O(Ri³): ζ ≈ Ri − Δ Ri² + (1.5 Δ² − 0.5 c1) Ri³.
- Compose to O(Ri²):
  - f_m(Ri) = φ_m(ζ(Ri)) ≈ 1 + s1_m Ri + s2_m Ri²
  - f_h(Ri) = φ_h(ζ(Ri)) ≈ 1 + s1_h Ri + s2_h Ri²
  - Coefficients:
    - s1_m = a_m, s1_h = a_h
    - s2_m = b_m − a_m Δ, s2_h = b_h − a_h Δ

Example (α=0.5, β=16): a_m = a_h = 8, b_m = b_h = 96, Δ = −8
- s2_m = 96 − 8·(−8) = 160 → f_m(Ri) ≈ 1 + 8 Ri + 160 Ri²
- s2_h = 96 − 8·(−8) = 160 → f_h(Ri) ≈ 1 + 8 Ri + 160 Ri²
Validity: accurate for small Ri; clip or switch to rational form for Ri ≳ 0.1.

2) Pade [1/1] rational approximation (wider Ri, still near‑neutral)
Fit a stable rational form to match the series up to O(Ri²):
- f(Ri) ≈ (1 + p Ri) / (1 − q Ri), where:
  - p − q = s1
  - q (p + q) = s2
- Solve q = 0.5 (−s1 + sqrt(s1² + 4 s2)), p = s1 + q. Apply separately to m,h.

Example (α=0.5, β=16): s1 = 8, s2 = 160 → q ≈ 9.2665, p ≈ 17.2665
- f_m(Ri) ≈ (1 + 17.27 Ri) / (1 − 9.27 Ri)
- f_h(Ri) ≈ (1 + 17.27 Ri) / (1 − 9.27 Ri)
Notes:
- Pole at Ri ≈ 1/q ≈ 0.108: use only below that or blend to a capped form.
- For stronger stability or larger Ri, prefer ζ inversion or a capped surrogate.

Optional “power” alternative (match s1,s2):
- f(Ri) ≈ (1 + q Ri)^(−e), with e = −s1/q and 0.5 e(e+1) q² = s2.
- For α=0.5, β=16, q = (s1² − 2 s2)/s1 = −32, e = 0.25.
- Growth is rapid; pole at Ri = −1/q = 0.03125 → too restrictive in practice.

Recommended procedure
- Compute (a_m, b_m, a_h, b_h, Δ, c1) from α,β.
- Build series closures (s1, s2) for m,h.
- If Ri ≤ 0.05–0.1, use series; else use Pade [1/1] with pole guard (e.g., cap Ri < 0.8/q).
- For full‑range fidelity: use exact ζ(Ri) inversion (series seed + 1–2 Newton steps) and evaluate φ_{m,h}(ζ).

Implementation hint
- See profiles.py helpers ri_closure_series and ri_closure_pade for ready callables that return f_m(Ri), f_h(Ri) directly.

## 27. Cross-Reference
For full derivation details of \(d^{2}Ri_g/d\zeta^{2}\) and inversion series see: “Richardson number curvature.md”. This note emphasizes discretization, mapping, and surrogate blending.
