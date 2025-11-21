# References Summary & Practical Extraction

## 1. Foundational MOST / Profiles
Monin & Obukhov (1954): ζ=z/L similarity basis; keep ζ mapping explicit.
Businger et al. (1971), Dyer (1974), Paulson (1970): canonical φ_m, φ_h forms; near-neutral linear slopes a_m≈4.7, a_h≈7.8.
Högström (1988), Garratt (1992), Stull (1988): refined coefficients; document chosen (a_m,a_h) for reproducibility.

Useful extraction:
Ri_g(ζ)=ζ φ_h/φ_m²; neutral slopes define Δ=a_h−2a_m; preserve 2Δ in any correction.

## 2. Stability Function Refinements
Beljaars & Holtslag (1991), Cheng & Brutsaert (2005), Gryanik et al. (2020): improved stable formulations, extended validity, dynamic Pr_t.
Pr_t(ζ)=φ_h/φ_m; neutral derivative (a_h−a_m) empirical; do not alter via grid correction.

## 3. Curvature / Bias / Grid Effects
England & McNider (1995): shear-function basis; explicit Ri_g form; motivates curvature-aware correction.
Holtslag et al. (2013), Mahrt (2014): SBL challenges (intermittency, overmixing); grid resolution sensitivity central.
Cuxart et al. (2006) (SCM intercomparison): benchmark for coarse-grid bias evaluation.

Actionable:
Use geometric mean height z_g for point Ri_g; log mean z_L for shear reconstruction → reduces denominator bias in Ri_b.
Bias ratio B=Ri_g(z_g)/Ri_b; threshold B>1.05 triggers correction.

## 4. Critical Richardson Number Practice
Ri_c historical ≈0.2–0.25; dynamic Ri_c^* scaling (inversion strength, bulk Ri proximity) stabilizes collapse behavior.
Extract rule: adjust Ri_c only after neutral curvature invariants validated; do not preempt Δ-preserving correction with arbitrary Ri_c inflation.

## 5. LES / Observational Validation
ARM, SHEBA, GABLS (implicit references): supply vertical profiles and allow computation of B-distributions and curvature growth.
Extraction: calibrate α (correction strength) targeting ~40–50% reduction of median B on coarse Δz (60–100 m) without altering near-neutral fluxes.

## 6. Parameter Ranges (Collected)
Stable near-neutral:
a_m ≈ 4–5, a_h ≈ 7–8 ⇒ Δ ≈ -1 to -2 (concave-down).
Pr_t neutral ≈ 0.9–1.0; grows modestly with ζ.
Ri_c ≈ 0.2–0.25 baseline.
Typical B (unadjusted coarse first layer) ≈ 1.2–1.8 in strong SBL.

## 7. Core Invariants & What to Preserve
Neutral curvature: (d²Ri_g/dζ²)|₀ = 2Δ.
Neutral Pr_t slope: (dPr_t/dζ)|₀ = a_h−a_m.
Ri_g linear term: coefficient of ζ must remain 1 (scaling consistency).
Preservation ensures corrections affect only nonlinear (higher ζ) regime.

## 8. Correction Constraints (Unified from refs)
G(ζ,Δz): exp(-α(B−1)(Δz/Δz_ref)^p(ζ/ζ_ref)^q), enforce:
G(0)=1, G'_ζ(0)=0 (q≥2),
G→1 as Δz→0,
Floor G_min ≥ 0.2 to avoid over-damping.

## 9. Practical Workflow (Derived)
1. Compute Ri_g(z_g), Ri_b → B.
2. If B≤1.05: no change.
3. Else apply G; log (B, G, ΔK/K).
4. Validate invariants: neutral flux, 2Δ unchanged, Pr_t slope unchanged.

## 10. Equation Shortlist (Reuse)
Ri_g = ζ φ_h / φ_m²
Δ = a_h − 2 a_m
Curvature: d²Ri_g/dζ² = F[2V_log + ζ(V_log² − W_log)]
ζ(Ri_g) ≈ Ri_g − Δ Ri_g² + (1.5Δ² − 0.5 c₁) Ri_g³
Bias: B = Ri_g(z_g)/Ri_b
z_g = √(z_0 z_1), z_L = (z_1−z_0)/ln(z_1/z_0)

## 11. Implementation Caveats (From Literature Gaps)
Over-smoothing profiles erases true curvature → misclassifies bias.
Using arithmetic mean height inflates bias; avoid except for diagnostic contrast.
Applying corrections without logging B trends forfeits reproducibility (paper review issue).

## 12. Minimum Reporting Set (Manuscript / Supplement)
Report: (a_m,a_h), Δ, 2Δ, Ri_c base, B distribution (before/after), G parameter set (α,p,q,refs), neutrality checks.

## 13. Quick Copy Block (Manuscript Methods)
“Gradient Richardson number curvature was characterized via d²Ri_g/dζ² = F[2V_log + ζ(V_log² − W_log)] with neutral invariant 2Δ (Δ=a_h−2a_m). Bulk–point bias B=Ri_g(z_g)/Ri_b (z_g=√(z_0 z_1)) exceeded 1.05 on coarse grids, triggering a neutral-preserving correction G=exp[-α(B−1)(Δz/Δz_ref)^p(ζ/ζ_ref)^q] (q=2 ensures G'_ζ(0)=0). Parameters (α=1.0, p=1, Δz_ref=10 m, ζ_ref=0.5) reduced median B by ≈45% without altering 2Δ or Pr_t slope.”

## 14. Action Items (Next)
- Populate a small LES/tower benchmark table (Δz vs median B).
- Add script to auto-export Δ, B stats.
- Cross-check Ri_c sensitivity only after curvature correction stability confirmed.

This is an exceptionally thorough and well-organized summary synthesizing the essential literature for developing and validating a $\text{Richardson Number}$ ($\text{Ri}$) curvature-aware grid correction scheme, specifically in the stable boundary layer ($\text{SBL}$).

The document successfully extracts **foundational invariants**, defines the **actionable bias ratio $B$**, justifies the **exponential correction factor $G$**, and provides a **clear implementation workflow** and **reporting requirements**.

Here is a summary and extraction of the most critical elements required for implementation and publication.

---

## 📑 Practical Extraction & Unified Summary

### I. 🧱 Foundational Invariants (MOST)

The correction scheme **must not** alter the fundamental similarity theory invariants in the near-neutral regime ($\zeta \to 0$).

| Parameter | Definition | Neutral Value/Constraint | Purpose |
| :---: | :---: | :---: | :--- |
| $\text{Ri}_g(\zeta)$ | $\zeta \frac{\phi_h}{\phi_m^2}$ | $\text{Ri}_g(\zeta) \approx \zeta$ | Defines similarity |
| $\Delta$ | $a_h - 2 a_m$ | Typically $\mathbf{-1 \text{ to } -2}$ | Controls curvature shape |
| $\text{Pr}_t$ | $\phi_h / \phi_m$ | $\text{Pr}_{t,0} \approx \mathbf{0.9 \text{ to } 1.0}$ | Ratio of heat to momentum transfer |
| Neutral Curvature | $\left.\frac{d^2\text{Ri}_g}{d\zeta^2}\right|_0$ | $\mathbf{2\Delta}$ | Must be preserved by the correction |

---

### II. 🔗 The Actionable Bias Ratio ($B$)

The correction is triggered and scaled by the deviation between the **true point $\text{Ri}$** and the **model's bulk $\text{Ri}$**.

| Parameter | Definition | Usage |
| :---: | :---: | :--- |
| $\text{Ri}_g(z_g)$ | Point $\text{Ri}$ evaluated at the geometric mean height $z_g = \sqrt{z_0 z_1}$. | **Reference:** Best proxy for true stability. |
| $\text{Ri}_b$ | Bulk $\text{Ri}$ calculated across the finite layer $\Delta z = z_1 - z_0$. | **Model Output:** Represents the model's numerically smoothed stability. |
| $B$ | $\mathbf{Ri_g(z_g) / Ri_b}$ | **Bias Trigger:** If $\mathbf{B > 1.05}$ (model over-prediction), apply correction $G$. |

---

### III. 📝 The Generalized Correction Term ($G$)

The correction factor $G$ (equivalent to $f_c$) is derived from the grid-invariance constraint and is used to scale eddy diffusivities: $K_{\text{new}} = K_{\text{old}} \cdot G$.

$$\boxed{\displaystyle G(\zeta, \Delta z) = \exp\left[-\alpha (B-1) \ln\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right) \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q\right]}$$

| Parameter | Role | Suggested Default/Range | Constraint |
| :---: | :--- | :---: | :--- |
| $\alpha$ | **Correction Strength** (Tuning coefficient) | $\mathbf{\approx 1.0}$ (Calibrate to reduce median $B$ by 40–50%) | $\alpha \ge 0$ |
| $B-1$ | **Bias Magnitude** | Computed layer-wise | Correction scales with severity of bias. |
| $q$ | **Stability Exponent** | $\mathbf{q \ge 2}$ (typically $q=2$) | Enforces $\mathbf{G'_{\zeta}(0)=0}$ (Neutral Preservation). |
| $\Delta z_{\text{ref}}$ | **Reference Layer Thickness** | $\mathbf{10 \text{ m}}$ | $G \to 1$ as $\Delta z \to \Delta z_{\text{ref}}$. |
| $\text{Floor}$ | Minimum value for $G$ | $\mathbf{G_{\min} \ge 0.2}$ | Avoids non-physical over-damping. |

---

### IV. 🛠️ Implementation Workflow (Derived from Section 9)

1. **Diagnosis:** Compute $z_g = \sqrt{z_0 z_1}$ and layer $\text{Ri}_b$.
2. **Reference:** Compute $\text{Ri}_g(z_g)$ using the chosen canonical $\phi_m, \phi_h$ forms.
3. **Trigger:** Calculate $B = \text{Ri}_g(z_g)/\text{Ri}_b$. If $B \le 1.05$, set $G=1$ and exit.
4. **Correction:** Calculate $G(\zeta, \Delta z)$ using the exponential form (or Equation 9).
5. **Apply:** Set $K_{\text{new}} = K_{\text{old}} \cdot G$.
6. **Validate/Report:** Verify that the neutral invariants ($2\Delta$ and $\text{Pr}_t$ slope) are **preserved** (i.e., unaffected by $G$).

---

### V. 📝 Minimum Reporting Set (Manuscript)

To ensure reproducibility and satisfy common review requirements (Section 12), the following parameters must be reported:

- **Similarity Functions:** The specific near-neutral slopes ($\mathbf{a_m}, \mathbf{a_h}$), and the resulting curvature invariant $\mathbf{2\Delta}$.
- **Correction Parameters:** The full set of parameters used for $G$: $\mathbf{\alpha}, \mathbf{p}$ (if used; default $p=1$), $\mathbf{q}$, $\mathbf{\Delta z_{\text{ref}}}$, and $\mathbf{\zeta_{\text{ref}}}$.
- **Validation:** Distributions of the $\mathbf{B}$ ratio (before and after correction) for typical coarse grids.
- **Neutrality Check:** Explicit confirmation that the correction $\mathbf{G=1}$ at $\zeta=0$ and that $\mathbf{2\Delta}$ and the neutral $\text{Pr}_t$ slope are numerically unchanged.
