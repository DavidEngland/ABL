# Guest Lecture Outline — MOST, Stable BL, and Curvature-Aware Corrections

Audience: Graduate students and scientists (STEM); goal is intuition → equations → implications → open problems.

Time plan (60–90 min)
- 0–5 min: Motivation and objectives
- 5–20 min: MOST primer (intuitive → formal)
- 20–35 min: Gradient Richardson number Ri_g and curvature
- 35–50 min: Grid sensitivity and curvature-aware correction
- 50–65 min: Results summary and validation
- 65–80 min: Open questions, future work, discussion
- 80–90 min: Live demo (optional) and Q&A

Learning objectives
- Understand the role of similarity variables and Obukhov length L in the surface layer.
- See how Ri_g(ζ) curvature explains coarse-grid bias in stable BLs.
- Learn neutral-curvature invariance (2Δ) and why preserving it matters.
- Review practical correction strategies and where validation is still needed.

---

## 0. Motivation (Slide 1–2)
- Arctic and nocturnal stable layers are grid-sensitive; coarse Δz misclassifies stability and biases turbulent fluxes.
- We need a resolution-aware fix that preserves near-neutral physics and reduces curvature-induced bias.

Key claim (working result)
- A neutral-curvature–preserving modifier reduces coarse-grid curvature error by 40%+ without ad hoc diffusion floors.

---

## 1. MOST Primer (Slide 3–6)

Physical intuition
- Near-surface exchange driven by shear and buoyancy; dimensional groups collapse variability when scaled properly.

Core definitions
- Obukhov length: \(L=-\frac{u_*^3\theta}{\kappa g\,\overline{w'\theta'}}\) (sign sets stable/unstable).
- Dimensionless height: \(\zeta=z/L\) or \(\zeta=z/L(z)\) if local.
- Stability functions:
  - \(\phi_m=\frac{\kappa z}{u_*}\frac{\partial U}{\partial z}\), \(\phi_h=\frac{\kappa z}{\theta_*}\frac{\partial \theta}{\partial z}\), \(\theta_*=-\overline{w'\theta'}/u_*\).
  - Example stable branch: \(\phi_{m,h}=(1-\beta_{m,h}\zeta)^{-\alpha_{m,h}}\) with domain guard \(1-\beta\zeta>0\).

Gradient Richardson number (MOST form)
- $Ri_g(\zeta)=\zeta\,\phi_h/\phi_m^2=\zeta F(\zeta)$.

Notes for students
- L can vary with height and time in real SBLs; start with bulk L and add local-L mapping when needed.

---

## 2. Curvature of Ri_g and Why It Matters (Slide 7–10)

Compact curvature formula
- $ \dfrac{d^2 Ri_g}{d\zeta^2}=F\big[2V_{\log}+\zeta\,(V_{\log}^2-W_{\log})\big]$, with $F=\phi_h/\phi_m^2$, $V_{\log}=(\phi_h'/\phi_h)-2(\phi_m'/\phi_m)$, $W_{\log}=dV_{\log}/d\zeta$.

Neutral curvature invariant
- $\Delta=\alpha_h\beta_h-2\alpha_m\beta_m$, and $\left.\dfrac{d^2 Ri_g}{d\zeta^2}\right|_0=2\Delta$.
- Interpretation: sign sets initial concavity; magnitude sets strength of early departure from linearity.

Practical impact
- For Δ<0 (typical SBL), Ri_g bends down quickly; layer-averaging then underestimates stability at the lowest level → overmixing.

NEW (Slide visual cue)
- Show a single panel with three thin curves sharing the same initial slope (tangent to Ri_g = ζ):
  - Δ = 0: straight reference line (linear).
  - Δ < 0: concave-down (typical SBL) — highlight early drop below ζ line.
  - Δ > 0: concave-up (rare) — rises above ζ line.
- Spoken emphasis: “Δ sets the *rate of first departure* from neutrality; preserving 2Δ anchors physics at ζ → 0.”

Optional quick plot (Python)
```python
zeta = np.linspace(0,0.25,200)
def ri(z,a): return z + a*z*z   # quadratic near-neutral sketch only
plt.plot(zeta, ri(zeta, 0.0), label='Δ=0')
plt.plot(zeta, ri(zeta,-0.8), label='Δ<0 (concave-down)')
plt.plot(zeta, ri(zeta, 0.6), label='Δ>0 (concave-up)')
plt.plot(zeta, zeta, 'k--', lw=0.7, label='Ri_g = ζ')
plt.xlabel('ζ'); plt.ylabel('Ri_g'); plt.legend()
```

Speaker note
- “Everything downstream (bias, correction) originates from this curvature contrast right after neutrality.” 

---

## 3. Grid Sensitivity and Correction (Slide 11–14)

Observation
- Coarse Δz (10–100 m) smears near-surface curvature; bulk Ri_b < point Ri_g at \(z_g=\sqrt{z_0 z_1}\).

NEW clarification: geometric mean height
- Use \(z_g\) because for power-law / logarithmic-like vertical structure, integrals over [z0,z1] are best represented by exp(average ln z) → geometric mean.
- Concavity logic: If \(Ri_g(z)\) is concave-down over [z0,z1], then by Jensen:
  \[
  Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_g)
  \]
  ⇒ \(Ri_b\) underestimates true local stability ⇒ turbulent mixing coefficients too large ⇒ overmixing.

Strategy
- Preserve neutral curvature (2Δ), adjust only tail behavior (ζ>0) with a grid-aware modifier:
  - \(f_c(\zeta,\Delta z)=\exp[-D (\zeta/\zeta_r)(\Delta z/\Delta z_r)]\), choose exponents so \(V_{\log}(0)\) unchanged.

Alternatives
- Q‑SBL (quadratic surrogate) for ζ≤0.2–0.3; Ri-based closures via series + Newton; regularized power law to avoid poles.

Simplified correction framing (slide-friendly)
Goal
- Adjust \(K_{m,h}\) on coarse grids without altering neutral entry curvature (2Δ).

Introduce generic damping factor
\[
K_{m,h}^* = K_{m,h} \times G(\zeta,\Delta z)
\]
Constraints
1. \(G(0, \Delta z)=1\) (preserve 2Δ).
2. \(\partial_\zeta G|_{\zeta=0}=0\) (do not perturb first derivative → keeps neutral Taylor series).
3. \(G \rightarrow 1\) as \(\Delta z \rightarrow 0\) (grid convergence).
4. \(G\) monotone non-increasing in ζ for fixed coarse Δz (only damp tails).
Minimal functional template
\[
G(\zeta,\Delta z)=\exp\!\left[-D\left(\frac{\Delta z}{\Delta z_r}\right)^p \left(\frac{\zeta}{\zeta_r}\right)^q\right],\quad p,q>0
\]
Choose q ≥ 1 so \(\partial_\zeta G|_{0}=0\); tune D,p,q from target bias reduction (e.g., 40% curvature error cut at Δz = 60–100 m).

Speaker wording
- “Not modifying neutrality; only suppressing excessive tail influence introduced by coarse vertical averaging.”

Bias diagnostic (show on slide)
\[
B = \frac{Ri_g(z_g)}{Ri_b};\quad B>1 \text{ signals curvature-induced underestimation.}
\]

---

## 4. Results Snapshot (Slide 15–17)

- Neutral curvature preserved within <5% across tested grids.
- Curvature error reduced by 40%+ at Δz=60–100 m (idealized cases).
- ζ(Ri) inversion: series seed + 1 Newton step achieves machine-precision equality with low cost.
- Diagnostics: amplification ratio A(ζ), inflection height, omission metric for variable L(z).

Caveats
- Parameter transfer from USL fits can worsen SBL curvature; prefer SBL-calibrated or surrogate forms.
- Guard poles: restrict ζ<~0.7/β or switch to surrogate.

---

## 5. Dynamic Critical Richardson Number & Practical Options (talking points)
- Motivation: fixed Ri_c (≈0.25) misses hysteresis, intermittent turbulence and inversion-strength dependence — propose dynamic Ri_c* that adapts to local inversion strength, shear, and turbulence memory.
- Two operational choices when Ri exceeds threshold:
  1. Modify mixing-length l → l* = l · L_mod(Ri,Ri_c*). Reduces eddy scale directly (preferred for McNider-style slope/terrain work).
  2. Modify diffusivity K → K* = K · F_mod(Ri,Ri_c*). Simpler multiplier (preferred for Biazar-style air-quality/model-integration).
- Canonical dynamic Ri_c* prototype (lecture-friendly):
  - Ri_c* = Ri_c0 + α_inv · min(Γ/Γ_ref,1) + β_mem · TKE_rel, where Γ = lapse/inversion strength, TKE_rel ∈ [0,1].
- Recommended default actions:
  - If intermittent (TKE low, strong inversion): reduce l first, then apply K damping.
  - If operational cost dominant: apply K multiplier (exponential in Ri/Ri_c*).

## 6. Jensen, Bulk Ri vs Gradient Ri — Estimation Techniques (slide bullets)
- Jensen: for concave-down Ri_g(z), layer average Ri_b < Ri_g(z_g) where z_g = √(z0 z1). Use this to explain coarse-grid underestimation bias.
- Representative heights:
  - Use geometric mean z_g for point evaluations of log/power-law profiles.
  - Use logarithmic mean z_L = (z2−z1)/ln(z2/z1) when matching ΔU exactly.
- Recommended finite-difference estimators:
  - Centered gradient (interior): (U_{k+1}-U_{k-1})/(z_{k+1}-z_{k-1}).
  - First layer: forward difference with geometric/log mean for height.
  - Bulk Ri_b: trapezoid or Simpson on Ri_g(z) if full profile available.
- Quick pseudocode (slide-ready):
```python
# compute Ri_b and Ri_g at geometric mean
z_g = sqrt(z0*z1)
Ri_g_zg = z_g/L * phi_h(z_g/L) / phi_m(z_g/L)**2
Ri_b = (g/theta_ref)*(theta1-theta0)*(z1-z0) / ((U1-U0)**2)
B = Ri_g_zg / Ri_b  # Jensen bias >1 indicates underestimation
```

## 7. Roles — McNider & Biazar (lecture note)
- McNider: lead on dynamic Ri_c derivation, slope/terrain modifications, mixing-length based interventions and collapse/LLJ studies.
- Biazar: lead on K-based multiplier design, operational integration (WRF/CMAQ), urban/air-quality validation and remote-sensing assimilation.
- Joint tasks: tuning Ri_c* function, shared validation on tower/LES cases.

---

## 8. Live Demo (5–10 min, optional)

- Demo A: ζ(Ri) inversion accuracy (series + Newton) and pole guard behavior.
- Demo B: Curvature profile vs Δz grids; show reduction with neutral-preserving modifier.
- Demo C: Quick variable-L(z) mapping and omission metric E_omit thresholding.

Refined execution plan
- Demo B (Core): Single figure with three curves for one case:
  1. Fine-grid reference \(Ri_g^{fine}(z)\).
  2. Coarse-grid reconstructed (layer-averaged) profile (biased).
  3. Corrected coarse \(Ri_g^*(z)\) after applying \(G\).
  Add inset: table with (Δz, B_before, B_after).
- Demo C (If time): Plot \(E_{\text{omit}}(z)\) for a case with variable \(L(z)\); gray band where \(E_{\text{omit}}<0.05\) (safe shortcut region).
- Keep code hidden; narrate physical interpretation (“We recover local curvature signature without touching 2Δ.”).

Slide micro-labels
- “Anchor: neutrality (2Δ)”
- “Problem: concave-down ⇒ Ri_b deficit”
- “Fix: tail damping G(ζ,Δz)”

Time guard
- Abort Demo C if Demo B + questions exceed 6 min.

---

## Appendix A — Copilot Q&A: Canonical derivation & drop‑in curvature correction (for slides)

- Why a correction term appears (short): turbulent eddies sample a finite vertical extent ℓ(z); if Ri(z) has curvature the eddy-averaged stability differs from the point value and the closure must include a second‑order correction proportional to ℓ^2.

- Kernel-average result (one slide):
  - Effective averaged stability: ⟨F(Ri)⟩_W = F(Ri) + (M2/2)[ F'(Ri) Ri'' + F''(Ri) (Ri')^2 ] + O(ℓ^3)
  - With M2 ≈ α ℓ^2 (α depends on kernel; use α∈[1/12, 1/2]).

- Drop‑in diffusivity correction (single equation for slides/code):
  - K_eff(z) = K_0(z) · ⟨F(Ri)⟩_W
  - ≈ K_0(z) [ F(Ri) + (α ℓ^2 / 2) ( F'(Ri) Ri'' + F''(Ri) (Ri')^2 ) ]

- Practical expressions using MOST (one bullet each):
  - Ri = ζ r(ζ) with r = φ_h / φ_m^2 and ζ = z/L (local-L terms add chain-rule corrections).
  - Ri' = (1/L) [ r + ζ r' ], Ri'' = (1/L^2) [ 2 r' + ζ r'' ] for locally-constant L.
  - Mixing length estimate: ℓ ≈ κ z / φ_m(ζ).

- Implementation notes (short):
  1. Evaluate r, r', r'' from chosen φ_m/φ_h (analytically if power-law, or numerically).
  2. Choose kernel constant α (default 1/12 or 1/2); compute ℓ and M2 = α ℓ^2.
  3. Compute F, F', F'' (F = φ_h / φ_m^2).
  4. Compute curvature terms Ri', Ri'' and assemble K_eff via formula above.
  5. Use K_eff in place of K_0·F for flux computations or as multiplicative correction G = K_eff / (K_0 F).

- Slide-friendly phrasing for the lecture:
  - "Eddies average stability; curvature yields a second‑order correction ∝ ℓ^2. Use K_eff = K_0·(F + correction) as a principled, non‑ad hoc fix for coarse Δz."
  - Provide the one-line drop-in formula and a 3-line pseudocode snippet on the slide.

- Demo tip: implement formula with symbolic/automatic-differentiation for r', r'' (SymPy or small finite difference) and show before/after K profiles for Demo B.

## Appendix B — Short code pseudocode (for handout)
```python
# quick reference pseudocode
# inputs: z, L (or L(z)), phi_m(zeta), phi_h(zeta), K0(z)
zeta = z / L
phi_m = phi_m_func(zeta); phi_h = phi_h_func(zeta)
r = phi_h / (phi_m**2)
# r', r'': analytic or small-h central difference in zeta
r1 = d_dzeta(r, zeta); r2 = d2_dzeta(r, zeta)
Ri = zeta * r
Ri_p = (r + zeta * r1)/L
Ri_pp = (2*r1 + zeta * r2)/(L**2)
F = r
F1 = dF_dRi(F, Ri)   # or chain via zeta derivatives
F2 = d2F_dRi2(F, Ri)
ell = kappa * z / phi_m
alpha = 1.0/12.0  # kernel constant
K_eff = K0 * ( F + 0.5 * alpha * ell**2 * ( F1*Ri_pp + F2*(Ri_p**2) ) )
# use K_eff in flux calculation
```

## Quick answer: φ‑agnostic curvature corrections (slide / FAQ)

- Yes — curvature corrections can be applied without knowing the model's internal φ(ζ) form by operating on diagnosed quantities (Ri_g, Ri_b, z_g) and using neutral‑preserving modifiers or simple surrogates.
- Minimal φ‑agnostic recipe (for slides):
  1. Diagnose point Ri_g at the geometric mean z_g (or compute Ri_b and Ri_g(z_g)).  
  2. Compute bias ratio B = Ri_g(z_g)/Ri_b. If B ≲ 1.05 do nothing. If B > threshold (e.g., 1.1), apply correction.
  3. Apply multiplicative damping to diffusivities: K* = K · G(ζ,Δz) with G from Section 3. Choose q≥2 so G′(0)=0 to preserve 2Δ.
  4. Optionally apply mixing length reduction: l* = l / (1 + a_l (Ri/Ri_c*)^n) when you cannot modify K directly.
- Excel / spreadsheet quick formulas (one cell per quantity):
  - z_g = SQRT(z0*z1)
  - Ri_g_zg = (z_g / L) * phi_h(z_g/L) / (phi_m(z_g/L)^2)  // if φ unavailable compute Ri_g_zg from local gradients instead
  - Ri_b = (g/theta_ref)*(theta1 - theta0)*(z1 - z0)/( (U1-U0)^2 )
  - B = Ri_g_zg / Ri_b
  - G = EXP( -D * (Δz/Δz_ref)^p * (ζ/ζ_ref)^q )
  - K_star = K * G
