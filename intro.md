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
- \(Ri_g(\zeta)=\zeta\,\phi_h/\phi_m^2=\zeta F(\zeta)\).

Notes for students
- L can vary with height and time in real SBLs; start with bulk L and add local-L mapping when needed.

---

## 2. Curvature of Ri_g and Why It Matters (Slide 7–10)

Compact curvature formula
- \(\frac{d^2 Ri_g}{d\zeta^2}=F[2V_{\log}+\zeta(V_{\log}^2-W_{\log})]\),
  with \(F=\phi_h/\phi_m^2\), \(V_{\log}=(\phi_h'/\phi_h)-2(\phi_m'/\phi_m)\), \(W_{\log}=dV_{\log}/d\zeta\).

Neutral curvature invariant
- \(\Delta=\alpha_h\beta_h-2\alpha_m\beta_m\), \( \left.\frac{d^2 Ri_g}{d\zeta^2}\right|_0=2\Delta\).
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

## 5. Open Questions and Areas Needing Work (Slide 18–21)

- Calibration of D (tail strength):
  - Functional form vs. Δz, Ri, and target amplification A(ζ)? Site dependence vs. universal scaling?
- Variable L(z):
  - When is constant-L mapping sufficient (E_omit thresholds) in real nocturnal layers?
- Ri-only closures:
  - Pade [1/1] vs. ζ inversion tradeoffs for operational stability and range; blending strategy near poles.
- Critical Richardson number:
  - Fixed 0.25 vs. dynamic Rc; reconcile with observed intermittent turbulence at Ri>1.
- ODE route:
  - Directly solving reduced ODEs for Ri with lower-boundary curvature constraints—benefits vs. cost?
- Data coverage:
  - Extreme stability regimes (polar, over ice) for robust calibration; how to handle sparse segments.

Discussion prompts (with possible answers)
- “Should Rc be fixed?” → Likely context-dependent; normalize s=Ri/Rc and let Rc vary with inversion strength.
- “Why preserve 2Δ?” → Anchors neutral physics; prevents “fixing” the tail by breaking the entry condition.
- “Q‑SBL vs power-law?” → Q‑SBL is safer near poles; power-law gives analytic leverage—blend pragmatically.

---

## 6. Future Directions (Slide 22–24)

Near-term
- 1D slope-flow notebook (katabatic/anabatic) for rapid testing of closures and Δz effects.
- LES/tower cross-validation (GABLS, ARM NSA) with shared diagnostics and bias tables.
- Julia acceleration path for large sensitivity sweeps and ERA5 climatology subsetting.

Medium-term
- Urban megacity Ri_g mapping (remote sensing blend) → grid-dependence thresholds.
- Planetary scaling demonstrations (Mars/Titan): reuse curvature logic with different L and thermodynamics.

Tooling and collaboration
- GitHub-first (issues/PRs/branches); Markdown/LaTeX equations; no Word math.
- Reproducible notebooks with pinned environments; tagged releases for student projects.

---

## 7. Live Demo (5–10 min, optional)

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

## 8. Takeaways (Slide 25)
- Preserve neutral curvature (2Δ); correct tails, not the physics at ζ→0.
- Use geometric mean height \(z_g\) and analytic curvature for consistent first-layer behavior.
- Diagnostics enable principled Δz choices and adaptive refinement flags.
- Grid Damping Factor G: satisfies invariance at ζ=0 while targeting curvature-driven bias aloft in first layer.

Reading (short list for class)
- Stull (1988), Högström (1988), Beljaars & Holtslag (1991)
- England & McNider (1995), Cheng & Brutsaert (2005)
- Holtslag et al. (2013), Cuxart et al. (2006)

---

# Introduction to Monin–Obukhov Surface-Layer Theory and Richardson Number Curvature (SBL Focus)

Audience: STEM background (math / physics / engineering) new to atmospheric boundary-layer similarity theory (MOST), with emphasis on the Stable Boundary Layer (SBL).  
Goal: Build from physical intuition → key definitions → why curvature of the gradient Richardson number matters on coarse vertical grids → motivation for curvature-aware correction.

---

## 1. Physical Setting: Surface Layer and Boundary Layer Regimes

Atmospheric flow near the ground (first few 10–200 m) exchanges momentum and heat with the surface through turbulent eddies.  
Two canonical regimes:

- Unstable / convective (USL): surface heating, vigorous mixing.
- Stable (SBL): surface cooling, turbulence suppressed, strong vertical gradients.

MOST (Monin–Obukhov Similarity Theory) provides dimensionless relationships for mean gradients in the **(constant-flux) surface layer**, typically lowest ~10% of the full boundary layer depth. In practice, SBL model grids can be coarse enough that near-surface curvature in stability metrics is smeared, biasing mixing parameterizations.

---

## 2. Fundamental Variables and Fluxes

Let:
- \(u, v\): horizontal wind components; speed \(U = (u^2+v^2)^{1/2}\).
- \(\theta\) (or \(\theta_v\)): potential (virtual) temperature.
- Surface turbulent (Reynolds) fluxes:
  - Momentum: \(\tau = -\rho \overline{u'w'} = \rho u_*^2\) → defines friction velocity \(u_*\).
  - Heat: \(H = -\rho c_p \overline{w' \theta'}\) (stable SBL: typically negative at night).
- Turbulent kinematic heat flux: \(\overline{w'\theta'} = H / (\rho c_p)\).

Von Kármán constant: \(\kappa \approx 0.4\).

---

## 3. Obukhov Length \(L\): Definition and Physical Meaning

Canonical (bulk) form:
\[
L = -\frac{u_*^3\,\theta_{\text{ref}}}{\kappa g\,\overline{w' \theta'}}.
\]

Interpretation:
- Sign:  
  - Unstable (surface heating): \(\overline{w'\theta'}>0\) ⇒ \(L<0\).  
  - Stable (surface cooling): \(\overline{w'\theta'}<0\) ⇒ \(L>0\).
- Magnitude: ratio of mechanical (shear) production to buoyancy (stability) effects; smaller \(|L|\) → stronger buoyancy influence.

### 3.1 Local vs Bulk and Height / Time Variability
While classical MOST uses a single (layer-constant) \(L\), in real SBLs:
- Fluxes can vary with height → **local** Obukhov length \(L(z)\).
- Surface forcing evolves (nighttime cooling, transient clouds) → time-dependent \(L(t)\).
- For strongly stratified shallow layers, \(L(z)\) variation can matter when mapping curvature from dimensionless height to physical height.

We denote \(\zeta = z / L\) (bulk) or \(\zeta = z / L(z)\) (local). Small \(\zeta\) (near-neutral), larger \(\zeta\) (stronger stratification for stable).

---

## 4. Similarity Functions and Dimensionless Gradients

MOST postulates that non-dimensionalized mean gradients depend only on \(\zeta\).  
Define stability (or “correction”) functions for momentum (m) and heat (h):

\[
\phi_m(\zeta) = \frac{\kappa z}{u_*}\frac{\partial U}{\partial z},\qquad
\phi_h(\zeta) = \frac{\kappa z}{\theta_*}\frac{\partial \theta}{\partial z},
\quad \theta_* = -\frac{\overline{w'\theta'}}{u_*}.
\]

In neutral conditions (\(\zeta=0\)): \(\phi_m=\phi_h=1\).  
Stable example family (power-law form with a finite “pole”):
\[
\phi_m = (1 - \beta_m \zeta)^{-\alpha_m},\quad
\phi_h = (1 - \beta_h \zeta)^{-\alpha_h},\quad (1 - \beta_{m,h}\zeta) > 0.
\]

Near-neutral Taylor expansion (generic smooth \(\phi\)):
\[
\phi(\zeta) \approx 1 + a\zeta + b\zeta^2 + \dots
\]
with \(a=\alpha\beta\), \(b=\tfrac12\alpha(\alpha+1)\beta^2\) for the power-law example.

---

## 5. Richardson Numbers: Gradient vs Bulk

### 5.1 Gradient Richardson Number (local)
\[
Ri_g = \frac{(g/\theta)\, \partial \theta / \partial z}{\left(\partial U / \partial z\right)^2}.
\]
Using MOST gradients:
\[
Ri_g(\zeta) = \zeta\frac{\phi_h}{\phi_m^2} = \zeta F(\zeta),\quad F(\zeta) = \frac{\phi_h}{\phi_m^2}.
\]

### 5.2 Bulk Richardson Number (layer averaged)
Across first layer \(z_0\)–\(z_1\):
\[
Ri_b = \frac{g}{\theta}\frac{(\theta_1 - \theta_0)(z_1 - z_0)}{(U_1 - U_0)^2}.
\]
On coarse vertical grids \(Ri_b\) underestimates local (point) stability when \(Ri_g(\zeta)\) is strongly **concave-down** near the surface (common in SBL), leading to **overmixing**.

---

## 6. Why Curvature of \(Ri_g(\zeta)\) Matters

Near the surface:
\[
Ri_g(\zeta) = \zeta + \Delta \zeta^2 + \tfrac12(\Delta^2 + c_1)\zeta^3 + \dots
\]
where
\[
\Delta = \alpha_h\beta_h - 2\alpha_m\beta_m,\qquad
c_1 = \alpha_h\beta_h^2 - 2\alpha_m\beta_m^2.
\]

Second derivative (“curvature”) at neutral limit:
\[
\left.\frac{d^2 Ri_g}{d\zeta^2}\right|_{0} = 2\Delta.
\]

If \(\Delta < 0\) (typical stable sets), \(Ri_g\) bends downward rapidly: averaging over a thick layer lowers the apparent stability ⇒ excessive turbulent diffusion.

Full curvature expression (generic differentiable \(\phi\)):
\[
\frac{d^2Ri_g}{d\zeta^2} = F\left[2V_{\log} + \zeta\left(V_{\log}^2 - W_{\log}\right)\right],
\]
with
\[
V_{\log} = \frac{\phi_h'}{\phi_h} - 2\frac{\phi_m'}{\phi_m},\quad
W_{\log} = \frac{dV_{\log}}{d\zeta}.
\]

---

## 7. Parameter Transfer Caution (USL → SBL)

Many empirical \((\alpha,\beta)\) fits originate from **unstable** or weakly stable, near-neutral data (USL dominance). Directly applying those to **strong SBL**:
- Exaggerates curvature (large \(|\Delta|\)).
- Reduces usable \(\zeta\) range (pole near \(\zeta = 1/\beta\)).
- Increases sensitivity to vertical resolution.

Recommendation:  
(1) Fit or adjust coefficients using stable segments (low-level nocturnal data, LES).  
(2) Consider **quadratic surrogate (Q‑SBL)** forms for \(\zeta\) up to ~0.2–0.3 to avoid pole artifacts.

---

## 8. Variable Obukhov Length \(L(z)\) Mapping

If \(L\) varies with height:
\[
\zeta(z) = \frac{z}{L(z)},\quad
\frac{d\zeta}{dz} = \frac{L - z L'}{L^2},\quad
\frac{d^2\zeta}{dz^2} = -\frac{2L'}{L^2} - \frac{zL''}{L^2} + \frac{2z(L')^2}{L^3}.
\]
Map curvature:
\[
\frac{d^2 Ri_g}{dz^2} = \left(\frac{d\zeta}{dz}\right)^2 \frac{d^2 Ri_g}{d\zeta^2} + \frac{d^2\zeta}{dz^2}\frac{d Ri_g}{d\zeta},
\quad
\frac{d Ri_g}{d\zeta} = F(1 + \zeta V_{\log}).
\]

Diagnostic omission metric (decide if constant \(L\) shortcut is adequate):
\[
E_{\text{omit}} = \left|\frac{(d^2\zeta/dz^2)(dRi_g/d\zeta)}{(d\zeta/dz)^2(d^2Ri_g/d\zeta^2)}\right|.
\]
If \(E_{\text{omit}} < 0.05\), use simpler \(\frac{d^2Ri_g}{dz^2} \approx \frac{1}{L^2} \frac{d^2Ri_g}{d\zeta^2}\).

---

## 9. Curvature-Aware Correction (Concept)

Problem: Coarse first-layer thickness \(\Delta z\) → \(Ri_b < Ri_g(z_g)\) (with \(z_g = \sqrt{z_0 z_1}\), geometric mean height) → overestimated turbulent mixing coefficients \(K_m,K_h\).

Solution strategy:
1. Evaluate analytic \(Ri_g\) & curvature at \(z_g\).
2. Estimate layer bias \(B = Ri_g(z_g)/Ri_b\).
3. Apply a damping factor to eddy diffusivities:
   \[
   K_{m,h}^{*} = K_{m,h} \, G(\Delta z, \Delta),\quad
   G < 1 \text{ when curvature magnitude and }\Delta z \text{ are large}.
   \]
4. Preserve neutral curvature \(2\Delta\) (do not alter near-neutral baseline).
5. Optionally embed grid spacing in \(\phi\) via multiplicative tail modifier \(f_c(\zeta,\Delta z) = \exp[-D (\zeta/\zeta_r)(\Delta z / \Delta z_r)]\) choosing exponents to keep \(V_{\log}(0)\) unchanged.

---

## 10. Practical Workflow (Stable Case)

| Step | Action | Output |
|------|--------|--------|
| 1 | Gather tower / LES: \(u,v,\theta,\overline{w'\theta'}, u_*\) | Profiles, flux |
| 2 | Compute (bulk or local) \(L\) and \(\zeta=z/L\) | \(\zeta\) array |
| 3 | Fit or choose \(\alpha_{m,h}, \beta_{m,h}\) (or quadratic surrogate) | Parameter set |
| 4 | Compute \(\Delta, c_1\); record \(2\Delta\) | Neutral curvature |
| 5 | Evaluate \(Ri_g(\zeta)\), curvature via formula | Reference analytic fields |
| 6 | Compare with bulk \(Ri_b\) (layer 0–1) → bias \(B\) | Bias metric |
| 7 | If coarse: apply correction \(G\) or tail modifier | Adjusted \(K_{m,h}\) |
| 8 | Map to height curvature if \(L(z)\) variable | \(d^2Ri_g/dz^2\) |
| 9 | Diagnostics: \(B, 2\Delta, E_{\text{omit}},\) inflection (if any) | QC / metadata |

---

## 11. Typical Stable Parameter Ranges and Coefficients

Approximate literature near-neutral ranges (site dependent):
- \(\alpha_{m,h} \approx 0.45\text{–}0.55\), \(\beta_{m,h} \approx 14\text{–}16\).
Derived (Q‑SBL) coefficients:
- \(a = \alpha\beta \approx 6.3\text{–}8.8\).
- \(b = 0.5\alpha(\alpha+1)\beta^2 \approx 64\text{–}110\).
Neutral curvature coefficient:
\[
2\Delta = 2(\alpha_h\beta_h - 2\alpha_m\beta_m)\ \text{(often negative in stable sets)}.
\]

---

## 12. Minimal Numerical Core (Pseudocode)

```python
def phi_power(zeta, alpha, beta):
    d = 1 - beta*zeta
    if d <= 1e-8: raise ValueError("ζ beyond domain")
    return d**(-alpha)

def rig_and_curv(zeta, am,bm,ah,bh):
    pm = phi_power(zeta, am, bm)
    ph = phi_power(zeta, ah, bh)
    F  = ph / (pm*pm)
    V  = (ah*bh)/(1 - bh*zeta) - 2*(am*bm)/(1 - bm*zeta)
    W  = (ah*bh*bh)/(1 - bh*zeta)**2 - 2*(am*bm*bm)/(1 - bm*zeta)**2
    Ri = zeta * F
    curv = F*(2*V + zeta*(V*V - W))
    return Ri, curv

# Neutral coefficients
Delta = ah*bh - 2*am*bm
c1    = ah*bh*bh - 2*am*bm*bm
```

Variable \(L(z)\) mapping (if required):
\[
\frac{d^2 Ri_g}{dz^2} \approx \frac{1}{L^2} \frac{d^2 Ri_g}{d\zeta^2} \quad\text{(if }E_{\text{omit}} \text{ small)}.
\]

---

## 13. Glossary (Condensed)

| Term | Meaning |
|------|---------|
| MOST | Monin–Obukhov Similarity Theory |
| \(L\) | Obukhov length (stability scaling) |
| \(\zeta=z/L\) | Dimensionless height |
| \(\phi_{m,h}\) | Dimensionless gradient (momentum / heat) |
| \(Ri_g\) | Gradient Richardson number (point stability) |
| \(Ri_b\) | Bulk Richardson number (layer difference) |
| \(\Delta\) | Neutral curvature half-coefficient (since curvature = \(2\Delta\)) |
| \(V_{\log}, W_{\log}\) | Log-derivative combination and its derivative |
| Q‑SBL | Quadratic stable-layer surrogate |
| \(E_{\text{omit}}\) | Metric for deciding if variable \(L\) effects can be ignored |

---

## 14. Core References (Foundational & Stable BL)

1. Monin & Obukhov (1954) – Original similarity framework.  
2. Businger et al. (1971), Dyer (1974) – Classical empirical functions.  
3. Högström (1988) – Review / stable fits.  
4. Beljaars & Holtslag (1991) – Stable extensions / surrogate forms.  
5. Cheng & Brutsaert (2005) – Monotonic, pole‑free formulations.  
6. McNider et al. (2012) – Resolution issues in stable boundary layers.

---

## 15. Next Step (Optional)

If desired, proceed to:  
“Coupling curvature-aware adjustment into the implicit vertical diffusion (TDMA / Crank–Nicolson) scheme” — boundary-condition handling and stability impacts.

---

## 16. Summary

For stable boundary layers on coarse grids, uncorrected bulk Richardson number underestimates true near-surface stability due to concave-down \(Ri_g(\zeta)\). Analytic curvature (through \(\Delta\) and higher terms) enables bias-aware eddy diffusivity adjustment while preserving neutral physics. Careful parameter selection (SBL-focused), optional quadratic surrogates, and variable-\(L\) diagnostics reduce grid sensitivity and improve physical fidelity of nocturnal cooling and turbulence suppression.

<!-- End revised introductory overview -->
