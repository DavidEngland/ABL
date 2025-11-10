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
