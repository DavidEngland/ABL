# Gradient and Bulk Richardson Numbers: ζ Mapping and Curvature Basis

## 1. Definitions
Gradient:
$$Ri_g=\frac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2}.$$
Bulk (layer z₀→z₁):
$$Ri_b=\frac{g}{\theta}\frac{\Delta\theta\,\Delta z}{(\Delta U)^2}.$$

## 2. MOST Relation
$$Ri_g(\zeta)=\zeta\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=\zeta F(\zeta),\qquad \zeta=z/L.$$

## 3. Near-Neutral Series
Let
$$\phi_{m,h}=1+a_{m,h}\zeta+b_{m,h}\zeta^2+O(\zeta^3).$$
Then
$$Ri_g=\zeta + \Delta\zeta^2 + \tfrac12(\Delta^2+c_1)\zeta^3+O(\zeta^4),$$
$$\Delta=a_h-2a_m,\quad c_1=b_h-2b_m.$$

## 4. Inversion ζ(Ri)
$$\zeta=Ri_g - \Delta Ri_g^2 + \Big(\tfrac32\Delta^2 - \tfrac12 c_1\Big)Ri_g^3+O(Ri_g^4).$$
Seed for Newton refinement when evaluating φ at given Ri.

## 5. Curvature
Log derivatives:
$$V_{\log}=\frac{\phi_h'}{\phi_h}-2\frac{\phi_m'}{\phi_m},\quad W_{\log}=V_{\log}'.$$
Curvature:
$$\frac{d^2Ri_g}{d\zeta^2}=F\big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\big],\quad F=\frac{\phi_h}{\phi_m^2}.$$
Neutral:
$$\left.\frac{d^2Ri_g}{d\zeta^2}\right|_0=2\Delta.$$

## 6. Bulk vs Point Bias
Concave-down ($\Delta<0$) ⇒
$$Ri_b < Ri_g(z_g),\qquad z_g=\sqrt{z_0 z_1},\quad B=\frac{Ri_g(z_g)}{Ri_b}>1.$$

## 7. Correction Principle
Grid damping $G(\zeta,\Delta z)$ with $G(0,\Delta z)=1$, $\partial_\zeta G|_0=0$ preserves 2Δ while reducing $B$ at coarse Δz.

## 8. Critical Richardson Number
Fixed $Ri_c$ vs dynamic $Ri_c^*$ informed by curvature growth (magnitude of ζ(V_{\log}^2-W_{\log})) and inversion strength.

## 9. Key Identities
$$Pr_t=\frac{\phi_h}{\phi_m},\quad f_m(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))^2},\quad f_h(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\phi_h(\zeta(Ri_g))}.$$

## 10. Generic Scalar Closure

For any scalar $q$:
\[
f_q(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\,\phi_q(\zeta(Ri_g))},\quad K_q=l_m^2 S f_q.
\]
Schmidt number:
\[
Sc_t^{(q)}=\frac{\phi_q}{\phi_m}.
\]
Series (near-neutral):
\[
f_q\approx 1+a_q Ri_g+(b_q-a_q\Delta+2a_m a_q)Ri_g^2.
\]

Concise algorithm (near-neutral):
1. Compute Ri_g.
2. ζ ≈ Ri_g − Δ Ri_g² (cubic if needed).
3. Evaluate φ_m, φ_h; obtain K_m,K_h.

## 11. Numerical Estimation of Ri_g and Ri_b

Given discrete $z_k, U_k, \theta_k$:

**Point gradient (centered):**
$$
\partial U/\partial z \Big|_{z_k} \approx (U_{k+1} - U_{k-1})/(z_{k+1} - z_{k-1}).
$$

**Bulk Ri_b (layer $[z_0,z_1]$):**
- Definition: $Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz$.
- Trapezoid: $Ri_b \approx \frac{1}{2}[Ri_g(z_0) + Ri_g(z_1)]$.
- Simpson (3-pt): $Ri_b \approx \frac{1}{6}[Ri_g(z_0) + 4Ri_g(z_g) + Ri_g(z_1)]$.

**Representative heights:**
- $z_g = \sqrt{z_0 z_1}$ (geometric mean, midpoint in $\ln z$).
- $z_L = (z_1 - z_0)/\ln(z_1/z_0)$ (logarithmic mean, exact for $\Delta U$ in log wind).
- $z_a = (z_0 + z_1)/2$ (arithmetic mean, biases high for log profiles).

Use $z_g$ for point $Ri_g$ evaluation; use $z_L$ for exact layer-averaged gradient matching.

## 12. Practical Estimation Techniques & Jensen reminder
- Representative heights:
  - z_g = √(z0 z1) (geometric mean) → use for evaluating Ri_g for log/power-law profiles.
  - z_L = (z1-z0)/ln(z1/z0) (log mean) → use when matching ΔU exactly.
- Finite-difference estimates:
  - Centered (interior): (U_{k+1}-U_{k-1})/(z_{k+1}-z_{k-1})
  - First-layer forward difference: (U1-U0)/(z1-z0) with z_rep = z_g or z_L
- Bulk vs gradient correction workflow:
  1. Compute Ri_g(z_g).
  2. Compute Ri_b (bulk formula or integral).
  3. Compute B = Ri_g(z_g)/Ri_b. If B>1.1, flag for curvature-aware correction.
- Numerical integration: prefer Simpson/trapezoid on Ri_g(z) over bulk formula for curved profiles.

## 13. Quick decision tree
- If B ≤ 1.05 → no correction.
- If 1.05 < B ≤ 1.3 → mild correction (K multiplier with small γ).
- If B > 1.3 and strong inversion (Γ large) → apply mixing-length reduction + K damping; consider raising Ri_c*.

## 14. Mixed Concavity and Inflection Handling
If \(d^2Ri_g/d\zeta^2\) changes sign inside a layer (inflection at \(\zeta_{\text{inf}}\)):
- Split layer at \(z_{\text{inf}}=\zeta_{\text{inf}}L\).
- Apply bias logic (Ri_b < Ri_g midpoint) only to concave-down segment.
- Concave-up segment may yield \(Ri_b > Ri_g\) locally.
Report:
\[
Ri_b = \frac{(z_{\text{inf}}-z_0)Ri_{b1} + (z_1-z_{\text{inf}})Ri_{b2}}{\Delta z},\quad
B_1=\frac{Ri_g(z_{g1})}{Ri_{b1}},\ B_2=\frac{Ri_g(z_{g2})}{Ri_{b2}}.
\]
Use damping only below \(\zeta_{\text{inf}}\).

## X. φ‑Agnostic Diagnosis & Correction (practical cookbook)

When φ(ζ) is unknown inside a model, use these robust, model‑agnostic steps:

1) Representative heights
- z_g = sqrt(z0*z1)  (use for Ri_g point evaluation)
- z_L = (z1 - z0)/ln(z1/z0)  (use if you need exact ΔU reconstruction)

2) Finite-difference estimators (use available levels)
- centered shear at z_k:
  U_z = (U_{k+1} - U_{k-1})/(z_{k+1} - z_{k-1})
  theta_z = (TH_{k+1} - TH_{k-1})/(z_{k+1} - z_{k-1})
- point Ri_g:
  Ri_g = (g/theta_k) * theta_z / (U_z**2)

3) Bulk Ri_b (two-level):
- Ri_b = (g/theta_ref) * (TH1 - TH0) * (z1 - z0) / ((U1 - U0)**2)

4) Bias check and correction (spreadsheet formula friendly)
- B = Ri_g(z_g) / Ri_b
- If B ≤ 1.05 → no change
- Else apply K modifier:
  - G = EXP( -D*(Δz/Δz_ref)^p * (ζ/ζ_ref)^q )
  - K_new = K_old * G

5) Safe default surrogate (if you must produce f_m(Ri) or f_h(Ri))
- Exponential Ri closure (pole-free):
  f_m(Ri) = exp( -γ_m * Ri / Ri_c* )
  f_h(Ri) = exp( -γ_h * Ri / Ri_c* )
  suggested γ_m ≈ 1.8, γ_h ≈ 1.5; Ri_c* dynamic or 0.25 default

6) Minimal pseudocode
```python
# φ-agnostic correction pseudocode
z_g = sqrt(z0*z1)
Ri_g_zg = compute_point_Ri(z_g)  # centered diffs or interpolation
Ri_b = compute_bulk_Ri(z0,z1)
B = Ri_g_zg / Ri_b
if B > 1.1:
    G = exp(-D*(dz/10.0)**p * (zeta/zeta_ref)**q)
    K_star = K * G
else:
    K_star = K

# fallback Ri-based closure
f_m = exp(-gamma_m * Ri / Ric_star)
```

Notes
- Always verify neutral‑curvature preservation by testing near ζ→0 that your G does not change first derivative (G′(0)=0).
- For Excel: use LOG(), SQRT(), EXP() and user‑defined constants in header cells to allow easy tuning.

## 13. Mixed Concavity Handling
If \(\tfrac{d^2Ri_g}{d\zeta^2}\) changes sign:
- Split layer at inflection \(\zeta_{\text{inf}}\).
- Apply bias logic separately to concave-down and concave-up segments.
- Recombine weighted averages.

---

## 14. φ‑Agnostic Surrogate
When φ-functions are unknown:
- Use exponential Ri closures:
  \[
  f_m(Ri) = e^{-\gamma_m Ri / Ri_c^*}, \quad f_h(Ri) = e^{-\gamma_h Ri / Ri_c^*}
  \]
- Suggested: \(\gamma_m \approx 1.8, \gamma_h \approx 1.5\).
- \(Ri_c^*\) dynamic or default 0.25.

---

## 15. Minimal Pseudocode
```python
z_g = sqrt(z0*z1)
Ri_g_zg = compute_point_Ri(z_g)
Ri_b = compute_bulk_Ri(z0,z1)
B = Ri_g_zg / Ri_b

if B > 1.1:
    G = exp(-D*(dz/dz_ref)**p * (zeta/zeta_ref)**q)
    K_star = K * G
else:
    K_star = K

# φ-agnostic fallback
f_m = exp(-gamma_m * Ri / Ric_star)
```

---

# Curvature in ζ versus z

Curvature in ζ is natural in MOST because the similarity functions are defined in terms of the non-dimensional height \(\zeta=z/L\). If \(L\) is treated as locally constant across a thin layer, then derivatives transform simply: \(\frac{d}{dz}=\frac{1}{L}\frac{d}{d\zeta}\) and \(\frac{d^2}{dz^2}=\frac{1}{L^2}\frac{d^2}{d\zeta^2}\). In that case, the sign and relative magnitude of curvature are preserved between \(z\) and \(\zeta\). When \(L\) varies with height, curvature in \(z\) picks up extra terms involving \(dL/dz\), so \(\zeta\)-based diagnostics cleanly separate profile shape from coordinate effects.

---

## Near-neutral coefficients from Businger–Dyer

For the classic Businger–Dyer (BD) stable formulations near \(\zeta\to 0^+\):
\[
\phi_m = 1 + a_m \zeta,\quad \phi_h = 1 + a_h \zeta,
\]
with commonly used values \(a_m\approx 4.7\) and \(a_h\approx 7.8\). If you retain only the linear terms (i.e., \(b_m=b_h=0\)), then
\[
\Delta = a_h - 2a_m = 7.8 - 9.4 = -1.6,\qquad c_1 = b_h - 2b_m = 0.
\]
Implications:

- **Neutral curvature:** \(\left.\frac{d^2 Ri_g}{d\zeta^2}\right|_0 = 2\Delta = -3.2\) → concave-down.
- **Bias direction:** concave-down implies \(Ri_b < Ri_g(z_g)\) and a bulk-versus-point bias factor \(B>1\).
- **ζ inversion (cubic):**
  \[
  \zeta \approx Ri_g - \Delta Ri_g^2 + \Big(\tfrac32\Delta^2 - \tfrac12 c_1\Big)Ri_g^3
  = Ri_g + 1.6\,Ri_g^2 + 3.84\,Ri_g^3.
  \]

For unstable BD, \(\phi\) functions are non-linear (e.g., \(\phi_m=(1-16\zeta)^{-1/4}\), \(\phi_h=(1-16\zeta)^{-1/2}\)). A Taylor expansion about \(\zeta=0^-\) yields finite linear coefficients as well:
\[
\phi_m \approx 1 + 4\zeta + 10\zeta^2 + \dots,\quad
\phi_h \approx 1 + 8\zeta + 48\zeta^2 + \dots,
\]
so near-neutral on the unstable side you’d have \(a_m\approx 4\), \(a_h\approx 8\), giving \(\Delta\approx 0\) (specifically \(\Delta = 8 - 2\cdot 4 = 0\)), i.e., weak curvature in the immediate neutral limit and rapidly increasing nonlinearity at larger \(|\zeta|\). If you prefer other empirical constants (e.g., 5 and 5 for modified BD), update \(a_m,a_h\) and recompute \(\Delta,c_1\) accordingly.

---

## Special case — identical linear φ (φ_h = φ_m = 1 + β ζ)

Assume
\[
\phi_h(\zeta)=\phi_m(\zeta)=1+\beta\zeta,\qquad \beta=4.7.
\]

Then
- F = φ_h/φ_m^2 = 1/(1+\beta\zeta).
- Closed form:
  \[
  Ri_g(\zeta)=\frac{\zeta}{1+\beta\zeta}.
  \]
- Near‑neutral Taylor series:
  \[
  Ri_g(\zeta)=\zeta - \beta\zeta^2 + \beta^2\zeta^3 + O(\zeta^4).
  \]
  Hence the near‑neutral coefficients are a_m=a_h=β and
  \[
  \Delta = a_h - 2a_m = -\beta = -4.7,\qquad 2\Delta=-9.4.
  \]
- ζ(Ri) inversion (exact rational form and series):
  \[
  \zeta=\frac{Ri}{1-\beta Ri} = Ri + \beta Ri^2 + \beta^2 Ri^3 + O(Ri^4).
  \]
- Turbulent Prandtl: Pr_t = φ_h/φ_m = 1 (unit Prandtl in this special case).

Implications (one line)
- This symmetric linear choice yields strong negative neutral curvature (2Δ = −9.4), so Jensen bias is significant for coarse Δz; apply curvature‑aware damping G(ζ,Δz) or the Q‑SBL surrogate as described in the main text.

## Curvature mapping: ζ to z

- **If L is uniform in the layer:**
  - **Label:** Derivative scaling
  - \[
    \frac{d^2 Ri_g}{dz^2} = \frac{1}{L^2}\frac{d^2 Ri_g}{d\zeta^2}
    \]
  - **Result:** Same concavity and bias logic; only magnitude rescales by \(1/L^2\).

- **If L varies with z:**
  - **Label:** Extra terms
  - \[
    \frac{d^2 Ri_g}{dz^2} = \frac{1}{L^2}\frac{d^2 Ri_g}{d\zeta^2}
    - 2\frac{1}{L^3}\frac{dL}{dz}\frac{d Ri_g}{d\zeta}
    - \text{terms with } \frac{d^2 L}{dz^2}
    \]
  - **Result:** Curvature in \(z\) combines profile shape and stability variation; using \(\zeta\) isolates the MOST shape. Practically, use \(\zeta\)-curvature for diagnostics and treat \(L(z)\) variability via layer splitting or effective \(L\).

---

## Practical guidance

- **Near-neutral diagnostics:** Use BD-derived \(a_m,a_h\) to compute \(\Delta\) and neutral curvature \(2\Delta\). This sets the expected sign of the bulk vs point bias.
- **Bias correction:** Apply damping \(G(\zeta,\Delta z)\) designed so \(G(0)=1\) and \(G'(0)=0\), preserving neutral curvature while reducing coarse-grid bias.
- **Representative height:** Evaluate point \(Ri_g\) at \(z_g=\sqrt{z_0z_1}\) and use Simpson/trapezoid integration for \(Ri_b\) when curvature is non-negligible.
- **ζ inversion for closures:** Use \(\zeta(Ri_g)\) from the series as a seed for Newton to evaluate \(\phi_m,\phi_h\) robustly when Ri is the control variable.

---

## Quick check with BD-stable

- **Coefficients:** \(a_m=4.7,\ a_h=7.8 \Rightarrow \Delta=-1.6\).
- **Neutral curvature:** \(2\Delta=-3.2\) → concave-down; expect \(B>1\).
- **ζ seed:** \(\zeta \approx Ri_g + 1.6\,Ri_g^2 + 3.84\,Ri_g^3\).
- **Action:** Prefer Simpson for \(Ri_b\); if \(B>1.05\), apply mild damping; if \(B>1.3\) with strong inversion, reduce mixing length and consider raising \(Ri_c^*\).

If you’re using a different BD constant set, share them and I’ll plug them in to give you the updated \(\Delta\), curvature, and ζ-inversion coefficients.

### 6A. Log‑mean (z_L) sandwiched between geometric and arithmetic means

Recall the standard ordering for positive heights z0<z1:
- geometric mean: z_g = sqrt(z0*z1)
- logarithmic mean: z_L = (z1 - z0) / ln(z1/z0)
- arithmetic mean: z_a = (z0 + z1)/2

These satisfy
z_g ≤ z_L ≤ z_a (strict inequalities when z0 ≠ z1).
Because Ri_g typically decreases with z in the stable surface layer, this implies the pointwise ordering
Ri_g(z_g) ≥ Ri_g(z_L) ≥ Ri_g(z_a).

Combining Jensen's inequality (for concave Ri_g)
Ri_b = (1/Δz) ∫_{z0}^{z1} Ri_g(z) dz ≤ Ri_g(z_a),
we obtain the useful sandwich
Ri_b ≤ Ri_g(z_a) ≤ Ri_g(z_L) ≤ Ri_g(z_g).

Practical takeaway:
- If U(z) is log-like, the log-mean z_L is the exact height that makes a point gradient reproduce ΔU; using z_L for the shear estimate reduces denominator bias in Ri_b.
- If the profile is closer to a power-law/log shape, z_g remains a good representative for evaluating Ri_g point values (e.g., curvature diagnostics).
- For robust Ri estimation:
  1. Use z_L to compute ΔU (or use log-mean when reconstructing shear from a log wind law).
  2. Evaluate Ri_g at z_g (or z_L) depending on which quantity you wish to represent (point gradient vs reproducing ΔU exactly).
  3. Compute Ri_b using Simpson/trapezoid on Ri_g(z) when a fine profile is available; otherwise compute Ri_b from observed Δθ and ΔU but use z_L for ΔU reconstruction when the wind is near-log.

## Log‑mean justification & practical Ri_b estimation

Summary (quick)
- Geometric mean z_g = √(z0 z1) is the midpoint in ln z (useful for evaluating point Ri_g when profiles are log/power‑law).
- Logarithmic mean z_L = (z1−z0)/ln(z1/z0) is the unique height where a point gradient reproduces the exact layer ΔU for a log wind law; it sits between z_g and arithmetic mean z_a:
  z_g ≤ z_L ≤ z_a.
- Use z_L when your denominator (shear) should match a log law (exact ΔU). Use z_g when evaluating pointwise Ri_g or curvature diagnostics.

Short justification of ordering
- ln is concave on (0,∞). By Jensen for the integral average over [z0,z1],
  (1/Δz) ∫_{z0}^{z1} ln z dz ≤ ln( (1/Δz) ∫_{z0}^{z1} z dz ) = ln z_a,
  so ln z_L ≤ ln z_a ⇒ z_L ≤ z_a.
- The geometric mean satisfies ln z_g = ½(ln z0 + ln z1). For smooth log-like profiles and thin layers the integral mean (ln z_L) lies between the midpoint-in-log (ln z_g) and the arithmetic log average; standard mean inequality gives z_g ≤ z_L. (Reference: inequality chain for classical means: geometric ≤ logarithmic ≤ arithmetic.)

When to use which mean (practical rule)
- If U(z) ≈ (u*/κ) ln(z/z0)+C (log law) then compute shear from z_L — it reproduces ΔU exactly.
- If you need a representative height for MOST φ evaluation or curvature diagnostics, use z_g — it's the midpoint in ln z and minimizes log-space quadrature error.
- If layer is very thin (z1/z0 → 1) all means converge and choice is immaterial.

Practical recipes to estimate bulk Ri near the surface

Notation: layer [z0,z1], Δz = z1−z0, surface fluxes u_*, θ_* (if available), mean θ̄ representative.

A. Direct (observations / profiles)
- ΔU = U1 − U0; Δθ = θ1 − θ0
- Ri_b (direct) = (g/θ̄) · (Δθ · Δz) / (ΔU)^2
  - Use centered differences for interior layers when data at 3+ levels are available.

B. Log‑law reconstruction from surface fluxes (use when only surface fluxes available)
- Assume log wind & temperature increments:
  ΔU_log = (u*/κ) ln(z1/z0)
  Δθ_log = (θ*/κ) ln(z1/z0)   (use carefully: θ* = −w'θ'/u*; sign convention)
- Then
  Ri_b (from surface fluxes) = (g/θ̄) · [ (θ*/κ) ln(z1/z0) · Δz ] / [ (u*/κ ln(z1/z0))^2 ]
  which simplifies to
  Ri_b = (g/θ̄) · (θ* · Δz) / [ u*^2 · ln(z1/z0) ].
- Use absolute/consistent sign for θ* (stable nights θ* typically >0 when using the definition above).

C. Use z_L for shear when profile is log‑like (exact matching)
- Compute ΔU from a point shear at z_L:
  (∂U/∂z)|_{z_L} = (u*/κ z_L)  ⇒ ΔU ≈ (∂U/∂z)|_{z_L} · Δz  (for small Δz)
  or directly use ΔU_log = (u*/κ) ln(z1/z0) which is exact for log-law.

D. Simpson / Trapezoid on Ri_g when Ri_g(z) is known (preferred if full profile available)
- Compute Ri_g(z) at high resolution inside layer (from gradients or analytic form).
- Trapezoid: Ri_b ≈ 0.5 [Ri_g(z0) + Ri_g(z1)].
- Simpson (3‑point using z_g): Ri_b ≈ (1/6)[Ri_g(z0) + 4 Ri_g(z_g) + Ri_g(z1)].
- These numeric integrals are more accurate than the Δθ/ΔU bulk formula when curvature is important.

E. Bias comparator & decision metric
- Compute Ri_g at z_g (point evaluator) and Ri_b (from any of A–D).
- Bias ratio B = Ri_g(z_g) / Ri_b.
  - If B ≈ 1 → negligible curvature bias.
  - If B > 1.05 → curvature-induced underestimation; consider correction (see curvature-aware modifiers).
  - If B < 0.95 → layer concave-up behavior; treat accordingly or split layer.

F. φ‑agnostic fallback (no φ knowledge, minimal data)
- If only surface fluxes and one upper-level measurement available:
  1. Compute ΔU from log law (z_L or ln formula) using u_*.
  2. Compute Δθ from θ* if available, otherwise use measured Δθ.
  3. Form Ri_b with the expression in B.
- Document uncertainty: use ensemble of z_rep choices (z_g, z_L) to bound Ri_b.

Implementation tips & caveats
- Use z_L when denominator (shear estimate) must be consistent with log-law surface schemes (e.g., when u_* is provided by surface layer parameterization).
- Use z_g for MOST‑based φ evaluations and Jensen bias checks — it aligns with log-space midpoint and curvature diagnostics.
- Be careful with roughness heights z0, z0_h and displacement d; use displaced geometrical/log means when d ≠ 0 (z replaced by z−d).
- For noisy tower data, smooth lightly before differentiation; prefer Simpson/trapezoid on Ri_g if fine resolution data exist.
- Preserve sign conventions for θ* and Ri (stable Ri>0 using consistent θ* sign).

Short worked algebraic reminder (for copy/paste)
- If u_*, θ* known and log law assumed:
  ΔU = (u*/κ) ln(z1/z0)
  Δθ = (θ*/κ) ln(z1/z0)
  Ri_b = (g/θ̄) · [ (θ*/κ) ln(z1/z0) · Δz ] / [ (u*/κ ln(z1/z0))^2 ]
       = (g/θ̄) · (θ* Δz) / [ u*^2 ln(z1/z0) ].