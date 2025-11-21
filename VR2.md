
$$
f_s(Ri) = \exp\left(- \frac{\gamma Ri}{Ri_c}\right)

$$
f_{c}(\Delta z) \;=\; e^{D\left(\frac{\gamma}{Ri_c}\right)\;Ri\left(1-\frac{\Delta z_r}{\Delta z}\right)}
$$

$$
f_s*f_c = e^{-\frac{\gamma}{Ri_c}Ri} * e^{D\left(\frac{\gamma}{Ri_c}\right)Ri\left(1 - \frac{\Delta z_r}{\Delta z}\right)}
$$

$$
f_{is}=f_s*f_c = e^{-\frac{\gamma}{Ri_c}Ri\left[\frac{(1-D)\Delta z + D \Delta z_r}{\Delta z} \right]}
$$

### Reconciliation: McNider–Biazar D versus Bias B

Original bulk-driven ODE (paper draft form):
\[
\frac{d\ln f_c}{d\ln \Delta z} = - D \left(\frac{Ri_b}{Ri_{\text{ref}}}\right)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q
\]
with layer bulk Richardson number \(Ri_b\), reference \(Ri_{\text{ref}}\) (e.g. 0.25), exponent \(q\ge1\).

Bias-based form (current framework):
\[
\frac{d\ln f_c}{d\ln \Delta z} = - \alpha (B-1)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q,\qquad
B=\frac{Ri_g(z_g)}{Ri_b},\ z_g=\sqrt{z_0 z_1}.
\]

Parameter mapping (equate RHS):
\[
\alpha (B-1) = D \frac{Ri_b}{Ri_{\text{ref}}} \;\Rightarrow\;
D = \alpha (B-1)\frac{Ri_{\text{ref}}}{Ri_b}, \qquad
B-1 = \frac{D}{\alpha}\frac{Ri_b}{Ri_{\text{ref}}}.
\]

Implication:
Since typically \(Ri_b < Ri_g(z_g)\) in concave-down stable layers (\(B>1\)), using \(Ri_b\) directly underestimates needed correction amplitude by roughly a factor \(1/B\). Driving the ODE with \((B-1)\) captures curvature-induced bulk vs point deficit explicitly.

Numeric example:
Layer: \(Ri_g(z_g)=0.60,\ Ri_b=0.30 \Rightarrow B=2.0,\ B-1=1.0\).
Take \(\alpha=1.0,\ Ri_{\text{ref}}=0.25\).
Then \(D = 1.0 \cdot 1.0 \cdot (0.25/0.30) \approx 0.83\).
If manuscript reports \(D=1.0\), equivalent \(\alpha \approx 1.2\) by inversion.

Recommended operational ODE:
\[
\frac{d\ln f_c}{d\ln \Delta z} = - \alpha (B-1)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q,\quad
f_c(\Delta z,\zeta)=\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}.
\]

Practical workflow:
1. Compute \(Ri_b\) (bulk) and \(Ri_g(z_g)\) (gradient at geometric mean).
2. Form \(B\). If \(B \le B_{\text{thr}}\) (≈1.05) set \(f_c=1\).
3. Else apply exponential or power-law form with driver \((B-1)\).

Code snippet (bulk-to-bias conversion):
```python
# Given D (bulk form) from paper and measured Ri_b, Ri_ref, alpha target:
# Recover alpha or translate D to bias-based alpha
alpha_equiv = D * (Ri_b / Ri_ref) / (B - 1.0)

# Preferred direct usage:
exponent = -alpha * (B - 1.0) * (dz/dz_ref)**p * (zeta/zeta_ref)**q
f_c = math.exp(exponent)
```

Summary:
- Bulk-driven D hides curvature strength when \(Ri_b\) is suppressed.
- Bias form with \(B-1\) is dimensionless, directly tied to Jensen concavity effect, and transfers across regimes.
- Report both \(B\) and chosen \(\alpha\) for reproducibility; provide conversion via \(D = \alpha (B-1)(Ri_{\text{ref}}/Ri_b)\).

## Improved Mathematical Reconstruction and Critique of the McNider–Biazar D Formulation

### 1. Original handwritten form (as entered)
\[
f_c(\Delta z) = e^{D\left(\frac{\gamma}{Ri_c}\right) Ri \left(1 - \frac{\Delta z_r}{\Delta z}\right)}
\]
Issues:
- Sign: positive exponent increases mixing for larger Ri; for stability damping a negative sign is normally required.
- Units: \((\gamma/Ri_c)Ri\) is dimensionless only if \(\gamma\) itself is dimensionless (OK), but \(D\) must be dimensionless too; confirm.
- Grid factor: \(1 - \Delta z_r/\Delta z\) → 0 at fine resolution (good) but saturates (→1) at coarse resolution; lacks ζ (stability-depth coupling).
- Driver: uses Ri (unspecified point or bulk); if Ri_b is used, curvature bias is embedded and weakens response relative to Ri_g.

### 2. Corrected “linear tail-extension” form (keeping their structure)
Damping (physically consistent) version:
\[
f_c^{(lin)}(\Delta z) = \exp\left[-\,D \left(\frac{\gamma}{Ri_c}\right) Ri \left(1 - \frac{\Delta z_r}{\Delta z}\right)\right].
\]
Properties:
- \(f_c^{(lin)} \to 1\) as \(\Delta z \to \Delta z_r\).
- Monotone decrease with \(\Delta z\) for \(D,\gamma,Ri>0\).
Limitations:
- No neutral-preserving slope guarantee (\(\partial_\zeta f_c|_{\zeta=0} \neq 0\) unless Ri→0).
- No curvature awareness (Ri'' absent).
- No explicit bias correction; purely grid-scalar.

### 3. Bias / ODE-driven derivation (recommended)
Starting local invariance condition (grid independence of effective stability):
\[
\frac{d\ln f_c}{d\ln \Delta z} = - \alpha (B - 1)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q.
\]
Integrates to (power law):
\[
f_c^{(pl)}(\Delta z,\zeta) = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}.
\]
Exponential surrogate (smoother tail):
\[
f_c^{(exp)}(\Delta z,\zeta) = \exp\left[-\alpha (B-1)\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q\right].
\]
Neutral preservation: choose \(q\ge 2\) so \(\partial_\zeta f_c|_0=0\). Grid convergence: \(f_c\to1\) as \(\Delta z\to0\).

### 4. Curvature-aware extension for D
Observed: D constant fails when curvature varies strongly with height. Define:
\[
D_{\text{eff}}(z) = D_0 + M \frac{|Ri''(z)|}{|Ri''(z)| + C},
\]
with:
- \(Ri''(z)\): second derivative (height or ζ-based, consistent scaling).
- \(C\): curvature soft-threshold (prevents runaway at small |Ri''|).
- Bound: \(D_{\min} \le D_{\text{eff}} \le D_{\max}\) (e.g. 0.2 ≤ D_eff ≤ 0.7).
Use \(Ri''\) from analytic MOST curvature if φ known; else numerical second difference with smoothing.

### 5. Reconciliation (verification)
Earlier mapping:
\[
\alpha (B-1) = D \frac{Ri_b}{Ri_{\text{ref}}} \Rightarrow D = \alpha (B-1)\frac{Ri_{\text{ref}}}{Ri_b}.
\]
Still valid if their ODE used \(Ri_b/Ri_{\text{ref}}\) as driver. If they reported D from bulk-Ri fits, recover \(\alpha\) via:
\[
\alpha = D \frac{Ri_b}{Ri_{\text{ref}}(B-1)}.
\]
Caution: use time/height-averaged Ri_b and Ri_g(z_g) over representative stable intervals; instant noisy values distort mapping.

### 6. Critique Summary
| Aspect | Handwritten form | Recommended |
|--------|------------------|-------------|
| Sign | Growth (mixing increase) | Damping (negative exponent) |
| Driver | Ri (ambiguous; likely bulk) | (B - 1) explicit bias |
| Neutral invariance | Not guaranteed | Enforced (q ≥ 2) |
| Curvature use | None | Optional D_eff(Ri'') |
| Grid/ζ scaling | Linear (1 - Δz_r/Δz) only | Power/exponential with ζ coupling |
| Calibration clarity | Weak (single D) | Multi-parameter (α, p, q, bounds) |

### 7. Path Forward (detailed steps)
1. Archive current formulations (keep original fc text for provenance).
2. Implement bias-based exponential form \(f_c^{(exp)}\) with initial parameters: \(\alpha=1.0, p=1, q=2, \Delta z_{\text{ref}}=10\text{ m}, \zeta_{\text{ref}}=0.5\).
3. Compute per-layer B = Ri_g(z_g)/Ri_b; apply fc only if B > B_thresh (e.g. 1.05).
4. Add curvature diagnostic: compute Ri'' analytically (if φ known) or numerically; derive D_eff; replace α by α·(D_eff/D_ref) if needed.
5. Bound fc: fc = max(fc, fc_min) with fc_min ≈ 0.2.
6. Validation loop:
   - Metrics: bias reduction (ΔB), neutral curvature error (<5%), temperature profile RMSE vs fine grid, inversion height bias.
   - Sensitivity: vary α, p, q, and curvature parameters (D_0,M,C).
7. Reconcile historical D values: produce table converting published D to α for key cases.
8. Reporting: log (Ri_b, Ri_g, B, Ri'', D_eff, fc) each timestep for evaluation.
9. Optional refinement: switch from Ri_b to Ri_g(z_g) in curvature zone to avoid bulk underestimation bias feeding back into D.
10. Transition plan: after stable performance, test dynamic Ri_c*; ensure fc and Ri_c* changes do not conflict (log both).

### 8. Practical Implementation Snippet (concise)
```python
# inputs: z0,z1,Ri_b,Ri_g_zg,zeta,alpha,p,q,dz,dz_ref,zeta_ref,B_thresh,fc_min
B = Ri_g_zg / Ri_b if Ri_b > 1e-6 else 1.0
if B <= B_thresh:
    fc = 1.0
else:
    # optional curvature weighting
    D_eff = D0 + M * (abs(Ri_ddz)/(abs(Ri_ddz)+C))  # Ri_ddz precomputed
    alpha_eff = alpha * (D_eff / D_ref)
    exponent = -alpha_eff * (B - 1) * (dz/dz_ref)**p * (zeta/zeta_ref)**q
    fc = max(math.exp(exponent), fc_min)
K_new = K_old * fc
```

### 9. Where We Can Help (targeted)
- Supply analytic Ri'' for chosen φ set to avoid noisy numerical second differences.
- Provide automated blending of bias-based fc with curvature-weighted D_eff.
- Build validation scripts: bias statistics, curvature preservation, profile RMSE.
- Convert legacy D values to α and populate a reconciliation table for manuscript appendix.
- Advise on dynamic Ri_c* integration (sequence: stabilize fc first, then Ri_c*).

### 10. Key Risks and Mitigations
| Risk | Mitigation |
|------|------------|
| Noisy Ri'' inflates D_eff | Apply small centered smoothing window; cap derivative |
| Over-damping (fc too small) | Enforce fc_min; monitor flux sensitivity |
| Misuse of bulk Ri for driver | Replace with (B - 1) bias; retain mapping for comparison only |
| Neutral distortion | Verify \(\partial_\zeta fc|_0 \approx 0\) numerically each run |

### 11. Final Summary (one paragraph)
The original D-based correction as entered lacks a damping sign, explicit curvature coupling, and a bias-normalized driver; adopting a bias-triggered, neutral-preserving exponential (or power-law) form with optional curvature-weighted D_eff provides a controlled, dimensionless, and verifiable path that reduces grid-induced stability underestimation while safeguarding near-neutral invariants (2Δ). Mapping between historical D and current α allows continuity; curvature diagnostics refine spatial targeting; a structured validation loop consolidates performance before dynamic Ri_c* extension.
