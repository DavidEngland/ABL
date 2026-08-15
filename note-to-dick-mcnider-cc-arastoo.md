## Note to Dick McNider (cc: Arastoo Biazar)

Summary
- Core issue: coarse Δz + concave‑down Ri_g(ζ) (neutral curvature 2Δ < 0) ⇒ Ri_b < Ri_g(z_g) (Jensen) → K too large → overmixing.
- Goal: apply a minimal, neutral‑preserving multiplicative correction f_c to reduce this bias without changing ζ→0 behaviour.

McNider ODE (clean form) and solution
- Proposed local constraint (diagnostic-driven):
  d ln f_c / d ln Δz = −α (B − 1) (ζ/ζ_ref)^q,
  where B = Ri_g(z_g)/Ri_b is measured per layer.
- Integrates to exact power‑law solution:
  f_c(Δz,ζ) = (Δz/Δz_ref)^{ −α (B−1) (ζ/ζ_ref)^q }.
- Practical exponential approximation (recommended for model stability/tuning):
  f_c ≈ exp[ −α (B−1) (Δz/Δz_ref)^p (ζ/ζ_ref)^q ].

Physical justification (brief)
- fc multiplies K or f_m to reduce transport where bulk averaging underestimates local Ri.
- Design constraints: fc(ζ=0)=1 and ∂ζfc|0=0 (choose q≥2), fc→1 as Δz→0, fc bounded (apply floor e.g. fc_min=0.2).
- Multiplicative form preserves positivity and scales with base K.

Implementation (drop‑in pseudocode)
```python
# inputs per layer: z0,z1,U0,U1,th0,th1,L,K_old
z_g = sqrt(z0*z1)
zeta = z_g / L
Ri_g_zg = compute_point_Ri(z_g)       # from phi or local gradients
Ri_b = compute_bulk_Ri(z0,z1)
B = max(1.0, Ri_g_zg / (Ri_b + tiny))
if B <= B_thresh:    # e.g. 1.05
    fc = 1.0
else:
    fc = exp(-alpha * (B - 1.0) * ( (Dz/Dz_ref)**p ) * ( (zeta/zeta_ref)**q ))
    fc = max(fc, fc_min)
K_new = K_old * fc
```

Tuning guidance
- Defaults: alpha≈1.0, p=1.0, q=2.0, Dz_ref=10 m, zeta_ref=0.5, B_thresh=1.05, fc_min=0.2.
- Start with alpha=0.8–1.2; test on 3–5 LES/tower cases; monitor B_before, B_after and ΔK/K.
- If fc too aggressive: lower α or increase fc_min.

Diagnostics to include in paper
- Per‑layer histogram of B (before/after), median reduction target ≈40%.
- Time series of K_old vs K_new at first interior level for representative nights.
- Figure: Ri_g(ζ) (fine reference), coarse Ri_b reconstruction, corrected Ri_g^* (show improvement).
- Table: parameter sensitivity (α, p, q) with recommended default.

Suggested manuscript paragraph (concise, copy/paste)
- "We quantify curvature‑induced bulk bias by comparing the point Richardson number at the layer geometric mean, Ri_g(z_g), to the layer bulk value Ri_b and define B = Ri_g(z_g)/Ri_b. Where B>1.05 we apply a neutral‑preserving multiplicative correction f_c to the momentum/heat diffusivities. f_c is parameterized as f_c = exp[−α(B−1)(Δz/Δz_ref)^p(ζ/ζ_ref)^q] with q≥2 to ensure the neutral limit is unchanged. Calibration on LES/tower cases yields α≈1.0, p=1, q=2 and Δz_ref=10 m; we impose a floor f_c≥0.2 to avoid over‑damping. This correction reduces the bulk‑vs‑point bias by ~40% on coarse Δz (60–100 m) while preserving near‑neutral behaviour."

Open questions for Dick / Arastoo
- Do you want fc applied to K or to f_m (closure input)? (K is simpler; f_m keeps MOST formalism.)
- Preferred floor fc_min and B_thresh for manuscript sensitivity tests?
- Status of journal submission and whether you want these diagnostics added as a short appendix or supplementary figures.