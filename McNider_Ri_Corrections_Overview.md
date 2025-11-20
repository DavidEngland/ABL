# McNider — Practical Overview: Richardson‑based Grid Corrections (concise)

Purpose
- Explain a minimal, neutral‑preserving correction workflow that reduces coarse‑grid bulk-vs-point Ri bias and limits excessive mixing in the SBL.

Core idea
- Compare point gradient Ri_g at geometric mean z_g to bulk Ri_b. If Ri_b underestimates Ri_g (B>1), reduce the apparent stability used by the closure via a multiplicative factor fc(B,Δz,ζ) that →1 as ζ→0 and Δz→0.

1. Diagnostics (per layer)
- z_g = sqrt(z0*z1)
- Ri_g(z_g) — from local gradients or MOST φ if available
- Ri_b = (g/θ_ref)·(Δθ·Δz)/(ΔU^2)
- Bias: B = Ri_g(z_g) / Ri_b

2. Decision
- If B ≤ B_thresh (default 1.05) → no correction
- Else apply fc

3. Two recommended fc templates (both preserve fc→1 at small ζ; choose q≥2 to make dfc/dζ|0≈0)

A) Exponential (smooth, robust)
- fc = exp( − α * (B − 1) * (Δz/Δz_ref)^p * (ζ/ζ_ref)^q )
- defaults: α=1.0, p=1.0, q=2, Δz_ref=10 m, ζ_ref=0.5

B) Rational (bounded, milder tail)
- fc = 1 / (1 + α * max(0,B−1) * (Δz/Δz_ref)^p * (ζ/ζ_ref)^q )
- defaults: α=0.8, p=1.0, q=2

Notes
- Use ζ = z_g / L (if L unknown, use ζ proxy or omit ζ-term and rely on Δz scaling).
- Choose q≥2 to ensure neutral‑preserving slope near ζ=0.
- Set hard caps: fc_min ≥ 0.2 to avoid over‑mixing.

4. Where to apply (operational guidance)
- K‑based models (flux computed from K at first interior level): apply fc to K_m and K_h computed using fm: K* = K · fc.
- Surface‑prescribed flux models: do NOT change surface flux directly; apply fc to interior K (first interior layer) to limit vertical intrusion.
- Mixing‑length alternative: l* = l · (1 + β * (B−1) · (Δz/Δz_ref)^p)^{-1} — reduces eddy scale rather than K.

5. Minimal insertion pseudocode (drop‑in)
```python
# inputs: z0,z1,U0,U1,th0,th1,phi_m,phi_h,K_old,z,L,Delta_z
z_g = sqrt(z0*z1)
Ri_g_zg = compute_Ri_point(z_g, profile or phi_m,phi_h, L)
Ri_b = (g/theta_ref) * (th1-th0) * (z1-z0) / ((U1-U0)**2)
B = Ri_g_zg / (Ri_b if Ri_b>0 else tiny)
if B > B_thresh:
    # choose template A or B
    fc = exp( -alpha * (B-1) * (Delta_z/Delta_z_ref)**p * ( (z_g/L)/zeta_ref )**q )
    fc = max(fc, fc_min)
else:
    fc = 1.0
K_new = K_old * fc
# use K_new in diffusion update
```

6. Tuning & validation checklist
- Unit tests:
  - fc(ζ→0)=1 within numerical tolerance.
  - dfc/dζ|_{0} ≈ 0 (verify numerically).
- Case tests:
  - Fine reference (LES or Δz≈10 m) → compute B_ref ≈1.
  - Coarse baseline (Δz=60–100 m) → B_before; after applying fc ensure B_after reduced by target (e.g., 30–50% reduction toward 1).
- Monitor: surface flux change, inversion height, K percent change, and numerical stability (no oscillations).
- Logging: per timestep save B, fc, K_old/K_new for QA and tuning.

7. Default quick parameters (starting point)
- B_thresh = 1.05
- Exponential: α=1.0, p=1.0, q=2, Δz_ref=10 m, ζ_ref=0.5, fc_min=0.2
- Rational: α=0.8, same p,q and refs, fc_min=0.25

8. Practical remarks
- If φ-functions are trusted, compute Ri_g via MOST; otherwise use finite differences (centered/forward) and geometric/log mean for shear reconstruction.
- For mixed-concavity (inflection inside layer), split layer or compute fc only on concave‑down sublayer.
- Document and archive tuning cases; prefer mild correction and iterate.

References & contacts
- Short derivation and diagnostics live in Ri.md and Curvature_Demo.ipynb.
- For tuned implementation and examples, see notebooks/Curvature_Demo.ipynb and the synthetic NetCDF demo.
- Contact: David E. England — provide case specifics and I will generate tuned fc suggestions.
