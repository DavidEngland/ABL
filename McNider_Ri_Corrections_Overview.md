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

GitHub Copilot

Step‑by‑step physical justification (concise)

1. Origin of the bias
   - MOST: Ri_g(ζ) = ζ·F(ζ) with F=φ_h/φ_m^2. For typical SBL φ, V_log ≡ (φ_h'/φ_h) − 2(φ_m'/φ_m) < 0 near ζ→0 so d²Ri_g/dζ²|0 = 2Δ < 0 — Ri_g is concave‑down.
   - Jensen: for a concave function f, average(f) < f(average). Thus layer average Ri_b = (1/Δz)∫ Ri_g dz < Ri_g(z_g) where z_g = √(z0 z1). Physically the curved profile lifts the point value above the layer mean → bulk underestimates local stability.
   - Denominator (shear) bias: if shear is reconstructed using arithmetic/log-approximations instead of the log‑mean z_L, ΔU tends to be overestimated (shear underestimated), further lowering Ri_b.

2. Why a multiplicative correction fc is appropriate
   - Closure uses K computed from φ (or f_m(Ri)). If Ri_b underestimates true Ri_g, the computed K is too large. Multiplying K by fc<1 reduces K and increases effective stability where coarse averaging caused spurious mixing.
   - Multiplicative form is physically minimal: it rescales turbulent transport while leaving base MOST structure intact and is easy to localize (apply only in layers where bias detected).

Design constraints for fc (physical)

- Preserve neutral physics: fc(ζ→0)=1 and ∂ζfc|0=0 so neutral slope/curvature (2Δ) unchanged.
- Grid consistency: fc→1 as Δz→0 (no change on fine grids).
- Locality/monotonicity: fc decreases with ζ and with coarsening Δz (only damp tails).
- Boundedness: avoid fc→0 (set fc_min ≳ 0.2) to prevent numerical starvation of mixing.

Simple, interpretable templates (pick one)

- Exponential (smooth, tunable tail)
  fc = exp[ − α·(B−1)·(Δz/Δz_ref)^p·(ζ/ζ_ref)^q ]
  - Physics: fractional reduction ∝ bias magnitude (B−1) scaled by coarse‑grid severity.
  - Choose q≥2 so ∂ζfc|0=0; α tunes strength.
- Rational (milder)
  fc = 1 / (1 + α·max(0,B−1)·(Δz/Δz_ref)^p·(ζ/ζ_ref)^q)
  - Physics: bounded, less aggressive at large bias.

Practical implementation notes

- Use B = Ri_g(z_g)/Ri_b as trigger and amplitude; compute Ri_g(z_g) from local gradients or MOST φ if available.
- Apply fc to K_m,K_h (or to mixing length l) in the diffusion step — do not change surface‑prescribed fluxes.
- Cap fc to avoid numerical issues; log diagnostics (B, fc, ΔK/K).
- Tune α,p,q,Δz_ref using LES/tower cases so target reduction (e.g., 30–50% of bias) is achieved without harming neutral behavior.

Validation checklist (short)

- Unit test: fc(ζ→0)=1 and ∂ζfc|0≈0.
- Case tests: coarse vs fine grids → B_before → B_after reduced toward 1.
- Monitor: surface fluxes, inversion height, no spurious oscillations.

If you want, I can:

- add a one‑cell notebook snippet that computes B and fc for the synthetic profile and shows K_old vs K_new.

### Cleaned McNider ODE and practical closed form

Notation
- fc(Δz, ζ): multiplicative correction (applied to K or f_m) for layer thickness Δz and nondimensional height ζ = z/L.
- B = Ri_g(z_g) / Ri_b : observed bias ratio for the layer (z_g = √(z0 z1)).
- Δz_ref, ζ_ref : reference scales (e.g., 10 m and 0.5).
- α ≥ 0 : tuning coefficient (strength), q ≥ 1 : ζ‑exponent to enforce neutral preservation (q≥2 recommended).

Design goal (McNider): make fc→1 as Δz→0 and reduce fc monotonically with Δz when the bias B>1.

Clean ODE (logarithmic form, simple and integrable)
\[
\frac{d\ln f_c}{d\ln\Delta z} \;=\; -\,\alpha\,(B-1)\left(\frac{\zeta}{\zeta_{\rm ref}}\right)^q .
\]

Integrating with boundary condition f_c(\Delta z_{\rm ref})=1 gives the closed form
\[
f_c(\Delta z,\zeta) \;=\; \left(\frac{\Delta z}{\Delta z_{\rm ref}}\right)^{-\alpha\,(B-1)\,(\zeta/\zeta_{\rm ref})^q}
\;=\; \exp\!\Big[-\alpha\,(B-1)\,(\tfrac{\zeta}{\zeta_{\rm ref}}^q)\,\ln\!\tfrac{\Delta z}{\Delta z_{\rm ref}}\Big].
\]

Remarks and practical defaults
- This form satisfies f_c→1 as Δz→Δz_ref (choose Δz_ref small, e.g. 10 m) and tends to 1 as Δz→0.
- Choosing q≥2 ensures ∂f_c/∂ζ|_{ζ→0}=0 (neutral‑preserving slope).  
- For small departures you can approximate with the exponential template used elsewhere:
  \[
  f_c \approx \exp\big[-\alpha'(B-1)\,(\tfrac{\Delta z}{\Delta z_{\rm ref}})^p\,(\tfrac{\zeta}{\zeta_{\rm ref}})^q\big]
  \]
  which is equivalent in spirit; pick the version that best fits tuning/implementation constraints.
- Suggested starting coefficients: α≈0.8–1.2, Δz_ref=10 m, ζ_ref=0.5, q=2. Apply floor f_c ≥ 0.2 to avoid excessive reduction.

Simple pragmatic estimator (no ODE solve required)
- Compute B for the layer. If B ≤ 1.05 → f_c = 1. Else compute f_c from the closed form above and apply K_new = K_old · f_c.
