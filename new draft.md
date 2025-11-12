# Surface Roughness, Geometric Mean Heights, and Drag Coefficient Bias
Title: Geometric Mean Heights in Logarithmic Boundary Layers: Implications for Bulk Transfer Coefficients and Richardson Number Reconstruction

Target Journal: Quarterly Journal of the Royal Meteorological Society or Journal of Geophysical Research: Atmospheres

Abstract
Coarse vertical discretization in the stable surface layer introduces systematic bias when arithmetic mean heights are used to represent logarithmic profiles. Using convexity of the logarithm and midpoint accuracy in log-coordinates, we (i) prove that arithmetic mean heights bias neutral drag coefficients C_D and bulk Richardson number Ri_b, (ii) show that the geometric mean height z_g = √(z₁z₂) is the midpoint in ln z and minimizes layer-reconstruction error for log/power-law structures, and (iii) quantify additional bias from heat–momentum roughness contrasts through kB⁻¹ = ln(z₀/z₀_h). We provide a practical workflow for calibrating z₀, z₀_h from tower and remote sensing, and deliver SBL-focused diagnostics and algorithms suitable for operations.

1. Introduction
- Problem: Layer-mean approximations in z introduce bias for log-like profiles typical of the SBL, degrading bulk transfers and stability inference.
- Key idea: Work in s = ln(z − d) where surface-layer mean states are near-linear; then the representative height is the geometric mean in z, not the arithmetic mean.
- Contributions:
  1) Bias sign and magnitude via Jensen’s inequality for ln.
  2) Representative height theory: geometric mean (midpoint in ln z) vs logarithmic mean (exact for gradients).
  3) SBL implications with distinct roughness for heat and momentum.

2. Theory
2.1 Preliminaries and notation
- Displacement height: d (urban/canopy); define ζ̃ = (z − d).
- Neutral/log wind: U(z) = (u*/κ) ln(ζ̃/z₀) + B; temperature analog with z₀_h.
- Geometric mean: z_g = √(z₁z₂); in displaced coordinates: ζ̃_g = √(ζ̃₁ζ̃₂).
- Logarithmic mean: L(ζ̃₁, ζ̃₂) = (ζ̃₂ − ζ̃₁)/ln(ζ̃₂/ζ̃₁).

2.2 Arithmetic vs geometric mean: bias via Jensen
- For concave ln(·), Jensen gives:
  ln((ζ̃₁ + ζ̃₂)/2) ≥ (ln ζ̃₁ + ln ζ̃₂)/2 = ln(√(ζ̃₁ζ̃₂)).
- Neutral drag at a “representative” height ζ̃_r:
  C_D(ζ̃_r) = [κ / ln(ζ̃_r/z₀)]².
- Using ζ̃_r = (ζ̃₁ + ζ̃₂)/2 yields ln denominator ≥ ln(ζ̃_g/z₀) ⇒ C_D(arith) ≤ C_D(geom).
  Therefore arithmetic-mean height underestimates drag and surface stress.

Bias ratio:
- B_CD = C_D(ζ̃_geom)/C_D(ζ̃_arith) = [ ln((ζ̃_arith)/z₀) / ln((ζ̃_geom)/z₀) ]² > 1.

2.3 Representative height for gradients: the logarithmic mean
- Exact layer wind increment for log layer:
  ΔU = (u*/κ) ln(ζ̃₂/ζ̃₁).
- First-order reconstruction using a point gradient:
  ΔU ≈ (∂U/∂z)|_{ζ̃_r} Δz = (u*/κ) Δz / ζ̃_r.
- Matching exact increment gives ζ̃_r = L(ζ̃₁, ζ̃₂) (logarithmic mean).
- Relation to geometric mean: For thin layers (r = ζ̃₂/ζ̃₁ → 1),
  L = ζ̃_g × [(r − 1)/(√r ln r)] = ζ̃_g × [1 + O((ln r)²)].
  Thus ζ̃_g is a second-order accurate proxy for ζ̃_r in thin layers and is the exact midpoint in ln z.

2.4 Bulk Richardson number reconstruction
- Definition across [z₁, z₂]:
  Ri_b = (g/θ̄) (Δθ Δz) / (ΔU)².
- With log wind and near-linear θ in ln z, evaluating point gradients at ζ̃_g minimizes second-order error in ln-space midpoint rule, while using ζ̃_arith yields a consistent positive bias in Ri_b (weaker shear estimate).
- Practical rule:
  - Use ζ̃_r = L for single-interval gradient matching (exact ΔU).
  - Use ζ̃_g for midpoint-based multi-level reconstructions, cubic-accurate in ln z for smooth profiles.

2.5 Heat–momentum roughness contrast (kB⁻¹)
- kB⁻¹ = ln(z₀/z₀_h) > 0 in most SBL cases (z₀_h < z₀).
- Neutral transfer coefficients:
  C_D(z) = [κ/ln(ζ̃/z₀)]²,    C_H(z) = [κ/ln(ζ̃/z₀_h)]·[κ/ln(ζ̃/z₀)].
- At identical representative height, C_H/C_D ≈ ln(ζ̃/z₀)/ln(ζ̃/z₀_h) > 1 if z₀_h < z₀.
- Consequence: Using inconsistent representative heights between momentum and heat exacerbates bias in Ri_b and in θ_* retrievals, especially in SBL with large kB⁻¹ (sea ice/snow, polar night).

3. Practical workflow (SBL-focused)
3.1 Inputs
- Multi-level tower: U(z), θ(z), z; roughness estimates or priors for z₀, z₀_h; displacement d (if canopy/urban).
- Remote sensing: Doppler lidar (U), microwave radiometer/temperature profiler (θ) for cross-checks.

3.2 Steps
1) Displacement correction: ζ̃_i = z_i − d (d from canopy metrics or fit).
2) Choose interval [z₁, z₂] (first layer or arbitrary pair).
3) Compute representative heights:
   - ζ̃_geom = √(ζ̃₁ζ̃₂) (midpoint in ln z),
   - ζ̃_logm = (ζ̃₂ − ζ̃₁)/ln(ζ̃₂/ζ̃₁) (exact for gradient matching).
4) Bulk coefficients:
   - C_D^geom = [κ/ln(ζ̃_geom/z₀)]²,
   - C_H^geom = κ² / [ln(ζ̃_geom/z₀) ln(ζ̃_geom/z₀_h)].
5) Ri_b with exact ΔU, Δθ:
   - ΔU = (u*/κ) ln(ζ̃₂/ζ̃₁)  (if u* available) or ΔU = |U₂ − U₁|,
   - Δθ = θ₂ − θ₁,
   - Ri_b = (g/θ̄)(Δθ Δz)/(ΔU)².
6) Bias diagnostics:
   - B_CD = [ln(ζ̃_arith/z₀)/ln(ζ̃_geom/z₀)]²,
   - B_Ri ≈ [Δz/ζ̃_arith]/ln(ζ̃₂/ζ̃₁)]² (first-order shear error proxy).
7) kB⁻¹ calibration:
   - Fit z₀, z₀_h by minimizing mismatch between observed and model transfer coefficients over stable segments (ζ > 0.05),
   - Enforce physical priors: z₀_h ≤ z₀; regularize by roughness similarity class.

3.3 Robust choices
- Prefer ζ̃_logm for single-interval ΔU reconstructions (exact).
- Prefer ζ̃_geom when using midpoint quadrature in ln z (multi-level stacks; second-order accurate, simple).
- Always use the same representative height for momentum and heat; differences in C_D vs C_H should arise only from z₀ vs z₀_h.

4. Numerical examples and expected magnitudes
- Example 1 (neutral, sea ice): z₀ = 1 mm, z₀_h = 0.1 mm, z₁=10 m, z₂=30 m, d=0:
  - ζ̃_arith=20 m, ζ̃_geom≈17.32 m, ζ̃_logm≈19.10 m,
  - B_CD≈[ln(20/0.001)/ln(17.32/0.001)]²≈[9.90/9.76]²≈1.03 → 3% drag underestimate using arithmetic mean.
- Example 2 (urban canopy, d=8 m, z₁=20 m, z₂=40 m, z₀=0.5 m):
  - ζ̃_arith=22 m, ζ̃_geom≈√(12·32)=19.6 m, ζ̃_logm≈(32−12)/ln(32/12)=20.7 m,
  - B_CD≈[ln(22/0.5)/ln(19.6/0.5)]²≈[3.78/3.67]²≈1.06 → 6% bias.

5. Validation plan
- Towers: CASES-99, Cabauw, ARM NSA (stable nights). Metrics: stress bias, Ri_b bias, retrieval consistency of z₀, z₀_h.
- LES: Neutral and weakly stable ensembles with known z₀, z₀_h. Recover bias scaling vs layer aspect ratio r = z₂/z₁.
- Polar focus: SHEBA/MOSAiC sites to probe large kB⁻¹ and thin roughness sublayers.
- Success criteria: |Bias(C_D)| < 3% for r ≤ 3 using ζ̃_geom; Ri_b RMSE reduction ≥ 20% vs arithmetic-mean baseline.

6. Extension hooks
- Urban/canopy displacement d:
  - Use ζ̃ = z − d in all formulas; estimate d from morphology or inertial-layer fits.
- Sea-ice/snow cool-skin:
  - Adjust z₀_h via cool-skin models; propagate to C_H via ln(ζ̃/z₀_h).
- Multi-level drag profiles:
  - Define C_D(z_k+½) on staggered ln z grid with ζ̃_geom at interfaces; ensures second-order accuracy.
- SBL coupling:
  - Blend with curvature-aware SBL scheme by evaluating φ functions at ζ̃_geom, ensuring consistent neutral curvature constraints.

7. Algorithm sketches
7.1 Representative height and coefficients (Python-like)
```python
def repr_heights(z1, z2, d=0.0):
    zt1, zt2 = z1 - d, z2 - d
    z_geom = (zt1 * zt2) ** 0.5
    z_logm = (zt2 - zt1) / math.log(zt2 / zt1)
    return z_geom, z_logm

def coeffs(z_rep, z0, z0h, kappa=0.4):
    Lm = math.log(z_rep / z0)
    Lh = math.log(z_rep / z0h)
    CD = (kappa / Lm) ** 2
    CH = kappa**2 / (Lm * Lh)
    return CD, CH
```

7.2 Bias diagnostics
```python
def cd_bias(z1, z2, z0, d=0.0):
    zt1, zt2 = z1 - d, z2 - d
    z_arith = 0.5 * (zt1 + zt2)
    z_geom = (zt1 * zt2) ** 0.5
    return (math.log(z_arith / z0) / math.log(z_geom / z0)) ** 2
```

8. Proof sketches (appendix)
- A1. Jensen bias:
  - For concave ln, ln(E[ζ̃]) ≥ E[ln ζ̃]. Substitute E over {ζ̃₁, ζ̃₂} with equal weights to establish C_D(arith) ≤ C_D(geom).
- A2. Logarithmic mean for exact ΔU:
  - Solve ΔU = ∫(u*/κ)(dz/ζ̃) = (u*/κ) ln(ζ̃₂/ζ̃₁) = (u*/κ)(Δz/ζ̃_r) ⇒ ζ̃_r = L.
- A3. Midpoint optimality in ln z:
  - Let s = ln ζ̃, U(s) = As + B. Midpoint rule is exact; s_m = (s₁+s₂)/2 ⇒ ζ̃_m = exp(s_m) = √(ζ̃₁ζ̃₂).
  - For smooth deviations (power-law z^p), midpoint error is O((Δs)³), minimized at s_m.

9. Implications and recommendations
- Use ζ̃_logm for exact single-interval gradient matching; use ζ̃_geom for multi-level midpoint/quadrature in ln z.
- Always apply the same representative height for momentum and heat; differences enter only via z₀ vs z₀_h (kB⁻¹).
- Replace arithmetic-mean representative heights in codes handling bulk transfers, Ri_b reconstruction, and stability diagnostics.

10. Reproducibility
- Provide scripts to compute geometric/logarithmic means, coefficients, and biases; pair with tower/LES example notebooks.
- Unit tests: exactness of ΔU with ζ̃_logm; Jensen inequality checks for B_CD > 1 under varying r.

Keywords: geometric mean, logarithmic mean, displacement height, roughness length, kB⁻¹, bulk Richardson number, drag coefficient, stable boundary layer
