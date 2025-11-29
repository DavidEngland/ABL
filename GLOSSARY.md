# ABL Project Glossary — concise definitions

This glossary records the canonical names, symbols, and one-line meanings used across the ABL repository.

- Obukhov length (L)  
  Definition: L = −u_*^3 θ / (κ g <w'θ'>). Scale where shear and buoyancy balance. Units: m.

- Dimensionless height (ζ)  
  Definition: ζ = z / L. Positive (ζ>0) indicates stable conditions.

- Gradient Richardson number (Ri_g)  
  Definition: Ri_g(z) = (g/θ) (∂θ/∂z) / (∂U/∂z)^2. Local stability diagnostic (unitless).

- Bulk Richardson number (Ri_b)  
  Definition (layer z0→z1): Ri_b = (g/θ̄) (Δθ · Δz) / (ΔU)^2. Layer-averaged stability (unitless).

- Geometric mean height (z_g)  
  Definition: z_g = √(z0 z1). Recommended representative height for log-like profiles and curvature diagnostics.

- Logarithmic mean height (z_L)  
  Definition: z_L = (z1 − z0)/ln(z1/z0). Exact representative height for ΔU in a log wind law.

- Bias Ratio / Richardson Number Bias (B) — primary project term  
  Preferred name: Bias Ratio (B) — also called "Richardson Number Bias".  
  Formula: B = Ri_g(z_g) / Ri_b (unitless).  
  Interpretation: B > 1 indicates concave‑down Ri_g over the layer (typical SBL) and that the bulk value Ri_b underestimates local stability; larger B → larger coarse-grid bias. Use B as the trigger/metric for curvature-aware corrections and for reporting bias-reduction performance.

- Neutral curvature invariant (Δ)  
  Definition: Δ = a_h − 2 a_m (for linear-stable φ_m/h: φ ≈ 1 + a ζ). Neutral curvature: (d^2Ri_g/dζ^2)|_0 = 2Δ. Sign of Δ sets concavity near neutrality.

- Curvature (d^2Ri_g/dζ^2)  
  Definition: Second derivative of Ri_g with respect to ζ. Controls Jensen/bulk vs. point bias and inflection behavior.

- Grid damping factor (G)  
  Definition: Multiplicative modifier G(ζ, Δz) applied to K (or φ) to reduce curvature-induced bias while preserving neutral behavior (G(0,Δz)=1 and ∂G/∂ζ|_0=0).

- Eddy diffusivities (K_m, K_h)  
  Definition: Momentum and heat diffusivities used by first-order closures. Units: m^2 s^-1.

- Critical Richardson number (Ri_c, Ri_c*)  
  Definition: Threshold used by closures for turbulence collapse; Ri_c* denotes a dynamic, state-dependent critical value.

- Turbulent kinetic energy (TKE)  
  Definition: 0.5*(u'^2+v'^2+w'^2) — used for turbulence memory and intermittency diagnostics.

Notes:
- All Ri-type diagnostics are unitless.
- When reporting B, always state the layer bounds (z0→z1) and the method used to compute Ri_b (bulk formula vs numerical integral).
- Recommended notation in manuscripts: use "Bias Ratio (B)" on first mention, then B thereafter; include formula and representative heights in a caption or glossary entry.

