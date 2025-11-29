# ABL Canonical Glossary

Brief, canonical definitions and formulas used across the ABL repository.

- Obukhov length (L)  
  Definition: L = −u_*^3 θ / (κ g ⟨w'θ'⟩). Scale where shear and buoyancy balance. Units: m.
  
  **Sign convention and physical meaning:**
  - **L > 0** (stable): Surface cooling (⟨w'θ'⟩ < 0); buoyancy opposes mixing; Arctic nights, nocturnal BL.
  - **L < 0** (unstable): Surface heating (⟨w'θ'⟩ > 0); buoyancy enhances mixing; convective daytime BL.
  - **Sign in equations:** Use $L$ directly (preserves sign) in ζ = z/L. Do NOT use |L| unless explicitly justified (loses regime information).
  
  **Branch selection for φ functions:**
  - For L > 0 (ζ > 0): Use **stable-branch** φ (linear or Beljaars-Holtslag forms).
  - For L < 0 (ζ < 0): Use **unstable-branch** φ (power-law $(1-β\zeta)^{-α}$ forms).
  - **Critical error:** Applying unstable power-law to ζ > 0 produces unphysical poles and divergent curvature.

- Dimensionless height (ζ)  
  Definition: ζ = z / L. Positive (ζ>0) indicates stable conditions.

- Geometric mean height (z_g)  
  Definition: z_g = √(z0 · z1). Representative height for log-like profiles; use for evaluating Ri_g point values.

- Logarithmic mean height (z_L)  
  Definition: z_L = (z1 − z0) / ln(z1 / z0). Exact representative height for reproducing ΔU for a log wind law.

- Gradient Richardson number (Ri_g)  
  Definition: Ri_g(z) = (g/θ) (∂θ/∂z) / (∂U/∂z)^2. Local (point) stability diagnostic (unitless).

- Bulk Richardson number (Ri_b)  
  Definition: Ri_b = (g/θ̄) (Δθ · Δz) / (ΔU)^2. Layer-averaged stability (unitless).

- Bias Ratio / Richardson Number Bias (B)  
  Definition: B = Ri_g(z_g) / Ri_b. Interpretation: B > 1 → Ri_b underestimates local stability (concave‑down Ri_g over layer); used to trigger curvature-aware corrections.

- Neutral curvature invariant (Δ)  
  Definition (linear-stable φ): Δ = a_h − 2 a_m. Neutral curvature: (d^2Ri_g/dζ^2)|_0 = 2Δ. Sign of Δ sets initial concavity of Ri_g near neutrality.

- Curvature (d^2Ri_g/dζ^2)  
  Compact formula: d^2Ri_g/dζ^2 = F [ 2 V_log + ζ (V_log^2 − W_log) ], with F = φ_h / φ_m^2, V_log = φ_h'/φ_h − 2 φ_m'/φ_m, W_log = dV_log/dζ. Controls bulk vs point bias and inflection behavior.

- Grid damping factor (G)  
  Definition: multiplicative modifier G(ζ, Δz) applied to K (or φ) to reduce curvature-induced bias while preserving neutral behavior (G(0,Δz) = 1 and ∂G/∂ζ|_0 = 0). Typical template: G = exp[−D (Δz/Δz_ref)^p (ζ/ζ_ref)^q] with q ≥ 2.

- Eddy diffusivities (K_m, K_h)  
  Definition: Momentum and heat diffusivities used by first-order closures. Units: m^2 s^-1.

- Critical Richardson number (Ri_c, Ri_c*)  
  Definition: Threshold used by closures for turbulence collapse; Ri_c* denotes a dynamic, state-dependent critical value (learned or parameterized).

- Turbulent kinetic energy (TKE)  
  Definition: 0.5·(u'^2 + v'^2 + w'^2). Used for turbulence memory and intermittency diagnostics.

- Omission metric (E_omit)  
  Definition: measures validity of constant-L approximation when mapping ζ-curvature to z-curvature; if E_omit ≪ 1 the simple scaling d^2/dz^2 ≈ (1/L^2) d^2/dζ^2 is acceptable.

Usage notes
- When reporting B, always state the layer bounds (z0→z1) and method used to compute Ri_b (bulk formula vs numerical integral).  
- Use z_g for point Ri_g evaluation; use z_L when reconstructing ΔU for a log wind law.  
- Preserve neutral curvature (2Δ) when applying any correction (G or φ modifier).

Short citation
- Suggested first mention in manuscripts: "Bias Ratio (B) = Ri_g(z_g)/Ri_b (see glossary)." Subsequent mentions: B.

