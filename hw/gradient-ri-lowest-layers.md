GitHub Copilot

Short answer — recommended practice
- Prefer evaluating point Ri_g for the lowest model layer at the geometric-mean height z_g = sqrt(z0*z1). That minimizes log-space bias for MOST/log-like profiles.
- Use the log-mean z_L = (z1−z0)/ln(z1/z0) when reconstructing shear ΔU for the bulk denominator (it exactly reproduces ΔU for a log law). This avoids denominator bias that inflates Ri_b low.
- If surface fluxes (u*, θ*) are available, compute Ri_g(z_g) from MOST (φ_m, φ_h at ζ=z_g/L) — best option.
- If only discrete values exist, use a forward difference for the lowest-layer gradients but evaluate Ri at z_g and compute ΔU using z_L. Smooth noisy gradients before dividing.

Steps (compact)
1. compute z_g = sqrt(z0*z1) and z_L = (z1−z0)/ln(z1/z0).
2. if u*,θ* known:
   - ζ_g = z_g/L; evaluate φ_m(ζ_g), φ_h(ζ_g)
   - Ri_g(z_g) = ζ_g * φ_h(ζ_g) / φ_m(ζ_g)^2
3. else:
   - dθ/dz ≈ (θ1−θ0)/Δz (forward); ΔU = use log-law or (U1−U0) but reconstruct shear with z_L
   - Ri_g(z_g) ≈ (g/θ_g) * (dθ/dz) / ( (ΔU/Δz_using_zL)^2 )
4. guard: if Ri_b very small or noisy, low-pass smooth before computing B = Ri_g/Ri_b.

One-line pseudocode
```python
z_g = sqrt(z0*z1); z_L = (z1-z0)/log(z1/z0)
if have_u_star:
    zeta = z_g/L; Ri_g = zeta * phi_h(zeta) / phi_m(zeta)**2
else:
    dth = (th1-th0)/Dz; dU = (U1-U0)/Dz_using_zL
    Ri_g = (g/theta_g) * dth / (dU**2)
```

Practical tips
- Use q‑preserving smoothing (small window) to reduce noise in finite differences.
- Prefer geometric mean for point diagnostics and log-mean for shear reconstruction — using arithmetic mean risks additional bias.