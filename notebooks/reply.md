Dick —

Thanks for the read and clear summary.  Short answer: implement a multiplicative correction fc to the stability function fm that depends on the observed bulk-vs-point bias B (and optionally on Δz).  When the finite-difference (bulk) Ri grows relative to the point Ri at z_g (B>1), reduce fm so the stability function has a longer tail (more mixing) at large Ri.  If Ri_g is linear (B≈1) fc→1 and no correction is applied.

Minimal pseudocode (drop into the diffusion routine prior to K calculation):

    # compute bulk and point Ri for the layer
    Ri_b = compute_bulk_Ri(z0,z1,theta,U)
    z_g = sqrt(z0*z1)
    Ri_g = compute_point_Ri(z_g, profile)   # use local gradients or analytic estimate
    B = Ri_g / Ri_b if Ri_b>0 else 1.0

    # baseline correction - simple rational form
    alpha = 1.0            # tuning parameter (start ~1)
    fc = 1.0 / (1.0 + alpha * max(0.0, B - 1.0))

    # optional: make correction stronger for coarser layers
    # Dz_ref = 10.0  # m reference
    # p = 1.0
    # fc = exp( -alpha * (B - 1.0) * (Dz / Dz_ref)**p )

    # apply correction to stability function used to compute K
    fm_corrected = fm * fc
    K_m = compute_Km(fm_corrected, ...)

Tuning notes
- Start with alpha=1.0; if correction is too strong reduce alpha (0.3–0.7), if too weak increase (1.5–2.0).
- Use fc exponential (fc = exp(−alpha*(B−1))) if you prefer strictly positive smooth decay.
- Preserve neutral behavior by ensuring fc( B≈1 ) ≈ 1 and dfc/dζ|_{0} ≈ 0 (use q≥2 if making ζ-dependent).
- Log results for a few nights: report B distribution and percent change in K at z_g to pick alpha.

Why this is simple and robust
- Directly targets the observed bias B you described.
- No iterative inversion required.
- If Ri_g is linear with height (no curvature), B≈1 and fc≈1 so no change.
- Easily implemented and tested on one-night cases (ARM/SHEBA/LES) before operational adoption.

If you want I will:
- produce a one-file patch that implements fc in your demo notebook and prints B, fc and K change for the synthetic profile, or
- create a short note with recommended alpha/p values from a few example nights.

— Dave (drafted reply)