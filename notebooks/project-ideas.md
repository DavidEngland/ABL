Project ideas for short Jupyter notebooks

Below are compact, high‑impact notebook projects you can finish quickly and use as teaching examples, reproducible tests, or manuscript figures. Each item includes goal, key steps, expected outputs, difficulty, and data/code hints — ready to turn into a runnable notebook cell outline.

---

1. Bias diagnostic: Ri_g vs Ri_b for common MOST closures

• Goal: Compare point and bulk Ri for log‑linear, BD (USL), exponential, and polylog closures.
• Key steps: implement φm, φh; compute Ri_g(ζ); integrate Ri_b over layer [ζ0,ζ1]; compute B(ζ).
• Expected outputs: table and plot of Ri_g, Ri_b, B versus ζ; small‑Δζ asymptotic check.
• Difficulty: easy.
• Hints: vectorized NumPy; Matplotlib plots; include analytic closed forms for log‑linear.


---

2. f_c sensitivity map (McNider power law)

• Goal: Visualize correction factor f_c(Δz, ζ) sensitivity to α and q.
• Key steps: implement f_c formula; sweep Δz/Δz_ref on log scale and ζ on linear scale; surface and contour plots.
• Expected outputs: contour plot of f_c; recommended α,q ranges marked; export PNG for manuscript.
• Difficulty: easy–moderate.
• Hints: use logspace, pcolormesh, and interactive sliders (ipywidgets) for α,q.


---

3. Riemann–Stieltjes IBP demo and discrete estimator

• Goal: Demonstrate IBP identity numerically and build robust discrete estimator for ∫ z dRi_g.
• Key steps: derive discrete Stieltjes sum; compare continuous and discrete for coarse Δz; add TV‑regularized smoothing for noisy Ri_g.
• Expected outputs: diagnostic table of IBP terms; error vs grid spacing plot; notebook function ready for model integration.
• Difficulty: moderate.
• Hints: use scipy.signal for smoothing; show convergence plots.


---

4. ADM worked example with explicit G(z)

• Goal: Produce closed‑form ε0, ε1 for G(z)=G0 (constant) and demonstrate truncation error.
• Key steps: compute ε0 via double integral; evaluate ε1 = L^{-1}[-c ε0^{3/2}]; compare series to Newton iteration.
• Expected outputs: symbolic expressions (SymPy) and numeric comparisons; plot ε(z) approximations.
• Difficulty: moderate.
• Hints: use SymPy for integrals and to generate LaTeX‑ready equations.


---

5. Polylog closure exploration: s=1 to s=3

• Goal: Explore φ = exp(a Li_s(b ζ)) behavior across s and map curvature regimes.
• Key steps: implement Li_s via mpmath or scipy.special; compute Ri_g and Ri_g’’; identify Δ = q-2p analogue using coefficient mapping.
• Expected outputs: family plot of φ, Ri_g, Ri_g’’; classification table (concave up/down).
• Difficulty: moderate.
• Hints: include asymptotic series near ζ→0 and ζ→large; show where analytic continuation would introduce complex parts.


---

6. Calibration against tower/LES case (toy)

• Goal: Fit α, q, γ to match LES/tower Ri and flux profiles (toy synthetic dataset or small public case).
• Key steps: generate or load reference Ri_g/Ri_b; minimize error between corrected K·f_c flux and reference; show sensitivity.
• Expected outputs: best‑fit parameters; residual maps; recommended α,q per regime.
• Difficulty: moderate–hard.
• Hints: use scipy.optimize.least_squares; bootstrap to estimate uncertainty.


---

7. Grid refinement study: apply f_c in simple 1D model

• Goal: Embed f_c in a simple 1D diffusion model and show grid‑independence of integrated fluxes.
• Key steps: implement 1D vertical diffusion timestepper with K_old and K_new = K_old·f_c; run for multiple Δz and compare layer‑integrated fluxes and temperature tendency.
• Expected outputs: convergence plots of fluxes vs Δz with/without f_c; short animation of profile evolution.
• Difficulty: moderate.
• Hints: explicit/implicit time stepping examples; CFL discussion.


---

8. Teaching notebook: intuitive visualizations and exercises

• Goal: Produce a short interactive lesson for students illustrating curvature bias and correction.
• Key steps: narrative cells, minimal code, interactive sliders for ζ, Δz, a_m/a_h; include exercises (compute B, tune α).
• Expected outputs: shareable Binder/Colab link; succinct exercises with solutions.
• Difficulty: easy.
• Hints: use ipywidgets and markdown cells; export as HTML for teaching.


---

Quick implementation outline (one‑cell scaffold)

• Imports: numpy, scipy, mpmath/sympy (optional), matplotlib, ipywidgets.
• Functions: phi_m, phi_h, Ri_g(ζ), Ri_b(z0,z1), B(ζ), f_c(Δz,ζ).
• Plot utilities: plot_family(ζ), contour_fc(Δz,ζ).
• Save figures: high‑dpi PNG for manuscript.


---

If you want, I’ll generate a ready‑to‑run notebook for one of these projects (pick the project number). I can include executable code, plots, and LaTeX‑ready output cells you can paste into your manuscript.