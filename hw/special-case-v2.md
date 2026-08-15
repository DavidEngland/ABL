# Homework: Special case Ri_g(ζ) = ζ / (1 + β ζ) — derivations and diagnostics

Learning objectives
- Derive Taylor series, first/second derivatives, and curvature for the rational Ri_g.
- Compute the analytic layer‑average (Ri_b as an integral) for a finite layer and compare to point Ri_g at a representative height.
- Relate the analytic integral result to the common finite‑difference/bulk formula Ri_b_fd = (g/θ) Δθ Δz / (ΔU)² and identify sources of bias.
- Quantify the bias ratio B = Ri_g(z_rep)/Ri_b (choose z_rep = geometric mean) and obtain leading‑order asymptotics for thin layers.
- Perform simple numeric experiments and plots to illustrate bias dependence on β and Δz, and test a curvature‑preserving correction prototype.

Problem set

A. Symbolic derivations (analytic)
1. Starting from
   \[
   Ri_g(\zeta) = \frac{\zeta}{1+\beta\zeta},
   \]
   derive the Taylor series up to O(ζ³). Show the first and second derivatives:
   \[
   Ri_g'(\zeta)=\frac{1}{(1+\beta\zeta)^2},\qquad
   Ri_g''(\zeta)=-\frac{2\beta}{(1+\beta\zeta)^3}.
   \]
   Comment on sign of curvature for β>0.

2. For a finite nondimensional layer ζ∈[ζ0,ζ1] (ζ0>0), derive the layer average
   \[
   Ri_b=\frac{1}{\Delta\zeta}\int_{\zeta_0}^{\zeta_1}\frac{\xi}{1+\beta\xi}\,d\xi
   \]
   in closed form. (Hint: perform the integral explicitly; result involves ln(1+βξ).)

3. Define the geometric‑mean representative ζ_g = √(ζ0 ζ1). Write expression for point value Ri_g(ζ_g). Define bias ratio
   \[
   B=\frac{Ri_g(\zeta_g)}{Ri_b}.
   \]

B. Small-layer asymptotics (analytical)
4. Let ζ0 = ζ̄ − Δζ/2, ζ1 = ζ̄ + Δζ/2 with Δζ ≪ ζ̄. Expand Ri_b and Ri_g(ζ_g) in powers of Δζ and obtain leading term for B−1 (show the bias scales with Δζ²/ζ̄ terms). Interpret result: how does bias depend on layer thickness and curvature?

C. Finite‑difference (practical) comparison
5. Starting from sampled values at two z‑levels (physical heights z0, z1), the common finite‑difference (bulk) estimator is
   \[
   Ri_{b,\text{fd}}=\frac{g}{\theta_{\rm ref}}\frac{\Delta\theta\ \Delta z}{(\Delta U)^2}.
   \]
   (i) For the analytic Ri_g form, derive the corresponding Δθ and ΔU across the layer in terms of ζ mapping (assume L known or L=1 for nondimensional demo). (ii) Compare Ri_{b,fd} to the integral Ri_b from part A. Identify algebraic discrepancy terms and show how shear reconstruction (use of log‑mean z_L vs arithmetic mean) affects the denominator.

D. Numerical exercises (coding / plotting)
6. For β = 1, 4.7, 10 and several representative layers (ζ0, ζ1) corresponding to Δz = 5, 10, 50, 100 m after scaling by a chosen L (e.g., L=50 m):
   - Compute Ri_g(ζ_g), Ri_b (analytic integral), Ri_{b,fd} (finite difference with ΔU computed either directly from U(z)=seed profile or from log‑law using u*).
   - Tabulate B = Ri_g(ζ_g)/Ri_b and B_fd = Ri_g(ζ_g)/Ri_{b,fd}.
   - Plot B vs Δz for each β; plot Ri profiles (fine reference) and layer averages before/after a simple correction.

7. Implement a neutral‑preserving multiplicative damping
   \[
   G(\zeta,\Delta z) = \exp\!\Big[-D\Big(\frac{\Delta z}{\Delta z_{\rm ref}}\Big)^p\Big(\frac{\zeta}{\zeta_{\rm ref}}\Big)^q\Big],
   \]
   choose p=1, q=2, Δz_ref=10 m, ζ_ref=0.5 and test D values (0.3, 0.6, 1.0). Show how applying G to K (or equivalently to Ri_g for demonstration) changes Ri_b and B for coarse Δz.

E. Short essay (1 page)
8. Discuss practical implications: when is analytic integral preferable to bulk finite‑difference? When must you prefer log‑mean shear reconstruction (z_L) vs geometric mean (z_g)? Suggest model implementation guidance (what to compute in model per time step to detect and reduce curvature bias automatically).

Hints & short answers (instructor hints; students should attempt before reading)
- Integral (part A.2) result:
  \[
  Ri_b=\frac{1}{\Delta\zeta\,\beta^2}\Big(\beta\zeta_1-\beta\zeta_0 - \ln\frac{1+\beta\zeta_1}{1+\beta\zeta_0}\Big).
  \]
  (Set ζ0=0 as limiting case only with care — avoid starting layer at zero in practical exercises.)

- Small‑layer expansion (part B): for symmetric small Δζ about ζ̄,
  \[
  Ri_b \approx Ri_g(\zetā) - \tfrac{1}{24}\,Ri_g''(\zetā)\,(\Delta\zeta)^2 + O(\Delta\zeta^4),
  \]
  so bias ∼ −(Ri_g''/24)Δζ². Since Ri_g''<0 for β>0, Ri_b < Ri_g ⇒ B>1. This shows leading dependence ∝Δζ².

- For forward finite‑difference (part C) if ΔU is reconstructed using arithmetic mean height, denominator is biased high for log‑like wind → Ri_{b,fd} biased low. Use z_L = (z1−z0)/ln(z1/z0) for exact log-law shear reconstruction.

- Numeric sanity checks: for small Δz the integral and finite‑difference converge; for large Δz they diverge and B increases.

Suggested deliverables
- Short pdf or notebook with derivations, numeric tables, and three figures: (i) Ri profiles fine vs coarse, (ii) B vs Δz for chosen β, (iii) effect of G correction on B.
- One‑paragraph conclusion with recommended default choices for operational models (compute Ri_g at z_g, use Simpson integration for Ri_b when profile known, otherwise use log‑mean for shear in ΔU).

Minimal example code pseudocode (student starter)
```python
# filepath: hw/example_calc.py  (not required; for students)
import numpy as np
def Ri_g(zeta,beta): return zeta/(1+beta*zeta)
def Ri_b_integral(z0,z1,beta):
    bz = lambda xi: xi/(1+beta*xi)
    # analytic form: (1/((z1-z0)*beta**2))*(beta*z1 - beta*z0 - log((1+beta*z1)/(1+beta*z0)))
    return (1.0/((z1-z0)*beta**2))*(beta*z1 - beta*z0 - np.log((1+beta*z1)/(1+beta*z0)))
# numeric test
beta=4.7; z0=0.1; z1=2.0
z_g = np.sqrt(z0*z1)
print("Ri_g(zg)", Ri_g(z_g,beta))
print("Ri_b", Ri_b_integral(z0,z1,beta))
```

Instructor notes (brief)
- Emphasize that the integral derivation is more informative than the black‑box bulk formula because it isolates curvature contributions and enables small‑Δ asymptotics; comparing analytic Ri_b to Ri_{b,fd} reveals where shear‑denominator reconstruction errors dominate (use z_L vs arithmetic mean).
- Encourage students to use both analytic and numeric approaches: algebraic series + small‑Δ expansion gives insight; plots show practical magnitude.

References / further reading
- Jensen inequality & representative heights: geometric vs log mean discussion in the course notes.
- England & McNider (1995), Businger et al. (1971) for MOST forms.
- GABLS intercomparison for LES benchmarks (practical validation).