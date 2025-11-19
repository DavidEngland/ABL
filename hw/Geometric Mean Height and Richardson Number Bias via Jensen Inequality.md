# Homework: Geometric Mean Height & Richardson Number Bias (Jensen's Inequality)

Course: Boundary-Layer Meteorology / Advanced Atmospheric Physics  
Level: Undergraduate and above

---

## Problem statement (brief)
Atmospheric models compute a bulk Richardson number across a layer,
$$
Ri_b = \frac{g}{\theta}\frac{\Delta\theta\,\Delta z}{(\Delta U)^2},
$$
while the local (gradient) Richardson number is
$$
Ri_g(z)=\frac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2}.
$$
Show that when \(Ri_g(z)\) is concave-down over \([z_0,z_1]\) (typical SBL), the layer average underestimates the point value at the geometric-mean height \(z_g=\sqrt{z_0 z_1}\):
$$
Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_g).
$$

---

## Part A — Jensen's inequality & concavity (warm-up)

A1. State Jensen's inequality for a concave function \(f\):
For concave \(f\) on \([a,b]\),
$$
f\!\left(\frac{1}{b-a}\int_a^b x\,dx\right) \ge \frac{1}{b-a}\int_a^b f(x)\,dx.
$$
Concavity (twice differentiable): \(f''(x)\le 0\) on \([a,b]\).

A2. Show \(\ln z\) is concave for \(z>0\):
$$
\frac{d^2}{dz^2}\ln z = -\frac{1}{z^2}<0,\quad z>0.
$$

---

## Part B — Why the geometric mean?

B1. Logarithmic mean for log wind:
For \(U(z)=\frac{u_*}{\kappa}\ln(z/z_0)+C\),
the exact layer-averaged gradient equals the point gradient at the logarithmic mean
$$
\bar z_{\ln}=\frac{z_1-z_0}{\ln(z_1)-\ln(z_0)}.
$$

B2. Geometric mean is midpoint in log space:
$$
\ln z_g=\tfrac12(\ln z_0+\ln z_1)\quad\Rightarrow\quad z_g=\sqrt{z_0 z_1}.
$$

B3. For thin layers (\(z_1/z_0\to1\)) the log-mean and geometric mean coincide to \(O((\ln(z_1/z_0))^2)\); hence \(z_g\) is a practical representative height for log-like profiles.

---

## Part C — Concave-down \(Ri_g\) and bias proof

C1. Apply Jensen to concave \(Ri_g(z)\):
Since \(Ri_g\) concave → by Jensen
$$
\frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz \le Ri_g\!\left(\frac{1}{\Delta z}\int_{z_0}^{z_1} z\,dz\right) = Ri_g\!\left(\frac{z_0+z_1}{2}\right).
$$
So \(Ri_b < Ri_g\) at the arithmetic mean.

C2. Use concavity of \(\ln z\) to compare means:
Because \(\ln\) is concave,
$$
\ln z_g = \tfrac12(\ln z_0+\ln z_1) > \ln\!\left(\tfrac{z_0+z_1}{2}\right),
$$
hence \(z_g < (z_0+z_1)/2\). For typical SBL where \(Ri_g\) decreases with height (near-surface), this implies
$$
Ri_g(z_g) > Ri_g\!\left(\tfrac{z_0+z_1}{2}\right),
$$
so \(Ri_b < Ri_g(z_g)\).

C3. Alternative, cleaner route: change variable \(s=\ln z\). If \(\tilde Ri(s)=Ri_g(e^s)\) is concave in \(s\), apply Jensen on \(s\)-interval to get
$$
\frac{1}{\ln(z_1/z_0)}\int_{\ln z_0}^{\ln z_1}\tilde Ri(s)\,ds < \tilde Ri\!\left(\tfrac{\ln z_0+\ln z_1}{2}\right)=Ri_g(z_g).
$$
For thin layers the log-weighted average approximates the arithmetic layer average, yielding \(Ri_b \lesssim Ri_g(z_g)\).

---

## Part D — Numerical verification (analytic quadratic example)

Assume \(Ri_g(z)=c_1 z + c_2 z^2\) with \(c_1>0, c_2<0\).
Example: \(z_0=10\) m, \(z_1=100\) m, \(c_1=0.01\), \(c_2=-0.0001\).
Compute:
- \(z_g=\sqrt{z_0 z_1}\).
- \(Ri_g(z_g)=c_1 z_g + c_2 z_g^2\).
- Analytic integral:
$$
Ri_b=\frac{1}{\Delta z}\int_{z_0}^{z_1}(c_1 z + c_2 z^2)\,dz
=\frac{1}{\Delta z}\left[\tfrac{c_1}{2}(z_1^2-z_0^2)+\tfrac{c_2}{3}(z_1^3-z_0^3)\right].
$$

Short Python snippet to reproduce:
```python
# compute example values
import math
z0, z1 = 10.0, 100.0
c1, c2 = 0.01, -0.0001
zg = math.sqrt(z0*z1)
Rig_zg = c1*zg + c2*zg*zg
Dz = z1 - z0
Rb = (0.5*c1*(z1*z1 - z0*z0) + (c2/3.0)*(z1**3 - z0**3))/Dz
print("z_g =", zg)
print("Ri_g(z_g) =", Rig_zg)
print("Ri_b =", Rb)
```
You should observe `Ri_b < Ri_g(z_g)` for the example.

---

## Part E — Interpretation and practical correction

E1. Why underestimation leads to overmixing: models map \(Ri\) → eddy diffusivities \(K_m,K_h\) (decreasing functions of \(Ri\)). If \(Ri_b\) underestimates true local stability, computed \(K\) are too large → excessive vertical mixing.

E2. Using arithmetic mean height worsens the bias (since \(z_g<z_{\text{arith}}\) and \(Ri_g\) typically decreases with height in SBL).

E3. Practical correction (no grid change): evaluate representative Ri at \(z_g\) (or use log-mean) when computing stability functions; or apply a curvature-aware multiplicative correction \(G(\zeta,\Delta z)\) to \(K\) (ensure \(G(0)=1\), \(G'(0)=0\)).

---

## Optional challenge (derivation hint)

Taylor expand \(Ri_g\) about \(z_g\):
$$
Ri_g(z)\approx Ri_g(z_g) + Ri_g'(z_g)(z-z_g) + \tfrac12 Ri_g''(z_g)(z-z_g)^2.
$$
Integrate over \([z_0,z_1]\); the linear term integrates to zero about the symmetric point \(z_g\), leaving a quadratic bias term proportional to \(Ri_g''(z_g)\). Use this to construct a leading-order correction factor \(G\).

---

## Instructor Hints, References, and Student Resources

Brief purpose
- Provide students practical tools to reproduce examples, extend the exercise, and run small projects using tower / LES / reanalysis data.

Recommended datasets (public)
- ARM/SGP (https://www.arm.gov) — tower + sondes; good for stable-night case studies.
- SHEBA archive (Polar Data) — Arctic stable BL examples.
- CASES-99 (Kansas) — intermittent turbulence / LLJ cases.
- Cabauw (KNMI) — long-term tower profiles.
- GABLS LES cases (https://climate.northwestern.edu/gabls/) — idealized high-resolution truth.
- ERA5 (Copernicus) — coarse reanalysis for synthetic aggregation tests.

Visual aids & figures to assign
- Plot Ri_g(z) fine vs coarse aggregation with shaded layer [z0,z1] and mark z_g.
- Show bias ratio B vs Δz for a set of cases (line plot).
- Show Taylor-series seed inversion (ζ seed → Newton) convergence plot.
- Simple diagram: geometric vs arithmetic vs log means on ln(z) axis.

Possible student projects (short)
1. Numeric verification: reproduce Part D example and sweep z0,z1 to map B(Δz).
2. Tower/remote-sensing case: compute Ri_g profiles from ARM or Cabauw and quantify B distribution for nighttime hours.
3. LES aggregation: aggregate LES profiles to coarse Δz and test curvature-aware multiplicative correction G.
4. Dynamic Ri_c test: propose and test a simple Ri_c*(Γ) form on one CASES-99 night.

Suggested deliverables
- Short report (2–4 pages): method, plots, discussion.
- Notebook (Jupyter) with code to reproduce figures.
- Optional: small script to compute B per time step from provided profile.

Grading rubric (out of 100)
- Correct derivation and explanation (20)
- Numerical implementation correctness (30)
- Figures clarity and interpretation (20)
- Reproducibility (notebook runs end-to-end) (15)
- Novel extension or thoughtful discussion (15)

Quick project timeline (suggested for 1–2 week assignment)
- Day 1–2: data ingestion + basic diagnostics (Ri_g, Ri_b).
- Day 3–4: core plots and verification of Jensen inequality numerically.
- Day 5–7: extension (aggregation sweep or short case study).
- Day 8–10: polish figures and notebook; prepare short report.

Pseudocode — minimal Python (for notebook)
```python
# compute representative heights and bias B
import numpy as np

def z_geom(z0,z1): return np.sqrt(z0*z1)
def z_logm(z0,z1): return (z1-z0)/np.log(z1/z0)

def ri_point(z,theta,U,g=9.81):
    # simple centered difference estimates (assumes arrays, index k)
    # implement robustly in notebook with boundary checks
    dz = z[1:]-z[:-1]            # interface dz
    dtheta = theta[1:]-theta[:-1]
    dU = U[1:]-U[:-1]
    Ri_b = (g/theta[:-1]) * (dtheta*dz) / (dU**2)
    return Ri_b

# example: analytic quadratic case from hw
z0,z1 = 10.0,100.0
zg = z_geom(z0,z1)
# compute Ri_g(zg) and analytic Ri_b per homework formula
```

Pseudocode — minimal Julia (for variety)
```julia
# Julia snippet
z_geom(z0,z1) = sqrt(z0*z1)
function analytic_Ri_b(z0,z1,c1,c2)
    dz = z1 - z0
    I = 0.5*c1*(z1^2 - z0^2) + (c2/3.0)*(z1^3 - z0^3)
    return I/dz
end
```

Notebook structure (recommended)
- Cell 1 (markdown): Problem statement & objectives.
- Cell 2 (code): imports and helper functions.
- Cell 3 (code): analytic quadratic verification (Part D).
- Cell 4 (code): ingestion of sample data (or synthetic profile).
- Cell 5 (code): compute Ri_g profiles and Ri_b, compute B.
- Cell 6 (plots): visualizations (Ri profiles, B vs Δz).
- Cell 7 (markdown): conclusions and extension notes.

Prompts for automated assistance (useful for students)
- "Compute Ri_g and Ri_b for provided profile and return bias ratio distribution."
- "Aggregate a high-resolution profile to coarser grids by averaging in z and report B(Δz)."
- "Plot Ri_g(z) at fine resolution and show trapezoid and Simpson approximations for Ri_b."

Quick references (papers & docs)
- Businger et al., 1971 — flux-profile relationships (classic).
- Holtslag & de Bruijn, 1991 — stable parameterizations & blends.
- GABLS intercomparison papers — LES benchmarks.
- Mahrt 2014 review — stably stratified BL behavior overview.
- England & McNider (1995) — Ri-based stability functions (historical).

Software & tool suggestions
- Python: numpy, scipy, xarray, matplotlib, netCDF4, pandas.
- Julia: DataFrames, Plots, NetCDF.jl.
- Notebooks: JupyterLab, Google Colab (if students lack local setups).
- Version control: git + GitHub classroom or private repo for submissions.

Accessibility & reproducibility tips
- Provide a small synthetic dataset (10–20 profiles) in repo for students to start.
- Include exact commands to install environment (requirements.txt / Pipfile / Project.toml).
- Use seeded random generators for synthetic tests to ensure reproducibility.

Assessment aids for instructors
- Provide a reference notebook with complete solution (hidden or available after submission).
- Use unit tests (pytest / Julia Test) to check key functions: z_geom, analytic_Ri_b, bias threshold detection.

Further reading & online resources
- ARM data tutorials: https://www.arm.gov/capabilities/data
- GABLS: https://climate.northwestern.edu/gabls/
- ERA5 access and examples: https://cds.climate.copernicus.eu/
- Example Jupyter notebooks (repo): add link to course repo once created.

---

## Appendix C — Additional Problems & Mini‑Projects

C1. Simpson‑rule correction and practical K adjustment (Problem)

Goal: quantify the Simpson approximation error for Ri_b and derive a simple multiplicative correction to K based on the Simpson-integral bias.

Task:
1. Given a smooth Ri_g(z) on [z0,z1], compute:
   - Simpson integral Ri_b^S = (1/6)[Ri_g(z0) + 4 Ri_g(z_g) + Ri_g(z1)].
   - Point value Ri_g(z_g).
   - Simpson bias S_err = Ri_g(z_g)/Ri_b^S.
2. Show for a quadratic expansion around z_g:
   Ri_g(z) ≈ Ri_g(z_g) + (1/2) Ri_g''(z_g)(z − z_g)^2,
   and derive analytic Ri_b^S and S_err to leading order in (Δz)^2.
3. Propose a practical correction to eddy diffusivity:
   K_corrected = K · S_err^α (choose α = 1 as baseline, discuss α tuning).
4. Numerical exercise:
   - Use the quadratic example from Part D (c1,c2). Compute Ri_b^S, Ri_g(z_g), S_err for Δz=90 m.
   - Apply K_corrected to a simple K profile (e.g., K0 = κ u_* z / φ_m) and show percent change at z_g.

Short solution hint:
- For quadratic Ri_g(z) = A + B (z − z_g)^2, Simpson integral over symmetric interval [z_g−h, z_g+h] yields exact integral: Ri_b^S = A + (B h^2)/3, while Ri_g(z_g) = A. Hence S_err = A / (A + B h^2/3) ≈ 1 − (B h^2)/(3A) to leading order.
- Translate B = (1/2) Ri_g''(z_g) and h = Δz/2 to express S_err in terms of Ri_g'' and Δz.

Pseudocode (Python):
```python
# Simpson-based correction
zg = math.sqrt(z0*z1)
Ri_zg = compute_Ri(zg)
Ri0, Ri1 = compute_Ri(z0), compute_Ri(z1)
Ri_S = (Ri0 + 4*Ri_zg + Ri1)/6.0
S_err = Ri_zg / Ri_S
K_corr = K * S_err  # baseline alpha=1
```

Grading notes (5–10 points)
- Derivation correctness (Taylor → Simpson): 4
- Numerical implementation & plot: 4
- Discussion of α tuning and limitations: 2

---

C2. Dynamic critical Richardson number Ri_c* — design & exploratory project

Motivation: fixed Ri_c (≈0.25) fails to capture hysteresis and turbulence persistence. Students design, implement, and evaluate candidate Ri_c*(t,z).

Part tasks:
1. Propose 2–3 functional forms for Ri_c* (lecture seed examples):
   - Linear inversion form:
     Ri_c* = Ri_c0 + a1 · clamp(Γ/Γ_ref,0,1)
   - TKE-memory form:
     Ri_c* = Ri_c0 + a2 · (1 − TKE/TKE_ref)
   - Combined form:
     Ri_c* = Ri_c0 + a1·min(Γ/Γ_ref,1) + a2·(1 − TKE/TKE_ref)
   Explain physical meaning of Γ (inversion lapse) and TKE memory.

2. Implementation exercise:
   - Compute Ri_c*(t) for a chosen night from tower or LES using running-window TKE and inversion estimate over the first 200 m.
   - Replace fixed Ri_c in a simple Ri-based f_m(Ri) = exp(−γ Ri / Ri_c*) and compare K, surface fluxes and turbulence flags (Ri>Ri_c*).

3. Calibration & evaluation:
   - Define objective metrics: percent reduction in flux RMSE vs observations, fraction of timesteps where Ri>Ri_c* but turbulence persisted (false positives).
   - Tune coefficients (a1,a2,Γ_ref,TKE_ref,γ) on training nights and test on held-out nights.

4. Project variants:
   - Multi-site: compare Ri_c* behavior at ARM SGP vs SHEBA (polar).
   - Model-integration demo: apply K multiplier path (Biazar style) in a 1D column and assess surface temperature bias.

Hints & diagnostics:
- Estimate inversion strength Γ as θ(z_top) − θ(z_base) over a small window (e.g., 10–50 m) above the surface.
- Use exponential smoothing for TKE memory: TKE_rel(t)=α TKE(t−1)+(1−α)TKE(t).
- Keep Ri_c* bounded (e.g., [0.1,1.5]) to avoid pathological closures.

Pseudocode (outline):
```python
# dynamic Ri_c* prototype
Ri_c0 = 0.25
Gamma_ref = 0.01  # K/m typical strong inversion
TKE_ref = 0.1     # m2/s2
a1, a2 = 0.5, 0.5

Gamma = (theta[z_top] - theta[z_base]) / (z_top - z_base)
TKE_rel = running_mean(TKE, window=10)
Ri_c_star = Ri_c0 + a1 * min(Gamma / Gamma_ref, 1.0) + a2 * (1.0 - min(TKE_rel / TKE_ref, 1.0))
Ri_c_star = np.clip(Ri_c_star, 0.1, 1.5)
```

Assessment rubric (project, out of 100)
- Clear definition & reasoning for Ri_c* form(s): 20
- Correct implementation & reproducible notebook: 30
- Calibration procedure & validation metrics: 20
- Interpretation and physical discussion: 20
- Code quality & documentation: 10

Suggested datasets and timeline
- Use one night from ARM SGP (training) and one night from SHEBA (test) — 2 week mini‑project.
- For LES-based testing use GABLS outputs to measure "truth" turbulence persistence.

---

C3. Combined capstone mini‑project (optional, advanced)

Combine C1 and C2: use Simpson-based correction to reduce instantaneous bias and implement dynamic Ri_c* to control regime switching. Evaluate combined impact on surface flux RMSE and the bias ratio B across a suite of nights (report).

Deliverable:
- Notebook, 4 figures (Ri profiles, B vs Δz, K before/after, time series of Ri_c*), short write‑up (2 pages).

Grading addendum:
- Combined project judged on improvement over baseline (uncorrected) and clarity of causal attribution (how Simpson correction vs Ri_c* each contributed).
