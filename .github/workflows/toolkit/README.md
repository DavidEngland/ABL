Richardson number curvature toolkit specifications

This toolkit turns curvature-aware theory for gradient and bulk Richardson numbers into a reproducible, testable, and teachable package. It’s built to be solo‑developable, pedagogically annotated, and ready for validation with UAH collaborators.

---

Scope and objectives

• Problem focus: Systematic bias in Richardson number diagnostics arising from vertical curvature, discretization, and stability transitions in the atmospheric boundary layer (ABL).
• Core deliverables: Curvature-aware corrections, bias diagnostics, reference tables, and reproducible workflows for students and researchers.
• Design principles: Analytic rigor, error control, reproducibility, pedagogy-first docs, modular APIs, CI-tested, publishable notebooks.

---

Core features and algorithms

Richardson numbers and curvature-aware forms

• Gradient Richardson number:Ri_g(z)=\frac{\frac{g}{\theta_0}\,\frac{\partial \theta}{\partial z}}{\left(\frac{\partial U}{\partial z}\right)^2+\left(\frac{\partial V}{\partial z}\right)^2}

• Bulk Richardson number (layer \(z_1\) to \(z_2\)):Ri_b=\frac{\frac{g}{\theta_0}\,(\theta_2-\theta_1)\,(z_2-z_1)}{(U_2-U_1)^2+(V_2-V_1)^2}

• Curvature-aware gradient approximation (second‑order, error‑controlled):\left(\frac{\partial \phi}{\partial z}\right)_{z_0}\approx \frac{\phi_{+}-\phi_{-}}{2\Delta z}\;-\;\frac{\phi^\prime^\prime(z_0)}{6}\,\Delta z
with local curvature proxy via three-point stencil:\phi^\prime^\prime(z_0)\approx \frac{\phi_{+}-2\phi_0+\phi_{-}}{\Delta z^2}

• Curvature-corrected Ri\(_g\) estimator:\widehat{Ri_g}(z_0)=\frac{\frac{g}{\theta_0}\left[\frac{\theta_{+}-\theta_{-}}{2\Delta z}-\frac{\theta_{+}-2\theta_0+\theta_{-}}{6\Delta z}\right]}{\left[\frac{U_{+}-U_{-}}{2\Delta z}-\frac{U_{+}-2U_0+U_{-}}{6\Delta z}\right]^2+\left[\frac{V_{+}-V_{-}}{2\Delta z}-\frac{V_{+}-2V_0+V_{-}}{6\Delta z}\right]^2}

• Bias ratio diagnostics (gradient vs. bulk):\mathcal{B}=\frac{\widehat{Ri_g}(z_0)}{Ri_b(z_1,z_2)}
with stratified interpretation (e.g., stable/neutral/unstable) and confidence intervals via bootstrap over profiles.
• Correction factor ODE (curvature-aware smoothing across stability transitions):• Formulation:\frac{dC}{dz}=\alpha\,\kappa(z)\,C\;-\;\beta\,\sigma_{\text{shear}}(z)\,(C-1)
where \(\kappa(z)\) is curvature magnitude (combined from \(\theta\), \(U\), \(V\)), \(\sigma_{\text{shear}}(z)\) is shear scale, \(\alpha,\beta>0\) tunable.
• Application: Use \(C(z)\) to adjust Ri estimates:Ri_g^{(corr)}(z)=C(z)\,\widehat{Ri_g}(z),\quad Ri_b^{(corr)}=C_{\text{layer}}\;Ri_b

• Numerics: Stable explicit/implicit integrators with monotonicity constraints and parameter identifiability checks.

Stability-aware reference tables and regimes

• Regime bins: Stable (\(Ri>0.25\)), near-neutral (\(|Ri|\le 0.05\)), unstable (\(Ri<-0.0\)) with configurable thresholds.
• Canonical tables: Error-controlled entries for typical ABL profiles (SBL, CBL, UBL), with parameter sweeps over \(\Delta z\), roughness, and curvature proxies.
• Visualization: Phase plots of \((Ri_b, Ri_g)\), bias ratio contours, and transition ladders across height.

---

Repository structure

richardson-curvature-toolkit/
├─ .github/
│  └─ workflows/ci.yml
├─ src/rct/
│  ├─ core/
│  │  ├─ derivatives.py
│  │  ├─ ri_estimators.py
│  │  ├─ curvature.py
│  │  └─ correction_ode.py
│  ├─ diagnostics/
│  │  ├─ bias.py
│  │  ├─ stability.py
│  │  └─ bootstrap.py
│  ├─ data/
│  │  ├─ loaders.py
│  │  └─ schemas.py
│  ├─ viz/
│  │  ├─ tables.py
│  │  └─ plots.py
│  ├─ ml/
│  │  ├─ surrogates.py
│  │  └─ constraints.py
│  └─ utils/
│     ├─ config.py
│     └─ logging.py
├─ notebooks/
│  ├─ 01_intro_richardson.ipynb
│  ├─ 02_curvature_corrections.ipynb
│  ├─ 03_bias_diagnostics.ipynb
│  ├─ 04_reference_tables.ipynb
│  └─ 05_validation_uah.ipynb
├─ tests/
│  ├─ test_derivatives.py
│  ├─ test_ri_estimators.py
│  ├─ test_bias.py
│  ├─ test_correction_ode.py
│  └─ test_viz.py
├─ examples/
│  ├─ quickstart.py
│  └─ pipeline_cli.py
├─ pyproject.toml
├─ rct_config.yaml
├─ README.md
└─ LICENSE

---

APIs, interfaces, and data

Core Python APIs

• Derivatives• Function: central_with_curvature(phi_minus, phi_0, phi_plus, dz)
• Returns: First derivative estimate and curvature proxy.
• Notes: Error-controlled second-order; optional TVD limiter.

• Curvature• Function: curvature_proxy(theta_triplet, U_triplet, V_triplet, dz)
• Returns: Combined curvature magnitude for ODE driving.

• Estimators• Function: ri_gradient(theta_triplet, U_triplet, V_triplet, dz, theta0)
• Function: ri_bulk(theta1, theta2, U1, U2, V1, V2, z1, z2, theta0)
• Function: ri_gradient_corrected(..., params)
• Returns: Ri estimate plus uncertainty metadata.

• Diagnostics• Function: bias_ratio(ri_g, ri_b)
• Function: bootstrap_bias(profiles, n=1000, seed=None)
• Function: stability_regime(ri, thresholds)

• Correction ODE• Class: CorrectionODE(alpha, beta, integrator="rk23")
• Method: solve(z_grid, kappa_z, shear_z, C0=1.0)
• Output: C(z) monotone-constrained solution.

• Viz and Tables• Function: make_reference_table(params_grid, regimes)
• Function: plot_bias_phase(ri_b_array, ri_g_array)
• Function: plot_C_profile(z, C)

CLI interface

• Entry point: rct (via console_scripts)
• Commands:• rct compute: Load profile CSV/NetCDF, output Ri estimates + corrections.
• rct diagnose: Generate bias reports and stability regime summaries.
• rct table: Build reference tables over parameter grids.
• rct validate: Run notebook-driven validations.

Data expectations and schema

• Input schema (columnar):• Required: z, theta, U, V
• Optional: theta0, metadata (station, time, instrument)

• Units: z in meters, theta in K, U,V in m/s; enforce via validators.
• Formats: CSV for quickstart, NetCDF for operational datasets.
• Config: Global parameters in rct_config.yaml (thresholds, integrator options, bootstrap sizes).

---

Validation, testing, and CI

Unit and property tests

• Derivatives: Second-order convergence on analytic profiles; limiter monotonicity.
• Ri estimators: Consistency between gradient and bulk under affine profiles; boundedness under noise.
• Bias diagnostics: Bootstrap confidence intervals, regime classification stability.
• Correction ODE: Positivity/monotonicity, parameter sensitivity, integrator equivalence.

Synthetic and real-data validation

• Synthetic profiles:• Label: Affine, quadratic, cubic temperature and wind profiles with known Ri.
• Runs: Parameter sweeps over \(\Delta z\), noise levels, curvature.

• Real datasets:• Label: Tower/Lidar/Radio soundings; initial validation notebooks structured to plug UAH datasets.

• Metrics:• Bias: Mean relative error between corrected and baseline Ri.
• Uncertainty: Bootstrap CI width, regime classification accuracy.
• Robustness: Performance across instrument resolutions.

CI/CD

• CI workflow: Lint (ruff), type check (mypy), tests (pytest + coverage), notebook execution (nbval/pytest‑nb), packaging (build).
• Badges: Coverage, CI status, PyPI readiness (optional).
• Releases: Semantic versioning with changelog; GitHub Actions for tagged releases.

---

Documentation and pedagogy

• README.md: Problem overview, quickstart, figures that tell the bias story.
• Docs site (optional mkdocs):• Concepts: Richardson numbers, curvature, discretization bias, stability regimes.
• Math pages: Step‑by‑step derivations with visual analogies and error bounds.
• API reference: Autogenerated from docstrings.
• Tutorials: Notebook series integrated with datasets and exercises.

• Teaching assets:• Exercises: Calibrate corrections on synthetic vs. real profiles.
• Visual analogies: “Curvature as lens distortion” for gradients; phase portraits for bias ratios.

---

Development roadmap

Milestone 0 — bootstrap (2–3 weeks)

• Setup: Repo, packaging, CI, config, baseline notebooks.
• Core math: Curvature-aware gradient and bulk estimators.
• Validation: Synthetic profiles; initial plots.

Milestone 1 — diagnostics and ODE (3–5 weeks)

• Bias modules: Ratios, bootstrap CIs, regime logic.
• Correction ODE: Implement, test monotonicity, parameter sweeps.
• Reference tables: Generate canonical entries, save as CSV/PNG.

Milestone 2 — CLI, docs, pedagogy (3–4 weeks)

• CLI: compute, diagnose, table.
• Docs: Tutorials, API pages, figures.
• Examples: quickstart.py, classroom-ready notebooks.

Milestone 3 — data integration and UAH validation (4–6 weeks)

• Loaders: NetCDF/CSV pipelines, unit enforcement.
• UAH validation: Plug tower/LiDAR datasets; comparative reports.
• Release: v0.1 with changelog, DOI via Zenodo (optional).

Milestone 4 — ML constraints and surrogates (optional, 4–8 weeks)

• Surrogates: PINN/physics‑aware regressors for \(C(z)\) or direct Ri mapping.
• Constraints: Enforce physical bounds and regime consistency.
• Benchmarking: Cross‑validation against analytic baselines.

---

Example usage

from rct.core.ri_estimators import ri_gradient_corrected, ri_bulk
from rct.core.curvature import curvature_proxy
from rct.core.correction_ode import CorrectionODE

# Triplets around z0 with uniform dz

theta_minus, theta_0, theta_plus = 300.2, 300.5, 300.9
U_minus, U_0, U_plus = 2.1, 2.6, 3.4
V_minus, V_0, V_plus = 0.4, 0.6, 1.0
dz = 5.0
theta0 = 300.0

ri_g_corr = ri_gradient_corrected(
    theta_triplet=(theta_minus, theta_0, theta_plus),
    U_triplet=(U_minus, U_0, U_plus),
    V_triplet=(V_minus, V_0, V_plus),
    dz=dz,
    theta0=theta0,
    params={"limiter": "tvd", "curvature_weight": 1.0}
)

ri_b = ri_bulk(
    theta1=300.2, theta2=300.9,
    U1=2.1, U2=3.4,
    V1=0.4, V2=1.0,
    z1=0.0, z2=10.0,
    theta0=theta0
)

kappa = curvature_proxy(
    (theta_minus, theta_0, theta_plus),
    (U_minus, U_0, U_plus),
    (V_minus, V_0, V_plus),
    dz
)

ode = CorrectionODE(alpha=0.3, beta=0.7, integrator="rk23")
z_grid = [0, 5, 10, 15, 20]
shear_z = [0.5, 0.6, 0.9, 0.8, 0.7]
C = ode.solve(z_grid, kappa_z=[kappa]*len(z_grid), shear_z=shear_z, C0=1.0)

ri_g_final = C[1] * ri_g_corr  # apply at the corresponding z index

---

GitHub and VS Code setup tips

• Dev environment: uv or pip-tools for locked deps; VS Code tasks for lint/test/notebooks; .vscode/settings.json for black/ruff/mypy integration.
• Pre-commit hooks: Format, lint, static types, notebook strip outputs.
• Issue templates: “Bug”, “Feature”, “Validation dataset” with checklists.
• Branching: main (stable), dev (integration), feature branches per module.
• Visibility: README figures, badges, minimal examples; use Discussions for pedagogy and validation threads.

---

If you want, I can generate the initial repo scaffolding (files and stubs) and a first notebook outlining the curvature derivations and bias plots—you’ll have a working prototype in under an hour.
