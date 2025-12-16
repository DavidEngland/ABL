# Richardson Curvature Toolkit — Implementation Guide

**Version:** 0.0.1-alpha  
**Status:** Milestone 0 (Bootstrap Phase)  
**Created:** December 16, 2025  
**Target Completion:** Week of December 23–30, 2025

---

## Quick Start for Developers

### Clone & Setup

```bash
cd /Users/davidengland/Documents/GitHub/ABL
python -m venv .venv
source .venv/bin/activate
pip install numpy scipy matplotlib pydantic click pytest
pip install -e .  # Install ABL/src in development mode
```

### Module Structure

```
src/rct/
├── __init__.py                # Package initialization, version
├── core/                      # Core Richardson number algorithms
│   ├── __init__.py
│   ├── derivatives.py         # ✅ DONE: Central differences, curvature
│   ├── ri_estimators.py       # ✅ DONE: Ri_g, Ri_b, bias_ratio
│   ├── profiles.py            # ⏳ TODO: Extract from code/profiles.py
│   ├── curvature.py           # ⏳ TODO: Curvature proxy, phase analysis
│   └── correction_ode.py      # ⏳ TODO: Correction factor ODE solver
├── diagnostics/               # Analysis & characterization
│   ├── __init__.py
│   ├── bias.py                # ⏳ TODO: Bias diagnostics, z_g
│   ├── stability.py           # ⏳ TODO: Regime classification
│   └── bootstrap.py           # ⏳ TODO: Bootstrap confidence intervals
├── data/                      # Data I/O & validation
│   ├── __init__.py
│   ├── loaders.py             # ⏳ TODO: CSV, NetCDF loaders
│   └── schemas.py             # ⏳ TODO: Pydantic models
├── viz/                       # Visualization & tables
│   ├── __init__.py
│   ├── plots.py               # ⏳ TODO: Matplotlib wrappers
│   └── tables.py              # ⏳ TODO: Reference table generation
├── utils/                     # Support utilities
│   ├── __init__.py
│   ├── config.py              # ⏳ TODO: YAML config loading
│   └── logging.py             # ⏳ TODO: Logging setup
└── ml/                        # Machine learning (Milestone 4)
    └── __init__.py
```

---

## Completed Components

### 1. `core/derivatives.py` — ✅ DONE

**Status:** Production-ready with tests

**What it does:**
- `central_with_curvature(phi_minus, phi_0, phi_plus, dz, limiter="tvd")` — 2nd-order central difference with optional TVD limiter
- `second_derivative(phi_minus, phi_0, phi_plus, dz)` — Curvature estimation
- `numerical_jacobian(func, x, h=1e-5)` — Jacobian matrix via finite differences
- `richardson_extrapolation(func, x, h0, order, n_steps)` — High-precision derivative via Romberg
- Vectorized: `gradient_array(f, z)`, `curvature_array(f, z)` — Full profile operations

**Usage example:**
```python
from rct.core.derivatives import central_with_curvature

# Three-point stencil: θ at z-Δz, z, z+Δz
theta = [300.2, 300.5, 300.9]
dz = 5.0

first_deriv, curvature = central_with_curvature(theta[0], theta[1], theta[2], dz)
print(f"dθ/dz ≈ {first_deriv:.4f} K/m")
print(f"d²θ/dz² ≈ {curvature:.6f} K/m²")
```

**Tests needed:** `tests/test_derivatives.py`
- Convergence on polynomial functions (O(h²) accuracy)
- TVD limiter monotonicity
- Richardson extrapolation error reduction

---

### 2. `core/ri_estimators.py` — ✅ DONE

**Status:** Production-ready with docstrings

**What it does:**
- `ri_gradient(theta_triplet, U_triplet, V_triplet, dz, theta0)` — Point-wise gradient Richardson
- `ri_bulk(theta1, theta2, U1, U2, V1, V2, z1, z2, theta0)` — Layer-averaged Richardson
- `bias_ratio(ri_g, ri_b, min_denominator=1e-6)` — Bias ratio B = Ri_g / Ri_b
- `estimate_uncertainty(ri_estimate, derivative_step, noise_level)` — Error metadata
- Vectorized: `ri_gradient_array(theta, U, V, z, theta0, g)` — Full profile

**Usage example:**
```python
from rct.core.ri_estimators import ri_gradient, ri_bulk, bias_ratio

# Gradient Richardson at z=10m (triplet from model/tower)
rig = ri_gradient(
    theta_minus=300.0, theta_0=300.2, theta_plus=300.5,
    U_minus=2.0, U_0=2.5, U_plus=3.0,
    V_minus=0.1, V_0=0.2, V_plus=0.3,
    dz=10.0, theta0=300.0
)
print(f"Ri_g(z=10m) = {rig:.4f}")

# Bulk Richardson over layer 10-100m
rib = ri_bulk(
    theta1=300.0, theta2=302.0,
    U1=2.0, U2=5.0,
    V1=0.1, V2=0.5,
    z1=10.0, z2=100.0,
    theta0=300.0
)
print(f"Ri_b(10-100m) = {rib:.4f}")

# Bias diagnostic
B = bias_ratio(rig, rib)
print(f"Bias ratio B = {B:.2f}")  # B > 1: coarse grid underestimates
```

**Tests needed:** `tests/test_ri_estimators.py`
- Consistency: Ri_g ≈ Ri_b on uniform profiles
- Boundedness: No NaN/inf on edge cases (calm wind, large noise)
- Physical: Unstable, neutral, stable regimes reproduce expected signs

---

### 3. `rct_config.yaml` — ✅ DONE

**Status:** Default configuration ready for customization

**Key sections:**
- `stability_thresholds`: Ri values for regime classification
- `derivatives`: Numerical step size, limiter choice
- `correction_ode`: ODE solver parameters (α, β, integrator)
- `bootstrap`: Bootstrap CI settings (n_samples, confidence_level)
- `visualization`: Matplotlib defaults (DPI, colormap)
- `reference_tables`: Parameter grid for canonical tables
- `data`: Units, required/optional columns, validation

**How to use:**
```python
import yaml

with open('rct_config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Access threshold for critical Ri
Ri_c = config['stability_thresholds']['supercritical']
```

---

## In-Progress / TODO Components

### Priority 1 — **Core APIs** (Week 1-2)

#### `core/profiles.py` — ⏳ EXTRACT FROM `code/profiles.py`

**Work:**
- Move 589-line `code/profiles.py` into `src/rct/core/profiles.py`
- Ensure imports work (`numpy`, `math`)
- Public API:
  - `make_profile(tag: str, pars: dict) → (phi_m, phi_h)` callables
  - `zeta_from_ri_gradient(ri_target, profile_tag, pars) → float` (Newton inversion)
  - `profile_catalog() → dict` with canonical IDs: BD71, BH91, CB05, GF97
- Add docstring with example usage

**Estimated time:** 1 hour (move + polish)  
**Test:** Unit tests in `tests/test_profiles.py` (convergence of Newton solver)

---

#### `core/curvature.py` — ⏳ IMPLEMENT

**Work:**
- Implement `curvature_proxy(theta_triplet, U_triplet, V_triplet, dz) → float`
  - Combine curvature from θ, U, V for ODE driving
  - Normalize by 2Δ (neutral invariant) for regime-independent scaling
- Implement `compute_neutral_curvature(phi_m_params, phi_h_params) → float`
  - Use log-derivative lemma (see `code/log-derivatives lemma.md`)
  - Return 2Δ = α_h β_h - 2 α_m β_m
- Import helpers from `derivatives.py` and `profiles.py`

**Estimated time:** 2–3 hours  
**Reference:** `code/log-derivatives lemma.md`, `Curvature of the Gradient Richardson Number.tex`

---

#### `core/correction_ode.py` — ⏳ IMPLEMENT

**Work:**
- Implement `CorrectionODE` class:
  - `__init__(alpha, beta, integrator="rk45")` — parameters
  - `solve(z_grid, kappa_z, shear_z, C0=1.0) → C_array` — ODE integration
- ODE to solve: dC/dz = α κ(z) C - β σ_shear(z) (C - 1)
- Use scipy.integrate.solve_ivp with RK45/RK23
- Apply monotonicity constraint: 0.5 ≤ C(z) ≤ 2.0

**Estimated time:** 3–4 hours  
**Reference:** Toolkit spec, `McNider_Ri_Corrections_Overview.md`

---

### Priority 2 — **Diagnostics & Analysis** (Week 2-3)

#### `diagnostics/bias.py` — ⏳ IMPLEMENT

**Work:**
- `bias_ratio_diagnostic(ri_g, ri_b) → dict` — metadata wrapper
  - Return: B, B-1, confidence interval placeholder
- `geometric_mean_height(z1, z2) → float` — $z_g = \sqrt{z_1 z_2}$
- Integration with bootstrap (next item)

**Estimated time:** 1 hour  
**Reference:** `HW_Jensen_Geometric_Mean.tex`, `CANONICAL_GLOSSARY.md`

---

#### `diagnostics/stability.py` — ⏳ IMPLEMENT

**Work:**
- `classify_regime(ri, thresholds=None) → str` — return "unstable" | "neutral" | "stable" | "supercritical"
- `regime_thresholds(config=None) → dict` — load from config or default
- Use thresholds from `rct_config.yaml`

**Estimated time:** 0.5 hour

---

#### `diagnostics/bootstrap.py` — ⏳ IMPLEMENT

**Work:**
- `bootstrap_bias(profiles, n_samples=1000, seed=42, method="percentile") → dict`
  - Input: List of (z, theta, U, V) tuples
  - Resample with replacement, compute bias for each sample
  - Return: CI bounds, percentiles, error bars
- Use `random.seed(seed)` for reproducibility

**Estimated time:** 2 hours  
**Reference:** `rct_config.yaml` bootstrap settings

---

### Priority 3 — **Data & Visualization** (Week 3-4)

#### `data/loaders.py` & `schemas.py` — ⏳ IMPLEMENT

**Work:**
- `load_csv(filepath) → dict` with columns: z, theta, U, V, optional V
- `load_netcdf(filepath, vars=None) → dict` via xarray
- Pydantic `ProfileData` class for validation
- Unit conversion helpers

**Estimated time:** 2–3 hours

---

#### `viz/plots.py` & `tables.py` — ⏳ IMPLEMENT

**Work:**
- `plot_bias_phase(ri_b_array, ri_g_array)` — 2D scatter, regime coloring
- `plot_ri_profile(z, ri_g, ri_b=None)` — line plot with z axis
- `plot_curvature(z, curvature)` — d²Ri_g/dζ² profile
- `make_reference_table(params_grid)` — parameter sweep CSV output

**Estimated time:** 3–4 hours

---

### Priority 4 — **CLI & Testing** (Week 4-5)

#### `cli.py` — ⏳ IMPLEMENT

**Work:**
- Use Click framework
- Subcommands: `rct compute`, `rct diagnose`, `rct table`, `rct validate`
- Wire into `pyproject.toml` console_scripts

**Estimated time:** 2–3 hours

---

#### `tests/` — ⏳ IMPLEMENT

**Work:**
- `test_derivatives.py` — Convergence, limiter monotonicity
- `test_ri_estimators.py` — Consistency, edge cases
- `test_bias.py` — Bootstrap CI coverage
- `test_correction_ode.py` — Monotonicity, energy conservation
- `test_profiles.py` — Newton convergence, canonical IDs
- `test_viz.py` — Plot generation (no regression, just no errors)

**Estimated time:** 8–10 hours total

---

## Development Workflow

### 1. Check Status
```bash
cd /Users/davidengland/Documents/GitHub/ABL
git status
```

### 2. Create Feature Branch
```bash
git checkout -b feature/rct-core-curvature
```

### 3. Implement (with TDD)
```python
# Write test first
def test_curvature_proxy():
    result = curvature_proxy((300.2, 300.5, 300.9), (2.0, 2.5, 3.0), ...)
    assert 0.5 < result < 2.0  # Sanity check

# Implement function
def curvature_proxy(...):
    ...

# Run tests
pytest tests/test_curvature.py -v
```

### 4. Commit & Document
```bash
git add src/rct/core/curvature.py tests/test_curvature.py
git commit -m "Implement curvature_proxy with neutral invariant preservation"
```

### 5. Update Audit
- Move completed items from TODO to DONE in [TOOLKIT_AUDIT.md](.github/workflows/TOOLKIT_AUDIT.md)
- Run: `python -m pytest tests/ --cov=src/rct` and note coverage

---

## Next Actions (Immediate)

1. **This week:**
   - [ ] Extract & refactor `code/profiles.py` → `src/rct/core/profiles.py`
   - [ ] Implement `core/curvature.py` with neutral invariant logic
   - [ ] Implement `core/correction_ode.py` with scipy solver
   - [ ] Write first batch of tests

2. **Validate:**
   - [ ] `python -c "from rct.core import ri_gradient; print(ri_gradient(...))"` works
   - [ ] First Jupyter notebook runs without errors

3. **Document:**
   - [ ] Update main README.md with toolkit overview
   - [ ] Create `docs/` folder structure

---

## References

- **Spec:** [toolkit/README.md](.github/workflows/toolkit/README.md)
- **Audit:** [TOOLKIT_AUDIT.md](.github/workflows/TOOLKIT_AUDIT.md)
- **Theory:** 
  - `Curvature of the Gradient Richardson Number.tex`
  - `code/log-derivatives lemma.md`
  - `McNider_Ri_Corrections_Overview.md`
- **Data/Validation:**
  - `data/GABLS.md` — LES validation cases
  - `config/schema_complete.sql` — Database schema

---

## Contact & Questions

- **Lead:** David England (davidengland@uah.edu)
- **Theory:** Dick McNider, Arastoo Pour-Biazar (UAH)
- **Issues:** Create GitHub issue or Slack message

---

**Last Updated:** December 16, 2025  
**Next Review:** December 23, 2025
