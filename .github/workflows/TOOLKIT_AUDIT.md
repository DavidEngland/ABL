# Richardson Curvature Toolkit — Implementation Audit

**Date:** December 16, 2025  
**Status:** Planning Phase (Milestone 0)

---

## Summary

The toolkit specifications in [toolkit/README.md](.github/workflows/toolkit/README.md) define a comprehensive curvature-aware Richardson number diagnostic and correction framework. This audit maps what exists across the ABL workspace against the planned repository structure and identifies implementation priorities.

**TL;DR:**
- ✅ **Theory documented:** ~15 papers/notes on curvature, bias, dynamic Ri_c
- ✅ **Profile generation:** `profiles.py` implements 8 MOST families + Newton inversion
- ✅ **Derivatives & formulas:** Scattered across `.md` files, not yet unified in Python
- ⚠️ **Core modules:** Missing organized `src/rct/` structure with unit tests
- ⚠️ **Diagnostics:** Bias calculation exists piecemeal; no bootstrap, regime logic yet
- ⚠️ **Correction ODE:** Pseudocode in docs; not yet implemented
- ⚠️ **CLI:** No entry points yet
- ⚠️ **Notebooks:** Web tools exist; validation notebooks not yet written

---

## Existing Assets by Toolkit Component

### 1. **Core Math & Derivatives**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `Ri_g` definition & formula | ✅ | `Ri2.md`, `Ri.md`, `scratch.md`, `web/most.js` | Documented; JS impls exist |
| `Ri_b` definition | ✅ | `Ri2.md`, `HW_Jensen_Geometric_Mean.tex` | Well-documented |
| Curvature formula (∂²Ri_g/∂ζ²) | ✅ | `Curvature of the Gradient Richardson Number.tex`, `Physical Interpretation...tex` | Complete derivation |
| Log-derivative lemma | ✅ | `code/log-derivatives lemma.md` | Model-agnostic compact form |
| Central difference (O(h²)) | ⚠️ | `code/profiles.py` (helper `_central_diff`) | Implemented but not exposed in core API |
| TVD limiter option | ❌ | Mentioned in spec; not found | Needs implementation |

**Action:** Unify derivative code into `src/rct/core/derivatives.py` with public API.

---

### 2. **Profile Generation & MOST Families**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| BD_PL (Businger–Dyer) | ✅ | `code/profiles.py` | Complete |
| BD_CLASSIC (piecewise) | ✅ | `code/profiles.py` | Complete |
| HOG88 (Högström linear) | ✅ | `code/profiles.py` | Complete |
| QSBL (quadratic truncation) | ✅ | `code/profiles.py` | Complete |
| CB (Cheng–Brutsaert) | ✅ | `code/profiles.py` | Complete |
| RPL (regularized power-law) | ✅ | `code/profiles.py` | Complete |
| VEXP (variable exponent) | ✅ | `code/profiles.py` | Complete |
| DTP (dynamic Prandtl) | ✅ | `code/profiles.py` | Complete |
| URC (Ri-based direct) | ✅ | `code/profiles.py` | Complete |
| Canonical IDs (BD71, BH91, CB05, GF97) | ⚠️ | `code/profiles.py` line ~500+ | Exists; needs verification & docs |
| Newton inversion (ζ from Ri_g) | ✅ | `code/profiles.py` | Complete |
| Series seed (ζ₀) | ✅ | `code/profiles.py` | Complete |

**Action:** Move `profiles.py` → `src/rct/core/profiles.py`; document canonical ID registry.

---

### 3. **Richardson Number Estimators**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `ri_gradient(theta_triplet, ...)` | ⚠️ | `web/most.js` (JS); scattered in `.py` | Not unified Python API |
| `ri_bulk(theta1, theta2, ...)` | ⚠️ | `web/most.js` (JS); scattered in `.py` | Not unified Python API |
| `ri_gradient_corrected(...)` | ⚠️ | `code/profiles.py` partial | Needs explicit correction factor integration |
| Error bounds / uncertainty metadata | ❌ | Not found | Spec calls for metadata; not implemented |

**Action:** Create `src/rct/core/ri_estimators.py` with NumPy-vectorized implementations.

---

### 4. **Curvature & Bias Diagnostics**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `curvature_proxy(...)` | ⚠️ | `web/most.js` partial; `code/` scattered | Not unified; needs z-domain + ζ-domain variants |
| `bias_ratio(ri_g, ri_b)` | ⚠️ | `code/profiles.py` line ~300+ | Exists; needs extraction to core module |
| `bootstrap_bias(profiles, n=1000)` | ❌ | Not found | High priority; missing |
| `stability_regime(ri, thresholds)` | ⚠️ | Thresholds in docs; not coded | Needs implementation |
| Neutral curvature (2Δ) invariant | ✅ | `code/profiles.py`, docs | Well-documented |
| Geometric mean height (z_g) | ✅ | `HW_Jensen_Geometric_Mean.tex`, `CANONICAL_GLOSSARY.md` | Well-documented |

**Action:** Create `src/rct/diagnostics/bias.py` and `stability.py`; implement bootstrap.

---

### 5. **Correction ODE**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| ODE formulation (dC/dz = ...) | ✅ | Toolkit spec, McNider papers | Theory complete |
| Pseudocode | ✅ | `Curvature-Aware...Grid-Dependent...md` | Detailed walkthrough |
| Monotonicity constraints | ⚠️ | Mentioned; not coded | Needs implementation |
| Explicit/implicit integrators | ❌ | Not found | Needs RK23, RK45 wrappers |
| Class `CorrectionODE` with `.solve()` | ❌ | Not found | Priority implementation |

**Action:** Create `src/rct/core/correction_ode.py` with scipy integrators.

---

### 6. **Reference Tables & Visualization**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `make_reference_table(params_grid, ...)` | ❌ | Not found | High value; needs implementation |
| `plot_bias_phase(ri_b, ri_g)` | ⚠️ | `web/most2.html`, `web/curvature_standalone.html` | Web-based; needs Python/matplotlib version |
| `plot_C_profile(z, C)` | ⚠️ | Web tools partial | Needs Python version |
| Phase plots (Ri_b, Ri_g) | ✅ | Embedded in web tools | Visual implemented; needs extraction |
| Bias ratio contours | ✅ | Embedded in web tools | Visual implemented; needs extraction |
| Transition ladders | ⚠️ | Conceptual; not visualized | Needs implementation |

**Action:** Create `src/rct/viz/tables.py` and `plots.py` using matplotlib.

---

### 7. **Data Handling & Schema**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| Column schema (z, theta, U, V) | ✅ | `config/schema_complete.sql`, `data/` | Well-defined |
| Unit enforcement (validators) | ⚠️ | `config/` partial | Needs Python pydantic models |
| CSV loader | ⚠️ | `implementations/data_pipeline.py` | Exists; needs integration |
| NetCDF loader | ❌ | Not found | Priority; add xarray wrapper |
| Config system (rct_config.yaml) | ❌ | Not found | Needs creation |

**Action:** Create `src/rct/data/loaders.py` and `schemas.py`; add `rct_config.yaml`.

---

### 8. **CLI Interface**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `rct compute` | ❌ | Not found | Needs Click/argparse entry point |
| `rct diagnose` | ❌ | Not found | Needs Click/argparse entry point |
| `rct table` | ❌ | Not found | Needs Click/argparse entry point |
| `rct validate` | ❌ | Not found | Needs Click/argparse entry point |
| Console scripts in pyproject.toml | ❌ | Not found | Needs setup |

**Action:** Create `src/rct/cli.py` with Click; configure pyproject.toml.

---

### 9. **Testing**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `test_derivatives.py` | ❌ | Not found | Convergence tests needed |
| `test_ri_estimators.py` | ❌ | Not found | Consistency, boundedness tests |
| `test_bias.py` | ❌ | Not found | Bootstrap CI, regime logic |
| `test_correction_ode.py` | ❌ | Not found | Monotonicity, integrator tests |
| `test_viz.py` | ❌ | Not found | Plot generation, no regression |
| Synthetic profile suite | ⚠️ | `implementations/toy_sc_m.py` | Partial; needs expansion |
| CI/CD (pytest, coverage, mypy) | ⚠️ | `.github/workflows/` | Exists; needs RCT-specific config |

**Action:** Create `tests/` directory with full suite; integrate pytest + coverage.

---

### 10. **Documentation & Pedagogy**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `README.md` (main) | ⚠️ | `/ABL/README.md` | Exists; needs RCT-specific version |
| Concept pages (math tutorials) | ⚠️ | Scattered across `.md`, `.tex` | High-quality content; needs organization |
| API reference (auto-docstring) | ❌ | Not yet; will use Sphinx | Needs setup |
| Jupyter tutorials | ⚠️ | `notebooks/` exist; not RCT-focused | Need 5 main notebooks per spec |
| Teaching exercises | ⚠️ | `hw/` folder; not integrated | Needs adaptation |
| Visual analogies | ✅ | Embedded in web tools | Good; needs LaTeX rendering fix |

**Action:** Create `docs/` with mkdocs structure; write 5-notebook tutorial series.

---

### 11. **Packaging & Deployment**

| Spec Item | Status | Location(s) | Notes |
|-----------|--------|-------------|-------|
| `pyproject.toml` | ⚠️ | Root; exists but not RCT-configured | Needs `[tool.poetry]`, console scripts |
| `setup.py` | ⚠️ | Not found; recommend pyproject only | Modern approach: pyproject + Poetry |
| Semantic versioning | ❌ | Not tracked; needs CHANGELOG | Create v0.0.1-alpha |
| Zenodo/DOI integration | ❌ | Not found | Optional; defer to v0.2+ |
| PyPI readiness | ❌ | Not yet; needs build system | Defer to v0.1 release |

**Action:** Update pyproject.toml; create CHANGELOG.md for RCT.

---

### 12. **Related Existing Work**

| Asset | Type | Value | Link(s) |
|-------|------|-------|--------|
| GABLS validation framework | 📚 | High | `data/GABLS.md` |
| ML integration strategy | 📚 | Medium | `ML/ML_Curvature_Corrections_Guide.md` |
| Dynamic Ri_c formulation | 📚 | High | `Dynamic Critical Richardson Number.md` |
| McNider reconciliation | 📚 | High | `McNider_Ri_Corrections_Overview.md` |
| Bias teaching homework | 📚 | Medium | `HW_Jensen_Geometric_Mean.tex` |
| Web curvature tools | 🔧 | Medium | `web/curvature_standalone.html`, `most2.html` |
| Profile catalog (Biazar) | 📚 | High | `code/profiles.py` |
| Database schema | 🔧 | Medium | `config/schema_complete.sql` |

---

## Implementation Roadmap

### Milestone 0 — **Bootstrap & Core Math** (2–3 weeks)

**Deliverables:**
1. ✅ Repo structure: `src/rct/core/`, `tests/`, `notebooks/`, `docs/`
2. ✅ `src/rct/core/derivatives.py` — unified derivative & curvature API
3. ✅ `src/rct/core/ri_estimators.py` — `ri_gradient()`, `ri_bulk()`, metadata
4. ✅ `src/rct/core/profiles.py` — move & refactor from `code/profiles.py`
5. ✅ `src/rct/core/__init__.py` — version, imports
6. ✅ `tests/test_derivatives.py` — convergence tests on synthetic profiles
7. ✅ `pyproject.toml` — modern packaging
8. ✅ `notebooks/01_intro_richardson.ipynb` — theory + definitions
9. ✅ `rct_config.yaml` — thresholds, integrator settings

**Success Metrics:**
- `python -c "from rct.core import ri_gradient"` works
- `pytest tests/test_derivatives.py -v` passes (O(h²) convergence proven)
- Notebook 01 runs without errors

---

### Milestone 1 — **Diagnostics & ODE** (3–5 weeks)

**Deliverables:**
1. ✅ `src/rct/diagnostics/bias.py` — `bias_ratio()`, confidence intervals
2. ✅ `src/rct/diagnostics/stability.py` — regime classification
3. ✅ `src/rct/diagnostics/bootstrap.py` — `bootstrap_bias()` with seed control
4. ✅ `src/rct/core/correction_ode.py` — `CorrectionODE` class, monotonicity guards
5. ✅ `tests/test_bias.py`, `test_correction_ode.py` — property tests
6. ✅ `notebooks/02_curvature_corrections.ipynb` — ODE solver tutorial
7. ✅ `notebooks/03_bias_diagnostics.ipynb` — bootstrap & regimes
8. ✅ Reference table generation (`src/rct/viz/tables.py` basic version)

**Success Metrics:**
- Bootstrap CI width ≤ 5% Ri for typical profiles
- `CorrectionODE.solve()` monotone + energy-conserving
- Notebooks 02, 03 executable end-to-end

---

### Milestone 2 — **CLI & Docs** (2–3 weeks)

**Deliverables:**
1. ✅ `src/rct/cli.py` — `click` entry points (`compute`, `diagnose`, `table`, `validate`)
2. ✅ `src/rct/data/loaders.py` — CSV, NetCDF I/O
3. ✅ `src/rct/viz/plots.py` — matplotlib wrappers for bias phase plots
4. ✅ `docs/` structure — mkdocs or Sphinx
5. ✅ `examples/quickstart.py` — minimal working example
6. ✅ `CHANGELOG.md` (RCT-specific)
7. ✅ API reference auto-generated

**Success Metrics:**
- `rct compute sample.csv --output results.json` runs
- `rct diagnose --help` displays all subcommands
- Docs deploy locally without warnings

---

### Milestone 3 — **Data Integration & UAH Validation** (4–6 weeks)

**Deliverables:**
1. ✅ ARM NSA + GABLS LES loaders (CSV templates provided)
2. ✅ `notebooks/04_reference_tables.ipynb` — generate canonical tables
3. ✅ `notebooks/05_validation_uah.ipynb` — plug-and-play for tower data
4. ✅ Validation report templates (HTML + PDF)
5. ✅ v0.1 release with DOI

**Success Metrics:**
- Reference tables saved to `data/canonical_tables.csv`
- UAH validation notebook produces bias corrections for real tower data
- GitHub release + Zenodo DOI

---

### Milestone 4 — **ML Constraints & Surrogates** (Optional, 4–8 weeks)

**Deliverables:**
1. ⏳ `src/rct/ml/surrogates.py` — PINN/physics-aware regressors
2. ⏳ `src/rct/ml/constraints.py` — bound enforcement, regime consistency
3. ⏳ Benchmarking notebook
4. ⏳ v0.2 release

---

## Recommended Priority Sequence

**Phase 1 (Week 1–2):** *Organize & Stabilize*
- [ ] Create `richardson-curvature-toolkit/` directory structure
- [ ] Move `code/profiles.py` → `src/rct/core/profiles.py`
- [ ] Extract & unify derivative helpers into `core/derivatives.py`
- [ ] Write `pyproject.toml`, `rct_config.yaml`
- [ ] Set up pytest, pre-commit hooks

**Phase 2 (Week 2–3):** *Core APIs*
- [ ] `ri_estimators.py` (vectorized)
- [ ] `curvature.py` (unified z & ζ variants)
- [ ] `test_derivatives.py` (convergence proofs)
- [ ] First Jupyter notebook (intro + plots)

**Phase 3 (Week 3–4):** *Diagnostics*
- [ ] `bias.py`, `bootstrap.py`, `stability.py`
- [ ] ODE solver skeleton (choose scipy integrator)
- [ ] Notebooks 02–03

**Phase 4 (Week 4–5):** *Deployment*
- [ ] CLI entry points
- [ ] Data loaders (CSV, NetCDF)
- [ ] Examples & quickstart
- [ ] Basic docs

**Phase 5 (Week 5–6+):** *Validation & Release*
- [ ] UAH tower data integration
- [ ] Reference table generation
- [ ] Release v0.1 with changelog

---

## Next Actions (Immediate)

1. **Clarify:** Do you want:
   - A. Standalone `richardson-curvature-toolkit/` repo (cleaner; separate from ABL)?
   - B. Integrate into existing `/ABL` workspace (more coupled; easier testing)?

2. **Choose integrator:** Scipy (RK23/RK45)? JAX (autodiff + performance)? Custom (pedagogical)?

3. **Assign ownership:** 
   - Theory & math: McNider/Biazar review?
   - Code implementation: Your lead?
   - Validation data: UAH collaboration?

4. **Select first win:** Which Milestone 0 component to code first?

---

## Notes

- All `.md` and `.tex` files in `/ABL` are high-value **source material**; keep as `docs/references/`.
- `web/most.js` has solid D3.js visualizations; can be ported to Python for consistency.
- `profiles.py` is **production-ready**; just needs API polish & tests.
- MathJax rendering issue in web tools is a **separate, lower-priority** fix (use unicode fallback for now).

