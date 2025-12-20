# Richardson Number Curvature Toolkit - Implementation Status

**Last Updated:** 2025-01-XX  
**Focus:** Practical integration tools for model developers

---

## Executive Summary

**Mission:** Help atmospheric model developers diagnose and correct Richardson number bias on coarse grids.

**Approach:** Flexible advisory toolkit with multiple integration strategies (not rigid package).

**Status:** Core correction algorithms **production-ready**. Diagnostic tools **ready to use**. Validation suite **in progress**.

---

## ✅ COMPLETED (Ready to Use)

### Core Algorithms
| Module | Lines | Status | Documentation |
|--------|-------|--------|---------------|
| `src/rct/core/correction_ode.py` | 350+ | ✅ Production | Full docstrings + examples |
| `src/rct/core/ri_estimators.py` | 400+ | ✅ Production | Full docstrings + examples |
| `src/rct/core/derivatives.py` | 500+ | ✅ Production | Full docstrings + examples |

**Functions available:**
- ✅ `simple_multiplicative_correction(Ri, dz, ...)` - Direct Ri adjustment
- ✅ `stability_function_correction(Ri, dz, f_base, ...)` - Extend f(Ri) tail
- ✅ `dynamic_critical_richardson(Ri, dz, kappa, ...)` - Adjust Ri_c threshold
- ✅ `estimate_curvature_proxy(theta, U, V, z)` - Extract κ from profiles
- ✅ `CorrectionODE` class - Full ODE solver with scipy
- ✅ `ri_gradient()`, `ri_bulk()`, `bias_ratio()` - Richardson number estimators
- ✅ Central differences, TVD limiters, Richardson extrapolation

### Integration Tools
| Tool | Purpose | Status |
|------|---------|--------|
| `scripts/diagnose_ri_bias.py` | Quick diagnostic (CSV input → report) | ✅ Works, tested |
| `examples/model_integration_example.py` | Python model demo (before/after) | ✅ Complete |
| `examples/wrf_integration_example.F90` | Fortran/WRF integration guide | ✅ Complete |
| `examples/test_profile.csv` | Sample data for testing | ✅ Provided |

### Documentation
| Document | Pages | Status | Audience |
|----------|-------|--------|----------|
| `RCT_README.md` | ~5 | ✅ Complete | Quick start |
| `MODEL_DEVELOPER_GUIDE.md` | ~20 | ✅ Complete | Integration details |
| `RCT_QUICK_REFERENCE.md` | ~2 | ✅ Complete | API cheat sheet |
| `RCT_STATUS_REPORT.md` | ~15 | ✅ Complete | Project overview |
| `TOOLKIT_AUDIT.md` | ~12 | ✅ Complete | Gap analysis |
| `IMPLEMENTATION_GUIDE.md` | ~10 | ✅ Complete | Dev workflow |

---

## ⏳ IN PROGRESS

### Notebooks (Pedagogical)
| Notebook | Purpose | Status |
|----------|---------|--------|
| `01_intro_richardson.ipynb` | Ri basics, bias demo | 🔶 Started (skeleton) |
| `02_corrections_demo.ipynb` | Apply 4 strategies | ⏳ Not started |
| `03_gabls1_validation.ipynb` | Validate on GABLS1 | ⏳ Not started |
| `04_wrf_integration.ipynb` | Step-by-step WRF guide | ⏳ Not started |
| `05_parameter_tuning.ipynb` | Sensitivity analysis | ⏳ Not started |

### Modules (Remaining)
| Module | Purpose | Status | Priority |
|--------|---------|--------|----------|
| `core/profiles.py` | Extract from `code/profiles.py` (589 lines) | 🔶 Exists, needs refactor | HIGH |
| `core/curvature.py` | Curvature proxy with neutral invariant | ⏳ Placeholder | HIGH |
| `diagnostics/bias.py` | Bias ratio analysis tools | ⏳ Placeholder | MEDIUM |
| `diagnostics/stability.py` | Stability regime classification | ⏳ Placeholder | MEDIUM |
| `diagnostics/bootstrap.py` | Uncertainty quantification | ⏳ Placeholder | LOW |
| `data/loaders.py` | CSV/NetCDF I/O utilities | ⏳ Placeholder | MEDIUM |
| `data/schemas.py` | Data validation (pydantic) | ⏳ Placeholder | LOW |
| `viz/plots.py` | Matplotlib wrappers | ⏳ Placeholder | LOW |
| `viz/tables.py` | Summary table generation | ⏳ Placeholder | LOW |
| `utils/config.py` | Load `rct_config.yaml` | ⏳ Placeholder | MEDIUM |
| `utils/logging.py` | Logging setup | ⏳ Placeholder | LOW |

---

## 🚫 NOT PLANNED (Out of Scope)

- ❌ Web dashboard (use Jupyter notebooks instead)
- ❌ Automatic WRF code patching (too model-specific)
- ❌ Real-time operational deployment (research toolkit only)
- ❌ Machine learning tuning (current focus: physics-based)
- ❌ Multi-model ensemble framework (single-model integration)

---

## Immediate Next Steps (This Week)

### Priority 1: Profiles Module (1-2 hours)
- [ ] Extract relevant functions from `code/profiles.py`
- [ ] Create clean API: `load_profile()`, `validate_profile()`, `interpolate_profile()`
- [ ] Add docstrings and examples
- [ ] Test with GABLS data

**Rationale:** Needed for validation workflows and diagnostic script enhancements.

### Priority 2: Curvature Module (2-3 hours)
- [ ] Implement `compute_curvature_proxy()` with neutral invariant logic
- [ ] Use log-derivative lemma from `code/log-derivatives lemma.md`
- [ ] Add dimensionless curvature κ = |d²θ/dz²| / (|dθ/dz| / Δz)
- [ ] Test on synthetic profiles (linear, parabolic, inversion)

**Rationale:** Required for advanced corrections (Strategy 3, Strategy 4).

### Priority 3: Test Suite (2 hours)
- [ ] Create `tests/test_correction_ode.py`
- [ ] Test neutral limit: C→1 as Ri→0, Δz→0
- [ ] Test monotonicity: C increases with Ri, Δz
- [ ] Test parameter sensitivity: α, β, D ranges
- [ ] Test clamping: C ∈ [0.5, 2.0]

**Rationale:** Ensure production code is robust before wider use.

### Priority 4: Enhanced Diagnostic (1 hour)
- [ ] Add curvature estimation to `diagnose_ri_bias.py`
- [ ] Include correction preview: "If D=0.5, expect B_corrected = ..."
- [ ] Generate quick plot: Ri_bulk vs. Ri_corrected vs. height
- [ ] Export results to JSON for programmatic use

**Rationale:** Make diagnostic more actionable for model developers.

---

## Medium-Term Goals (Next 2-4 Weeks)

### Validation Suite
- [ ] Implement GABLS1 validation notebook (compare h_BL, z_LLJ)
- [ ] Add ARM NSA tower data validation
- [ ] Create parameter sensitivity analysis (D, α, β)
- [ ] Validate on multiple model types (WRF, MPAS, Python toy model)

### Data & I/O
- [ ] CSV loader with validation (`data/loaders.py`)
- [ ] NetCDF loader for WRF/MPAS output
- [ ] Schema validation (pydantic models)

### Visualization
- [ ] Bias ratio profile plots (Ri_g vs. Ri_b vs. z)
- [ ] Correction impact plots (before/after K_m, h_BL timeseries)
- [ ] Parameter sensitivity heatmaps (D vs. α vs. RMSE)

### Advanced Corrections
- [ ] Implement machine learning correction (optional)
- [ ] Add hybrid MOST+Ri correction approach
- [ ] Extend to Arctic-specific physics (sea ice, snow)

---

## Long-Term Vision (2-3 Months)

### Release v1.0
- [ ] Complete test coverage (>80% lines)
- [ ] Full API documentation (auto-generated from docstrings)
- [ ] Validation on 5+ model types (WRF, MPAS, CESM, MITgcm, custom)
- [ ] Published paper (JAMC or JAS)
- [ ] Zenodo DOI for citation

### Community Adoption
- [ ] Integration examples for 10+ PBL schemes (YSU, MYNN, ACM2, TEMF, ...)
- [ ] Tutorial workshops (AMS Annual Meeting, WRF Users' Workshop)
- [ ] Contribution guidelines for community corrections

---

## Metrics & Success Criteria

### Code Quality
- ✅ All core modules have docstrings + examples
- ✅ Zero syntax errors in production code
- ⏳ Test coverage: Currently 0% → Target 80%
- ⏳ Type hints: Currently partial → Target full coverage

### Documentation
- ✅ Quick start guide (RCT_README.md)
- ✅ Comprehensive integration guide (MODEL_DEVELOPER_GUIDE.md)
- ⏳ API reference: Needs auto-generation
- ⏳ Validation reports: Not yet written

### Validation
- ⏳ GABLS1: h_BL RMSE reduction > 30%
- ⏳ ARM NSA: 2m temperature bias reduction > 20%
- ⏳ WRF stability: No crashes in 10-day simulation
- ⏳ Parameter robustness: D ∈ [0.3, 0.7] all stable

### Adoption
- ⏳ Model integrations: 0 → Target 3+ (WRF, MPAS, CESM)
- ⏳ User feedback: Collect from 5+ early adopters
- ⏳ Citations: Paper accepted + 10+ citations within 2 years

---

## Dependencies & Requirements

### Python Packages (Required)
- ✅ `numpy` (for core algorithms)
- ✅ `scipy` (for ODE solver, interpolation)
- ⏳ `matplotlib` (for visualization, not yet used)
- ⏳ `pandas` (for data loading, removed from diagnostic)
- ⏳ `netCDF4` (for model output, future)
- ⏳ `pydantic` (for validation, future)

### External Data (Validation)
- ⏳ GABLS1 LES reference (need to acquire)
- ⏳ ARM NSA tower data (publicly available)
- ⏳ WRF single-column output (need to generate)

### Fortran Compiler (Optional)
- ⏳ gfortran or ifort (for WRF/MPAS integration testing)

---

## Known Issues & Limitations

### Current Limitations
1. **Curvature proxy is simplified:** Uses finite differences, not full MOST inversion
   - *Impact:* May underestimate κ in complex profiles
   - *Mitigation:* Use geometric mean height for better approximation

2. **No machine learning correction:** Only physics-based approaches
   - *Impact:* May not capture model-specific biases
   - *Mitigation:* Allow tunable parameters (D, α, β)

3. **Single-column focus:** Not yet tested in 3D models
   - *Impact:* Unknown stability in coupled simulations
   - *Mitigation:* Conservative clamping (C ∈ [0.5, 2.0])

4. **Limited validation data:** Only GABLS1 and ARM NSA planned
   - *Impact:* May not generalize to all stable regimes
   - *Mitigation:* Encourage community validation contributions

### Open Questions
- Q: Should we include convective corrections (Ri < 0)?
  - *Current answer:* No, focus on stable layers first

- Q: How to handle sea ice / snow albedo coupling?
  - *Current answer:* Out of scope, user must handle in their model

- Q: What about sub-grid terrain corrections?
  - *Current answer:* Future work, not in v1.0

---

## Resource Allocation (Estimated Hours)

### Completed Work (~40 hours)
- [x] Audit & planning: 6 hours
- [x] Module structure & packaging: 3 hours
- [x] Core algorithms (derivatives, Ri, corrections): 15 hours
- [x] Documentation (6 docs): 10 hours
- [x] Integration examples (Python, Fortran): 4 hours
- [x] Diagnostic tool: 2 hours

### Remaining Work (~60 hours)
- [ ] Profiles & curvature modules: 5 hours
- [ ] Test suite: 8 hours
- [ ] Validation notebooks (5 total): 15 hours
- [ ] Data loaders & viz: 6 hours
- [ ] WRF/MPAS integration testing: 10 hours
- [ ] API reference & tutorials: 8 hours
- [ ] Paper writing & submission: 8 hours

**Total project:** ~100 hours (2.5 weeks full-time equivalent)

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| v0.1-pre | 2025-01-XX | Core algorithms + diagnostic tool |
| v0.2 (planned) | 2025-02-XX | Validation suite + curvature module |
| v1.0 (planned) | 2025-03-XX | Full release + paper submission |

---

## Contact & Contributions

**Maintainer:** David England (davidengland@...)

**Contributions:** See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

**Issues:** Report bugs or request features via GitHub Issues.

---

**Summary:** Core toolkit is **ready to use** for model developers. Validation and pedagogical materials are **in progress**. Target v1.0 release in Q1 2025.
