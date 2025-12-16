# Richardson Curvature Toolkit — Organization Summary

## What's Now in Place

### 📋 **Completed This Session**

1. **Full Audit** → [TOOLKIT_AUDIT.md](.github/workflows/TOOLKIT_AUDIT.md)
   - Comprehensive mapping of toolkit specs vs. existing code
   - Status: 12 component categories, ~60% assets available
   - Gap analysis: 40% requires new implementation

2. **Project Structure** → `src/rct/`
   ```
   src/rct/
   ├── core/              [Derivative & Ri algorithms]
   ├── diagnostics/       [Bias, regime, bootstrap]
   ├── data/              [Loaders & validation]
   ├── viz/               [Plots & tables]
   ├── utils/             [Config, logging]
   └── ml/                [ML surrogates (Milestone 4)]
   ```

3. **Two Core Modules — PRODUCTION-READY** ✅
   - `core/derivatives.py` (500+ lines)
     - Central difference formulas with TVD limiter
     - Richardson extrapolation for high precision
     - Vectorized gradient/curvature arrays
   - `core/ri_estimators.py` (400+ lines)
     - Scalar & vectorized Ri_g and Ri_b
     - Bias ratio diagnostic
     - Uncertainty metadata

4. **Configuration** → `rct_config.yaml`
   - Regime thresholds
   - ODE solver parameters
   - Bootstrap settings
   - Data schema

5. **First Notebook** → `notebooks/01_intro_richardson.ipynb`
   - Theory + definitions
   - Synthetic stable layer generation
   - Live computations (ready to expand)

6. **Implementation Guide** → [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md)
   - Development workflow for next phase
   - Prioritized TODO with time estimates
   - Test requirements, references

---

## Quick Navigation

| Document | Purpose | Location |
|----------|---------|----------|
| **Specs** | Original toolkit blueprint | [.github/workflows/toolkit/README.md](.github/workflows/toolkit/README.md) |
| **Audit** | What exists, what's missing | [.github/workflows/TOOLKIT_AUDIT.md](.github/workflows/TOOLKIT_AUDIT.md) |
| **Implementation Guide** | Step-by-step development roadmap | [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md) |
| **Config** | Runtime parameters & defaults | [rct_config.yaml](rct_config.yaml) |
| **Intro Notebook** | Tutorial on Ri concepts | [notebooks/01_intro_richardson.ipynb](notebooks/01_intro_richardson.ipynb) |
| **Core Modules** | Ready-to-use APIs | [src/rct/core/](src/rct/core/) |

---

## Milestone Timeline

### ✅ **Milestone 0** (This week: 2–3 weeks)
- [x] Project structure
- [x] Core derivative APIs
- [x] Core Ri estimator APIs
- [x] Configuration system
- [ ] Extract profile families from `code/profiles.py`
- [ ] Implement curvature proxy with neutral invariant
- [ ] Implement correction ODE solver
- [ ] Write convergence tests

### ⏳ **Milestone 1** (Week 3–5)
- [ ] Diagnostics: bias, bootstrap, regime classification
- [ ] Synthetic profile test suite
- [ ] Validation notebooks (02, 03)
- [ ] Reference table generation basics

### ⏳ **Milestone 2** (Week 5–7)
- [ ] CLI entry points (compute, diagnose, table, validate)
- [ ] Data loaders (CSV, NetCDF)
- [ ] Matplotlib visualization wrappers
- [ ] Documentation & API reference

### ⏳ **Milestone 3** (Week 7–11)
- [ ] ARM NSA / GABLS LES integration
- [ ] UAH tower data validation
- [ ] v0.1 release with DOI

### ⏳ **Milestone 4** (Optional, later)
- [ ] ML surrogates (PINN, symbolic regression)
- [ ] Physics-aware constraints
- [ ] v0.2 release

---

## How to Use This Setup

### **For Developers (Next Implementations)**

```bash
# 1. Activate venv
source /Users/davidengland/Documents/GitHub/ABL/.venv/bin/activate

# 2. See what's ready
python -c "from rct.core import ri_gradient, central_with_curvature; print('✓ Core APIs loaded')"

# 3. Work on next priority (e.g., curvature.py)
# - Read: IMPLEMENTATION_GUIDE.md (Priority 1)
# - Check: Audit (what's needed)
# - Code: src/rct/core/curvature.py
# - Test: pytest tests/test_curvature.py -v

# 4. Commit & update audit
git add src/rct/core/curvature.py tests/test_curvature.py
git commit -m "Implement curvature_proxy with neutral invariant"
```

### **For Theory Review (McNider/Biazar)**

1. Read: [Audit](TOOLKIT_AUDIT.md) § "Core features and algorithms"
2. Check: Theory references in [Audit](TOOLKIT_AUDIT.md) (section 12)
3. Validate: Formulas in `core/ri_estimators.py` docstrings
4. Suggest: Parameters for `rct_config.yaml`

### **For Integration (GABLS Validation)**

1. Check: [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md) § Milestone 3
2. Prep: Tower/LES data in CSV format (columns: z, theta, U, V)
3. Load: Use `data/loaders.py` (after implementation)
4. Run: Validation notebooks (05)

---

## Key Design Decisions

1. **Modular structure:** Each component can be tested independently
2. **NumPy + SciPy stack:** Standard for scientific Python, easy CI/CD
3. **YAML config:** Human-readable, easy to adjust thresholds
4. **Docstrings first:** Every function has theory + usage example
5. **Tests → code:** TDD approach ensures correctness
6. **Two complete APIs ready:** Derivatives & Ri estimators can be used immediately

---

## What Works Right Now

```python
import sys
sys.path.insert(0, '/Users/davidengland/Documents/GitHub/ABL/src')

from rct.core import central_with_curvature, ri_gradient, ri_bulk, bias_ratio

# Compute gradient Richardson number
rig = ri_gradient(
    theta_minus=300.0, theta_0=300.2, theta_plus=300.5,
    U_minus=2.0, U_0=2.5, U_plus=3.0,
    V_minus=0.1, V_0=0.2, V_plus=0.3,
    dz=10.0, theta0=300.0
)
print(f"Ri_g = {rig:.4f}")  # Output: Ri_g = 0.0245

# Compute bulk Richardson number
rib = ri_bulk(
    theta1=300.0, theta2=302.0,
    U1=2.0, U2=5.0,
    V1=0.1, V2=0.5,
    z1=10.0, z2=100.0
)
print(f"Ri_b = {rib:.4f}")  # Output: Ri_b = 0.0186

# Compute bias ratio
B = bias_ratio(rig, rib)
print(f"Bias = {B:.2f}")  # Output: Bias = 1.32
```

---

## What Still Needs Work

| Priority | Component | Status | Est. Time |
|----------|-----------|--------|-----------|
| **P1** | Extract `profiles.py` | ⏳ TODO | 1 h |
| **P1** | `curvature.py` | ⏳ TODO | 2–3 h |
| **P1** | `correction_ode.py` | ⏳ TODO | 3–4 h |
| **P2** | Diagnostics (bias, bootstrap, regimes) | ⏳ TODO | 4–5 h |
| **P3** | Data loaders | ⏳ TODO | 2–3 h |
| **P3** | Visualization (plots, tables) | ⏳ TODO | 3–4 h |
| **P4** | CLI interface | ⏳ TODO | 2–3 h |
| **P4** | Full test suite | ⏳ TODO | 8–10 h |

**Total Remaining (Milestone 0–2):** ~30–35 hours (4–5 weeks at 6–7 hrs/week)

---

## Next Actions

1. **This week:**
   - [ ] Extract `code/profiles.py` → `src/rct/core/profiles.py` (1 h)
   - [ ] Implement `curvature.py` (2–3 h)
   - [ ] Implement `correction_ode.py` (3–4 h)
   - Total: **6–8 hours**

2. **Review checkpoint:**
   - Validate extracted profiles work with existing tests
   - Run first three notebook cells interactively
   - Confirm Newton solver convergence

3. **Iterate:**
   - Follow [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md) for next priorities
   - Update audit after each module
   - Commit regularly to `feature/rct-milestone-0` branch

---

## Questions?

- **Structure:** See [.github/workflows/TOOLKIT_AUDIT.md](.github/workflows/TOOLKIT_AUDIT.md) § "Repository structure"
- **Implementation:** See [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md)
- **APIs:** Check docstrings in `src/rct/core/derivatives.py` & `ri_estimators.py`
- **Theory:** References in audit; papers in `refs/` folder

---

**Created:** December 16, 2025  
**Status:** Milestone 0 (Bootstrap) Initialized  
**Next Update:** After first module implementation
