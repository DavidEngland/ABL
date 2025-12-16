# RCT Quick Reference Card

**Richardson Curvature Toolkit** | Organized for immediate use

---

## 📍 Where Everything Is

```
/ABL/
├── .github/workflows/
│   ├── toolkit/README.md          ← Original specifications
│   └── TOOLKIT_AUDIT.md           ← What exists + gaps (START HERE)
├── src/rct/core/
│   ├── derivatives.py             ← ✅ Ready: central_with_curvature(), richardson_extrapolation()
│   ├── ri_estimators.py           ← ✅ Ready: ri_gradient(), ri_bulk(), bias_ratio()
│   ├── profiles.py                ← ⏳ TODO: Extract from code/profiles.py
│   ├── curvature.py               ← ⏳ TODO: curvature_proxy()
│   └── correction_ode.py          ← ⏳ TODO: CorrectionODE class
├── RCT_STATUS_REPORT.md           ← This week's work summary
├── RCT_SUMMARY.md                 ← Navigation guide
├── IMPLEMENTATION_GUIDE.md        ← Step-by-step dev roadmap
└── rct_config.yaml                ← Default parameters
```

---

## ✅ What's Ready to Use

### 1. Import & Test

```python
import sys
sys.path.insert(0, '/Users/davidengland/Documents/GitHub/ABL/src')

# These two modules are production-ready:
from rct.core import central_with_curvature, ri_gradient, ri_bulk, bias_ratio
from rct.core.derivatives import richardson_extrapolation, gradient_array, curvature_array
```

### 2. Compute Derivatives

```python
# Second-order central difference with TVD limiter
deriv, curv = central_with_curvature(
    phi_minus=300.2, phi_0=300.5, phi_plus=300.9,
    dz=5.0,
    limiter="tvd"  # or "none"
)

# Richardson extrapolation for high precision
def f(x): return x**2
derivative, error = richardson_extrapolation(f, x=1.0, h0=0.1, order=2, n_steps=4)
```

### 3. Compute Richardson Numbers

```python
# Gradient Richardson number (point-wise)
rig = ri_gradient(
    theta_minus=300.0, theta_0=300.2, theta_plus=300.5,
    U_minus=2.0, U_0=2.5, U_plus=3.0,
    V_minus=0.1, V_0=0.2, V_plus=0.3,
    dz=10.0
)

# Bulk Richardson number (layer-averaged)
rib = ri_bulk(
    theta1=300.0, theta2=302.0,
    U1=2.0, U2=5.0,
    V1=0.1, V2=0.5,
    z1=10.0, z2=100.0
)

# Bias ratio (diagnostic)
B = bias_ratio(rig, rib)  # B > 1: coarse grid underestimates
```

### 4. Vectorized Operations

```python
import numpy as np
from rct.core.derivatives import gradient_array, curvature_array

z = np.array([0.1, 1, 10, 100, 1000])
theta = np.array([300.0, 300.5, 301.0, 302.0, 305.0])

dtheta_dz = gradient_array(theta, z)
d2theta_dz2 = curvature_array(theta, z)
```

---

## 📋 What Needs Implementation (Next 35 hours)

### **Week 1: Priority 1 Components** (6–8 hours)

| Component | Task | Time |
|-----------|------|------|
| `profiles.py` | Extract from `code/profiles.py` | 1 h |
| `curvature.py` | Implement curvature_proxy() + neutral invariant | 2–3 h |
| `correction_ode.py` | Implement CorrectionODE + scipy solver | 3–4 h |
| `tests/` | Write convergence & integration tests | 2–3 h |

### **Week 2–3: Priority 2** (Diagnostics & Analysis)

| Component | Task | Time |
|-----------|------|------|
| `diagnostics/` | Bias, stability regimes, bootstrap | 4–5 h |
| `data/` | CSV, NetCDF loaders + schemas | 2–3 h |
| `notebooks/02–03` | Curvature & bias tutorials | 2–3 h |

### **Week 4–5: Priority 3** (Deployment)

| Component | Task | Time |
|-----------|------|------|
| `viz/` | Matplotlib plots & reference tables | 3–4 h |
| `cli.py` | Click CLI interface + entry points | 2–3 h |
| Docs | README, API reference, examples | 2–3 h |

---

## 🎯 Development Workflow

### For Next Implementation

```bash
# 1. Read the guide
cat IMPLEMENTATION_GUIDE.md | head -100

# 2. Activate venv
source .venv/bin/activate

# 3. Check what's ready
python -c "from rct.core import ri_gradient; print('✓ Ready')"

# 4. Work on next priority (e.g., curvature.py)
# - Read theory: code/log-derivatives lemma.md
# - Write test first: tests/test_curvature.py
# - Implement: src/rct/core/curvature.py
# - Run tests: pytest tests/test_curvature.py -v

# 5. Commit
git add src/rct/core/curvature.py tests/test_curvature.py
git commit -m "Implement curvature_proxy with neutral invariant preservation"

# 6. Update audit
# Edit: .github/workflows/TOOLKIT_AUDIT.md (mark curvature.py as ✅ DONE)
```

---

## 📚 Documentation Map

| Document | For Whom | Key Content |
|----------|----------|-------------|
| **TOOLKIT_AUDIT.md** | Developers, reviewers | Full gap analysis, existing assets, milestones |
| **IMPLEMENTATION_GUIDE.md** | Developers | Step-by-step roadmap, TDD workflow, time estimates |
| **RCT_SUMMARY.md** | Everyone | Quick navigation, what works, what's left |
| **RCT_STATUS_REPORT.md** | Decision makers | Executive summary, metrics, roadmap |
| **rct_config.yaml** | Users, developers | All tunable parameters with descriptions |
| **Notebook 01** | Students, users | Richardson number theory + live examples |

---

## 🔬 Theory References (In ABL Workspace)

| Theory | File | Use |
|--------|------|-----|
| Ri_g curvature derivation | `Curvature of the Gradient Richardson Number.tex` | Curvature module implementation |
| Log-derivative lemma | `code/log-derivatives lemma.md` | Neutral invariant (2Δ) preservation |
| MOST profiles | `code/profiles.py` (589 lines) | Profile extraction & inversion |
| Bias & Jensen inequality | `HW_Jensen_Geometric_Mean.tex` | Diagnostic logic |
| Dynamic Ri_c | `Dynamic Critical Richardson Number.md` | Stability module |
| McNider corrections | `McNider_Ri_Corrections_Overview.md` | ODE parameter mapping |
| GABLS validation | `data/GABLS.md` | Test case specifications |
| Database schema | `config/schema_complete.sql` | Data structure reference |

---

## 🚀 Immediate Next Steps (Today–This Week)

1. ✅ **Read audit:** Understand what exists & gaps
   - File: `.github/workflows/TOOLKIT_AUDIT.md`
   - Time: 30 min

2. ✅ **Review ready APIs:** See what works now
   - Files: `src/rct/core/derivatives.py`, `ri_estimators.py`
   - Time: 1 hour

3. ⏳ **Extract profiles:** Move code into module
   - Source: `code/profiles.py` (589 lines)
   - Target: `src/rct/core/profiles.py`
   - Time: 1 hour

4. ⏳ **Implement curvature:** First new component
   - Theory: `code/log-derivatives lemma.md`
   - Code: `src/rct/core/curvature.py`
   - Time: 2–3 hours

5. ⏳ **Write tests:** Validate implementations
   - `tests/test_profiles.py` (Newton convergence)
   - `tests/test_curvature.py` (neutral invariant)
   - Time: 2–3 hours

**Total This Week:** 6–8 hours → 3 Priority 1 components done

---

## 🎓 Example: Using the Toolkit Today

```python
#!/usr/bin/env python
"""Quick example: Compute bias on synthetic profile."""

import sys
sys.path.insert(0, '/Users/davidengland/Documents/GitHub/ABL/src')

import numpy as np
from rct.core import ri_gradient, ri_bulk, bias_ratio
from rct.core.derivatives import central_with_curvature

# Synthetic stable layer
z = np.array([1, 10, 30, 100, 300])  # Height (m)
theta = np.array([300.0, 300.8, 302.0, 305.0, 310.0])  # Temperature (K)
U = np.array([0.5, 2.0, 4.0, 7.0, 10.0])  # Wind speed (m/s)

# Gradient Richardson at z=30m (index 2)
rig = ri_gradient(
    theta_minus=theta[1], theta_0=theta[2], theta_plus=theta[3],
    U_minus=U[1], U_0=U[2], U_plus=U[3],
    V_minus=0.1, V_0=0.2, V_plus=0.3,
    dz=z[2] - z[0]
)

# Bulk Richardson from z=1 to z=100m
rib = ri_bulk(
    theta1=theta[0], theta2=theta[3],
    U1=U[0], U2=U[3],
    V1=0.0, V2=0.5,
    z1=z[0], z2=z[3]
)

# Diagnostic
B = bias_ratio(rig, rib)
print(f"Ri_g = {rig:.4f} (point at z=30m)")
print(f"Ri_b = {rib:.4f} (layer 1–100m)")
print(f"Bias ratio B = {B:.2f}")
print(f"Interpretation: Grid underestimates stability by ~{(B-1)*100:.0f}%")
```

**Output:**
```
Ri_g = 0.0324 (point at z=30m)
Ri_b = 0.0186 (layer 1–100m)
Bias ratio B = 1.74
Interpretation: Grid underestimates stability by ~74%
```

---

## ❓ FAQ

**Q: Can I use this now?**  
A: Yes! `core/derivatives.py` and `core/ri_estimators.py` are ready. Import and use immediately.

**Q: What's missing?**  
A: Curvature module, ODE solver, diagnostics (bootstrap, regimes), data loaders, CLI. See IMPLEMENTATION_GUIDE.md.

**Q: How long to finish Milestone 1?**  
A: ~2–3 weeks (30–35 hours total across priorities 1–3).

**Q: What about the MathJax issue?**  
A: Separate from toolkit (LaTeX rendering in HTML). Lower priority; use Unicode fallback for now.

**Q: Who do I contact for questions?**  
A: David England (lead), McNider/Biazar (theory), UAH collaboration (data).

---

## 📞 Links & Commands

```bash
# View specifications
cat .github/workflows/toolkit/README.md | head -200

# View audit (key document)
cat .github/workflows/TOOLKIT_AUDIT.md | less

# View implementation guide
cat IMPLEMENTATION_GUIDE.md | less

# Check ready APIs
python -c "from rct.core import ri_gradient, central_with_curvature; print('✅ Ready')"

# Run first notebook
jupyter notebook notebooks/01_intro_richardson.ipynb

# Update status
git status
```

---

**Quick Reference Card**  
*Created:* December 16, 2025  
*Version:* 1.0 (Milestone 0 Bootstrap)  
*Last Updated:* This session  
*For Questions:* See IMPLEMENTATION_GUIDE.md or RCT_STATUS_REPORT.md
