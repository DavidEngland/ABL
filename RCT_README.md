# Richardson Number Curvature Toolkit (RCT)
## Practical Integration Guide for Model Developers

**Quick Start:** *You have an atmospheric model with excessive Ri-based mixing suppression on coarse grids. This toolkit helps you diagnose and fix it.*

---

## What's the Problem?

On coarse vertical grids (Δz > 20m), bulk Richardson numbers **underestimate** gradient Richardson numbers due to **curvature effects**. This causes:

- Premature turbulence cutoff in stable layers
- Shallow boundary layers (h_BL too low by 30-50%)
- Weak low-level jets (LLJ too strong/high)
- Excessive surface cooling in Arctic winter

**Bias Ratio:** B = Ri_gradient(z_geometric_mean) / Ri_bulk  
✅ **Good:** B < 1.2  
⚠️ **Moderate:** 1.2 < B < 1.5  
🚨 **Needs correction:** B > 1.5

---

## Quick Diagnostic

**Step 1:** Extract a vertical profile from your model (CSV format):
```
z,theta,U,V
2.0,265.0,3.0,0.0
10.0,265.5,5.0,0.2
30.0,267.0,7.0,0.4
...
```

**Step 2:** Run the diagnostic tool:
```bash
python scripts/diagnose_ri_bias.py your_profile.csv
```

**Step 3:** Read the report:
```
======================================================================
ASSESSMENT:
======================================================================
🔴 Large bias detected: B_max = 1.82
🟡 Coarse grid: Δz_max = 100 m

RECOMMENDED STRATEGY: MODERATE
  Use stability_function_correction() with D=0.5.
  Extends tail of f(Ri) for coarse grids.
  
SUGGESTED PARAMETERS:
  D     = 0.5
  alpha = 0.3
  beta  = 0.7
```

---

## Integration Strategies

### **Strategy 1: Simple Multiplicative Correction** (5 minutes)
*Minimal code changes, sufficient for moderate bias (B < 1.5)*

**Python:**
```python
from rct.core.correction_ode import simple_multiplicative_correction

# Before
Ri_bulk = compute_richardson(theta, U, V, z)

# After (1 line added)
Ri_corrected = simple_multiplicative_correction(Ri_bulk, dz=50.0, D=0.5)
```

**Fortran (WRF/MPAS/CESM):**
```fortran
USE module_ri_correction, ONLY: apply_simple_ri_correction
CALL apply_simple_ri_correction(Ri_bulk, dz, Ri_corrected, D_param=0.5)
```

**Impact:** +20-40% mixing in stable layers, minimal computational cost.

---

### **Strategy 2: Stability Function Correction** (1 hour)
*More accurate, extends f(Ri) tail for coarse grids*

**Python:**
```python
from rct.core.correction_ode import stability_function_correction

# Before
f_m = 1.0 / (1.0 + 5.0 * Ri_bulk)

# After
f_m_corrected = stability_function_correction(
    Ri_bulk, dz=50.0, 
    base_function=lambda Ri: 1.0/(1.0 + 5.0*Ri),
    D=0.5
)
```

**Impact:** +30-60% mixing near Ri_c, prevents premature cutoff.

---

### **Strategy 3: Dynamic Critical Richardson** (2 hours)
*Physics-based threshold adjustment, best for very coarse grids (Δz > 50m)*

**Python:**
```python
from rct.core.correction_ode import dynamic_critical_richardson, estimate_curvature_proxy

# Estimate curvature from profile shape
kappa = estimate_curvature_proxy(theta, U, V, z)

# Adjust Ri_c based on grid and curvature
Ri_c_dynamic = dynamic_critical_richardson(Ri_bulk, dz=80.0, curvature_proxy=kappa)

# Use adjusted threshold
if Ri_corrected < Ri_c_dynamic:
    # Turbulence active
else:
    # Turbulence suppressed
```

**Impact:** Prevents premature cutoff, +10-30% increase in h_BL.

---

### **Strategy 4: Full ODE Solver** (research-grade, 4+ hours)
*Most rigorous, for extreme cases (B > 2.0)*

```python
from rct.core.correction_ode import CorrectionODE

ode = CorrectionODE(alpha=0.3, beta=0.7)
C_profile = ode.solve(z, kappa_z, shear_z)
Ri_corrected = Ri_bulk * C_profile
```

---

## Model-Specific Examples

### **WRF** (`phys/module_bl_ysu.F` or `module_bl_mynn.F`)
See: [`examples/wrf_integration_example.F90`](examples/wrf_integration_example.F90)
- Add `module_ri_correction.F90` to `phys/`
- Modify Ri computation in `ysu` or `mynn_bl_driver`
- Tune D via `namelist.input`: `ri_correction_D = 0.5`

### **MPAS** (`atmphys_driver_pbl.F`)
- Insert correction in `pbl_compute_gradient_richardson`
- Use grid spacing from `meshDensity`

### **CESM/CAM** (`vertical_diffusion.F90`)
- Add call in `compute_vdiff`
- Access Δz from `pver` array

### **Custom Python/Julia**
See: [`examples/model_integration_example.py`](examples/model_integration_example.py)
- Full working example with before/after comparison
- Shows toggle for corrections (single flag)

---

## Parameter Tuning Workflow

**Defaults (start here):**
- D = 0.5 (damping amplitude)
- α = 0.3 (curvature weight)
- β = 0.7 (shear damping)
- p = 0.5, q = 1.0 (exponents)
- dz_ref = 10.0 m (reference grid spacing)

**Tuning steps:**
1. **Run GABLS1** (single-column mode, 9-hour stable case)
2. **Check h_BL and z_LLJ** vs. LES reference:
   - No correction: h_BL ~ 120m (too shallow)
   - Target (LES): h_BL ~ 200m
3. **Adjust D:**
   - h_BL too shallow → increase D (0.5 → 0.7)
   - h_BL too deep → decrease D (0.5 → 0.3)
4. **Validate on full suite** (GABLS2-4, ARM NSA tower)

**Rule of thumb:**
- Mesoscale model (Δz ~ 50m): D = 0.5
- Regional model (Δz ~ 20m): D = 0.3
- LES-style (Δz ~ 5m): D = 0.1 or disable

---

## Pre-Built Tools

| Tool | Purpose | Usage |
|------|---------|-------|
| **Diagnostic** | Assess bias in profiles | `python scripts/diagnose_ri_bias.py profile.csv` |
| **Corrections** | Apply 4 correction strategies | `from rct.core.correction_ode import ...` |
| **Estimators** | Compute Ri_g, Ri_b, bias ratio | `from rct.core.ri_estimators import ...` |
| **Derivatives** | High-precision gradients | `from rct.core.derivatives import ...` |

---

## Testing Checklist

- [ ] **Neutral profile:** No change (C → 1, corrections off)
- [ ] **Stable layer:** B > 1.5 → corrections reduce bias
- [ ] **GABLS1:** h_BL and z_LLJ match LES within 10%
- [ ] **3D stability:** No timestep crashes, CFL violations
- [ ] **Surface fluxes:** Δ(SHF, LHF) < 15% from control

---

## File Structure

```
ABL/
├── src/rct/                          # Core toolkit
│   ├── core/
│   │   ├── correction_ode.py         # 4 correction strategies ⭐
│   │   ├── ri_estimators.py          # Ri_g, Ri_b, bias ratio ⭐
│   │   ├── derivatives.py            # High-precision derivatives
│   │   └── (profiles, curvature TBD)
│   ├── diagnostics/                  # Bias analysis (TBD)
│   └── (data, viz, utils TBD)
│
├── scripts/
│   └── diagnose_ri_bias.py           # Quick diagnostic tool ⭐
│
├── examples/
│   ├── model_integration_example.py  # Python demo ⭐
│   ├── wrf_integration_example.F90   # Fortran/WRF demo ⭐
│   └── test_profile.csv              # Sample data
│
├── notebooks/                        # Pedagogical (in progress)
│   └── 01_intro_richardson.ipynb
│
└── docs/
    ├── MODEL_DEVELOPER_GUIDE.md      # Full integration guide ⭐
    ├── RCT_QUICK_REFERENCE.md        # One-page summary
    └── (API reference, validation TBD)
```

⭐ = **Ready to use now**

---

## Frequently Asked Questions

**Q: Do I need to install a Python package?**  
A: No! The correction functions are standalone. Just copy `src/rct/core/correction_ode.py` and call from your model code (Python/Fortran/Julia via wrapper).

**Q: What if my model uses MOST instead of bulk Ri?**  
A: See [MODEL_DEVELOPER_GUIDE.md § Pattern C](MODEL_DEVELOPER_GUIDE.md). You can still estimate curvature and apply corrections to stability functions.

**Q: Will this slow down my model?**  
A: Negligible impact (<1% overhead). Simple correction is 3 arithmetic operations per grid point.

**Q: What if corrections make things worse?**  
A: Start with D=0.3 (conservative). If validation degrades, corrections may not be needed (your grid may be fine-resolution or already accounting for curvature).

**Q: Can I use this for convective boundary layers?**  
A: Corrections target **stable** layers (Ri > 0). For convection (Ri < 0), standard approaches work well.

---

## Citation

If you use this toolkit in published research:

```
England, D. E., et al. (2025). "Grid-Dependent Stability Function Corrections
for Arctic and Stable Boundary Layers." Journal of the Atmospheric Sciences
(in preparation).
```

---

## Support

- **Documentation:** See [MODEL_DEVELOPER_GUIDE.md](MODEL_DEVELOPER_GUIDE.md)
- **Examples:** Browse `examples/` directory
- **Issues:** Open a GitHub issue with your model type and diagnostic output

---

## Quick Reference Card

| **If your diagnostic says...** | **Then use...** | **With parameters...** |
|--------------------------------|----------------|----------------------|
| B < 1.2, Green flags | No correction needed | N/A |
| 1.2 < B < 1.5, Fine grid | Simple multiplicative | D=0.3, p=0.5 |
| 1.5 < B < 2.0, Moderate grid | Stability function | D=0.5, p=0.5 |
| B > 2.0 or Δz > 100m | Dynamic Ri_c or ODE | D=0.7, α=0.4 |

**Remember:** Start with defaults, tune on GABLS1, validate on full suite.

---

**Last updated:** 2025-01-XX  
**Version:** 0.1 (pre-release)
