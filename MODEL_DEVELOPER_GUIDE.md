# Richardson Number Curvature Corrections for Model Developers

**Purpose:** Help atmospheric modelers diagnose and correct systematic Ri bias in their boundary layer schemes

**Target Audience:** Developers of WRF, MPAS, AM4, CESM, regional models, LES codes, or custom ABL schemes

**Philosophy:** We don't know your code structure. This guide provides flexible tools you can adapt.

---

## Quick Diagnostic: Does Your Model Need Curvature Corrections?

### Red Flags (You Probably Need This)

1. **Coarse vertical resolution** in stable boundary layer:
   - Lowest model level > 10 m AGL
   - First layer thickness Δz > 50 m in SBL
   - Stretched vertical grid (Δz grows rapidly with height)

2. **Excessive mixing in stable conditions:**
   - Nighttime boundary layer too warm/well-mixed
   - LLJ (low-level jet) too weak or missing
   - Arctic/polar simulations struggle with persistent stable layers

3. **Bulk Richardson number used for stability:**
   - Using Ri_b = (g/θ) Δθ Δz / (ΔU)² to determine mixing
   - Critical Richardson number Ri_c ≈ 0.25 as turbulence cutoff
   - Flux-gradient relationships depend on layer-averaged Ri

4. **Grid sensitivity:**
   - Finer vertical resolution → different SBL structure
   - Increasing Δz in stable layers → unrealistically strong mixing

### Green Flags (You're Probably Fine)

- High vertical resolution (Δz < 5 m in SBL)
- LES or DNS (resolved turbulence)
- Using gradient Richardson number Ri_g directly from high-res gradients
- Already implementing McNider-style or other curvature corrections

---

## Assessment Workflow

### Step 1: Extract Profiles from Your Model

**What we need:**
```python
# At one timestep, one location (stable conditions preferred)
data = {
    'z': [z0, z1, z2, ..., zn],        # Height AGL (m)
    'theta': [θ0, θ1, θ2, ..., θn],    # Potential temperature (K)
    'U': [u0, u1, u2, ..., un],        # Zonal wind (m/s)
    'V': [v0, v1, v2, ..., vn],        # Meridional wind (m/s) [optional]
    'L': Obukhov_length,                # Monin-Obukhov length (m) [if available]
    'z0': roughness_length,             # Surface roughness (m)
}
```

**Example output formats:**
- CSV: `z,theta,U,V` (simplest)
- NetCDF: Standard CF conventions
- JSON: Nested structure
- Python dict: Direct use in diagnostic scripts

### Step 2: Run Diagnostic Tool

```python
from rct.diagnostics import assess_model_bias

# Load your profile
profile = load_profile('model_output.csv')  # Or your custom loader

# Assess bias
assessment = assess_model_bias(
    z=profile['z'],
    theta=profile['theta'],
    U=profile['U'],
    V=profile.get('V', None),
    layer_indices=[(0,5), (5,10), (10,15)],  # Define layers to check
    verbose=True
)

print(assessment['summary'])
# Output:
#   Layer 0-5 (Δz=50m):   Bias ratio B = 1.85 → 85% underestimate
#   Layer 5-10 (Δz=100m): Bias ratio B = 2.31 → 131% underestimate
#   Recommendation: Apply corrections to all layers
```

### Step 3: Understand Your Model's Ri Calculation

**Common patterns we'll help you identify:**

#### Pattern A: Bulk Ri with Stability Function
```fortran
! Your code might look like:
Ri_bulk = (g/theta_ref) * (theta(k+1) - theta(k)) * dz / wind_shear**2
IF (Ri_bulk < Ri_crit) THEN
    K_m = f_m(Ri_bulk) * K_neutral  ! Stability function
    K_h = f_h(Ri_bulk) * K_neutral
END IF
```

**→ Correction strategy:** Adjust Ri_bulk or modify f_m/f_h

#### Pattern B: Gradient Ri from Finite Differences
```python
# Your code might look like:
dtheta_dz = (theta[k+1] - theta[k-1]) / (2*dz)
du_dz = (u[k+1] - u[k-1]) / (2*dz)
Ri_g = (g/theta0) * dtheta_dz / du_dz**2
```

**→ Correction strategy:** Adjust finite difference stencil or add curvature term

#### Pattern C: MOST-Based (Monin-Obukhov)
```python
# Your code might look like:
zeta = z / L_obukhov
phi_m = stability_function_momentum(zeta)
phi_h = stability_function_heat(zeta)
Ri_MOST = zeta * phi_h / phi_m**2
```

**→ Correction strategy:** Modify φ_m, φ_h, or add grid-dependent tail extension

#### Pattern D: TKE-Based Closure (e.g., Mellor-Yamada)
```fortran
! Mixing length limited by Ri:
l_mix = MIN(kappa*z, l_max / (1 + alpha*Ri)**beta)
```

**→ Correction strategy:** Adjust Ri in mixing length formula or modify α, β

---

## Implementation Strategies

### Strategy 1: Correct Ri Before Use (Simplest)

**Idea:** Compute corrected Ri_bulk that accounts for curvature, then feed into existing stability functions.

```python
# In your model's stability calculation subroutine:

# BEFORE (original code):
Ri_bulk = compute_bulk_richardson(theta, U, V, z)
K_m = stability_function(Ri_bulk)

# AFTER (with correction):
Ri_bulk_raw = compute_bulk_richardson(theta, U, V, z)
Ri_bulk_corrected = apply_curvature_correction(
    Ri_bulk_raw, 
    dz=z[k+1] - z[k],
    L_obukhov=L,
    curvature_params={'alpha': 0.3, 'beta': 0.7}
)
K_m = stability_function(Ri_bulk_corrected)
```

**Pros:**
- Minimal code changes
- Easy to A/B test (toggle correction on/off)
- Works with any downstream stability function

**Cons:**
- Doesn't fundamentally fix discretization issues
- May need tuning of α, β parameters

---

### Strategy 2: Modify Stability Functions (More Accurate)

**Idea:** Extend f_m(Ri) and f_h(Ri) with grid-dependent tail for large Ri.

```python
def stability_function_corrected(Ri, dz, dz_ref=10.0):
    """
    Grid-aware stability function with extended tail.
    
    Standard form: f(Ri) = 1 / (1 + α Ri)
    Corrected form: f(Ri, Δz) = f(Ri) * C(Ri, Δz)
    
    where C(Ri, Δz) = exp(-D * (Δz/Δz_ref)^p * (Ri/Ri_ref)^q)
    """
    # Base stability function (your existing form)
    f_base = 1.0 / (1.0 + 5.0 * Ri)
    
    # Grid correction factor
    D = 0.5  # Damping amplitude (tune to your model)
    p = 0.5  # Grid sensitivity exponent
    q = 1.0  # Ri sensitivity exponent
    Ri_ref = 0.25
    
    C = np.exp(-D * (dz/dz_ref)**p * (Ri/Ri_ref)**q)
    
    return f_base * C
```

**Pros:**
- Physics-based correction
- Preserves neutral limit (C→1 as Ri→0)
- Can be validated against LES/tower data

**Cons:**
- Requires tuning D, p, q
- More complex code modification

---

### Strategy 3: Dynamic Critical Richardson Number (Advanced)

**Idea:** Adjust Ri_c threshold based on grid resolution and local curvature.

```python
def dynamic_Ri_critical(Ri_bulk, dz, curvature_proxy, Ri_c0=0.25):
    """
    Compute grid- and curvature-aware critical Richardson number.
    
    Standard: Turbulence cutoff at Ri = Ri_c0 = 0.25
    Dynamic:  Ri_c* = Ri_c0 * [1 + γ * (B - 1)]
    
    where B = bias ratio (≈ 1 + curvature effect)
    """
    # Estimate bias from grid spacing and curvature
    B = estimate_bias_ratio(dz, curvature_proxy)
    
    gamma = 0.5  # Sensitivity parameter (tune to data)
    Ri_c_star = Ri_c0 * (1.0 + gamma * (B - 1.0))
    
    # Clamp to reasonable range
    return np.clip(Ri_c_star, 0.15, 0.5)

# In your turbulence cutoff logic:
Ri_c = dynamic_Ri_critical(Ri_bulk, dz, curvature)
if Ri_bulk < Ri_c:
    # Turbulence active
    K_m = compute_mixing(...)
else:
    # Laminar regime
    K_m = K_laminar
```

**Pros:**
- Physically intuitive (higher threshold for coarse grids)
- Preserves turbulence in realistic stable layers
- Can prevent premature cutoff

**Cons:**
- Needs curvature proxy (extra computation)
- Requires careful tuning

---

## Curvature Proxy Calculation

**Problem:** How to estimate curvature d²Ri_g/dζ² without knowing full MOST parameters?

**Solution:** Use local profile shape:

```python
def estimate_curvature_proxy(theta, U, V, z):
    """
    Estimate combined curvature from profile second derivatives.
    
    Returns dimensionless proxy κ ∈ [0, ∞) where:
      κ ≈ 0: Nearly linear profiles (small curvature, low bias)
      κ > 1: Strong curvature (large bias expected)
    """
    # Compute second derivatives
    d2theta = np.gradient(np.gradient(theta, z), z)
    d2U = np.gradient(np.gradient(U, z), z)
    
    # Normalize by layer scales
    theta_scale = np.std(theta)
    U_scale = np.std(U)
    z_scale = z[-1] - z[0]
    
    # Combined curvature magnitude
    kappa_theta = np.abs(d2theta) * z_scale**2 / theta_scale
    kappa_U = np.abs(d2U) * z_scale**2 / U_scale
    
    return np.mean(kappa_theta + kappa_U)
```

**Use cases:**
- Trigger corrections only when κ > threshold
- Scale correction amplitude by κ
- Diagnostic output for model tuning

---

## Testing & Validation Framework

### Test 1: Neutral Profile (No Correction Needed)

```python
# Generate neutral profile: Ri ≈ 0 everywhere
z = np.linspace(1, 100, 20)
theta = 300.0 + 0.001 * z  # Tiny lapse
U = 5.0 * np.ones_like(z)  # Constant wind

# Apply correction
Ri_corrected = apply_correction(theta, U, z)

# EXPECTED: Ri_corrected ≈ Ri_original (no change)
assert np.allclose(Ri_corrected, Ri_original, rtol=0.05)
```

### Test 2: Stable Layer (Correction Active)

```python
# Generate stable layer with known bias
z = np.array([1, 10, 50, 100, 200])
theta = np.array([300.0, 300.5, 302.0, 304.0, 307.0])  # Strong inversion
U = np.array([2.0, 3.0, 6.0, 8.0, 10.0])  # Shear layer

# Compute bias ratio
B = compute_bias_ratio(theta, U, z, layer=(1, 4))  # Layer 10-200m

# EXPECTED: B > 1.5 (significant underestimate)
assert B > 1.5

# Apply correction
Ri_corrected = apply_correction(theta, U, z)

# EXPECTED: Corrected Ri closer to gradient Ri at geometric mean height
z_g = np.sqrt(z[1] * z[4])  # Geometric mean
Ri_g_target = compute_Ri_gradient_at(z_g, theta, U, z)
Ri_b_corrected = Ri_corrected[2]  # Middle layer

assert np.abs(Ri_b_corrected - Ri_g_target) / Ri_g_target < 0.2  # Within 20%
```

### Test 3: LES/Tower Validation (Real Data)

```python
# Load reference data (GABLS, ARM NSA, CASES-99, etc.)
ref = load_reference('GABLS1_hour9.nc')

# Run your model (offline single-column mode)
model_output = run_model_offline(
    initial_profile=ref['initial'],
    forcing=ref['forcing'],
    timesteps=ref['time']
)

# Compare key metrics
metrics = {
    'h_BL': boundary_layer_height,
    'LLJ_height': low_level_jet_peak_height,
    'LLJ_speed': low_level_jet_peak_speed,
    'surface_flux': sensible_heat_flux,
}

for metric in metrics:
    model_val = model_output[metric]
    ref_val = ref[metric]
    error = abs(model_val - ref_val) / ref_val
    print(f"{metric}: {error*100:.1f}% error")
    assert error < 0.3  # Within 30% (adjust threshold as needed)
```

---

## Common Model Architectures & Integration Points

### WRF (Weather Research & Forecasting)

**File to modify:** `phys/module_bl_ysu.F` or `module_bl_mynn.F`

**Location:** Subroutine computing Ri and stability functions
```fortran
! Look for lines like:
CALL SFCLAYINIT(...)
! or
Ri = ALOG((Z(K+1)-Z0)/(Z(K)-Z0)) * ...

! Add correction call:
CALL COMPUTE_RI_CORRECTION(Ri, dz, L_obukhov, Ri_corrected)
```

**Namelist control:**
```fortran
&physics
  ri_curvature_correction = .true.
  ri_corr_alpha = 0.3
  ri_corr_beta = 0.7
/
```

---

### MPAS (Model for Prediction Across Scales)

**File to modify:** `src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F`

**Key subroutine:** `driver_pbl` → calls to `pbl_compute_gradient_richardson`

**Integration:**
```fortran
! After computing gradient Richardson:
call compute_curvature_proxy(theta, u, v, mesh%z, kappa)
call apply_grid_correction(Ri_gradient, dz, kappa, Ri_corrected)

! Use Ri_corrected downstream
```

---

### CESM/CAM (Community Earth System Model)

**File to modify:** `src/physics/cam/vertical_diffusion.F90`

**Subroutine:** `compute_vdiff` → Ri calculation before eddy diffusivity

**Integration:**
```fortran
! After line computing local Richardson number:
ri_local = ... ! Original computation

! Add correction module:
if (ri_curvature_correction_enabled) then
    call apply_curvature_correction(ri_local, dz(k), L_obukhov, ri_local)
endif
```

---

### Custom Python/Julia Models

**Example Python integration:**

```python
class BoundaryLayerScheme:
    def __init__(self, enable_curvature_correction=False):
        self.enable_correction = enable_curvature_correction
        if enable_curvature_correction:
            from rct.core import RiCurvatureCorrector
            self.corrector = RiCurvatureCorrector(alpha=0.3, beta=0.7)
    
    def compute_mixing(self, theta, U, V, z):
        # Original Ri calculation
        Ri_bulk = self._compute_bulk_ri(theta, U, V, z)
        
        # Apply correction if enabled
        if self.enable_correction:
            dz = np.diff(z)
            Ri_bulk = self.corrector.correct(Ri_bulk, dz, theta, U, V, z)
        
        # Downstream mixing calculation
        K_m = self._stability_function(Ri_bulk)
        return K_m
```

---

## Parameter Tuning Guide

### Default Parameters (Start Here)

```python
params = {
    'alpha': 0.3,        # Curvature weight in correction ODE
    'beta': 0.7,         # Shear damping weight
    'D': 0.5,            # Damping amplitude for stability function
    'p': 0.5,            # Grid sensitivity exponent (Δz^p)
    'q': 1.0,            # Ri sensitivity exponent (Ri^q)
    'Ri_ref': 0.25,      # Reference Richardson number
    'dz_ref': 10.0,      # Reference grid spacing (m)
}
```

### Tuning Workflow

1. **Start with defaults** → Run GABLS1 case
2. **If too much mixing** (h_BL too high, LLJ too weak):
   - Increase α (0.3 → 0.5)
   - Increase D (0.5 → 0.8)
3. **If too little mixing** (h_BL too low, unrealistic stratification):
   - Decrease α (0.3 → 0.1)
   - Decrease D (0.5 → 0.3)
4. **If grid-sensitive** (results change a lot with Δz):
   - Increase p (0.5 → 0.8) for stronger grid dependence
5. **Validate** → Run full suite (GABLS1-4, ARM NSA, CASES-99)

---

## Diagnostic Outputs (Add to Your Model)

**Recommended new output variables:**

```python
output_vars = {
    'Ri_bulk_raw': 'Bulk Richardson number (uncorrected)',
    'Ri_bulk_corrected': 'Bulk Richardson number (curvature-corrected)',
    'bias_ratio': 'B = Ri_g / Ri_b diagnostic',
    'curvature_proxy': 'Estimated curvature magnitude κ',
    'correction_factor': 'C(z) applied to stability functions',
    'Ri_critical_dynamic': 'Grid-aware critical Richardson number',
}
```

**Use for:**
- Debugging correction implementation
- Tuning parameters
- Publication figures showing correction impact

---

## Pre-Built Tools We Provide

### 1. Diagnostic Script

```bash
python -m rct.diagnose model_profile.csv --output report.html
```

Generates:
- Bias ratio plots by layer
- Recommendations (correction needed? which layers?)
- Comparison to MOST theory

### 2. Correction Module (Fortran/Python/Julia)

```python
from rct.core import CurvatureCorrector

corrector = CurvatureCorrector(alpha=0.3, beta=0.7)
Ri_corrected = corrector.apply(Ri_raw, dz, L_obukhov, curvature_proxy)
```

Fortran wrapper:
```fortran
CALL rct_correct_ri(Ri_raw, dz, L, kappa, Ri_corrected)
```

### 3. Validation Suite

```python
from rct.validation import run_validation

results = run_validation(
    your_model_function,
    test_cases=['GABLS1', 'GABLS3', 'ARM_NSA'],
    metrics=['h_BL', 'LLJ_speed', 'surface_flux']
)
print(results.summary())
```

---

## Quick Start Checklist

- [ ] Extract 1 stable profile from your model (CSV: z, theta, U, V)
- [ ] Run diagnostic: `python -m rct.diagnose profile.csv`
- [ ] Read bias ratio → Is B > 1.3? (correction recommended)
- [ ] Identify your model's Ri calculation pattern (A/B/C/D above)
- [ ] Choose implementation strategy (1/2/3)
- [ ] Add correction to 1 test case (GABLS1 or tower case)
- [ ] Compare before/after: h_BL, LLJ, surface flux
- [ ] Tune α, β if needed
- [ ] Validate on full test suite
- [ ] Add diagnostic outputs to production runs

---

## Support & Collaboration

**We provide:**
- Code review of your integration
- Parameter tuning assistance
- Custom corrections for unusual model structures
- Co-authorship on validation papers

**Contact:**
- David England (davidengland@uah.edu)
- Dick McNider, Arastoo Pour-Biazar (UAH)

**Repository:**
- GitHub: `ABL/src/rct/` (flexible modules)
- Examples: `ABL/examples/model_integration/`

**Citation:**
```
England et al. (2025). Curvature-Aware Richardson Number Corrections
for Atmospheric Boundary Layer Models. In preparation.
```

---

**Last Updated:** December 16, 2025  
**Version:** 1.0 (Model Developer Focus)
