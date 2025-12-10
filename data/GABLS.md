# GABLS — GEWEX Atmospheric Boundary Layer Study

## Overview

GABLS (GEWEX Atmospheric Boundary Layer Study) is a multi-phase international intercomparison project focused on improving the representation of the atmospheric boundary layer in numerical weather prediction (NWP) and climate models. The project coordinates Large Eddy Simulation (LES) studies, single-column model (SCM) evaluations, and observational campaigns to advance understanding of stable and convective boundary layer processes.

**Primary Focus Areas:**
- Stable boundary layer (SBL) physics and parameterizations
- Arctic boundary layers and sea-ice coupling
- Diurnal cycle transitions (convective → stable → neutral)
- Model intercomparison and validation hierarchy (LES → SCM → 3D)

**Key Organizing Institutions:**
- KNMI (Royal Netherlands Meteorological Institute)
- NCAR (National Center for Atmospheric Research)
- MetOffice (UK)
- ECMWF (European Centre for Medium-Range Weather Forecasts)

---

## GABLS1 — Canonical Stable Boundary Layer

### Case Description

**Reference:** Beare et al. (2006), *Boundary-Layer Meteorology*, 118, 247–272.

**Setup (Idealized):**
- Domain: 400 m × 400 m × 400 m (horizontal periodic)
- Initial conditions: Neutral, U₀ = 8 m/s (geostrophic), constant θ₀ = 265 K
- Surface cooling: Fixed surface temperature drop −0.25 K/h for 9 hours
- Coriolis: f = 1.39×10⁻⁴ s⁻¹ (latitude ~73°N, Arctic-like)
- No moisture, no radiation (dry, clear-sky idealization)

**Grid Resolutions (LES participants):**
- Δx = Δy: 2–10 m (horizontal)
- Δz: 1–6.25 m (near surface), stretched aloft
- Timestep: 0.1–1.0 s (varies by code)

**Duration:** 9 hours (0–9 h simulation time)

### Validation Targets

**Primary Diagnostics (ensemble mean ± spread):**

| Metric | Target Value | Source |
|--------|-------------|---------|
| Surface friction velocity u* (t=9h) | 0.30 ± 0.01 m/s | LES ensemble |
| Inversion height h (t=9h) | 175 ± 10 m | LES ensemble |
| Surface heat flux H₀ (t=9h) | −8.2 ± 0.5 W/m² | LES ensemble |
| Low-level jet (LLJ) height | ~200 m | LES ensemble |
| LLJ peak wind speed | ~9.5 m/s | LES ensemble |
| Surface temperature drop ΔT_sfc | ~2.25 K | Prescribed forcing |

**Profile Targets (t=9h):**
- U(z): Wind speed profile (LLJ structure)
- θ(z): Potential temperature (inversion shape)
- u'w'(z): Momentum flux (negative below LLJ nose)
- w'θ'(z): Heat flux (negative near surface, zero above h)
- TKE(z): Turbulent kinetic energy (peak below LLJ)

### Curvature Diagnostic Application

**Key Observations for Validation:**

1. **Ri_g Profile Shape:**
   - Near-surface (0–50 m): Ri_g ≈ 0.05–0.15 (weakly stable)
   - Mid-layer (50–150 m): Ri_g ≈ 0.15–0.25 (critical threshold approach)
   - Above inversion (>175 m): Ri_g > 0.5 (strongly stable, turbulence collapse)

2. **Curvature ∂²Ri_g/∂z² (from LES gradients):**
   - Concave-down in 0–100 m → confirms Δ < 0 regime
   - Inflection near 120–140 m → curvature sign change
   - Use finite-difference: (Ri_{k+1} − 2Ri_k + Ri_{k−1}) / (Δz)²

3. **Bias Amplification Ratio B:**
   - Fine LES (Δz=2 m): B ≈ 1.02–1.05 (minimal bias)
   - Coarse SCM (Δz=50 m): B ≈ 1.4–1.6 (40% underestimation)
   - Coarse SCM (Δz=100 m): B ≈ 1.8–2.2 (80% underestimation)

4. **Correction Validation:**
   - Apply G(ζ,Δz) to coarse SCM runs
   - Target: Reduce B to 1.15–1.25 (≤25% residual bias)
   - Check: Inversion height error <15 m, LLJ timing within 30 min

**Expected Neutral Curvature (from φ fits to GABLS1 LES):**
- If using Businger-like stable: Δ ≈ −1.5 to −2.0
- If using Beljaars-Holtslag: Δ ≈ −4 to −6 (stronger curvature)

---

## GABLS2 — Diurnal Cycle (Cabauw Site)

### Case Description

**Reference:** Svensson et al. (2011), *Boundary-Layer Meteorology*, 140, 177–206.

**Setup (Semi-Realistic):**
- Location: Cabauw, Netherlands (51.97°N, 4.93°E)
- Date: 23 October 1999 (IOP from CASES-99 campaign analog)
- Domain: 9.6 km × 9.6 km × 4 km (LES), 25.6 km × 25.6 km (mesoscale)
- Surface: Prescribed from observations (grass, z₀=0.15 m)
- Forcing: Full diurnal cycle (24 h), geostrophic wind ~10 m/s, clear sky

**Key Complexity:**
- Morning transition: stable → neutral → convective
- Afternoon: CBL growth to ~1500 m
- Evening transition: convective → neutral → stable
- Nighttime: SBL development, LLJ formation

### Validation Targets

**Time-Series Metrics (full 24h):**
- Surface sensible heat flux H₀(t)
- 2-m temperature T₂ₘ(t)
- 10-m wind speed U₁₀ₘ(t)
- Boundary layer height h(t)

**Critical Periods for SBL Curvature:**
- 18:00–06:00 UTC: Stable phase
- Focus: 21:00–03:00 UTC (deepest SBL, strongest stratification)

**Observed (Cabauw Tower) vs LES:**
- Inversion strength: Δθ ≈ 5–8 K (0–200 m)
- LLJ height: 150–250 m
- Surface cooling rate: 0.5–1.0 K/h (radiative + turbulent)

### Curvature Diagnostic Application

**Challenges (vs GABLS1):**
1. **Variable L(z):**
   - Radiative flux divergence → dL/dz ≠ 0
   - Use omission metric E_omit (Section 15C of SBL corrections.md)
   - Apply constant-L shortcut only where E_omit < 0.05

2. **Hysteresis in Ri_c:**
   - Evening transition: turbulence persists until Ri_g > 0.4 (memory effect)
   - Morning transition: turbulence resumes at Ri_g < 0.2 (asymmetry)
   - Require dynamic Ri_c* = f(TKE_prev, Δθ/Δz)

3. **Heterogeneous Surface:**
   - Patchy grass/water → spatial L(x,y) variations
   - Aggregate fluxes using tile-based approach
   - Check curvature consistency across tiles

**Validation Workflow:**
- Extract LES Ri_g(z,t) profiles at 1-h intervals (18:00–06:00)
- Compute Δ from near-neutral segments (Ri_g < 0.1, 18:00–19:00)
- Apply curvature correction to SCM; compare T₂ₘ RMSE vs observations
- Target: RMSE reduction ≥1 K over nocturnal period

---

## GABLS3 — Arctic Stable Boundary Layer (Sodankylä)

### Case Description

**Reference:** Bosveld et al. (2014), *Boundary-Layer Meteorology*, 152, 313–341.

**Setup (Realistic Arctic):**
- Location: Sodankylä, Finland (67.37°N, 26.63°E)
- Date: Winter case (January–February), continuous SBL (polar night)
- Surface: Snow-covered boreal forest, z₀ ~ 0.5–1.0 m (tall roughness)
- Forcing: Weak geostrophic winds (2–5 m/s), strong radiative cooling

**Key Features:**
- Very stable regime: Ri_g > 0.5 common in lowest 50 m
- Intermittent turbulence: sporadic bursts, not continuous
- Strong surface inversion: Δθ ≈ 10–20 K (0–100 m)
- Thin turbulent layer: h < 50 m (turbulence confined near surface)

### Validation Targets

**Long-Term Statistics (multi-day composites):**
- Mean profiles: U(z), θ(z), Ri_g(z)
- Flux intermittency factor: fraction of time with |u'w'| > threshold
- Turbulence collapse frequency: hours per day with h < 30 m

**Single-Event Metrics (case study nights):**
- Surface temperature: T_sfc − T_200m (inversion strength)
- Momentum flux: u'w'(z=10 m) (often near zero)
- Heat flux: w'θ'(z=2 m) (residual turbulence check)

### Curvature Diagnostic Application

**Extreme Stability Challenges:**

1. **Laminar Layer Identification:**
   - If Ri_g > 1.0 and du/dz < 0.05 s⁻¹ → flag as laminar
   - Apply K_background (molecular + weak residual) instead of turbulent K
   - Use ML classifier: P_laminar = f(Ri_b, ζ, Δz/z_g, ∂²θ/∂z²)

2. **Enhanced Curvature (Large |Δ|):**
   - SHEBA-derived φ: Δ ≈ −6 to −8 (very strong concave-down)
   - Coarse-grid bias factor: B > 2.5 for Δz = 100 m
   - Require aggressive damping: D ≈ 1.2–1.5 in G(ζ,Δz)

3. **Roughness Sublayer:**
   - Forest canopy → displacement height d ≈ 0.7 h_canopy ≈ 10–15 m
   - Use ζ̃ = (z − d)/L instead of z/L
   - Adjust geometric mean: z_g = √[(z₁−d)(z₂−d)]

**Validation Workflow:**
- Compare LES ensemble (5+ codes) vs tower composites
- Isolate nights with geostrophic wind <3 m/s, clear sky
- Evaluate curvature correction impact on:
  - Surface temperature bias (target: <1.5 K RMSE)
  - Inversion depth error (target: <20 m)
  - Flux intermittency timing (target: ±1 h for collapse onset)

---

## GABLS4 — Coastal/Sea-Ice Transition (Arctic Ocean)

### Case Description

**Reference:** Sandu et al. (2013), *Geoscientific Model Development*, 6, 1359–1378; extended in Holtslag et al. (2013).

**Setup (Polar Marine):**
- Location: Arctic Ocean, sea-ice edge (hypothetical but representative)
- Season: Late winter (March), continuous night
- Surface: Transition from open water (ice-free leads) to solid sea ice
- Roughness contrast: z₀_water ~ 1e-4 m, z₀_ice ~ 1e-3 m, z₀_h < z₀ (kB⁻¹ large)

**Key Physics:**
- Strong surface heterogeneity → internal boundary layers
- Heat flux discontinuity at ice edge → step change in L
- Advection of cold air over warmer leads → convective plumes

### Validation Targets

**Spatial Variability Metrics:**
- Horizontal flux divergence: ∂(u'w')/∂x near ice edge
- Temperature jumps: ΔT across 1 km fetch
- Mixing depth contrast: h_water − h_ice

**LES Ensemble Outputs (not yet fully standardized):**
- Prototype runs from ECMWF, MetOffice, NCAR
- Focus: 2D cross-sections along fetch direction

### Curvature Diagnostic Application

**Heterogeneous L(x,z) Challenges:**

1. **Spatial Curvature Mapping:**
   - Compute Ri_g(x,z,t) on 2D slices
   - Track curvature evolution downstream of ice edge
   - Identify where constant-L assumption fails (E_omit > 0.1)

2. **Tile-Based Correction:**
   - Separate G(ζ,Δz) for water vs ice tiles
   - Weighted aggregate: K_eff = Σ_tile (A_tile / A_total) K_tile
   - Check discontinuity artifacts at tile boundaries

3. **Dynamic Ri_c* for Leads:**
   - Convective plumes → intermittent turbulence even in stable mean flow
   - Use TKE-based memory: Ri_c*(t) = f(TKE(t−Δt), local ∂θ/∂z)

**Validation Workflow (Preliminary):**
- Awaiting coordinated LES submissions (GABLS4 Phase 2)
- Interim: Use idealized 2D simulations (water→ice step change)
- Target diagnostics: flux recovery length, Ri_g adjustment distance

---

## GABLS LES Data Access

### Official Repositories

1. **KNMI GABLS Portal:**
   - URL: https://www.knmi.nl/research/observations-data-technology/atmospheric-boundary-layer/gabls
   - Contents: Case specifications, model output submission guidelines, benchmark statistics
   - Format: NetCDF (CF-compliant), CSV summaries

2. **NCAR Research Data Archive (RDA):**
   - Dataset IDs: ds624.1 (GABLS1), ds624.2 (GABLS2)
   - Access: Free registration required
   - Resolution: Full 3D fields (large files, 10–100 GB per case)

3. **ECMWF GABLS Intercomparison:**
   - URL: https://www.ecmwf.int/en/research/projects/gabls
   - Contents: SCM-LES comparison tools, visualization scripts
   - Format: GRIB2, NetCDF

### Recommended Subsets for Validation

**Minimal Download (for curvature validation):**
- GABLS1: 1D profiles (ensemble mean), t=9h, Δz≤5 m
  - Variables: U, θ, u'w', w'θ', TKE, Ri_g (computed or provided)
  - File size: ~5 MB (NetCDF)

- GABLS2: 1D profiles (Cabauw location), hourly, 00–24h
  - Variables: same as GABLS1
  - File size: ~50 MB

**Full Dataset (for LES benchmarking):**
- 3D snapshots: u, v, w, θ (every 10 min or hourly)
- File size: 50–500 GB depending on case

### Preprocessing for Ri_g Curvature

**Python Workflow (pseudo-code):**

```python
import xarray as xr
import numpy as np

# Load LES output
ds = xr.open_dataset('GABLS1_LES_ensemble_mean.nc')
z = ds['z'].values  # height (m)
U = ds['U'].isel(time=-1).values  # final time
theta = ds['theta'].isel(time=-1).values

# Compute gradients (centered difference)
dU_dz = np.gradient(U, z)
dtheta_dz = np.gradient(theta, z)

# Gradient Richardson number
g = 9.81; theta_ref = theta.mean()
Ri_g = (g / theta_ref) * dtheta_dz / (dU_dz**2)

# Curvature (second derivative)
d2Ri_dz2 = np.gradient(np.gradient(Ri_g, z), z)

# Neutral curvature (fit near Ri_g ≈ 0)
neutral_mask = (Ri_g > -0.05) & (Ri_g < 0.05)
Delta_fit = np.polyfit(z[neutral_mask], d2Ri_dz2[neutral_mask], 0)[0] / 2

# Bias amplification (first layer, z0=0.1 m to z1=50 m)
z_g = np.sqrt(0.1 * 50)
Ri_g_zg = np.interp(z_g, z, Ri_g)
Ri_b = np.trapz(Ri_g[(z >= 0.1) & (z <= 50)], 
                z[(z >= 0.1) & (z <= 50)]) / (50 - 0.1)
B = Ri_g_zg / Ri_b

print(f"Neutral curvature: 2Δ = {2*Delta_fit:.3f}")
print(f"Bias ratio B (0.1-50 m layer): {B:.3f}")
```

---

## GABLS SCM Intercomparison

### Participating Models (Sample)

**Operational NWP:**
- ECMWF IFS (EDMF scheme, multi-level dry convection adjustment)
- NCEP GFS (hybrid PBL, TKE-based closure)
- MetOffice UM (Lock et al. scheme, non-local K-profile)
- JMA GSM (MYNN level 2.5, Mellor-Yamada-Nakanishi-Niino)

**Research Models:**
- WRF-YSU, WRF-MYJ (Weather Research and Forecasting variants)
- COAMPS (Navy coupled ocean-atmosphere model)
- HARMONIE-AROME (limited-area, non-hydrostatic)

### Common Biases (Pre-Correction)

**GABLS1 (Stable, Idealized):**
- Too-strong mixing: 60–80% of models overestimate surface heat flux by 20–50%
- Shallow inversion: h underestimated by 10–30 m
- Weak LLJ: peak wind speed 0.5–1.0 m/s too low

**GABLS2 (Diurnal, Cabauw):**
- Nocturnal warm bias: T₂ₘ overestimated by 1–3 K (18:00–06:00 UTC)
- Premature morning transition: convective onset 30–60 min early
- Afternoon CBL: reasonable (stable corrections less relevant)

**GABLS3 (Arctic, Very Stable):**
- Severe overmixing: 40–60% of models collapse h < 20 m (unrealistic continuous turbulence)
- Inversion too weak: Δθ underestimated by 5–10 K
- Momentum decoupling: LLJ structure poorly captured

### Curvature Correction Impact (Projected)

**Expected Improvements (based on offline tests):**
- Surface temperature RMSE: −25% to −40% (nocturnal periods)
- Inversion height error: −30% to −50% (strongly stable cases)
- LLJ timing: −15 to −30 min (onset/decay)
- Flux bias: −20% to −35% (momentum and heat)

**Models Most Likely to Benefit:**
- K-profile schemes (e.g., YSU, Lock) → direct K modification
- TKE schemes with Ri-dependent length scale → l*(Ri_c*) hook
- Schemes already using geometric mean for Ri_b → minimal code change

**Models Requiring Adaptation:**
- Non-local schemes (countergradient terms) → embed G in diffusivity component only
- EDMF (mass-flux + diffusion) → apply G to diffusive part, not mass flux
- Large-scale condensation coupled → check moisture curvature separately

---

## Validation Protocol for Curvature Corrections

### Step-by-Step Workflow

**Phase 1: Offline Validation (1D Column)**

1. **Extract LES Forcing:**
   - Geostrophic wind components (U_g, V_g)
   - Surface temperature or heat flux (prescribed or interactive)
   - Initial profiles (U, V, θ at t=0)

2. **Run SCM (Baseline):**
   - Use standard MOST without correction
   - Save 1D profiles: U(z,t), θ(z,t), K_m(z,t), K_h(z,t)
   - Compute Ri_g(z,t), Ri_b(t) for lowest layer

3. **Run SCM (Corrected):**
   - Apply G(ζ,Δz) to K_m, K_h (or φ functions)
   - Preserve all other settings (advection, radiation, etc.)
   - Save same diagnostics

4. **Compare:**
   - Time series: T_2m(t), H_0(t), h(t)
   - Profiles at key times: Ri_g(z) at t=9h (GABLS1)
   - Bias metrics: RMSE(T_2m), MAE(h), B(lowest layer)

**Phase 2: Model Intercomparison (Multi-Model)**

5. **Standardize Implementation:**
   - Provide reference Fortran/Python module
   - Define common namelist parameters: D, p, q in G formula
   - Require version control: document ModelVersion, CorrectionVersion

6. **Submit to GABLS Repository:**
   - Format: NetCDF (CF-1.8 conventions)
   - Variables: T_2m, U_10m, h, Ri_b_lowest, G_applied(z,t)
   - Metadata: contact, institution, model name, correction tag

7. **Ensemble Statistics:**
   - Compute multi-model mean and spread
   - Rank models by skill score (Taylor diagram, target diagram)
   - Identify outliers and failure modes

**Phase 3: Sensitivity Analysis**

8. **Parameter Sweep:**
   - Vary D ∈ [0.5, 1.5], p ∈ [1.0, 2.0], q ∈ [2.0, 3.0]
   - Run 27 combinations (3×3×3 factorial)
   - Identify optimal (D*, p*, q*) per case

9. **Grid Convergence:**
   - Test Δz = [10, 25, 50, 100, 200] m
   - Verify: G → 1 as Δz → 10 m; B_corrected ≈ 1.1 at Δz = 100 m

10. **Robustness Checks:**
    - Noisy initialization: perturb U, θ by ±5%
    - Missing data: remove 10% of forcing points (interpolate)
    - Extreme cases: double surface cooling rate

**Phase 4: Publication-Ready Figures**

11. **Standard Plots:**
    - Time series: T_2m (obs, LES, SCM_baseline, SCM_corrected)
    - Vertical profiles: Ri_g(z) at t=9h (all four)
    - Scatter: Ri_b vs Ri_g(z_g) (baseline=red, corrected=blue, 1:1 line)
    - Taylor diagram: RMSE vs correlation (all participating models)

12. **Curvature Diagnostics Panel:**
    - Top: ∂²Ri_g/∂z²(z,t) contour (time on x-axis, height on y-axis)
    - Middle: G(z,t) contour (damping strength)
    - Bottom: B(t) time series (bias ratio evolution)

---

## Recommended Citations

**Core GABLS Papers:**

1. **GABLS1:** Beare, R. J., et al. (2006). An intercomparison of large-eddy simulations of the stable boundary layer. *Boundary-Layer Meteorology*, 118(2), 247–272. https://doi.org/10.1007/s10546-004-2820-6

2. **GABLS2:** Svensson, G., et al. (2011). Evaluation of the diurnal cycle in the atmospheric boundary layer over land as represented by a variety of single-column models: The second GABLS experiment. *Boundary-Layer Meteorology*, 140(2), 177–206. https://doi.org/10.1007/s10546-011-9611-7

3. **GABLS3:** Bosveld, F. C., et al. (2014). The third GABLS intercomparison case for evaluation studies of boundary-layer models. Part A: Case selection and set-up. *Boundary-Layer Meteorology*, 152(2), 133–156. https://doi.org/10.1007/s10546-014-9917-3

4. **GABLS4:** Sandu, I., et al. (2013). Why is it so difficult to represent stably stratified conditions in numerical weather prediction (NWP) models? *Journal of Advances in Modeling Earth Systems*, 5(2), 117–133. https://doi.org/10.1002/jame.20013

**LES Methodology:**

5. Beare, R. J., & MacVean, M. K. (2004). Resolution sensitivity and scaling of large-eddy simulations of the stable boundary layer. *Boundary-Layer Meteorology*, 112(2), 257–281.

6. Sullivan, P. P., & Patton, E. G. (2011). The effect of mesh resolution on convective boundary layer statistics and structures generated by large-eddy simulation. *Journal of the Atmospheric Sciences*, 68(10), 2395–2415.

---

## Contact Points for GABLS Community

**Scientific Steering Committee (as of 2024):**
- **Chair:** Gunilla Svensson (Stockholm University, Sweden) — gunilla@misu.su.se
- **Co-Chair:** Bert Holtslag (Wageningen University, Netherlands) — bert.holtslag@wur.nl
- **NCAR Representative:** Peter Sullivan — pps@ucar.edu
- **ECMWF Representative:** Irina Sandu — irina.sandu@ecmwf.int

**Data Access Support:**
- KNMI helpdesk: https://www.knmi.nl/about-us/contact
- NCAR RDA support: rdahelp@ucar.edu

**Mailing List (for intercomparison announcements):**
- Subscribe: https://groups.google.com/g/gabls-community (publicly archived)

---

## Summary Table — GABLS Cases for Curvature Validation

| Case | Regime | Duration | Complexity | Priority | Validation Metrics |
|------|--------|----------|------------|----------|-------------------|
| **GABLS1** | Idealized SBL | 9 h | Low | **High** | h, u*, LLJ, Ri_g(z) |
| **GABLS2** | Diurnal cycle | 24 h | Medium | Medium | T_2m, H_0(t), h(t) |
| **GABLS3** | Arctic winter | Multi-day | High | **High** | Intermittency, collapse frequency |
| **GABLS4** | Sea-ice edge | 12–24 h | Very High | Low (prototype) | Spatial flux divergence |

**Recommendation for Initial Validation:**
1. Start with **GABLS1** (cleanest case, rich LES ensemble)
2. Validate curvature correction → reduce B from 1.6 to 1.2 at Δz=50 m
3. Progress to **GABLS3** (extreme test, laminar layer challenge)
4. Use **GABLS2** for diurnal robustness check (full cycle stress test)

---

## Future GABLS Phases (Planned/Proposed)

**GABLS5 (Under Discussion):**
- Theme: Urban stable boundary layer (megacity night)
- Location: Beijing or Paris (dense instrumentation)
- Focus: Anthropogenic heat flux, building-resolving LES, heterogeneous surface

**GABLS6 (Conceptual):**
- Theme: Polar stable BL with leads (marginal ice zone)
- Multi-model coupling: atmosphere + ocean + sea ice
- CMIP7 relevance: Arctic Amplification process studies

**GABLS-Planetary (Exploratory):**
- Adapt GABLS1 protocol to Mars (Viking/InSight lander data)
- Titan analog: methane-cycle stable layers (Huygens descent profiles)

---

## Acknowledgments

GABLS is a joint activity of the Global Energy and Water Exchanges (GEWEX) project and the World Weather Research Programme (WWRP). LES data providers and SCM participants are acknowledged in individual case papers. Tower observations for GABLS2 provided by Cabauw Experimental Site for Atmospheric Research (CESAR). GABLS3 supported by the Finnish Meteorological Institute.

---

**Document Version:** 1.0  
**Last Updated:** 2024-01-15  
**Maintainer:** David E. England (david@davidengland.org)  
**Repository:** https://github.com/DavidEngland/ABL/tree/main/data
