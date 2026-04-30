For Earth’s polar regions, the key is that surface-atmosphere coupling is dominated by stable stratification, longwave radiation, and heterogeneity in ice/snow state, while some “midlatitude intuitions” are much less important.

**Most Important Factors**
1. Surface energy balance partition
- Net longwave loss, low sun angle, cloud phase, and snow/ice emissivity dominate near-surface temperature structure.

2. Stable boundary-layer structure and intermittency
- Very stable, shallow, intermittently turbulent layers strongly control momentum, heat, moisture, and tracer exchange.

3. Sea-ice and snow properties
- Ice concentration, leads/polynyas, snow depth, grain metamorphism, albedo, and thermal conductivity set lower-boundary forcing.

4. Surface heterogeneity and fetch transitions
- Ice-ocean-land transitions and roughness contrasts create strong local gradients and nonstationary flux footprints.

5. Moisture phase-change physics
- Sublimation/deposition, fog/diamond dust, mixed-phase clouds, and cloud-radiative feedbacks strongly affect stability.

6. Synoptic forcing and advection
- Warm-air intrusions, pressure-gradient events, and frontal passages can rapidly reset local ABL regimes.

7. Vertical resolution near the surface
- Coarse lowest layers miss sharp inversions and flux-gradient curvature, causing large bias in stable conditions.

8. Time-scale competition
- Turbulence, radiation, and chemistry can occur on similar timescales; this matters for ozone and other reactive tracers.

**Often Less Important (or Secondary) in Polar ABL Context**
1. Deep daytime convective mixing as a persistent control
- Important seasonally/episodically, but much less dominant than in lower latitudes.

2. Canonical neutral-layer assumptions as defaults
- Neutral assumptions are frequently wrong in polar night and cold-season conditions.

3. Constant/universal transfer coefficients
- Fixed coefficients usually underperform because regime shifts are frequent.

4. Single “global” tracer exponent
- Tracer behavior is often regime- and species-dependent; one exponent is usually too rigid.

5. Midlatitude roughness analogs without ice-state dependence
- Roughness and exchange over sea ice/snow differ substantially from temperate land/ocean assumptions.

6. Ignoring intermittency
- Treating turbulence as continuous misses key burst-driven transport in very stable layers.

## Use/No-Use Checklist for Ultraspherical Harmonics in Polar ABL

Use this as a pre-run gate. If most answers in a row are Yes, proceed.

| Decision item | Use ultraspherical now? | Why |
|---|---|---|
| Need compact representation across stability regimes (unstable to very stable) | Yes | Orthogonal modes separate low-order mean structure from high-order intermittency |
| Multi-station fusion (tower, buoy, reanalysis) needed | Yes | Fit local modal coefficients first, then map coefficients regionally/globally |
| You need interpretable diagnostics for Arctic amplification episodes | Yes | Mode-energy shifts are clearer than raw parameter drift |
| Strong surface heterogeneity (ice-ocean-land transitions) present | Yes | Modal representation handles nonstationary footprint effects better than single fixed coefficients |
| Very sparse observations with no regime diversity | Maybe | Start with low-order modes only and strong regularization |
| Dominant process is unresolved chemistry/microphysics with weak transport constraints | Not yet | Add process model first; transport-only modes can be misleading |
| Need fast surrogate for sensitivity/ensemble tests | Yes | Truncated mode sets create efficient reduced-order emulators |
| Goal is only a quick neutral-layer fit | Not necessary | Classical low-parameter MOST fit is often sufficient |

## Practical Implementation Guide (Code, Results, Graphics)

## Step 1: Baseline physics run (existing code)

1. Use Fortran baseline similarity functions from [code/most.f](code/most.f).
2. Keep a fixed benchmark parameter set (Dyer or Hogstrom) for comparison.
3. Export baseline profiles and residuals versus observations.

## Step 2: Julia prototyping layer (existing code)

1. Use [julia/MOSTProfiles.jl](julia/MOSTProfiles.jl) to generate and invert profile families.
2. Use [julia/SCMSkeleton.jl](julia/SCMSkeleton.jl) for single-column experiments.
3. Add ultraspherical residual fit in Julia first (fast iterate, easier plotting).

## Step 3: Ultraspherical closure insertion

1. Fit baseline MOST parameters for each regime.
2. Map stability coordinate to bounded space, then fit Gegenbauer coefficients.
3. Keep truncation order low at first (for example, N=2,4,6 cross-validation).
4. Reconstruct corrected transfer functions and compare skill against baseline.

## Step 4: Required output figures

1. Observed vs baseline MOST vs ultraspherical-corrected transfer function.
2. Mode coefficient spectra by regime class (unstable/neutral/stable/very stable).
3. Event timeline: low-order vs high-order modal energy during Arctic transitions.
4. Multi-site map of selected mode coefficients (or mode energy fraction).

## Step 5: Minimal success criteria

1. Lower test RMSE than baseline in at least one stable regime.
2. Improved behavior during intermittency events (not just mean conditions).
3. Coefficients remain physically interpretable and stable under resampling.

## FORTRAN and Julia Libraries: What to use

Below are practical options by language for immediate development.

| Task | Fortran options | Julia options |
|---|---|---|
| Orthogonal polynomials / special functions | self-implemented recurrence (recommended for portability), SLATEC-style routines if available in your stack | SpecialFunctions.jl for gamma/beta support; evaluate Gegenbauer via recurrence or package helper |
| Nonlinear fitting / inversion | MINPACK wrappers or custom Newton/Levenberg routines | LsqFit.jl, Optim.jl, NLsolve.jl |
| Linear algebra for modal solve | LAPACK/BLAS | LinearAlgebra (stdlib), IterativeSolvers.jl |
| NetCDF I/O for met data | netcdf-fortran | NCDatasets.jl |
| Dataframes and tabular workflow | custom arrays + CSV parse | DataFrames.jl, CSV.jl |
| Plotting and graphics | write CSV for external plotting | CairoMakie.jl or Plots.jl |

## Recommendation for fastest path in this repo

1. Prototype and visualize in Julia using [julia/MOSTProfiles.jl](julia/MOSTProfiles.jl).
2. Once stable, port the final low-order correction kernel to Fortran near [code/most.f](code/most.f).
3. Validate in SCM before any WRF/MYNN integration using [julia/SCMSkeleton.jl](julia/SCMSkeleton.jl) and then [code/module_bl_mynn.F90](code/module_bl_mynn.F90).

## 2-week practical sprint

1. Days 1-3: Baseline fit and residual diagnostics at one station.
2. Days 4-6: Ultraspherical residual fit and cross-validation.
3. Days 7-9: Add second station, compare coefficient transferability.
4. Days 10-12: Build four required figures and summary table.
5. Days 13-14: Freeze configuration and draft methods/results memo.

If desired, the next concrete step is a Julia script that ingests a station CSV and produces all four figures in one run.

## Ready-to-run script in this repo

Starter script created at [julia/ultraspherical_practical_run.jl](julia/ultraspherical_practical_run.jl).

Expected input CSV columns:

1. zeta
2. phi_obs

Suggested Julia package setup (one time):

1. using Pkg
2. Pkg.add(["CSV", "DataFrames", "LsqFit", "CairoMakie"])

Example run:

1. julia julia/ultraspherical_practical_run.jl data/station_phi_m.csv output/ultra_demo

Outputs:

1. output/ultra_demo_metrics.csv
2. output/ultra_demo_params.csv
3. output/ultra_demo_pred_test.csv
4. output/ultra_demo_comparison.png (if CairoMakie is installed)
