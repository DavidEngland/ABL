# HSNBL Ultraspherical Workflow Status (API + Julia)

## Purpose

This note records what has already been implemented and validated for applying ultraspherical methods to highly stable nocturnal boundary layer (HSNBL) data.

## Implemented Components

### 1. Practical ultraspherical fitter

- Script: julia/ultraspherical_practical_run.jl
- Required input columns: zeta, phi_obs
- Optional input column: time (used for blocked split)
- Outputs:
  - metrics CSV
  - parameter CSV
  - test prediction CSV
  - coefficient CSV
  - comparison plot PNG (if CairoMakie is available)

### 2. Synthetic validation mode

- Built into julia/ultraspherical_practical_run.jl
- Command pattern:
  - julia julia/ultraspherical_practical_run.jl --synthetic output/ultra_synth [noise_frac] [n_samples]
- Status: validated and running.

### 3. Tower preprocessor (CSV + API)

- Script: julia/preprocess_tower_to_ultra_input.jl
- Modes:
  - raw: consumes direct covariances and optional direct gradients
  - two-level: computes gradients from two profile levels
  - api-smear: fetches SmartSMEAR data directly, then preprocesses
- Output schema:
  - time
  - zeta
  - phi_obs
  - phi_m
  - phi_h
  - u_star
  - L
  - quality_pass

## Real-data Demonstration Completed

### Data source used

- SmartSMEAR (HYY station)

### Demonstration dataset generated

- output/hyy_raw_flux_jan2018.csv
- output/hyy_ultra_input_jan2018.csv

### Demonstration fit outputs generated

- output/hyy_ultra_fit_jan2018_metrics.csv
- output/hyy_ultra_fit_jan2018_params.csv
- output/hyy_ultra_fit_jan2018_pred_test.csv
- output/hyy_ultra_fit_jan2018_coeffs.csv
- output/hyy_ultra_fit_jan2018_comparison.png

### Demonstration performance snapshot

From output/hyy_ultra_fit_jan2018_metrics.csv:

- MOST test RMSE: 1.7778
- MOST+ULTRA test RMSE: 0.7994

Interpretation:

- The ultraspherical residual correction substantially improved held-out test error in this first real-data run.

## Notes on Physics Fidelity

This first real-data run is a workflow proof and not final publication-grade calibration. Main reasons:

- Wind-gradient proxy used from two levels with practical assumptions.
- Displacement height was fixed in the demonstration.
- Additional strict quality filtering can be added (night-only windows, flux quality classes, precipitation/rime filtering, fetch filtering).

## Recommended Current Command Patterns

### A. Preprocess from local CSV (two-level mode)

julia julia/preprocess_tower_to_ultra_input.jl input.csv output/prepped.csv 24.0 0.0 --mode=two-level --z1=16.8 --z2=24.0 --stable-only --phi=phi_m

### B. Preprocess directly from SmartSMEAR API

julia julia/preprocess_tower_to_ultra_input.jl HYY output/prepped_api.csv 24.0 0.0 --mode=api-smear --profile-mode=two-level --from=2018-01-01T00:00:00.000 --to=2018-02-01T00:00:00.000 --interval=30 --aggregation=ARITHMETIC --quality=ANY --tv-uw=HYY_EDDY233.u_star --tv-vw=HYY_EDDY233.u_star --tv-wthetav=HYY_EDDY233.u_star --tv-thetav=HYY_EDDY233.u_star --tv-u1=HYY_META.WSU168 --tv-u2=HYY_EDDY233.U --tv-theta1=HYY_META.T168 --tv-theta2=HYY_EDDY233.av_t --z1=16.8 --z2=24.0 --phi=phi_m

Important:

- The tablevariable selections above illustrate structure only.
- Use physically consistent variable selections (momentum flux and buoyancy terms from compatible instruments/processing streams) before final inference.

### C. Fit ultraspherical model

julia julia/ultraspherical_practical_run.jl output/prepped_api.csv output/hyy_ultra_run

## Current Bottom Line

- End-to-end API-to-fit workflow is operational.
- Real data can now be ingested, transformed into zeta and phi_obs, fit with ultraspherical correction, and plotted.
- Next gains will come from stricter variable pairing and QC, not from missing software pieces.
