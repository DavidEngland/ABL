# HSNBL Ultraspherical Next Actions

## Goal

Move from workflow proof to robust, publishable HSNBL parameter estimation and diagnostics.

## What Can Be Done Next

### 1. Build station presets for API ingest

Add named presets (for example HYY_eddy24_profile168) so users only set:

- station label
- date range
- z and d
- phi target

Outcome:

- One-line repeatable commands
- Fewer mapping errors in tablevariable flags

### 2. Add metadata-assisted variable auto-discovery

Use SmartSMEAR variable metadata endpoint to propose candidate variables for:

- momentum terms
- buoyancy terms
- wind and temperature profile levels

Outcome:

- Reduced manual setup
- Better portability to additional Finland/EU sites

### 3. Implement stricter HSNBL quality control

Recommended additional filters:

- stable-only: zeta > 0
- nighttime subset (local time windows)
- minimum u_star threshold
- precipitation/surface wetness exclusion
- quality-class filtering where available (for example Qc_* fields)
- outlier rejection in phi_m and phi_h tails

Outcome:

- Lower noise in fitted coefficients
- More defensible physical interpretation

### 4. Add bootstrapped uncertainty for parameters

Run blocked bootstrap resampling in zeta-time blocks and export:

- confidence intervals for baseline MOST parameters
- confidence intervals for ultraspherical coefficients
- spread of test RMSE and MAE

Outcome:

- Quantified uncertainty for manuscript reporting

### 5. Add standard diagnostics panels

In addition to the current comparison plot, generate:

- residual vs zeta
- binned bias vs zeta
- train vs test distribution checks
- coefficient spectrum decay with mode index

Outcome:

- Faster model vetting
- Better communication in papers and talks

### 6. Add cross-site transfer tests

Train on one period/site and test on another:

- seasonal transfer
- site transfer
- year-to-year transfer

Outcome:

- Demonstrates generalization, not just fit quality

## Suggested Implementation Order

1. Station presets
2. QC profile and flags
3. Diagnostics panels
4. Bootstrap uncertainty
5. Cross-site transfer experiments

## Minimal Publishable Pipeline (Target)

1. API ingest with preset and fixed variable map
2. Deterministic QC profile for HSNBL
3. Blocked-split fit and metrics
4. Bootstrap confidence intervals
5. Standard figures and tables exported automatically

## Example Near-Term Runbook

### Step A. Ingest and preprocess

julia julia/preprocess_tower_to_ultra_input.jl HYY output/hyy_qc.csv 24.0 0.0 --mode=api-smear --profile-mode=two-level --from=2018-01-01T00:00:00.000 --to=2018-03-01T00:00:00.000 --interval=30 --aggregation=ARITHMETIC --quality=ANY --tv-uw=<TABLE.VAR> --tv-vw=<TABLE.VAR> --tv-wthetav=<TABLE.VAR> --tv-thetav=<TABLE.VAR> --tv-u1=<TABLE.VAR> --tv-u2=<TABLE.VAR> --tv-theta1=<TABLE.VAR> --tv-theta2=<TABLE.VAR> --z1=16.8 --z2=24.0 --stable-only --phi=phi_m

### Step B. Fit ultraspherical model

julia julia/ultraspherical_practical_run.jl output/hyy_qc.csv output/hyy_ultra_qc

### Step C. Archive outputs for analysis

Keep together:

- metrics
- parameters
- coefficients
- test predictions
- figures

## Decision Criteria for "Ready to Publish"

- Stable improvement over MOST in blocked test RMSE across multiple months
- Physically plausible baseline parameters and residual structure
- Coefficient behavior robust under QC and bootstrap resampling
- Comparable results across at least two independent periods or sites
