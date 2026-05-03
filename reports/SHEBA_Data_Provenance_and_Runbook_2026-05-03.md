# SHEBA Data Provenance and Runbook (2026-05-03)

## Purpose

This document captures where SHEBA data was found, what was downloaded, how it was preprocessed for ultraspherical fitting, and exactly how to rerun the workflow.

## Public Data Source (Confirmed)

Primary archive root:
https://psl.noaa.gov/arctic/sheba/netcdf/

Primary dataset directory used:
https://psl.noaa.gov/arctic/sheba/netcdf/PerssonDatasets/

Files used for this run:
- https://psl.noaa.gov/arctic/sheba/netcdf/PerssonDatasets/main_file6_hd.txt
- https://psl.noaa.gov/arctic/sheba/netcdf/PerssonDatasets/readme_ASFG3.0.txt

Notes:
- Access is anonymous HTTP (no account/login required).
- Data are tab-delimited text with 9999-style missing values.
- The readme documents column definitions for the ASFG tower products.

## Local Retained Raw Snapshot

Saved local copies for reproducibility:
- data/external/sheba/raw/main_file6_hd.txt
- data/external/sheba/raw/readme_ASFG3.0.txt
- data/external/sheba/raw/SHA256SUMS.txt

SHA-256 checksums:
- main_file6_hd.txt: ca97864a7ab372d148b11ba652a74c9ccab30b36b9aa963ed506bad1b51a34e8
- readme_ASFG3.0.txt: f198e9c96c38f68b09e2c3fdd00e91cfb696dd47a6a0dd84b30aedb15e9a6459

## Preprocessing Implementation

Script:
- julia/preprocess_sheba_main.jl

Input:
- main_file6_hd.txt (hourly ASFG data)

Output schema:
- time, zeta, phi_obs

Core formulas:
- Obukhov length:
  L = -(u*^3 * T_K) / (kappa * g * H), with H = hs / (rho_cp)
- Stable coordinate:
  zeta = z_ref / L, where z_ref = sqrt(2.5 * 10.0)
- Two-level momentum stability function:
  phi_m = (kappa * z_ref / u*) * (ws10 - ws2.5) / (10.0 - 2.5)

QC filters applied:
- u* >= 0.05 m/s
- |hs| >= 2.0 W/m2
- ws10 >= ws2.5
- 0 < zeta <= 10
- 0 < phi_obs <= 30

Column normalization caveat handled:
- CSV.jl normalizenames can map u* -> u_ and ws2.5 -> ws2_5
- Script includes robust candidate matching for column names.

## Executed SHEBA Preprocessing Run

Run folder:
- output/runs/SHEBA/20260503_sheba_grachev/

Input CSV produced:
- output/runs/SHEBA/20260503_sheba_grachev/input/sheba_input.csv

Observed preprocessing counts:
- total rows loaded: 8112
- missing u*: 2074
- missing hs: 6
- missing wind/T: 392
- unstable (zeta <= 0): 1468
- QC failures: 1906
- good stable rows: 2266

Quick file check:
- line count of sheba_input.csv: 2267 (header + 2266 records)

## SHEBA Ultraspherical Fit Run

Driver:
- julia/sheba_ultra.jl

Run outputs prefix:
- output/runs/SHEBA/20260503_sheba_grachev/fit/sheba_ultra_grachev

Primary report:
- output/runs/SHEBA/20260503_sheba_grachev/fit/sheba_ultra_grachev_report.md

Metrics from run:
- Grachev baseline test RMSE: 0.3410727610284458
- Grachev+ULTRA test RMSE: 0.3301766276026332
- Relative RMSE gain: 3.19%

Fitted parameters:
- a = 2.49844027746819
- b = 0.6766125007389536
- alpha_xi = 1.7427953628298036
- lambda_star = 0.75
- ridge = 0.0001
- n_ultra = 5
- xi_mode = log
- split_mode = blocked

## Exact Rerun Commands

From repo root:

1) Preprocess SHEBA main file
julia julia/preprocess_sheba_main.jl \
  output/runs/SHEBA/$(date +%Y%m%d)_sheba_grachev/input/sheba_input.csv

2) Fit SHEBA ultraspherical model
julia julia/sheba_ultra.jl \
  output/runs/SHEBA/$(date +%Y%m%d)_sheba_grachev/input/sheba_input.csv \
  output/runs/SHEBA/$(date +%Y%m%d)_sheba_grachev/fit/sheba_ultra_grachev \
  SHEBA

## Additional Public SHEBA Files Worth Future Use

Potential next-step files in the same public directory:
- prof_file_all6_ed_hd.txt (multi-level profile data, better gradient estimates)
- sheba_composite_data.txt
- ASFG_data_10min/ (higher-frequency data, potentially better for HSNBL structure)

## Practical Retention Guidance

- Keep raw snapshots under data/external/sheba/raw with checksums.
- Keep each run immutable under output/runs/SHEBA/<run_id>/.
- If moving to a dedicated repo, copy this file first, then preserve:
  - data/external/sheba/raw/
  - julia/preprocess_sheba_main.jl
  - julia/sheba_ultra.jl
  - output/runs/SHEBA/20260503_sheba_grachev/
