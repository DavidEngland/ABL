# Ultraspherical (Gegenbauer) Status Report for PBL and HSNBL (2026-05-03)

## Executive Status

Current state: active and productive.

What is already demonstrated in this workspace:
- End-to-end preprocessing and fitting pipeline for observed data.
- Multiple baseline modes and xi-map strategies implemented.
- Strong improvement over baseline MOST in several observed cases.
- Public SHEBA ingestion path established and repeatable.

## Implementation Status (Code)

Core scripts:
- julia/ultraspherical_practical_run.jl
- julia/sheba_ultra.jl
- julia/preprocess_tower_to_ultra_input.jl
- julia/preprocess_sheba_main.jl

Notable capabilities now present:
- Baselines: dyer47, linear-fit, ultra-only, most-free, grachev
- Xi-map options: tanh and log
- Blocked split support for chronological holdout
- Hyperparameter search over alpha_xi, lambda_star, ridge, polynomial order
- Unified artifact exports for metrics, params, predictions, formula, model, report, and plots

## Site and Benchmark Results Snapshot

### 1) SHEBA (public NOAA PSL data)

Run:
- output/runs/SHEBA/20260503_sheba_grachev/fit/sheba_ultra_grachev_report.md

Metrics:
- Grachev test RMSE: 0.3410727610284458
- Grachev+ULTRA test RMSE: 0.3301766276026332
- Relative gain: 3.19%

Interpretation:
- Improvement is real but moderate.
- Stable-data-only SHEBA run is functioning and reproducible.

### 2) SMEAR I (Varrio) observed run family

Representative run (linear-fit baseline):
- output/runs/SMEARI/20260503_varrio_mo_linear/fit/varrio_ultra_linearfit_report.md

Metrics (linear-fit):
- MOST test RMSE: 0.8342888174337197
- MOST+ULTRA test RMSE: 0.6769012227541137
- Relative gain: 18.86%

Metrics (dyer47):
- MOST test RMSE: 2.8381327719620715
- MOST+ULTRA test RMSE: 0.9949474093467571

Interpretation:
- Ultraspherical correction is high-value for Varrio conditions.
- Baseline choice strongly affects raw MOST error; correction remains beneficial.

### 3) Hyytiala (January 2018)

Report:
- output/hyy_ultra_fit_jan2018_reported_report.md

Metrics:
- MOST test RMSE: 1.0761533626117377
- MOST+ULTRA test RMSE: 0.6376367216687656
- Relative gain: 40.75%

Interpretation:
- Largest observed relative gain among current documented runs.

### 4) Synthetic validation

Files:
- output/ultra_synth_metrics.csv
- output/ultra_synth_params.csv

Metrics:
- MOST test RMSE: 0.015667122542349088
- MOST+ULTRA test RMSE: 0.009857180354991649

Interpretation:
- Method recovers known structure under controlled noise.

## HSNBL-Focused Assessment

Strengths already shown:
- Grachev baseline is implemented for very stable regimes.
- Log xi-map reduces saturation issues when zeta is large.
- Stable-only SHEBA workflow now works from public raw data to final fit report.

Current limitations:
- SHEBA currently uses two-level gradient proxy from main_file6_hd.txt.
- Improvement at SHEBA is smaller than SMEAR-site gains, suggesting either lower signal-to-noise in the proxy or under-resolved structure.
- Some runs show very large fitted alpha_xi in non-SHEBA modes, which may signal over-compression tuning under specific baselines and data partitions.

Near-term technical priorities:
- Test profile-rich SHEBA source prof_file_all6_ed_hd.txt to improve phi_m gradient fidelity.
- Test 10-minute ASFG files for stronger HSNBL signal.
- Add repeated blocked-fold validation to reduce single-split sensitivity.

## Reproducibility and Data Retention Status

SHEBA provenance and runbook:
- reports/SHEBA_Data_Provenance_and_Runbook_2026-05-03.md

Local retained public raw data with checksums:
- data/external/sheba/raw/main_file6_hd.txt
- data/external/sheba/raw/readme_ASFG3.0.txt
- data/external/sheba/raw/SHA256SUMS.txt

## Should You Spawn a Dedicated Repo?

Recommendation: yes, if you want cleaner collaboration, easier paper support, and less coupling to broad ABL notes.

Decision criteria (current state):
- Enough code exists now to justify isolation (preprocess + fit + reporting).
- Enough run artifacts exist to define baseline CI regression checks.
- Scope is already coherent: ultraspherical closures for PBL/HSNBL.

Suggested target structure:
- src/ (Julia drivers and shared model utilities)
- data/ (small examples only; large data ignored with documented fetch scripts)
- runs/ (optional local output, not committed)
- docs/ (theory notes, runbook, benchmark status)
- tests/ (synthetic recovery tests and smoke tests)

Minimum migration set from this workspace:
- julia/ultraspherical_practical_run.jl
- julia/sheba_ultra.jl
- julia/preprocess_tower_to_ultra_input.jl
- julia/preprocess_sheba_main.jl
- reports/SHEBA_Data_Provenance_and_Runbook_2026-05-03.md
- output/ultra_synth_metrics.csv (or regenerate via CI)

## Practical Next Actions (Time-Bounded)

1) Consolidate one canonical benchmark matrix CSV across SHEBA, SMEAR I, Hyytiala, synthetic.
2) Add one script to run all benchmark cases and emit a single summary markdown.
3) Add one unit/integration test that checks synthetic RMSE threshold stability.
4) Run SHEBA with profile-rich source and compare against current two-level result.

Overall status: the Gegenbauer-ultraspherical analysis for PBL and HSNBL is beyond prototype stage and ready for disciplined packaging into a focused repository.
