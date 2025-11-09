# Summary (Remote Sensing 2024: High-Resolution Ri_g in Megacity Boundary Layer)

## Instruments & Validation
- Doppler wind lidar (WS vertical res: 20–25 m from 50–3000 m; blind below 50 m).
- Microwave radiometer (temperature/humidity, vertical res 25–40 m in BL).
- 325 m tower (reference, 15 levels).
- Validation metrics (tower vs remote sensing):
  - Temperature RMSE ≤ 1.66 K (R ≥ 0.97).
  - Relative humidity RMSE ≤ 7.9% (R ≈ 0.93–0.96).
  - Wind speed RMSE ≤ 1.45 m/s (R ≈ 0.81–0.9 increasing with height).
  → Confirms high-fidelity input for Ri_g and supports use of analytic curvature forms with observed gradients.

## Observed Vertical Structure (Jan winter period)
- Surface (0–100 m): frequent strong temperature inversion → stable layer (Ri_g > 0.25).
- Active turbulence layer: ~100–500 m daytime (Ri_g < 0; strongest around 300 m).
- Nighttime neutral: ~200–400 m (0 < Ri_g < 0.25).
- Stable above ~500 m; peak stability near ~1000 m (Ri_g > 0.25).
- Thin neutral/weakly unstable layer near ~1500 m (shear/directional changes).
- Additional inversion zone 600–2000 m linked to warm/cold advection (WD southwest events).

## Resolution Sensitivity (3×3 matrix: Spatial 25/50/100 m × Temporal 1 min/30 min/1 h)
- Coarser temporal (1 h) overestimates turbulence intensity & spatial extent (bluer patches).
- Coarser spatial (100 m) overestimates onset height of turbulence by ~50% near surface.
- RMSE vs benchmark (25 m / 1 min) > 1; correlation R < 0.8 for other resolutions.
- Differences grow with altitude: masking by clouds/water vapor & accumulated gradient smoothing.
→ Direct evidence of grid-dependent Ri_g distortion; supports curvature-aware BC and Δz selection criteria.

## Critical Richardson Number
- Operational threshold used: Ri_g = 0.25 for turbulence cessation.
- Literature acknowledges Rc range 0.21–0.25 and persistence of turbulence until Ri_g > 1 if previously active.
→ Suggest dynamic Rc parameterization or normalization s = Ri/Ri_c in closures.

## Temperature Inversion Impacts
- Near-surface inversion drives early stability, reduces effective turbulent depth.
- Inversion aloft (advection type) modulates shear & stratification layering → potential secondary inflection in curvature profile.

## Modeling / Current Study Extraction
1. Use analytic curvature neutral limit 2(α_hβ_h − 2α_mβ_m) to calibrate against observed near-neutral range 200–400 m.
2. Validate curvature sign transitions around 100 m (inversion) and ~500 m (stability onset).
3. Grid dependence experiment: replicate resolution matrix numerically; compare curvature BC vs standard Ri finite difference (expect reduced overestimation at coarse Δz).
4. Introduce dynamic Ri_c selection influenced by inversion strength (e.g., lapse rate metrics).
5. Parameter fitting: prioritize layers 100–500 m for α_m, α_h, β_m, β_h (active turbulence) and 500–1000 m for stability constraints.
6. Mask heights with strong inversion (Δθ_v > 0, sustained) when estimating power-law φ functions to avoid bias.
7. Reporting set: (i) neutral curvature at ~250 m; (ii) first curvature sign change height; (iii) RMSE improvement across synthetic coarse grids.
8. Error analysis: quantify resolution-induced bias factor B_res = Ri_coarse / Ri_fine; test curvature BC reduction B_res' < B_res.

## Proposed Figures (cross-link)
- Comparison: Observed Ri_g profile vs analytic curvature-guided reconstructed profile (fine vs coarse Δz).
- Resolution bias curves: turbulence layer thickness vs Δz & Δt.
- Inversion effect: curvature enhancement factor C(ζ) overlay on temperature lapse/inversion frequency.

## Data Utilization Strategy
- Use validated remote sensing gradients as ground truth for evaluating inversion of ζ(Ri_g).
- Derive empirical Δz thresholds where coarse-grid misclassification (stable vs neutral) exceeds tolerance; calibrate curvature BC.
- Leverage observed active layer thickness (≈400 m) for selecting vertical extent of power-law validity before singular guards.

## Action Items
- Extract tower & lidar gradient pairs for near-neutral curvature numeric verification.
- Implement synthetic aggregation (25 → 50 → 100 m and 1 → 30 → 60 min) to measure error decline with curvature BC.
- Evaluate dynamic Rc: compare fixed 0.25 vs lapse-rate-informed Rc* in Ri-based closure performance.# Summary (Remote Sensing 2024: High-Resolution Ri_g in Megacity Boundary Layer)

## Instruments & Validation
- Doppler wind lidar (WS vertical res: 20–25 m from 50–3000 m; blind below 50 m).
- Microwave radiometer (temperature/humidity, vertical res 25–40 m in BL).
- 325 m tower (reference, 15 levels).
- Validation metrics (tower vs remote sensing):
  - Temperature RMSE ≤ 1.66 K (R ≥ 0.97).
  - Relative humidity RMSE ≤ 7.9% (R ≈ 0.93–0.96).
  - Wind speed RMSE ≤ 1.45 m/s (R ≈ 0.81–0.9 increasing with height).
  → Confirms high-fidelity input for Ri_g and supports use of analytic curvature forms with observed gradients.

## Observed Vertical Structure (Jan winter period)
- Surface (0–100 m): frequent strong temperature inversion → stable layer (Ri_g > 0.25).
- Active turbulence layer: ~100–500 m daytime (Ri_g < 0; strongest around 300 m).
- Nighttime neutral: ~200–400 m (0 < Ri_g < 0.25).
- Stable above ~500 m; peak stability near ~1000 m (Ri_g > 0.25).
- Thin neutral/weakly unstable layer near ~1500 m (shear/directional changes).
- Additional inversion zone 600–2000 m linked to warm/cold advection (WD southwest events).

## Resolution Sensitivity (3×3 matrix: Spatial 25/50/100 m × Temporal 1 min/30 min/1 h)
- Coarser temporal (1 h) overestimates turbulence intensity & spatial extent (bluer patches).
- Coarser spatial (100 m) overestimates onset height of turbulence by ~50% near surface.
- RMSE vs benchmark (25 m / 1 min) > 1; correlation R < 0.8 for other resolutions.
- Differences grow with altitude: masking by clouds/water vapor & accumulated gradient smoothing.
→ Direct evidence of grid-dependent Ri_g distortion; supports curvature-aware BC and Δz selection criteria.

## Critical Richardson Number
- Operational threshold used: Ri_g = 0.25 for turbulence cessation.
- Literature acknowledges Rc range 0.21–0.25 and persistence of turbulence until Ri_g > 1 if previously active.
→ Suggest dynamic Rc parameterization or normalization s = Ri/Ri_c in closures.

## Temperature Inversion Impacts
- Near-surface inversion drives early stability, reduces effective turbulent depth.
- Inversion aloft (advection type) modulates shear & stratification layering → potential secondary inflection in curvature profile.

## Modeling / Current Study Extraction
1. Use analytic curvature neutral limit 2(α_hβ_h − 2α_mβ_m) to calibrate against observed near-neutral range 200–400 m.
2. Validate curvature sign transitions around 100 m (inversion) and ~500 m (stability onset).
3. Grid dependence experiment: replicate resolution matrix numerically; compare curvature BC vs standard Ri finite difference (expect reduced overestimation at coarse Δz).
4. Introduce dynamic Ri_c selection influenced by inversion strength (e.g., lapse rate metrics).
5. Parameter fitting: prioritize layers 100–500 m for α_m, α_h, β_m, β_h (active turbulence) and 500–1000 m for stability constraints.
6. Mask heights with strong inversion (Δθ_v > 0, sustained) when estimating power-law φ functions to avoid bias.
7. Reporting set: (i) neutral curvature at ~250 m; (ii) first curvature sign change height; (iii) RMSE improvement across synthetic coarse grids.
8. Error analysis: quantify resolution-induced bias factor B_res = Ri_coarse / Ri_fine; test curvature BC reduction B_res' < B_res.

## Proposed Figures (cross-link)
- Comparison: Observed Ri_g profile vs analytic curvature-guided reconstructed profile (fine vs coarse Δz).
- Resolution bias curves: turbulence layer thickness vs Δz & Δt.
- Inversion effect: curvature enhancement factor C(ζ) overlay on temperature lapse/inversion frequency.

## Data Utilization Strategy
- Use validated remote sensing gradients as ground truth for evaluating inversion of ζ(Ri_g).
- Derive empirical Δz thresholds where coarse-grid misclassification (stable vs neutral) exceeds tolerance; calibrate curvature BC.
- Leverage observed active layer thickness (≈400 m) for selecting vertical extent of power-law validity before singular guards.

## Action Items
- Extract tower & lidar gradient pairs for near-neutral curvature numeric verification.
- Implement synthetic aggregation (25 → 50 → 100 m and 1 → 30 → 60 min) to measure error decline with curvature BC.
- Evaluate dynamic Rc: compare fixed 0.25 vs lapse-rate-informed Rc* in Ri-based closure performance.