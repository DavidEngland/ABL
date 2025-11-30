Step-by-step
1. Specify objective and success metrics (what Ri_c* must predict, tolerances, and operational constraints).
2. Prepare datasets: synthetic MOST/LES profiles + labeled events (collapse/onset) + tower/LES observational folds. Define train/val/test splits by site/time to avoid leakage.
3. Feature engineering: physics-first features (Γ=dθ/dz, S=|dU/dz|, ζ, Δz, z_g, u_*, θ_*, prior flags, time-of-day, season), normalized/anchored scalings, and engineered interaction terms.
4. Symbolic regression setup (PySR or equivalent): choose operator set, complexity penalty, bounds, and physics constraints (monotonicity, bounds on outputs). Use custom loss combining event/timing metrics and regularization.
5. Training protocol: multi-stage (coarse search → refine), cross‑validation by site, include Pareto frontier selection (error vs complexity). Store candidate formulas + metrics.
6. Filter and enforce constraints: clamp expressions to [Ri_min,Ri_max], apply monotonicity checks (numerical tests), and calibrate coefficients via bounded regression (if needed).
7. Validation: evaluate regression error (MAE/RMSE on Ri_c*), event-timing metrics (precision/recall ± tolerance), closed-loop test in column model (stability, energy conservation), and robustness to Δz.
8. Select finalists: prefer simple interpretable expressions meeting safety constraints; produce fallback (fixed Ri_c) and LUT/ONNX exports.
9. Deployment recipe: test inference speed, integrate as callable function in model (with guard rails), logging schema (ModelVersion, DatasetVersion, Thresholds).
10. Reproducibility: containerize training environment, seed RNGs, record hyperparams, publish artifact bundle (formula, tests, small reference dataset).

### [Symbolic Regression Ri_c.md](file:///Users/davidengland/Documents/GitHub/ABL/ML/Symbolic%20Regression%20Ri_c.md)

Guidance: implement symbolic‑regression pipeline to derive an interpretable dynamic critical Richardson number Ri_c*(state). Focus on physics, constraints, and safe deployment.

````markdown
# Symbolic Regression Workflow — Dynamic Critical Richardson Number (Ri_c*)

Goal
- Discover a compact, interpretable formula Ri_c* = f(state) that improves onset/cessation timing of turbulence vs fixed Ri_c (e.g., 0.25).
- Constraints: Ri_c* ∈ [0.15, 2.0], monotonic where physics demands, robust across sites and Δz.

1) Define targets & metrics
- Targets:
  - Primary: Ri_c* (continuous) estimated from labeled events (onset/cessation) or optimization on LES/tower windows.
  - Secondary: boolean event labels derived from turbulence diagnostics (TKE>thr).
- Metrics:
  - Regression: MAE, RMSE on Ri_c* (validation / held-site).
  - Event: precision, recall, F1 for predicted collapse/restore within ±Δt window (e.g., 10 min).
  - Operational: closed-loop test — change in surface bias, energy consistency, no solver crashes.
- Success criteria (example): MAE_Ri_c* < 0.05 and event F1 increase ≥20% vs baseline.

2) Data curation
- Sources: LES ensembles (GABLS-type), tower obs (SHEBA, CASES-99, ARM), synthetic MOST param sweeps.
- Labels:
  - Define turbulence ON/OFF using TKE, vertical flux variance, or energy threshold; produce Ri_c* per event as the Ri where event flips in reanalysis.
  - Aggregate per-layer samples (include Δ