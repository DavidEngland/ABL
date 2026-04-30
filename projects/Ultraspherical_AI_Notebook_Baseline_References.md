# Ultraspherical ABL: Baseline References for External AI Notebook Platform

## Purpose

Provide a minimal, high-quality reference set and reproducible notebook scaffold so an external AI notebook system can generate useful first-pass analyses without drifting from your scientific framing.

## Baseline concept references

1. Monin-Obukhov similarity foundation
- Businger et al. (1971)
- Dyer (1974)
- Hegstrom/Hogstrom reanalyses

2. Stable ABL and intermittency constraints
- Grachev et al. stable boundary-layer studies
- McNider and related dynamic critical Richardson literature

3. Orthogonal polynomial and spectral methods
- Classical Gegenbauer/ultraspherical references
- Legendre/Chebyshev Sturm-Liouville texts

4. Arctic amplification process studies
- Surface-based amplification diagnostics
- Boundary-layer coupling/decoupling mechanisms

5. Planetary boundary-layer references
- Mars ABL observational/model papers
- Titan/Venus near-surface or lower-atmosphere transport references

## Notebook baseline structure

Notebook sections (required):

1. Data intake and QC
- Station metadata, units, missingness, flags

2. Baseline MOST implementation
- Compute $\zeta$, $Ri_g$, $\phi_m$, $\phi_h$

3. Ultraspherical representation
- Fit modal coefficients to stability bins

4. Fractional-dimension diagnostic
- At least one estimator with uncertainty

5. Tracer extension
- Humidity as first tracer case

6. Validation and diagnostics
- Residuals, bias by stability class, event-based skill

7. Reproducibility block
- Config snapshot, package versions, random seeds

## Required plots

1. Observed vs modeled transfer functions by stability bin.
2. Modal energy distribution across regimes.
3. Effective dimension vs stability index.
4. Arctic event composite showing mode shifts.

## Required tables

1. Parameter estimates and confidence intervals.
2. Skill metrics for baseline MOST vs ultraspherical model.
3. Sensitivity to truncation order $N$.

## Metadata contract for external AI system

Each run should emit:

1. Data source identifiers and time coverage.
2. Versioned model spec (equations, priors, constraints).
3. Exact solver settings and convergence checks.
4. Structured result summary in JSON.

## Minimal JSON schema (example)

```json
{
  "experiment_id": "string",
  "site_id": "string",
  "time_range": "string",
  "baseline_model": "MOST",
  "spectral_model": "Gegenbauer",
  "truncation_order": 6,
  "parameters": {
    "alpha_m": 0.25,
    "alpha_h": 0.50,
    "b_m": 16.0,
    "b_h": 16.0
  },
  "effective_dimension": {
    "method": "structure_function",
    "value": 2.31,
    "ci95": [2.10, 2.52]
  },
  "metrics": {
    "rmse_baseline": 0.0,
    "rmse_spectral": 0.0,
    "delta_rmse_percent": 0.0
  }
}
```

## First-pass benchmark protocol

1. Run one Arctic station month.
2. Compare baseline vs ultraspherical model with same forcing/QC.
3. Repeat for two additional stations.
4. Accept framework only if improvement generalizes across sites.

## Notes

- Keep the first notebook intentionally small and auditable.
- Scale complexity only after cross-site robustness is demonstrated.