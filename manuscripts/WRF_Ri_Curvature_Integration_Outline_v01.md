# WRF Integration Manuscript Outline (Ri Curvature + Dynamic Ric)

Status: Outline / pre-draft  
Date: 2026-04-22  
Target journals: Monthly Weather Review (primary), JAMC (secondary)

## Working Title

Grid-Aware Richardson-Number Curvature Corrections and Dynamic Critical Richardson Number for Stable Boundary Layers in WRF

## Scope and WRF-Fit Statement

This paper is scoped as an implementation and evaluation study for WRF PBL physics (initially YSU, optional MYNN extension), with emphasis on:

- physically consistent stable-layer corrections
- minimal and restart-safe code changes
- namelist-controlled activation
- reproducible SCM and 3D validation

The implementation anchor is the example module in examples/wrf_integration_example.F90.

## Candidate Authors

- David E. England
- Richard T. McNider
- Arastoo P. Biazar

## 1. Abstract (150-250 words target)

One paragraph each for:

- problem: coarse-grid over-mixing in stable ABL
- method: neutral-preserving damping + dynamic Ric + hysteresis
- platform: WRF integration points and computational cost
- validation: GABLS1 + tower/field + 3D WRF case
- impact: improved h_BL, LLJ timing, and flux errors

## 2. Introduction

### 2.1 Operational Problem

- Stable ABL biases in NWP and Arctic/continental night conditions.
- Why fixed Ric=0.25 is often too rigid.

### 2.2 Grid-Curvature Mechanism

- Concave-down Ri_g and Jensen/bulk-vs-point bias framing.
- Why coarse vertical spacing amplifies over-mixing.

### 2.3 Objectives

- O1: WRF-ready correction that preserves neutral behavior.
- O2: Dynamic Ric* for regime-adaptive transitions.
- O3: Demonstrate skill gains without numerical instability.

## 3. Theory and Formulation

### 3.1 Baseline Relations

- Ri_b, Ri_g, representative-height rationale (z_g and/or z_L as needed).

### 3.2 Neutral-Preserving Stable-Only Damping

- fc = exp[-D (dz/dz_ref)^p (zeta/zeta_ref)^q max(0, B-B_th)]
- Constraints:
- fc -> 1 as zeta -> 0
- fc -> 1 as dz -> small
- fc in [fc_min, 1]

### 3.3 Dynamic Critical Richardson Number

- Ric* as function of Gamma, shear S, TKE memory, curvature proxy kappa.
- Clamped bounds and physical interpretation.

### 3.4 Hysteresis Logic

- suppress threshold: Ri > 1.5 Ric*
- restart threshold: Ri < 0.5 Ric*

## 4. WRF Implementation Details

### 4.1 Files and Schemes

- YSU primary insertion point (module_bl_ysu.F).
- MYNN optional pathway via stability function correction.

### 4.2 Minimal Patch Philosophy

- One new module/use block + local call sites.
- No interface break for existing PBL APIs.
- Optional activation by namelist switches.

### 4.3 Namelist/Runtime Controls (table)

- apply_ri_correction
- use_dynamic_ric
- D_param, p, q, fc_min, B_threshold
- Ric weights and bounds

### 4.4 Numerical Safeguards

- floors for shear and Ri denominators
- bounds for fc and Ric*
- stable-only gating
- no added iterative solver in strong-stable branch

### 4.5 Computational Cost Reporting

- wall-clock overhead (%) in SCM and 3D test
- number/fraction of timesteps using dynamic transitions

## 5. Experimental Design

### 5.1 WRF-SCM Benchmarks

- GABLS1 baseline and corrected runs.
- Optional GABLS3 diurnal/transition test.

### 5.2 Observational/Field Comparison

- SHEBA and/or ARM SGP stable cases (availability-based).

### 5.3 3D WRF Case

- Nocturnal LLJ or Arctic stable event.
- Domain, vertical grid, forcing, physics suite documented.

### 5.4 Sensitivity Matrix

- D in [0.5, 0.8, 1.2]
- with/without dynamic Ric*
- with/without hysteresis

## 6. Evaluation Metrics

- Boundary-layer height h_BL bias and RMSE
- LLJ peak speed/height/timing
- Surface momentum and heat flux error
- Near-surface T2 and wind10 bias
- Regime classification skill (active/intermittent/suppressed)
- Stability/robustness: CFL violations, crashes, outliers

## 7. Results

### 7.1 SCM Skill Changes

- Main table: baseline vs corrected metrics.

### 7.2 Mechanistic Diagnostics

- B ratio and fc distributions by stability class.
- Ric* distribution and transition frequencies.

### 7.3 3D Impacts

- Forecast-relevant fields, not only internal turbulence metrics.

### 7.4 Cost-Benefit

- accuracy gains per compute overhead.

## 8. Discussion

- Why gains occur (curvature correction and adaptive transitions).
- Regimes with strongest/weakest benefit.
- Transferability to MPAS/CMAQ-style diffusion frameworks.
- Limitations: parameter identifiability and case dependence.

## 9. Conclusions

- 3-5 concise findings tied to operational value.
- Recommended default settings for first WRF adoption.

## 10. Reproducibility and WRF-Guideline Checklist

- Provide exact WRF version/commit and compile options.
- Keep corrections optional and backward compatible by default.
- Document all namelist changes and defaults.
- Provide SCM and 3D run scripts and post-processing scripts.
- Include conservative defaults that avoid destabilizing integration.
- Report negative cases (where correction is neutral or degrades skill).

## Figure Plan (initial)

- Fig 1: Schematic of Ri_g curvature and bulk-vs-point bias.
- Fig 2: fc response surface vs (dz, stability).
- Fig 3: WRF insertion points (YSU and MYNN pathways).
- Fig 4: SCM time series (h_BL, LLJ, fluxes).
- Fig 5: Regime transition diagram with Ric* and hysteresis.
- Fig 6: 3D case maps/time-height diagnostics.

## Table Plan (initial)

- Table 1: Parameter defaults and tuned ranges.
- Table 2: SCM benchmark metrics.
- Table 3: 3D case skill and runtime overhead.

## BibTeX Plan

Primary bib file for manuscript build:

- manuscripts/references.bib (existing)

Additional bib sources to merge/check keys from:

- grid.bib
- FOUNDATIONAL.bib
- drafts/references.bib

## Seed References (to cite in draft)

- EnglandMcNider1995BLM
- Businger_1971 or businger1971flux
- BeljaarsHoltslag1991JAMC or beljaars1991flux
- Holtslag_2013 or holtslag2013stable
- Cuxart2006BLM
- Bosveld2014BLM
- Cassano2001MWR
- Grachev_2012b
- Audouin2021JAMES

Note: Key names differ across bib files. Standardize into manuscripts/references.bib before submission formatting.

## Writing Allocation (suggested)

- England: Sections 3-4, first-pass Results text.
- McNider: Introduction framing + operational implications.
- Biazar: Validation sections and 3D case interpretation.

## Next Draft Target

- v02 should include filled Sections 4-7 and full in-text citations.
