# WRF Integration Manuscript Outline (Ri Curvature + Dynamic Ric)

Status: Outline / pre-draft  
Date: 2026-04-22  
Target journals: Monthly Weather Review (primary), JAMC (secondary)

## Working Title

Grid-Aware Richardson-Number Curvature Corrections and Dynamic Critical Richardson Number for Stable Boundary Layers in WRF

## Scope and WRF-Fit Statement

This paper is scoped as an implementation and evaluation study for WRF PBL physics (initially YSU, optional MYNN extension), with emphasis on:

- physically consistent stable-layer corrections
- explicit Richardson-number input handling (local gradient vs bulk)
- minimal and restart-safe code changes
- namelist-controlled activation
- reproducible SCM and 3D validation

The implementation anchor is the example module in examples/wrf_integration_example.F90.
Reusable profile helpers are prototyped in examples/module_most_profile_utils.F90.

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
- Explicit bulk-to-effective-gradient mapping when only Ri_b is available:
- Ri_eff = B_ratio * Ri_b, with B_ratio constrained and diagnosed per layer.
- API clarity via Ri-input-kind tagging on each stability-function call.

### 3.2 Neutral-Preserving Stable-Only Damping

- fc = exp[-D (dz/dz_ref)^p (zeta/zeta_ref)^q max(0, B-B_th)]
- Constraints:
- fc -> 1 as zeta -> 0
- fc -> 1 as dz -> small
- fc in [fc_min, 1]
- Stable-only gating: no damping is applied in unstable or neutral branches.

### 3.3 Dynamic Critical Richardson Number

- Ric* as function of Gamma, shear S, TKE memory, curvature proxy kappa.
- Clamped bounds and physical interpretation.

### 3.4 Hysteresis Logic

- suppress threshold: Ri > 1.5 Ric*
- restart threshold: Ri < 0.5 Ric*

### 3.5 Similarity-Function Consistency

- Distinguish MOST similarity functions phi_m(zeta), phi_h(zeta) from Ri-based closure functions f_m(Ri), f_h(Ri).
- Stable f_m, f_h options to document and compare:
- algebraic Ri forms (linear/quadratic)
- MOST-consistent inversion path (invert Ri_g -> zeta, then evaluate phi).
- Clarify that BD71 unstable radicals are used only for zeta < 0 and are not analytically continued into operational SBL closures.

### 3.6 Backward Derivation Track (CBC -> Gegenbauer/Legendre -> Ri Forms)

- Start from unstable MOST forms with explicit coefficients:
- phi_h = (1 - b_h zeta)^(-1/2), phi_m = (1 - b_m zeta)^(-1/4).
- Show near-neutral expansions in eta = -zeta > 0:
- heat branch via central-binomial coefficients (CBC),
- momentum branch via Gegenbauer C_n^(1/4) equatorial values.
- Connect CBC coefficients to even Legendre values P_{2n}(0) and P_n(cos(theta)) representation.
- Document the degenerate identity when b_m = b_h:
- phi_h = phi_m^2 and Ri_g = zeta.
- Document non-degenerate corrections when b_m != b_h and implications for Ri inversion.
- State parameter-aware critical limits:
- UBL momentum limit: Ri_c,UBL,m = -1/b_m,
- UBL heat limit: Ri_c,UBL,h = -1/b_h,
- SBL linear limit: Ri_c,SBL = 1/beta (or component-wise 1/beta_m, 1/beta_h).

## 4. WRF Implementation Details

### 4.1 Files and Schemes

- YSU primary insertion point (module_bl_ysu.F).
- MYNN optional pathway via stability function correction.
- Utility-module pathway for clean reuse:
- phys/module_ri_mapping_utils.F90 (Ri_b -> Ri_g-effective and safeguards)
- phys/module_most_profile_utils.F90 (phi profiles, Ri_g(zeta), zeta inversion)
- examples/module_cbc_legendre_most.F90 (CBC/Gegenbauer recurrences + Ri_c limits)
- examples/driver_cbc_gegenbauer_errors.F90 (tabulated exact-vs-series errors)

### 4.2 Minimal Patch Philosophy

- One new module/use block + local call sites.
- No interface break for existing PBL APIs.
- Optional activation by namelist switches.
- Default behavior remains backward compatible when new switches are off.
- Keep scalar safeguards local (bounds, floors) and avoid branch-dependent side effects.

### 4.3 Namelist/Runtime Controls (table)

- apply_ri_correction
- use_dynamic_ric
- ri_input_kind (gradient|bulk)
- stable_form_selector (LINEAR|QUADRATIC|MOST_BD71)
- D_param, p, q, fc_min, B_threshold
- Ric weights and bounds

### 4.4 Numerical Safeguards

- floors for shear and Ri denominators
- bounds for fc and Ric*
- stable-only gating
- iterative inversion only when MOST_BD71 is selected; algebraic otherwise
- bounded Newton iteration with fallback and clipping for robustness

### 4.5 Computational Cost Reporting

- wall-clock overhead (%) in SCM and 3D test
- number/fraction of timesteps using dynamic transitions
- number/fraction of timesteps using inversion path (MOST_BD71)

### 4.6 Numerical Verification Artifacts (pre-results appendix)

- Include driver-generated tables for phi_h and phi_m approximation errors vs zeta and truncation order N.
- Report convergence envelope and practical switch threshold for series-vs-direct evaluation.
- Include one table using b_m=b_h=16 (canonical) and one with b_m!=b_h (Businger set) to demonstrate non-degenerate behavior.

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
- Ri-input-kind sensitivity (native Ri_g vs mapped Ri_b)
- stable-form sensitivity (LINEAR, QUADRATIC, MOST_BD71)
- polynomial truncation sensitivity for unstable evaluation (N = 4, 8, 12, 16)

## 6. Evaluation Metrics

- Boundary-layer height h_BL bias and RMSE
- LLJ peak speed/height/timing
- Surface momentum and heat flux error
- Near-surface T2 and wind10 bias
- Regime classification skill (active/intermittent/suppressed)
- Stability/robustness: CFL violations, crashes, outliers
- Closure-consistency diagnostics: Ri-input-kind usage and mapped-vs-native Ri differences

## 7. Results

### 7.1 SCM Skill Changes

- Main table: baseline vs corrected metrics.

### 7.2 Mechanistic Diagnostics

- B ratio and fc distributions by stability class.
- Ric* distribution and transition frequencies.
- Frequency and impact of MOST_BD71 inversion branch use.

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
- Fig 2a: CBC/Gegenbauer truncation error vs zeta and N (heat and momentum).
- Fig 3: WRF insertion points (YSU and MYNN pathways).
- Fig 3a: Decision tree for Ri-input-kind handling and stable-form selector.
- Fig 3b: Backward derivation map (MOST -> CBC/Gegenbauer -> Legendre P_n(cos theta) -> Ri_g/Ri_c/Pr).
- Fig 4: SCM time series (h_BL, LLJ, fluxes).
- Fig 5: Regime transition diagram with Ric* and hysteresis.
- Fig 6: 3D case maps/time-height diagnostics.

## Table Plan (initial)

- Table 1: Parameter defaults and tuned ranges.
- Table 1a: Ri-input-kind and stable-form options with default safeguards.
- Table 1b: Parameter-dependent critical limits (Ri_c,UBL,m, Ri_c,UBL,h, Ri_c,SBL) for BD71/Businger/Hogstrom sets.
- Table 2: SCM benchmark metrics.
- Table 3: 3D case skill and runtime overhead.
- Table A1 (appendix): exact vs CBC/Gegenbauer relative errors from driver output.

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

- v02 should include filled Sections 4-7, full in-text citations, and an appendix-ready
  numerical verification package from examples/driver_cbc_gegenbauer_errors.F90.
