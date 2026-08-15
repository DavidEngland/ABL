# Reviewer Comments: Summary and Response Plan (McNider & Biazar)

## 1) Executive Summary of Reviewer Themes
- Scope/novelty: Clarify what is new beyond MOST and standard φ-functions.
- Physical interpretation: Make curvature of Ri_g and Ri-based closures intuitive.
- Mathematical rigor: Show exact forms, limits (neutral/critical), and domains.
- Validation: Add quantitative tests (LES/field), error bars, and sensitivity.
- Grid dependence: Demonstrate boundary-condition utility at coarse Δz.
- Notation/consistency: Distinguish φ, ψ_m, Ri_g vs bulk Ri; equal-β vs general case.
- Reproducibility: Provide data provenance, code, and parameter identifiability.
- Figures/readability: Improve axis choices, neutral limit checks, and legends.
- Related work: Position vs Businger–Dyer, Katul et al., Gryanik et al.

## 2) Authors’ Responses (What was improved)
- Theory: Provided closed-form curvature
  d²Ri_g/dζ² = F(ζ)[2V_log + ζ(V_log² − W_log)], with expanded self-contained form and neutral limit 2(α_hβ_h − 2α_mβ_m).
- Ri-based closures: Derived f_m(Ri), f_h(Ri): power laws for equal-β; cubic/algebraic inversion for integer-exponent cases; near-neutral series for general parameters.
- Bias analysis: Quantified geometric-mean vs arithmetic-height bias in C_D and implications for Ri.
- Critical Ri: Defined Ri_c linkage via F(ζ_c); gave normalized forms in s = Ri/Ri_c.
- HS utility: Outlined HS-assisted series/inversion with error control and branch handling.

## 3) Remaining Gaps → Next-Paper Action Items
- Validation suite:
  - LES: stable/unstable runs; compute Ri_g curvature vs analytic; RMSE, bias, R².
  - Field: canopy/urban/open terrain; stratification bins; bootstrap CIs.
- Grid dependence: Controlled Δz experiments; show reduction in near-surface Ri bias using curvature BC.
- Identifiability: Profile-fitting experiments for (α_m,β_m,α_h,β_h,Ri_c); profile-only vs flux-assisted; sloppiness analysis.
- Robustness: Domain guards (1−βζ>0), critical behavior when 2α_m≈α_h, parameter regularization.
- Reproducibility: Open scripts, config files, and data links; seed control; unit tests.
- Figures: Neutral-limit overlays, curvature sign-change heights, sensitivity heatmaps.

## 4) Plan for Contribution in New Paper
- Thesis: A curvature-aware boundary condition and Ri-based shear closures materially reduce grid-dependent errors and improve stability diagnostics in the ASL.
- Core contributions:
  1) Exact curvature formula with neutral/critical behavior and domains.
  2) Ri-based shear closures with equal-β closed forms and cubic inversion cases.
  3) Grid-dependence mitigation via curvature BC at the first level.
  4) HS-assisted series/inversion with error bounds for fast/stable evaluation.
- Evidence: LES + field validation, Δz sweeps, identifiability and robustness analyses.
- Artifacts: Code, datasets, notebooks, and parameter tables.

## 5) Action Matrix
- Theory: finalize cubic/degree-n inversions and edge-case analysis (2α_m≈α_h).
- Data/Experiments: assemble LES cases; curate field datasets; define QC.
- Metrics: RMSE/bias for φ_m, φ_h, Ri_g curvature; Δ bias vs Δz; AIC/BIC for fits.
- Reproducibility: repo structure, tests, containers; data DOIs.
- Writing: figures (curvature, s=Ri/Ri_c power laws, Δz sensitivity), methods narrative.

# Response to Reviewer Comments (Mathematical Clarifications)

## 1. Novelty Clarification
We supply: (i) closed-form curvature \( \partial_{\zeta}^2 Ri_g\) for general power-law φ; (ii) neutral / inflection coefficients (\(\Delta,c_1\)); (iii) ζ–Ri inversion series + single-step Newton eliminating iterative L–ζ loops; (iv) controlled error bounds (binomial / HS) for truncation.

## 2. Derivation Transparency
Key identity:
\[
Ri_g=\zeta (1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{2\alpha_m},\quad
\partial_{\zeta}^2 Ri_g = F\Big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\Big].
\]
All higher derivatives expressible by polynomials in \(V_{\log},W_{\log}, dW_{\log}/d\zeta\), supplied in supplementary (Section 12 main ms).

## 3. Boundary Condition Rationale
Finite-difference surface layering missing curvature term yields O(Δz) bias in Ri_g; imposing analytic second derivative restores O(Δz^2) accuracy. Neutral coefficient is a simple linear combination of fitted φ parameters, enabling direct plausibility checks.

## 4. Ri-Based Closure Advantage
Eliminates transformation back to ζ each sub-iteration: ζ(Ri) series:
\[
\zeta = Ri_g - \Delta Ri_g^2 + \left(\tfrac32\Delta^2-\tfrac12 c_1\right) Ri_g^3 + O(Ri_g^4).
\]
Practical: residual after one Newton step ~ O(Ri_g^4).

## 5. Sensitivity / Identifiability
Curvature neutral coefficient \(2\Delta\) collapses sensitivity to *one* scalar; simultaneous small |Δ| and small |c₁| produce ill-conditioning (flagged numerically). Provide Jacobians for inversion.

## 6. Singular Behavior Handling
Algebraic blow-up exponent \(\alpha_h+1\) (or \(\alpha_m+1\) if momentum singular first) informs ζ upper cut. Domain guard: if \(1-\beta_{(\cdot)}\zeta < 0.1\) switch to asymptotic outer/inner composite.

## 7. Error Bounds
Remainder in φ series:
\[
|R_{N+1}| \le \frac{(\beta\zeta)^{N+1}}{(N+1)(1-\beta\zeta)}.
\]
Propagated curvature truncation < ε ensures Ri bias contribution < O(ε ζ). We tabulate N for target ε=10^{-4}.

## 8. Figures (Planned)
F1 Analytic vs finite-difference curvature (LES).  
F2 ζ–Ri inversion residual distribution.  
F3 Neutral curvature histogram across sites.  
F4 Sensitivity of curvature to Δ near zero crossing.  
F5 Δz bias reduction (baseline vs analytic BC).  

## 9. Reproducibility
Module will expose: curvature(ζ, params), invert_Ri(Ri, params), error_bound(ζ,N,params). Versioned notebooks demonstrate matches for test parameter sets.

## 10. Summary to Reviewers
The curvature expression is *not* a cosmetic derivative; it is the enabling analytic ingredient for consistent Ri-based closure on coarse grids, with quantifiable accuracy and stability controls.

# Response to Reviewer Comments — References to Equations (21), (22), and Parameters

Citations
- Eq. (21): see Appendix paste target in JAMC_format_Grid_Dependence_V15_10-31-2025.md.
- Eq. (22): see Appendix paste target in JAMC_format_Grid_Dependence_V15_10-31-2025.md.
- Parameter table: see Appendix table in JAMC_format_Grid_Dependence_V15_10-31-2025.md.

Once pasted, cross-reference here with equation numbers and parameter names to support the responses.

# Response to Reviewer Comments (Placeholder)

- We added a multi‑profile curvature analysis (BD, QSBL, CB‑type, DTP). All use the same closed curvature form via V_log and W_log.
- New reporting includes neutral curvature (2Δ), inflection height, ζ‑range, and physical curvature via 1/L².
- Grid‑dependence section repeats the tests for each profile with consistent metrics.
