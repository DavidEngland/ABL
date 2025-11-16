# Publication Topics and Strategy — McNider, Biazar & England

Strategic focus: leverage curvature-aware MOST corrections, grid-dependence diagnostics, and Richardson number closures for stable/polar boundary layers. Build from foundational theory → operational implementation → climate/planetary extensions.

---

## Tier 1: Core Methodology Papers (Near-Term, 6–12 months)

### 1A. Grid-Dependent Curvature Correction for Stable Boundary Layers
**Title (draft):** "Curvature-Aware Stability Function Corrections for Coarse-Grid Atmospheric Models: Preserving Neutral Invariance in the Stable Boundary Layer"

**Target Journal:** *Boundary-Layer Meteorology* (primary) or *Journal of the Atmospheric Sciences*

**Core Contribution:**
- Analytic derivation of ∂²Ri_g/∂ζ² and neutral curvature invariant Δ = α_h β_h − 2α_m β_m.
- Grid damping factor G(ζ,Δz) preserving 2Δ while reducing coarse-grid bias.
- Validation: tower (ARM NSA, SHEBA) + LES (GABLS suite) demonstrating 40%+ curvature error reduction.
- Implementation: drop-in modifier for operational K-closure schemes.

**Why BLM:** Authoritative venue for MOST refinements; direct readership overlap with model developers.

**Extension Hooks:**
- Variable L(z) omission metric E_omit.
- Geometric mean height justification for layer reconstruction.
- Comparison with Beljaars–Holtslag, Cheng–Brutsaert forms.

---

### 1B. Richardson Number Series Inversion and Pole-Free Closures
**Title (draft):** "Direct Richardson Number Closures via Near-Neutral Series Expansion and Rational Approximation: Eliminating Iterative Obukhov Length Solvers"

**Target Journal:** *Monthly Weather Review* (operational focus) or *Journal of Applied Meteorology and Climatology*

**Core Contribution:**
- ζ(Ri) series (quadratic + cubic) with Newton refinement achieving machine precision in 1–2 steps.
- Padé [1/1] rational form for f_m(Ri), f_h(Ri) avoiding hard poles.
- Quadratic SBL (Q-SBL) surrogate matched on neutral coefficients (a_m, b_m, a_h, b_h, Δ, c_1).
- Cost analysis: eliminate ζ iteration overhead in time-stepping loops.

**Why MWR/JAMC:** Emphasizes operational utility; algorithm details appreciated by forecasting community.

**Extension Hooks:**
- Hybrid ζ/Ri switching criteria.
- Dynamic turbulent Prandtl Pr_t(Ri) = 1 + a₁Ri + a₂Ri².
- Regularized power-law (RPL) and variable-exponent (VEXP) forms.

---

### 1C. Surface Roughness, Geometric Mean Heights, and Drag Coefficient Bias
**Title (draft):** "Geometric Mean Heights in Logarithmic Boundary Layers: Implications for Bulk Transfer Coefficients and Richardson Number Reconstruction"

**Target Journal:** *Quarterly Journal of the Royal Meteorological Society* or *Journal of Geophysical Research: Atmospheres*

**Core Contribution:**
- Mathematical proof: arithmetic mean height biases C_D and Ri_b when profiles are logarithmic.
- z_g = √(z₁z₂) minimizes layer-reconstruction error for power-law / log structures.
- Heat vs momentum roughness (z₀_h ≠ z₀): kB⁻¹ implications for stable nights, polar regions.
- Practical workflow: calibrating z₀, z₀_h from tower + remote sensing; curvature consistency checks.

**Why QJRMS/JGR-A:** International scope; surface-layer parameterization audience includes climate modelers.

**Extension Hooks:**
- Urban / canopy displacement height d treatment.
- Sea-ice / snow surface renewal and cool-skin adjustments.
- Multi-level drag coefficient vertical profiles.

---

### 1D. Dynamic Critical Richardson Number and Hybrid MOST/Ri Closures
**Title (draft):** "Adaptive Regime Transitions in Stable Boundary Layers: A Dynamic Critical Richardson Number Framework for Hybrid MOST/Ri Closures"

**Target Journal:** *Monthly Weather Review* or *Journal of Applied Meteorology and Climatology*

**Core Contribution:**
- Dynamic $Ri_c^*$ formulation informed by inversion strength, shear, TKE memory
- Seamless MOST→Ri transition with curvature-aware blending
- Eliminates iterative L solver overhead in strong stability (>40% cost reduction)
- Validation: SHEBA, ARM SGP, GABLS LES; regime classification accuracy >85%

**Why MWR/JAMC:** Operational focus; algorithm suitable for operational NWP

**Extension Hooks:**
- Hysteresis for intermittent turbulence
- Machine-learned blend functions
- Coupling with TKE-based schemes (IFS, ECMWF)

---

## Tier 2: Application & Validation Papers (12–24 months)

### 2A. Arctic / Polar Stable Boundary Layer Performance
**Title (draft):** "Curvature-Corrected MOST Functions for Arctic Winter: Reducing Mixing Biases in Persistent Stable Layers"

**Target Journal:** *Journal of Geophysical Research: Atmospheres* or *Climate Dynamics*

**Core Contribution:**
- Apply 1A framework to SHEBA, Barrow/Utqiaġvik, Polarstern tower data.
- Document systematic Ri_b < Ri_g(z_g) bias in 100 m first layers → overmixing → warm-biased inversion.
- Show curvature correction reduces surface temperature RMSE by 1–2 K over multi-day stable episodes.
- Link to Arctic Amplification: improved inversion representation → better longwave trapping, albedo feedback timing.

**Why JGR-A/Climate Dynamics:** High-impact polar process focus; climate relevance draws citations.

**Extension Hooks:**
- Sea-ice concentration / melt pond influence on z₀, z₀_h.
- Permafrost carbon feedback sensitivity to SBL mixing depth.
- AO/NAO phase modulation via improved polar BL.

---

### 2B. Megacity Urban Boundary Layer (Remote Sensing Integration)
**Title (draft):** "High-Resolution Gradient Richardson Number Profiles from Lidar-Radiometer Fusion: Grid-Dependence Diagnostics for Urban Stable Layers"

**Target Journal:** *Atmospheric Measurement Techniques* (instrumentation emphasis) or *Boundary-Layer Meteorology* (process emphasis)

**Core Contribution:**
- Validate geometric mean + curvature framework using Remote Sensing 2024 megacity dataset (325 m tower + Doppler lidar + microwave radiometer).
- Resolution matrix (25/50/100 m spatial × 1 min/30 min/1 h temporal): document Ri_g bias amplification.
- Show curvature-aware BC reduces coarse-grid turbulence overestimation by 35–50%.
- Identify inversion-layer secondary inflection points in ∂²Ri_g/∂z².

**Why AMT:** Method paper for instrument synergy and retrieval validation; citable by remote sensing community.

**Extension Hooks:**
- Urban heat island inversion structure and curvature signatures.
- Dynamic Ri_c selection informed by lapse-rate metrics.
- Machine learning surrogate for real-time curvature estimation.

---

### 2C. Nocturnal Low-Level Jets and Decoupling
**Title (draft):** "Curvature-Driven Stability Transitions in Nocturnal Low-Level Jet Formation: A Case Study Analysis with Curvature-Aware Closures"

**Target Journal:** *Journal of the Atmospheric Sciences* or *Boundary-Layer Meteorology*

**Core Contribution:**
- Apply curvature diagnostics to classic Great Plains LLJ cases (ARM SGP, Kansas).
- Show jet nose height correlates with inflection ζ_inf in Ri_g(ζ) where ∂²Ri_g/∂ζ² = 0.
- Demonstrate that preserving 2Δ improves timing of turbulence collapse → better LLJ onset prediction.
- Compare standard MOST, Q-SBL, and curvature-corrected MOST in 1D column mode.

**Why JAS/BLM:** Process-oriented; links classic SBL feature to new diagnostic framework.

**Extension Hooks:**
- Ageostrophic wind component evolution and curvature feedback.
- Shear instability triggering: Ri_g crossing 0.25 vs dynamic threshold.
- Sensitivity to land-surface heterogeneity and advection.

---

## Tier 3: Extended Framework & Theory (18–36 months)

### 3A. Variable Obukhov Length L(z) Mapping and Omission Metrics
**Title (draft):** "Height-Dependent Obukhov Length in Strongly Stratified Layers: Omission Metrics and Curvature Mapping"

**Target Journal:** *Boundary-Layer Meteorology* or *Journal of Fluid Mechanics* (if focusing on mathematical rigor)

**Core Contribution:**
- Full chain-rule derivation: ∂²Ri_g/∂z² = (dζ/dz)² ∂²Ri_g/∂ζ² + (d²ζ/dz²) ∂Ri_g/∂ζ.
- Omission metric E_omit quantifies when constant-L shortcut is valid.
- Structured L(z) forms: affine, power-law, exponential → analytic ζ'(z), ζ''(z).
- Application: very stable SBL with strong flux divergence (radiative cooling gradient).

**Why BLM/JFM:** Formal mathematical contribution; extends MOST applicability range.

**Extension Hooks:**
- Flux divergence parameterizations and L(z) closure.
- LES-derived L(z) profiles as benchmark.
- Optimal layer thickness Δz(z) informed by E_omit.

---

### 3B. Unified Multi-Profile Pack (Comparative Assessment)
**Title (draft):** "A Unified Multi-Profile Framework for Monin–Obukhov Stability Functions: Comparative Curvature, Inversion Accuracy, and Operational Robustness"

**Target Journal:** *Geoscientific Model Development* (GMD) or *Monthly Weather Review*

**Core Contribution:**
- Catalog: power-law (Businger–Dyer), Cheng–Brutsaert, Beljaars–Holtslag, Gryanik (SHEBA), Q-SBL, RPL, VEXP, DTP, NLM, URC.
- Benchmark metrics: neutral curvature 2Δ, inflection height ζ_inf, Ri→ζ inversion accuracy, pole proximity guard.
- Standardized test suite: GABLS LES, ARM tower, synthetic profiles.
- Recommendations by use case: operational NWP, climate long integration, polar extremes, urban heterogeneity.

**Why GMD:** Model intercomparison focus; code/data archival requirements ensure reproducibility.

**Extension Hooks:**
- Machine-learning-assisted profile selection (XGBoost on meteorological predictors).
- Hybrid blending strategies (e.g., MOST below, local scaling aloft).
- Sensitivity to vertical discretization (sigma vs height coordinates).

---

### 3C. Dynamic Turbulent Prandtl and Richardson Asymptote
**Title (draft):** "Dynamic Turbulent Prandtl Number Closures for Stable Layers: Reconciling Near-Neutral Behavior with Critical Richardson Asymptotes"

**Target Journal:** *Journal of Fluid Mechanics* or *Physics of Fluids*

**Core Contribution:**
- Pr_t(Ri) = 1 + a₁Ri + a₂Ri² matched to observed heat/momentum flux ratio evolution.
- Asymptotic Ri_crit behavior: continuous vs sharp cutoff; observational evidence from SHEBA, CASES-99.
- Curvature impact: Pr_t slope enters V_log → affects 2Δ calibration and inflection location.
- Thermodynamic consistency: energy dissipation constraints on Pr_t(Ri) functional form.

**Why JFM/PoF:** Fundamental turbulence theory; rigorous derivation audience.

**Extension Hooks:**
- Higher-moment closures (variance, TKE) coupling.
- Intermittency and regime transitions (turbulent ↔ wave-like).
- Laboratory analogue experiments (stratified wind tunnel).

---

## Tier 4: Planetary & Climate Extensions (24–48 months)

### 4A. Mars & Titan Boundary Layers
**Title (draft):** "Planetary Adaptation of Curvature-Aware MOST: Mars Dust Devil Convection and Titan Methane Cycle Boundary Layers"

**Target Journal:** *Icarus* or *Planetary and Space Science*

**Core Contribution:**
- Scale MOST framework: replace (g, R, c_p, θ_v) with Mars CO₂ and Titan N₂–CH₄ values.
- Compute L, ζ, Ri_g curvature for observed Mars lander (InSight) and Titan Huygens descent profiles.
- Show curvature diagnostics correlate with dust devil tracks (Mars) and methane condensation onset (Titan).
- Propose curvature-based triggers for subgrid convection schemes in GCMs.

**Why Icarus:** Premier planetary atmospheres journal; interdisciplinary readership.

**Extension Hooks:**
- Venus dense-atmosphere SBL analogue.
- Gas giants: use retrieved N², S from Juno/Cassini for J=N²/S² curvature analogue.
- Exoplanet habitability: boundary layer mixing and surface–atmosphere decoupling.

---

### 4B. Arctic Amplification Feedback Attribution
**Title (draft):** "Boundary Layer Mixing Biases and Arctic Amplification: Quantifying the Role of Stable Layer Curvature in Climate Model Projections"

**Target Journal:** *Nature Climate Change* (high-impact target) or *Geophysical Research Letters*

**Core Contribution:**
- CMIP6 model survey: diagnose Ri_b vs Ri_g(z_g) bias in polar winter monthly means.
- Show correlation between excessive SBL mixing (low 2Δ awareness) and underestimated amplification factor.
- Offline curvature correction experiment: 1.5–2.5% amplification factor recovery in idealized runs.
- Policy relevance: improved SBL → better sea-ice melt timing → adaptation timeline shifts.

**Why NCC/GRL:** Climate policy impact; short format emphasizes key result.

**Extension Hooks:**
- Permafrost thaw rate sensitivity to SBL depth.
- AO/NAO teleconnection strength and curvature-corrected polar heating.
- Cloud–radiation interaction modulation via improved inversion structure.

---

### 4C. Rotating Cores & MHD Geodynamo (Exploratory)
**Title (draft):** "Curvature-Classifier Diagnostics for Stratified Rotating Flows: Application to Earth's Outer Core Quasi-Geostrophic Turbulence"

**Target Journal:** *Journal of Fluid Mechanics* (theory) or *Geophysical Journal International* (Earth application)

**Core Contribution:**
- Adapt curvature framework to (Ro, E, Λ) parameter space (Rossby, Ekman, Elsasser numbers).
- Define J = N²/S² analogue for core stratification vs shear; compute "curvature" of stability metric.
- Show regime classification (Coriolis-dominated, magnetically-suppressed, turbulent) via curvature signatures.
- Validate against numerical dynamo simulations and geomagnetic secular variation constraints.

**Why JFM/GJI:** Novel application domain; showcases universality of curvature diagnostics.

**Extension Hooks:**
- Stellar convection zones (solar/stellar astrophysics).
- Laboratory rotating annulus experiments.
- Coupling to inverse problem (core flow inference from geomagnetic observations).

---

## Tier 5: Methodological & Tool Papers (Ongoing Support)

### 5A. Open-Source Toolkit Release
**Title (draft):** "ABL-Curvature: An Open-Source Python/Julia Toolkit for Curvature-Aware Monin–Obukhov Similarity Diagnostics"

**Target Journal:** *Geoscientific Model Development* or *Journal of Open Source Software* (JOSS)

**Core Contribution:**
- GitHub repository: https://github.com/DavidEngland/ABL
- Module suite: profiles.py (φ functions, ζ↔Ri inversion), curvature.py (diagnostics), correction.py (G factor, Q-SBL).
- Documentation: tutorial notebooks (Jupyter), test suite (pytest), continuous integration (GitHub Actions).
- Archival DOI via Zenodo; versioned releases for reproducibility.

**Why GMD/JOSS:** Software citation credit; community adoption facilitated.

**Extension Hooks:**
- Fortran/C++ bindings for legacy model integration.
- Cloud-optimized NetCDF I/O for ERA5 / CMIP climatologies.
- Web interface (Streamlit/Dash) for interactive curvature exploration.

---

### 5B. Graduate Student Training Module
**Title (draft):** "Teaching Stable Boundary Layer Physics via Curvature Diagnostics: A Modular Curriculum for Graduate Atmospheric Science"

**Target Journal:** *Bulletin of the American Meteorological Society* (education section) or *Journal of Geoscience Education*

**Core Contribution:**
- Lecture outline (as developed in intro.md).
- Hands-on exercises: ζ inversion, curvature calculation, grid-bias quantification.
- Assessment rubric: problem sets + mini-project (analyze local tower data).
- Supplementary materials: slides (LaTeX Beamer), Jupyter notebooks, solution keys.

**Why BAMS:** High visibility; education community engagement.

**Extension Hooks:**
- Online course adaptation (Coursera, edX).
- Flipped classroom implementation reports.
- Interdisciplinary variants (engineering, planetary science).

---

## Strategic Priorities & Sequencing

### Phase 1 (Year 1): Establish Core Framework
1. **1A (BLM):** Curvature correction methodology → foundational citation.
2. **1B (MWR):** Ri inversion → operational utility demonstration.
3. **1C (QJRMS):** Surface roughness → closes geometric mean justification loop.

### Phase 2 (Year 2): Validate & Extend
4. **2A (JGR-A):** Arctic application → climate relevance.
5. **2B (AMT/BLM):** Urban remote sensing → instrumentation synergy.
6. **5A (GMD/JOSS):** Toolkit release → community adoption.

### Phase 3 (Year 3+): Broaden Impact
7. **3B (GMD):** Multi-profile comparison → comprehensive reference.
8. **4B (NCC/GRL):** Arctic Amplification → high-impact policy connection.
9. **4A (Icarus):** Planetary → showcase universality.

### Opportunistic (As Data/Collaborations Emerge)
- **2C (JAS/BLM):** LLJ case studies.
- **3A (BLM/JFM):** Variable L(z) theory.
- **3C (JFM):** Dynamic Pr_t.
- **4C (JFM/GJI):** MHD geodynamo (if collaboration with geophysics group).

---

## Collaboration Strategy with McNider & Biazar

### McNider Strengths & Interests
- 45-year SBL focus; limit-cycle mixing, slope flows, iterative L solvers.
- Strong MOST foundation; appreciates rigorous math paired with operational pragmatism.
- Policy/impact orientation (Arctic, air quality, wind energy).

**Ideal Roles:**
- Co-PI on NSF proposals linking curvature work to climate/Arctic themes (2A, 4B).
- Senior author on foundational BLM papers (1A, 3A) providing historical context + field campaign connections.
- Manuscript reviews emphasizing physical intuition and mesoscale model integration.

### Biazar Strengths & Interests
- Air quality modeling (urban BL, chemistry coupling).
- Computational efficiency; grid-dependence expertise.
- Remote sensing retrieval algorithms.

**Ideal Roles:**
- Lead author on urban application (2B) given UAH megacity focus.
- Algorithm implementation specialist for operational NWP (1B).
- Co-author on toolkit paper (5A) for testing/validation infrastructure.

### Joint Authorship Templates
- **Core theory (1A, 1C, 3A):** England (lead derivation), McNider (context/validation), Biazar (numerical experiments).
- **Applications (2A, 2B, 2C):** Biazar or England (lead), McNider (interpretation), shared data analysis.
- **Planetary/climate (4A, 4B):** England (lead), McNider (climate feedback framing), Biazar (sensitivity tests).

---

## Funding Opportunities (Aligned with Publication Plan)

### NSF Atmospheric & Geospace Sciences
- **Program:** Physical and Dynamic Meteorology; Polar Programs
- **Proposal Theme:** "Curvature-Aware Boundary Layer Closures for Arctic Climate Models"
- **Budget:** Graduate student (2 years), field campaign participation (MOSAIC follow-on), computing.
- **Deliverables:** Papers 1A, 2A, 3A, 4B.

### DOE Atmospheric System Research (ASR)
- **Program:** Boundary Layer & Turbulence
- **Proposal Theme:** "Grid-Dependence Diagnostics for ARM Site Observations: Curvature-Based SBL Corrections"
- **Budget:** Postdoc (18 months), ARM data access/analysis, LES validation.
- **Deliverables:** Papers 1A, 2A, 2C, 5A.

### NASA Planetary Atmospheres
- **Program:** Solar System Workings
- **Proposal Theme:** "MOST Adaptation for Mars InSight and Titan: Boundary Layer Curvature Diagnostics"
- **Budget:** Co-I support, Mars/Titan GCM runs, student travel to DPS meeting.
- **Deliverables:** Paper 4A.

### Industry / Operational Partners
- **NOAA/NCEP:** Transition curvature correction (1B) into GFS/RAP → tech report + internal seminar.
- **ECMWF:** Contribute Q-SBL option to IFS boundary layer scheme → workshop participation + GMD note.
- **Wind Energy (NREL):** Improve LLJ forecasts (2C) for wind farm siting → white paper + site study.

---

## Success Metrics & Timeline

### Year 1 Milestones
- Submit 1A to BLM (Q2).
- Draft 1B and 1C (Q3–Q4).
- NSF proposal submitted (Q4).
- Toolkit alpha release (Q4).

### Year 2 Milestones
- 1A published; 1B, 1C submitted.
- 2A drafted (Arctic application).
- Toolkit v1.0 + GMD submission (5A).
- Present at AMS Annual Meeting (poster + oral).

### Year 3 Milestones
- 2A, 2B published.
- 4B submitted (Arctic Amplification).
- Planetary paper (4A) drafted.
- Invited seminar circuit (NCAR, ECMWF, UAH, GFDL).

### Long-Term (5 years)
- 8–12 publications in Tiers 1–4.
- Toolkit adopted by ≥2 operational centers.
- Graduate course module (5B) taught at ≥3 universities.
- Follow-on proposals (DOE, NASA) funded.

---

## Risk Mitigation & Contingencies

### Risk 1: Data Access / Validation Bottlenecks
**Mitigation:** Prioritize public datasets (ARM, SHEBA archive, GABLS LES); coordinate early with site PIs.

### Risk 2: Journal Rejection / Major Revisions
**Mitigation:** Pre-submission reviews by McNider/Biazar + external collaborator (e.g., Gert-Jan Steeneveld, Wageningen); prepare fallback journals.

### Risk 3: Computational Resource Limits
**Mitigation:** Use NSF XSEDE allocation for LES; negotiate UAH ESSC cluster time; leverage cloud credits (AWS, Google Earth Engine) for climatology.

### Risk 4: Collaboration Dynamics / Authorship Disputes
**Mitigation:** Contribution matrix drafted at proposal stage; quarterly authorship review meetings; CRediT taxonomy adopted.

### Risk 5: Scope Creep (Planetary, MHD Extensions)
**Mitigation:** Clearly delineate Tier 4/5 as "opportunistic"; focus 80% effort on Tiers 1–2 until core framework published.

---

## Summary Table: Top 5 Near-Term Targets

| # | Title (Short) | Journal | Lead | Timeline | Funding Link |
|---|---------------|---------|------|----------|--------------|
| 1A | Curvature Correction (SBL) | BLM | England | 6 mo | NSF-AGS |
| 1B | Ri Inversion & Closures | MWR | England | 9 mo | DOE-ASR |
| 1C | Geometric Mean Heights | QJRMS | England | 12 mo | NSF-AGS |
| 2A | Arctic Validation | JGR-A | Biazar/England | 18 mo | NSF-Polar |
| 5A | Toolkit Release | GMD/JOSS | England | 15 mo | DOE-ASR |

---

## Next Actions

1. **Draft 1A manuscript:** Target 6000 words + 8 figures; circulate to McNider/Biazar by [date].
2. **Refine NSF proposal:** Integrate 1A, 2A, 3A into 3-year plan; budget justification for tower data + LES.
3. **Schedule authorship meeting:** Zoom call to assign roles for 1B, 1C; CRediT taxonomy discussion.
4. **Toolkit sprint:** Consolidate profiles.py, curvature.py, add docstrings, write tutorial notebook (Jupyter).
5. **Conference abstracts:** AMS 2026 (oral: 1A results), EGU 2026 (poster: 2B urban application).

---

**Document Status:** Living strategy; update quarterly as papers progress and new opportunities emerge.

**Contact:** David E. England (david@davidengland.org), Richard T. McNider, Arastoo P. Biazar

**Repository:** https://github.com/DavidEngland/ABL (public; see CONTRIBUTING.md for collaboration guidelines)

---
<!-- End publication topics document -->
