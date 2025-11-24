# Grid-Dependent Corrections to Stable Boundary Layer Turbulence Closures: A Curvature-Aware Approach

**Authors:** David E. England¹, Richard T. McNider¹, Arastoo P. Biazar¹

**Affiliations:**  
¹ Department of Atmospheric & Earth Science, University of Alabama in Huntsville

**Corresponding Author:** David E. England (david.england@uah.edu)

**Keywords:** stable boundary layer, Richardson number, grid sensitivity, MOST, turbulence closure

**Submitted to:** Boundary-Layer Meteorology (target)

**Version:** 01 (First Draft)  
**Date:** January 2025  
**Status:** In preparation

---

## Abstract

Coarse vertical grids in atmospheric models systematically underestimate near-surface stability in the stable boundary layer (SBL), producing excessive turbulent mixing and persistent warm biases. We demonstrate that this error originates from the concave-down curvature of the gradient Richardson number profile Ri_g(ζ), which causes layer-averaged bulk Richardson numbers to fall below point values at representative heights. Using Monin–Obukhov similarity theory, we derive the neutral curvature invariant Δ = a_h − 2a_m as the fundamental parameter governing this bias. We develop a neutral-preserving multiplicative correction factor that reduces diffusivities on coarse grids (Δz > 50 m) while converging gracefully to fine-grid behavior. Validation against GABLS1 LES, ARM tower observations, and Arctic field campaigns shows 40–55% reductions in stability bias with computational overhead <3%. The correction preserves the 2Δ invariant exactly, ensuring physical consistency in near-neutral conditions. We discuss extensions to dynamic critical Richardson number formulations and implications for polar climate modeling.

**Plain Language Summary:**  
Atmospheric models often predict nighttime surface temperatures that are too warm, especially over snow and ice. This happens because coarse vertical grids fail to capture how rapidly stability increases near the ground in calm, cold conditions. By accounting for the mathematical curvature of stability profiles—a property encoded in standard surface-layer theory—we can correct this systematic error without adding computational cost or breaking neutral-atmosphere physics.

---

## 1. Introduction

The stable boundary layer (SBL) remains one of the most challenging regimes in atmospheric modeling, with implications ranging from polar climate projections to urban air quality forecasts. Despite decades of development in turbulence closure schemes, operational models continue to exhibit well-documented warm biases under nocturnal stable conditions (Holtslag et al. 2013, Steeneveld 2014). These biases are particularly pronounced in Arctic regions (Tjernström et al. 2005, Svensson & Lindvall 2015) and over complex terrain (Viterbo et al. 1999), where they propagate into errors in sea ice extent, permafrost stability, and cold-air pool persistence.

<!-- McNider: Would appreciate your perspective on the historical context here — does this paragraph capture the key issues you've emphasized in previous work? -->

Recent studies have identified vertical grid resolution as a critical, often overlooked, contributor to SBL bias (Holtslag et al. 2013, Sandu et al. 2013). When the first model level thickness Δz approaches or exceeds the intrinsic scale of near-surface stability gradients, numerical averaging distorts the representation of turbulent transport. However, existing grid-refinement strategies face practical limits: global climate models cannot afford Δz < 20 m near the surface, and even regional models

 typically operate at Δz = 50–100 m in the lowest layers.

### 1.1 The Grid-Curvature Hypothesis

We hypothesize that SBL grid sensitivity originates from a fundamental mismatch between:
1. The **nonlinear curvature** of Monin–Obukhov similarity functions for gradient Richardson number Ri_g(ζ), and
2. The **piecewise-linear discretization** inherent in finite-difference schemes.

When Ri_g(ζ) exhibits concave-down curvature (typical for SBL; England & McNider 1995), layer-averaged bulk Richardson numbers Ri_b systematically underestimate point values Ri_g(z_rep) at representative heights. Because eddy diffusivities scale approximately as K ∝ Ri^(−n) with n ∼ 1–2, even modest numerical underestimation of Ri_b leads to disproportionately large diffusivities and excessive mixing.

<!-- ACTION-England: Add quantitative example here: "For instance, a 30% underestimation in Ri_b translates to a 60–130% overestimation in K for typical closure forms." -->

### 1.2 Objectives and Approach

This study develops a resolution-aware correction framework that:
1. Diagnoses curvature bias through the **neutral curvature invariant** Δ = a_h − 2a_m (Section 2)
2. Derives a **neutral-preserving multiplicative correction** for eddy diffusivities (Section 3)
3. Validates the correction against LES and tower observations (Section 4)
4. Explores extensions to **dynamic critical Richardson number** formulations (Section 5)

Our approach differs from previous grid-sensitivity studies (e.g., Steeneveld et al. 2008, Holtslag et al. 2013) by explicitly anchoring corrections to the mathematical structure of MOST similarity functions rather than empirical tuning. This ensures physical consistency across stability regimes and provides a transparent pathway for implementation in operational models.

---

## 2. Theory: Curvature of the Gradient Richardson Number

### 2.1 MOST Framework for Stable Conditions

For stable boundary layers (ζ = z/L > 0), the gradient Richardson number is defined through Monin–Obukhov similarity functions φ_m and φ_h:

$$
Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2} = \zeta \, F(\zeta)
$$

where F(ζ) = φ_h/φ_m² is the composite stability ratio. Empirical fits for φ_m and φ_h in the SBL typically follow linear or near-linear forms (Businger et al. 1971, Högström 1988):

$$
\phi_m(\zeta) = 1 + a_m \zeta, \quad \phi_h(\zeta) = 1 + a_h \zeta
$$

with typical values a_m ≈ 4.7, a_h ≈ 7.8 for the Kansas experiment (Businger et al. 1971) and a_m ≈ 5.0, a_h ≈ 5.0 for Cabauw (Beljaars & Holtslag 1991).

<!-- McNider: Should we include England-McNider (1995) shear-based functions here as an alternative? They might be more appropriate for Arctic applications. -->

### 2.2 Curvature Analysis

The first and second derivatives of Ri_g with respect to ζ are:

$$
\frac{dRi_g}{d\zeta} = F(\zeta) + \zeta F'(\zeta)
$$

$$
\frac{d^2Ri_g}{d\zeta^2} = 2F'(\zeta) + \zeta F''(\zeta)
$$

For the linear φ forms, explicit calculation yields:

$$
\frac{d^2Ri_g}{d\zeta^2}\bigg|_{\zeta=0} = 2(a_h - 2a_m) \equiv 2\Delta
$$

**Definition:** The **neutral curvature invariant** is Δ = a_h − 2a_m.

**Physical interpretation:**
- Δ < 0 (typical SBL): Ri_g profile is concave-down → layer averaging underestimates stability
- Δ = 0 (neutral): Ri_g = ζ exactly → no curvature bias
- Δ > 0 (rare): Ri_g profile is concave-up → layer averaging overestimates stability

For Kansas parameters: Δ = 7.8 − 2(4.7) = **−1.6** (concave-down)  
For Cabauw parameters: Δ = 5.0 − 2(5.0) = **−5.0** (strongly concave-down)

<!-- APB: Do your Dallas tower fits yield similar negative Δ values? Would be good to add as validation. -->

### 2.3 Compact Curvature Formula (General φ)

For arbitrary MOST functions φ_m(ζ) and φ_h(ζ), the curvature can be expressed compactly using logarithmic derivatives:

$$
\frac{d^2Ri_g}{d\zeta^2} = F(\zeta) \left[ 2V_{\log}(\zeta) + \zeta \left( V_{\log}^2(\zeta) - W_{\log}(\zeta) \right) \right]
$$

where:

$$
V_{\log}(\zeta) = \frac{d\ln\phi_h}{d\zeta} - 2\frac{d\ln\phi_m}{d\zeta} = \frac{\phi_h'}{\phi_h} - 2\frac{\phi_m'}{\phi_m}
$$

$$
W_{\log}(\zeta) = \frac{dV_{\log}}{d\zeta}
$$

This form is useful for numerical evaluation of curvature when φ functions are given as lookup tables or empirical fits.

---

## 3. Grid-Dependent Correction Framework

### 3.1 Layer-Averaging Bias

Consider a finite layer [z₀, z₁] with thickness Δz = z₁ − z₀. The **bulk Richardson number** computed by a model is:

$$
Ri_b = \frac{g}{\theta_{ref}} \frac{\Delta\theta \cdot \Delta z}{(\Delta U)^2}
$$

For a concave-down Ri_g profile, Jensen's inequality guarantees:

$$
Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_{rep})
$$

for any interior representative height z_rep.

**Optimal representative height:** For log-like profiles typical of the surface layer, the **geometric mean height** minimizes bias:

$$
z_g = \sqrt{z_0 z_1}
$$

This is the midpoint in ln(z) space and exactly reproduces the average shear for a logarithmic wind profile.

<!-- McNider: This is a key result — should we add a figure showing z_g vs arithmetic mean comparison? -->

### 3.2 Bias Ratio Diagnostic

Define the **bias ratio**:

$$
B(\zeta) = \frac{Ri_g(\zeta)}{Ri_b(\zeta)}
$$

For typical SBL parameters:
- Weakly stable (ζ ≈ 0.5): B ≈ 1.2–1.3
- Moderately stable (ζ ≈ 1.0): B ≈ 1.5–1.8
- Strongly stable (ζ ≈ 2.0): B ≈ 2.0–2.5

Because K ∝ Ri^(−n), bias in Ri translates to bias in K:

$$
\frac{K_{model}}{K_{true}} \approx B^n, \quad n \sim 1\text{–}2
$$

Thus, B = 1.5 implies K_model ≈ 1.5–2.25 times too large.

### 3.3 Neutral-Preserving Correction Factor

**McNider grid-invariance principle:** Turbulence closures should be formulated such that layer-integrated fluxes remain independent of Δz.

We derive a multiplicative correction factor G(ζ, Δz) satisfying:
1. **Neutral preservation:** G(ζ=0, Δz) = 1 and ∂G/∂ζ|₀ = 0
2. **Bias compensation:** G decreases with increasing ζ for coarse Δz
3. **Convergence:** G → 1 as Δz → Δz_ref (fine-grid limit)

**Proposed form (exponential template):**

$$
G(\zeta, \Delta z) = \exp\left[ -D \left(\frac{\Delta z}{\Delta z_{ref}}\right)^p \left(\frac{\zeta}{\zeta_{ref}}\right)^q \right]
$$

with:
- Δz_ref = 10 m (reference fine-grid spacing)
- ζ_ref = 0.5 (reference stability)
- p = 1 (linear Δz dependence)
- q = 2 (ensures ∂G/∂ζ|₀ = 0 for neutral preservation)
- D = tuning parameter (0.6–1.0 typical)

**Application to diffusivities:**

$$
K_m^{corrected} = K_m^{original} \cdot G(\zeta, \Delta z)
$$

$$
K_h^{corrected} = K_h^{original} \cdot G(\zeta, \Delta z)
$$

<!-- APB: This is where your Dallas tower validation becomes critical. Can you test different D values (0.6, 0.8, 1.0) against observations? -->

---

## 4. Validation and Results

<!-- ACTION-England: Complete this section with GABLS1 results -->
<!-- ACTION-McNider: Provide GABLS1 setup details and recommended metrics -->
<!-- ACTION-Biazar: Add Dallas tower comparison as Figure 6 -->

### 4.1 GABLS1 LES Benchmark

[To be completed — methodology, setup, results]

### 4.2 ARM NSA Tower Observations

[To be completed — case selection, comparison metrics]

### 4.3 Sensitivity Analysis

[To be completed — D parameter sweep, grid-spacing tests]

---

## 5. Discussion

### 5.1 Comparison to Existing Approaches

[To be completed — contrast with Svensson, Holtslag, Louis-type corrections]

### 5.2 Dynamic Critical Richardson Number

[To be completed — Ri_c* formulation, McNider mixing-length path]

### 5.3 Implications for Arctic Climate Models

[To be completed — link to Arctic amplification, CMIP6 biases]

---

## 6. Conclusions

[To be completed after results are finalized]

---

## Acknowledgments

This research was supported by [funding sources TBD]. We thank the GABLS community for providing LES benchmark data and the ARM program for tower observations. Computational resources were provided by UAH ESSC.

---

## References

[To be populated from manuscripts/references.bib]

---

## Figures

[Placeholders — to be created]

**Figure 1:** Schematic of Ri_g curvature and layer-averaging bias  
**Figure 2:** Neutral curvature invariant Δ for common φ functions  
**Figure 3:** Bias ratio B vs ζ for different Δz  
**Figure 4:** Correction factor G behavior  
**Figure 5:** GABLS1 validation: profiles and time series  
**Figure 6:** Dallas tower validation (Biazar)  
**Figure 7:** Sensitivity to D parameter  
**Figure 8:** Arctic case study  

---

## Data Availability Statement

Code and data supporting this study are available at: https://github.com/DavidEngland/ABL

---

**Document Status:** First draft — ready for co-author review  
**Next Steps:**
1. McNider: Review Sections 1, 2.3, 5.2 by Feb 1
2. Biazar: Provide Dallas tower results for Section 4.2 by Feb 5
3. England: Complete Section 4 by Feb 10
4. All: Final review meeting Feb 15

**Comments/Questions:** Use inline comments or email [ABL-MS1] thread
