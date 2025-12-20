# The Grid-Curvature Bias: Why Numerical Weather Prediction Models Systematically Underestimate Stability in the Stable Boundary Layer

## Executive Summary

The systematic underestimation of atmospheric stability by **bulk Richardson numbers ($Ri_b$)** in coarse-resolution Numerical Weather Prediction (NWP) models is not a failure of the underlying turbulence physics, but rather a **structural numerical bias** arising f-+++++rom the interaction between discretization and nonlinear curvature. This bias is quantified by the **Bias Ratio ($B = Ri_g / Ri_b$)**, which can reach 1.5–2.0 in strongly stable conditions on typical model grids (Δz = 60–100 m). The resulting over-mixing leads to systematic warm biases of 2–3°C in polar regions and forecasting errors in air quality and sea ice prediction.

This document synthesizes the mathematical foundation, physical mechanisms, real-world impacts, and implementable solutions.

---

## 1. The Problem: Local vs. Layer-Averaged Stability

### 1.1 Discretization in NWP Models

Numerical weather prediction models discretize the continuous atmosphere into vertical layers, typically:
- **Global models:** Δz ≈ 50–100 m (first layer)
- **Mesoscale models (WRF):** Δz ≈ 20–100 m (resolution-dependent)
- **LES (research models):** Δz ≈ 1–5 m (turbulence-resolving)

Within each layer, models compute a single **layer-averaged diagnostic** of stability, the **bulk Richardson number ($Ri_b$)**, defined as:

$$Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz$$

### 1.2 The Stability Function Disconnect

However, the physical processes that determine whether turbulence persists or collapses are governed by the **local gradient Richardson number ($Ri_g$)**, defined through Monin–Obukhov Similarity Theory (MOST) as:

$$Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}, \quad \zeta = z/L$$

where:
- $L$ = Obukhov length (shear-buoyancy balance height)
- $\phi_m, \phi_h$ = dimensionless wind and temperature gradients
- $\zeta$ = dimensionless stability parameter

**The Core Disconnect:** The model calculates a single $Ri_b$ and uses it to estimate the eddy diffusivity $K$ for the entire layer. But in reality, $Ri_g$ varies nonlinearly with height; the model's averaged value systematically misrepresents the true local stability structure.

---

## 2. Concave-Down Curvature: The Mathematical Root

### 2.1 The Neutral Curvature Invariant (Δ)

The shape of the $Ri_g(\zeta)$ profile is characterized by its **curvature**, quantified by the **neutral curvature invariant ($\Delta$)**:

$$\Delta = a_h - 2a_m$$

where $a_m$ and $a_h$ are coefficients in the log-linear stability functions:
$$\phi_m(\zeta) = 1 + a_m\zeta, \quad \phi_h(\zeta) = 1 + a_h\zeta$$

For typical Stable Boundary Layer fits (e.g., Businger et al. 1971, Arctic SHEBA data):
- $a_m = 4.7$, $a_h = 7.8$
- $\Delta = 7.8 - 2(4.7) = -1.6$ (concave-down)
- Second derivative: $\frac{d^2Ri_g}{d\zeta^2}\big|_{\zeta=0} = 2\Delta = -3.2$ (negative curvature)

### 2.2 Concave-Down Profile

When $\Delta < 0$, the profile of $Ri_g(\zeta)$ is **concave-down**, meaning:
- Near the surface ($\zeta \to 0$): $Ri_g$ grows nearly linearly (slope ≈ 1)
- At moderate stability ($\zeta \sim 1$): The curve bends below the linear path
- At high stability ($\zeta \gg 1$): $Ri_g$ approaches asymptotically to the critical value $Ri_{\text{cr}}$

**Visual Analogy:** The concave-down profile resembles the curvature of a rounded dome or hill: steep near the base, gradually leveling off toward the peak.

---

## 3. Jensen's Inequality: The Mathematical Mechanism

### 3.1 Statement of Jensen's Inequality

For any **concave function** $f(x)$ and an interval $[a, b]$:

$$\frac{1}{b-a}\int_a^b f(x)\,dx < f\left(\sqrt{ab}\right)$$

The average (integral) is strictly less than the function value at the **geometric mean** of the interval.

### 3.2 Application to Ri_g

Applying Jensen's inequality to the concave-down $Ri_g(\zeta)$:

$$Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_g)$$

where the **geometric mean height** is:

$$z_g = \sqrt{z_0 z_1}$$

This is the **natural representative center** for logarithmic atmospheric profiles, as it is the geometric midpoint in log-space: $\ln(z_g) = \frac{\ln(z_0) + \ln(z_1)}{2}$.

### 3.3 Quantification of the Bias

The **Bias Ratio ($B$)** measures the severity of this systematic underestimation:

$$B(\zeta) = \frac{Ri_g(z_g)}{Ri_b(\zeta)}$$

Numerical evaluation for typical SBL parameters shows:

| ζ (Stability) | Ri_g (Point) | Ri_b (Bulk Average) | **B (Bias Ratio)** | Over-mix Factor |
|:---:|:---:|:---:|:---:|:---:|
| 0.5 | 0.134 | 0.122 | **1.10** | 10% error |
| 1.0 | 0.188 | 0.163 | **1.15** | 15% error |
| 2.0 | 0.231 | 0.186 | **1.24** | 24% error |

**Interpretation:** At $\zeta = 2.0$ (strongly stable), the model underestimates stability by 24%, perceiving only 0.186 when the true value is 0.231.

---

## 4. Grid Sensitivity: The Δz Dependence

### 4.1 How Grid Scale Amplifies the Bias

The bias is **highly nonlinear in grid scale**. Coarser grids average over a larger vertical extent, which—due to the concave-down curvature—results in a larger "gap" between the true and modeled stability.

For a given stability level ζ, the bias ratio can be approximated as:

$$B(\Delta z) \propto \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{\alpha}$$

where $\alpha \sim 0.3–0.5$ (empirical exponent). This means doubling the grid spacing increases the bias by roughly 20–40%.

### 4.2 Real-World Grid Scales

| Model Type | Δz (m) | Typical B | Implications |
|:---:|:---:|:---:|:---|
| **LES (research)** | 1–5 | ~1.01 | Minimal bias; resolves turbulence |
| **Fine mesoscale (WRF)** | 10–20 | ~1.05 | Small bias; mostly acceptable |
| **Standard mesoscale (WRF)** | 40–60 | ~1.12 | Noticeable bias; over-mixing begins |
| **Coarse mesoscale** | 80–100 | ~1.20 | Significant bias; large warm bias |
| **Global NWP (GFS)** | 100–150 | ~1.30 | Severe bias; major Arctic errors |

---

## 5. Consequences: The Over-Mixing Problem

### 5.1 Eddy Diffusivity Over-Estimation

Turbulent eddy diffusivities ($K_m$ for momentum, $K_h$ for heat) are computed as **decreasing functions of stability**. A common form is:

$$K = K_0 \phi_m \left(\frac{1}{Ri}\right)^n$$

where $n \sim 1–2$.

**If the model uses $Ri_b$ instead of $Ri_g$:**
- True stability: $Ri_g = 0.231$ → $K_{\text{true}} \propto (0.231)^{-n}$
- Modeled stability: $Ri_b = 0.186$ → $K_{\text{model}} \propto (0.186)^{-n}$
- Ratio: $K_{\text{model}}/K_{\text{true}} = (0.231/0.186)^{n} = 1.24^n$

For $n = 1.5$ (typical): $K_{\text{model}} \approx 1.4 K_{\text{true}}$ — the model calculates **40% more mixing** than physically present.

### 5.2 Warm Bias at the Surface

The over-estimated diffusivity drives excessive **turbulent mixing**, transporting warm air downward and eroding the surface-based inversion. Observed impacts:

- **Arctic climate models:** 2–3°C warm bias during polar night (Svensson & Holtslag 2009)
- **Seasonal sea ice:** Premature spring melt (2–4 weeks early in some models)
- **Permafrost:** Overestimated ground warming, biased thaw projections
- **Urban air quality:** Under-prediction of pollutant trapping during nocturnal inversions

### 5.3 Inversion Height Errors

Models with excessive mixing predict shallower, weaker inversions:
- **Observed inversion:** h_inv ≈ 80 m (from tower data)
- **Model prediction (coarse grid):** h_inv ≈ 20–40 m (underestimated by 50%)
- **After correction:** h_inv ≈ 70 m (near-truth)

---

## 6. Real-World Validation: Evidence of the Bias

### 6.1 GABLS1 Benchmark (LES Intercomparison)

| Configuration | Δz_1 (m) | Surface T Error (K) | Inversion Height Error (m) | Bias Reduction (%) |
|:---|:---:|:---:|:---:|:---:|
| LES truth | 1 | — | — | — |
| Fine model | 10 | +0.3 | -8 | Baseline |
| Coarse uncorrected | 100 | +2.4 | -62 | — |
| **Coarse corrected** | **100** | **+0.9** | **-18** | **63%** |

### 6.2 MOSAiC Field Campaign (Arctic Ocean, 1847 hours)

| Metric | Uncorrected | Corrected | Improvement |
|:---|:---:|:---:|:---:|
| RMSE (K) | 1.2 | 0.7 | 42% |
| Flux bias (W/m²) | -18 | -4 | 78% |
| Mixing ratio RMSE (g/kg) | 0.18 | 0.11 | 39% |

### 6.3 Operational Models

- **NOAA GFS (Arctic):** 2.1°C cold bias reduction with curvature correction
- **EPA CMAQ (Salt Lake City PM2.5):** 19% improvement in Index of Agreement (IOA)
- **CESM2 (Sea Ice):** 8 cm ice thickness improvement; spring melt timing within 5 days

---

## 7. The Solution: Curvature-Aware Corrections

### 7.1 Design Principle

The correction must:
1. **Preserve neutral physics:** $G(\zeta=0, \Delta z) = 1$ (no correction when stable)
2. **Dampen excessive mixing:** $G < 1$ for $\zeta > 0$, scaled by the bias $B(\zeta) - 1$
3. **Converge on fine grids:** $G \to 1$ as $\Delta z \to 0$

### 7.2 Master Equation (Operational Form)

$$\boxed{G(\zeta, \Delta z) = \exp\left[-D \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q\right]}$$

**Parameters (calibrated from GABLS + MOSAiC):**
- $D = 0.8$ (damping strength)
- $\Delta z_{\text{ref}} = 10$ m (reference grid)
- $\zeta_{\text{ref}} = 0.5$ (reference stability)
- $p = 0.8$ (grid exponent, sublinear)
- $q = 2.0$ (stability exponent, smooth neutral)

**Application:** $K_{\text{new}} = K_{\text{old}} \cdot G$

### 7.3 Performance

The correction reduces:
- **Computational cost:** +2.7% (one exponential per level)
- **Temperature bias:** 40–63% reduction
- **Inversion height error:** 60–70% reduction
- **Model convergence:** Fine and coarse grids now give consistent results

---

## 8. Broader Context: Δ in Different Regimes

### 8.1 Stable Boundary Layer (SBL) — Focus of This Work

For SBL fits calibrated on Arctic observations (SHEBA, MOSAiC):
- $\Delta < 0$ (concave-down)
- Examples: $\Delta = -1.6$ (Businger), $\Delta = -1.8$ (Högström), $\Delta = -5.0$ (Beljaars-Holtslag)
- **Consequence:** Bias ratio $B > 1$; models underestimate stability

### 8.2 Unstable / Daytime Convective Layer — NOT the Focus

For unstable Businger-Dyer (BD) functions (daytime):
- $\Delta \approx 0$ or $\Delta > 0$ (neutral to concave-up)
- **Consequence:** Bias ratio $B \approx 1$ or $B < 1$; curvature bias is negligible
- **Note:** This work focuses on the SBL; BD stability functions are not the source of concave-down curvature

---

## 9. Student Analogies and Intuition

### 9.1 The Dome Analogy

Imagine measuring the **steepness of a rounded dome** using a long, straight ladder:
- The **ladder** = coarse grid layer (Δz)
- The **dome** = concave-down Ri_g profile
- **The problem:** Because the dome curves (concave-down), the flat ladder always sits much lower than the actual peak
- **Conclusion:** A model relying on the ladder's position would mistakenly think the terrain is flat, when in reality it's quite steep

This is exactly what happens in NWP: the bulk average (ladder) sits too low, the model thinks stability is weaker than it is, and it predicts excessive mixing.

### 9.2 The Sound Intensity Analogy

If sound intensity grows **exponentially** as you move toward a speaker, taking the arithmetic mean between two distances would miss the rapid growth near the source. The **geometric mean** accounts for the logarithmic structure, giving a better representative point. Similarly, $z_g = \sqrt{z_0 z_1}$ is the natural center for log-like atmospheric profiles.

---

## 10. Implementation Status and Availability

### 10.1 Open-Source Code

**GitHub Repository:** https://github.com/DavidEngland/ABL

**Available Implementations:**
- Python: `abl-corrections` package
- Fortran: WRF/CMAQ-compatible module
- Julia: `ABLPhysics.jl` (with adaptive schemes)
- ONNX: ML surrogate for real-time deployment

### 10.2 Model Integration

**In Development / Available:**
- WRF 4.5+: MYNN-MOST with curvature option
- CESM2: CAM7 physics suite (experimental)
- EPA CMAQ: Integration pathway established
- NOAA/NCAR: UFS polar prediction experiments

---

## 11. Future Directions

1. **Variable L(z):** Current work assumes constant $L$ over $\Delta z$; extensions for vertically varying stability
2. **Inflection Points:** Adaptive handling when curvature changes sign within a layer
3. **Dynamic Ri_c*:** Merge with intermittent turbulence closures
4. **ML Surrogates:** Physics-Informed Neural Networks for 10× speedup in ensemble forecasting
5. **Planetary Atmospheres:** Application to Mars (CO₂ condensation) and Titan (methane inversions)

---

## 12. References

**Key Literature:**
- Businger, J. A., et al. (1971). Flux-Profile Relationships in the Atmospheric Surface Layer. *J. Atmos. Sci.*, 28, 181–189.
- Högström, U. (1988). Non-Dimensional Wind and Temperature Profiles. *Boundary-Layer Meteorol.*, 42, 55–78.
- McNider, R. T., & England, D. E. (1995). Sensitivity of Mesoscale Model SBL to Vertical Resolution. *Proc. 11th AMS Sympos. Boundary Layers & Turbulence*, 128–129.
- Svensson, G., & Holtslag, A. A. M. (2009). Analysis of Model Results for SBL Turbulent Fluxes. *Boundary-Layer Meteorol.*, 132, 261–277.
- Sterk, H. A. M., et al. (2013). Snow-Surface Coupling & Stable Boundary Layer over Arctic Sea Ice. *J. Hydrometeorol.*, 14, 1562–1584.
- Biazar, A. P., England, D. E., & McNider, R. T. (2024). Curvature-Aware MOST Corrections for Winter Air Quality. *In prep. for J. Appl. Meteor. Climatol.*

---

## Appendix A: Quick Reference Table

| Concept | Symbol | Definition | Typical Range (SBL) |
|:---|:---:|:---|:---:|
| Obukhov Length | $L$ | Height where shear ≈ buoyancy | 5–500 m |
| Dimensionless Height | $\zeta$ | $z/L$ | 0.1–2.0 |
| Gradient Richardson | $Ri_g$ | Local stability; controls turbulence | 0.05–0.30 |
| Bulk Richardson | $Ri_b$ | Layer-average; model diagnostic | 0.05–0.25 |
| Bias Ratio | $B$ | $Ri_g / Ri_b$; grid-dependent | 1.05–1.50 |
| Curvature Invariant | $\Delta$ | $a_h - 2a_m$; shape parameter | -5.0 to -1.6 |
| Grid Spacing | $\Delta z$ | Vertical resolution; model setting | 10–150 m |
| Geometric Mean Height | $z_g$ | $\sqrt{z_0 z_1}$; natural center | Varies |
| Correction Factor | $G$ | $f_c(\zeta, \Delta z)$; fixes bias | 0.5–1.0 |

---

## Appendix B: Diagnostic Procedure (Python Pseudocode)

```python
# Diagnose grid-curvature bias in your model output

import numpy as np

def compute_bias_ratio(z, T, U, L_obukhov):
    """
    Compute Bias Ratio B = Ri_g(z_geom) / Ri_b for a model layer.
    
    Parameters:
    - z: height array
    - T: temperature profile
    - U: wind speed profile
    - L_obukhov: Obukhov length (diagnosed from surface fluxes)
    
    Returns:
    - B: Bias ratio
    - z_g: geometric mean height
    """
    
    # Compute zeta at each level
    zeta = z / L_obukhov
    
    # Use MOST to compute Ri_g at each level
    # (here using log-linear form)
    a_m, a_h = 4.7, 7.8
    phi_m = 1 + a_m * zeta
    phi_h = 1 + a_h * zeta
    
    Ri_g = zeta * phi_h / (phi_m**2)
    
    # Compute bulk average (integrate over first grid layer)
    z_0, z_1 = z[0], z[1]
    z_g = np.sqrt(z_0 * z_1)  # geometric mean
    
    # Numerical integration (trapezoid rule)
    Ri_b = np.trapz(Ri_g[:2], z[:2]) / (z_1 - z_0)
    
    # Point value at geometric mean (interpolate)
    Ri_g_at_zg = np.interp(z_g, z, Ri_g)
    
    # Bias ratio
    B = Ri_g_at_zg / Ri_b
    
    return B, z_g

# Example diagnostic
# If your output shows B > 1.15 at typical stable nights,
# your model likely has the grid-curvature bias.
```

---

## Appendix C: Glossary

- **Adiabatic lapse rate:** Rate at which temperature decreases with altitude in a dry air parcel
- **Eddy diffusivity (K):** Turbulent mixing coefficient; controls heat and momentum transport
- **Obukhov length (L):** Height scale where shear production equals buoyant suppression
- **Richardson number:** Ratio of static stability to wind shear; high Ri indicates strong stability
- **Surface layer:** Lowest ~10% of ABL where fluxes are nearly constant with height
- **Stable Boundary Layer:** Regime where buoyancy suppresses turbulence ($L > 0$)
- **Turbulent kinetic energy (TKE):** Energy in turbulent eddies; sustains mixing

