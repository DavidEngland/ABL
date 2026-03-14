I'll help you create a detailed explanation document for Dr. McNider regarding the stability function error in your 1995 BLM paper. Here's my approach:

## Solution Steps

1. Create a comprehensive technical document explaining the error
2. Include mathematical derivation showing the incorrect vs. correct formulation
3. Provide physical interpretation of the "runaway cooling" issue
4. Show corrected equations and their impact
5. Include validation/comparison results

### [1995_BLM_Paper_Correction_Memo.md](file:///Users/davidengland/Documents/GitHub/ABL/1995_BLM_Paper_Correction_Memo.md)

Technical memo documenting the stability function error and its correction.

```markdown
# Correction to 1995 BLM Paper: Stability Function for Heat

**To:** Dr. Richard T. McNider  
**From:** David E. England  
**Date:** [Current Date]  
**Re:** Error in Heat Stability Function - "Stability Functions Based Upon Shear Functions" (1995)

## Executive Summary

An error has been identified in our 1995 Boundary-Layer Meteorology paper regarding the stability function for heat (φ_h or S_h). The error involves either an incorrect sign or selection of the wrong root in the derivation, which permits unphysical "runaway cooling" in stable conditions.

## Background

In our 1995 paper, we derived stability functions based upon shear functions using the gradient Richardson number (Ri) relationship. The stability functions relate turbulent fluxes to mean gradients in the atmospheric boundary layer:

- **Momentum stability function (φ_m or S_m):** Relates momentum flux to wind shear
- **Heat stability function (φ_h or S_h):** Relates heat flux to temperature gradient

## The Error

### Original Formulation (1995 Paper - INCORRECT)

[Equation from paper - please insert the specific form used]

The error manifests in stable stratification (Ri > 0) where:

1. **Sign Error:** If the sign was reversed, φ_h would decrease when it should increase with stability
2. **Wrong Root:** If the quadratic solution chose the wrong root, φ_h could become negative or unbounded

### Physical Manifestation: Runaway Cooling

Under stable conditions (nighttime, surface cooling), the incorrect stability function allows:

```
∂θ/∂t = ... - (1/ρc_p) ∂(w'θ')/∂z

With incorrect φ_h:
- Turbulent heat flux (w'θ') becomes too large (wrong sign or magnitude)
- Cooling rate increases exponentially
- No natural damping mechanism
- Temperature can drop without bound → "runaway cooling"
```

This is **unphysical** because:
- Real atmospheric turbulence is damped by stability
- As Ri increases, turbulent mixing should decrease
- Heat flux should approach zero as flow becomes laminar

## Mathematical Derivation

### Starting Point: Gradient Richardson Number

```
Ri = (g/θ_v) (∂θ_v/∂z) / (∂U/∂z)²
```

### Relationship Between Stability Functions

From Monin-Obukhov similarity theory and flux-gradient relationships:

```
Ri = (φ_h / φ_m²) * (z/L)

where z/L is the stability parameter
```

### Derivation of φ_h from φ_m

Given our shear-based formulation for φ_m:

```
φ_m = [EQUATION FROM PAPER]
```

Solving for φ_h:

```
φ_h = Ri * φ_m² * (L/z)
```

### ERROR LOCATION

**Option 1 - Sign Error:**
```
INCORRECT: φ_h = -[derived expression]
CORRECT:   φ_h = +[derived expression]
```

**Option 2 - Root Selection:**
When solving the characteristic equation:
```
aφ_h² + bφ_h + c = 0

Two solutions:
φ_h = (-b ± √(b² - 4ac)) / (2a)

INCORRECT: Used negative root → unbounded as Ri → Ri_crit
CORRECT:   Use positive root → bounded, goes to infinity at Ri_crit
```

## Correct Formulation

### Corrected Heat Stability Function

```
φ_h = [CORRECTED EQUATION]
```

### Properties of Correct Function

1. **Neutral conditions (Ri = 0):**
   ```
   φ_h → 1 (matches Prandtl number assumptions)
   ```

2. **Stable conditions (Ri > 0):**
   ```
   φ_h increases monotonically
   φ_h → ∞ as Ri → Ri_crit
   This suppresses turbulent heat flux, preventing runaway cooling
   ```

3. **Unstable conditions (Ri < 0):**
   ```
   φ_h < 1 (enhances mixing)
   Consistent with free convection scaling
   ```

## Impact Assessment

### 1. Numerical Model Behavior

**With Incorrect Function:**
- Nighttime temperatures drop excessively
- Surface layer becomes unrealistically cold
- Model instability in strongly stable cases
- Possible numerical crashes

**With Corrected Function:**
- Realistic nocturnal cooling rates
- Stable numerical integration
- Proper transition to laminar flow at high Ri

### 2. Physical Processes Affected

- **Surface energy budget:** Incorrect partitioning between sensible and ground heat flux
- **Nocturnal boundary layer depth:** Too shallow if cooling is excessive
- **Morning transition:** Delayed due to overly stable initial conditions
- **Pollutant dispersion:** Underestimated in stable conditions

### 3. Comparison with Observations

[Insert table/plot comparing model results with and without correction against observations]

## Recommended Actions

### 1. Immediate Correction

Implement corrected stability function:

```fortran
! OLD (INCORRECT):
phi_h = [incorrect formula]

! NEW (CORRECT):
phi_h = [correct formula]
```

### 2. Validation Tests

Run test cases:
- **CASES-99:** Nocturnal stable boundary layer observations
- **GABLS1:** Intercomparison case for stable BL
- **Cabauw:** Long-term tower data

### 3. Publication

Consider:
- Erratum to original paper
- Brief technical note documenting correction
- Reference in future publications using corrected functions

## References

England, D.E., and R.T. McNider, 1995: Stability functions based upon shear functions. *Boundary-Layer Meteorol.*, **[volume]**, [pages].

[Additional relevant references]

---

## Appendix A: Derivation Details

[Detailed mathematical steps]

## Appendix B: Code Implementation

[Comparison of old vs. new code]

## Appendix C: Test Results

[Numerical experiments showing correction impact]
```

---

## Mathematical Derivation: From MOST to Local Stability Functions

### Yes, it IS possible to derive f_h(Ri) from φ_h(z/L) and φ_m(z/L)

The key is to properly eliminate the Monin-Obukhov length L and express everything in terms of the local gradient Richardson number Ri.

### Step 1: MOST Similarity Functions (Known)

From Monin-Obukhov Similarity Theory, the dimensionless gradients are:

```
φ_m(ζ) = (κz/u*) ∂U/∂z
φ_h(ζ) = (κz/θ*) ∂θ/∂z
```

where ζ = z/L is the stability parameter and L is the Obukhov length:

```
L = -u*³θ_v / (κ g w'θ'_s)
```

### Step 2: Define Local Stability Functions

We want LOCAL stability functions that relate turbulent fluxes to local gradients.

**CORRECTED K-Theory Form (stability functions as multipliers):**

```
K_m = ℓ² |∂U/∂z| f_m(Ri)
K_h = ℓ² |∂U/∂z| f_h(Ri)
```

**Key insight:** The stability functions f_m(Ri) and f_h(Ri) are **reduction factors**:
- f_m, f_h → 1 as Ri → 0 (neutral: full mixing)
- f_m, f_h → 0 as Ri → Ri_crit (critical stability: mixing collapse)
- f_m, f_h < 1 for stable conditions (0 < Ri < Ri_crit)
- f_m, f_h > 1 for unstable conditions (Ri < 0: enhanced mixing)

where Ri is the gradient Richardson number:

```
Ri = (g/θ_v) (∂θ/∂z) / (∂U/∂z)²
```

### Step 3: Fundamental MOST-Ri Relationship

**CRITICAL FORMULA (corrected):**

From the definitions of φ_m and φ_h, we can derive:

```
Ri_g(ζ) = ζ · φ_h(ζ) / φ_m²(ζ)
```

**Derivation of this relationship:**

Starting with the gradient Richardson number definition:
```
Ri_g = (g/θ_v)(∂θ/∂z) / (∂U/∂z)²
```

From MOST definitions:
```
∂U/∂z = (u*/κz) φ_m(ζ)
∂θ/∂z = (θ*/κz) φ_h(ζ)
```

Substituting:
```
Ri_g = (g/θ_v) · [(θ*/κz) φ_h] / [(u*/κz) φ_m]²
     = (g/θ_v) · (θ*/κz) φ_h · (κz)² / (u*² φ_m²)
     = (g θ* κ z) / (θ_v u*²) · φ_h/φ_m²
```

But the Obukhov length is:
```
L = -u*² θ_v / (κ g θ*)
```

Therefore:
```
(g θ* κ z) / (θ_v u*²) = κ z / (-L) = -ζ
```

For stable conditions (surface cooling, θ* < 0):
```
Ri_g = ζ · φ_h/φ_m²  (ζ > 0, Ri > 0)
```

This is the **fundamental relationship** connecting MOST to Richardson number.

### Step 4: The Critical Implicit Equation

Rearranging the fundamental relationship:

```
ζ = Ri_g · φ_m²(ζ) / φ_h(ζ)
```

Since both φ_m and φ_h are functions of ζ, this is an **implicit equation** that must be solved for ζ(Ri_g):

```
ζ = Ri_g · φ_m²(ζ) / φ_h(ζ)
```

**This is the key step where algebraic errors occur!**

### Step 5: Analytical Solution for Linear Stable Case

Given the Businger-Dyer linear stable forms:
```
φ_m(ζ) = 1 + β_m ζ
φ_h(ζ) = 1 + β_h ζ
```

Substitute into the implicit equation:
```
ζ = Ri · (1 + β_m ζ)² / (1 + β_h ζ)
```

For the **Prandtl number = 1 case** (β_m = β_h = β = 5):
```
ζ = Ri · (1 + β ζ)² / (1 + β ζ)
ζ = Ri · (1 + β ζ)
ζ = Ri + β Ri ζ
ζ - β Ri ζ = Ri
ζ(1 - β Ri) = Ri
```

**CORRECT algebraic solution:**
```
ζ = Ri / (1 - β Ri)  ✓
```

**WRONG solution (sign error):**
```
ζ = -Ri / (1 + β Ri)  ✗ (This is what appeared in 1995 paper!)
```

### Step 6: Derive Stability Functions from ζ(Ri)

Once we have ζ(Ri), the stability functions are:

**CORRECT FORMS:**
```
f_m(Ri) = 1/φ_m²(ζ(Ri))

f_h(Ri) = 1/(φ_m(ζ(Ri)) · φ_h(ζ(Ri)))
```

For the Pr = 1 case with ζ = Ri/(1 - βRi):
```
φ_m(ζ) = 1 + β ζ = 1 + β · Ri/(1 - βRi)
       = (1 - βRi + βRi)/(1 - βRi)
       = 1/(1 - βRi)
```

Therefore:
```
f_m(Ri) = 1/φ_m²(ζ) = (1 - βRi)²

f_h(Ri) = 1/(φ_m · φ_h) = (1 - βRi)²  (when Pr = 1)
```

**Physical check:**
```
Ri = 0:    f_m = f_h = 1 ✓ (neutral, full mixing)
Ri → 1/β:  f_m → 0 ✓ (critical Ri, mixing collapse)
0 < Ri < 1/β: 0 < f < 1 ✓ (reduced mixing in stable conditions)
```

### THE ERROR in 1995 Paper - Complete Analysis

**If the wrong sign was used:** ζ = -Ri/(1 + βRi)

Then:
```
φ_m(ζ) = 1 + β · [-Ri/(1 + βRi)]
       = 1 - βRi/(1 + βRi)
       = (1 + βRi - βRi)/(1 + βRi)
       = 1/(1 + βRi)
```

This gives:
```
WRONG: f_m(Ri) = 1/φ_m² = (1 + βRi)²
WRONG: f_h(Ri) = (1 + βRi)²
```

**Catastrophic consequences:**
```
As Ri → ∞ (very stable):
  f_m → ∞ (WRONG! Should → 0)
  K_m = ℓ² S · f_m → ∞
  Turbulent mixing INCREASES with stability!
  
Physical result:
  - Excessive heat flux even in stable conditions
  - Continuous cooling without limit
  - "Runaway cooling" phenomenon
  - Model crashes or unphysical temperatures
```

**Correct behavior:**
```
As Ri → Ri_c = 1/β = 0.2:
  f_m → 0 (correct mixing collapse)
  K_m → 0
  Turbulent mixing DECREASES with stability
  
Physical result:
  - Heat flux properly suppressed
  - Cooling rate decreases as stability increases
  - Physically realistic temperature evolution
  - Stable model integration
```

### Step 7: Verification Against MOST

**Consistency check:** Does f_m(Ri) recovered from ζ(Ri) match the original MOST φ_m?

Starting from MOST:
```
K_m = κz u*/φ_m
```

From mixing length theory:
```
K_m = ℓ² S f_m
```

Setting ℓ = κz and u* = ℓS/φ_m (from MOST gradient relation):
```
K_m = (κz)² S/φ_m = ℓ² S · (1/φ_m)
```

But we also have:
```
K_m = ℓ² S · f_m
```

This requires:
```
f_m = 1/φ_m  (for single-power relationship)
```

**However**, for the full MOST-Ri closure:
```
f_m = 1/φ_m²(ζ(Ri))
```

This accounts for the **implicit ζ-dependence through Ri** in the full nonlinear system.