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

