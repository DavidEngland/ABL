# Email Draft: Recent Mixing Length & Grid Spacing Findings

**To:** Dick McNider, Arastoo Biazzar  
**From:** David England  
**Subject:** Recent Findings: Mixing Length Parameterization, Grid Spacing, and Displacement Height  
**Date:** March 2025

---

## Body

Dear Dick and Arastoo,

I hope this note finds you well. I wanted to share some recent developments in our surface-layer parameterization work that I think will interest you, particularly around mixing length behavior relative to vertical grid spacing, the neutral-case linearity constraint, and integration with displacement height.

### 1. Linear Mixing Length in Neutral Conditions

A key anchor point we've established is that in neutral conditions ($\zeta \approx 0$), the mixing length must remain **strictly linear in effective height**:

$$l_m = \kappa(z - d)$$

where $\kappa = 0.40$ is the von Kármán constant and $d$ is displacement height. This Prandtl form is not just historical convenience—it's a *fundamental constraint* derived from the surface-layer momentum balance under MOST.

**Why this matters for grid design:** The mixing length must scale with $\Delta z$ in the vertical, not with grid cell aspect ratio. This means:
- On **coarse grids** (e.g., $\Delta z = 100$ m at $z = 10$ m), the mixing length is correctly *small* relative to the layer depth
- On **fine grids** (e.g., $\Delta z = 1$ m), the mixing length remains proportional to the smallest resolved scale, preserving near-surface structure
- A single $z_0$ and $d$ should work across different grid spacings—the Richardson number and stability classification do the rest

### 2. Three-Term Reciprocal Mixing Length

One of the most useful findings is the **three-term reciprocal form** for stability correction:

$$\frac{1}{l_m} = \frac{1}{\kappa(z - d)} + \frac{1}{\lambda_{\text{asym}}} + \mathbb{1}_{L>0} \cdot \frac{1}{c_s |L|}$$

where:
- **First term:** Prandtl near-surface structure
- **Second term:** Asymptotic limit (Blackadar saturation; ~0.1 BL height or fixed scale)
- **Third term:** Obukhov length suppression (active only when $L > 0$, i.e., stable)

This is not new (Nieuwstadt 1984; England & McNider 1995), but the implementation is clean: the three reciprocals act *independently*, making it transparent which physics dominates at each height and stability. In the neutral limit, the third term vanishes, recovering the linear Prandtl form. In strong stability ($L \ll z$), the third term dominates, dramatically suppressing $l_m$.

**Grid implications:** The reciprocal sum naturally handles the crossover from grid-scale ($\Delta z$) dominance at coarse resolution to Obukhov-scale dominance in the SBL. No ad hoc damping needed.

### 3. Displacement Height as a Unifying Constraint

A critical integration point we've worked out is the **five-point consistency contract** for displacement height:

1. **Profile argument:** $U(z) = (u_*/\kappa) \ln\left(\frac{z-d}{z_0}\right) + \cdots$
2. **Shear from profile:** $\partial U/\partial z = u_*/[\kappa(z-d)]$
3. **Richardson number:** Uses $(z-d)$ as the local scale height
4. **Mixing-length arguments:** $l_m = \kappa(z-d)$ and ultimately $\psi_m(\zeta = (z-d)/L)$
5. **Obukhov-length consistency:** $L = u_*^2 T_0 / (\kappa g \theta_*)$ in all contexts

When all five points use the **same** $(z-d)$, the scheme is internally consistent. Mismatching them (e.g., using $z$ in the profile but $(z-d)$ in mixing length) creates artificial curvature and Richardson-number bias.

### 4. Recent Physics Corrections

While implementing these modules, we've made two important corrections to the flux framework:

**a) MOST constant-flux relation:**  
The local gradient Richardson number is:
$$Ri_g = \zeta \frac{\phi_h}{\phi_m^2}$$
(not $\zeta \phi_m \phi_h$, which is dimensionally wrong). This is now consistent across the richardson.md, profiles.md, and mixing-length implementations.

**b) Heat-flux and Obukhov-length sign convention:**  
We unified the sign convention for sensible heat flux and Obukhov length to match the canonical definitions:
- $H_s = -\rho c_p u_* \theta_*$ (explicit negative sign)
- $L = u_*^2 T_0 / (\kappa g \theta_*)$ (derived from the canonical form $L = -u_*^3 T_0 / (\kappa g \langle w'\theta' \rangle)$)
- **Sign meaning:** $\theta_* < 0$ ⟺ $H_s > 0$ (daytime upward heating); $\theta_* > 0$ ⟺ $H_s < 0$ (nighttime cooling)

This is crucial for the iterative flux solver, where L drives the stability functions.

### 5. Practical Implementation

We've documented all of this in a suite of modular reference files under `param/core/`:
- **mixing_length.md:** Prandtl → Blackadar → three-term reciprocal progression, with examples
- **displacement.md:** Definition, canonical values, and consistency contract  
- **richardson.md:** Local and bulk Richardson numbers, corrected MOST relation, curvature bias fix
- **gradients.md:** Exact log-ratio finite-difference shear (avoids bias on coarse grids)
- **profiles.md:** Full MOST profiles with corrected ψ arguments
- **fluxes.md:** Friction velocity and temperature scale inversion (iterative solver detailed)
- **regimes.md:** Five-regime classification and regime-specific closures

Each module includes theory, pseudocode, numerical examples, and footnoted references to the bibliography.

### 6. Next Steps & Questions

Some areas where I'd appreciate your input:

1. **Asymptotic length $\lambda_{\text{asym}}$:** Should we use a fraction of BL height ($\lambda = 0.1 h$) in the unstable/neutral branch, or prefer the Obukhov form in all stability classes? The paper suggests a hybrid, but I'd value your intuition.

2. **Thermal roughness ($z_{0h}$):** The temperature profile currently uses $z_{0h}$ as an independent parameter. Do you have strong constraints (from tower data or land-surface schemes) on how $z_{0h}/z_0$ should vary? This affects temperature-scale estimates and bulk heat coefficients.

3. **Rimodal behavior near $Ri_c$:** In strongly stable conditions near the critical Richardson number (~0.2), tower data sometimes shows oscillations between mixing regimes. The current scheme applies gentle damping. Have you observed this, and do you favor damping or a bistability correction?

I'll attach the core modules so you can review the full derivations and sign conventions.

---

## Attachments

- `param/core/mixing_length.md` — Prandtl form, Blackadar composite, three-term reciprocal with stability correction  
- `param/core/displacement.md` — Definition, consistency contract, geometric estimates  
- `param/core/richardson.md` — Gradient and bulk Richardson numbers; MOST relation correction; curvature bias  
- `param/core/gradients.md` — Finite-difference shear methods (exact log-ratio preferred)  
- `param/core/profiles.md` — Wind and temperature profiles under MOST  
- `param/core/fluxes.md` — Surface-flux definitions, inversion algorithm, corrected Obukhov-length update  
- `param/core/regimes.md` — Five-regime stability classification  

All files are cross-referenced and use consistent notation (von Kármán $\kappa = 0.40$; critical $Ri_c = 0.2$; canonical glossary alignment).

---

## Questions for Discussion

Please let me know if you'd like me to:
- Walk through the three-term reciprocal derivation in detail
- Compare grid-spacing scaling on CASES99, SHEBA, or ARM tower data
- Discuss the thermal roughness parameterization further
- Share results from any test cases with these modules

Looking forward to your thoughts.

Best,  
David

---

*This work is part of the grid-aware surface-layer parameterization effort for Paper 2. All modules are under version control and ready for collaborative refinement.*
