# Displacement Height: Definition, Valid Domain, Consistent Usage

**File:** `param/core/displacement.md`  
**Purpose:** Define and justify displacement height $d$; ensure consistent usage across schemeinterface.  
**Related:** `drafts/log-expansions.md` (theoretical expansion and pitfalls), `SCAFFOLDING.md` §8 (consistency contract)

---

## 1. Definition and Physical Meaning

### 1.1 What is $d$?

The **displacement height** is the height above ground at which the **zero-plane displacement** occurs due to surface roughness (buildings, vegetation). It is the **effective surface** from which the logarithmic wind profile emerges.

**Sign convention:**
- $d = 0$ for aerodynamically smooth surfaces (e.g., water, sand, crops < 0.5 m)
- $d > 0$ for rough surfaces (forests: $d \approx 0.6$–$0.8 \times h_{\text{canopy}}$; cities: $d \approx 0.5$–$0.8 \times h_{\text{building}}$)

### 1.2 Why it matters

The log profile is written as:
$$U(z) = \frac{u_*}{\kappa} \ln\left(\frac{z - d}{z_0}\right)$$

If $d$ is **omitted**, the effective roughness $z_0$ becomes a fudge factor that artificially absorbs displacement, leading to:
- **Overestimated** surface stress at low levels
- **Correct** stress at mid-levels by accident
- **Wrong parametrization** — the $z_0$ learned from one site doesn't transfer to another

### 1.3 Canonical values

| Surface | Typical $d$ (m) | Notes |
|---------|---|---|
| Water, sand | 0.0 | Smooth; no displacement |
| Grass, short crops | 0.0–0.1 | $\approx 0.1 \times h_{\text{veg}}$ |
| Tall grass, field crops | 0.1–0.5 | $\approx 0.6 \times h_{\text{veg}}$–ish |
| Young forest | 2–5 | $\approx 0.6$–$0.8 \times h$ |
| Mature forest | 5–25 | $0.6$–$0.8 \times h$; $h \sim 30$–50 m |
| Urban / city | 5–30 | $0.5$–$0.8 \times h_{\text{building}}$ |

---

## 2. Integration with Mixing Length

### 2.1 Consistent usage in surface-layer scheme

The mixing length must use the **same** effective height as the wind profile:
$$l_m = \kappa(z - d)$$

**Consistency is critical.** If the profile uses $(z-d)$ but the mixing length ignores $d$, the two will diverge above low heights.

### 2.2 Shear computation

Vertical wind shear from the log profile:
$$\frac{\partial U}{\partial z} = \frac{u_*}{\kappa(z - d)}$$

This appears in the **finite-difference form** as the **exact log-ratio shear:**
$$\left(\frac{\partial U}{\partial z}\right)_{\text{exact}} = \frac{u_*}{\kappa \Delta z} \ln\left(\frac{z_{k+1} - d}{z_k - d}\right)$$

(See `gradients.md` for details.)

---

## 3. Displacement Height from Geometric Data

### 3.1 Vegetation-based estimate (Thom 1971)

For plant canopies:
$$d = (0.6 \text{ to } 0.8) \times h_{\text{canopy}}$$

More refined (depends on leaf-area index LAI and canopy structure):
$$d \approx h_{\text{canopy}} \left(1 - \frac{z_0}{h_{\text{canopy}}}\right)$$

where the zero-plane roughness emerges from momentum absorption within the canopy.

### 3.2 Building / urban canopy (Oke 1987)

For city / built-up areas:
$$d \approx (0.5 \text{ to } 0.8) \times \bar{h}_{\text{building}}$$

Refined models (Macdonald et al. 1998) account for building-width-to-height ratios; typically more complex.

### 3.3 Blending across surfaces (heterogeneous terrain)

**TODO:** If the surface has mixed land use (e.g., 60% forest, 40% grassland), compute $d$ as an *area-weighted* mean or *aerodynamic blending*. This is non-trivial in rough terrain.

---

## 4. Displacement Height and Obukhov Length

### 4.1 Stability correction

The Obukhov length should follow the canonical definition

$$L = -\frac{u_*^3 T_0}{\kappa g\,\overline{w'\theta'}} = \frac{u_*^2 T_0}{\kappa g\,\theta_*}$$

and the dimensionless height should use effective height when displacement is active:

$$\zeta = \frac{z_{\text{eff}}}{L} = \frac{z-d}{L}$$

**Do not** hide displacement implicitly inside $z$ in $\zeta$; using $z/L$ instead of $(z-d)/L$ causes regime misclassification near rough canopies.

### 4.2 Gradient Richardson number

Similarly, $Ri_g$ should use shear computed from $(z-d)$:

$$Ri_g = \frac{g}{T_0} \frac{1}{u_*^2} \left(\frac{\partial \theta}{\partial z}\right) \left(\frac{\partial U}{\partial z}\right)^{-2}$$

where $\partial U/\partial z$ is evaluated at the effective height.

---

## 5. Consistency Contract

(Extracted from `SCAFFOLDING.md` §8)

### 5.1 The Five-Point Contract

Ensure **all uses of **$z$ in the surface-layer scheme** use the same definition and interpretation:

1. **Wind profile:** $U(z) = (u_*/\kappa) \ln[(z-d)/z_0]$ uses $z - d$
2. **Temperature profile:** $\Theta(z) = \Theta_s + (\theta_*/\kappa) \ln[(z-d)/z_{0h}] + \psi_h(\zeta)$ uses $z - d$
3. **Mixing length:** $l_m = \kappa(z - d)$ uses $z - d$
4. **Shear:** $\partial U/\partial z = (u_*/\kappa)/(z-d) \cdot \phi_m(\zeta)$ uses $z - d$
5. **Richardson number:** $Ri_g$, $\zeta = (z-d)/L$ both use $z - d$ (effective height)

**Violation example:** If (1)–(4) use $(z-d)$ but (5) uses bare $z$, then:
- Shear is computed for effective height $(z-d)$
- But $Ri_g$ is computed using $z/L$, not $(z-d)/L$
- **Result:** Incoherent stability classification, spurious regime transitions

### 5.2 Check: Neutral consistency

Under $\zeta = 0$ (neutral), the scheme should recover:
- $K_M = \kappa u_* (z - d)$ (MOST form)
- $l_m = \kappa(z - d)$ (Prandtl)
- $K_M = l_m^2 |\partial U/\partial z| = [\kappa(z-d)]^2 \cdot [u_*/\kappa(z-d)] = \kappa u_* (z-d)$ ✓

If the formulas don't close, check: *did you use $(z-d)$ everywhere?*

---

## 6. Special Cases

### 6.1 Over water (no displacement)

Set $d = 0.0$ everywhere. The logarithmic profile becomes:
$$U(z) = \frac{u_*}{\kappa} \ln\left(\frac{z}{z_0}\right)$$

**Caution:** Check that your code has no singularities when $d = 0$ (e.g., `if z - d <= 0` guards).

### 6.2 Very rough surfaces ($d \approx z$)

If displacement height approaches the model level height, the effective height becomes very small:
$$z_{\text{eff}} = z - d \ll \text{relevant scales}$$

The log profile, mixing length, and stability functions all become ill-conditioned (mixing length → 0). In practice:
- **Model top level should be well above the geometry.** If $z_{\text{max}} < 1.5 d$, increase model domain.
- **RSL (Roughness SubLayer) corrections** become necessary (see `SCAFFOLDING.md` §10.2).

### 6.3 Time-varying displacement height

In some applications (e.g., seasonal leaf-area index variation), $d$ varies. If $d(t)$ changes:
- Recompute effective height $z_{\text{eff}} = z - d(t)$ everywhere
- Re-evaluate shear, mixing length, and Richardson number
- Check: do $Ri_g$ and regime classification change smoothly, or do discontinuities arise?

---

## 7. Implementation Notes

### 7.1 Input specification

The scheme's interface (§3 of `SCAFFOLDING.md`) includes `d` as an **input parameter**:

```python
def surface_layer_fluxes(z_k, z_k1, dz, d, z0m, z0h, U_k, U_k1, Th_k, Th_k1, 
                         Th_bar, g, kappa, nu, ...):
    """
    Args:
        d : displacement height (m)
        ...
    """
    z_eff_k = z_k - d
    z_eff_k1 = z_k1 - d
    
    # Check validity
    if z_eff_k <= 0 or z_eff_k1 <= 0:
        raise ValueError(f"Model levels too close to surface: z_eff = {z_eff_k}, {z_eff_k1}")
    
    # Use z_eff in all subsequent calculations
    ...
```

### 7.2 Validation: diagnose $d$ from observations

**Reverse problem:** Given observed $U(z_1)$, $U(z_2)$, and assuming known $u_*$:
$$\frac{d}{z_1} = 1 - \frac{\kappa(U_2 - U_1)/u_*}{\ln(z_2/z_1)}$$

If the scheme's predicted $d$ doesn't match observations to ~±20%, suspect:
- $z_0$ or $z_{0h}$ incorrect
- Stability function $\psi_m$ not appropriate for site
- Displacement height actually varies with stability

---

## 8. TODOs & Open Items

- [ ] **Blending law for heterogeneous surfaces:** How to combine $d$ from multiple land classes in a grid cell?
- [ ] **Seasonal variation:** If LAI drives $d$ in vegetation models, develop $d(t)$ profile
- [ ] **Urban morphology:** Implement Macdonald (1998) or similar for more sophisticated city geometry
- [ ] **RSL corrections:** For very large $d$, add roughness-sublayer form to profiles (see `SCAFFOLDING.md` §10.2)
- [ ] **Validation:** Compare scheme-predicted surface stress (from $u_*$) to tower-observed $\tau$ for a range of surfaces

---

## References

Verified BibTeX keys:

- `Garratt1992Book` in `refs/master_bibliography.bib`
- `Stull1988Book` in `refs/master_bibliography.bib`
- `KaimalFinnigan1994Book` in `refs/master_bibliography.bib`

TODO (missing from current refs/*.bib):

- Add BibTeX entry for Macdonald et al. (1998) roughness-length urban formula
- Add BibTeX entry for Oke (1987) Boundary Layer Climates
- Add BibTeX entry for Thom (1971) canopy momentum absorption
