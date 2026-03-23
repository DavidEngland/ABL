# Richardson Number ($Ri_g$, $Ri_b$) & Curvature Correction

**File:** `param/core/richardson.md`  
**Purpose:** Compute gradient Richardson number, bulk Richardson number, and grid-aware curvature corrections.  
**Related:** `Candidate Ri-Based Stability Functions.md`, `Critical Richardson Number.md`, `SCAFFOLDING.md` §4 (regime classifier)

Notation follows `CANONICAL_GLOSSARY.md` (local: $Ri_g$, bulk: $Ri_b$).

---

## 1. Gradient Richardson Number ($Ri_g$)

### 1.1 Definition

The **gradient Richardson number** measures the relative importance of stratification to shear:

$$Ri_g = \frac{g}{T_0} \frac{\partial \Theta/\partial z}{(\partial U/\partial z)^2}$$

where:
- $g = 9.81\,\text{m}\,\text{s}^{-2}$ (gravity)
- $T_0 \approx 288\,\text{K}$ (reference temperature)
- $\partial \Theta/\partial z$ = vertical potential temperature gradient (K m$^{-1}$)
- $\partial U/\partial z$ = vertical wind shear (s$^{-1}$)

**Interpretation:**
- $Ri_g < 0$: buoyancy-driven (unstable); convection overpowers shear
- $Ri_g \approx 0$: neutral; shear dominates buoyancy effects
- $0 < Ri_g < Ri_c$: weak stable; turbulence persists
- $Ri_g > Ri_c$ (critical): laminar; turbulence collapses

### 1.2 Connection to Obukhov parameter $\zeta$

Under MOST with constant-flux assumptions:
$$\zeta = \frac{z_{\mathrm{eff}}}{L} = Ri_g\,\frac{\phi_m^2}{\phi_h}$$

where $\phi_m$, $\phi_h$ are stability functions (dimensionless wind and temperature shear ratios). So:

$$Ri_g = \zeta\,\frac{\phi_h}{\phi_m^2}$$

This is **not an identity**; it emerges under the assumption of constant flux. In the surface layer, $Ri_g$ and $\zeta$ are related but distinct: $Ri_g$ is **local and diagnostic**, while $\zeta$ is **integrated** over the constant-flux layer depth.

---

## 2. The Curvature Bias Correction

### 2.1 Finite-difference shear error

On a coarse grid ($\Delta z$ large compared to the turbulent length scale $l_m$), the **finite-difference shear underestimates** the true gradient:

$$\left(\frac{\partial U}{\partial z}\right)_{\text{FD}} < \left(\frac{\partial U}{\partial z}\right)_{\text{exact}}$$

**Why?** The linear FD interpolates across the log profile, which is concave. The secant (FD) lies below the tangent (exact derivative). More precisely:

$$\left(\frac{\partial U}{\partial z}\right)_{\text{FD}} = \left(\frac{\partial U}{\partial z}\right)_{\text{exact}} \left(1 - \frac{\Delta z}{2 z_{\text{eff}}} + O(\Delta z^2)\right)$$

### 2.2 Cascading Richardson number bias

Since $(\partial U/\partial z)^2$ appears in the denominator of $Ri_g$:

$$Ri_g \propto \frac{1}{(\partial U/\partial z)^2}$$

An underestimated shear squaring makes $Ri_g$ **overestimated**:

$$Ri_{g,\text{FD}} = Ri_{g,\text{exact}} \cdot \frac{1}{(1 - \Delta z/(2z_{\text{eff}}))^2} \approx Ri_{g,\text{exact}} \cdot \left(1 + \frac{\Delta z}{z_{\text{eff}}} + O(\Delta z^2)\right)$$

### 2.3 The correction

Define the **bias factor**:

$$\text{bias\_factor} = \left(1 - \frac{\Delta z}{2 z_{\text{eff}}}\right)^2$$

or equivalently:

$$\text{bias\_factor} = \left(\frac{r}{1-r}\right)^2, \quad r = \frac{\Delta z}{2(z_{\text{eff}} + \Delta z/2)}$$

Apply in regime classifier (§4 of `SCAFFOLDING.md`):

$$Ri_{g,\text{corr}} = Ri_{g,\text{FD}} \times \text{bias\_factor}$$

**Constraint:** Ensure $0 < \text{bias\_factor} \le 1.0$.

### 2.4 Implementation

```python
def curvature_correction_factor(z_k, dz, d):
    """
    Compute the finite-difference shear curvature correction.
    
    Args:
        z_k : lower model level height (m)
        dz : vertical grid spacing (m)
        d : displacement height (m)
    
    Returns:
        bias_factor : correction factor; Ri_g_corrected = Ri_g_FD * bias_factor
    """
    z_eff = z_k - d
    if z_eff <= 0:
        return 1.0
    
    ratio = dz / (2.0 * z_eff)
    if ratio >= 1.0:
        # Grid spacing too large; return conservative factor
        return 1.0
    
    bias_factor = (1.0 - ratio) ** 2
    return min(bias_factor, 1.0)

def richardson_gradient_corrected(d_U_dz_FD, d_Theta_dz, g=9.81, T0=288.0, 
                                  z_k=None, dz=None, d=None):
    """
    Compute grid-corrected gradient Richardson number.
    
    Args:
        d_U_dz_FD : finite-difference wind shear (s^-1)
        d_Theta_dz : potential temperature gradient (K m^-1)
        g : gravity (m s^-2)
        T0 : reference temperature (K)
        z_k, dz, d : optional; if provided, apply curvature correction
    
    Returns:
        Ri_g, Ri_g_corr : uncorrected and corrected Richardson numbers
    """
    if abs(d_U_dz_FD) < 1e-10:
        return 0.0, 0.0
    
    shear_sq = d_U_dz_FD ** 2
    Ri_g = (g / T0) * d_Theta_dz / shear_sq
    
    if z_k is not None and dz is not None and d is not None:
        bias = curvature_correction_factor(z_k, dz, d)
        Ri_g_corr = Ri_g * bias
    else:
        Ri_g_corr = Ri_g
    
    return Ri_g, Ri_g_corr
```

### 2.5 Numerical example

Given:
- $z_k = 5\,\text{m}$, $\Delta z = 10\,\text{m}$, $d = 0\,\text{m}$
- $\partial U/\partial z = 0.1\,\text{s}^{-1}$ (exact)
- FD underestimation: suppose $(\partial U/\partial z)_{\text{FD}} = 0.08\,\text{s}^{-1}$ (20% low)
- $\partial \Theta/\partial z = 0.01\,\text{K}\,\text{m}^{-1}$ (stable)

**Uncorrected:**
$$Ri_{g,\text{FD}} = \frac{9.81}{288} \times \frac{0.01}{0.08^2} = 0.0340 \times 1.5625 = 0.0531$$

**Bias factor:**
$$\text{bias\_factor} = \left(1 - \frac{10}{2 \times 5}\right)^2 = (1 - 1)^2 = 0$$

Oops! $\Delta z = 2 z_k$ means the grid is essentially at the surface; correction breaks down. More realistic:

- $z_k = 20\,\text{m}$, $\Delta z = 5\,\text{m}$, $d = 1\,\text{m}$
- $z_{\text{eff}} = 19\,\text{m}$
- $\text{bias\_factor} = (1 - 5/38)^2 = (0.868)^2 = 0.753$

$$Ri_{g,\text{corr}} = 0.0531 \times 0.753 = 0.0400$$

(Correction lowers the overestimated $Ri_g$ back toward the true value.)

---

## 3. Bulk Richardson Number ($Ri_b$)

### 3.1 Definition

The **bulk Richardson number** compares integrated stability over a column to integrated shear:

$$Ri_b = \frac{g (z_2 - z_1)}{\bar{T}} \frac{\Delta \Theta}{\Delta U^2}$$

where $\Delta \Theta = \Theta_2 - \Theta_1$ and $\Delta U = U_2 - U_1$ are finite differences over a large layer.

### 3.2 When to use

- **Column stability:** Does the whole PBL remain turbulent?
- **Bulk estimates:** If you don't have continuous profiles, $Ri_b$ gives a quick assessment of $Ri_c$ breach
- **Caution:** $Ri_b$ can be misleading if shear is confined to a thin layer (e.g., jet); the bulk form misses vertical structure

### 3.3 Relationship to $Ri_g$

In the limit $\Delta z \to 0$:
$$Ri_b \to Ri_g$$

So $Ri_b$ is a **layer-integrated analog** of $Ri_g$. If the layer is deep and gradients smooth, $Ri_b \approx \text{(some average of local } Ri_g \text{)}$.

### 3.4 Implementation

```python
def richardson_bulk(z_2, z_1, Theta_2, Theta_1, U_2, U_1, Theta_bar, g=9.81):
    """
    Compute bulk Richardson number over a layer.
    
    Args:
        z_2, z_1 : upper and lower layer heights (m)
        Theta_2, Theta_1 : potential temperatures (K)
        U_2, U_1 : wind speeds (m/s)
        Theta_bar : reference potential temperature (K)
        g : gravity (m s^-2)
    
    Returns:
        Ri_b : bulk Richardson number
    """
    dz = z_2 - z_1
    dTheta = Theta_2 - Theta_1
    dU = U_2 - U_1
    
    if abs(dU) < 1e-6:
        return float('inf') if dTheta > 0 else -float('inf')
    
    Ri_b = (g * dz / Theta_bar) * (dTheta / dU**2)
    return Ri_b
```

---

## 4. Critical Richardson Number ($Ri_c$)

### 4.1 Definition

The **critical Richardson number** is the threshold above which turbulence is **suppressed by stratification alone, independent of shear magnitude**. In the surface layer, $Ri_c \approx 0.2$ to $0.25$ (England & McNider 1995):

$$Ri_c = \frac{1}{\beta}$$

where $\beta = 5$ is an empirical parameter (Businger et al. 1971, Paulson 1970).

### 4.2 Physical interpretation

**Near $Ri_g = Ri_c$:**
- Thermal stratification exactly balances mechanical shearing
- Turbulent overturning becomes marginally unstable
- For $Ri_g > Ri_c$, buoyancy suppresses all eddy growth (laminar limit, $K_M \to \nu$)

### 4.3 Dynamic vs. fixed $Ri_c$

**Currently:** Use fixed $Ri_c = 0.2$ (from $\beta = 5$).

**Future (Paper 2 open question):** $Ri_c$ might depend on:
- Reynolds number $Re = u_* z / \nu$ (higher $Re$ → higher Ri_c?)
- Stratification history (hysteresis?)
- Height in PBL

**For now:** Assume $Ri_c(z) = 1/\beta = 0.2$ everywhere. Update here once Paper 2 research is done.

### 4.4 Implementation

```python
def critical_richardson_number(beta=5.0):
    """Return critical Richardson number."""
    return 1.0 / beta

Ri_c = 0.2  # Global default
```

---

## 5. Regime Detection via Richardson Number

(See also `SCAFFOLDING.md` §4 and `regimes.md`)

### 5.1 Thresholds

| Regime | $Ri_g$ range | Physical meaning |
|--------|---|---|
| **Unstable** | $Ri_g < -0.05$ | Strong convection; $Ri_g$ significantly negative |
| **Near-neutral** | $-0.05 \le Ri_g \le 0.05$ | MOST near neutral; $\|\zeta\| < 0.05$ approximately |
| **Stable (weak)** | $0.05 < Ri_g < 0.9 \, Ri_c$ | Turbulence persists; stratification growing |
| **Stable (strong)** | $0.9 \, Ri_c \le Ri_g < Ri_c$ | Approaching collapse; eddy diffusivity declining rapidly |
| **Laminar** | $Ri_g \ge Ri_c$ | Turbulence collapsed; $K_M = \nu$ (molecular) |

### 5.2 Hysteresis and bimodality

**Warning:** Near $Ri_g \approx Ri_c$, the model can exhibit **flip-flop behavior**:
- If $Ri_g$ touches $Ri_c$ from below, $K_M$ suddenly collapses
- Reduced mixing allows temperature to drop further, pushing $Ri_g$ higher
- Model oscillates between laminar and weakly-turbulent states

**Mitigation:**
- Apply **damping** within $0.9 Ri_c < Ri_g < 1.1 Ri_c$ (see `SCAFFOLDING.md` §8, near-$Ri_c$ guard)
- Use **hysteresis filter:** require $\Delta Ri_g / \Delta t$ to exceed a threshold before flipping regimes
- **Validate** with tower data to see if bimodality is real or numerical artifact

---

## 6. Summary Table

| Quantity | Symbol | Formula | Use |
|----------|--------|---------|-----|
| Gradient Richardson number | $Ri_g$ | $(g/T_0)(\partial\Theta/\partial z)/(\partial U/\partial z)^2$ | Local stability; regime ID |
| Grid-corrected Richardson number | $Ri_{g,\text{corr}}$ | $Ri_g \times (1 - \Delta z / 2z_{\text{eff}})^2$ | Remove FD shear bias |
| Bulk Richardson number | $Ri_b$ | $(g \Delta z / \Theta_b)(\Delta\Theta / \Delta U^2)$ | Layer-integrated; coarse grid |
| Critical Richardson number | $Ri_c$ | $1/\beta = 0.2$ | Turbulence collapse threshold |
| Obukhov parameter | $\zeta$ | $z_{\text{eff}} / L = Ri_g\,\phi_m^2/\phi_h$ | Surface-layer stability; MOST |

---

## 7. TODOs & Open Items

- [ ] **Dynamic $Ri_c$:** Explore dependence on $Re$ or other flow parameters; implement if justified
- [ ] **Hysteresis:** Implement state-machine filter for regime transitions near $Ri_c$
- [ ] **Bimodality validation:** Compare model flip-flop rate to observations (tower data)
- [ ] **Higher-order Richardson number:** Smear $Ri_g$ over a vertical stencil to reduce oscillations?
- [ ] **Bulk-to-local mapping:** When running on coarse grid, how to infer $Ri_g$ from $Ri_b$?

---

## References

Verified BibTeX keys:

- `BusingerEtAl1971` in `refs/master_bibliography.bib`
- `EnglandMcNider1995BLM` in `refs/master_bibliography.bib`
- `Paulson1970` in `refs/master_bibliography.bib`
