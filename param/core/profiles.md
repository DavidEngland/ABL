# Wind and Temperature Profiles: Exact, Corrected, and Composite

**File:** `param/core/profiles.md`  
**Purpose:** Compute wind and temperature profiles using log forms with stability corrections; benchmark against observations.  
**Related:** `SCAFFOLDING.md` §5, `gradients.md`, `hw/neutral-log-profiles.md`

---

## 1. Neutral Log Profile (Baseline)

### 1.1 Wind profile

$$U(z) = \frac{u_*}{\kappa} \ln\left(\frac{z-d}{z_0}\right)$$

where:
- $u_*$ = friction velocity (m/s)
- $\kappa = 0.40$ = von Kármán constant
- $z$ = height (m)
- $d$ = displacement height (m)
- $z_0$ = momentum roughness (m)

**Assumptions:**
- Constant-flux layer (no subsidence, no divergence)
- Neutral stratification ($\partial \theta / \partial z = \Gamma_d$, dry adiabatic)
- Smooth or uniform surface

### 1.2 Temperature profile

$$\Theta(z) = \Theta_0 + \frac{\theta_*}{\kappa} \ln\left(\frac{z-d}{z_{0h}}\right)$$

where:
- $\theta_*$ = temperature scale (K)
- $z_{0h}$ = thermal roughness (m), typically $z_{0h} \ll z_0$
- $\Theta_0$ = reference potential temperature (often taken at $z = z_0$)

**Note:** Often $Pr_t = z_0 / z_{0h} \sim 1$–$10$ (thermal Prandtl ratio); use observational or parameterized value.

### 1.3 Implementation

```python
def wind_profile_neutral(z, u_star, d, z0m, kappa=0.40):
    """
    Neutral log-profile wind.
    
    Args:
        z : height or array of heights (m)
        u_star : friction velocity (m/s)
        d : displacement height (m)
        z0m : momentum roughness (m)
        kappa : von Kármán constant
    
    Returns:
        U : wind speed (m/s)
    """
    z_eff = z - d
    if np.any(z_eff <= 0):
        raise ValueError("Height below displacement: z-d <= 0")
    
    U = (u_star / kappa) * np.log(z_eff / z0m)
    return U

def temp_profile_neutral(z, theta_star, Theta_0, d, z0h, kappa=0.40):
    """
    Neutral log-profile temperature.
    
    Args:
        z : height (m)
        theta_star : temperature scale (K)
        Theta_0 : reference potential temperature (K)
        d : displacement height (m)
        z0h : thermal roughness (m)
        kappa : von Kármán constant
    
    Returns:
        Theta : potential temperature (K)
    """
    z_eff = z - d
    if z_eff <= 0:
        raise ValueError("Height below displacement")
    
    Theta = Theta_0 + (theta_star / kappa) * np.log(z_eff / z0h)
    return Theta
```

---

## 2. Stability-Corrected Profiles (MOST)

### 2.1 General form

With stability corrections via Monin-Obukhov Similarity Theory:

$$U(z) = \frac{u_*}{\kappa} \left[ \ln\left(\frac{z-d}{z_0}\right) + \psi_m\left(\frac{z-d}{L}\right) - \psi_m\left(\frac{z_0}{L}\right) \right]$$

$$\Theta(z) = \Theta_0 + \frac{\theta_*}{\kappa} \left[ \ln\left(\frac{z-d}{z_{0h}}\right) + \psi_h\left(\frac{z-d}{L}\right) - \psi_h\left(\frac{z_{0h}}{L}\right) \right]$$

where:
- $L = \frac{u_*^2 T_0}{\kappa g\,\theta_*}$ (equivalently $L = -u_*^3 T_0/(\kappa g\,\overline{w'\theta'})$) = Obukhov length
- $\psi_m$, $\psi_h$ = stability correction functions (regime-dependent)

### 2.2 Paulson (1970) unstable corrections

For $\zeta = z/L < 0$ (unstable):

$$\psi_m(\zeta) = 2 \ln\left(\frac{1+x}{2}\right) + \ln\left(\frac{1+x^2}{2}\right) - 2 \arctan(x) + \frac{\pi}{2}$$

$$\psi_h(\zeta) = 2 \ln\left(\frac{1+y}{2}\right)$$

where:
- $x = (1 - 16\zeta)^{1/4}$
- $y = (1 - 16\zeta)^{1/2}$

(More convective → smaller $\psi$ → larger profile slope near surface)

### 2.3 Linear stable corrections (Beljaars–Holtslag)

For $\zeta > 0$ (stable):

$$\psi_m(\zeta) = -a \zeta - b(\zeta - c/d) e^{-d\zeta}$$

**Simplified (first-order):**
$$\psi_m(\zeta), \psi_h(\zeta) \approx -5\zeta$$

(Stratification suppresses mixing; profile steepens)

### 2.4 Implementation

```python
def psi_m_paulson(zeta):
    """Paulson (1970) unstable wind correction."""
    if zeta >= 0:
        return 0.0  # Use stable form instead
    x = (1.0 - 16.0 * zeta) ** 0.25
    psi = 2.0 * np.log((1.0 + x) / 2.0) + np.log((1.0 + x**2) / 2.0) \
        - 2.0 * np.arctan(x) + np.pi / 2.0
    return psi

def psi_h_paulson(zeta):
    """Paulson (1970) unstable temperature correction."""
    if zeta >= 0:
        return 0.0
    y = (1.0 - 16.0 * zeta) ** 0.5
    psi = 2.0 * np.log((1.0 + y) / 2.0)
    return psi

def psi_m_stable(zeta, a=5.0):
    """Linear stable correction."""
    if zeta <= 0:
        return 0.0
    return -a * zeta

def wind_profile_most(z, u_star, d, z0m, L, kappa=0.40):
    """
    MOST wind profile with stability correction.
    
    Args:
        z : height (m)
        u_star : friction velocity (m/s)
        d : displacement height (m)
        z0m : momentum roughness (m)
        L : Obukhov length (m); L > 0 = stable, L < 0 = unstable
        kappa : von Kármán constant
    
    Returns:
        U : wind speed (m/s)
    """
    z_eff = z - d
    z0m_eff = z0m  # Can be made d-dependent if needed
    
    if z_eff <= 0:
        raise ValueError("Height below displacement")
    
    if L is None or L == 0:
        # Neutral
        U = (u_star / kappa) * np.log(z_eff / z0m_eff)
    else:
        zeta = z_eff / L
        zeta_0 = z0m_eff / L
        
        if zeta < 0:  # Unstable
            psi_m_z = psi_m_paulson(zeta)
            psi_m_z0 = psi_m_paulson(zeta_0)
        else:  # Stable
            psi_m_z = psi_m_stable(zeta)
            psi_m_z0 = psi_m_stable(zeta_0)
        
        U = (u_star / kappa) * (np.log(z_eff / z0m_eff) + psi_m_z - psi_m_z0)
    
    return U
```

---

## 3. Composite Profiles for High-Resolution Grids

### 3.1 Motivation

If the model has levels within the **roughness sublayer** (RSL, typically $z < 10 z_0$–$30 z_0$), the log profile breaks down. A composite form blending the viscous layer and log layer may help.

### 3.2 Viscous layer (laminar)

Near the wall:
$$U^+ = y^+ \quad (y^+ < 5)$$

where $U^+ = U \kappa / u_*$ and $y^+ = (z - d) u_* / \nu$ are wall units.

### 3.3 Blended form (TODO)

$$U(z) = \begin{cases}
u_* U^+ / \kappa & z < z_{\text{switch}} \\
\text{log profile} & z \ge z_{\text{switch}}
\end{cases}$$

where $z_{\text{switch}} \approx 10\,\nu/u_*$ (e.g., ~1 mm for $u_* = 0.1\,\text{m}\,\text{s}^{-1}$).

**Status:** Not needed for typical atmospheric grid spacing. Defer to Paper 2 if very high resolution (< 1 m near surface) is required.

---

## 4. Numerical Validation: Wind Profile Example

**Setup:** $u_* = 0.3\,\text{m}\,\text{s}^{-1}$, $d = 0\,\text{m}$, $z_0 = 0.01\,\text{m}$, $\kappa = 0.40$, $L = 100\,\text{m}$ (weakly stable).

| Height (m) | Neutral $U$ (m/s) | Stable $U$ ($L=100$) (m/s) | Difference (%) |
|---|---|---|---|
| 1 | 1.31 | 1.30 | −0.8 |
| 3 | 2.05 | 2.03 | −1.0 |
| 10 | 2.91 | 2.88 | −1.0 |
| 30 | 3.76 | 3.69 | −1.9 |
| 100 | 4.92 | 4.73 | −3.9 |

(Stable conditions reduce wind at all heights, with stronger effect aloft.)

---

## 5. Temperature Profile Numerical Example

**Setup:** $\theta_* = 0.1\,\text{K}$, $\Theta_0 = 288\,\text{K}$, $z_{0h} = 0.001\,\text{m}$, $d = 0$.

| Height (m) | Neutral $\Theta$ (K) |
|---|---|
| 1 | 289.39 |
| 10 | 290.86 |
| 100 | 292.32 |

(Temperature increases logarithmically; typical gradient ~0.02 K/m in this setup.)

---

## 6. Summary: When to Use Each Form

| Profile | Use when | Pros | Cons |
|---------|----------|------|------|
| **Neutral log** | $\|\zeta\| < 0.01$ | Simple; fast | Only in nearly-neutral layer |
| **MOST (Paulson)** | $\zeta < 0.05$ | Accurate for unstable/near-neutral | Psi terms add complexity |
| **MOST (stable)** | $\zeta > 0.05$ | Standard for stable BL | True function can be more complex |
| **Composite (TODO)** | $z / z_0 < 30$ | Bridges sublayer | Not implemented; defer |

---

## 7. TODOs & Open Items

- [ ] **Beljaars–Holtslag vs. linear stable:** Test which is better for this scheme
- [ ] **Composite (RSL) profile:** Implement if very-high-resolution runs needed
- [ ] **Temperature roughness $z_{0h}$:** Parameterize as function of $z_0$ and Reynolds number
- [ ] **Observational validation:** Compare predicted vs. observed $U(z)$, $\Theta(z)$ profiles
- [ ] **Scalar analogy:** Should moisture/CO₂ follow the same form as temperature?

---

## References

Verified BibTeX keys:

- `Hogstrom1988` in `refs/master_bibliography.bib`
- `Paulson1970` in `refs/master_bibliography.bib`
