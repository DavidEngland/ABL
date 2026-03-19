# Mixing Length: Prandtl, Blackadar, Stability-Modified

**File:** `param/core/mixing_length.md`  
**Purpose:** Compute turbulent mixing length $l_m$ across the neutral-to-stable spectrum; bridge surface-layer K-theory to PBL.  
**Related:** `hw/mixingLength.md` (pedagogical problems), `SCAFFOLDING.md` §6 (theory background)

---

## 1. Prandtl Mixing Length (Neutral Surface Layer)

### 1.1 Definition

In the neutral log layer, the mixing length is proportional to distance from the surface:

$$l_m = \kappa(z - d)$$

where:
- $\kappa = 0.40$ (von Kármán constant)
- $z$ = height above ground
- $d$ = displacement height (for vegetation, typically 0.6–0.8 × vegetation height)

**Derivation:** The eddy-diffusivity form $u_* = l_m |\partial U/\partial z|$ combined with the neutral MOST shear $\partial U/\partial z = u_*/(\kappa z)$ gives:
$$l_m = \kappa(z - d)$$
(see `hw/mixingLength.md` Problem 1 for full derivation)

### 1.2 Limits

- **Near surface** ($z \approx d$): $l_m \to 0$ (no turbulence at the boundary)
- **Log-layer** ($d < z < 0.1 h$): $l_m$ grows linearly
- **Free atmosphere** ($z > 0.1 h$): $l_m$ unbounded growth is unphysical

### 1.3 Implementation

```python
def prandtl_mixing_length(z, d, kappa=0.40):
    """
    Prandtl mixing length in neutral conditions.
    
    Args:
        z : height above ground (m)
        d : displacement height (m)
        kappa : von Kármán constant (default 0.40)
    
    Returns:
        l_m : mixing length (m)
    """
    z_eff = z - d
    if z_eff <= 0:
        return 0.0
    return kappa * z_eff
```

---

## 2. Blackadar Composite Mixing Length

### 2.1 Motivation

The Prandtl form $l_m = \kappa(z - d)$ overestimates mixing length above ~0.1 BL height. The asymptotic value should be limited by the larger-scale flow (BL height, Obukhov length scale).

Blackadar (1962) proposed a **reciprocal form** that blends the near-surface and asymptotic regimes:

$$\frac{1}{l_m} = \frac{1}{\kappa(z - d)} + \frac{1}{\lambda}$$

### 2.2 Rearranged form

$$l_m = \frac{\kappa(z - d)}{1 + \kappa(z - d) / \lambda}$$

**Limiting behavior:**
- $z \to d$: $l_m \to \kappa(z - d)$ (recovers Prandtl)
- $z \gg \lambda/\kappa$: $l_m \to \lambda$ (asymptotic saturation)

### 2.3 Asymptotic mixing length $\lambda$

Common choices:

| Source | Formula | Typical value | When to use |
|--------|---------|---|---|
| **Blackadar (1962)** | Integral over BL velocity | ~50–80 m | Research; requires $h$ input |
| **Fraction of BL height** | $\lambda = c_\lambda h$ with $c_\lambda \approx 0.1$ | $0.1 h$ | Daytime ML; simple if $h$ available |
| **Obukhov length** (stable) | $\lambda = c_s \|L\|$ with $c_s \approx 0.5$ | ~10–100 m (SBL) | Stable/nocturnal; self-contained |
| **Neutral limiting** | $\lambda = \max(c_\lambda h, 50\text{ m})$ | — | Hybrid; fallback if $h$ undefined |

**For this scheme:** We recommend a hybrid approach:
- **Unstable/neutral** ($\zeta < 0.05$): $\lambda = c_\lambda h$ where $h$ is the BL height passed from the PBL scheme
- **Stable** ($\zeta > 0.05$): $\lambda = c_s |L|$ where $L$ is the Obukhov length
- **Fallback** (if $h$ or $L$ unavailable): $\lambda = 100\,\text{m}$ (fixed default)

### 2.4 Implementation

```python
def asymptotic_mixing_length(regime, h_BL, L, c_lambda=0.1, c_s=0.5, lambda_default=100.0):
    """
    Compute asymptotic mixing length lambda based on regime.
    
    Args:
        regime : str, one of 'unstable', 'neutral', 'stable'
        h_BL : BL height (m); required for unstable/neutral
        L : Obukhov length (m); required for stable (L > 0)
        c_lambda : BL fraction coefficient (default 0.1)
        c_s : stable Obukhov coefficient (default 0.5)
        lambda_default : fallback if h_BL or L undefined (default 100 m)
    
    Returns:
        lambda : asymptotic mixing length (m)
    """
    if regime in ['unstable', 'neutral']:
        if h_BL is not None and h_BL > 0:
            return c_lambda * h_BL
    elif regime == 'stable':
        if L is not None and L > 0:
            return c_s * L
    return lambda_default

def blackadar_mixing_length(z, d, lambda_asym, kappa=0.40):
    """
    Blackadar composite mixing length.
    
    Args:
        z : height (m)
        d : displacement height (m)
        lambda_asym : asymptotic mixing length (m)
        kappa : von Kármán constant (default 0.40)
    
    Returns:
        l_m : mixing length (m)
    """
    z_eff = z - d
    if z_eff <= 0 or lambda_asym <= 0:
        return 0.0
    
    denom = 1.0 + (kappa * z_eff) / lambda_asym
    return (kappa * z_eff) / denom
```

### 2.5 Numerical example

With $u_* = 0.3\,\text{m}\,\text{s}^{-1}$, $d = 1\,\text{m}$, $\lambda = 100\,\text{m}$, $\kappa = 0.40$:

| $z$ (m) | Prandtl $l_m$ (m) | Blackadar $l_m$ (m) | Ratio |
|---|---|---|---|
| 5 | 1.6 | 1.56 | 0.97 |
| 10 | 3.6 | 3.43 | 0.95 |
| 20 | 7.6 | 6.89 | 0.91 |
| 50 | 19.6 | 13.79 | 0.70 |
| 100 | 39.6 | 19.52 | 0.49 |

Transition height: $z_t \approx \lambda / \kappa = 100 / 0.40 = 250\,\text{m}$ (above this level, Blackadar asymptotes).

---

## 3. Stability Correction: Three-Term Reciprocal Form

### 3.1 General form

Under stable conditions, the mixing length is further suppressed by stratification. The **three-term reciprocal form** (Nieuwstadt 1984, England & McNider 1995) is:

$$\frac{1}{l_m} = \frac{1}{\kappa(z - d)} + \frac{1}{\lambda} + \frac{1}{c_s |L|}$$

where $L$ is the Obukhov length (only active if $L > 0$, i.e., stable).

**Rearranged:**
$$l_m = \frac{1}{\frac{1}{\kappa(z-d)} + \frac{1}{\lambda} + \frac{\mathbb{1}_{L>0}}{c_s |L|}}$$

### 3.2 Alternative: multiplicative $f_m$ correction

In the unstable/neutral branch, apply the stability function $f_m$ as a damping factor:

$$l_m^{(\text{eff})} = l_m^{(\text{Blackadar})} \sqrt{f_m(Ri_g)}$$

where $f_m$ is computed from the appropriate MOST formula (Paulson for unstable, linear for stable; see `stability_functions/` module).

**Advantage:** $f_m$ already incorporates regime-specific stability physics.  
**Disadvantage:** Requires computing $Ri_g$ before mixing length (iteration coupling).

### 3.3 Recommended approach

**For this scheme:** Use a hybrid:
- **Unstable** ($\zeta < -0.05$): Multiplicative $l_m^{(\text{eff})} = l_m \sqrt{f_m}$ where $f_m = (1 - 16\zeta)^{-1/4}$ (Paulson)
- **Neutral** ($-0.05 \le \zeta \le 0.05$): Blackadar $l_m$ without additional correction
- **Stable** ($\zeta > 0.05$): Three-term reciprocal with Obukhov term

### 3.4 Implementation

```python
def stability_corrected_mixing_length(z, d, L, lambda_asym, zeta, regime, 
                                      kappa=0.40, c_s=0.5):
    """
    Compute stability-corrected mixing length.
    
    Args:
        z : height (m)
        d : displacement height (m)
        L : Obukhov length (m)
        lambda_asym : asymptotic mixing length (m)
        zeta : normalized height z/L (can be None for neutral)
        regime : str, one of 'unstable', 'neutral', 'stable'
        kappa : von Kármán constant
        c_s : stable Obukhov coefficient
    
    Returns:
        l_m : effective mixing length (m)
    """
    z_eff = z - d
    if z_eff <= 0:
        return 0.0
    
    # Blackadar base
    inv_lm = 1.0 / (kappa * z_eff) + 1.0 / lambda_asym
    
    # Regime-specific correction
    if regime == 'unstable':
        # Multiplicative: apply phi_m^(-1/2) = (1-16*zeta)^(1/8)
        # TODO: Import f_m_paulson or compute phi_m here
        phi_m = (1.0 - 16.0 * zeta) ** (-0.25) if zeta is not None else 1.0
        inv_lm_corrected = inv_lm / (phi_m ** 0.5)  # divide by sqrt(phi_m) = multiply by 1/sqrt
    elif regime == 'stable' and L is not None and L > 0:
        # Add Obukhov term
        inv_lm_corrected = inv_lm + 1.0 / (c_s * abs(L))
    else:
        # neutral
        inv_lm_corrected = inv_lm
    
    l_m = 1.0 / inv_lm_corrected
    return max(l_m, 0.0)
```

### 3.5 Numerical examples

**Example 1: Neutral, $z = 10\,\text{m}$, $d = 1\,\text{m}$, $\lambda = 100\,\text{m}$**
$$l_m = \frac{0.40 \times 9}{1 + 0.40 \times 9 / 100} = \frac{3.6}{1.036} = 3.47\,\text{m}$$

**Example 2: Stable with $L = 50\,\text{m}$, same heights**
$$\frac{1}{l_m} = \frac{1}{3.6} + \frac{1}{100} + \frac{1}{0.5 \times 50} = 0.278 + 0.01 + 0.04 = 0.328$$
$$l_m = 3.05\,\text{m}$$
(Note: $l_m$ reduced by ~12% due to stability)

**Example 3: Very stable with $L = 10\,\text{m}$**
$$\frac{1}{l_m} = 0.278 + 0.01 + 0.2 = 0.488 \Rightarrow l_m = 2.05\,\text{m}$$
(Further reduction; $l_m$ now ~42% of neutral)

---

## 4. Eddy Diffusivity: Mixing-Length Form

### 4.1 Definition

$$K_M = l_m^2 \left| \frac{\partial U}{\partial z} \right|$$

With $l_m$ = stability-corrected mixing length (from §3) and shear = exact log-ratio (see `gradients.md`).

### 4.2 Consistency with MOST

In the neutral surface layer with $l_m = \kappa(z-d)$ and $\partial U/\partial z = u_*/[\kappa(z-d)]$:

$$K_M = [\kappa(z-d)]^2 \cdot \frac{u_*}{\kappa(z-d)} = \kappa u_* (z-d) = \kappa u_* z_{\text{eff}} \cdot f_m(0)$$

**✓ Matches MOST ($f_m = 1$ at $\zeta = 0$)**

### 4.3 Why Prandtl form matters

The MOST form $K_M = \kappa u_* z$ is valid **only in the surface layer** (constant-flux layer). Above, the flux diverges because of entrainment and PBL-scale subsidence. The mixing-length form $K_M = l_m^2 |\partial U/\partial z|$ extends smoothly to the PBL interior because $l_m$ smoothly transitions to $\lambda$, not to unbounded $\kappa z$.

### 4.4 Implementation

```python
def eddy_diffusivity_mixing_length(l_m, shear_exact, nu_mol=1.5e-5):
    """
    Compute eddy diffusivity from mixing length and shear.
    
    Args:
        l_m : mixing length (m)
        shear_exact : exact vertical wind shear (s^-1)
        nu_mol : molecular viscosity (m^2 s^-1, default 1.5e-5)
    
    Returns:
        K_M : eddy diffusivity (m^2 s^-1)
    """
    K_M_turb = l_m**2 * abs(shear_exact)
    K_M = max(K_M_turb, nu_mol)  # Apply molecular floor
    return K_M
```

---

## 5. Diagnostic Quantities

### 5.1 Grid-to-Mixing-Length Ratio

$$\frac{\Delta z}{l_m}$$

This ratio indicates **grid resolution relative to the turbulent length scale**:
- $\Delta z / l_m < 1$: grid resolves eddy structure (good)
- $\Delta z / l_m \sim 1$–10: transition; FD errors moderate
- $\Delta z / l_m > 10$: grid too coarse; layer-average bias significant

When this ratio is large, use **exact log-ratio shear** (see `gradients.md`) rather than simple FD.

### 5.2 Mixing-Length Scale Separation

$$\frac{l_m}{\Delta z}$$

For validation: if model mixes stably, can compute what mixing length **would have** been needed to give the observed mixing:

$$l_m^{(\text{diagnosed})} = \sqrt{\frac{K_M}{|\partial U/\partial z|}}$$

Compare to predicted $l_m$ from this module.

---

## 6. Summary Table: Mixing Length Hierarchy

| Regime | Expression | When to use | Notes |
|--------|-----------|---|---|
| **Neutral, near-surface** | $l_m = \kappa(z-d)$ | $z < 50\,\text{m}$, $\zeta \approx 0$ | Prandtl; simple |
| **Neutral, PBL** | $l_m = \frac{\kappa(z-d)}{1 + \kappa(z-d)/\lambda}$ | $50\,\text{m} < z < h$ | Blackadar; requires $\lambda$ or $h$ |
| **Unstable** | $l_m^{(\text{eff})} = l_m \sqrt{f_m}$ where $f_m = (1-16\zeta)^{-1/4}$ | $\zeta < -0.05$ | Paulson ψ; $l_m$ amplified by buoyancy |
| **Stable (mild)** | $l_m = \frac{1}{\frac{1}{\kappa(z-d)} + \frac{1}{\lambda}}$ | $0 < \zeta < 0.1$ | Blackadar + light stability |
| **Stable (strong)** | $l_m = \frac{1}{\frac{1}{\kappa(z-d)} + \frac{1}{\lambda} + \frac{1}{c_s\|L\|}}$ | $\zeta > 0.1$, $L < 50\,\text{m}$ | Three-term; Obukhov-limited |

---

## 7. TODOs & Open Items

- [ ] **Blackadar integral formula:** Implement full $\lambda = \int_0^{z_{\max}} \|U(z)\| dz / \|U\|_{\max}$ if requested
- [ ] **TKE coupling:** If using TKE-based closure, adapt to $K_M = c_\mu l_m \sqrt{e}$ form
- [ ] **Dynamic $Ri_c$:** If critical Richardson number is made flow-dependent (Paper 2 §11 Q1), update three-term form
- [ ] **RSL (roughness sublayer):** Add height-dependent corrections for $z < z_0 \times 10$ (currently ignored)
- [ ] **Table comparison:** Run against McNider 1DBLM mixing-length profile to validate $\lambda$ choice
- [ ] **Validation:** Compare diagnostics $l_m^{(\text{diagnosed})}$ vs. predicted against tower data (±15% target)

---

## References

Verified BibTeX keys:

- `k1962` in `refs/complete_bibliography.bib` (Blackadar 1962)
- `Paulson1970` in `refs/master_bibliography.bib`
- `EnglandMcNider1995BLM` in `refs/master_bibliography.bib`
- `Nieuwstadt1984` in `refs/pbl.bib`

Note: keep citations keyed to BibTeX entries above to avoid free-typed reference drift.
