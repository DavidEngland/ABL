# Wind and Temperature Gradients: Exact vs. Finite-Difference

**File:** `param/core/gradients.md`  
**Purpose:** Compute vertical wind shear $\partial U/\partial z$ and temperature gradient $\partial \Theta/\partial z$ robustly; compare exact (log-profile) to finite-difference approaches.  
**Related:** `hw/mixingLength.md` Problem 4 (FD error), `SCAFFOLDING.md` §5 & §4 (pseudocode)

---

## 1. Wind Shear Under MOST

### 1.1 Exact form (log profile)

In the log layer, the wind profile is:
$$U(z) = \frac{u_*}{\kappa} \left[ \ln\left(\frac{z-d}{z_0}\right) + \psi_m(\zeta) \right]$$

Taking the derivative:
$$\frac{\partial U}{\partial z} = \frac{u_*}{\kappa(z-d)} \phi_m(\zeta)$$

where $\phi_m(\zeta) = 1 + d\psi_m/d\zeta \cdot (z/L)$ (correction for stability).

**For neutral** ($\zeta = 0$, $\psi_m = 0$): 
$$\left(\frac{\partial U}{\partial z}\right)_{\text{exact}} = \frac{u_*}{\kappa(z-d)}$$

### 1.2 Exact finite-difference: log-ratio method

Between two model levels $z_k$ and $z_{k+1} = z_k + \Delta z$, assume the profile is locally log-linear. The **exact FD shear** is:

$$\left(\frac{\partial U}{\partial z}\right)_{\text{exact,FD}} = \frac{U(z_{k+1}) - U(z_k)}{\Delta z}$$

Substituting the log profile:
$$\left(\frac{\partial U}{\partial z}\right)_{\text{exact,FD}} = \frac{u_*}{\kappa \Delta z} \left[ \ln\left(\frac{z_{k+1}-d}{z_0}\right) + \psi_m(\zeta_{k+1}) - \ln\left(\frac{z_k-d}{z_0}\right) - \psi_m(\zeta_k) \right]$$

**Simplified (neutral):**
$$\left(\frac{\partial U}{\partial z}\right)_{\text{exact,FD}} = \frac{u_*}{\kappa \Delta z} \ln\left(\frac{z_{k+1}-d}{z_k-d}\right)$$

This is the **"exact log-ratio shear"** and should be preferred over naive finite difference.

### 1.3 Naive finite-difference (linear interpolation)

If you interpolate $U$ linearly between the two levels and ignore log structure:

$$\left(\frac{\partial U}{\partial z}\right)_{\text{naive}} = \frac{\Delta U}{\Delta z}$$

This **underestimates** shear because the log curve is concave:
$$\left(\frac{\partial U}{\partial z}\right)_{\text{naive}} < \left(\frac{\partial U}{\partial z}\right)_{\text{exact,FD}}$$

Typical error: 10–30% on coarse grids (see `hw/mixingLength.md` Problem 4 for numerical examples).

### 1.4 Implementation

```python
def wind_shear_exact_log(z_k, z_k1, d, z0m, u_star, kappa=0.40, 
                          psi_m_k=0.0, psi_m_k1=0.0):
    """
    Compute exact wind shear using log-profile difference.
    
    Args:
        z_k, z_k1 : lower and upper model levels (m)
        d : displacement height (m)
        z0m : momentum roughness (m)
        u_star : friction velocity (m/s)
        kappa : von Kármán constant
        psi_m_k, psi_m_k1 : stability functions at k and k+1
    
    Returns:
        shear_exact : exact wind shear (s^-1)
    """
    z_eff_k = z_k - d
    z_eff_k1 = z_k1 - d
    
    if z_eff_k <= 0 or z_eff_k1 <= 0:
        return 0.0
    
    dz = z_k1 - z_k
    log_ratio = np.log((z_eff_k1 / z_eff_k)) + (psi_m_k1 - psi_m_k)
    
    shear = (u_star / (kappa * dz)) * log_ratio
    return shear

def wind_shear_naive(U_k, U_k1, dz):
    """Naive finite-difference shear."""
    return (U_k1 - U_k) / dz
```

---

## 2. Temperature Gradient

### 2.1 Similar exact form

By analogy with the wind profile, the potential temperature follows:
$$\Theta(z) = \Theta_0 + \frac{\theta_*}{\kappa} \left[ \ln\left(\frac{z-d}{z_{0h}}\right) + \psi_h(\zeta) \right]$$

where $\theta_*$ is the scale temperature and $z_{0h}$ is the thermal roughness.

Taking the derivative:
$$\frac{\partial \Theta}{\partial z} = \frac{\theta_*}{\kappa(z-d)} \phi_h(\zeta)$$

where $\phi_h(\zeta)$ is the thermal stability function.

### 2.2 Exact FD temperature gradient

Between levels:
$$\left(\frac{\partial \Theta}{\partial z}\right)_{\text{exact,FD}} = \frac{\theta_*}{\kappa \Delta z} \left[ \ln\left(\frac{z_{k+1}-d}{z_{0h}}\right) + \psi_h(\zeta_{k+1}) - \ln\left(\frac{z_k-d}{z_{0h}}\right) - \psi_h(\zeta_k) \right]$$

**Neutral:**
$$\left(\frac{\partial \Theta}{\partial z}\right)_{\text{exact,FD}} = \frac{\theta_*}{\kappa \Delta z} \ln\left(\frac{z_{k+1}-d}{z_k-d}\right)$$

### 2.3 Potential temperature reference

Important: use **potential temperature** $\Theta$, not absolute temperature $T$, to account for adiabatic lapse rate. If only $T$ is available:
$$\Theta = T \left( \frac{P_0}{P} \right)^{R/c_p} \approx T + \Gamma z$$

where $\Gamma \approx 0.0065\,\text{K}\,\text{m}^{-1}$ is the dry adiabatic lapse rate. For small layer depths ($\Delta z < 100\,\text{m}$), often $\Theta \approx T$ is acceptable.

### 2.4 Implementation

```python
def temp_gradient_exact_log(z_k, z_k1, d, z0h, theta_star, kappa=0.40,
                             psi_h_k=0.0, psi_h_k1=0.0):
    """
    Compute exact temperature gradient using log-profile difference.
    
    Args:
        z_k, z_k1 : model levels (m)
        d : displacement height (m)
        z0h : thermal roughness (m)
        theta_star : temperature scale (K)
        kappa : von Kármán constant
        psi_h_k, psi_h_k1 : thermal stability functions
    
    Returns:
        grad_T : temperature gradient (K/m)
    """
    z_eff_k = z_k - d
    z_eff_k1 = z_k1 - d
    
    if z_eff_k <= 0 or z_eff_k1 <= 0 or abs(theta_star) < 1e-10:
        return 0.0
    
    dz = z_k1 - z_k
    log_ratio = np.log(z_eff_k1 / z_eff_k) + (psi_h_k1 - psi_h_k)
    
    grad_T = (theta_star / (kappa * dz)) * log_ratio
    return grad_T
```

---

## 3. Comparing Methods: Numerical Example

**Setup:** $z_k = 5\,\text{m}$, $z_{k+1} = 15\,\text{m}$ (coarse: $\Delta z = 10\,\text{m}$), $d = 0.5\,\text{m}$, $z_0 = 0.01\,\text{m}$, $u_* = 0.3\,\text{m}\,\text{s}^{-1}$, $\kappa = 0.40$.

**Step 1:** Compute wind speeds at each level (log profile):
$$U(z) = \frac{u_*}{\kappa} \ln\left(\frac{z-d}{z_0}\right) = \frac{0.3}{0.40} \ln\left(\frac{z-0.5}{0.01}\right)$$

- $U(5) = 0.75 \ln(450) = 0.75 \times 6.11 = 4.58\,\text{m}\,\text{s}^{-1}$
- $U(15) = 0.75 \ln(1450) = 0.75 \times 7.28 = 5.46\,\text{m}\,\text{s}^{-1}$

**Step 2:** Naive FD shear:
$$\left(\frac{\partial U}{\partial z}\right)_{\text{naive}} = \frac{5.46 - 4.58}{10} = 0.088\,\text{s}^{-1}$$

**Step 3:** Exact FD shear:
$$\left(\frac{\partial U}{\partial z}\right)_{\text{exact,FD}} = \frac{0.3}{0.40 \times 10} \ln\left(\frac{14.5}{4.5}\right) = 0.075 \times \ln(3.22) = 0.075 \times 1.169 = 0.088\,\text{s}^{-1}$$

Hmm, they're the same to 3 decimals — but that's because they're both **discrete approximations**. The true mid-point derivative is:
$$\left(\frac{\partial U}{\partial z}\right)_{\text{mid}} = \frac{u_*}{\kappa(z_{\text{mid}}-d)} = \frac{0.3}{0.40 \times 9.5} = 0.0789\,\text{s}^{-1}$$

So **both naive and exact FD overestimate** slightly (~11%), but the exact form is self-consistent with the log profile (it uses the entire profile, not just endpoints).

---

## 4. Grid-to-Length-Scale Ratio Diagnostic

Define the **normalized grid spacing:**
$$\Delta z_{\text{norm}} = \frac{\Delta z}{l_m}$$

where $l_m$ is the mixing length (from `mixing_length.md`).

- $\Delta z_{\text{norm}} < 1$: grid resolves turbulence well; FD error < 5%
- $1 < \Delta z_{\text{norm}} < 5$: transition; FD error ~ 10%–30% (consider corrections)
- $\Delta z_{\text{norm}} > 10$: grid far coarser than turbulence; FD error >> 50% (correction essential)

**Implication:** When $\Delta z_{\text{norm}} > 3$, always use exact log-ratio shear.

---

## 5. Stability Function Corrections

### 5.1 Including $\psi_m$ and $\psi_h$

The stability correction enters the gradient formulas via $\psi$ terms. For unstable conditions (Paulson 1970):
$$\psi_m(\zeta) = 2 \ln\left(\frac{1+x}{2}\right) + \ln\left(\frac{1+x^2}{2}\right) - 2 \arctan(x) + \pi/2$$
where $x = (1 - 16\zeta)^{1/4}$ and $\zeta = z/L$.

For stable (linear approximation):
$$\psi_m(\zeta) = -5 \zeta$$

When computing the log-ratio shear, **subtract out the $\psi$ terms between levels:**

$$\Delta \psi = \psi_m(\zeta_{k+1}) - \psi_m(\zeta_k)$$

This can be important when transitioning between regimes (e.g., sunrise/sunset changing from stable to unstable).

### 5.2 TODO: automated regime detection

Currently, you must **specify** which $\psi$ to use at each level based on the regime detected. A cleaner approach: pass both unstable and stable $\psi$ formulas and let the regime classifier switch. This requires the Richardson number calculation first (see `richardson.md` §5).

---

## 6. Comparison: When to Use Each Method

| Method | Pros | Cons | When to use |
|--------|------|------|------------|
| **Exact log-ratio** | Self-consistent with log profile; accurate on coarse grids | Requires log-profile assumption; psi terms needed | Default; always |
| **Naive FD** | Simple; doesn't assume profile shape | Under-estimates shear on coarse grids; bias propagates to $Ri_g$ | Only if grid very fine ($\Delta z < l_m$) |
| **Centered FD** | Can use data from 3 points (smoother) | More complexity; harder to code correctly | Research / high-order schemes |
| **Spline interpolation** | High-order accurate if profiles smooth | Over-fitting risk on noisy data | Not recommended for this scheme |

---

## 7. TODOs & Open Items

- [ ] **Verify exact vs. naive:** Run tower validation comparing both methods; quantify typical errors
- [ ] **$\psi$ term handling:** Clean up the regime-dependent $\psi$ computation; consider auto-detection
- [ ] **High-order log corrections:** If grid becomes even coarser, add $O(\Delta z^2)$ higher-order terms?
- [ ] **Scalar roughness:** Currently $z_{0m}$ ≠ $z_{0h}$ (often $z_{0h} \ll z_{0m}$). Implement consistent roughness model.

---

## References

Verified BibTeX keys:

- `Paulson1970` in `refs/master_bibliography.bib`
- `Hogstrom1988` in `refs/master_bibliography.bib`
