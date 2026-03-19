# Surface Fluxes: $u_*$, $\theta_*$, $\tau$, $H_s$ — Derivations & Iterative Solution

**File:** `param/core/fluxes.md`  
**Purpose:** Compute surface friction velocity, temperature scale, stresses, and sensible heat flux; iteratively invert for consistency.  
**Related:** `SCAFFOLDING.md` §7–§8 (iteration & displacement contract), `hw/neutral-log-profiles.md` Problem 5

---

## 1. Surface Friction Velocity ($u_*$)

### 1.1 Definition

The **friction velocity** is related to surface wind stress $\tau_s$ by:

$$u_* = \sqrt{\frac{\tau_s}{\rho}}$$

where, for 1D shear,

$$\tau_s = -\rho\,\overline{u'w'} = \rho u_*^2$$

and $\rho$ is air density.

In the constant-flux layer, $\tau$ is approximately uniform with height:
$$\tau(z) \approx \tau_s \quad (0 < z < z_{\text{inversion}})$$

### 1.2 Inverse of the wind profile

Given observed wind at two heights $z_k$ and $z_{k+1}$, we can **invert** the log profile to find $u_*$.

Starting from:
$$U(z_k) = \frac{u_*}{\kappa} \left[ \ln\left(\frac{z_k - d}{z_0}\right) + \psi_m(\zeta_k) \right]$$

$$U(z_{k+1}) = \frac{u_*}{\kappa} \left[ \ln\left(\frac{z_{k+1} - d}{z_0}\right) + \psi_m(\zeta_{k+1}) \right]$$

Subtracting (to eliminate $z_0$ and $\psi_m(z_0)$ terms):

$$U(z_{k+1}) - U(z_k) = \frac{u_*}{\kappa} \left[ \ln\left(\frac{z_{k+1} - d}{z_k - d}\right) + \psi_m(\zeta_{k+1}) - \psi_m(\zeta_k) \right]$$

Rearranging:
$$u_* = \frac{\kappa [U(z_{k+1}) - U(z_k)]}{\ln\left(\frac{z_{k+1} - d}{z_k - d}\right) + \psi_m(\zeta_{k+1}) - \psi_m(\zeta_k)}$$

### 1.3 Implementation

```python
def friction_velocity_from_winds(U_k, U_k1, z_k, z_k1, d, z0m, L=None, kappa=0.40):
    """
    Invert wind profile to compute friction velocity.
    
    Args:
        U_k, U_k1 : wind speeds at two levels (m/s)
        z_k, z_k1 : heights (m)
        d : displacement height (m)
        z0m : momentum roughness (m)
        L : Obukhov length (m); if None, assume neutral
        kappa : von Kármán constant
    
    Returns:
        u_star : friction velocity (m/s)
    """
    z_eff_k = z_k - d
    z_eff_k1 = z_k1 - d
    
    if z_eff_k <= 0 or z_eff_k1 <= 0:
        raise ValueError("Heights too close to displacement")
    
    dU = U_k1 - U_k
    dz = z_k1 - z_k
    
    # Log-ratio term
    log_term = np.log(z_eff_k1 / z_eff_k)
    
    # Stability correction (if L provided)
    if L is not None and L != 0:
        zeta_k = z_eff_k / L
        zeta_k1 = z_eff_k1 / L
        
        if zeta_k < 0:  # Unstable branch
            psi_m_k = psi_m_paulson(zeta_k)
            psi_m_k1 = psi_m_paulson(zeta_k1)
        else:  # Stable branch
            psi_m_k = psi_m_stable(zeta_k)
            psi_m_k1 = psi_m_stable(zeta_k1)
        
        denom = log_term + (psi_m_k1 - psi_m_k)
    else:
        denom = log_term
    
    u_star = (kappa * dU) / denom
    
    if u_star < 0:
        raise ValueError(f"Negative friction velocity: {u_star}")
    
    return u_star
```

### 1.4 Typical values

| Surface / Condition | $u_*$ (m/s) | Notes |
|---|---|---|
| Smooth water, calm | 0.05–0.10 | Very weak wind stress |
| Light wind over grass | 0.10–0.20 | 10 m/s mean wind |
| Moderate wind, rougher surface | 0.20–0.40 | 10 m/s wind, rough surface |
| Strong wind, complex terrain | 0.40–1.0+ | High-wind events or urban areas |

---

## 2. Temperature Scale ($\theta_*$)

### 2.1 Definition

The **temperature scale** represents the cross-correlation between temperature and vertical velocity fluctuations:

$$\theta_* = -\frac{\overline{w'\theta'}}{u_*}$$

$$\theta_* = -\frac{H_s}{\rho c_p u_*}$$

where $H_s = \rho c_p\,\overline{w'\theta'}$ is the sensible heat flux.

Inversely, by analogy with the wind:
$$\Theta(z) = \Theta_0 + \frac{\theta_*}{\kappa} \left[ \ln\left(\frac{z-d}{z_{0h}}\right) + \psi_h\left(\frac{z-d}{L}\right) \right]$$

### 2.2 Inverse of the temperature profile

Given $\Theta_k$, $\Theta_{k+1}$, and assuming constant heat flux:

$$\theta_* = \frac{\kappa [\Theta(z_{k+1}) - \Theta(z_k)]}{\ln\left(\frac{z_{k+1}-d}{z_k-d}\right) + \psi_h(\zeta_{k+1}) - \psi_h(\zeta_k)}$$

### 2.3 Implementation

```python
def temperature_scale_from_profile(Theta_k, Theta_k1, z_k, z_k1, d, z0h, L=None, kappa=0.40):
    """
    Invert temperature profile to compute temperature scale.
    
    Args:
        Theta_k, Theta_k1 : potential temperatures (K)
        z_k, z_k1 : heights (m)
        d : displacement height (m)
        z0h : thermal roughness (m)
        L : Obukhov length (m)
        kappa : von Kármán constant
    
    Returns:
        theta_star : temperature scale (K)
    """
    z_eff_k = z_k - d
    z_eff_k1 = z_k1 - d
    
    dTheta = Theta_k1 - Theta_k
    dz = z_k1 - z_k
    
    log_term = np.log(z_eff_k1 / z_eff_k)
    
    if L is not None and L != 0:
        zeta_k = z_eff_k / L
        zeta_k1 = z_eff_k1 / L
        
        if zeta_k < 0:
            psi_h_k = psi_h_paulson(zeta_k)
            psi_h_k1 = psi_h_paulson(zeta_k1)
        else:
            psi_h_k = psi_h_stable(zeta_k)  # TODO: implement stable form
            psi_h_k1 = psi_h_stable(zeta_k1)
        
        denom = log_term + (psi_h_k1 - psi_h_k)
    else:
        denom = log_term
    
    theta_star = (kappa * dTheta) / denom
    
    return theta_star
```

### 2.4 Sign convention

- $\theta_* < 0$: upward sensible heat flux ($H_s > 0$), typical daytime heating
- $\theta_* > 0$: downward sensible heat flux ($H_s < 0$), typical nighttime cooling
- $|\theta_*| \sim 0.01$–$0.1\,\text{K}$ (typical magnitudes)

**Note:** $\theta_*$ can be singular ($\to \infty$) when $dTheta \to 0$ (neutral). Use guard: if $|\theta_*| > 1\,\text{K}$, flag as suspicious.

---

## 3. Surface Stresses and Fluxes

### 3.1 Momentum flux (stress)

$$\tau = \rho u_*^2$$

**Wind stress components:**
$$\tau_x = -\rho u_* \overline{u'w'}  \approx -\rho u_*^2 \frac{u_k}{\|U_k\|}$$
$$\tau_y = -\rho u_* \overline{v'w'}  \approx -\rho u_*^2 \frac{v_k}{\|U_k\|}$$

(Projects the friction velocity onto the wind direction.)

### 3.2 Sensible heat flux

$$H_s = -\rho c_p u_* \theta_*$$

where $c_p = 1005\,\text{J}\,\text{kg}^{-1}\,\text{K}^{-1}$ (specific heat at constant pressure).

**Sign convention:**
- $H_s > 0$: upward heat flux from surface to air (daytime heating)
- $H_s < 0$: downward heat flux from air to surface (nighttime cooling)

### 3.3 Bulk coefficients

Often written in bulk form:
$$H_s = \rho c_p C_H |\vec{U}| (\Theta_s - \Theta_k)$$
$$\tau = \rho C_D |\vec{U}|^2$$

where $C_D, C_H$ are the **drag and heat transfer coefficients**. Inverting:
$$C_D = \frac{u_*^2}{|\vec{U}|^2}$$

$$C_H = -\frac{u_* \theta_*}{|\vec{U}| (\Theta_s - \Theta_k)}$$

(Typically $C_D \sim 1$–$3 \times 10^{-3}$ and $C_H \sim 1$–$2 \times 10^{-3}$, depending on stability and roughness.)

### 3.4 Implementation

```python
def compute_surface_fluxes(u_star, theta_star, U_wind, rho=1.2, cp=1005.0):
    """
    Compute surface momentum and heat fluxes.
    
    Args:
        u_star : friction velocity (m/s)
        theta_star : temperature scale (K)
        U_wind : wind speed (m/s)
        rho : air density (kg/m³)
        cp : specific heat (J/kg/K)
    
    Returns:
        tau : wind stress (N/m²)
        H_s : sensible heat flux (W/m²)
    """
    tau = rho * u_star**2
    H_s = -rho * cp * u_star * theta_star
    
    return tau, H_s
```

---

## 4. Iterative Inversion: Coupling $u_*$, $\theta_*$, and $L$

### 4.1 The coupling problem

$u_*$ depends on the stability function $\psi_m$, which depends on $L$, which depends on $\theta_*$, which depends on... $u_*$ again. **Circular dependency.**

Resolution: **iterative inversion** (see `SCAFFOLDING.md` §8 for pseudocode).

### 4.2 Algorithm outline

1. **Initial guess:** $u_{*,0} = \sqrt{\kappa (U_k - U_{k+1}) / (\ln(z_{k+1}/z_k) + \delta)}$ (neutral or weakly corrected)
2. **Iterate** $n = 1, 2, \ldots, N_{\max}$ (typically $N_{\max} = 50$):
   - Compute $L_n$ from current $u_{*,n-1}$ and $\theta_{*,n-1}$
   - Update $u_{*,n}$ using inversion formula with $L_n$
   - Update $\theta_{*,n}$ using inversion formula with $L_n$
   - Check convergence: if $| u_{*,n} - u_{*,n-1} | < \epsilon_u$ and $| \theta_{*,n} - \theta_{*,n-1} | < \epsilon_\theta$, **exit**
3. **Near-$Ri_c$ guard:** If $Ri_g > 0.9 Ri_c$ during iteration, apply damping (see `regimes.md` §4)

### 4.3 Implementation

```python
def invert_fluxes_iterative(U_k, U_k1, Theta_k, Theta_k1, z_k, z_k1, d, 
                             z0m, z0h, g, T0, rho, cp, kappa=0.40,
                             max_iter=50, tol_u=1e-4, tol_theta=1e-6, Ri_c=0.2):
    """
    Iteratively invert for u*, theta*, L.
    
    Args:
        U_k, U_k1 : wind speeds (m/s)
        Theta_k, Theta_k1 : potential temperatures (K)
        z_k, z_k1 : heights (m)
        d : displacement height (m)
        z0m, z0h : roughness lengths (m)
        g : gravity (m/s²)
        T0 : reference temperature (K)
        rho : air density (kg/m³)
        cp : specific heat (J/kg/K)
        kappa : von Kármán constant
        max_iter : maximum iterations
        tol_u, tol_theta : convergence tolerances
        Ri_c : critical Richardson number
    
    Returns:
        u_star, theta_star, L : final estimates
        converged : bool
    """
    z_eff_k = z_k - d
    
    # Initial guess (neutral)
    u_star = friction_velocity_from_winds(U_k, U_k1, z_k, z_k1, d, z0m, L=None)
    theta_star = temperature_scale_from_profile(Theta_k, Theta_k1, z_k, z_k1, d, z0h, L=None)
    
    for it in range(max_iter):
        # Compute Obukhov length using canonical definition:
        # L = -u_*^3 * T0 / (kappa * g * <w'theta'>)
        # with <w'theta'> = -u_* * theta_star -> L = u_*^2 * T0 / (kappa * g * theta_star)
        if abs(theta_star) < 1e-10:
            L = float('inf')  # Neutral approximation
        else:
            L = (u_star**2 * T0) / (kappa * g * theta_star)
        
        # Compute Ri_g to check for bimodality risk
        shear = (u_star / (kappa * z_eff_k))  # Approximate; could use full form
        dTheta_dz = theta_star / (kappa * z_eff_k)
        Ri_g = (g / T0) * dTheta_dz / (shear**2 + 1e-10)
        
        # Near-Ri_c damping (TODO: implement full guard)
        if Ri_g > 0.9 * Ri_c and Ri_g < Ri_c:
            damp_factor = 1.0 - (Ri_g - 0.9*Ri_c) / (0.1*Ri_c)
        else:
            damp_factor = 1.0
        
        # Invert for new u*, theta*
        u_star_new = friction_velocity_from_winds(U_k, U_k1, z_k, z_k1, d, z0m, L=L)
        theta_star_new = temperature_scale_from_profile(Theta_k, Theta_k1, z_k, z_k1, d, z0h, L=L)
        
        # Damp update
        u_star_new = damp_factor * u_star_new + (1 - damp_factor) * u_star
        theta_star_new = damp_factor * theta_star_new + (1 - damp_factor) * theta_star
        
        # Check convergence
        du = abs(u_star_new - u_star)
        dtheta = abs(theta_star_new - theta_star)
        
        u_star = u_star_new
        theta_star = theta_star_new
        
        if du < tol_u and dtheta < tol_theta:
            return u_star, theta_star, L, True
    
    # Max iterations reached
    return u_star, theta_star, L, False
```

### 4.4 Convergence issues

- **Divergence:** If gradients are inverted (e.g., model error), iteration may diverge (e.g., $u_* \to 0$ or $\infty$)
  - **Fix:** Add guards on physically reasonable ranges, e.g., $0.01 < u_* < 2.0\,\text{m}\,\text{s}^{-1}$
- **Slow convergence near neutral:** If $|\theta_*| \ll 0.01\,\text{K}$, the iteration can be sluggish
  - **Fix:** Use a more aggressive step or line-search

---

## 5. Free Convection Limit

### 5.1 When wind is very weak

If $U_k, U_{k+1} \to 0$ (near-calm conditions), the shear-based inversions can become singular or ill-conditioned.

**Free-convection regime:** Use $u_*$ from convective scaling instead:
$$u_* = (w_* h)^{1/3}$$

where $w_* = [gH_s h / (\rho c_p T_0)]^{1/3}$ is the convective velocity scale and $h$ is PBL height.

**TODO:** For now, add a guard: if $|U_k - U_{k+1}| < 0.1\,\text{m}\,\text{s}^{-1}$, warn or use a default (e.g., $u_* = 0.05\,\text{m}\,\text{s}^{-1}$).

---

## 6. TODOs & Open Items

- [ ] **Thermal roughness parameterization:** How does $z_{0h}$ depend on $z_0$, $Re$, or other parameters?
- [ ] **Bulk coefficient validation:** Compare $C_D, C_H$ from scheme to observational databases
- [ ] **Free-convection handling:** Implement proper $u_*$ scaling for calm conditions
- [ ] **Iteration robustness:** Add more aggressive guards and fallbacks if convergence fails
- [ ] **Latent heat flux $L_E$:** Extend to moisture by analogy with $H_s$

---

## References

Verified BibTeX keys:

- `Hogstrom1988` in `refs/master_bibliography.bib`
- `Paulson1970` in `refs/master_bibliography.bib`
