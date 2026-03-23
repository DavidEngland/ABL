# Regime Classification: Neutral, Unstable, Near-Neutral, Stable, Laminar

**File:** `param/core/regimes.md`  
**Purpose:** Classify flow into stability regime based on Richardson number; select appropriate parameterization.  
**Related:** `richardson.md` (Ri_g computation), `SCAFFOLDING.md` §4 (pseudocode)

Notation follows `CANONICAL_GLOSSARY.md` ($Ri_g$, $Ri_b$, $\zeta = z/L$ with $z \equiv z_{\text{eff}}$ when displacement is used).

---

## 1. Regime Classes

### 1.1 Hierarchical classification

```
Unstable        : Ri_g < -0.05
Near-neutral    : -0.05 <= Ri_g <= +0.05
Weakly stable   : 0.05 < Ri_g < 0.9*Ri_c
Strongly stable : 0.9*Ri_c <= Ri_g < Ri_c
Laminar         : Ri_g >= Ri_c
```

### 1.2 Regime-specific parameter choices

| Regime | $Ri_g$ range | $\phi_m$ form | $K_M$ behavior | Notes |
|--------|---|---|---|---|
| **Unstable** | $< -0.05$ | Paulson $(1-16\zeta)^{-1/4}$ | Amplified by buoyancy | Convection dominates |
| **Near-neutral** | $[-0.05, +0.05]$ | Linear $\phi_m \approx 1$ | Shear-driven | Uncertain; often log-linear |
| **Weak stable** | $[0.05, 0.9 Ri_c)$ | Linear $1 + 5\zeta$ or Holtslag | Modest reduction | Turbulence persists |
| **Strong stable** | $[0.9 Ri_c, Ri_c)$ | Rapid $K_M$ decline | Transitional collapse | Bimodal risk |
| **Laminar** | $\ge Ri_c$ | $K_M = \nu$ | Molecular diffusion only | Fully suppressed turbulence |

---

## 2. Decision Tree for Regime Classification

### 2.1 Flow chart

```
1. Compute Ri_g from local shear & stratification
   ├─ Input: ∂U/∂z, ∂Θ/∂z
   └─ Output: Ri_g (uncorrected)

2. Apply grid-curvature correction
   ├─ Input: Δz, z_k, d
   ├─ Compute bias_factor = (1 - Δz/(2z_eff))²
   └─ Output: Ri_g,corr = Ri_g · bias_factor

3. Check against regime thresholds
   ├─ if Ri_g,corr < -0.05  → UNSTABLE
   ├─ elif -0.05 ≤ Ri_g,corr ≤ +0.05  → NEAR_NEUTRAL
   ├─ elif Ri_g,corr < 0.9·Ri_c  → WEAKLY_STABLE
   ├─ elif Ri_g,corr < Ri_c  → STRONGLY_STABLE
   ├─ elif Ri_g,corr ≥ Ri_c  → LAMINAR
   └─ else  → ERROR (shouldn't happen)

4. Select parameterization
   └─ Based on regime, choose ψ, f_m, K_M floor, etc.
```

### 2.2 Implementation

```python
def classify_regime(Ri_g_corrected, Ri_c=0.2, 
                    thres_unstable=-0.05, thres_neutral=0.05,
                    thres_weak_stable=None):
    """
    Classify flow stability regime.
    
    Args:
        Ri_g_corrected : grid-corrected Richardson number
        Ri_c : critical Richardson number (default 0.2)
        thres_unstable : boundary neutral/unstable (default -0.05)
        thres_neutral : boundary neutral/stable (default +0.05)
        thres_weak_stable : boundary weak/strong stable (default 0.9·Ri_c)
    
    Returns:
        regime : str, one of 'unstable', 'near_neutral', 'weak_stable', 'strong_stable', 'laminar'
    """
    if thres_weak_stable is None:
        thres_weak_stable = 0.9 * Ri_c
    
    if Ri_g_corrected < thres_unstable:
        return 'unstable'
    elif Ri_g_corrected <= thres_neutral:
        return 'near_neutral'
    elif Ri_g_corrected < thres_weak_stable:
        return 'weak_stable'
    elif Ri_g_corrected < Ri_c:
        return 'strong_stable'
    else:
        return 'laminar'
```

---

## 3. Regime-Specific Parameterizations

### 3.1 UNSTABLE

**Characteristics:** Buoyancy-driven turbulence; convection overpowers shear.

**Wind stability function (Paulson 1970):**
$$\phi_m(\zeta) = (1 - 16\zeta)^{-1/4} \quad (\zeta < 0)$$

**Temperature stability function:**
$$\phi_h(\zeta) = (1 - 16\zeta)^{-1/2}$$

**Mixing length:** $l_m = \kappa(z-d) / \phi_m(\zeta) > \kappa(z-d)$ (enhanced by buoyancy)

**Eddy diffusivity:** $K_M = l_m^2 |\partial U/\partial z|$ (large; vigorous mixing)

**Notes:**
- $Ri_g$ can be very negative (e.g., $-1.0$ during strong surface heating)
- Profiles exhibit strong near-surface gradients (plume-like structure)
- This regime is **well-understood** from laboratory and field studies

### 3.2 NEAR-NEUTRAL

**Characteristics:** Marginal stability; uncertain parameterization.

**Empirical approach:** Use a **log-linear blend**:
$$\frac{\partial U}{\partial z} = \frac{u_*}{\kappa z_{\text{eff}}} (1 + c_N \zeta)$$
$$\frac{\partial \Theta}{\partial z} = \frac{\theta_*}{\kappa z_{\text{eff}}} (1 + c_N \zeta)$$

where $c_N \sim 1$–$5$ (observed range; often left as tuning parameter).

**Problem:** Near $\zeta = 0$, MOST exhibits a **discontinuity** or **corner** in slope (Högström 1988, Grachev et al.). A smooth transition is physically uncertain.

**Pragmatic choice:** Use $\phi_m = 1 + 5\zeta$ (linear stable) even for near-neutral unstable (small $|\zeta| < 0.05$). Introduces ~10% error but avoids the singularity.

**Mixing length:** Continue with Blackadar form; no $\phi_m$ damping ($\phi_m \approx 1$).

**Status:** Paper 2 open question: is there observational justification for a specific form near $\zeta = 0$?

### 3.3 WEAKLY STABLE

**Characteristics:** Stratification suppresses turbulence but shear maintains eddies.

**Stability function (linear):**
$$\phi_m(\zeta) = 1 + 5\zeta \quad (\zeta > 0)$$
$$\phi_h(\zeta) = (1 + 5\zeta) P_r^{-1/2} \quad (\text{with thermal Prandtl ratio})$$

Alternatively, more sophisticated (Beljaars–Holtslag, Cheng & Brutsaert):
$$\phi_m(\zeta) = 1 + 5\zeta + a(\zeta - c/d) e^{-d\zeta}$$

**Mixing length:** Three-term reciprocal form (see `mixing_length.md` §3):
$$\frac{1}{l_m} = \frac{1}{\kappa(z-d)} + \frac{1}{\lambda} + \frac{1}{c_s |L|}$$

**Eddy diffusivity:** $K_M = l_m^2 |\partial U/\partial z|$ (reduced vs. neutral)

**Typical scenario:** Clear sky, nighttime over land; wind backed toward geostrophic.

### 3.4 STRONGLY STABLE (near $Ri_c$)

**Characteristics:** Approaching turbulence collapse; strong bimodality risk.

**Mitigation:**
1. Apply **near-$Ri_c$ damping** filter (§5 below)
2. Use enhanced stability function that smoothly approaches zero as $Ri_g \to Ri_c$:

$$f_m = (1 - Ri_g / Ri_c)^2 \quad (\text{England \& McNider 1995})$$

3. Monitor for regime oscillations

**Mixing length:** Aggressively suppressed by the three-term reciprocal form.

**Notes:** Bimodal behavior (flip-flop between turbulent and laminar) is poorly understood; tower data often show it, but unclear if it's real or instrumental artifact.

### 3.5 LAMINAR

**Characteristics:** No turbulence; only molecular diffusion.

**Parameterization:** 
$$K_M = \nu \quad (\text{molecular viscosity})$$

**Temperature gradient:** Can become very steep (strong inversion).

**Wind shear:** Decouples from shear-driven turbulence; geostrophic balance dominates.

**Typical scenario:** Very cold surface (sea ice, high mountains); stable layer may extend 100+ m.

---

## 4. Regime Transition: Damping Near $Ri_c$

### 4.1 The flip-flop problem

Near $Ri_g = Ri_c$, sudden collapse of $K_M$ can cause:
- Temperature to drop rapidly (due to reduced mixing)
- $\partial \Theta/\partial z$ to increase sharply
- $Ri_g$ to jump above $Ri_c$
- Turbulence to vanish entirely
- Subsequently, reduced mixing allows re-stratification
- Cycle repeats (bimodal oscillation)

### 4.2 Damping formulation

Introduce a **smooth transition function** in the range $[0.9 Ri_c, Ri_c]$:

$$f_{\text{damp}} = \max\left( 0, 1 - \frac{Ri_g - 0.9 Ri_c}{0.1 Ri_c} \right)$$

Reduce the iteration step in flux inversion (§8 of SCAFFOLDING.md) when $f_{\text{damp}} < 1$:

```python
if 0.9 * Ri_c <= Ri_g_corr < Ri_c:
    damping_factor = 1.0 - (Ri_g_corr - 0.9*Ri_c) / (0.1*Ri_c)
    # Reduce iteration step:
    du_star_next = damping_factor * du_star_proposed + (1 - damping_factor) * du_star_prev
```

### 4.3 Hysteresis option

**Alternative:** Introduce hysteresis — require $Ri_g$ to rise above $1.1 Ri_c$ to flip to LAMINAR, and drop below $0.8 Ri_c$ to flip back to STRONGLY_STABLE. This mimics observed bimodality but adds memory.

**Status:** TODO — validate against tower data to see if damping or hysteresis is justified.

---

## 5. Summary: Which Regime Am I In?

**Decision table:**

| Scenario | Example conditions | Regime | Typical $K_M$ (m²/s) |
|----------|---|---|---|
| Sunny afternoon, active convection | $T_s > T_2$, $\partial \Theta/\partial z < 0$ | Unstable | 0.5–2.0 |
| Overcast midday, weak wind | $\partial \Theta/\partial z \approx \Gamma_d$, $u_* \sim 0.2$ | Near-neutral | 0.01–0.1 |
| Evening, light wind, clear sky | $\partial \Theta/\partial z > \Gamma_d$ gently, $u_* \sim 0.1$ | Weak stable | 0.001–0.01 |
| Night, strong inversion, calm | $\partial \Theta/\partial z > 0.01$ K/m, $u_* < 0.05$ | Strong stable | 1e-5–1e-4 |
| Arctic winter, sea ice | $\partial \Theta/\partial z > 0.1$ K/m, $Ri_g >> Ri_c$ | Laminar | $\sim \nu = 1.5$e-5$ |

---

## 6. TODOs & Open Items

- [ ] **Near-neutral form:** Resolve discontinuity at $\zeta = 0$; validate against observations
- [ ] **Hysteresis vs. damping:** Implement both; run GABLS1 cases to see which is more stable
- [ ] **Bimodality in tower data:** Examine real datasets near $Ri_c$; is flip-flop real or artifact?
- [ ] **Regime transition smoothness:** Ensure no spurious discontinuities or spikes when switching parameterizations
- [ ] **Cloud-topped BL:** How to classify regime when surface is unstable but clouds suppress convection?

---

## References

Verified BibTeX keys:

- `EnglandMcNider1995BLM` in `refs/master_bibliography.bib`
- `Grachev2005` in `refs/pbl.bib`
- `Hogstrom1988` in `refs/master_bibliography.bib`
