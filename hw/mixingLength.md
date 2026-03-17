In the context of Monin-Obukhov Similarity Theory (MOST), the shear function $\phi_m(\zeta)$ acts as a scaling factor that "stretches" or "shrinks" the effective mixing length relative to its neutral state.

### 1. The Mathematical Derivation
We start by defining the mixing length $l$ in terms of the momentum flux and the vertical wind gradient:
$$u_*^2 = l^2 \left| \frac{\partial \overline{u}}{\partial z} \right|^2 \implies u_* = l \frac{\partial \overline{u}}{\partial z}$$

From the MOST definition of the dimensionless wind shear:
$$\frac{\partial \overline{u}}{\partial z} = \frac{u_*}{\kappa z} \phi_m(\zeta)$$

Substituting this into the mixing length equation:
$$u_* = l \left( \frac{u_*}{\kappa z} \phi_m(\zeta) \right)$$

Solving for $l$ gives the **stability-corrected mixing length** ($l_\zeta$):
$$l_\zeta = \frac{\kappa z}{\phi_m(\zeta)}$$

---

### 2. Physical Influence of $\phi_m$

The shear function $\phi_m$ essentially modulates the efficiency of the "mechanical stirring" at height $z$:

* **Unstable Conditions ($\zeta < 0$):**
    Here, $\phi_m < 1$. Consequently, $l_\zeta > \kappa z$.
    **Interpretation:** Buoyancy assists the vertical movement of air parcels, effectively allowing eddies to be "larger" and more efficient at transporting momentum than they would be by shear alone.

* **Neutral Conditions ($\zeta = 0$):**
    Here, $\phi_m = 1$. Consequently, $l_\zeta = \kappa z$.
    **Interpretation:** Turbulence is purely mechanical; the size of the eddies is limited only by the geometric distance to the boundary.

* **Stable Conditions ($\zeta > 0$):**
    Here, $\phi_m > 1$. Consequently, $l_\zeta < \kappa z$.
    **Interpretation:** Negative buoyancy (stratification) acts as a cap, suppressing vertical motion. Eddies are "flattened" or shrunken, making them less efficient at transporting momentum.



---

### 3. Impact on your Numerical Model
When you implement your black box routine, the mixing length becomes the "knob" that controls the eddy diffusivity $K_m$. Since $K_m \approx u_* l_\zeta$, the $\phi_m$ function is what prevents your model from over-mixing in stable conditions (like a cold night over sea ice) and ensures it captures the vigorous mixing of a convective land slab.



---

## Problem Set: Mixing Length Theory & Numerics

**Context:** The Prandtl mixing length $l_m = \kappa(z - d)$ is the **foundation** of surface-layer turbulence closure. In this problem set, you will:
1. Derive $l_m$ from first principles (neutral layer)
2. Extend it with Blackadar's composite form to account for BL height
3. Apply stability corrections via MOST
4. Compute eddy diffusivity $K_M$ and see how grid spacing affects it

**Given constants:** $\kappa = 0.40$, $\nu = 1.5 \times 10^{-5}\,\text{m}^2\,\text{s}^{-1}$, $g = 9.81\,\text{m}\,\text{s}^{-2}$, $T_0 = 288\,\text{K}$.

---

## Problem 1: Prandtl Mixing Length in the Neutral Boundary Layer (25 pts)

**Background:** In neutral conditions (zero buoyancy, $\zeta = 0$), the wind profile is purely logarithmic:
$$U(z) = \frac{u_*}{\kappa}\ln\left(\frac{z}{z_0}\right)$$

The shear is:
$$\frac{\partial U}{\partial z} = \frac{u_*}{\kappa z}$$

Turbulent momentum flux is carried by eddies of size $l_m$. The Prandtl hypothesis states:
$$u_* = l_m \left|\frac{\partial U}{\partial z}\right|$$

**a) (8 pts)** Substitute the neutral shear expression into Prandtl's eddy-diffusivity form and solve for the mixing length $l_m(z, d)$ in terms of $z$, $d$ (displacement height if present), and $\kappa$. Show that $l_m = \kappa(z-d)$.

**b) (8 pts)** The eddy diffusivity is defined as $K_M = l_m^2 |\partial U/\partial z|$. Substitute your result from (a) and the neutral shear to show that:
$$K_M = \kappa u_* (z - d)$$
(This is equivalent to MOST; good sign!)

**c) (9 pts)** **Numerical example:** Suppose $u_* = 0.3\,\text{m}\,\text{s}^{-1}$, $d = 1\,\text{m}$ (short grass), and $\kappa = 0.40$.
- At $z = 5\,\text{m}$: compute $l_m(5)$ and $K_M(5)$.
- At $z = 20\,\text{m}$: compute $l_m(20)$ and $K_M(20)$.
- Comment: Over what height range is $l_m$ still "close to" the naïve $\kappa z$ estimate (within 10%)? Why might this matter for a model?

---

## Problem 2: Blackadar Composite Mixing Length (30 pts)

**Background:** The Prandtl form $l_m = \kappa(z - d)$ grows without bound, which is unphysical in the PBL. Blackadar (1962) introduced an **asymptotic mixing length** $\lambda$ (related to the BL height) via the reciprocal form:
$$\frac{1}{l_m} = \frac{1}{\kappa(z - d)} + \frac{1}{\lambda}$$

This ensures $l_m \to \lambda$ as $z \to \infty$, preventing unbounded growth.

**a) (10 pts)** **Derivation:** Rearrange the Blackadar reciprocal formula to solve for $l_m$ explicitly:
$$l_m = \frac{\kappa(z - d)}{1 + \kappa(z - d)/\lambda}$$

Verify the two limits:
- $z \to d$ (near surface): $l_m \to ?$ (should recover Prandtl)
- $z \gg \lambda/\kappa$ (free atmosphere): $l_m \to ?$ (should approach $\lambda$)

**b) (12 pts)** **Numerical:** Use $u_* = 0.3\,\text{m}\,\text{s}^{-1}$, $d = 1\,\text{m}$, and assume $\lambda = 0.1 \times h$ where $h = 1000\,\text{m}$ is the PBL height (so $\lambda = 100\,\text{m}$).

| Height $z$ (m) | Prandtl $l_m = \kappa(z-d)$ (m) | Blackadar $l_m$ (m) | Ratio Blackadar/Prandtl |
|---|---|---|---|
| 5 | — | — | — |
| 20 | — | — | — |
| 50 | — | — | — |
| 100 | — | — | — |
| 300 | — | — | — |

Fill in the table. At what height does Blackadar $l_m$ deviate from Prandtl by $> 10\%$? By $> 50\%$?

**c) (8 pts)** The **transition height** is roughly $z_t \approx \lambda/\kappa$. With $\lambda = 100\,\text{m}$ and $\kappa = 0.4$, estimate $z_t$. Does your table from (b) show the transition occurring near this height?

---

## Problem 3: Stability Correction to Mixing Length (28 pts)

**Background:** The stability-corrected mixing length (from the introductory section) is:
$$l_\zeta = \frac{\kappa z}{\phi_m(\zeta)}$$

For this problem, use the **Beljaars–Holtslag** stability function:
- **Unstable** ($\zeta < 0$): $\phi_m = (1 - 16\zeta)^{-1/4}$ (Paulson 1970)
- **Stable** ($\zeta > 0$): $\phi_m = 1 + 5\zeta$ (linear, simple)
- **Neutral** ($\zeta = 0$): $\phi_m = 1$

Recall that $\zeta = (gz/T_0) \cdot ((\theta_i - \theta_{i+1})/u_*^2)$ (gradient Richardson number approximation).

**a) (8 pts)** **Conceptual:** 
- If $\phi_m = (1 - 16\zeta)^{-1/4}$ with $\zeta = -0.01$ (unstable), compute $\phi_m$.
- Then compute $l_\zeta / (\kappa z)$ at that $\zeta$.
- Is $l_\zeta$ larger or smaller than $\kappa z$? Why does this make physical sense?

**b) (10 pts)** **Three scenarios:** At $z = 10\,\text{m}$, $d = 1\,\text{m}$, $\kappa = 0.40$, compute $l_\zeta$ for:

| Scenario | $\zeta$ | $\phi_m$ | $l_\zeta$ (m) |
|---|---|---|---|
| Strong unstable (convection) | $-0.5$ | — | — |
| Neutral | $0.0$ | — | — |
| Strong stable (SBL) | $+0.1$ | — | — |

Interpret: in which scenario is the mixing length largest? Smallest?

**c) (10 pts)** **Application to $K_M$:**
From Problem 1(b), we know $K_M = l_m^2 |\partial U/\partial z|$ where $|\partial U/\partial z| = u_*/(\kappa z) \cdot \phi_m(\zeta)$ (MOST form).

Substituting $l_\zeta = \kappa z / \phi_m(\zeta)$:
$$K_M = \left(\frac{\kappa z}{\phi_m(\zeta)}\right)^2 \cdot \frac{u_*}{\kappa z} \cdot \phi_m(\zeta)$$

Simplify this to show that:
$$K_M = \kappa u_* z$$

**Key insight:** $K_M$ **does not depend on $\zeta$** when written this way. Where did the stability correction go? (Hint: look at the shear term $\partial U/\partial z$ — how does it encode stability?)

---

## Problem 4: Finite Difference Error and Grid-Aware Mixing Length (32 pts)

**Background:** In a numerical model, we don't compute $\partial U/\partial z$ exactly. Instead, we use finite difference (FD):
$$\left(\frac{\partial U}{\partial z}\right)_{\text{FD}} = \frac{\Delta U}{\Delta z}$$

For two model levels $z_k = 5\,\text{m}$ and $z_{k+1} = 10\,\text{m}$ (so $\Delta z = 5\,\text{m}$), the exact log-profile shear is:
$$\left(\frac{\partial U}{\partial z}\right)_{\text{exact}} = \frac{u_*}{\kappa z_{\text{mid}}}$$

where $z_{\text{mid}} = (z_k + z_{k+1})/2 = 7.5\,\text{m}$.

But the FD approximation (using the midpoint rule) gives:
$$\left(\frac{\partial U}{\partial z}\right)_{\text{FD}} \approx \frac{u_*}{\kappa \Delta z} \ln\left(\frac{z_{k+1}}{z_k}\right)$$

**a) (10 pts)** With $z_k = 5\,\text{m}$, $z_{k+1} = 10\,\text{m}$, $u_* = 0.3\,\text{m}\,\text{s}^{-1}$, $\kappa = 0.40$:

- Compute $(\partial U/\partial z)_{\text{exact}}$ using $z_{\text{mid}} = 7.5\,\text{m}$
- Compute $(\partial U/\partial z)_{\text{FD}}$ using the log-ratio formula
- What is the **percent difference**? (Exact is "true"; FD is your approximation.)

**b) (8 pts)** Now suppose you use the **naive FD** instead (linear interpolation without log):
$$\left(\frac{\partial U}{\partial z}\right)_{\text{naive}} = \frac{U(z_{k+1}) - U(z_k)}{\Delta z}$$

With the log profile $U(z) = (u_*/\kappa) \ln(z/z_0)$ and $z_0 = 0.01\,\text{m}$:
- Compute $U(5)$ and $U(10)$.
- Compute the naive FD shear.
- Compare to your exact result from (a). Which FD is better?

**c) (14 pts)** **The mixing length trap:** Suppose you compute the mixing length using Prandtl, $l_m = 0.40 \times (5 - 0) = 2.0\,\text{m}$. Then:
$$K_{M,\text{FD}} = l_m^2 \times (\partial U/\partial z)_{\text{FD}}$$

- Compute $K_{M,\text{FD}}$ using the log-ratio FD shear from (a).
- Compute $K_{M,\text{exact}} = \kappa u_* z_{\text{mid}} = 0.40 \times 0.3 \times 7.5$ (using MOST).
- What is the **percent error** in estimated $K_M$?
- **Physical interpretation:** If your model uses $K_{M,\text{FD}}$ to diffuse temperature, and the true diffusivity is $K_{M,\text{exact}}$, will you over-mix or under-mix heat?

**d) (Bonus, 5 pts)** The **grid-to-mixing-length ratio** is $\Delta z / l_m = 5 / 2 = 2.5$. As this ratio grows (coarser grid), the FD error gets worse. What grid spacing $\Delta z$ would you need (at $z_k = 5\,\text{m}$) to achieve $\Delta z / l_m = 1.0$?

---

## Problem 5: Blackadar + Stability + Grid (35 pts)

**Integration Challenge:** Combine all three: Blackadar mixing length, stability correction, and **exact shear** to compute a grid-aware $K_M$.

**Setup:**
- Surface conditions: $u_* = 0.25\,\text{m}\,\text{s}^{-1}$, $\theta_* = 0.05\,\text{K}$ (slight heating), $T_0 = 288\,\text{K}$, $d = 0.5\,\text{m}$, $\lambda = 150\,\text{m}$ (BL height effects)
- Model levels: $z_k = 3\,\text{m}$, $z_{k+1} = 8\,\text{m}$ (coarse: $\Delta z = 5\,\text{m}$)
- Stability: $\zeta_k = -0.02$ (weakly unstable)
- For this level use $\phi_m(\zeta = -0.02) = (1 - 16 \times (-0.02))^{-1/4} = (1.32)^{-1/4} \approx 0.928$

**a) (8 pts)** Compute the Blackadar mixing length at $z_k = 3\,\text{m}$:
$$l_m = \frac{\kappa(z - d)}{1 + \kappa(z - d)/\lambda}$$

**b) (8 pts)** Compute the stability-corrected effective mixing length:
$$l_{\text{eff}} = \frac{l_m}{\phi_m(\zeta)} = \frac{l_m}{0.928}$$

**c) (10 pts)** Compute the **exact log-ratio shear** (no Taylor approximation):
$$(\partial U/\partial z)_{\text{exact}} = \frac{u_*}{\kappa \Delta z} \ln\left(\frac{z_{k+1}}{z_k}\right) = \frac{0.25}{0.40 \times 5} \ln\left(\frac{8}{3}\right)$$

(Compute $\ln(8/3) \approx 0.981$.)

**d) (9 pts)** Compute the grid-aware $K_M$:
$$K_M = l_{\text{eff}}^2 \times (\partial U/\partial z)_{\text{exact}}$$

Compare to the naïve MOST form $K_M = \kappa u_* z_k = 0.40 \times 0.25 \times 3 = 0.30\,\text{m}^2\,\text{s}^{-1}$. Are they close? Why or why not?

---

**Summary:** By working through this sequence, you have:
1. **Derived** Prandtl's $l_m = \kappa(z-d)$ from momentum conservation
2. **Extended** it with Blackadar's composite form to respect BL structure
3. **Corrected** it with MOST stability functions ($\phi_m$)
4. **Computed** eddy diffusivity and seen how finite-difference errors propagate
5. **Seen** that exact (log-ratio) shear estimates eliminate the dominant FD bias

This foundation underlies the new grid-aware surface-layer scheme in `param/SCAFFOLDING.md` §6.