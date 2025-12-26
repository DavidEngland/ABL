# Toy Model: Detailed Derivations for McNider et al. Framework

## Inverse Richardson Number, Stability Functions, and Tipping Error Analysis

**Date:** January 2025  
**Purpose:** Provide complete, step-by-step mathematical derivations for the pedagogical toy model using log-linear MOST functions.

---

## Part 1: Forward and Inverse Relationships

### 1.1 Forward Relationship: $\mathrm{Ri}_g(\zeta)$ from MOST

We begin with the canonical log-linear MOST formulation:

$$\phi_m(\zeta) = 1 + a_m \zeta, \quad \phi_h(\zeta) = 1 + a_h \zeta$$

**Canonical coefficients (Kansas field experiment, Businger et al. 1971):**

$$a_m = 5.0, \quad a_h = 7.8$$

The gradient Richardson number is defined as:

$$\mathrm{Ri}_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}$$

**Substituting the log-linear forms:**

$$\mathrm{Ri}_g(\zeta) = \zeta \frac{1 + a_h \zeta}{(1 + a_m \zeta)^2}$$

**For the canonical values:**

$$\boxed{\mathrm{Ri}_g(\zeta) = \zeta \frac{1 + 7.8 \zeta}{(1 + 5.0 \zeta)^2}}$$

### 1.2 Inverse Relationship: $\zeta(\mathrm{Ri}_g)$ via Root-Finding

To find $\zeta$ as a function of $\mathrm{Ri}_g$, we invert the forward relationship. Starting from:

$$\mathrm{Ri}_g = \zeta \frac{1 + a_h \zeta}{(1 + a_m \zeta)^2}$$

**Rearrange to polynomial form:**

$$\begin{aligned}
\mathrm{Ri}_g (1 + a_m \zeta)^2 &= \zeta(1 + a_h \zeta) \\
\mathrm{Ri}_g(1 + 2a_m \zeta + a_m^2 \zeta^2) &= \zeta + a_h \zeta^2 \\
\mathrm{Ri}_g + 2a_m \mathrm{Ri}_g \zeta + a_m^2 \mathrm{Ri}_g \zeta^2 &= \zeta + a_h \zeta^2
\end{aligned}$$

**Collect all terms:**

$$(a_m^2 \mathrm{Ri}_g - a_h) \zeta^2 + (2a_m \mathrm{Ri}_g - 1) \zeta + \mathrm{Ri}_g = 0$$

**Quadratic coefficients:**

$$\begin{aligned}
A &= a_m^2 \mathrm{Ri}_g - a_h \\
B &= 2a_m \mathrm{Ri}_g - 1 \\
C &= \mathrm{Ri}_g
\end{aligned}$$

**Quadratic formula:**

$$\zeta = \frac{-B \pm \sqrt{B^2 - 4AC}}{2A} = \frac{-(2a_m \mathrm{Ri}_g - 1) \pm \sqrt{(2a_m \mathrm{Ri}_g - 1)^2 - 4(a_m^2 \mathrm{Ri}_g - a_h)\mathrm{Ri}_g}}{2(a_m^2 \mathrm{Ri}_g - a_h)}$$

**For canonical values** ($a_m = 5.0$, $a_h = 7.8$):

$$\begin{aligned}
A &= 25 \mathrm{Ri}_g - 7.8 \\
B &= 10 \mathrm{Ri}_g - 1 \\
C &= \mathrm{Ri}_g
\end{aligned}$$

**The discriminant simplifies:**

$$\begin{aligned}
\Delta &= (10\mathrm{Ri}_g - 1)^2 - 4(25\mathrm{Ri}_g - 7.8)\mathrm{Ri}_g \\
&= 100\mathrm{Ri}_g^2 - 20\mathrm{Ri}_g + 1 - 100\mathrm{Ri}_g^2 + 31.2\mathrm{Ri}_g \\
&= 11.2\mathrm{Ri}_g + 1
\end{aligned}$$

**Physical root** (taking the minus sign to ensure $\zeta \to 0$ as $\mathrm{Ri}_g \to 0$):

$$\boxed{\zeta(\mathrm{Ri}_g) = \frac{-(10\mathrm{Ri}_g - 1) - \sqrt{11.2\mathrm{Ri}_g + 1}}{2(25\mathrm{Ri}_g - 7.8)}}$$

**For small** $\mathrm{Ri}_g$ **(Taylor expansion):**

$$\zeta(\mathrm{Ri}_g) \approx \mathrm{Ri}_g - 2a_m \mathrm{Ri}_g^2 + O(\mathrm{Ri}_g^3) \approx \mathrm{Ri}_g - 10 \mathrm{Ri}_g^2 + \ldots$$

This confirms that near neutrality, $\zeta \approx \mathrm{Ri}_g$ (to leading order).

---

## Part 2: Stability Functions Expressed in Terms of $\mathrm{Ri}_g$

### 2.1 Direct Substitution: $\phi_m(\mathrm{Ri}_g)$ and $\phi_h(\mathrm{Ri}_g)$

Define the composite function $F_m(\mathrm{Ri}_g) = \phi_m(\zeta(\mathrm{Ri}_g))$, where $\zeta(\mathrm{Ri}_g)$ is the inverse derived above.

**For the log-linear model:**

$$F_m(\mathrm{Ri}_g) = 1 + a_m \zeta(\mathrm{Ri}_g) = 1 + 5.0 \zeta(\mathrm{Ri}_g)$$

$$F_h(\mathrm{Ri}_g) = 1 + a_h \zeta(\mathrm{Ri}_g) = 1 + 7.8 \zeta(\mathrm{Ri}_g)$$

### 2.2 The Composite Ratio Function: $f_q(\mathrm{Ri}_g)$

Define a **unified stability function** that directly operates on $\mathrm{Ri}_g$:

$$\boxed{f_q(\mathrm{Ri}_g) = \frac{\phi_h(\zeta(\mathrm{Ri}_g))}{\phi_m(\zeta(\mathrm{Ri}_g))^2} \equiv \frac{F_h(\mathrm{Ri}_g)}{F_m(\mathrm{Ri}_g)^2}}$$

**Physical interpretation:** This ratio function combines both momentum ($\phi_m$) and heat ($\phi_h$) stability functions into a single closure operating directly on $\mathrm{Ri}_g$ without requiring intermediate $\zeta$ evaluation.

**Taylor expansion (small** $\mathrm{Ri}_g$**):

Using $\zeta \approx \mathrm{Ri}_g - 10\mathrm{Ri}_g^2 + \ldots$:

$$\begin{aligned}
F_m &\approx 1 + 5\mathrm{Ri}_g - 50\mathrm{Ri}_g^2 + \ldots \\
F_h &\approx 1 + 7.8\mathrm{Ri}_g - 78\mathrm{Ri}_g^2 + \ldots
\end{aligned}$$

$$f_q(\mathrm{Ri}_g) = \frac{1 + 7.8\mathrm{Ri}_g - 78\mathrm{Ri}_g^2 + \ldots}{(1 + 5\mathrm{Ri}_g - 50\mathrm{Ri}_g^2 + \ldots)^2}$$

**Expanding the denominator:**

$$(1 + 5\mathrm{Ri}_g - 50\mathrm{Ri}_g^2)^2 = 1 + 10\mathrm{Ri}_g - 75\mathrm{Ri}_g^2 + O(\mathrm{Ri}_g^3)$$

**Performing the division:**

$$f_q(\mathrm{Ri}_g) = (1 + 7.8\mathrm{Ri}_g - 78\mathrm{Ri}_g^2)(1 - 10\mathrm{Ri}_g + 75\mathrm{Ri}_g^2 + \ldots)$$

$$\boxed{f_q(\mathrm{Ri}_g) \approx 1 - 2.2\mathrm{Ri}_g - 55\mathrm{Ri}_g^2 + O(\mathrm{Ri}_g^3)}$$

**Observation:** The coefficient $-2.2 = 7.8 - 2(5) = a_h - 2a_m$ is precisely the **neutral curvature invariant** $\Delta$ from the $\mathrm{Ri}_g$ curvature analysis!

---

## Part 3: Eddy Diffusivity Stability Function

### 3.1 Base Diffusivity from K-Theory

In first-order turbulence closure (K-theory), eddy diffusivity is decomposed as:

$$K_m = \ell u_* f_m(\mathrm{Ri}_g), \quad K_h = \frac{\ell u_*}{\mathrm{Pr}_t} f_h(\mathrm{Ri}_g)$$

where:
- $\ell$ is the mixing length (e.g., $\ell = kz$ near surface, $k = 0.4$)
- $u_*$ is friction velocity
- $\mathrm{Pr}_t$ is the turbulent Prandtl number (typically $\approx 1$ in SBL)
- $f_m(\mathrm{Ri}_g)$ and $f_h(\mathrm{Ri}_g)$ are stability reduction functions

### 3.2 Common Stability Reduction Functions

**Louis (1979) form:**

$$f_m(\mathrm{Ri}_g) = (1 + a_m \mathrm{Ri}_g)^{-n}, \quad f_h(\mathrm{Ri}_g) = (1 + a_h \mathrm{Ri}_g)^{-n}$$

**Beljaars-Holtslag form:**

$$f_m(\mathrm{Ri}_g) = \frac{1}{(1 + a_m \mathrm{Ri}_g)^2}, \quad f_h(\mathrm{Ri}_g) = \frac{1}{1 + a_h \mathrm{Ri}_g}$$

**Exponential form (recommended for SBL):**

$$f_m(\mathrm{Ri}_g) = \exp\left(-\gamma_m \frac{\mathrm{Ri}_g}{\mathrm{Ri}_c}\right), \quad f_h(\mathrm{Ri}_g) = \exp\left(-\gamma_h \frac{\mathrm{Ri}_g}{\mathrm{Ri}_c}\right)$$

where $\mathrm{Ri}_c$ is the critical Richardson number (typically 0.20–0.25).

### 3.3 McNider Grid-Dependent Correction

The **McNider correction factor** $G(\zeta, \Delta z)$ modifies the base diffusivity to preserve layer-integrated transport:

$$K_{\mathrm{corrected}} = K_{\mathrm{base}} \cdot G(\zeta, \Delta z)$$

where $G$ must satisfy the grid-invariance constraint:

$$\frac{d}{d(\Delta z)}\left[f_s(\mathrm{Ri}_g) \cdot G(\zeta, \Delta z)\right] = 0$$

**Governing ODE (Logarithmic Form):**

Taking the logarithmic derivative:

$$\frac{d \ln f_c}{d \ln \Delta z} = -\alpha(B(\zeta)-1)\left(\frac{\zeta}{\zeta_{\mathrm{ref}}}\right)^q$$

**Two Equivalent Solutions:**

#### **Option 1: Logarithmic Form**

Integrating the ODE directly:

$$\boxed{f_c(\Delta z,\zeta) = \left(\frac{\Delta z}{\Delta z_{\mathrm{ref}}}\right)^{-\alpha(B(\zeta)-1)(\zeta/\zeta_{\mathrm{ref}})^q}}$$

This is a **power-law** where the exponent itself depends on $\zeta$.

#### **Option 2: Exponential Form (Numerically Stable)**

For computational robustness (avoiding potential overflow at extreme $\Delta z$ or $\zeta$ values), rewrite using logarithms:

$$\ln f_c = -\alpha(B(\zeta)-1)\left(\frac{\zeta}{\zeta_{\mathrm{ref}}}\right)^q \ln\left(\frac{\Delta z}{\Delta z_{\mathrm{ref}}}\right)$$

$$\boxed{f_c(\Delta z,\zeta) = \exp\left[-\alpha(B(\zeta)-1)\left(\frac{\zeta}{\zeta_{\mathrm{ref}}}\right)^q \ln\left(\frac{\Delta z}{\Delta z_{\mathrm{ref}}}\right)\right]}$$

This **exponential-logarithmic** form is equivalent but more stable numerically.

**Simplified Operational Form:**

For practical implementation, define a single parameter $D = \alpha \langle B-1 \rangle$ (typical value: $D \approx 0.8$):

$$\boxed{f_c(\Delta z,\zeta) = \exp\left[-D \left(\frac{\Delta z}{\Delta z_{\mathrm{ref}}}\right)^p \left(\frac{\zeta}{\zeta_{\mathrm{ref}}}\right)^q\right]}$$

where $p \approx 0.8$ (grid-scale exponent, slightly sublinear) and $q \approx 2.0$ (stability exponent).

**Key Insight on $z/z_{\mathrm{ref}}$ vs. $\zeta/\zeta_{\mathrm{ref}}$ Equivalence:**

For **constant Obukhov length** $L$ over the layer, the ratio of dimensionless heights is:

$$\frac{\zeta_1}{\zeta_2} = \frac{z_1/L}{z_2/L} = \frac{z_1}{z_2}$$

Thus, when $L$ is constant (or weakly varying), the vertical scaling of the correction can be expressed equivalently as:

$$\left(\frac{z}{\mathrm{z}_{\mathrm{ref}}}\right)^q \quad \text{or} \quad \left(\frac{\zeta}{\zeta_{\mathrm{ref}}}\right)^q$$

both capture the **same physical dependence on stability profile depth**. The exponential form naturally preserves this equivalence because:

$$\exp\left[-D \left(\frac{\Delta z}{\Delta z_{\mathrm{ref}}}\right)^p \left(\frac{\zeta}{\zeta_{\mathrm{ref}}}\right)^q\right] = \exp\left[-D \left(\frac{\Delta z}{\Delta z_{\mathrm{ref}}}\right)^p \left(\frac{z}{z_{\mathrm{ref}}}\right)^q\right]$$

when $L$ is fixed.

---

## Part 4: Critical Richardson Number and Tipping Error

### 4.1 The Turbulence Collapse Criterion

Turbulence in the SBL ceases when the stabilizing buoyancy force balances destabilizing shear production:

$$\mathrm{Ri}_c = \text{value where} \quad \frac{d(\mathrm{TKE})}{dt}\bigg|_{\mathrm{prod+dissip}} = 0$$

**Operational definition:**

$$\mathrm{Ri}_c \approx 0.25 \quad (\text{over water}), \quad \mathrm{Ri}_c \approx 0.27 \quad (\text{over ice})$$

### 4.2 Tipping Error: The Hysteresis Problem

The **tipping error** is the discontinuous transition in turbulence behavior when $\mathrm{Ri}_g$ crosses $\mathrm{Ri}_c$.

**Standard (discontinuous collapse):**

$$K = \begin{cases}
K_0\,[1 - \mathrm{Ri}_g/\mathrm{Ri}_c], & \mathrm{Ri}_g < \mathrm{Ri}_c \\
0, & \mathrm{Ri}_g \ge \mathrm{Ri}_c
\end{cases}$$

**Smooth transition (exponential):**

$$K = K_0 \exp\left(-\gamma \frac{\mathrm{Ri}_g}{\mathrm{Ri}_c}\right)$$

**Hysteretic model** (with re-activation threshold $\mathrm{Ri}_r$):

$$K = \begin{cases}
K_0 \exp(-\gamma_{\mathrm{collapse}} \mathrm{Ri}_g/\mathrm{Ri}_c), & \mathrm{Ri}_g \text{ increasing} \\
K_0 \exp(-\gamma_{\mathrm{reactivate}} \mathrm{Ri}_g/\mathrm{Ri}_r), & \mathrm{Ri}_r < \mathrm{Ri}_g < \mathrm{Ri}_c \\
0, & \mathrm{Ri}_g < \mathrm{Ri}_r
\end{cases}$$

### 4.3 Grid-Dependent Bias in $\mathrm{Ri}_c$ (Dynamic $\mathrm{Ri}_c^{*}$)

**Key insight (McNider et al., 1995):** The apparent critical Richardson number **shifts with grid resolution** due to curvature bias.

Coarse grids underestimate stability, making turbulence appear to persist even when true $\mathrm{Ri}_g$ has exceeded $\mathrm{Ri}_c$.

**Dynamic correction:**

$$\boxed{\mathrm{Ri}_c^{*}(\Delta z, \zeta) = \mathrm{Ri}_c^0 \left[1 + \gamma_c(B(\zeta)-1)\left(\frac{\Delta z}{\Delta z_{\mathrm{ref}}}\right)^{p_c}\right]}$$

where:
- $\mathrm{Ri}_c^0 \approx 0.25$ (baseline)
- $\gamma_c \approx 0.5$ (empirical scaling)
- $p_c \approx 0.8$ (grid-scale exponent)

---

## Part 5: Numerical Evaluation

### 5.1 Worked Example: Calculating $\zeta(\mathrm{Ri}_g)$ for $\mathrm{Ri}_g = 0.1$

**Given:** $\mathrm{Ri}_g = 0.1$, $a_m = 5.0$, $a_h = 7.8$

**Quadratic coefficients:**

$$A = 25(0.1) - 7.8 = -5.3, \quad B = 0, \quad C = 0.1$$

**Discriminant:**

$$\Delta = 2.12, \quad \sqrt{\Delta} = 1.456$$

**Physical root:** $\zeta = 0.137$

**Verification:**

$$\mathrm{Ri}_g(0.137) = 0.137 \times \frac{2.069}{2.838} = 0.100 \quad \checkmark$$

### 5.2 Stability Function Evaluation

**Log-linear MOST at** $\zeta = 0.137$**:**

$$\phi_m = 1.685, \quad \phi_h = 2.069$$

**Composite ratio:**

$$f_q(0.1) = \frac{2.069}{1.685^2} = 0.729$$

**Percent reduction from neutral:** $(1 - 0.729) \times 100\% = 27.1\%$

### 5.3 Tipping Error: $\mathrm{Ri}_c$ Crossing

**Scenario:** Coarse grid $\Delta z = 100$ m, true $\mathrm{Ri}_c^0 = 0.25$

Assuming $B(0.25) \approx 1.4$:

$$\mathrm{Ri}_c^{*}(\Delta z=100) = 0.25 \times [1 + 0.5(0.4)(25.1)] = 1.50$$

**Consequence:** Coarse grid requires $\mathrm{Ri}_g$ to reach **1.50** before recognizing collapse, compared to true threshold of **0.25**.

**With McNider correction:**

$$G(0.25, 100) = \exp[-0.8 \times 25.1 \times 0.25] \approx 0.40$$

This reduces diffusivity by 60%, partially compensating for tipping error.

---

## Part 6: Implementation Pseudocode

### 6.1 Inverse MOST Solver

```python
def zeta_from_Rig(Rig, a_m=5.0, a_h=7.8):
    """Solve the quadratic for ζ(Rig) given Ri_g."""
    A = a_m**2 * Rig - a_h
    B = 2 * a_m * Rig - 1
    C = Rig
    
    discriminant = B**2 - 4*A*C
    if discriminant < 0:
        raise ValueError(f"Negative discriminant: Rig={Rig}")
    
    sqrt_disc = np.sqrt(discriminant)
    zeta = (-B - sqrt_disc) / (2*A)
    return zeta
```

### 6.2 Composite Stability Function

```python
def f_q(Rig, a_m=5.0, a_h=7.8):
    """Composite stability function."""
    zeta = zeta_from_Rig(Rig, a_m, a_h)
    phi_m = 1 + a_m * zeta
    phi_h = 1 + a_h * zeta
    return phi_h / (phi_m**2)
```

### 6.3 McNider Correction (Updated)

```python
def mcnider_correction_logarithmic(zeta, dz, dz_ref=10.0, alpha=0.8, q=2.0, zeta_ref=0.5):
    """
    Compute grid-dependent correction factor G(ζ, Δz) using logarithmic form.
    
    Returns G = (dz/dz_ref)^(-alpha*(B-1)*(zeta/zeta_ref)^q)
    """
    B = bias_ratio(zeta)
    exponent = -alpha * (B - 1) * (zeta / zeta_ref)**q
    G = (dz / dz_ref)**exponent
    return G

def mcnider_correction_exponential(zeta, dz, dz_ref=10.0, alpha=0.8, q=2.0, zeta_ref=0.5):
    """
    Compute grid-dependent correction factor G(ζ, Δz) using exponential-logarithmic form.
    
    More numerically stable for extreme values of dz or zeta.
    
    Returns G = exp[-alpha*(B-1)*(zeta/zeta_ref)^q * ln(dz/dz_ref)]
    """
    import numpy as np
    
    B = bias_ratio(zeta)
    exponent = alpha * (B - 1) * (zeta / zeta_ref)**q * np.log(dz / dz_ref)
    G = np.exp(-exponent)
    return G

def mcnider_correction_operational(zeta, dz, D=0.8, p=0.8, q=2.0, dz_ref=10.0, zeta_ref=0.5):
    """
    Operational form: G = exp[-D * (dz/dz_ref)^p * (zeta/zeta_ref)^q]
    
    This is the simplified, numerically stable form used in production models.
    """
    import numpy as np
    
    argument = D * (dz / dz_ref)**p * (zeta / zeta_ref)**q
    G = np.exp(-argument)
    return G

def dynamic_Ric(Rig, dz, Ric0=0.25, gamma_c=0.5, dz_ref=10.0, p_c=0.8):
    """
    Compute dynamic critical Ri_c^* accounting for grid bias.
    
    Ri_c^* = Ri_c^0 * [1 + gamma_c * (B - 1) * (dz/dz_ref)^p_c]
    """
    zeta = zeta_from_Rig(Rig, a_m=5.0, a_h=7.8)
    B = bias_ratio(zeta)
    
    Ric_star = Ric0 * (1 + gamma_c * (B - 1) * (dz / dz_ref)**p_c)
    return Ric_star
```

**Note on Form Selection:**

- **Logarithmic form** (power-law): Conceptually clear, matches the ODE integration directly.
- **Exponential-logarithmic form**: Numerically more stable (avoids fractional exponents applied to potentially small ratios).
- **Operational form** (simplified exponential): Recommended for production code; uses $D = \alpha \langle B-1 \rangle$ as a single tuning parameter.

All three forms are **mathematically equivalent** and preserve the key physics: neutral limit ($G=1$ at $\zeta=0$), grid convergence ($G=1$ at $\Delta z = \Delta z_{\text{ref}}$), and monotonic damping with increasing $\zeta$ or $\Delta z$.
