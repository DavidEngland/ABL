## 💡 $\text{McNider}$ Correction: Physical Justification

The $\text{McNider}$ correction factor ($f_c$) is not an arbitrary fit but is mathematically derived from the physical requirement that **integrated turbulent transport** must be independent of the model's vertical grid resolution ($\Delta z$).

### 1. The Governing Constraint

The core principle is expressed by the differential constraint:
$$\frac{d}{d(\Delta z)} \left( f_s(\text{Ri}) \cdot f_c(\text{Ri}, \Delta z) \right) = 0$$

This enforces that the $\mathbf{product}$ of the **reference stability function** ($f_s$) and the **correction factor** ($f_c$) must remain locally **invariant** with respect to changes in layer thickness $\Delta z$. This ensures that the overall effect of turbulence closure remains consistent when refining the grid.

---

### 2. The Logarithmic ODE

Applying the logarithmic derivative to the constraint yields the governing Ordinary Differential Equation (ODE):

$$\frac{d \ln f_c}{d \ln \Delta z} = -\alpha (B-1) \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q$$

This ODE explicitly states that the **relative change in the correction factor** (LHS) is proportional to the **stability bias** ($B-1$) and scaled by the $\mathbf{local~stability}$ ($\zeta$).

---

### 3. The Closed-Form Solution

Integrating this ODE from the reference boundary condition $f_c(\Delta z_{\text{ref}}) = 1$ naturally results in the **power-law closed-form solution**:

$$\boxed{\displaystyle f_c(\Delta z, \zeta) = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}}$$

#### ➡️ Relationship to the Exponential Form

The power-law solution is mathematically equivalent to the exponential form when utilizing the property $x^y = e^{y \ln x}$:

$$f_c \approx \exp\left[-\alpha (B-1) \ln\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right) \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q\right]$$

The standard **practical exponential form** used in models (often replacing $\ln(\Delta z/\Delta z_{\text{ref}})$ with a simple power $\left(\Delta z/\Delta z_{\text{ref}}\right)^p$) is a **generalized empirical simplification** of this derived closed-form solution, but its origin is firmly rooted in the ODE expressing grid-invariance.
