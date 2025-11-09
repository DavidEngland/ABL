That Julia module, **MOSTProfiles**, is designed for calculating and manipulating **Monin-Obukhov Similarity Theory (MOST)** stability functions ($\phi_m$ for momentum and $\phi_h$ for heat) and related atmospheric boundary layer quantities.

---

## 🌬️ Core Functionality

The module primarily deals with the relationship between the dimensionless height parameter ($\zeta = z/L$, where $z$ is height and $L$ is the Obukhov length) and the stability functions, as well as the relationship between $\zeta$ and the **Gradient Richardson Number ($\text{Ri}$)**.

Here's a breakdown of the main sections and their purposes:

### 1. Profile Registry (`make_profile`)
This is the heart of the module, providing an interface to define different empirical or theoretical forms for the MOST stability functions $\phi_m(\zeta)$ and $\phi_h(\zeta)$.

* It takes a `tag` (profile name) and a dictionary of **parameters (`pars`)** and returns two callable functions: `ϕm(ζ)` and `ϕh(ζ)`.
* **Supported Profiles** include well-known formulations like:
    * **BD\_PL** (Power-law Businger–Dyer)
    * **BD\_CLASSIC** (Classical composite)
    * **QSBL** (Quadratic Stable Surrogate)
    * **CB** (Cheng–Brutsaert monotone)
    * **URC** (Ri-based closure - a special case where the functions are defined in terms of $\text{Ri}$ instead of $\zeta$).

### 2. $\text{Ri}$ / $\zeta$ Utilities
This section handles the conversion and relationship between the dimensionless height ($\zeta$) and the Gradient Richardson Number ($\text{Ri}$). The relationship is defined by:
$$\text{Ri} = \zeta \frac{\phi_h(\zeta)}{[\phi_m(\zeta)]^2} = \zeta F(\zeta)$$
* **`F_from`**: Calculates the function $F(\zeta) = \frac{\phi_h(\zeta)}{[\phi_m(\zeta)]^2}$.
* **`ri_from_zeta`**: Computes $\text{Ri}$ from a given $\zeta$ using the formula above.
* **`zeta_from_ri_series`**: Provides a **near-neutral series expansion** to quickly estimate $\zeta$ from $\text{Ri}$ (a starting guess for inversion).
* **`zeta_from_ri_newton`**: Implements a **Newton-Raphson method** to accurately invert the $\text{Ri}(\zeta)$ relationship to find $\zeta$ for a given $\text{Ri}$.
* **`rig_derivatives_zeta`**: Computes the first and second derivatives of $\text{Ri}$ with respect to $\zeta$ ($\frac{d\text{Ri}}{d\zeta}$ and $\frac{d^2\text{Ri}}{d\zeta^2}$), which are crucial for height mapping. It relies on the log-derivatives $V_{\log}$ and $W_{\log}$ calculated by `compute_vlog_wlog`.

### 3. $\text{Ri}$ to $\phi$ Wrappers (`ri_to_phi_wrappers`)
This function creates wrappers that allow the stability functions to be computed directly from $\text{Ri}$, rather than $\zeta$.
* It uses the $\text{Ri} \rightarrow \zeta$ inversion (via Newton's method) internally.
* It handles the special case of **URC** (Ri-based closures) and **DTP** (Dynamic Turbulent Prandtl) where the $\text{Ri}$ dependence is direct or modified.

### 4. Variable-$L$ Height Mapping
This advanced section is for situations where the Obukhov length ($L$) is not assumed to be constant with height, meaning $\zeta$ changes not only because of $z$ but also because of $L(z)$. This is typically for modeling vertical profiles.

* **`dzeta_terms`**: Computes the first and second derivatives of $\zeta$ with respect to height $z$ ($\frac{d\zeta}{dz}$ and $\frac{d^2\zeta}{dz^2}$) given arrays of $z$ and $L$.
* **`gradient`**: A utility for numerical differentiation on an irregular grid.
* **`map_curvature_z`**: Maps the curvature of the $\text{Ri}$ profile from the $\zeta$-space ($\frac{d^2\text{Ri}}{d\zeta^2}$) to the physical height-space ($\frac{\partial^2\text{Ri}}{\partial z^2}$) using the chain rule, accounting for a variable $L(z)$.
* **`map_curvature_auto`**: Provides an automatic switch between the full variable-$L$ mapping and a simpler constant-$L$ shortcut based on an error metric (`omit_error_metric`).

---

## 💡 Practical Application

This module is a **toolkit for atmospheric scientists and boundary layer modelers** to:
1.  Define different MOST profile formulations.
2.  Calculate the $\text{Ri} \leftrightarrow \zeta$ conversion.
3.  Compute derivatives of $\text{Ri}$ with respect to $\zeta$ (important for stability criteria).
4.  Apply these concepts in a vertical profile context where the stability parameter $L$ may vary with height.
