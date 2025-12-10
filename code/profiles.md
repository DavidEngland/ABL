This module provides utility functions for working with **Monin-Obukhov Similarity Theory (MOST)** stability functions, $\phi_m(\zeta)$ (momentum) and $\phi_h(\zeta)$ (heat), and for transforming between the **nondimensional height** $\zeta = z/L$ and the **gradient Richardson number** $\text{Ri}_g$.

---

## 🌬️ Core Concepts

The module is built around two key equations that define the relationship between the profiles and the Richardson number:

1.  **MOST Stability Functions ($\phi_m, \phi_h$):** These functions describe how turbulence profiles deviate from the neutral (logarithmic) profile as atmospheric stability changes (represented by $\zeta$).
    * $\phi_m(\zeta) = \frac{k z}{u_*} \frac{\partial \bar{u}}{\partial z}$ (Momentum)
    * $\phi_h(\zeta) = \frac{k z}{\theta_*} \frac{\partial \bar{\theta}}{\partial z}$ (Heat)
    * *Where $\phi$ often follows a **power-law baseline**: $\phi(\zeta) = (1 - \beta \zeta)^{-\alpha}$*

2.  **Gradient Richardson Number ($\text{Ri}_g$):** This number relates the buoyancy (thermal stratification) to the shear (mechanical mixing) and is defined in terms of the stability functions:
    $$\text{Ri}_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2} = \zeta F(\zeta)$$


---

## 🛠️ Key Utilities and Implementation Details

### 1. Profile Generation (`make_profile`)

The `make_profile` function acts as a registry to generate parameterized callable functions for $\phi_m(\zeta)$ and $\phi_h(\zeta)$ based on various well-known atmospheric stability profile families.

| Family Tag | Description | Key Characteristic |
| :--- | :--- | :--- |
| **BD\_PL** | Businger–Dyer Power-Law | Near-neutral behavior, $\phi \propto (1 - \beta \zeta)^{-\alpha}$ |
| **BD\_CLASSIC** | Businger–Dyer Classical | Piecewise; unstable power-law, stable linear |
| **BH91** | Beljaars–Holtslag | Hybrid polynomial + power-law tail for stable conditions |
| **CB** | Cheng–Brutsaert Monotone | Uses $|\zeta|$ for monotone behavior in stable/unstable conditions |
| **URC** | Ri-based Closure | Functions are defined directly in terms of $\text{Ri}_g$, not $\zeta$ |
| **DTP** | Dynamic Turbulent Prandtl | $\phi_h$ is derived from $\phi_m$ via a $\text{Pr}_t(\text{Ri}_g)$ relationship |

### 2. $\text{Ri}_g \leftrightarrow \zeta$ Transformation

The transformation between $\text{Ri}_g$ and $\zeta$ is central, especially the inversion $\zeta \to \text{Ri}_g$.

* **$\zeta \to \text{Ri}_g$ (Direct):** This is straightforwardly computed using the definition: $\text{Ri}_g(\zeta) = \zeta F(\zeta)$.
* **$\text{Ri}_g \to \zeta$ (Inversion):** This is solved iteratively because $\text{Ri}_g(\zeta)$ is a non-linear function of $\zeta$.

    1.  **Initial Guess (Seed):** For small $|\text{Ri}_g|$, the inversion starts with a third-order **near-neutral series expansion** (truncated at $O(\text{Ri}^3)$):
        $$\zeta_0 \approx \text{Ri} - \Delta \text{Ri}^2 + (1.5\Delta^2 - 0.5 c_1) \text{Ri}^3$$
        *Where $\Delta = \alpha_h \beta_h - 2 \alpha_m \beta_m$ and $c_1 = \alpha_h \beta_h^2 - 2 \alpha_m \beta_m^2$ are neutral curvature coefficients derived from the power-law profiles.*

    2.  **Newton Refinement:** The initial guess $\zeta_0$ is refined using **Newton's method** to solve the root-finding problem $f(\zeta) = \zeta F(\zeta) - \text{Ri}_{\text{target}} = 0$.
        * The derivative needed for the Newton step is:
            $$f'(\zeta) = F(\zeta) (1 + \zeta V_{\log})$$
            *Where $V_{\log} = \frac{\phi_h'}{\phi_h} - 2\frac{\phi_m'}{\phi_m}$ is calculated using **central difference approximation** on the logarithm of the ratio $F(\zeta)$.*

### 3. Logarithmic Derivatives

The derivative calculations for the Newton iteration are performed using **logarithmic derivatives** and numerical differentiation (`_central_diff`), which is more stable than directly differentiating $F(\zeta)$.

* The key logarithmic derivative is $V_{\log} = \frac{d}{d\zeta} \left[ \ln(\phi_h) - 2 \ln(\phi_m) \right] = \frac{d}{d\zeta} \left[ \ln(F) \right]$.

### 4. Closures and Machine Learning Integration

The module also provides utilities for creating **Ri-based closures** (functions that predict $\phi$ directly from $\text{Ri}_g$, like the `ri_closure_series` and `ri_closure_pade`) and includes hooks for integrating **Machine Learning (ML) correction factors** $G(\zeta, \Delta z)$ via ONNX models to damp or modify the profiles based on context or scale.

---

I can provide the specific $\phi_m$ and $\phi_h$ equations for one of the listed profile families, such as **BD\_PL** or **BH91**, if you're interested.