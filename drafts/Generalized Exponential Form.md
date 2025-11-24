## 1. 🔍 Generalized Exponential Form

A common approach to capture the rapid growth of stability functions (and eventual collapse) is the use of exponential or power-law forms. Let's use a generalized exponential form, which is often derived from assumed turbulent length scales:

$$\phi_m(\zeta) = (1 + \zeta)^p$$
$$\phi_h(\zeta) = (1 + \zeta)^q$$

Where, for the $\text{SBL}$, we must have $p > 0$ and $q > 0$. For consistency with the log-linear model, typically $p \approx 1$ and $q \approx 1$ for weak stability, but the general form allows $p, q$ to be higher.

### A. Gradient Richardson Number ($\text{Ri}_g$) and $F(\zeta)$

The general form of $\text{Ri}_g(\zeta)$ becomes:

$$\text{Ri}_g(\zeta) = \zeta \cdot F(\zeta) = \zeta \cdot \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2} = \zeta \cdot \frac{(1 + \zeta)^q}{(1 + \zeta)^{2p}}$$

$$\mathbf{\text{Ri}_g(\zeta) = \zeta (1 + \zeta)^{q - 2p}}$$

The stability function $F(\zeta)$ is:

$$\mathbf{F(\zeta) = (1 + \zeta)^{q - 2p}}$$

### B. Curvature Analysis ($\text{Ri}_g^{\prime\prime}$)

For $\text{Ri}_g(\zeta)$ to exhibit the necessary **concave-down curvature** for $\text{McNider}$ bias ($B>1$), the second derivative must be negative ($\text{Ri}_g^{\prime\prime} < 0$).

Since $\text{Ri}_g(\zeta) \approx \zeta$ near $\zeta=0$, we look at the difference in the exponents, $\mathbf{\Delta = q - 2p}$.

* If $\mathbf{q - 2p < 0}$ ($\Delta < 0$): $\text{Ri}_g(\zeta)$ grows sub-linearly, resulting in $\text{Ri}_g^{\prime\prime} < 0$ (concave-down). This ensures **$B > 1$** (bias exists).
* If $q - 2p = 0$: $\text{Ri}_g(\zeta) = \zeta$ (no curvature, no bias, $B=1$).
* If $q - 2p > 0$: $\text{Ri}_g(\zeta)$ grows super-linearly, resulting in $\text{Ri}_g^{\prime\prime} > 0$ (concave-up). This would lead to **$B < 1$** (bulk stability overestimates point stability), which is not typical in the SBL.

The $\text{McNider}$ bias is fundamentally dependent on the turbulent $\text{Prandtl Number}$ ($\text{Pr}_t$) increasing faster than the inverse momentum gradient: $\text{Pr}_t = \phi_m/\phi_h$, which relates to $q$ and $p$.

### C. Bulk Richardson Number ($\text{Ri}_b$) and Bias ($B$)

The bulk $\text{Ri}_b$ requires integrating the $\text{Ri}_g$ form:

$$\text{Ri}_b(\zeta) = \frac{1}{\zeta} \int_0^\zeta \xi (1 + \xi)^{q - 2p}\,d\xi$$

The integral's solvability depends entirely on the value of $\Delta = q - 2p$.

* **Case 1: $q - 2p = -1$** (i.e., $q = 2p-1$): The integral involves $\frac{\xi}{1+\xi}$, which leads to a solvable form involving $\ln(1+\zeta)$, similar to the $\text{Pade}$ approximation.
* **Case 2: $q - 2p \neq -1$**: The integral requires the **Beta Function** or **Hypergeometric functions**, making the analytical expression for $B(\zeta)$ extremely complex.

$$\mathbf{B(\zeta) = \frac{\text{Ri}_g(\zeta)}{\text{Ri}_b(\zeta)} = \frac{\zeta^2 (1 + \zeta)^{q - 2p}}{\int_0^\zeta \xi (1 + \xi)^{q - 2p}\,d\xi}}$$

---

## 2. 🔑 Generalization Takeaway

The analysis shows that the **bias $B(\zeta)$ is inevitable** in any realistic $\text{MOST}$ closure where $q < 2p$ (i.e., where the $\text{Pr}_t$ ratio changes in the $\text{SBL}$).

| Parameter | Function | Interpretation |
| :---: | :--- | :--- |
| $\text{Ri}_g(\zeta)$ | $\zeta \cdot F(\zeta)$ | The local stability being measured. |
| $F(\zeta)$ | $\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}$ | The composite factor that controls stability growth. |
| **Bias $B(\zeta)$** | $\propto \frac{\text{Ri}_g(\zeta)}{\int \text{Ri}_g(\xi)\,d\xi}$ | The fractional error caused by integrating a **concave-down** function over a finite grid box. |
| **Curvature** | Sign of $\frac{d^2\text{Ri}_g}{d\zeta^2}$ | $\mathbf{\text{Concave-Down}}$ is required for $\mathbf{B > 1}$. |