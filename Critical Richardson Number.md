That's a profound observation and a point of significant debate in atmospheric turbulence theory!

You are addressing the **critical condition for the cessation of turbulence** in a stably stratified fluid, which is indeed linked to the asymptotic behavior of the Gradient Richardson number, $Ri_g$, as stability becomes infinite ($\zeta \to \infty$).

### Critical Richardson Number ($Ri_c$) and the Limit of $Ri_g(\zeta)$

Your statement that **the critical Richardson number, $Ri_c$, should be the limit of $Ri_g(\zeta)$ as $\zeta$ increases without bound** ($\zeta \to \infty$) is **theoretically correct** under the framework of Monin-Obukhov Similarity Theory (MOST) when applied to the stable boundary layer (SBL).

Here's the breakdown of why this limit is important, how your power-law forms relate to it, and why this is a source of controversy:

---

### 1. The Asymptotic Limit in MOST

In the stable atmospheric surface layer (SBL), the theoretical relationship for the gradient Richardson number is:
$$
Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)}
$$

For **very strong stability** ($\zeta \to \infty$), Monin-Obukhov similarity theory predicts that $\phi_m$ and $\phi_h$ grow linearly with $\zeta$, leading to the **$z$-less stratification** concept:
$$\phi_m(\zeta) \sim \zeta$$
$$\phi_h(\zeta) \sim \zeta$$

If we assume the ratio of the stability functions approaches a constant, $\phi_h(\zeta)/\phi_m(\zeta) \to a_{h}/a_{m} = C_R$ (where $C_R$ is a constant related to the turbulent Prandtl number in this limit):

$$
\lim_{\zeta \to \infty} Ri_g(\zeta) = \lim_{\zeta \to \infty} \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)} \frac{1}{\phi_m(\zeta)} \propto \lim_{\zeta \to \infty} \zeta \cdot C_R \cdot \frac{1}{\zeta} = C_R
$$

The result is a finite, constant value:
$$
\lim_{\zeta \to \infty} Ri_g(\zeta) = Ri_c
$$

This asymptotic limit, $Ri_c$, represents the maximum $Ri_g$ that can be sustained while maintaining **fully developed Kolmogorov turbulence**. Beyond this value, turbulence is suppressed, mixing ceases, and the $z$-less hypothesis breaks down.

---

### 2. Relation to Your Power-Law Forms

Your proposed power-law forms for strong stability are:
$$\phi_m(\zeta) = (1-\beta_m \zeta)^{-\alpha_m}$$
$$\phi_h(\zeta) = (1-\beta_h \zeta)^{-\alpha_h}$$

While these forms are primarily used for unstable or moderate stability, their asymptotic behavior for large positive $\zeta$ (strong stability, where the terms become large) would be:
$$\phi_m(\zeta) \sim (\beta_m \zeta)^{-\alpha_m}$$
$$\phi_h(\zeta) \sim (\beta_h \zeta)^{-\alpha_h}$$

Substituting these into $Ri_g$:
$$
Ri_g(\zeta) \propto \zeta \frac{(\beta_h \zeta)^{-\alpha_h}}{[(\beta_m \zeta)^{-\alpha_m}]^2} = \zeta \frac{(\beta_h)^{-\alpha_h}}{(\beta_m)^{-2\alpha_m}} \zeta^{(2\alpha_m - \alpha_h)}
$$
$$
Ri_g(\zeta) \propto \text{Constant} \cdot \zeta^{(1 + 2\alpha_m - \alpha_h)}
$$

For $Ri_g(\zeta)$ to approach a **finite, non-zero limit ($Ri_c$ as $\zeta \to \infty$)**, the exponent must be zero:
$$
1 + 2\alpha_m - \alpha_h = 0 \quad \implies \quad \alpha_h = 2\alpha_m + 1
$$

If your power-law forms are intended to cover the strong stability regime, the exponents must satisfy this constraint to recover a finite $Ri_c$.

---

### 3. The Physical/Empirical Controversy 🤯

While theory requires $Ri_g(\zeta) \to Ri_c$ as $\zeta \to \infty$, the actual value and significance of $Ri_c$ are debated:

* **Classical Linear Stability Theory:** The foundational theory for laminar shear flows (Miles, Howard) yields $Ri_c = \mathbf{0.25}$. Below this value, Kelvin-Helmholtz instabilities can grow.
* **Empirical/MOST Limit:** Measurements in the SBL often show that **Kolmogorov-type turbulence breaks down** when the flux Richardson number ($Rf$) or $Ri_g$ exceeds values in the range $\mathbf{0.20}$ to $\mathbf{0.25}$. This is often used as the practical limit for MOST applicability.
* **Observations of Strong Stability:** Some field data and recent research suggest that turbulence can survive, albeit intermittently and in highly anisotropic form, at $Ri$ values **much greater than 0.25** (even $Ri > 1$). In this regime, the traditional MOST framework breaks down, and the flow is no longer considered "fully turbulent."

**Conclusion:** Your analytical structure is a robust way to enforce physical constraints. By requiring the exponents of your stability functions to satisfy $\mathbf{\alpha_h = 2\alpha_m + 1}$, you ensure that the *parameterization itself* respects the theoretical **asymptotic limit** $Ri_g \to Ri_c$ required for consistency with the overall MOST framework.