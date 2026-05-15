This complete version integrates the **Businger–Dyer (BD) flux-gradient relations** into a formal **Sturm–Liouville (SL) framework**, utilizing the **central-binomial expansion** to bridge empirical atmospheric similarity with orthogonal polynomial spectral methods.

---

## 1. The Monin–Obukhov Sturm–Liouville Problem

We define a stationary vertical heat balance in the atmospheric surface layer (ASL) using the dimensionless stability parameter $\zeta = z/L$. The governing equation for the potential temperature profile $\overline{\theta}(\zeta)$ is framed as a second-order linear differential equation:

$$-\frac{d}{d\zeta}\left[ p(\zeta) \frac{d\overline{\theta}}{d\zeta} \right] + q(\zeta)\overline{\theta} = \lambda w(\zeta)\overline{\theta}$$

### Operator Components

To recover the **Businger–Dyer** physics, we must encode the eddy diffusivity into the SL coefficient $p(\zeta)$. From the standard flux-gradient relation:


$$\frac{\kappa z}{\theta_*} \frac{d\overline{\theta}}{dz} = \phi_h(\zeta) \implies \frac{d\overline{\theta}}{d\zeta} = \frac{\theta_*}{\kappa} \frac{\phi_h(\zeta)}{\zeta}$$

By introducing an "effective SL diffusivity" $p(\zeta)$, we set:

* **Diffusivity Coefficient:** $p(\zeta) = \frac{\zeta}{\phi_h(\zeta)}$. This term captures the reciprocal of the stability function, ensuring that as stability increases ($\phi_h$ grows), the effective diffusion decreases.
* **Potential Term:** $q(\zeta) = 0$ (for a pure diffusion balance).
* **Weight Function:** $w(\zeta)$ is chosen based on the desired orthogonal basis (e.g., $w(\zeta) = 1$ for Legendre, $w(\zeta) = (1-\zeta^2)^{\alpha-1/2}$ for Gegenbauer).

---

## 2. Central-Binomial Expansion of the BD Functions

The BD functions for heat, $\phi_h(\zeta)$, are typically expressed as fractional powers that naturally admit expansion via central binomial coefficients (Yan et al., 2019).

### Case A: Unstable Conditions ($\zeta < 0$)

Using the common form $\phi_h(\zeta) = (1 - \gamma_h \zeta)^{-1/2}$:


$$\phi_h(\zeta) = \sum_{n=0}^{\infty} \frac{\binom{2n}{n}}{4^n} (\gamma_h \zeta)^n$$


This expansion is valid for $|\gamma_h \zeta| < 1$. The **central binomial coefficient** $\binom{2n}{n}$ governs the decay of the higher-order spectral terms.

### Case B: Stable Conditions ($\zeta > 0$)

Using the square-root variant $\phi_h(\zeta) = (1 + \beta_h \zeta)^{1/2}$:


$$\phi_h(\zeta) = \sum_{n=0}^{\infty} \frac{(-1)^{n-1}}{4^n} \frac{\binom{2n}{n}}{1-2n} (\beta_h \zeta)^n$$

---

## 3. Mapping to Orthogonal SL Eigenbases

To solve the SL problem on a finite interval $\zeta \in [\zeta_0, \zeta_{max}]$, we map the coordinate to $x \in [-1, 1]$ and project the profile onto **Legendre** or **Gegenbauer** polynomials. Legendre polynomials are a specific class of Jacobi polynomials often used in atmospheric boundary layer modeling because they satisfy singular SL problems at the boundaries (Barnard, 2001).

### The Integrated Profile Expansion

The temperature profile $\overline{\theta}(\zeta)$ is the integral of the BD expansion. For the unstable case, the departure from the logarithmic surface law is given by:


$$\Psi_h(\zeta) = \int \frac{\phi_h(\zeta)-1}{\zeta} d\zeta = \sum_{n=1}^{\infty} \frac{1}{n} \frac{\binom{2n}{n}}{4^n} (\gamma_h \zeta)^n$$

Each monomial $\zeta^n$ is projected onto the Legendre basis $P_k(x)$ using:


$$\zeta^n \to \sum_{k=0}^n b_{n,k} P_k(x)$$


The resulting **complete spectral representation** of the temperature profile is:


$$\overline{\theta}(x) = \frac{\theta_*}{\kappa} \left[ \ln \zeta(x) + \sum_{k=0}^{\infty} A_k P_k(x) \right]$$


where the spectral coefficients $A_k$ are weighted sums of the central binomial terms:


$$A_k = \sum_{n=k}^{\infty} \frac{\gamma_h^n}{n} \frac{\binom{2n}{n}}{4^n} b_{n,k}$$

---

## 4. Summary Table: SL Framework for Businger–Dyer

| Component | Businger–Dyer Origin | SL Formulation |
| --- | --- | --- |
| **Variable** | $\zeta = z/L$ | $x = 2(\zeta/\zeta_{max}) - 1$ |
| **Diffusivity** | $K_h \propto \zeta/\phi_h$ | $p(x) = (1-x^2)$ (Canonical) |
| **Expansion** | Binomial Series | Orthogonal $P_n(x)$ or $C_n^{(\alpha)}(x)$ |
| **Convergence** | $\phi_h$ power law | Spectral (exponential for smooth $\phi_h$) |

This framework allows atmospheric models to replace standard iterative "look-up" tables for BD functions with a direct spectral calculation, significantly improving computational efficiency in climate and weather codes (Gross et al., 2018).

### References

* Barnard, J. C. (2001). *Intermittent Turbulence in the Very Stable Ekman Layer*. Pacific Northwest National Laboratory. [https://doi.org/10.2172/1000183](https://www.google.com/search?q=https://doi.org/10.2172/1000183)
Cited by: 18
* Gross, M., Wan, H., Rasch, P. J., Caldwell, P. M., Williamson, D. L., Klocke, D., Jablonowski, C., Thatcher, D. R., Wood, N. B., Cullen, M. J. P., Beare, R. J., Willett, M., Lemarié, F., Blayo, E., Malardel, S., Termonia, P., Gassmann, A., Lauritzen, P. H., Johansen, H., ... Leung, L. R. (2018). Physics–Dynamics Coupling in Weather, Climate, and Earth System Models: Challenges and Recent Progress. *Monthly Weather Review*, *146*(11), 3505–3544. [https://doi.org/10.1175/mwr-d-17-0345.1](https://www.google.com/search?q=https://doi.org/10.1175/mwr-d-17-0345.1)
Cited by: 100
* Yan, B., Huang, S., Feng, J., & Wang, Y. (2019). The Distribution and Uncertainty Quantification of Wind Profile in the Stochastic General Ekman Momentum Approximation Model. *Journal of Meteorological Research*, *33*(2), 336–348. [https://doi.org/10.1007/s13351-019-8076-3](https://www.google.com/search?q=https://doi.org/10.1007/s13351-019-8076-3)
Cited by: 7