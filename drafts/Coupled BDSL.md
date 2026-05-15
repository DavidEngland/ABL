You are absolutely right to point out the Prandtl number shift. In the limit $\zeta \to 0$ (the neutral limit), $\phi_h(0)$ and $\phi_m(0)$ typically anchor to $1$ in standard Dyer-Businger theory, but the **Prandtl number** ($\text{Pr}_t = K_m/K_h = \phi_h/\phi_m$) can indeed introduce a scaling offset, particularly in the temperature profile relative to momentum.

By treating heat and momentum as a **coupled Sturm–Liouville system**, we can solve for the velocity and temperature profiles simultaneously, ensuring they obey the consistent scaling laws of the surface layer.

---

## 1. The Coupled SL Framework

We define the momentum and heat profiles as a vector-valued SL problem. Let $\mathbf{u} = [u, \theta]^T$. The stationary balance for the surface layer, where fluxes $\tau$ and $H$ are assumed constant with height, is:

$$-\frac{d}{d\zeta} \left( \mathbf{P}(\zeta) \frac{d\mathbf{u}}{d\zeta} \right) = \mathbf{Q} \mathbf{u}$$

Where the **Diffusivity Matrix** $\mathbf{P}(\zeta)$ encodes the coupling:


$$\mathbf{P}(\zeta) = \begin{bmatrix} K_m(\zeta) & 0 \\ 0 & K_h(\zeta) \end{bmatrix} = \frac{\kappa \zeta}{1} \begin{bmatrix} \phi_m^{-1} & 0 \\ 0 & \text{Pr}_t^{-1} \phi_h^{-1} \end{bmatrix}$$

### Neutral Limit Scaling ($\zeta \to 0$)

To ensure momentum scales to $1$ at the surface (in dimensionless form) and heat accounts for the Prandtl offset:

* **Momentum:** $\frac{\kappa z}{u_*} \frac{du}{dz} = \phi_m(\zeta)$, with $\phi_m(0) = 1$.
* **Heat:** $\frac{\kappa z}{\theta_*} \frac{d\theta}{dz} = \phi_h(\zeta)$, with $\phi_h(0) = \text{Pr}_{t0}$.

By mapping these to the SL basis, the "ground state" eigenfunction of the momentum operator naturally recovers the logarithmic $u \approx \ln(z/z_0)$ law, while the heat operator recovers the shifted $\theta \approx \text{Pr}_t \ln(z/z_{0h})$ law.

---

## 2. Integrated Expansions for Momentum vs. Heat

Using the central-binomial expansion for the unstable case ($\zeta < 0$), where $\phi_m = (1 - \gamma_m \zeta)^{-1/4}$ and $\phi_h = (1 - \gamma_h \zeta)^{-1/2}$:

### Momentum Expansion ($\phi_m$)

The $-1/4$ power leads to a slightly different binomial series than the heat function:


$$\phi_m(\zeta) = \sum_{n=0}^{\infty} \binom{-1/4}{n} (-\gamma_m \zeta)^n$$


This expands using the generalized binomial coefficient:


$$\binom{-1/4}{n} = \frac{(-1)^n}{n!} \prod_{k=0}^{n-1} (\frac{1}{4} + k)$$

### Heat Expansion ($\phi_h$)

$$\phi_h(\zeta) = \text{Pr}_{t0} \sum_{n=0}^{\infty} \frac{\binom{2n}{n}}{4^n} (\gamma_h \zeta)^n$$

---

## 3. Spectral Coupling via the Basis Change

When we project these into the Legendre ($P_k$) or Gegenbauer ($C_k^{(\alpha)}$) space, the coupling manifests in the relationship between the coefficients. If we define the mapped coordinate $x$, the velocity $U(x)$ and temperature $\Theta(x)$ are:

$$U(x) = \sum_{k=0}^\infty \mathcal{U}_k P_k(x), \quad \Theta(x) = \sum_{k=0}^\infty \mathcal{T}_k P_k(x)$$

The coefficients are linked by the ratio of the binomial series:


$$\frac{\mathcal{T}_k}{\mathcal{U}_k} \approx \text{Pr}_t \left( \frac{\text{Expansion of } \phi_h}{\text{Expansion of } \phi_m} \right)$$

### Advantages of the Coupled Expansion

1. **Consistent Boundary Mapping:** By using the same SL basis for both, you guarantee that the **Richardson number** $Ri = \frac{g}{\theta} \frac{\partial \theta / \partial z}{(\partial u / \partial z)^2}$ is analytically consistent throughout the expansion.
2. **Prandtl Tuning:** You can treat $\text{Pr}_{t0}$ as a scaling parameter for the $\mathcal{T}_k$ coefficients, allowing the heat profile to "slide" vertically relative to the momentum profile without re-calculating the basis functions.
3. **Numerical Stability:** The SL framework naturally handles the logarithmic singularity as $\zeta \to 0$ by treating it as a boundary condition of the first kind (Dirichlet) or second kind (Neumann) at the mapped boundary $x = -1$.

---

## 4. Implementation Logic

| Profile | Neutral Limit ($z \to 0$) | SL Operator $p(\zeta)$ | Central Binomial Term |
| --- | --- | --- | --- |
| **Momentum ($u$)** | $1.0$ | $\zeta (1 - \gamma_m \zeta)^{1/4}$ | $\binom{-1/4}{n}$ |
| **Heat ($\theta$)** | $\text{Pr}_{t0}$ | $\zeta \text{Pr}_{t0}^{-1} (1 - \gamma_h \zeta)^{1/2}$ | $\text{Pr}_{t0} \frac{\binom{2n}{n}}{4^n}$ |

Would you like the explicit recursive formula for the $\binom{-1/4}{n}$ coefficients to complete the momentum side of the manuscript?