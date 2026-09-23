This dynamical system represents a single-column, zero-dimensional prognostic model of the Atmospheric Boundary Layer (ABL) Turbulent Kinetic Energy (TKE) budget under stable stratification:

$$\frac{dE}{dt} = \mathcal{P} - \mathcal{B} - \varepsilon$$

where $E$ is TKE ($\text{m}^2\text{s}^{-2}$), $\mathcal{P}$ is mechanical shear production, $\mathcal{B}$ is buoyant consumption (loss of TKE to negative buoyancy fluxes), and $\varepsilon$ is molecular viscous dissipation.

### Physical Components and Parameterizations

The model uses a standard 1.5-order closure scheme where vertical turbulent transport is parameterized via an eddy viscosity $K_m = l \sqrt{E+\delta}$:

| Term | Mathematical Form | Atmospheric Physics & Interpretation |
| --- | --- | --- |
| **Mechanical Shear Production ($\mathcal{P}$)** | $l \sqrt{E+\delta} \, S^2$ | Generation of TKE by mean vertical wind shear $S = \Vert{}\partial \bar{\mathbf{u}} / \partial z\Vert{}$. Regularized by $\delta > 0$ to prevent numerical singularities at zero TKE. |
| **Background Buoyant Loss** | $l \sqrt{E+\delta} \, \phi N^2$ | Standard linear damping of TKE by ambient atmospheric stability, governed by the squared Brunt-Väisälä frequency $N^2 = \frac{g}{\theta_0} \frac{\partial \theta}{\partial z}$ and background stability function $\phi$. |
| **Modified Buoyant Destruction ($\mathcal{B}_{\text{mod}}$)** | $l \sqrt{E+\delta} \, D_b(E; N^2)$ | Hump-shaped heat flux saturation function $D_b = c_b N^2 \frac{\alpha E}{(E+\alpha)^2}$ capturing non-linear wave-turbulence interactions and turbulent flux suppression in the Very Stable Boundary Layer (VSBL). |
| **Viscous Dissipation ($\varepsilon$)** | $\frac{E^{3/2}}{l}$ | Kolmogorov isotropic dissipation rate parameterized using the master turbulence mixing length $l$. |

### Boundary-Layer Physics of the Non-Monotone Closure

Standard Monin-Obukhov similarity closures (e.g., Businger-Dyer or Louis functions) assume monotonic buoyant damping ($d\mathcal{B}/dE > 0$). This forces a single continuous equilibrium, preventing models from capturing the sudden **turbulence collapse** observed during night-time transitions.

The modified closure introduces a non-monotonic efficiency peak at $E = \alpha$:

* **Sub-critical Regime ($E \le \alpha$):** As TKE decreases, negative buoyancy fluxes efficiently quench vertical velocity variance ($d D_b / dE > 0$).
* **Super-critical Regime ($E > \alpha$):** At higher TKE, energetic turbulent eddies break up localized stratification structures, suppressing the relative efficiency of buoyant destruction ($d D_b / dE < 0$).

### Meteorological Interpretation of the Bifurcation

The saddle-node fold pair and critical stratification $N^2_{\text{crit}}$ define the regime transition thresholds between two canonical stable boundary layer states:

* **Weakly Stable Boundary Layer (WSBL):** High-TKE, shear-dominated, continuously turbulent equilibrium.
* **Very Stable Boundary Layer (VSBL):** Low-TKE (or quasi-laminar), buoyancy-dominated regime where turbulence is suppressed, causing surface radiation decoupling and runaway nocturnal cooling.

For $N^2 > N^2_{\text{crit}}$, the atmosphere exhibits **bistability and hysteresis**: a strong wind shear event can kick a decoupled, laminar boundary layer into a fully turbulent WSBL state, but once stratification dominates, the system snaps back into decoupling at the lower fold point.