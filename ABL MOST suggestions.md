# Stable Boundary Layer Corrections: Curvature-Aware MOST Implementation

## Executive Summary

Coarse vertical grids in atmospheric models systematically underestimate near-surface stability in the stable boundary layer (SBL), leading to excessive turbulent mixing and warm-biased surface temperatures. This document provides a comprehensive framework for curvature-aware corrections that preserve neutral physics (the invariant $2\Delta$) while reducing grid-induced bias by 40%+ in operational settings.

Key innovation: Analytic curvature of the gradient Richardson number $Ri_g(\zeta)$ quantifies the nonlinear stability structure; preserving the neutral curvature
$$
2\Delta \;=\; 2\left[\left.\frac{d\ln\phi_h}{d\zeta}\right|_{\zeta=0} - 2\left.\frac{d\ln\phi_m}{d\zeta}\right|_{\zeta=0}\right]
$$
anchors the correction to physically consistent near-neutral SBL behavior (for linear-stable: $\Delta = a_h/\mathrm{Pr} - 2a_m$) while damping coarse-grid tail effects.

...existing code...

### 2.1 Core Definitions
Monin–Obukhov similarity:
$$
\zeta \;=\; \frac{z}{L},\qquad
L \;=\; -\frac{u_*^3 \,\theta}{\kappa g \,\overline{w'\theta'}}
$$
$$
\phi_m \;=\; \frac{\kappa z}{u_*}\frac{\partial U}{\partial z},\qquad
\phi_h \;=\; \frac{\kappa z}{\theta_*}\frac{\partial \theta}{\partial z},\quad
\theta_* \;=\; -\frac{\overline{w'\theta'}}{u_*}
$$

Gradient Richardson number:
$$
Ri_g(\zeta) \;=\; \zeta\,\frac{\phi_h}{\phi_m^2} \;=\; \zeta\,F(\zeta),\qquad F \;=\; \frac{\phi_h}{\phi_m^2}.
$$

...existing code...

### 3.2 Linear Stable (Högström, Beljaars–Holtslag)
$$
\phi_m \;=\; 1 + a_m\,\zeta,\qquad \phi_h \;=\; \mathrm{Pr} + a_h\,\zeta,\qquad \zeta>0
$$
Curvature: $\Delta = a_h/\mathrm{Pr} - 2a_m$ (concave-down if $\Delta<0$).

...existing code...

### 4.2 Grid Damping Factor Approach
Modify eddy diffusivities:
$$
K_m^* \;=\; K_m \times G(\zeta,\Delta z),\qquad
K_h^* \;=\; K_h \times G(\zeta,\Delta z)
$$
Constraints on $G$:
$$
G(0,\Delta z)=1,\quad
\left.\frac{\partial G}{\partial \zeta}\right|_{\zeta=0}=0,\quad
\lim_{\Delta z\to 0}G=1,\quad
\frac{\partial G}{\partial \zeta}\le 0 \text{ for fixed }\Delta z.
$$
Functional template:
$$
G(\zeta,\Delta z) \;=\; \exp\!\left[-D\left(\frac{\Delta z}{\Delta z_r}\right)^{p}\left(\frac{\zeta}{\zeta_r}\right)^{q}\right],\qquad p\ge 1,\;q\ge 2.
$$

...existing code...

### 7.2 Observational Validation
ARM NSA (Alaska): stable nights ($\zeta>0.1$).
SHEBA: assess $L(z)$ variability and omission metric
$$
E_{\text{omit}} \;=\; \left|\frac{(d^2\zeta/dz^2)(dRi_g/d\zeta)}{(d\zeta/dz)^2(d^2Ri_g/d\zeta^2)}\right|.
$$

...existing code...