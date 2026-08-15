# Appendix: Curvature-Based Sub-Grid Buffer for the Richardson Number

## A.1 Physical Diagnosis of the Nocturnal "Death Spiral"

The standard local gradient Richardson number is defined as:

$$
Ri = \frac{\frac{g}{\theta_0} \frac{\partial \theta}{\partial z}}{\left(\frac{\partial u}{\partial z}\right)^2}
$$

where $u(z)$ is the horizontal wind speed, $\theta(z)$ is the potential temperature, and $\theta_0$ is a constant reference temperature.

Under strongly stable stratification (e.g., nocturnal or polar regimes), the standard local formulation triggers a numerical pathology via the following sequence:

1. Longwave radiative cooling at the surface rapidly intensifies thermal stability, driving $\partial \theta/\partial z \gg 0$.
2. Vertical momentum mixing is suppressed, causing the horizontal velocity profile to decouple and forcing the local shear gradient to vanish ($\partial u/\partial z \to 0$) at localized structural nodes.
3. As the denominator in the gradient Richardson number collapses to zero, $Ri \to \infty$. Traditional planetary boundary layer (PBL) closures respond by dropping the turbulent exchange coefficients for momentum and heat ($K_m, K_h$) to absolute zero.

This creates an unphysical decoupling loop known as the **"Nocturnal Death Spiral."** Because the land surface is entirely insulated from the warmer atmospheric reservoir aloft, it undergoes runaway radiative cooling.

---

## A.2 Spectral Coefficients and Dimensional Analysis

To construct a non-local regularizing framework, we project $u(z)$ and $\theta(z)$ onto an orthogonal Chebyshev polynomial basis over the transformed canonical interval $x \in [-1,1]$:

$$\begin{equation} f(x) = \sum_{n=0}^{N} c_n T_n(x), \qquad f \in \{u,\theta\} \end{equation}$$

The lower-order spectral modes carry explicit structural definitions of the vertical profile:

* **$c_0$:** The spatial domain mean, yielding dimensions of $[c_0] = [f]$.
* **$c_1$:** The domain-wide linear gradient mode. In physical space, its dimension scales inversely with the layer depth $L$:

$$\begin{equation} [c_1] \sim \frac{[f]}{L} \end{equation}$$


* **$c_2$:** The profile curvature mode, representing secondary derivative bending structures. Its dimensions map to:

$$\begin{equation} [c_2] \sim \frac{[f]}{L^2} \end{equation}$$



To combine these modes within a dimensionally homogeneous closure, we introduce a static characteristic vertical length scale $z_c$ (e.g., the tower observation height or canopy interface height), which yields identical physical dimensions for both terms:

$$\begin{equation} z_c c_2 \sim \frac{[f]}{L} \end{equation}$$

---

## A.3 Definition of the Effective Richardson Number

We exploit this dimensional consistency to define a non-local, curvature-regularized effective Richardson number ($Ri_{\mathrm{eff}}$):

$$\begin{equation}
Ri_{\mathrm{eff}} = \frac{\dfrac{g}{\theta_0} \, \mathrm{sign}(c_{1,\theta}) \sqrt{c_{1,\theta}^2 + \beta_\theta (z_c c_{2,\theta})^2}} {c_{1,u}^2 + \beta_u (z_c c_{2,u})^2}\end{equation}$$

where:

* $c_{1,u}$ and $c_{2,u}$ are the gradient and curvature Chebyshev coefficients of the wind profile,
* $c_{1,\theta}$ and $c_{2,\theta}$ are the gradient and curvature Chebyshev coefficients of the potential temperature profile,
* $z_c$ is the domain characteristic vertical height scale,
* $\beta_u > 0$ and $\beta_\theta > 0$ are positive, dimensionless regularization scaling parameters.

The numerator in Equation $\ref{eq:Ri-eff}$ retains the dimensions of a buoyancy gradient ($[t^{-2}]$), while the denominator retains the dimensions of squared mechanical shear ($[t^{-2}]$), rendering $Ri_{\mathrm{eff}}$ strictly dimensionless.

---

## A.4 Constraint 1: Recovery of Classical MOST in Linear Profiles

Consider a structurally simple, uniform boundary layer whose vertical profiles are perfectly linear. In this regime, the higher-order structural bending terms vanish identically:

$$\begin{equation} c_{2,u} = 0, \qquad c_{2,\theta} = 0 \end{equation}$$

Substituting these conditions into the regularized closure equation $\ref{eq:Ri-eff}$ yields:

$$\begin{equation} Ri_{\mathrm{eff}} = \frac{\dfrac{g}{\theta_0} \, \mathrm{sign}(c_{1,\theta}) \sqrt{c_{1,\theta}^2 + 0}} {c_{1,u}^2 + 0} = \frac{\dfrac{g}{\theta_0} \, \mathrm{sign}(c_{1,\theta}) |c_{1,\theta}|}{c_{1,u}^2} = \frac{\dfrac{g}{\theta_0} c_{1,\theta}}{c_{1,u}^2} \end{equation}$$

Because $c_{1,\theta} \propto \partial \theta / \partial z$ and $c_{1,u} \propto \partial u / \partial z$, the expression collapses back to the standard local gradient definition:

$$\begin{equation} Ri_{\mathrm{eff}} \equiv Ri \end{equation}$$

Thus, the formulation is exactly consistent with classic Monin-Obukhov Similarity Theory (MOST) in unstratified, non-complex linear profiles.

---

## A.5 Constraint 2: Singularity Prevention for Vanishing Local Shear

Next, we evaluate the extreme stable limit where intense local stratification decouples the momentum field, driving the local linear gradient mode to zero:

$$\begin{equation} c_{1,u} \to 0 \end{equation}$$

Under standard local MOST, this condition forces a division-by-zero singularity ($Ri \to \infty$). Applying this limit to the regularized closure $\ref{eq:Ri-eff}$ yields:

$$\begin{equation} \lim_{c_{1,u}\to 0} Ri_{\mathrm{eff}} = \frac{\dfrac{g}{\theta_0} \, \mathrm{sign}(c_{1,\theta}) \sqrt{c_{1,\theta}^2 + \beta_\theta (z_c c_{2,\theta})^2}} {\beta_u (z_c c_{2,u})^2} \end{equation}$$

As long as the vertical layer retains structural macro-curvature ($c_{2,u} \neq 0$), the denominator is bound quadratically away from zero. Physically, a non-zero $c_{2,u}$ indicates that while the net linear slope may momentarily cross zero at a localized grid node, the profile itself is fundamentally non-linear (e.g., due to an internal gravity wave inflection or a low-level jet boundary).

By tuning the positive parameters $\beta_u$ and $\beta_\theta$, the model ensures that:

$$\begin{equation} Ri_{\mathrm{eff}} < Ri_{\mathrm{crit}} \approx 0.25 \end{equation}$$

This bounds the effective Richardson number safely below its critical threshold, preventing the sudden collapse of vertical flux exchange. The curvature mode $c_2$ thus acts as a non-local **sub-grid buffer** that maintains a continuous baseline of turbulent kinetic energy production, keeping $K_m, K_h > 0$ and completely eliminating the Nocturnal Death Spiral.

---

## A.6 Schematic Diagram of the Curvature Buffer Mechanism

```
Height (z)
  ▲
  │                                      / [Blue Profile: Linear]
  │                                     /  c2 ≈ 0
  │                                    /
  │                                   /
z_node ──────────────────────────────● <─── [Inflection Node]
  │                                ╱ │      ∂u/∂z → 0
  │                              ╱   │      Local Ri → ∞
  │                            ╱     │      Ri_eff stays finite via c2 > 0
  │                           │      │
  │                           │      │
  │                            ╲     │
  │                              ╲   │      [Red Profile: Curved]
  │                                ╲ │      c2 > 0
  └──────────────────────────────────┴────────────────► Wind Speed u(z)

```

* **The Blue Profile ($c_2 \approx 0$):** The profile is entirely linear, meaning the sub-grid terms fall away and the closure reduces cleanly to classical MOST.
* **The Red Profile ($c_2 > 0$):** Even though the local gradient hits zero at $z_{\text{node}}$ (yielding a vertical tangent), the global macro-curvature acts as an algebraic buffer that keeps $Ri_{\mathrm{eff}}$ stable and well below the critical threshold.

---

## A.7 Summary

The curvature-based effective Richardson number formulated in Equation $\ref{eq:Ri-eff}$ achieves three critical objectives:

* It systematically recovers the standard local gradient Richardson number in simplified linear profiles, preserving classic MOST logic when complex stratification is absent.
* It regularizes the system against singularities when local wind gradients collapse, provided the broader column retains structural curvature.
* It provides a dimensionally consistent, physically interpretable non-local sub-grid buffer that prevents numerical decoupling in stable layers.

This framework directly leverages the low-order Chebyshev coefficients $c_1$ (gradient) and $c_2$ (curvature) to stabilize boundary-layer closures under extreme atmospheric stability.