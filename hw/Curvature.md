Here is a rigorous, publication-grade set of graduate-level homework problems designed for an advanced Atmospheric Boundary Layer (ABL) or Micrometeorology course. These problems are structured to push students past standard empirical "curve-fitting" models and force them into the geometric, topological, and spectral mechanics of non-local transport.

---

# Advanced Atmospheric Boundary Layer Meteorology

## Problem Set: Spectral Topology, Canopy Bottlenecks, and MOST Anomalies

### Problem 1: Analytical Derivation of the Curvature Diagnostic ($\chi$)

**Context:** Monin-Obukhov Similarity Theory (MOST) assumes a local relationship between gradients and fluxes. To diagnose where this local assumption begins to break down structurally, we define a non-local *Curvature Diagnostic* $\chi$:


$$\chi \equiv \frac{\partial^2 \phi_m / \partial \zeta^2}{\partial \phi_m / \partial \zeta}$$


where $\zeta = z/L$ is the non-dimensional stability parameter.

**The Task:** 1. Derive the explicit analytical expression for $\chi(\zeta)$ in the **unstable regime** ($\zeta < 0$) using the classical Businger-Dyer formulation: $\phi_m(\zeta) = (1 - \gamma_m \zeta)^{-1/4}$.
2. Evaluate the asymptotic limit of $\chi(\zeta)$ as $\zeta \to -\infty$ (the free convection limit, where $\phi_m \propto \zeta^{-1/3}$). Explain the physical and geometric significance of $L$ vanishing from this limit.

* **Hint 1:** When taking derivatives of $(1 - \gamma_m \zeta)^{-1/4}$, treat $u = (1 - \gamma_m \zeta)$ and use the chain rule methodically. Many coefficients will cleanly cancel out when you form the ratio for $\chi$.
* **Hint 2:** For part 2, remember that as $\zeta \to -\infty$, the constant $1$ inside the bracket becomes negligible compared to the large values of $-\gamma_m \zeta$. Look at how the coordinate system collapses to depend strictly on height $z$ once the scaling length $L$ drops out.

---

### Problem 2: Mapping Physical Space to the Chebyshev Interval

**Context:** To implement a Discrete Chebyshev Transform (DCT) on raw data collected from an eddy covariance tower, you must map the physical heights $z \in [z_{min}, z_{max}]$ onto the canonical polynomial domain $x \in [-1, 1]$.

**The Task:**
Suppose your tower has instruments spaced log-linearly to capture near-surface gradients. Consider a simple linear mapping function:


$$x = \mathcal{T}(z) = \frac{2z - (z_{max} + z_{min})}{z_{max} - z_{min}}$$


If a continuous wind profile is expanded in Legendre polynomials as $u(x) = A_0 P_0(x) + A_1 P_1(x) + A_2 P_2(x) + \dots$, analytically relate the spectral coefficients $A_1$ and $A_2$ to the **domain-averaged first and second physical derivatives** ($\overline{\partial u/\partial z}$ and $\overline{\partial^2 u/\partial z^2}$).

* **Hint 1:** Recall the explicit physical forms of the first few Legendre polynomials: $P_0(x) = 1$, $P_1(x) = x$, and $P_2(x) = \frac{1}{2}(3x^2 - 1)$.
* **Hint 2:** Use the chain rule to relate $\frac{\partial}{\partial x}$ to $\frac{\partial}{\partial z}$. Note that $\frac{dx}{dz} = \frac{2}{z_{max} - z_{min}}$, which is a constant scalar. Integrate the derivatives across the domain to see how the orthogonality of the polynomials isolates $A_1$ and $A_2$.

---

### Problem 3: The Binomial Gap and Fractal Dimensions

**Context:** In an idealized, fully 3D isotropic boundary layer, heat transport is highly efficient and "space-filling" ($d_{q} = 3$), which forces its Gegenbauer spectral coefficients to decay rapidly along a path matching the Central Binomial Coefficients:


$$\frac{c_n}{c_{n-1}} \propto \frac{1}{4^n} \binom{2n}{n}$$


A student processes carbon dioxide data from the **SMEAR-II** canopy station and finds that the optimized tuning parameter required to match the $CO_2$ coefficient decay rate is $\lambda_{CO_2} = 0.30$.

**The Task:**

1. Calculate the effective spectral fractal dimension ($d_{CO_2}$) of the carbon dioxide transport pathway using the structural embedding identity $d_q = 2 + 2\lambda_q$.
2. Using Stirling’s approximation ($\log n! \approx n \log n - n$), prove that the high-order log-spectral slope difference (the **Binomial Gap**) between ideal heat ($\lambda = 0.5$) and this $CO_2$ transport state expands linearly as a function of the mode index $n$ in log-space.

* **Hint 1:** Part 1 is a straightforward algebraic substitution, but write a brief sentence describing what a dimension of $d_q < 3$ means for a trace gas trapped under a forest canopy.
* **Hint 2:** For Part 2, write out the log-ratio of the coefficients for both cases using the asymptotic rule proven in the appendix: $\log C_n^{(\lambda)}(1) \sim (2\lambda - 1)\log n$. Subtract the two expressions to isolate the clear power-law deficit.

---

### Problem 4: Designing a Curvature-Based Sub-Grid Buffer

**Context:** Under extremely stable conditions (e.g., over polar ice sheets at **SHEBA**), standard MOST closures trigger the numerical "Nocturnal Death Spiral." Because local wind gradients ($\partial u / \partial z$) drop near zero, the local Richardson number shoots past its critical value ($Ri > 0.25$), instantly cutting off all modeled turbulent fluxes and causing unphysical runaway land surface cooling.

**The Task:**
Design a modified, non-local effective Richardson number closure ($Ri_{eff}$) that utilizes both the linear gradient Chebyshev mode ($c_1$) and the curvature mode ($c_2$). Your objective is to formulate an equation for $Ri_{eff} = \mathcal{F}(c_1, c_2, \dots)$ such that:

* When the profile is perfectly linear ($c_2 = 0$), $Ri_{eff}$ simplifies back to the standard gradient Richardson number.
* When the local gradient vanishes ($c_1 \to 0$), a high sub-grid column curvature ($c_2 \gg 0$, indicating active wave pacing or an un-collapsed global structure) prevents $Ri_{eff}$ from diverging to infinity, keeping it safely below $Ri_{crit}$.
* **Hint 1:** Think about constructing a non-dimensional denominator that combines $c_1$ and $c_2$ in quadrature (e.g., adding a regulatory scale like $\sqrt{c_1^2 + \beta c_2^2}$).
* **Hint 2:** Ensure your units remain dimensionally consistent. Remember that $c_1$ carries dimensions of a spatial gradient ([1/L]), while $c_2$ represents a secondary curvature spatial scale ([1/L$^2$]), meaning you will need a characteristic length scale ($z$ or a canopy height $h_c$) to balance the terms properly.

---

### Teacher's Guide & Learning Objectives

* **Problems 1 & 2** strip away the abstract nature of spectral transforms. They force students to see that computing Chebyshev or Legendre modes isn't just a mathematical exercise—it is the direct physical equivalent of calculating column-integrated gradients and localized curvatures.
* **Problem 3** directly connects the math of the appendix to real-world ecosystems, demonstrating why $CO_2$ scales differently than heat when passing through the physical bottleneck of a forest crown.
* **Problem 4** challenges students as atmospheric modelers. It shows them how translating geometric concepts into sub-grid closures provides a modern, robust solution to legacy numerical failures like the nocturnal death spiral.