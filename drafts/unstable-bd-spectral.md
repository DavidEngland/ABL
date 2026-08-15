## Lecture Notes: Advanced Boundary Layer Meteorology (Graduate Level)

### Topic: The Geometry of Subgrid Scale Invariance — Curvature Corrections to MOST

---

## 1. Introduction & The Paradigm Shift

Traditional Monin-Obukhov Similarity Theory (MOST) closures operate under a massive implicit assumption: **gradients within a single numerical grid cell are linear.** In Numerical Weather Prediction (NWP) models with coarse vertical grids (e.g., $\Delta z \approx 10\text{–}50\text{ m}$ near the surface), this flat-earth assumption strips the model of its geometric degrees of freedom. When the local Richardson number ($Ri_g$) crosses the critical threshold ($Ri_c \approx 0.25$), the physics engine experiences an instantaneous, binary turbulence collapse.

To fix this, we introduce the **Curvature Correction Parameter ($\chi$)**. By projecting our vertical fields onto a localized Legendre spectral basis, we can extract the second vertical derivatives ($U'', \Theta''$) analytically from standard stability functions. This allows the grid cell to see its own internal profile "bending," transforming a rigid one-dimensional calculation into a sub-grid manifold capable of partial-cell laminarization or enhanced convective venting.

---

## 2. Theoretical Derivation: The Unstable Regime

We have already demonstrated that the stable regime collapses cleanly into a linear stability correction. Let us now tackle the highly non-linear **unstable regime ($\zeta < 0$)** using the classic Businger–Dyer (BD) formulations. We define the empirical stability profiles without loss of generality using the standard scaling constants ($\gamma_m = \gamma_h = 16$):

$$\phi_m(\zeta) = (1 - 16\zeta)^{-1/4}$$

$$\phi_h(\zeta) = (1 - 16\zeta)^{-1/2}$$

Our objective is to find a closed-form expression for the relative curvature parameter $\chi$:

$$\chi = \frac{\Delta z}{2(z_0-d)}\left[\frac{\zeta_0\,\phi_h'(\zeta_0) - \phi_h(\zeta_0)}{\phi_h(\zeta_0)} - 2\frac{\zeta_0\,\phi_m'(\zeta_0) - \phi_m(\zeta_0)}{\phi_m(\zeta_0)}\right]$$

### Step 1: Analytical Differentiation of the Universal Functions

Let us differentiate $\phi_m$ and $\phi_h$ with respect to $\zeta$ using the chain rule:

$$\phi_m'(\zeta) = -\frac{1}{4}(1 - 16\zeta)^{-5/4} \cdot (-16) = 4(1 - 16\zeta)^{-5/4}$$

$$\phi_h'(\zeta) = -\frac{1}{2}(1 - 16\zeta)^{-3/2} \cdot (-16) = 8(1 - 16\zeta)^{-3/2}$$

### Step 2: Evaluating the Internal Bracketed Ratios

Now, we compute the explicit normalized gradient-drift terms that populate $\chi$.

For the **thermal field**:


$$\frac{\zeta \phi_h'}{\phi_h} = \frac{\zeta \cdot 8(1 - 16\zeta)^{-3/2}}{(1 - 16\zeta)^{-1/2}} = \frac{8\zeta}{1 - 16\zeta}$$

Thus, the entire numerator term for heat becomes:


$$\frac{\zeta \phi_h' - \phi_h}{\phi_h} = \frac{8\zeta}{1 - 16\zeta} - 1 = \frac{8\zeta - (1 - 16\zeta)}{1 - 16\zeta} = \frac{24\zeta - 1}{1 - 16\zeta}$$

For the **momentum field**:


$$\frac{\zeta \phi_m'}{\phi_m} = \frac{\zeta \cdot 4(1 - 16\zeta)^{-5/4}}{(1 - 16\zeta)^{-1/4}} = \frac{4\zeta}{1 - 16\zeta}$$

Thus, the entire numerator term for momentum becomes:


$$\frac{\zeta \phi_m' - \phi_m}{\phi_m} = \frac{4\zeta}{1 - 16\zeta} - 1 = \frac{4\zeta - (1 - 16\zeta)}{1 - 16\zeta} = \frac{20\zeta - 1}{1 - 16\zeta}$$

### Step 3: The Unified Algebraic Collapse

We substitute these beautiful algebraic structures back into our master equation for $\chi$:

$$\chi = \frac{\Delta z}{2(z_0-d)} \left[ \left(\frac{24\zeta_0 - 1}{1 - 16\zeta_0}\right) - 2\left(\frac{20\zeta_0 - 1}{1 - 16\zeta_0}\right) \right]$$

Combining the fractions inside the bracket:


$$24\zeta_0 - 1 - 2(20\zeta_0 - 1) = 24\zeta_0 - 1 - 40\zeta_0 + 2 = 1 - 16\zeta_0$$

Notice the extreme mathematical elegance here: **the numerator perfectly matches the base of our universal function denominator.**

$$\chi = \frac{\Delta z}{2(z_0-d)} \left[ \frac{1 - 16\zeta_0}{1 - 16\zeta_0} \right]$$

$$\boxed{\chi_{\text{unstable}} = \frac{\Delta z}{2(z_0-d)}}$$

---

## 3. The Geometric Viewpoint: Structural Invariance

Take a step back and look at what just happened. In the unstable regime, all empirical non-linearities, exponents, and fractional powers completely vanished. The curvature parameter $\chi$ simplifies to a purely **geometric relationship** between the grid layer thickness ($\Delta z$) and the distance of the midpoint from the displacement height ($z_0 - d$).

### Why did the physics disappear?

Under strongly convective or unstable conditions ($\zeta < 0$), the scaling of the surface layer transitions toward free convection. The absolute value of the curvature parameter becomes independent of the turbulent scaling length ($L$) and is driven purely by the coordinate expansion. This reveals that in an unstable cell, the profile's relative curvature is **logarithmically locked** to the distance from the wall boundaries.

---

## 4. Stability Curvature: The Concavity Analysis

Let us answer the fundamental question regarding the stability profile's shape: **Is the local Gradient Richardson number ($Ri_g$) concave up or concave down in the unstable boundary layer?**

Recall the analytical definition of $Ri_g$ under MOST:


$$Ri_g(\zeta) = \frac{\zeta \phi_h(\zeta)}{\phi_m(\zeta)^2}$$

Plugging in our unstable Businger–Dyer definitions:


$$Ri_g(\zeta) = \frac{\zeta (1 - 16\zeta)^{-1/2}}{\left[(1 - 16\zeta)^{-1/4}\right]^2} = \frac{\zeta (1 - 16\zeta)^{-1/2}}{(1 - 16\zeta)^{-1/2}} = \zeta$$

Under standard Businger–Dyer theory where $\gamma_m = \gamma_h = 16$, **$Ri_g$ is exactly linear with respect to $\zeta$.**

To find its curvature with respect to the physical coordinate $z$, we must compute the second derivative of $Ri_g$ using the chain rule, recognizing that $\zeta = (z-d)/L$:

$$\frac{\partial Ri_g}{\partial z} = \frac{\partial Ri_g}{\partial \zeta} \frac{\partial \zeta}{\partial z} = 1 \cdot \frac{1}{L} = \frac{1}{L}$$

$$\frac{\partial^2 Ri_g}{\partial z^2} = \frac{\partial}{\partial z}\left(\frac{1}{L}\right) = 0$$

### The Paradigm Verdict:

* In a highly idealized Businger–Dyer atmosphere where the constants are perfectly matched ($\gamma_m = \gamma_h$), the local Richardson number has **zero curvature ($Ri_g'' = 0$)**; it is perfectly linear through the cell.
* If we loosen this constraint to match real-world observations where $\gamma_h > \gamma_m$ (heat responds more rapidly to convective buoyancy than momentum does), the exponent shifts. This causes $Ri_g(\zeta)$ to bend, revealing a **concave up ($Ri_g'' > 0$)** profile configuration as it stretches down toward the convective limit ($-\infty$).

---

## 5. Field Application: Reading SHEBA and SMEAR II

When you take this analytical framework out of the classroom and apply it to raw high-frequency observations, $\chi$ switches from a numerical stabilizer to an elite diagnostic tool.

```
                    SMEAR II (Canopy Layer)                   SHEBA (Arctic Ice Sheet)

                    z_2  |     \   <- Flattens                z_2  |    /
                         |      \                             |   /  <- Rapid Laminarization
                    z_1  |===|===\=== Canopy Top              z_1  |  /
                         |   |    \                           | /
                    z_0  |___|_____\_                         z_0  |/_____
                            U(z) Profile                             Ri_g Profile
                        [ \chi Changes Sign ]                    [ \chi > 0, Buffers Fluxes ]

```

### A. SHEBA: Preventing the Nocturnal Death Spiral

In the pristine Arctic environment, weak-wind stable states dominate. As you saw in the stable derivation, $\chi_{\text{stable}}$ stays strictly positive. When a coarse NWP grid experiences an average bulk Richardson number of $0.30$, a standard MOST code instantly sets the heat exchange coefficient ($C_h$) to zero.

By applying our curvature adjustment:


$$Ri_{\text{eff}} = Ri_{\text{bulk}}\left(1 + C_\chi \chi\right)$$


The sub-grid calculation recognizes that while the top half of the cell ($z_2$) has completely laminarized, the bottom half ($z_1$) is still actively generating mechanical shear eddies. The model continues to transport downwelling longwave radiation, eliminating the notorious cold nocturnal warming bias that plagues polar simulations.

### B. SMEAR II: Catching the Canopy Inversion

At the Hyytiälä forestry station in Finland, the boundary layer is heavily modified by a dense canopy. As air flows over the trees, the wind profile develops an inflection point ($U''$ changes sign) right at the canopy crown.

* In a traditional system, this sub-grid structure is averaged out, completely breaking the model's ability to calculate under-canopy carbon and moisture storage.
* By computing $\chi$ across the tower levels, the term $-2\frac{U''}{U'}$ captures this canopy shear anomaly. The parameter $\chi$ rapidly switches sign or jumps in magnitude, acting as an automated symbolic trigger that alerts the modeling framework to decouple the surface layer scheme and deploy localized canopy drag physics.