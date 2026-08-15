Vertical grid spacing ($\Delta z$) is a critical determinant of the **stability correction factor** ($G$ or $f_c$) because it dictates the magnitude of the **grid-curvature bias** inherent in coarse-resolution models. This bias arises because traditional stability functions within Monin–Obukhov Similarity Theory (MOST) exhibit **concave-down curvature** in the Stable Boundary Layer (SBL), where the neutral curvature invariant $\Delta < 0$.

The relationship between grid spacing and the stability correction factor is governed by the following principles:

### **1. The Mechanism of Grid-Curvature Bias**
When vertical grid spacing is coarse (typically $\Delta z > 20$ m), the model calculates a **bulk Richardson number** ($Ri_b$) by averaging gradients over a finite layer. Because the underlying stability function $Ri_g(\zeta)$ is concave-down, **Jensen's Inequality** dictates that this layer average will systematically **underestimate the true local stability** at the geometric mean height $z_g$. This underestimation leads to **excessive vertical mixing** and a characteristic **warm bias** in surface temperature forecasts.

### **2. Scaling with Resolution**
The correction factor is designed to be **resolution-aware**, ensuring that model results remain consistent across varying grid spacings.

*   **Fine-Grid Convergence:** As grid spacing decreases ($\Delta z \to 0$), the correction factor $G$ (or $f_c$) approaches **unity**, allowing the model to recover the original physical parameterization without modification.
*   **Coarse-Grid Damping:** As $\Delta z$ increases, the correction factor **decreases monotonically**, acting as a damping agent on eddy diffusivities ($K$). This reduction in $K$ counteracts the spurious over-mixing caused by the underestimation of stability on thick grids.

### **3. Functional Form and Feature Importance**
The correction factor is typically expressed in an exponential form that explicitly includes $\Delta z$ as a scaled input:

$$G(\zeta, \Delta z) = \exp\left[ -D \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q \right]$$

In this framework, **Machine Learning (ML)** analysis using Physics-Informed Neural Networks (PINNs) reveals that grid-scale sensitivity ($\Delta z$) accounts for approximately **38% of the relative importance** in determining the correction factor, second only to the stability parameter $\zeta$ (45%).

### **4. Impact on Mixing Tails**
Coarse grid spacing necessitates using **"longer-tailed"** stability functions, which allow turbulent exchange to persist at higher stability values to compensate for the loss of resolved structure. The correction factor $G$ effectively transforms a fundamental "short-tail" function into a resolution-dependent form that mimics these longer-tailed behaviors on coarse grids while maintaining the correct physical limits as resolution improves.

***

### **Analogy: The Winding Mountain Road**

The influence of grid spacing is like **measuring the length of a winding mountain road on a map**. 

If you use a **low-resolution map** with measurement points 10 miles apart (a coarse grid), you will draw straight lines between points that ignore every bend and switchback. This makes the road seem shorter and flatter than it truly is—analogous to **underestimating stability** in the atmosphere.

The **stability correction factor** acts like a GPS navigation system that knows the exact curvature of every turn. The larger the distance between your map points ($\Delta z$), the more the GPS must **correct your estimated travel time** (mixing rate) to ensure you don't naively assume a straight path when the road actually spirals upward through hairpin turns.

Just as GPS compensates for map coarseness by adding time for curves you can't see, the correction factor $G$ compensates for grid coarseness by reducing mixing coefficients to prevent the model from "taking shortcuts" through the stable layer that don't exist in reality.