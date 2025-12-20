The systematic underestimation of stability by **bulk Richardson numbers ($Ri_b$)** in Numerical Weather Prediction (NWP) models is a structural bias caused by the interaction between **coarse vertical resolution** and the **nonlinear curvature** of atmospheric stability functions.

This phenomenon is primarily explained by the following factors:

### 1. Local vs. Layer-Averaged Stability
In NWP models, the atmosphere is discretized into vertical grid layers ($\Delta z$), often 50–100 meters thick. Models calculate a **bulk Richardson number ($Ri_b$)**, which is a layer-averaged diagnostic of stability. However, the physical processes that govern the collapse or persistence of turbulence are driven by the **local gradient Richardson number ($Ri_g$)**.

### 2. Concave-Down Curvature ($\Delta < 0$)
In the Stable Boundary Layer (SBL), the profile of local stability ($Ri_g$) is typically **concave-down**. This means that stability increases very rapidly near the surface and then begins to level off. This curvature is quantified by the **neutral curvature invariant ($\Delta$)**; when $\Delta < 0$, the $Ri_g$ profile "bends" below a linear path as it departs from neutrality.

### 3. Jensen’s Inequality
The mathematical core of the bias is **Jensen’s Inequality**, which states that for any concave-down function, the average value over an interval is **mathematically guaranteed to be less than** the function's value at its representative center. Specifically, models find that:
$$Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_g)$$
where **$z_g$ is the geometric mean height** ($\sqrt{z_0 z_1}$), the natural midpoint for logarithmic atmospheric profiles.

### 4. The Grid Sensitivity Crisis
The discrepancy between the model's bulk average ($Ri_b$) and the true local stability at the center ($Ri_g$) is highly sensitive to grid resolution. On coarse grids, the model averages over a much larger "gap" in the curve, causing the **Bias Ratio ($B = Ri_g/Ri_b$)** to increase significantly. Numerical evaluations show that as stability increases, this bias ratio can nearly double, reaching values where the model identifies only half the true local stability.

### 5. Consequences: The Over-mixing Problem
Because the model's $Ri_b$ is systematically too low, the NWP closure scheme perceives the atmosphere as being less stable than it actually is. This leads to:
*   **Excessive Diffusion:** The model calculates eddy diffusivities ($K_m, K_h$) that are too high.
*   **Spurious Mixing:** The model "over-mixes" heat and pollutants.
*   **Warm Biases:** Over-mixing results in erroneous surface temperatures, such as the 2–3°C warm biases frequently observed in Arctic climate models during the polar night.

**Analogy**
Using a bulk Richardson number on a coarse grid is like trying to measure the slope of a **steeply rounded dome** by laying a long, straight ladder against it. Because the dome is curved (concave-down), the flat ladder will always sit much lower than the actual peak of the curve. A model relying on that ladder would mistakenly conclude the terrain is much flatter than it really is, leading it to predict a smooth walk where there is actually a steep climb.