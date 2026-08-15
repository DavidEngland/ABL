The **Random Forest model (RF_G)**, used to predict the eddy diffusivity correction factor $G$, is primarily driven by two key variables that together account for **85% of its predictive power**:

*   **$\zeta_{\text{geom}}$ (Stability):** This represents the dimensionless height at the geometric mean, which captures the stability of the atmospheric layer. In the sources, $\zeta$ is defined as $z/L$, where $L$ is the Obukhov length.
*   **$\Delta z$ (Grid Scale):** This is the vertical grid spacing or layer thickness.

The model uses these variables to learn the core physical hypothesis regarding how coarse grid scales distort the representation of stability in the Stable Boundary Layer. By using these inputs, the Random Forest acts as a high-performance, near-perfect surrogate ($R^2 = 0.9984$) for the complex analytical integrals of **Monin-Obukhov Similarity Theory (MOST)**.

This approach allows for the instant calculation of corrected eddy diffusivities ($K_{\text{new}} = K_{\text{old}} \cdot G_{\text{pred}}$) without the computational latency required by traditional numerical integration.

**Analogy**
Think of the Random Forest model as a **highly experienced navigator** who only needs two pieces of information to tell you how much to adjust your speed: how steep the hill is (**$\zeta$**) and how far apart the markers on your map are (**$\Delta z$**). Instead of doing complex calculus on the fly, the navigator has "memorized" the corrections needed for every combination of slope and map resolution.