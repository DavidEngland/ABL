 **Triple-Point Invariant-Convergence Hypothesis**. This framework is used to determine whether an observed collapse in a vertical profile reflects a genuine physical transition or is merely an artifact of the numerical grid or sensor spacing.

The triad consists of three independent physical transition trackers that are monitored for vertical convergence as the measurement or model resolution increases:

### 1. The Triad of Heights
*   **Diffusivity Threshold Height $z_K$:** The altitude where the eddy diffusivity for momentum reaches its minimum value or floor $K_m \to K_{\text{min}}$.
*   **TKE Floor Height $z_e$:** The altitude where Turbulent Kinetic Energy (TKE) reaches its minimum threshold $e \to e_{\text{min}}$.
*   **TKE Gradient Height $z_{e_z}$:** The altitude where the vertical gradient of TKE $\partial e / \partial z$ reaches a localized extremum.

### 2. The Convergence Hypothesis
The hypothesis posits that while these heights may appear separated on a coarse measurement grid, they will **collapse onto a single invariant altitude** during a turbulence extinction event, regardless of instrument separation $\Delta z$ or grid spacing.

The mathematical formulation for this convergence is expressed as:
$$C_{\text{TP}} = \frac{\Delta z_{\text{TP}}}{H_{\text{SBL}}} \to 0 \quad \text{as } \Delta z \to 0$$
where:
*   $\Delta z_{\text{TP}}$ is the absolute spatial spread between the three heights $\max_{i,j} |z_i - z_j|$.
*   $H_{\text{SBL}}$ is the total boundary layer height.

### 3. Purpose and Significance
*   **Bypassing Grid Artifacts:** Coarse resolution $\Delta z \gg L$ acts as a spatial low-pass filter that can shift perceived turbulence thresholds. Tracking this triad ensures that a perceived "fold" in the profile is an invariant physical transition.
*   **Scientific Evidence Hierarchy:** Resolving these discretization errors and demonstrating numerical convergence is classified as a **Level 3 prerequisite**. No quantitative physical parameter inference (such as estimating mixing lengths) can be reliably drawn from a model until this numerical convergence is demonstrated.