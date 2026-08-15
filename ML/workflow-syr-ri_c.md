This is a detailed, well-structured workflow for Symbolic Regression of the Dynamic Critical Richardson Number (Ri_c^*). The process logically covers defining objectives, data preparation, modeling, validation, and deployment.
Given the highly technical nature of the problem, a diagram illustrating the Atmospheric Boundary Layer (ABL) structure and the key Richardson Number components would greatly aid understanding of the features and targets.
🌬️ Dynamic Critical Richardson Number (Ri_c^*) Workflow
The goal is to discover a compact, interpretable formula Ri_c^* = f(\text{state}) that improves the timing of turbulence onset/cessation compared to a fixed critical value (e.g., Ri_c = 0.25).
1. Key Physical Concepts & Features
The symbolic regression model aims to represent the physics governing the transition between turbulent and laminar flow, primarily controlled by the balance of buoyancy (stability) and shear.
| Concept | Feature/Variable | Description |
|---|---|---|
| Local Stability | \zeta = z/L | The dimensionless height in Monin-Obukhov Similarity Theory (MOST). L is the Obukhov length. Preserves the sign of stability. |
| Bulk Stability | Ri_b | The Bulk Richardson Number, the ratio of buoyancy term to shear term, typically calculated across a layer \Delta z. |
| Buoyancy/Shear Terms | \Gamma = d\theta/dz, $S = | dU/dz |
| Memory/Context | \text{previous\_turb\_flag} | Captures hysteresis/lag effects, acknowledging that a turbulent state is harder to shut down than a laminar state is to initiate. |
2. Symbolic Regression Setup Highlights
The setup leverages PySR (or similar) but is carefully constrained to ensure the output is physically meaningful and practically useful.
 * Formula Domain: Ri_c^* \in [0.15, 2.0]. (Constraint enforced via clipping).
 * Loss Function: A composite loss is critical for balancing regression accuracy with operational performance:

 * Physical Constraints: Monotonicity constraints must be enforced via penalty terms on the loss function, for example, penalizing \partial Ri_c^* / \partial \Gamma < 0 (as increasing stability \Gamma should increase the required critical Richardson number).
3. Training & Validation Protocol
The protocol is designed for robustness and avoiding site-specific overfitting.
| Stage | Focus | Key Action |
|---|---|---|
| Data Split | Generalization | Site/Experiment Holdout for validation and testing to ensure the formula works across different meteorologies and instrumentation. |
| Stage 1 (Search) | Exploration | Broad symbolic search using PySR on the train set with relaxed complexity. Generate a diverse Pareto front. |
| Stage 2 (Refine) | Accuracy | Refine top N candidates by performing local optimization of coefficients (bounded least squares) to minimize extrapolation errors. |
| Stage 3 (Test) | Final Check | Evaluate the final formula on unseen test sites and perform the critical closed-loop column-model test. |
4. Closed-Loop Validation (Operational Test)
This is the most critical validation step, as it verifies the formula's utility in a real-world modeling context.
 * Integration: Replace the fixed Ri_c threshold in the vertical diffusion closure of a single-column model with the dynamic Ri_c^*(\text{state}).
 * Success Metrics: Assess the change in surface bias (e.g., T_{\text{surf}} error), energy consistency, and the frequency of false collapses (where turbulence is wrongly shut off).
The ultimate deliverable is a set of compact, interpretable formulas (Top-3) along with the code and documentation for immediate integration and testing.
