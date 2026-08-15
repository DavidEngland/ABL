Grid-Curvature Correction for Eddy Diffusivity in the Stable Boundary Layer: A Unified McNider-Biazar Approach
Author: Your Name
Affiliation: Department of Atmospheric Sciences, University
Date: December 22, 2025
Abstract
Accurately representing turbulent exchange in the Stable Boundary Layer (SBL) remains a central challenge in modeling, often leading to systematic warm biases on coarse grids. This study addresses the grid–curvature bias where the bulk Richardson number ($Ri_b$) systematically underestimates the true gradient Richardson number ($Ri_g$) due to the non-linear curvature of Monin–Obukhov similarity functions. We develop a unified framework combining the McNider grid-invariance principle with techniques inspired by the Biazar methodology (Adomian Decomposition Method) for analyzing base diffusivities. By explicitly diagnosing the curvature distortion via the ratio $B = Ri_g / Ri_b$, we derive a correction factor, $f_c$, that scales the eddy diffusivity $K$. This factor ensures the layer-integrated turbulent transport remains invariant under changes in vertical resolution, providing a mathematically transparent and physically grounded pathway toward resolution-aware turbulence closures.
1. Introduction
Accurately representing turbulent exchange in the Stable Boundary Layer (SBL) remains one of the central challenges in contemporary atmospheric modeling (Holtslag et al., 2013; Mahrt, 2014). The combination of weak turbulence, strong stratification, and intermittent mixing produces a system highly sensitive to vertical resolution. Numerical weather prediction (NWP) and climate models continue to exhibit well-documented warm near-surface biases in stably stratified regimes (Davy & Esau, 2014; Tjernström et al., 2005). Even when underlying physical schemes are nominally well-behaved, structural numerical distortions persist.
A key contributor to this distortion is the grid–curvature bias: the mismatch between the nonlinear curvature of Monin–Obukhov Similarity Theory (MOST) (Businger et al., 1971; Monin & Obukhov, 1954) and its piecewise-linear representation on discrete grids. When vertical spacing exceeds the intrinsic curvature scale of stability functions, numerical differencing underestimates $Ri_b$ relative to the true $Ri_g$. Because turbulent closures, such as those by Beljaars and Holtslag (1991), construct eddy diffusivities as steeply decreasing functions of stability, even modest underestimations of $Ri_b$ lead to excessive mixing.
This study builds a unified framework combining the McNider grid-invariance principle (England & McNider, 1995), which requires that layer-integrated transport be resolution-independent, and nonlinear decomposition techniques developed by Biazar and collaborators (England et al., 2025), which provide analytical means for solving the nonlinear closure equations that determine base diffusivities.
2. Curvature and the Richardson Number Bias
The gradient Richardson number is traditionally defined through MOST functions as:
$$Ri_g = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}, \quad \zeta = z/L$$
Empirically validated stability functions (Businger et al., 1971; Dyer, 1974; Högström, 1988) show that $Ri_g(\zeta)$ is strictly concave down over most of the stable regime. By Jensen's inequality, the layer-average $Ri_b$ computed over a grid cell $\Delta z$ underestimates the pointwise value:
$$Ri_b = \frac{1}{\Delta z} \int_{z_0}^{z_1} Ri_g(z) dz < Ri_g(z_*)$$
The magnitude of this distortion is captured by the curvature ratio $B(\zeta) = Ri_g(\zeta) / Ri_b(\zeta)$. Because eddy diffusivity $K$ often scales as $Ri_g^{-n}$ (Louis, 1979), these numerical artifacts are frequently misinterpreted as physical failures of the turbulence scheme.
3. The Curvature Correction Framework
Following the grid-invariance principle (England & McNider, 1995), we require that the product of the stability function and the correction factor remains invariant with respect to grid spacing. Using $B$ as the diagnostic, we derive a closed-form power-law solution:
$$f_c(\Delta z, \zeta) = \left( \frac{\Delta z}{\Delta z_{ref}} \right)^{-\alpha(B(\zeta)-1)(\zeta/\zeta_{ref})^q}$$
The corrected diffusivity $K_{new} = f_c K_{old}$ ensures that for coarser grids, diffusivity is reduced in proportion to the curvature-induced bias.
4. Biazar-Style Decomposition and Machine Learning Integration
The curvature correction depends on the base diffusivities $K_{old}$, often determined by nonlinear TKE equations. Adomian Decomposition Methods (ADM) (England et al., 2025) represent the solution as a rapidly convergent series, avoiding the sensitivities of iterative solvers in strongly stable regimes (Gryanik & Lüpkes, 2020).
Recent advances in machine learning have further expanded the potential for SBL parameterization. Studies have utilized random forests for surface roughness estimation (Duan & Takemi, 2021) and ensemble models for marine MABL prediction (Chen et al., 2025). Physics-Informed Neural Networks (PINNs) offer a way to bridge these data-driven approaches with the physical laws described by the McNider framework (Alamu et al., 2025).
5. Conclusion
By diagnosing curvature through the ratio $B$ and enforcing grid invariance, we provide a resolution-aware turbulence closure. This synthesis addresses long-standing deficiencies in MOST-based closures and provides a mathematically transparent pathway for modeling stable regimes from LES to global scales.
Works Cited
Alamu, R., et al. (2025). Physics-Informed Neural Networks for Climate Modeling. ResearchGate / Preprint.
Beljaars, A. C. M., & Holtslag, A. A. M. (1991). Flux Parameterization Over Land Surfaces for Atmospheric Models. Journal of Applied Meteorology, 30(4), 327-341.
Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971). Flux-Profile Relationships in the Atmospheric Surface Layer. Journal of the Atmospheric Sciences, 28(2), 181-189.
Chen, Y., et al. (2025). A Comprehensive Ensemble Model for Marine Atmospheric Boundary-Layer Prediction. Remote Sensing, 17.
Davy, R., & Esau, I. (2014). Global Climate Models' Bias in Surface Temperature Trends and Variability. Environmental Research Letters, 9(11), 114024.
Duan, G., & Takemi, T. (2021). Predicting Urban Surface Roughness Aerodynamic Parameters Using Random Forest. Journal of Applied Meteorology and Climatology, 60(7), 999-1018.
Dyer, A. J. (1974). A Review of Flux-Profile Relationships. Boundary-Layer Meteorology, 7(3), 363-372.
England, D. E., & McNider, R. T. (1995). Stability Functions Based upon Shear Functions. Boundary-Layer Meteorology, 74(1-2), 113-130.
England, D. E., McNider, R. T., & Biazar, A. P. (2025). Curvature-Aware MOST Toolkit for Stable Boundary Layer Corrections. GitHub Repository.
Gryanik, V. M., & Lüpkes, C. (2020). New Modified and Extended Stability Functions for the Stable Boundary Layer Based on SHEBA and GABLS1 Data. Journal of the Atmospheric Sciences, 77(8), 2687-2716.
Högström, U. (1988). Non-Dimensional Wind and Temperature Profiles in the Atmospheric Surface Layer: A Re-Evaluation. Boundary-Layer Meteorology, 42(1-2), 55-78.
Holtslag, A. A. M., et al. (2013). Stable Atmospheric Boundary Layers and Diurnal Cycles: Challenges for Weather and Climate Models. Bulletin of the American Meteorological Society, 94(11), 1691-1706.
Louis, J. F. (1979). A Parametric Model of Vertical Eddy Fluxes in the Atmosphere. Boundary-Layer Meteorology, 17(2), 187-202.
Mahrt, L. (2014). Stably Stratified Atmospheric Boundary Layers. Annual Review of Fluid Mechanics, 46, 23-45.
Monin, A. S., & Obukhov, A. M. (1954). Basic Laws of Turbulent Mixing in the Ground Layer of the Atmosphere. Trudy Geofiz. Inst. Akad. Nauk SSSR, 151, 163-187.
Tjernström, M., et al. (2005). Modelling the Arctic Boundary Layer: An Evaluation of Six ARCMIP Regional Scale Models. Boundary-Layer Meteorology, 117(2), 337-381.
