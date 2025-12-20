Research Proposal: The IMPECCABLE Project (Improving Model Parameterizations of the Continental Arctic Boundary Layer Ecosystem)

1.0 Scientific Rationale: The Arctic Modeling Imperative

The Arctic is warming at a rate two to three times that of the rest of the globe. This phenomenon, known as "Arctic Amplification," is driven by a complex series of feedback mechanisms that are initiated at the surface and unfold within the lowest part of the atmosphere, the atmospheric boundary layer (ABL). Processes controlling the surface energy balance (SEB), turbulence, and the formation of clouds and surface-based temperature inversions (SBIs) are fundamental to this warming, yet they remain poorly represented in the current generation of weather and climate models. This leads to persistent biases that undermine our ability to forecast Arctic weather, project regional climate change, and predict associated impacts like air quality. This proposal outlines the IMPECCABLE project, a targeted research initiative designed to address these critical model deficiencies by advancing our fundamental understanding of the physical processes governing the continental Arctic ABL.

1.1 The Challenge of Stable Boundary Layers (SBLs)

A defining characteristic of the high-latitude atmosphere, particularly during the long winter, is the frequent formation of strong SBIs. These stable boundary layers (SBLs) trap cold air near the ground and decouple the surface from the warmer atmosphere above. Numerical models, however, exhibit a persistent warm bias in these conditions, largely because they generate excessive vertical mixing that erroneously erodes the inversion. This issue stems from what can be called a "Grid Sensitivity Crisis": the coarse vertical resolution common in models is insufficient to resolve the sharp temperature and wind gradients near the surface.

**Unlike prior SBL studies focused on tundra or homogeneous surfaces, IMPECCABLE explicitly targets the coupled canopy–snow–ABL system and its role in triggering regime transitions between weakly and strongly stable states.** Boreal forest canopies fundamentally alter the radiation budget, insulate the snowpack, and modify turbulence structure in ways that existing single-layer schemes cannot represent. This canopy–atmosphere coupling is a first-order control on Arctic Amplification yet remains largely unparameterized in operational models.

This coarse resolution leads to a systematic underestimation of atmospheric stability. This is because models fail to resolve the vertical 'curvature' of the stability profile. As established by Monin-Obukhov Similarity Theory (MOST), stability departs from neutrality at a specific rate near the surface (a property defined by the neutral curvature invariant, 2Δ). By averaging over coarse vertical grid cells, models violate this physical constraint, leading to a bulk Richardson number that is systematically underestimated (an effect of Jensen's inequality on a concave function). This miscalculation causes the model's physics schemes to predict far more turbulence than exists in reality, which in turn transfers too much heat down to the surface, weakening the inversion and producing the characteristic warm bias in forecasts.

1.2 The Unresolved Impact of Local-Scale Heterogeneity

The challenge of modeling the Arctic ABL is compounded by local-scale phenomena that are too small to be resolved by regional or global models but have a dominant effect on local conditions. The 2019 pre-ALPACA field campaign in Fairbanks, Alaska, provides a stark example. During a period of clear skies and strong radiative cooling, a strong SBI was expected to form across the region, and indeed, one did develop at the local airport. However, at our measurement site just 3.5 km away, a completely different regime was observed.

A topographically-driven drainage flow, originating from the nearby Goldstream valley, channeled a persistent flow of air across the site. This local wind maintained a high level of mechanical turbulence, which continuously mixed the near-surface atmosphere. Consequently, instead of a strong, radiatively-controlled inversion, the site remained in a weakly stable state. This single local flow completely altered the SEB and stability regime. This demonstrates a critical breakdown in model fidelity, where unresolved topographic flows can create a completely different boundary layer regime just a few kilometers away, rendering regional forecasts fundamentally unreliable.

1.3 Compounding Uncertainties: Clouds and Forest Canopies

Beyond stability and local flows, models must also contend with the complex influences of clouds and vegetation. Arctic clouds are extremely frequent, with coverage exceeding 85% from May to October over sea ice. Their impact on the surface radiation budget is bimodal and profound, creating two distinct atmospheric states: a "radiatively clear" state with strong surface cooling (net longwave flux of approximately -40 W/m² during winter over sea ice) and an "opaquely cloudy" state where the cloud's thermal emission warms the surface, bringing the net flux near zero.

The vast sub-Arctic boreal forest introduces another layer of unresolved physics. A forest canopy dramatically alters the SEB by interacting with radiation, insulating the snow-covered ground, and modifying turbulence. Modeling studies using the Weather Research and Forecasting (WRF) model show that standard single-layer schemes fail to capture the near-surface temperature structure over forested terrain. More advanced two-layer schemes, such as Noah with Multi-Parameterization options (Noah-MP), perform better but still suffer from critical flaws. These schemes often impose artificial, non-physical limits on turbulence in very stable conditions—for example, by capping the stability parameter ζ—which prevents the model from simulating the very strong inversions that are routinely observed at forested sites like the Poker Flats Research Range in Alaska.

The IMPECCABLE project is designed to directly address these interconnected modeling deficiencies through a targeted research campaign and a systematic model development effort.

2.0 Project Objectives

The primary goal of the IMPECCABLE project is to resolve the persistent warm bias and excessive vertical mixing that plague model simulations of the sub-Arctic winter boundary layer. We will achieve this by generating a benchmark dataset of canopy-atmosphere interactions and developing improved physical parameterizations that correctly capture the transition between weakly and strongly stable regimes.

**Central Hypothesis:** Explicit representation of canopy-mediated turbulence suppression eliminates the need for artificial stability caps while improving inversion strength and persistence, particularly during radiatively-driven regime transitions between sustainable and collapsing turbulence.

To achieve this overarching goal, we will pursue four specific and measurable objectives:

1. **Quantify Canopy-Atmosphere Interactions:** To perform a comprehensive set of measurements that characterize the complete SEB, turbulent fluxes, and thermodynamic structure of the ABL above and within a boreal forest canopy during winter.

2. **Isolate Key Stability Regimes:** To use the collected data to analyze the transition between weakly and strongly stable regimes in a forested environment, testing the applicability of organizing frameworks such as the "minimum wind speed for sustainable turbulence" (MWST) theory. **MWST provides a physically grounded criterion separating sustainable turbulence from radiatively collapsing flow, making it a natural organizing framework for regime transitions in forested SBLs.**

3. **Evaluate and Refine Model Physics:** To leverage the high-resolution field data to evaluate the performance of existing surface layer and land surface schemes (e.g., Mellor-Yamada-Janjić [MYJ], Noah-MP in WRF) and to test physics-based modifications that better represent observed processes, such as canopy effects and turbulence suppression in very stable conditions.

4. **Develop Improved Parameterizations:** To synthesize the observational findings and modeling results into refined, physics-based parameterizations that can be implemented in meso-scale and single-column models to reduce systematic errors in Arctic forecasts.

These objectives form a cohesive research plan that directly links targeted field observations to the development of improved predictive tools.

3.0 Research Plan and Methodology

The IMPECCABLE research plan is structured into three interconnected Work Packages (WPs) that form a complete "observation-to-model" pipeline. This integrated approach ensures that our advances in process-level understanding are directly translated into tangible improvements in numerical modeling capabilities.

3.1 Work Package 1: The Sodankylä Winter Field Campaign

**Site Justification.** Building upon the insights gained from our campaigns in Alaska (pre-ALPACA and Poker Flats), the logical next step is to acquire a comprehensive dataset in an environment representative of the vast sub-Arctic boreal forest ecosystem. The Sodankylä Arctic Research Centre in Finland offers an ideal location. It is situated within the taiga biome and provides the logistical support necessary for a long-term winter deployment. This campaign will directly extend the foundational process-level understanding gained at Alaskan sites, allowing us to test the generality of our findings in the vast, circumpolar taiga biome.

**Instrumentation and Measurement Strategy.** A suite of high-precision instruments will be deployed to continuously monitor the key variables governing the SEB and ABL structure. The proposed instrumentation is summarized below.

| Instrument Type | Measured Variable(s) | Scientific Purpose |
|----------------|---------------------|-------------------|
| 4-Component Net Radiometer | Upwelling/Downwelling Shortwave & Longwave Radiation | To close the surface radiative energy budget. |
| 3D Sonic Anemometer | 3D Wind Components, Sonic Temperature | To calculate turbulent sensible heat flux using the eddy covariance method. |
| Vertical Temperature Array | High-resolution temperature profile | To characterize the strength and structure of the SBI above and within the canopy. |
| Snow Depth Sensor | Snow depth (m) | To quantify the insulating effect of the snowpack on the ground heat flux (G). |
| Soil Temperature Probes | Sub-surface temperature profile | To directly measure inputs for calculating the ground heat flux (G). |
| Ceilometer or Micro-Lidar | Cloud base height, Cloud optical depth | To differentiate between clear-sky and cloudy radiative regimes. |

3.2 Work Package 2: Process Analysis and Characterization

The data collected in WP1 will be systematically analyzed to advance our fundamental understanding of key ABL processes in a forested winter environment. The analysis will focus on three areas:

1. **Surface Energy Budget Closure:** We will analyze the full SEB, quantifying the relative contributions of net radiation (Rₙ), sensible heat flux (H), and ground flux (G). A primary goal is to assess how this delicate energy balance shifts between different stability regimes (e.g., weakly stable vs. strongly stable) and under different cloud conditions.

2. **Turbulence and Stability Analysis:** We will use the high-resolution data to investigate the relationship between wind speed, radiative cooling, and the strength of the temperature inversion (ΔT). We will specifically test the characteristic 'S' shape relationship observed at other sites and evaluate the performance of existing stability functions used in models.

3. **Canopy Impact Assessment:** This analysis will isolate the forest canopy's role in modulating the ABL. By comparing measurements made above the canopy with those made near the surface, we will quantify how the canopy alters the energy balance, temperature gradients, and turbulent mixing, providing critical data to validate the physics of two-layer models.

3.3 Work Package 3: Model Evaluation and Development

The insights and benchmark data from WP1 and WP2 will be used to systematically evaluate and improve numerical model parameterizations.

1. **Benchmark Evaluation:** The Sodankylä dataset will serve as a high-fidelity benchmark. We will run "offline" single-column versions of WRF's surface and land-surface schemes (e.g., Noah Land Surface Model [Noah-LSM], Noah-MP) forced with observed conditions. The model output, particularly the relationship between the temperature inversion (ΔT) and wind speed, will be compared directly against the field measurements.

2. **Implementation of Physics-Based Improvements:** Based on known deficiencies identified in prior work and confirmed by our analysis, we will implement and test specific model improvements. This will include modifications developed in foundational thesis work, such as removing the artificial cap on the stability parameter (ζ) and implementing a new stability function derived directly from observations (the "Poker Flats Research Range [PRR] stability function"). **The PRR stability function has been shown in prior single-site studies to reproduce observed ΔT–U relationships but has not yet been tested across canopy regimes or implemented in a community model.** IMPECCABLE will provide the first multi-site validation and operational implementation pathway.

3. **Meso-scale Model Runs:** The most successful modified schemes (e.g., the "modified Noah-MP [mMP]" version) will be implemented back into the full 3D WRF model. A case study simulation, analogous to the 4–8 December 2019 period in Alaska, will be conducted over the Sodankylä region. This will allow us to assess whether the improvements demonstrated in the offline model translate to a more realistic spatial representation of SBLs across the complex, forested landscape.

**Although WRF is used as the development platform, the resulting parameterizations are expressed in generic surface-layer form suitable for implementation in other weather and climate models** (e.g., CESM, ICON, Unified Model), ensuring broad community benefit.

Together, these three Work Packages form a robust and cohesive pipeline, moving seamlessly from fundamental observation to applied model improvement.

4.0 Expected Outcomes and Broader Impacts

By improving the representation of fundamental physics in weather and climate models, the IMPECCABLE project will enhance predictive capabilities for a region of global climatic importance. The project will generate new scientific knowledge and deliver practical tools to the research and forecasting communities, leading to more reliable projections of Arctic change and its global consequences.

The key deliverables and their scientific impacts include:

- **A unique benchmark dataset** characterizing the wintertime ABL and SEB over a sub-Arctic boreal forest, which will be made publicly available to the wider scientific community for future research and model validation.

- **Improved physical process understanding** of the complex interplay between radiative cooling, turbulence, and canopy effects in the formation, maintenance, and erosion of strong SBIs.

- **Validated, open-source modifications** to widely used community models (WRF), specifically targeting the known biases in the representation of SBLs over vegetated, snow-covered surfaces.

- **Reduced systematic model biases** in forecasting near-surface temperature and stability in the Arctic. This has direct consequences for operational applications, including:
  - Improving **air quality forecasts** by predicting the trapping of pollutants during inversions
  - Enhancing assessments of **permafrost vulnerability**
  - Providing more reliable **regional climate change projections** for a globally critical region
  - Improving **renewable energy forecasting** (wind power) in high-latitude regions where SBL persistence directly affects turbine performance and wake recovery

**Scientific Significance:** IMPECCABLE addresses a fundamental gap in Earth system modeling by coupling canopy-scale processes to meso-scale atmospheric dynamics, providing the first comprehensive framework for representing regime transitions in forested SBLs. This advances the theoretical foundation of turbulence parameterization while delivering immediate improvements to operational forecast systems serving Arctic communities and ecosystems.

