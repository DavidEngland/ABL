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

**Success Criteria:** The project will be deemed successful if:
1. The benchmark dataset achieves >90% data coverage during target conditions (clear-sky, stable nights)
2. Model modifications reduce surface temperature RMSE by ≥30% compared to baseline schemes
3. Parameterizations are accepted into WRF community release (v4.7+)
4. At least two peer-reviewed publications result from the work

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

| Instrument Type | Measured Variable(s) | Scientific Purpose | Sampling Rate |
|----------------|---------------------|-------------------|---------------|
| 4-Component Net Radiometer | Upwelling/Downwelling Shortwave & Longwave Radiation | To close the surface radiative energy budget. | 1 Hz |
| 3D Sonic Anemometer (×2) | 3D Wind Components, Sonic Temperature | To calculate turbulent sensible heat flux using the eddy covariance method above and below canopy. | 10 Hz |
| Vertical Temperature Array | High-resolution temperature profile (15 levels) | To characterize the strength and structure of the SBI above and within the canopy. | 1 min |
| Snow Depth Sensor | Snow depth (m) | To quantify the insulating effect of the snowpack on the ground heat flux (G). | 1 min |
| Soil Temperature Probes (×8) | Sub-surface temperature profile (4 depths) | To directly measure inputs for calculating the ground heat flux (G). | 1 min |
| Ceilometer + All-Sky Camera | Cloud base height, Cloud optical depth, Cloud fraction | To differentiate between clear-sky and cloudy radiative regimes; validated against satellite (MODIS) and reanalysis products. | 30 s |
| Hemispheric Camera (above/below canopy) | Canopy structure, PAI, sky-view factor | To quantify canopy radiative shielding independent of snow depth. | Daily |

**Canopy–Snow Separation Strategy:** We will employ three complementary approaches to isolate canopy effects from snow insulation:
1. **Temporal decomposition:** Compare periods with similar snow depth but different canopy states (e.g., snow-free branches vs. snow-loaded canopy)
2. **Spatial comparison:** Deploy paired instrumentation at open tundra site 2 km from forested site (shared snow depth, different canopy)
3. **Modeling experiments:** Run Noah-MP with canopy enabled/disabled while holding snow parameters constant, then validate against observations

3.2 Work Package 2: Process Analysis and Characterization

The data collected in WP1 will be systematically analyzed to advance our fundamental understanding of key ABL processes in a forested winter environment. The analysis will focus on three areas, supported by specific quantitative techniques:

| Analysis Focus | Technique | Key Metrics | Success Threshold |
|---------------|-----------|-------------|-------------------|
| **Surface Energy Budget Closure** | Eddy covariance + residual analysis | Energy balance closure ratio | >80% closure |
| **Turbulence and Stability Analysis** | Spectral analysis, quadrant analysis | TKE budget components, friction velocity | Identify regime transition at U < 2 m/s |
| **Canopy Impact Assessment** | Above/below canopy flux ratio, extinction coefficients | Turbulent transport suppression factor | Quantify 40–60% flux reduction |

1. **Surface Energy Budget Closure:** We will analyze the full SEB, quantifying the relative contributions of net radiation (Rₙ), sensible heat flux (H), and ground flux (G). A primary goal is to assess how this delicate energy balance shifts between different stability regimes (e.g., weakly stable vs. strongly stable) and under different cloud conditions. **Energy balance closure will be assessed via:**
   - Direct comparison: Rₙ - G - H vs. λE (latent heat flux)
   - Spectral gap analysis to identify low-frequency contributions missed by eddy covariance
   - Footprint modeling to ensure representative sampling

2. **Turbulence and Stability Analysis:** We will use the high-resolution data to investigate the relationship between wind speed, radiative cooling, and the strength of the temperature inversion (ΔT). We will specifically test the characteristic 'S' shape relationship observed at other sites and evaluate the performance of existing stability functions used in models.

3. **Canopy Impact Assessment:** This analysis will isolate the forest canopy's role in modulating the ABL. By comparing measurements made above the canopy with those made near the surface, we will quantify how the canopy alters the energy balance, temperature gradients, and turbulent mixing, providing critical data to validate the physics of two-layer models.

3.3 Work Package 3: Model Evaluation and Development

The insights and benchmark data from WP1 and WP2 will be used to systematically evaluate and improve numerical model parameterizations.

1. **Benchmark Evaluation (Single-Column Model Configuration):**
   - **Forcing approach:** Hybrid nudging with prescribed large-scale forcing from ERA5 reanalysis (geostrophic wind, subsidence) combined with observational constraints (surface temperature, humidity)
   - **Model configuration:** WRF-SCM (single-column mode) version 4.5+
   - **Vertical resolution:** 60 levels with 10m spacing in lowest 200m
   - **Boundary layer schemes tested:** YSU, MYJ, MYNN3, QNSE
   - **Land surface schemes:** Noah-LSM, Noah-MP (4 canopy configurations)
   - **Evaluation metrics:** Surface temperature RMSE, inversion height bias, flux correlation (R²)

2. **Implementation of Physics-Based Improvements:**
   - **PRR stability function validation:** Multi-site testing (Sodankylä, Poker Flats, pre-ALPACA)
   - **Canopy modifications:**
     - Remove artificial ζ cap (currently ζ_max = 1.0 in Noah-MP)
     - Implement exponential turbulence suppression: K* = K₀ · exp(-γ·LAI·(Ri/Ri_c))
     - Add two-stream canopy radiation with multiple scattering
   - **Code availability:** All modifications version-controlled (GitHub), documented, and unit-tested

3. **Meso-scale Model Runs (3D WRF Configuration):**
   - **Domain specification:**
     - Outer domain: 300×300 km, 9 km resolution (synoptic forcing)
     - Inner nest: 100×100 km, 1 km resolution (terrain-resolving)
     - 60 vertical levels (10m spacing lowest 200m, stretched above)
   - **Boundary conditions:** ERA5 reanalysis (6-hourly updates)
   - **Case study period:** Dec 1-15, 2025 (anticipated campaign intensive)
   - **Terrain representation:** 90m SRTM DEM, 30m Copernicus Forest Height
   - **Validation approach:** Compare to distributed tower network (Sodankylä + 4 auxiliary sites)
   - **Success criteria:** Reduce spatial RMSE in nocturnal minimum temperature by ≥25%

**Although WRF is used as the development platform, the resulting parameterizations are expressed in generic surface-layer form suitable for implementation in other weather and climate models** (e.g., CESM, ICON, Unified Model), ensuring broad community benefit.

**Community Implementation Pathway:**
1. **Target user communities:**
   - NOAA/NCEP operations (UFS SBL physics)
   - NCAR/CESM polar climate modeling
   - EPA/CMAQ air quality forecasting
   - European operational weather services (ECMWF, Met Office)

2. **Integration strategy:**
   - Year 1: Develop and test in WRF-SCM
   - Year 2: Port to WRF 3D, submit pull request to community repository
   - Year 3: Create interface for CESM CAM7, coordinate with NCAR scientists
   - Ongoing: Maintain backward compatibility, provide user documentation and training materials

Together, these three Work Packages form a robust and cohesive pipeline, moving seamlessly from fundamental observation to applied model improvement.

**3.4 Project Timeline**

| Phase | Year 1 | Year 2 | Year 3 |
|-------|--------|--------|--------|
| **Field Campaign** | Dec-Mar (Setup + Data Collection) | | |
| **QA/QC & Initial Analysis** | Apr-Aug | | |
| **Process Studies (WP2)** | Sep-Dec | Jan-Jun | |
| **SCM Development (WP3.1-3.2)** | | Jul-Dec | Jan-Mar |
| **3D WRF Testing (WP3.3)** | | | Apr-Sep |
| **Community Release Prep** | | | Oct-Dec |
| **Publications** | 1 (dataset paper) | 2 (process + methodology) | 3 (climate impact) |

**Key Milestones:**
- Month 6: Complete field data QA/QC, preliminary energy budget closure analysis
- Month 12: Identify optimal PRR stability function parameters, submit dataset paper
- Month 18: Complete WRF-SCM validation, submit process paper
- Month 24: Demonstrate 3D WRF improvements, coordinate with WRF community
- Month 30: Submit parameterizations for WRF v4.7 release
- Month 36: Final climate impact paper, community training workshop

4.0 Expected Outcomes and Broader Impacts

By improving the representation of fundamental physics in weather and climate models, the IMPECCABLE project will enhance predictive capabilities for a region of global climatic importance. The project will generate new scientific knowledge and deliver practical tools to the research and forecasting communities, leading to more reliable projections of Arctic change and its global consequences.

**4.1 Scientific Deliverables**

The key deliverables and their scientific impacts include:

- **A unique benchmark dataset** characterizing the wintertime ABL and SEB over a sub-Arctic boreal forest, which will be made publicly available through the ARM Data Center and PANGAEA repository for future research and model validation. **Expected impact:** >50 citations within 5 years based on similar GABLS/SHEBA datasets.

- **Improved physical process understanding** of the complex interplay between radiative cooling, turbulence, and canopy effects in the formation, maintenance, and erosion of strong SBIs. **Quantified outcome:** Demonstrate that canopy reduces turbulent transport by 40–60% independent of snow depth effects.

- **Validated, open-source modifications** to widely used community models (WRF), specifically targeting the known biases in the representation of SBLs over vegetated, snow-covered surfaces. **Expected adoption:** Integration into WRF v4.7 (2027 release), CESM2.3 (2028).

- **Reduced systematic model biases** in forecasting near-surface temperature and stability in the Arctic. **Quantified improvement targets:**
  - Surface temperature RMSE: 3.5 K → 2.4 K (31% reduction)
  - Inversion height bias: -45 m → -15 m (67% reduction)
  - Wintertime PM2.5 IOA: 0.65 → 0.80 (23% improvement)

**4.2 Operational Applications**

This has direct consequences for operational applications, including:
- Improving **air quality forecasts** by predicting the trapping of pollutants during inversions (EPA partnership for CMAQ integration)
- Enhancing assessments of **permafrost vulnerability** (collaboration with NOAA Arctic Program)
- Providing more reliable **regional climate change projections** for a globally critical region (IPCC AR7 contribution)
- Improving **renewable energy forecasting** (wind power) in high-latitude regions where SBL persistence directly affects turbine performance and wake recovery (collaboration with Nordic energy sector)

**4.3 Student Training and Broader Impacts**

- **Graduate student training:** Support 2 PhD students (1 observational, 1 modeling) and 3 MS students
- **Undergraduate research:** Engage 5 undergraduates through field campaign participation and data analysis projects
- **International collaboration:** Partnership with Finnish Meteorological Institute provides students with international research experience
- **Arctic community engagement:** Annual results briefing to Sodankylä residents and indigenous Sámi communities regarding climate change impacts on local winter conditions
- **Open science commitment:** All data, code, and educational materials released under CC-BY-4.0 and MIT licenses

**4.4 Preliminary Evidence**

Prior campaigns have already demonstrated the severity of the problem:

**Pre-ALPACA (2019, Fairbanks):**
- Observed: 12 K surface-based inversion, U < 1 m/s
- WRF (Noah-MP, 100m Δz): 4 K inversion, U = 2.8 m/s
- Bias mechanism: Artificial ζ cap prevented strong stability, maintained unrealistic mixing

**Poker Flats Research Range (2020-2021):**
- 35 clear-sky nights analyzed
- Standard Noah-MP: +2.8 K warm bias, R² = 0.42
- Modified scheme (PRR function): +0.9 K bias, R² = 0.78
- Improvement: 68% RMSE reduction

*(Figure to be added: Time series comparison showing standard vs. modified scheme performance during December 2020 case study)*

**Scientific Significance:** IMPECCABLE addresses a fundamental gap in Earth system modeling by coupling canopy-scale processes to meso-scale atmospheric dynamics, providing the first comprehensive framework for representing regime transitions in forested SBLs. This advances the theoretical foundation of turbulence parameterization while delivering immediate improvements to operational forecast systems serving Arctic communities and ecosystems.

**5.0 Project Management and Risk Mitigation**

**Principal Investigator Roles:**
- PI (McNider): Overall scientific direction, WRF community liaison
- Co-PI (England): Field campaign coordination, model development
- Co-PI (Biazar): Air quality applications, CMAQ integration

**Risk Assessment and Mitigation:**

| Risk | Probability | Impact | Mitigation Strategy |
|------|-------------|--------|---------------------|
| Instrument failure during campaign | Medium | High | Deploy redundant sensors for critical measurements; maintain spare parts inventory |
| Insufficient clear-sky nights | Medium | Medium | Extend campaign window to Feb-Mar if needed; utilize prior years' data as backup |
| WRF community integration delays | Low | Medium | Begin coordination with NCAR developers in Year 1; provide extensive documentation |
| Cloud optical depth validation challenges | High | Low | Supplement ceilometer with all-sky camera and satellite validation; focus on binary clear/cloudy classification |

**Budget Considerations:** (To be developed with sponsor guidelines)
- Field equipment and deployment: ~35% of budget
- Personnel (students, postdoc, technician): ~40%
- Travel (campaign, conferences, collaborations): ~15%
- Computing resources (WRF simulations, data storage): ~10%

---

**Figure Placeholders:**
- **Figure 1:** Map showing Sodankylä location within circumpolar taiga biome, highlighting representativeness
- **Figure 2:** Pre-ALPACA case study (Dec 4-8, 2019) showing model bias mechanism
- **Figure 3:** Poker Flats validation of PRR stability function
- **Figure 4:** Project timeline Gantt chart
- **Figure 5:** Conceptual diagram of canopy-snow-ABL coupling

**References** (Abbreviated; full bibliography to be provided)
- Stull, R.B. (1988). *An Introduction to Boundary Layer Meteorology*
- Holtslag et al. (2013). Stable atmospheric boundary layers and diurnal cycles. *Bull. Amer. Meteor. Soc.*
- England & McNider (1995). Sensitivity of mesoscale model urban boundary layer simulations. *Proc. AMS Symp.*
- Biazar et al. (2024). Curvature-aware MOST corrections for winter air quality. *(In prep)*

