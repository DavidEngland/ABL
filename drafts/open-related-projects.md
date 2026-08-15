Here is a clean, up-to-date catalog of open-source numerical models relevant to the planetary boundary layer (PBL), atmospheric sciences, Arctic research, and Earth system modeling. I divide them by purpose and note model type, programming language, grid structure, and what part of the PBL physics they handle best.

This list includes all major community models and many specialized ones.

⸻

🚀 1. Large-Eddy Simulation (LES) Models

LES models explicitly resolve boundary-layer turbulence and are widely used for PBL research and training ML closures.

📌 PALM
   •   Language: Fortran
   •   Focus: Urban canopy, Arctic BL, turbulence-resolving LES
   •   Strength: Very active development; GPU support; cloud and radiation modules
   •   License: GPL
   •   URL: “PALM model system”

📌 DALES (Dutch Atmospheric LES)
   •   Language: Fortran
   •   Focus: Cloud-topped BL, idealized convection, stratocumulus
   •   Strength: Clean physics; widely used in PBL closure research
   •   License: GPL

📌 UCLA-LES
   •   Language: Fortran
   •   Focus: Marine stratocumulus, Arctic cloud layers
   •   Strength: Classic model with extensive validation

📌 MicroHH
   •   Language: C++
   •   Focus: High-performance LES; very fast; ideal for turbulence studies
   •   Strength: GPU-enabled; strong for Arctic stable BL

📌 OpenFOAM (atmospheric LES setups)
   •   Language: C++
   •   Focus: General CFD + ABL and wind engineering
   •   Strength: Flexible solver library; PBL via custom boundary conditions
   •   License: GPL

⸻

🌡️ 2. Weather & Regional Climate Models (WRF-class)

These are limited-area NWP/regional climate models with parameterized PBL physics.

📌 WRF (Weather Research and Forecasting Model)
   •   Language: Fortran + C
   •   Strength: Most widely used research NWP model; multiple PBL schemes (MYNN, YSU, Shin-Hong)
   •   License: Open-source (Apache-like)
   •   Notes: Used extensively in Arctic research

📌 MPAS-Atmosphere
   •   Language: Fortran
   •   Grid: Unstructured Voronoi (variable resolution)
   •   Strength: Global but high-resolution Arctic nests without boundaries
   •   License: LGPL

📌 HARMONIE-AROME / HIRLAM community
   •   Language: Fortran
   •   Strength: High-resolution Arctic weather; partially open components
   •   License: Mixed (core partially open)

📌 BRAMS (Brazilian Regional Atmospheric Modeling System)
   •   Based on RAMS; open-source
   •   Advanced PBL and cloud microphysics

⸻

🌍 3. Global Climate / Earth System Models (Open Source)

For full ESM simulations with interactive components (ocean, sea ice, land, chemistry).

📌 CESM (Community Earth System Model)
   •   Language: Fortran
   •   Strength: Full ESM with CAM6 + PBL schemes, active support
   •   License: Open-source
   •   Arctic: Strong sea-ice component (CICE)

📌 E3SM / E3SM-MMF
   •   Mostly open (some components)
   •   Uses turbulence-resolving “super-parameterization” options

📌 ICON-ESM
   •   Language: Fortran
   •   Grid: Icosahedral
   •   License: Partially open (atmosphere open; land/chemistry mixed)

📌 OpenIFS (ECMWF)
   •   Language: Fortran
   •   Strength: Research version of ECMWF model
   •   License: Open but requires agreement

📌 FESOM / AWI-CM
   •   German/Arctic research community
   •   Sea ice + ocean + atmosphere coupling

⸻

❄️ 4. Sea Ice, Arctic System, and Cryosphere Models

Useful for PBL–sea ice interactions and Arctic research.

📌 CICE / CICE Consortium
   •   Language: Fortran
   •   Strength: State-of-art sea-ice dynamics & thermodynamics
   •   License: Open

📌 MITgcm (with sea ice + ocean + atmospheric coupling options)
   •   Language: Fortran
   •   Strength: Arctic ocean + sea ice + limited atmospheric physics

📌 FV3GFS for polar regions
   •   Open-source core of NOAA’s Unified Forecast System (UFS)
   •   Used for polar prediction experiments

📌 Icepack
   •   Thermodynamic component of CICE (standalone)

⸻

🏞️ 5. Land Surface Models (LSMs)

Useful for Arctic BL surface flux and permafrost coupling.

📌 CLM (Community Land Model)
   •   Part of CESM
   •   Includes Arctic snow, permafrost, soil column models

📌 Noah-MP
   •   Part of WRF, standalone as well
   •   Widely used for snow/soil/permafrost processes

📌 JULES
   •   UK community model, open-source, strong snow physics

⸻

🧪 6. Column Models & Process Models (Great for ML Integration)

Lightweight 1D models ideal for prototyping ML closures.

📌 SCM from WRF / CESM (single column modes)
   •   Use: test PBL schemes in controlled forcing conditions
   •   Open and easily modified

📌 PyCLES
   •   Language: Python + Fortran kernels
   •   Focus: idealized LES/SCM; perfect for ML experiments

📌 MIMICA
   •   Language: Fortran
   •   Focus: Cloud microphysics + turbulence; used for Arctic mixed-phase clouds
   •   Partially open

📌 Mixed-layer models (various open repos)
   •   Simple slab ocean + BL models used for conceptual Arctic studies

⸻

🧭 7. Ocean–Atmosphere Coupled Arctic System Models

Many combine sea ice, ocean, and boundary-layer physics.

📌 NEMO–CICE system
   •   Fully open components
   •   Used in polar research

📌 COAWST
   •   Coupled ocean–atmosphere–waves–sediment system
   •   WRF + ROMS + SWAN
   •   Strong for coastal Arctic work

⸻

🌐 8. Air Quality / Chemistry Transport Models

Many interact with PBL height, turbulence, and surface fluxes.

📌 CMAQ
   •   Language: Fortran/C
   •   Open-source US EPA
   •   PBL: imports turbulence from WRF; key for Arctic haze

📌 GEOS-Chem (with GCHP)
   •   Open-source global CTM
   •   Supports Arctic aerosol studies, blowing-snow chemistry

📌 WRF-Chem
   •   Popular in Arctic air pollution / wildfire smoke transport

⸻

🛠️ 9. Hybrid ML–Physics & Data-Assimilation Frameworks

These aren’t models themselves but allow ML inside physical models.

📌 OpenFOAM (CFD + physics plug-ins)

📌 PySDM / microphysics ML toolkits

📌 NASA UFS + ML operator replacements (open repos)

📌 Clima.jl (Julia-based climate model)
   •   Built from scratch to be ML/differentiation-friendly
   •   Fully open-source
   •   Includes PBL physics modules designed for neural closures

⸻

🧩 10. Emerging ML-native or hybrid models

These hybrid approaches are gaining traction:

📌 ClimODE / ClimSim (ML-ready large datasets for turbulence & BL physics)
   •   Datasets derived from high-end climate model + LES
   •   Open-source, excellent for ML parameterization development

📌 Graph Neural PDE solvers (JAX/Julia repos)
   •   Experimental but promising for Arctic coastal turbulence and sea-ice leads

⸻

🎯 Summary Table (Quick Scan)

Purpose	Best Open Models
LES (turbulence)	PALM, DALES, MicroHH, UCLA-LES, OpenFOAM
Regional weather	WRF, MPAS, BRAMS
Global ESM	CESM, E3SM (mostly), ICON (partly), OpenIFS
Arctic system	CICE, FESOM, NEMO, MITgcm
Land/permafrost	CLM, Noah-MP, JULES
Column/SCM	WRF-SCM, CESM-SCM, PyCLES
Atmospheric chemistry	CMAQ, GEOS-Chem, WRF-Chem
Coupled coastal/Arctic	COAWST
ML-ready frameworks	Clima.jl, ClimSim, PySDM


⸻

If you want, I can (a) compare these models for Arctic PBL research specifically, (b) pick a subset suitable for ML-based turbulence parameterization work, or (c) build a shortlist based on your goals (ML closure, forecasting, LES training data, etc.).