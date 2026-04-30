Data Access Request for Sodankyla ABL Analysis (IMPECCABLE): Momentum, Heat, and Tracers

Form Fields (suggested)
1. Full name: David E. England, PhD
2. Email: dee0001@uah.edu
3. Affiliation: University of Alabama in Huntsville (UAH)
4. Role/title: Research Engineer III (contracted)
5. Project: IMPECCABLE (Improving Model Parameterizations of the Continental Arctic Boundary Layer Ecosystem)
6. Purpose: Boundary-layer physics and Richardson-number correction development/validation
7. Site: Sodankyla, Finland
8. Time period requested: Winter seasons (please provide available years; preferred continuous multi-year coverage)
9. Intended use: Research and model development (non-commercial), publication and model validation
10. Output format requested: NetCDF preferred, CSV acceptable

Request Text
I am requesting access to Sodankyla atmospheric and surface datasets for research under the IMPECCABLE project.
My goal is to evaluate and improve Richardson-number-based boundary-layer stability corrections in Arctic/sub-Arctic winter conditions using a staged workflow: momentum calibration first, then heat transfer, then additional tracers.

Please provide (or advise availability for) the following data streams at the highest practical cadence, with QC flags and metadata:
1. Wind components and turbulence fluxes (u, v, w, u'w', w'theta', friction velocity where available).
2. Thermodynamic profiles (temperature, humidity, pressure; tower levels and/or profile products).
3. Surface energy budget terms (radiation components, sensible/latent heat fluxes, ground/snow heat flux where available).
4. Stability diagnostics or ingredients needed to compute them (Ri_b, Ri_g inputs, roughness/displacement metadata if available).
5. Tracer-relevant fields if available (water vapor products, CO2, CH4, or related concentration/flux observations).
6. Cloud/context observations useful for stability-regime tagging (e.g., ceilometer/cloud-base or all-sky products).

Intended outputs are:
1. Calibrated momentum and heat closures.
2. Tracer transfer-ratio fits using a shared momentum backbone.
3. Improved Ri_b-to-Ri_g effective mapping and validated Richardson corrections against observed winter regimes.

I will comply with all data-license, attribution, and citation requirements.
Please let me know the required acknowledgement text and any restrictions on redistribution or derivative products.

Thank you for your support.

Where to submit
1. FMI Open Data support/request channel: https://en.ilmatieteenlaitos.fi/open-data-support-request
2. FMI Download Observations portal: https://en.ilmatieteenlaitos.fi/download-observations
3. Sodankyla Geophysical Observatory data pages (archive/licence/contact): https://www.sgo.fi/Data/archive.php and https://www.sgo.fi/Data/licence.php

Attachment A - Requested Variables Table (priority, units, cadence, substitutes)

| Group | Variable(s) | Priority | Preferred Unit(s) | Minimum cadence | Acceptable substitutes / notes |
|---|---|---|---|---|---|
| Momentum core | U, V (or wind speed + direction) | P1 | m s^-1, deg | 10 min (1 min preferred) | If no U/V, speed+direction accepted; derive U/V |
| Momentum turbulence | w, u'w', v'w', friction velocity u* | P1 | m s^-1, m^2 s^-2 | 30 min flux blocks (10 Hz raw preferred) | If covariance fluxes unavailable, provide eddy diffusivity or TKE proxies |
| Heat core | Air temperature profile (multi-level) | P1 | degC or K | 10 min (1 min preferred) | Single-level accepted only as fallback; note measurement height |
| Moisture core | RH and/or specific humidity q | P1 | %, kg kg^-1 | 10 min | Dew point + pressure acceptable to derive q |
| Pressure | Surface pressure p | P1 | hPa or Pa | 10 min | Station pressure preferred over sea-level pressure |
| Stability ingredients | Virtual temperature or inputs for theta_v, layer heights z | P1 | K, m | 10 min | Dry theta accepted if humidity available to reconstruct theta_v |
| Radiation | SWdown, SWup, LWdown, LWup, net radiation | P1 | W m^-2 | 10 min | If only partial radiation set, provide available components with metadata |
| Surface exchange | Sensible heat flux H, latent heat flux LE | P1 | W m^-2 | 30 min | Bowen ratio products acceptable if H/LE not directly available |
| Ground/snow | Ground heat flux G, snow depth, snow temperature | P2 | W m^-2, m, degC | 10 min to hourly | Soil/snow temperature profile accepted if G unavailable |
| Soundings | PTU/radiosonde profiles (T, RH, wind) | P2 | K, %, m s^-1 | 00/12 UTC minimum | Additional launches strongly preferred for event nights |
| Clouds/context | Cloud base height, cloud fraction, ceilometer class | P2 | m, % | 10 min | All-sky image products acceptable for cloud-state tagging |
| Tracer - CO2 | CO2 concentration and/or flux | P2 | ppm, umol m^-2 s^-1 | 30 min (10 min preferred) | Tower-level concentration accepted if flux unavailable |
| Tracer - CH4 | CH4 concentration and/or flux | P2 | ppb, nmol m^-2 s^-1 | 30 min | Concentration-only acceptable for first-pass transfer diagnostics |
| Ancillary metadata | Station coordinates, sensor height, canopy height, displacement d, roughness z0 | P1 | m, deg | Static + change log | Required for Ri_b to Ri_g mapping and inter-site consistency |
| Data quality | QC flags, processing level, gap-fill indicators | P1 | code list | with each record | Mandatory for calibration and uncertainty estimation |

Suggested first stations/products from the Arctic Space Centre list for rapid start:
1. MET0002 (micrometeorological mast, wind + radiation + snow depth)
2. MET0006 (CO2 flux and turbulence)
3. VUO0002 (CO2 mast + radiation)
4. SUO0010 (peatland methane station, humidity/radiation/soil context)
5. LUO0015 (operative weather station baseline meteorology)
6. LUO0001 (PTU soundings)

Minimum deliverable package for Phase 1 (momentum -> heat):
1. U, V, T, RH, p at known heights
2. Radiation (at least SWdown and LWdown)
3. H and LE (or substitutes)
4. QC flags and complete station metadata
