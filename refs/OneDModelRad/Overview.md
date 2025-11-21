Planning explanation sections before summarizing model.

Created 8 todos

Wrapping up plan; marking all explanation tasks complete.

Completed: *Summarize model purpose and scope* (1/8)

**Overview**
- **Purpose:** 1-D (vertical column) boundary layer model (1DBLM) resolving surface layer → boundary layer thermodynamics, momentum, moisture, and radiation over time.
- **Regime Focus:** Stable nocturnal boundary layer with options for daytime (geostrophic wind driving mixing and insolation).
- **Core Equation Types:** Prognostic time stepping for potential temperature, specific humidity, wind components; diagnostic surface fluxes; vertical diffusion with a K-profile; shortwave radiative adjustments; cloud optical property parameterizations.

**Initialization**
- **Input Profiles:** Vertical arrays of height Z, potential temperature TH (T), wind components U,V, humidity Q, pressure P (read from an input file—trimmed lines suggest INPUT.DAT; geostrophic winds from BL.DAT).
- **BL.DAT Parameters:** Surface ground temperature `T0`, roughness length `Z0`, surface pressure `PRESS`, geostrophic winds `UG, VG`, latitude `LAT`, start time/date (`IMO, IDY, TIME`), time step `DELT`, total integration time `TTIME`, output intervals (`TINC`, `TINCFLUX`), physical constants (`VONKAR`, `G`), logicals for first-write.
- **Derived Grids:** Mid-level heights, layer thickness `dz`, perhaps `dzm` (midpoint spacing) computed in subroutine (truncated here but referenced).
- **Flags:** Choice of flux scheme `FLUXTYPE` (1: Businger; 2: England/McNider nocturnal + Businger daytime), stability function selection `STABFUNC`.

**Time Integration Loop**
- **Loop Counter:** `IT = IT + 1` each time step; runs until `TIME` reaches `TTIME`.
- **Hourly / Interval Output:** Controlled by `TINC`, `TEMPOUT`, and logical FIRST flags (handled in `store` and `save`).
- **Assimilation Hooks:** Placeholder commented out code for assimilating skin temperature tendency and insolation (MA_ADJ section).
- **Sequence (Typical):**
  1. Update time variables.
  2. Compute/adjust surface fluxes (USTAR, TSTAR, QSTAR, L).
  3. Radiation / insolation (SWNET or observed).
  4. Vertical diffusion (K-profile generation).
  5. Centered finite-difference integration of prognostics (e.g., temperature, momentum, moisture).
  6. Apply Coriolis tendencies (PGF term).
  7. Store outputs when interval met.

**Surface Flux & Stability**
- **Key Flux Variables:**
  - `USTAR` friction velocity.
  - `TSTAR` friction temperature scale.
  - `QSTAR` friction humidity scale.
  - `L` Monin–Obukhov length.
- **Schemes:**
  - England & McNider (quadratic Ri-based for stable nights) vs classical Businger similarity (log/psi functions).
  - Stability branching: uses `stability` and `ZOL = Z/L` logic for unstable (negative) vs stable (positive).
- **Computation Flow:**
  - Get aerodynamic temperature / surface specific humidity (`getsfcspechum`), then `getqstar`, `getustar`, `gettstar`, then `getl`.
  - Iterative coupling (OLDL, new L) ensures consistency; may re-enter flux calculation if convergence criteria not shown here.
- **Ri / Mixing:** Stable regime uses gradient Richardson number to diminish K coefficients; unstable uses convective scaling and boundary layer height.

**Boundary Layer Height**
- **Subroutine `blheight`:** For unstable conditions, sets or evolves BL height (e.g., initial guess `ZI = (0.33*USTAR)/|F|` if absent). This influences vertical diffusion limits and mixing length maxima.

**Vertical Diffusion (K-Profile)**
- **Unstable:** O’Brian polynomial form (O'Brian 1970) for eddy diffusivity rising to max at some fraction of BL height, then falling; background value above BL.
- **Stable:** Richardson-number dependent reduction. Computes shear (`SHEAR = sqrt(UDIFF^2 + VDIFF^2)`), mixing length, then K from stability functions chosen by `STABFUNC`.
- **Separate Km, Kh:** Sometimes unified (`C = Km = Kh` if specified) else computed individually.

**Radiation & Insolation**
- **Shortwave Net:** `SWNET` subroutine calculates downward solar flux using latitude, day-of-year, hour angle → cosine zenith → incoming flux; applies albedo and atmospheric absorption. Placeholder sets `XIS = 0.0` for testing (turning off SW input).
- **Gas Absorption:**
  - chou_routines_flcode.f implements Chou’s parameterizations for CO2 and O2 band absorption (bands 10 & 11), updating downward/upward flux reductions arrays `df`, `dfu`.
  - Pressure-scaled amounts and look-up tables (`cah`, `coa`) adjust flux profiles.
  - Cloud layer interaction adjusts flux fractions via `add_chou3` (albedo adjustments, cloud-reflected partition).
- **Cloud Optics (Shortwave):**
  - water_hu2_flcode.f provides water cloud optical properties using Mie parameterizations from Yong Hu.
  - `water_hu(ib)`: For spectral band `ib` (1–6 shortwave bands), obtains extinction coefficient `bext`, co-albedo `coalb`, asymmetry factor `asy` via:
    - `cldopt_l20` for effective droplet size ≤ 20 μm
    - `cldopt_g20` for > 20 μm (adds truncation correction `trc`, `trcasy`).
  - Outputs layer optical depth `tw(i) = dz(i)*plwc(i)*bext`, single scattering albedo `ww(i) = 1 - coalb`, and phase function expansion coefficients `www(i,1..4)` from powers of asymmetry factor for two-stream or multi-moment radiative solver.
- **Optical Depth ↔ Water Path:** `tau_wp` function converts between visible optical depth and water/ice path given particle size (different coefficients for water vs ice). Direction controlled by `IDIR`.

**Key Data Structures**
- **Cloud Common Blocks:** `/clouds/ pre, plwc, pde, piwc` effective radii and water contents; `/thick/ dz`; `/wat/ tw, ww, www`.
- **Radiation Arrays:** `dn, up` for flux downward/upward per level; `df` fractional reductions; `valb0` albedos; `dn1, up1` adjusted flux after absorption.
- **BL Variables:** `UG,V G` geostrophic wind; `USTAR,TSTAR,QSTAR,L,ZOL,ZI`; mixing lengths; shear.
- **Time Control:** `TIME`, `DELT`, `TTIME`, `TNEXT`, `TINC`, `TEMPOUT` from BL.DAT.

**Data Flow Summary**
- Initialization: read profiles → derive grid → set initial surface values.
- Each step: radiation (if daylight) supplies surface/column flux forcing; surface flux routine transforms near-surface gradients to turbulent flux scales (USTAR etc.); K-profile diffuses scalars/momentum; prognostic update changes state; outputs saved periodically.
- Cloud/radiative modules feed optical properties → flux attenuation → heating rates (implicitly if tendencies derived from flux divergence elsewhere).

**Running the Model**
- **Compile:** Likely via provided Makefile in root (contains Fortran sources). On macOS/zsh:
```sh
make
```
- **Run Executable:** (Assuming target name produced, e.g., `model`):
```sh
./model
```
- **Prepare Inputs:** Ensure profile file matches expected format and `PRESS` in BL.DAT equals surface pressure in surface data (comment notes).
- **Tuning:** Adjust `DELT` for stability (small for stiff nocturnal regime), `FLUXTYPE` for desired flux physics, `STABFUNC` for stability function shape, `TINC` for output density.
- **Radiation Activation:** Replace test line `XIS = 0.0` with computed `XIS` or assimilation variable to enable solar forcing.

**Extension Ideas**
- **Add Longwave Radiation:** Introduce surface energy balance improvements (currently focus on shortwave + turbulent).
- **Moist Processes:** Activate/extend cloud microphysics using `plwc`, `piwc` for latent heat tendencies.
- **Output Diagnostics:** Add budgets (TKE proxy, flux profiles) to `store`.
- **Adaptive BL Height:** Couple `ZI` evolution to entrainment flux if unstable daytime runs needed.

**Interpreting Results**
- **Primary Outputs:** Time series of near-surface temperature, fluxes (USTAR, sensible/latent), boundary layer height, vertical profiles.
- **Stable Nights:** Expect shallow mixing, low USTAR, strong decoupling; ensure `L` large and positive; check Richardson limiters in K-profile.
- **Cloud Impact:** Optical depth `tw` influences shortwave attenuation and heating; asymmetry factor shapes scattering contribution (phase expansion coefficients).

Let me know if you’d like a diagram, a variable map, or help enabling the radiation and cloud modifications in the main loop.