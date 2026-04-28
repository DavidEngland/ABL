# Stratified Slab Column: Snow, Ice, Land, and Under-Ice Water

**File:** `param/core/slab_column.md`  
**Purpose:** Extend the surface-layer parameterization downward into a minimal, physically consistent 1D slab column for snow, ice, land, and an optional thin under-ice water layer.  
**Related:** `param/SCAFFOLDING.md`, `param/core/fluxes.md`, `param/core/regimes.md`, `param/core/mixing_length.md`

---

## 1. Design Intent

The goal is a **simple but realistic** lower-boundary model that can be coupled to
the grid-aware surface-layer scheme.

This module must:

1. Allow **snow thickness to be zero but never negative**.
2. Support **layered materials** with different thermal properties.
3. Handle two lower-boundary modes with one framework:
   - **land mode**: snow over soil, or bare soil
   - **ice-ocean mode**: snow over ice, bare ice, and optional thin under-ice water
4. Keep the air-side exchange in the existing MOST-based framework.
5. Use a **small number of prognostic slabs** rather than trying to resolve every
   fuzzy interface explicitly.
6. Leave room for **tower-based adaptive calibration** without letting machine
   learning replace the physics.

This is a **column model extension**, not a full snowpack model, land-surface
model, or ocean mixed-layer model.

---

## 2. Conceptual Structure

The parameterization is split into two coupled pieces:

### 2.1 Air-side exchange

The air column computes momentum, heat, moisture, and optional scalar exchange
between the lowest atmospheric level and the surface skin using the existing
surface-layer architecture.

### 2.2 Subsurface slab column

Below the surface skin, a 1D vertical column carries heat and selected material
properties through a small set of slabs:

```text
Atmosphere
   |
   |  MOST-based fluxes: tau, H, LE, F_q
   v
Surface skin (diagnosed Ts)
   |
   |  conductive / diffusive interface flux
   v
Snow slabs (0 to Ns)
   |
   |  conductive interface flux
   v
Underlying substrate
   |-- land mode: soil slabs or deep-soil boundary
   `-- ice-ocean mode: ice slabs + optional thin under-ice water slab
```

If snow thickness is zero, the snow block disappears and the surface skin sits
on the underlying land or ice substrate directly.

---

## 3. Surface-Type State Machine

Use a small deterministic state machine to prevent impossible configurations:

| State | Description | Active slabs |
| --- | --- | --- |
| `bare_land` | no snow, land exposed | soil only |
| `snow_land` | snow over land | snow + soil |
| `bare_ice` | no snow, ice exposed | ice only |
| `snow_ice` | snow over ice | snow + ice |
| `ice_ocean` | ice over ocean with optional under-ice layer | ice + water |
| `snow_ice_ocean` | snow over sea ice | snow + ice + water |

The state is diagnosed each timestep from snow thickness, ice thickness, and the
selected lower-boundary mode.

---

## 4. Prognostic and Diagnostic Variables

### 4.1 Prognostic state

For each material slab `i`:

| Symbol | Meaning | Units |
| --- | --- | --- |
| `h_i` | slab thickness | m |
| `T_i` | slab temperature | K |
| `rho_i` | bulk density | kg m^-3 |
| `c_i` | heat capacity | J kg^-1 K^-1 |
| `k_i` | thermal conductivity | W m^-1 K^-1 |
| `H_i` | enthalpy per unit mass | J kg^-1 |

Optional, by mode:

| Symbol | Meaning | Applies to |
| --- | --- | --- |
| `theta_w` | liquid water content or wetness proxy | wet snow |
| `S_w` | salinity | under-ice water |
| `q_i` | additional scalar carried in substrate | only if physically justified |

### 4.2 Column-level state

| Symbol | Meaning | Units |
| --- | --- | --- |
| `h_s` | total snow thickness | m |
| `h_ice` | total ice thickness | m |
| `T_s` | surface skin temperature | K |
| `mode` | land or ice-ocean lower boundary | -- |
| `state` | surface-type state from Section 3 | -- |

### 4.3 Atmospheric forcing inputs

| Symbol | Meaning | Units |
| --- | --- | --- |
| `U_a` | lowest-level wind speed | m s^-1 |
| `theta_a` | lowest-level potential temperature | K |
| `q_a` | lowest-level humidity or atmospheric scalar | varies |
| `SW_dn` | downward shortwave radiation | W m^-2 |
| `LW_dn` | downward longwave radiation | W m^-2 |
| `P_snow` | snowfall rate converted to thickness tendency | m s^-1 |
| `P_rain` | rainfall rate if needed for melt budget | kg m^-2 s^-1 |

### 4.4 Primary outputs

| Symbol | Meaning | Units |
| --- | --- | --- |
| `tau` | surface momentum flux | N m^-2 |
| `H` | sensible heat flux | W m^-2 |
| `LE` | latent heat flux | W m^-2 |
| `G` | conductive ground or ice flux into substrate | W m^-2 |
| `F_q` | scalar flux at the air-surface interface | varies |
| `T_s` | surface skin temperature | K |
| `h_s` | updated snow thickness | m |

---

## 5. Governing Equations

### 5.1 Generic slab conservation law

For any conserved slab quantity `chi_i`:

$$
\frac{d}{dt}\left(\rho_i h_i \chi_i\right)
= J_{i-\frac{1}{2}} - J_{i+\frac{1}{2}} + S_i
$$

where `J` is the vertical flux at each interface and `S_i` is a source or sink.

### 5.2 Heat in enthalpy form

For energy, solve in enthalpy form rather than pure temperature form:

$$
\frac{d}{dt}\left(\rho_i h_i H_i\right)
= F_{i-\frac{1}{2}} - F_{i+\frac{1}{2}} + Q_i
$$

with conductive interface flux

$$
F_{i+\frac{1}{2}} = -k_{i+\frac{1}{2}} \frac{T_{i+1} - T_i}{\Delta z_{i+\frac{1}{2}}}
$$

and interface conductivity computed with a harmonic mean:

$$
k_{i+\frac{1}{2}} = \left(\frac{\delta_i}{k_i} + \frac{\delta_{i+1}}{k_{i+1}}\right)^{-1}
$$

This gives a practical way to represent fuzzy material transitions without
introducing many thin layers.

### 5.3 Enthalpy-temperature relation

For a minimal phase-change treatment, define `H(T)` piecewise:

$$
H(T) =
\begin{cases}
c_{ice}(T - T_f), & T < T_f \\
L_f \lambda + c_{mix}(T - T_f), & T = T_f \text{ with melt fraction } \lambda \\
L_f + c_{liq}(T - T_f), & T > T_f
\end{cases}
$$

In the first implementation, this can be approximated by an apparent heat
capacity over a narrow temperature interval near freezing.

### 5.4 Snow-thickness evolution

Snow thickness is prognostic but constrained:

$$
\frac{d h_s}{dt} = P_{snow} - M - S_{sub} - C_{comp}
$$

where:

- `P_snow` is snowfall accumulation
- `M` is melt converted to thickness loss
- `S_sub` is sublimation or evaporation loss
- `C_comp` is compaction and settling

Constraint:

$$
h_s^{n+1} = \max(0, h_s^n + \Delta t \, \dot{h}_s)
$$

Negative snow thickness is never allowed.

### 5.5 Surface energy balance

The surface skin temperature `T_s` is diagnosed from the net energy balance:

$$
R_n - H - LE - G - Q_m = 0
$$

with

$$
R_n = SW_{dn}(1 - \alpha) + LW_{dn} - \epsilon \sigma T_s^4
$$

where:

- `H` is sensible heat flux to the atmosphere
- `LE` is latent heat flux to the atmosphere
- `G` is conductive flux into the top slab
- `Q_m` is latent energy consumed by melt or freezing at the surface

This nonlinear balance is solved each timestep for `T_s`.

### 5.6 Air-surface heat and scalar exchange

The air-side fluxes use bulk transfer coefficients from the existing MOST
scheme:

$$
H = \rho_a c_p C_H U_a (\theta_a - T_s)
$$

$$
LE = \rho_a L_v C_E U_a (q_a - q_s)
$$

$$
F_q = \rho_a C_q U_a (q_a - q_s)
$$

Momentum remains in the air-side module:

$$
\tau = \rho_a C_D U_a^2 = \rho_a u_*^2
$$

### 5.7 Transport of generic `q` below the surface

Below the snow or ice, only include a substrate scalar if the physics is clear.
The default form is diffusive:

$$
J_q = -K_q \frac{\partial q}{\partial z}
$$

Recommended use cases:

1. water vapor in porous snow
2. salinity in the thin under-ice water slab

Do **not** treat arbitrary atmospheric tracers as universal snow-column scalars
unless a specific snow chemistry or deposition model is being added.

---

## 6. Lower-Boundary Options

### 6.1 Land mode

The simplest land lower boundary is one or more soil slabs over a prescribed deep
temperature `T_deep`:

$$
F_{bottom} = -k_{soil} \frac{T_{N} - T_{deep}}{\Delta z_{deep}}
$$

For the MVP, one soil slab plus a fixed deep-soil boundary is sufficient.

### 6.2 Ice-ocean mode

For sea ice, use one or more ice slabs and optionally a thin under-ice water slab.
The under-ice slab represents the fresh or stratified layer that can exist after
salt rejection during freezing.

The minimal under-ice heat flux closure is:

$$
F_{iw} = h_w \rho_w c_w \frac{dT_w}{dt} = F_{ice\to w} - F_{w\to deep}
$$

with a relaxation or exchange closure to deeper water:

$$
F_{w\to deep} = \rho_w c_w C_w (T_w - T_{deep,ocean})
$$

If salinity is carried in the thin under-ice slab:

$$
h_w \frac{dS_w}{dt} = J_{salt,ice\to w} - J_{salt,w\to deep}
$$

This is enough to represent a freshened or stratified layer without a full ocean
model.

---

## 7. Default Complexity Levels

### 7.1 MVP

1. diagnosed `T_s`
2. one snow slab
3. one substrate slab: soil or ice
4. fixed deep boundary
5. air-side sensible, latent, and scalar fluxes from existing MOST code

### 7.2 Recommended research version

1. two or three snow slabs
2. one or two ice or soil slabs
3. apparent heat capacity or enthalpy phase change
4. optional under-ice water slab
5. merge-split logic for thin or disappearing snow

### 7.3 Not in the first implementation

1. full snow metamorphism microphysics
2. brine-pocket sea-ice thermodynamics
3. full multilayer soil hydrology
4. arbitrary chemistry in every slab

---

## 8. Pseudocode

### 8.1 Main column update

```python
def advance_slab_column(column, forcing, air_state, dt):
    """
    Advance snow/substrate slabs by one timestep.
    """
    state = diagnose_surface_state(column)

    # 1. Get air-side transfer coefficients from existing surface-layer scheme.
    flux_coeffs = compute_surface_exchange_coeffs(
        air_state=air_state,
        surface_state=state,
    )

    # 2. Solve surface energy balance for skin temperature.
    T_s = solve_surface_skin_temperature(
        column=column,
        forcing=forcing,
        flux_coeffs=flux_coeffs,
        air_state=air_state,
        dt=dt,
    )

    # 3. Compute atmosphere-surface fluxes using diagnosed skin state.
    tau, H, LE, F_q = compute_air_surface_fluxes(
        T_s=T_s,
        air_state=air_state,
        forcing=forcing,
        coeffs=flux_coeffs,
    )

    # 4. Compute conductive or diffusive interface fluxes through slabs.
    interface_fluxes = compute_internal_interface_fluxes(
        slabs=column.slabs,
        T_s=T_s,
        lower_mode=column.mode,
    )

    # 5. Advance enthalpy of each slab.
    for slab in column.slabs:
        slab.H = update_enthalpy(
            slab=slab,
            interface_fluxes=interface_fluxes,
            dt=dt,
        )
        slab.T = invert_enthalpy_to_temperature(slab.H, slab)

    # 6. Update snow thickness and remove extinct snow safely.
    column.h_s = max(0.0, update_snow_thickness(column, forcing, H, LE, dt))
    if column.h_s == 0.0:
        remove_snow_slabs(column)

    # 7. Apply merge or split rules if a slab becomes too thin or too thick.
    regularize_layer_layout(column)

    # 8. Optionally update online-calibrated parameters.
    update_adaptive_parameters(column, air_state, forcing, observations=None)

    return {
        "T_s": T_s,
        "tau": tau,
        "H": H,
        "LE": LE,
        "F_q": F_q,
        "G": interface_fluxes["surface_to_top_slab"],
    }
```

### 8.2 Surface skin solver

```python
def solve_surface_skin_temperature(column, forcing, flux_coeffs, air_state, dt):
    """
    Solve R_n - H - LE - G - Q_m = 0 for T_s.
    Use Newton iteration or a bracketed scalar solve.
    """
    T_guess = column.T_surface_prev

    for _ in range(MAX_ITERS):
        R_n = net_radiation(T_guess, forcing, column.surface_optics)
        H = sensible_heat_flux(T_guess, air_state, flux_coeffs)
        LE = latent_heat_flux(T_guess, air_state, flux_coeffs)
        G = top_conductive_flux(T_guess, column.top_slab)
        Q_m = surface_phase_change_sink(T_guess, column)

        residual = R_n - H - LE - G - Q_m
        if abs(residual) < ENERGY_TOL:
            return T_guess

        dR_dT = energy_balance_jacobian(T_guess, column, forcing, air_state, flux_coeffs)
        T_guess = T_guess - residual / dR_dT

    return enforce_temperature_bounds(T_guess, column)
```

### 8.3 Under-ice water option

```python
def update_under_ice_layer(layer, lower_flux, deep_ocean_state, dt):
    """
    Minimal prognostic update for a thin stratified layer below ice.
    """
    exchange = layer.exchange_coeff * (layer.T - deep_ocean_state.T)
    dTdt = (lower_flux - exchange) / (layer.rho * layer.cp * layer.h)
    layer.T = layer.T + dt * dTdt

    if layer.has_salinity:
        salt_exchange = layer.salt_exchange_coeff * (layer.S - deep_ocean_state.S)
        dSdt = (layer.salt_source - salt_exchange) / layer.h
        layer.S = clamp_salinity(layer.S + dt * dSdt)
```

---

## 9. Machine-Learning-Assisted Calibration

The machine-learning role is **parameter estimation**, not replacement of the
governing equations.

At a fixed tower site, adaptive calibration may update slowly varying uncertain
parameters such as:

1. `z0m`, `z0h`, or roughness-length ratios
2. snow conductivity multiplier
3. compaction-rate coefficient
4. surface scalar resistance
5. under-ice exchange coefficient

Hard constraints must always be enforced:

1. energy conservation
2. nonnegative snow thickness
3. positive conductivity and diffusivity
4. bounded salinity and liquid fraction
5. physically admissible roughness parameters

Good first choices are recursive least squares, EnKF, or Bayesian parameter
updating. A neural surrogate is optional later, not required for the first
implementation.

---

## 10. Diagnostics and Validation Targets

Recommended diagnostics:

| Diagnostic | Purpose |
| --- | --- |
| `T_s`, `T_i` profiles | check thermal structure and inversion strength |
| `h_s`, `h_ice` | ensure physically reasonable layer evolution |
| `H`, `LE`, `G`, `R_n` | verify surface energy closure |
| `tau`, `u_*`, `C_H`, `C_E` | connect substrate response to air-side coupling |
| melt energy and sublimation | diagnose phase-change budgets |
| under-ice `T_w`, `S_w` | assess stratified water option |

Validation priorities:

1. synthetic tests with known conductive solutions
2. snow onset and disappearance transitions
3. stable nighttime land cases with tower data
4. Arctic snow-on-ice cases with strong inversions
5. optional under-ice layer sensitivity tests

---

## 11. Recommended Initial Implementation

The most practical first build is:

1. existing air-side MOST fluxes
2. one snow slab with apparent heat-capacity phase change
3. one underlying slab: soil for land, ice for sea ice
4. fixed deep boundary temperature
5. optional site-specific adaptive calibration of `z0h` and snow conductivity

That configuration is simple enough to implement quickly while still capturing:

1. insulation by snow
2. zero-snow transitions
3. different lower-boundary behavior for land and sea ice
4. persistent stable coupling relevant to tower validation

---

## 12. Open Design Questions

1. Should phase change be handled by apparent heat capacity or exact enthalpy inversion in the first implementation?
2. Should the under-ice water slab be enabled only in sea-ice mode, or also in lake-ice mode?
3. Which scalar, if any, is worth carrying below the snow in version 1 besides heat?
4. Should adaptive calibration run continuously, or only in an offline reanalysis mode during initial development?
