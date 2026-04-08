# SCMSkeleton.jl: Design Notes for Numerical Modelers

## Purpose

[SCMSkeleton.jl](SCMSkeleton.jl) is a minimal 1D atmospheric single-column model framework.
It is not meant to be scientifically complete. Its job is to provide a clear numerical structure
so new physics can be inserted without rewriting the rest of the model.

This is the right level of code when you want to answer questions like:

- Where does the turbulence closure plug in?
- Where should radiative cooling be applied?
- How do I replace prescribed surface fluxes with a surface scheme?
- How can I compare corrected vs. uncorrected diffusivity profiles in a controlled column model?

## State Variables

The current prognostic state is stored in `ColumnState`:

- `theta`: potential temperature profile
- `q`: specific humidity profile
- `u`, `v`: horizontal wind components

The lower boundary is represented separately by `SurfaceState` and slab metadata/state:

- `temperature`
- `sensible_flux`
- `latent_flux`
- `net_radiation`
- `ground_flux`

Plus the coupled slab components:

- `SurfaceSlabParameters`: material properties (depth, rho, cp, conductivity,
  albedo, emissivity, roughness lengths, moisture availability)
- `SurfaceSlabState`: prognostic `skin_temperature`, `deep_temperature`,
  and `liquid_fraction`
- `TowerSite`: observation geometry/metadata (`z_t`, `z_q`, `z_u`, lat/lon, terrain)

This split is deliberate. In a fuller SCM, the surface scheme often evolves on a different logic than the atmospheric column.

## Coordinate System: A Direct Answer

The grid is a **simple physical-z array**. There is no covariant/contravariant coordinate transform active in the current
solver. The exponential stretching is applied once at grid-construction time, and thereafter the physics works entirely in
physical metres.

The current implementation uses **cell-centered levels**:

$$\eta_k = \frac{k - 1/2}{N}, \quad k=1,\ldots,N$$

so the first prognostic level is above the surface (`z_1 > 0`) and `z=0` remains a boundary interface.

However, the **Jacobian** J = dz/deta is now stored explicitly in `grid.jacobian`:

$$J(\eta) = z_\text{top} \cdot s \cdot \frac{e^{s\eta}}{e^s - 1}$$

where $s$ is the stretch parameter and $\eta \in [0,1]$ is the uniform computational coordinate.

Numerically, the grid builder uses `expm1` and a small-$s$ branch:

- `|s| < 1e-8`: uniform limit, `z = z_top * eta`, `J = z_top`
- otherwise: `z = z_top * expm1(s*eta) / expm1(s)`

This avoids cancellation in the `s -> 0` limit.

This is where the covariant/contravariant terms **would first appear** in a terrain-following upgrade:

| Operation | Physical z (current) | η-coordinate upgrade |
|---|---|---|
| Vertical gradient | $\partial\phi/\partial z$ | $(1/J)\,\partial\phi/\partial\eta$ |
| Flux divergence | $\partial F/\partial z$ | $(1/J)\,\partial(JF)/\partial\eta$ |
| Laplacian | $\partial^2\phi/\partial z^2$ | $(1/J)\,\partial[(1/J)\,\partial\phi/\partial\eta]/\partial\eta$ |

Keeping `grid.jacobian` populated means a terrain-following rewrite only touches the diffusion functions,
not the timestepper or the closure interface.

## Data Containers

The model is divided into a few small containers.

### `ModelConfig`

Static parameters for one experiment:

- domain height,
- number of levels,
- timestep,
- run duration,
- initial thermodynamic structure,
- a few bulk constants.

### `Grid`

Stores the vertical geometry:

- `z`: cell-center heights,
- `dz`: representative layer thickness,
- `nz`: level count,
- `jacobian`: analytical `dz/deta` at cell centers.

The current grid is exponentially stretched toward the surface so the lowest layers are resolved more tightly.

### `Forcing`

External tendencies for one timestep:

- column tendencies for `theta`, `q`, `u`, `v`,
- optional prescribed surface heat/moisture fluxes,
- downwelling shortwave and longwave radiation,
- reference atmospheric state (`air_temperature_ref`, `specific_humidity_ref`,
  `wind_speed_ref`, `surface_pressure`),
- switch `prescribed_surface_fluxes` to select prescribed vs resolved coupling.

This is where synoptic forcing, radiative cooling, geostrophic relaxation, nudging, or idealized experiments should enter.

### `SimulationHistory`

Small diagnostic package recorded over time. It currently stores enough to plot a first run quickly without saving full 2D output arrays.

## Solver Sequence

The timestepper in [SCMSkeleton.jl](SCMSkeleton.jl) follows a simple sequence.

1. The closure returns vertical diffusivity profiles `km` and `kh`.
2. Diffusive tendencies are computed in flux-divergence form.
3. External tendencies from `Forcing` are added.
4. Surface fluxes are either prescribed or resolved from slab/tower state and then applied to the first model level.
5. Surface diagnostics are updated.
6. Model time advances by one step.

This is intentionally explicit. It makes the code easy to audit, but it also means stability constraints will become important when diffusivities or resolution increase.

## Explicit vs. Implicit Diffusion

Set `use_implicit = true` in `ModelConfig` to switch from forward-Euler to backward-Euler.

The explicit diffusion CFL constraint is:

$$\Delta t \leq \frac{\Delta z^2}{2 K}$$

In a polar ABL with $\Delta z = 1$ m and $K = 0.5$ m²/s, this forces $\Delta t \leq 1$ s, making a 12-hour simulation
require 43200 steps. With intermittent K bursts, the model crashes.

The implicit backward-Euler solve is **unconditionally stable**. The tridiagonal system assembled each step is:

$$(I - \Delta t \, L)\,\phi^{n+1} = \phi^n$$

where $L$ is the second-order diffusion operator. This is O(N) per variable per step via
`LinearAlgebra.Tridiagonal`. The extra cost over explicit is negligible; the numerical benefit is large.

### Crank-Nicolson (future)
A natural next upgrade is Crank-Nicolson:
$$(I - \tfrac{1}{2}\Delta t\,L)\,\phi^{n+1} = (I + \tfrac{1}{2}\Delta t\,L)\,\phi^n$$
which is second-order accurate in time and still unconditionally stable. It requires two tridiagonal solves
but eliminates the first-order temporal damping of backward-Euler.

## Closure Interface

The most important design choice is the abstract closure hook.

The model defines:

```julia
abstract type AbstractClosure end
```

and expects each closure to provide:

```julia
diffusivities(closure, model)
```

which must return two profiles:

- `km`: momentum diffusivity
- `kh`: scalar diffusivity

This keeps the closure logic separate from the time integrator. That matters because your research question is likely to change the closure much more often than the rest of the numerics.

## Example: Where a Ri-Based Closure Would Go

Suppose you want to add a Richardson-number-dependent closure with your curvature-aware correction.

You would:

1. Create a new subtype of `AbstractClosure` containing the needed parameters.
2. Diagnose shear and stratification from `model.state` and `model.grid`.
3. Compute `Ri`, `Ri_g`, or related stability diagnostics.
4. Form `km(z)` and `kh(z)` from the chosen closure formula.
5. Return those profiles through `diffusivities`.

Nothing in `step!` needs to change.

That separation is the main reason to prefer this skeleton over a single monolithic script.

## RiBasedClosure: The First Physics Closure

`RiBasedClosure(K0, Ri_c, Pr_t)` implements local Richardson-number-dependent diffusion using a Webb (1970) quadratic stability reduction:

$$K_m(z) = K_0 \cdot \max\left(0,\, 1 - \frac{Ri(z)}{Ri_c}\right)^2$$

$$K_h(z) = K_m(z) / Pr_t$$

where $Ri(z) = N^2/S^2$ is diagnosed by centered differences at each interior level.

This is the minimum closure needed to observe stable-case inversion growth and compare corrected vs.
uncorrected mixing. The extension point for your curvature-aware correction is inside
`diffusivities(::RiBasedClosure, model)` — compute the bias ratio $B = Ri_g / Ri_b$ and the $f_c$
correction factor there, then multiply `km` and `kh` by `fc` before returning.

## TKE Budget Diagnostics

`SimulationHistory` now records:
- `max_shear_production`: column maximum of $P_s = K_m \cdot S^2$ (m²/s³)
- `max_buoyant_destruction`: column maximum of $B = (g/\theta_0) \cdot K_h \cdot N^2$ (m²/s³)

In near-equilibrium turbulence: $P_s \approx B + \varepsilon$. When a closure shuts off $K_m$ too
aggressively, $P_s$ drops while $B$ stays elevated. The ratio $P_s/B$ falling far below unity is the
numerical fingerprint of runaway cooling — the exact failure mode you are targeting in polar ABL work.

The companion function `compute_tke_budget(model, km, kh)` returns the full column profiles so a notebook
can plot the vertical structure of the budget, not just the column maximum.

## Surface Flux Handling

Surface fluxes are currently prescribed through the `Forcing` object and applied directly to the first level.

That is numerically simple, but physically incomplete. A more realistic next step would be to replace this with a Monin-Obukhov or bulk-aerodynamic surface routine that uses:

- near-surface wind,
- surface temperature,
- roughness lengths,
- stability functions.

At that point, `SurfaceState` would likely need to hold additional variables such as roughness, skin temperature, and possibly soil or snow properties.

## Numerical Limitations

The current implementation is a scaffold, not a verified solver. Important limitations are:

- no monotonic advection scheme is included (diffusion-only transport skeleton),
- top and bottom diffusive fluxes are zero unless prescribed separately,
- moisture is not coupled to latent heating,
- no pressure or density evolution is included,
- no Coriolis term or momentum forcing is active,
- diagnostics are intentionally minimal.

For serious stable boundary-layer work, the first numerical upgrade should probably be an implicit or semi-implicit vertical diffusion solve.

## Suggested Development Roadmap

### Stage 1: Clean research prototype

- Add a stability-dependent closure.
- Add shear and Richardson number diagnostics.
- Add profile output for `km`, `kh`, and `Ri`.
- Add a simple plotting notebook.

### Stage 2: Physically better SCM

- Replace prescribed fluxes with a surface-layer routine.
- Add longwave cooling or idealized radiative forcing.
- Add Coriolis forcing and geostrophic relaxation.
- Add an implicit diffusion step.

### Stage 3: Comparison framework

- Run corrected vs. uncorrected closures.
- Compare inversion strength, jet height, and mixing depth.
- Export results to CSV or NetCDF.
- Use the same model core in a Pluto notebook for interactive parameter sweeps.

## Minimal Usage Pattern

```julia
include("SCMSkeleton.jl")
using .SCMSkeleton

# Explicit solver, constant K (simplest possible run)
config  = ModelConfig(z_top=500.0, nz=30, dt=20.0, t_end=6.0 * 3600.0)
closure = ConstantDiffusivityClosure(0.6, 0.4)
model, history = run_model(config; closure=closure, log_every=60)

# Implicit solver + Ri-based closure (research-grade stable ABL case)
config2  = ModelConfig(z_top=1000.0, nz=60, dt=60.0, t_end=12.0*3600.0, use_implicit=true)
closure2 = RiBasedClosure(1.0, 0.25, 1.3)   # K0=1, Ri_c=0.25, Pr_t=1.3
model2, h2 = run_model(config2; closure=closure2, log_every=120)
```

The `history` object from the second run provides a time series of `max_shear_production`
and `max_buoyant_destruction` ready for direct comparison against your expectations for
stable-case TKE budget balance.

## Bottom Line

[SCMSkeleton.jl](SCMSkeleton.jl) is best understood as a numerically clean insertion point for new physics. It gives you a stable code shape now, so future closure work does not turn into a structural rewrite.