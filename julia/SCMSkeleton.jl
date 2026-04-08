"""
    SCMSkeleton

Minimal one-dimensional single-column model scaffold for atmospheric boundary-layer
experiments. The module is intentionally simple: it separates grid generation,
state storage, forcing, closure logic, time stepping, and diagnostics so that a
more realistic turbulence scheme can be dropped in without restructuring the code.

This is a teaching and prototyping skeleton, not a production-ready SCM.
"""
module SCMSkeleton

using Printf
using LinearAlgebra   # Tridiagonal type used by the implicit diffusion solver

const KAPPA_VK = 0.4
const SIGMA_SB = 5.670374419e-8
const LV_WATER = 2.5e6

export Grid,
       ModelConfig,
       ColumnState,
       SurfaceState,
    SurfaceSlabParameters,
    SurfaceSlabState,
    TowerSite,
       Forcing,
       SCMModel,
       SimulationHistory,
       AbstractClosure,
       ConstantDiffusivityClosure,
       RiBasedClosure,
       CurvatureRiClosure,
       create_grid,
    default_surface_material,
    default_tower_site,
       initialize_model,
       zero_forcing,
       default_forcing,
    resolve_surface_fluxes!,
       step!,
       run_model,
       example_run

"""
    Grid

Vertical grid description for the 1D column.

Fields:
- `z`:        cell-center heights in meters
- `dz`:       representative layer thickness at each cell center (m)
- `nz`:       number of vertical levels
- `jacobian`: local stretch factor J = dz_physical/dη at each cell center.
              For the exponential map z(η) = z_top·[eˢη−1]/(eˢ−1),
              J = z_top·s·eˢη/(eˢ−1). On a uniform grid J = z_top.
              This is the coefficient that would appear in a covariant
              divergence operator (1/J)·∂/∂η in terrain-following coords.
              Storing J in the grid struct keeps metric terms available
              without recomputing the grid.
"""
struct Grid
    z::Vector{Float64}
    dz::Vector{Float64}
    nz::Int
    jacobian::Vector{Float64}
end

"""
    ModelConfig

Configuration parameters controlling grid size, timestep, initial thermodynamic
state, and a few bulk constants. These are deliberately small in number so that
the model remains easy to reason about while the physics is still being designed.
"""
struct ModelConfig
    z_top::Float64
    nz::Int
    dt::Float64
    t_end::Float64
    theta_surface::Float64
    lapse_rate::Float64
    q_surface::Float64
    geostrophic_wind::Float64
    rho_air::Float64
    cp_air::Float64
    use_implicit::Bool  # true → backward-Euler diffusion; false → explicit forward-Euler
end

"""
    ColumnState

Prognostic column variables. The current skeleton evolves potential temperature,
specific humidity, and horizontal wind components.
"""
mutable struct ColumnState
    theta::Vector{Float64}
    q::Vector{Float64}
    u::Vector{Float64}
    v::Vector{Float64}
end

"""
    SurfaceState

Minimal lower-boundary state. A more complete model would replace this with an
explicit land, snow, ice, or ocean surface scheme.
"""
mutable struct SurfaceState
    temperature::Float64
    sensible_flux::Float64
    latent_flux::Float64
    net_radiation::Float64
    ground_flux::Float64
end

"""
    SurfaceSlabParameters

Physical parameters for the lower boundary slab. This supports water, ice,
sea-ice, and user-defined molten/liquid surfaces through parameter choices.
"""
struct SurfaceSlabParameters
    material::Symbol
    depth::Float64
    rho::Float64
    cp::Float64
    k::Float64
    albedo::Float64
    emissivity::Float64
    z0m::Float64
    z0h::Float64
    moisture_availability::Float64
end

"""
    SurfaceSlabState

Prognostic slab state used by the coupled surface-energy balance.
"""
mutable struct SurfaceSlabState
    skin_temperature::Float64
    deep_temperature::Float64
    liquid_fraction::Float64
end

"""
    TowerSite

Minimal tower metadata and reference heights used by the bulk flux scheme.
"""
struct TowerSite
    name::String
    latitude::Float64
    longitude::Float64
    elevation::Float64
    terrain::String
    z_t::Float64
    z_q::Float64
    z_u::Float64
end

"""
    Forcing

External tendencies and lower-boundary fluxes applied over a single timestep.
Large-scale advection, radiative cooling, nudging, or prescribed fluxes should
be represented here.
"""
mutable struct Forcing
    theta_tendency::Vector{Float64}
    q_tendency::Vector{Float64}
    u_tendency::Vector{Float64}
    v_tendency::Vector{Float64}
    surface_heat_flux::Float64
    surface_moisture_flux::Float64
    sw_down::Float64
    lw_down::Float64
    air_temperature_ref::Float64
    specific_humidity_ref::Float64
    wind_speed_ref::Float64
    surface_pressure::Float64
    prescribed_surface_fluxes::Bool
end

"""
    SCMModel

Top-level container that carries configuration, grid, evolving state, surface
state, and the current model time.
"""
mutable struct SCMModel
    config::ModelConfig
    grid::Grid
    state::ColumnState
    surface::SurfaceState
    slab_params::SurfaceSlabParameters
    slab::SurfaceSlabState
    tower::TowerSite
    time::Float64
end

"""
    SimulationHistory

Small diagnostics bundle recorded during `run_model`. This is intentionally
limited to the most useful quantities for quick plotting and debugging.
"""
mutable struct SimulationHistory
    time::Vector{Float64}
    theta_surface::Vector{Float64}
    theta_first_level::Vector{Float64}
    q_first_level::Vector{Float64}
    boundary_layer_proxy::Vector{Float64}
    max_shear_production::Vector{Float64}    # column-max P_s = Kₘ·S²  (m²/s³)
    max_buoyant_destruction::Vector{Float64} # column-max B  = (g/θ₀)·Kₕ·N²  (m²/s³)
end

abstract type AbstractClosure end

"""
    ConstantDiffusivityClosure(km, kh)

Placeholder turbulence closure returning vertically uniform diffusivities for
momentum and scalars. This is useful for testing the solver wiring before a more
realistic stability-dependent closure is added.
"""
struct ConstantDiffusivityClosure <: AbstractClosure
    km::Float64
    kh::Float64
end

"""
    RiBasedClosure(K0, Ri_c, Pr_t)

Local Richardson-number-dependent K-profile closure.

The stability reduction follows a Webb (1970) style quadratic decay:

    f(Ri) = max(0, 1 − Ri / Ri_c)²

so that K → K0 in neutral conditions (Ri = 0) and K → 0 as Ri → Ri_c.
This is the simplest closure that respects the critical Richardson number.

Fields
──────
- `K0`:   neutral diffusivity (m²/s)
- `Ri_c`: critical Richardson number (default 0.25)
- `Pr_t`: turbulent Prandtl number Kₘ/Kₕ (default 1.0; set > 1 for stable ABL)

Extension point
───────────────
To add your curvature-aware grid correction, compute the bias ratio B = Ri_g/Ri_b
and the fc factor inside `diffusivities(::RiBasedClosure, model)`, then multiply
km and kh by fc before returning.
"""
struct RiBasedClosure <: AbstractClosure
    K0::Float64
    Ri_c::Float64
    Pr_t::Float64
end

RiBasedClosure(K0::Float64) = RiBasedClosure(K0, 0.25, 1.0)

"""
    CurvatureRiClosure(K0, Ri_c, Pr_t, beta)

Extension of `RiBasedClosure` that multiplies the stability-reduced diffusivity
by a grid-dependent correction factor $f_c$ derived from the curvature of the
potential-temperature profile. This is the direct implementation of the
curvature-aware bias correction described in England & McNider (in review).

Physics
───────
At a given level k the local Richardson-number bias is driven by grid-scale
curvature: the bulk Ri (Ri_b) computed over a thick layer systematically
under- or over-estimates the gradient Ri (Ri_g) when the θ and wind profiles
are nonlinear. The correction factor

    fc = clamp(1 + β · κ̃, 0.1, 2.0)

where κ̃ = (d²θ/dz²) · z / max(dθ/dz, ε) is a nondimensional curvature
index, adjusts K upward (fc > 1) when a sharpening inversion would otherwise
cause numerical decoupling, and downward (fc < 1) when the gradient is
already over-estimated.

Fields
──────
- `K0`:   neutral diffusivity (m²/s)
- `Ri_c`: critical Richardson number (default 0.25)
- `Pr_t`: turbulent Prandtl number Kₘ/Kₕ (default 1.0)
- `beta`: curvature sensitivity factor β. Start at 0.5; tune against GABLS1.

Extension point
───────────────
Replace the `curvature_bias` expression below with your Ri_g/Ri_b ratio
formulation when that derivation is finalised.
"""
struct CurvatureRiClosure <: AbstractClosure
    K0::Float64
    Ri_c::Float64
    Pr_t::Float64
    beta::Float64
end

CurvatureRiClosure(K0::Float64) = CurvatureRiClosure(K0, 0.25, 1.0, 0.5)

"""
    ModelConfig(; kwargs...)

Construct a default configuration for an idealized stable boundary-layer test.
All arguments are keyword-based so notebooks can expose them directly as sliders
or form inputs.
"""
function ModelConfig(; z_top::Float64=1000.0,
                       nz::Int=40,
                       dt::Float64=10.0,
                       t_end::Float64=12.0 * 3600.0,
                       theta_surface::Float64=265.0,
                       lapse_rate::Float64=0.01,
                       q_surface::Float64=1.5e-3,
                       geostrophic_wind::Float64=8.0,
                       rho_air::Float64=1.2,
                       cp_air::Float64=1004.0,
                       use_implicit::Bool=false)
    return ModelConfig(z_top, nz, dt, t_end, theta_surface, lapse_rate,
                       q_surface, geostrophic_wind, rho_air, cp_air, use_implicit)
end

"""
    default_tower_site()

Return a generic tower configuration for initial testing. Replace with actual
station metadata for observation-constrained runs.
"""
function default_tower_site()
    return TowerSite("Generic Tower", 0.0, 0.0, 0.0, "unknown", 2.0, 2.0, 10.0)
end

"""
    default_surface_material([material])

Return slab parameters for a named surface type.
Supported materials: `:water`, `:ice`, `:seaice`, `:molten_liquid`.
"""
function default_surface_material(material::Symbol=:water)
    if material == :water
        return SurfaceSlabParameters(:water, 1.0, 1000.0, 4180.0, 0.6, 0.07, 0.98, 2e-4, 2e-5, 1.0)
    elseif material == :ice
        return SurfaceSlabParameters(:ice, 1.0, 917.0, 2100.0, 2.2, 0.65, 0.99, 1e-3, 2e-4, 0.2)
    elseif material == :seaice
        return SurfaceSlabParameters(:seaice, 1.0, 930.0, 2500.0, 1.8, 0.50, 0.99, 8e-4, 2e-4, 0.5)
    elseif material == :molten_liquid
        return SurfaceSlabParameters(:molten_liquid, 0.5, 2500.0, 1200.0, 1.5, 0.12, 0.95, 5e-4, 5e-5, 0.7)
    end
    error("Unsupported surface material: $material")
end

"""
    create_grid(config; stretch=3.0)

Create a stretched vertical grid with finer resolution near the surface.

The grid uses **cell-centred levels**: N interfaces span η ∈ [0,1] uniformly, and
the prognostic levels sit at the cell centres

    ηᵢ = (i − ½) / N,    i = 1 … N

so that the lowest level z₁ > 0 (above the ground).  The surface z = 0 and the
domain top z = z_top are interface levels that appear only as boundary conditions,
not as prognostic nodes.

Physical heights come from the exponential map applied to the cell-centre η values:

    z(η) = z_top · [exp(s·η) − 1] / [exp(s) − 1],    s = stretch

  • s > 0 → fine near the surface, coarse aloft
  • s = 0 → uniform spacing (z₁ = z_top/(2N))

For numerical robustness near s = 0, the implementation uses `expm1` and
falls back to the exact uniform-grid limit when |s| is very small.

The Jacobian J = dz/dη is computed analytically and stored in `grid.jacobian`.
It answers the coordinate-system question directly: the grid is a simple z-array
in physical space, but the _stretch function_ induces a non-trivial metric.
In a terrain-following (η-coordinate) model, the operators ∂/∂z and ∇·F become
(1/J)·∂/∂η and (1/J)·∂(J·F)/∂η respectively.  The current solver works in
physical z directly so no covariant terms appear, but storing J ensures they
are available painlessly when a terrain-following upgrade is needed.
"""
function create_grid(config::ModelConfig; stretch::Float64=3.0)
    config.nz >= 2 || error("create_grid requires nz >= 2")
    # Cell-centred levels: ηᵢ = (i − ½)/N, so z₁ > 0 (surface is a BC, not a node)
    eta   = [(i - 0.5) / config.nz for i in 1:config.nz]
    if abs(stretch) < 1e-8
        # Uniform-grid limit of exp-stretch map: z = z_top * η and dz/dη = z_top.
        z = config.z_top .* eta
        jac = fill(config.z_top, config.nz)
    else
        # Use expm1 for robust evaluation when stretch is small.
        denom = expm1(stretch)
        z = config.z_top .* expm1.(stretch .* eta) ./ denom
        # Analytical Jacobian: dz/dη = z_top * s * exp(s*η) / (exp(s)-1)
        jac = config.z_top .* stretch .* exp.(stretch .* eta) ./ denom
    end
    dz    = zeros(config.nz)
    dz[1] = z[2] - z[1]
    for i in 2:config.nz-1
        dz[i] = 0.5 * (z[i + 1] - z[i - 1])
    end
    dz[end] = z[end] - z[end - 1]
    return Grid(z, dz, config.nz, jac)
end

"""
    initialize_column_state(config, grid)

Build a simple initial condition: weakly stable potential temperature,
exponentially decaying humidity, and geostrophic wind aligned with the `u`
component.
"""
function initialize_column_state(config::ModelConfig, grid::Grid)
    theta = [config.theta_surface + config.lapse_rate * height for height in grid.z]
    q = [config.q_surface * exp(-height / 1200.0) for height in grid.z]
    u = fill(config.geostrophic_wind, grid.nz)
    v = fill(0.0, grid.nz)
    return ColumnState(theta, q, u, v)
end

"""
    initialize_model([config])

Allocate and initialize the grid, prognostic state, and lower-boundary state.
This is the main constructor numerical experiments should call before stepping.
"""
function initialize_model(config::ModelConfig=ModelConfig();
                          surface_material::Symbol=:water,
                          tower::TowerSite=default_tower_site())
    grid = create_grid(config)
    state = initialize_column_state(config, grid)
    slab_params = default_surface_material(surface_material)
    deep_temperature = config.theta_surface - 1.0
    liquid_fraction = surface_material == :water ? 1.0 : 0.0
    slab = SurfaceSlabState(config.theta_surface, deep_temperature, liquid_fraction)
    surface = SurfaceState(config.theta_surface, 0.0, 0.0, 0.0, 0.0)
    return SCMModel(config, grid, state, surface, slab_params, slab, tower, 0.0)
end

"""
    zero_forcing(model)

Return a forcing object with all tendencies and fluxes set to zero. This is a
useful starting point for custom forcing functions.
"""
function zero_forcing(model::SCMModel)
    nz = model.grid.nz
    return Forcing(zeros(nz), zeros(nz), zeros(nz), zeros(nz),
                   0.0, 0.0,
                   0.0, 0.0,
                   NaN, NaN, NaN, 101325.0,
                   false)
end

"""
    default_forcing(model)

Provide a very simple stable-case forcing: weak column cooling plus a constant
negative surface heat flux. This exists only to make the skeleton runnable.
"""
function default_forcing(model::SCMModel)
    forcing = zero_forcing(model)
    hour = model.time / 3600.0
    cooling = hour <= 6.0 ? -2.0e-5 : -5.0e-6
    forcing.theta_tendency .= cooling
    # Simple synthetic radiation cycle used when observational forcing is absent.
    dayfrac = (model.time / 86400.0) % 1.0
    forcing.sw_down = max(0.0, 180.0 * sin(2pi * dayfrac))
    forcing.lw_down = 260.0
    forcing.air_temperature_ref = model.state.theta[1]
    forcing.specific_humidity_ref = model.state.q[1]
    forcing.wind_speed_ref = max(0.5, hypot(model.state.u[1], model.state.v[1]))
    forcing.surface_pressure = 101325.0
    forcing.prescribed_surface_fluxes = false
    return forcing
end

@inline clamp01(x::Float64) = max(0.0, min(1.0, x))

"""
    saturation_specific_humidity(T, p)

Saturation specific humidity over liquid water using a Tetens-type vapor
pressure relation. Inputs: T in K, p in Pa.
"""
function saturation_specific_humidity(T::Float64, p::Float64)
    es = 611.2 * exp((17.67 * (T - 273.15)) / (T - 29.65))
    return 0.622 * es / max(p - 0.378 * es, 1.0)
end

"""
    bulk_transfer_coefficients(model)

Return neutral bulk transfer coefficients (Cm, Ch) derived from tower
reference heights and slab roughness lengths.
"""
function bulk_transfer_coefficients(model::SCMModel)
    z_u = max(model.tower.z_u, model.slab_params.z0m * 1.1)
    z_t = max(model.tower.z_t, model.slab_params.z0h * 1.1)
    cm = (KAPPA_VK / log(z_u / model.slab_params.z0m))^2
    ch = (KAPPA_VK^2) / (log(z_u / model.slab_params.z0m) * log(z_t / model.slab_params.z0h))
    return cm, ch
end

"""
    resolve_surface_fluxes!(model, forcing)

Compute sensible/latent heat fluxes and slab radiative-conductive terms from
tower-height atmospheric forcing and slab state. If
`forcing.prescribed_surface_fluxes` is true, provided fluxes are used directly.
"""
function resolve_surface_fluxes!(model::SCMModel, forcing::Forcing)
    if forcing.prescribed_surface_fluxes
        model.surface.sensible_flux = forcing.surface_heat_flux
        model.surface.latent_flux = forcing.surface_moisture_flux
        model.surface.net_radiation = 0.0
        model.surface.ground_flux = 0.0
        return forcing.surface_heat_flux, forcing.surface_moisture_flux, 0.0, 0.0
    end

    cfg = model.config
    sp = model.slab_params
    slab = model.slab
    T_air = isfinite(forcing.air_temperature_ref) ? forcing.air_temperature_ref : model.state.theta[1]
    q_air = isfinite(forcing.specific_humidity_ref) ? forcing.specific_humidity_ref : model.state.q[1]
    U_ref = isfinite(forcing.wind_speed_ref) ? forcing.wind_speed_ref : hypot(model.state.u[1], model.state.v[1])
    U_ref = max(U_ref, 0.25)
    p_sfc = max(forcing.surface_pressure, 5.0e4)
    sw = max(forcing.sw_down, 0.0)
    lw = max(forcing.lw_down, 0.0)

    _, ch = bulk_transfer_coefficients(model)
    q_sat = saturation_specific_humidity(slab.skin_temperature, p_sfc)

    H = cfg.rho_air * cfg.cp_air * ch * U_ref * (slab.skin_temperature - T_air)
    LE = cfg.rho_air * LV_WATER * ch * U_ref * sp.moisture_availability * (q_sat - q_air)
    Rn = (1.0 - sp.albedo) * sw + sp.emissivity * (lw - SIGMA_SB * slab.skin_temperature^4)
    G = sp.k * (slab.deep_temperature - slab.skin_temperature) / max(sp.depth, 0.05)

    model.surface.sensible_flux = H
    model.surface.latent_flux = LE
    model.surface.net_radiation = Rn
    model.surface.ground_flux = G
    return H, LE, Rn, G
end

function diffusivities(::AbstractClosure, model::SCMModel)
    error("diffusivities is not implemented for this closure")
end

"""
    diffusivities(closure, model)

Return momentum and scalar diffusivity profiles for the supplied closure.
Closures should implement this method rather than modifying the timestepper.
"""
function diffusivities(closure::ConstantDiffusivityClosure, model::SCMModel)
    nz = model.grid.nz
    km = fill(closure.km, nz)
    kh = fill(closure.kh, nz)
    return km, kh
end

"""
    diffusivities(closure::RiBasedClosure, model)

Diagnose local gradient Richardson number Ri = N²/S² at each interior level
using centered differences, then return stability-reduced diffusivity profiles.

The Brunt-Väisälä frequency and wind shear use the physical z-spacing directly.
Future upgrades that work in the stretched η-coordinate would divide gradients
by grid.jacobian[k] at each level to recover the correct metric.
"""
function diffusivities(closure::RiBasedClosure, model::SCMModel)
    grid  = model.grid
    state = model.state
    nz    = grid.nz
    g     = 9.81
    theta_ref = sum(state.theta) / nz   # reference θ for N² normalization

    km = fill(1.0e-4, nz)   # background minimum (prevents total decoupling)
    kh = fill(1.0e-4, nz)

    for k in 2:nz-1
        dz_span = grid.z[k + 1] - grid.z[k - 1]               # two-cell span
        d_theta = state.theta[k + 1] - state.theta[k - 1]
        du      = state.u[k + 1]     - state.u[k - 1]
        dv      = state.v[k + 1]     - state.v[k - 1]
        N2      = (g / theta_ref) * d_theta / dz_span          # N²
        S2      = (du^2 + dv^2) / dz_span^2                    # S²
        Ri      = S2 > 1.0e-9 ? N2 / S2 : (N2 >= 0.0 ? 1.0e6 : -1.0e6)
        if Ri < closure.Ri_c
            f      = max(0.0, 1.0 - Ri / closure.Ri_c)^2
            km[k]  = closure.K0 * f
            kh[k]  = closure.K0 * f / closure.Pr_t
        end
    end
    return km, kh
end

"""
    diffusivities(closure::CurvatureRiClosure, model)

Diagnose Ri and second-order curvature of θ at each interior level, then
return diffusivity profiles modified by the curvature correction factor fc.

Algorithm
─────────
1. Central-difference first derivatives → N², S², Ri_g.
2. Second-order central-difference of θ using the non-uniform physical
   spacing (dz_up, dz_dn), consistent with the stretched-grid geometry.
3. Nondimensional curvature index κ̃ = (d²θ/dz²) · z / max(dθ/dz, ε).
4. fc = clamp(1 + β·κ̃, 0.1, 2.0).
5. K = K0 · f_stability(Ri) · fc.

The Jacobian `grid.jacobian[k]` is available for a future η-coordinate
upgrade; in physical z the spacing ratios already capture the stretching.
"""
function diffusivities(closure::CurvatureRiClosure, model::SCMModel)
    grid  = model.grid
    state = model.state
    nz    = grid.nz
    g     = 9.81
    theta_ref = sum(state.theta) / nz

    km = fill(1.0e-4, nz)   # background minimum prevents full decoupling
    kh = fill(1.0e-4, nz)

    for k in 2:nz-1
        # ── 1. First derivatives (two-cell central difference) ────────────
        dz_span = grid.z[k + 1] - grid.z[k - 1]
        d_theta = state.theta[k + 1] - state.theta[k - 1]
        du      = state.u[k + 1]     - state.u[k - 1]
        dv      = state.v[k + 1]     - state.v[k - 1]
        N2      = (g / theta_ref) * d_theta / dz_span
        S2      = (du^2 + dv^2) / dz_span^2
        Ri_g    = S2 > 1.0e-9 ? N2 / S2 : (N2 >= 0.0 ? 1.0e6 : -1.0e6)

        # ── 2. Second derivative of θ on non-uniform spacing ─────────────
        dz_up   = grid.z[k + 1] - grid.z[k]
        dz_dn   = grid.z[k]     - grid.z[k - 1]
        d2_theta = ((state.theta[k + 1] - state.theta[k]) / dz_up -
                    (state.theta[k]     - state.theta[k - 1]) / dz_dn) / grid.dz[k]

        # ── 3. Nondimensional curvature index κ̃ ──────────────────────────
        # Denominator guard: use |dθ/dz| evaluated at the two-cell span;
        # clamp to 1e-4 K/m so κ̃ remains finite in nearly-neutral layers.
        grad_theta_abs = abs(d_theta / dz_span)
        kappa_tilde    = d2_theta * grid.z[k] / max(grad_theta_abs, 1.0e-4)

        # ── 4. Correction factor fc ───────────────────────────────────────
        # fc > 1: curvature sharpening the inversion → allow more mixing
        # fc < 1: curvature flattening the gradient   → reduce mixing
        fc = max(0.1, min(2.0, 1.0 + closure.beta * kappa_tilde))

        # ── 5. Stability function + correction ───────────────────────────
        if Ri_g < closure.Ri_c
            f_stability = max(0.0, 1.0 - Ri_g / closure.Ri_c)^2
            km[k] = closure.K0 * f_stability * fc
            kh[k] = closure.K0 * f_stability * fc / closure.Pr_t
        end
    end
    return km, kh
end

"""
    diffusion_tendency(grid, field, kappa)

Compute the vertical diffusion tendency of a single prognostic variable using a
flux-divergence form on the staggered interfaces implied by the cell centers.
Top and bottom diffusive fluxes are set to zero here; explicit surface fluxes
are applied separately in `apply_surface_fluxes!`.
"""
function diffusion_tendency(grid::Grid, field::Vector{Float64}, kappa::Vector{Float64})
    tendency = zeros(grid.nz)
    flux = zeros(grid.nz + 1)
    for k in 2:grid.nz
        k_interface = 0.5 * (kappa[k] + kappa[k - 1])
        gradient = (field[k] - field[k - 1]) / (grid.z[k] - grid.z[k - 1])
        flux[k] = -k_interface * gradient
    end
    flux[1] = 0.0
    flux[end] = 0.0
    for k in 1:grid.nz
        tendency[k] = -(flux[k + 1] - flux[k]) / grid.dz[k]
    end
    return tendency
end

"""
    implicit_diffusion_step!(field, grid, kappa, dt)

Solve the backward-Euler (implicit) vertical diffusion equation in-place:

    (I − dt·L) ϕ_new = ϕ_old

where L is the second-order finite-difference diffusion operator with
zero-flux upper and lower boundaries.  The tridiagonal system is assembled
and solved in O(N) operations using `LinearAlgebra.Tridiagonal`.

Why this matters for the polar ABL
─────────────────────────────────
- The explicit CFL limit Δt ≤ Δz²/(2K) can require Δt < 1 s when dz ~ 1 m.
- K can burst intermittently; an explicit solver will crash at those moments.
- This solver is unconditionally stable: Δt is limited only by the advective
  CFL and the accuracy requirements of the forcing, not by diffusion stiffness.

Note: surface fluxes are still applied explicitly via `apply_surface_fluxes!`
after this call, which is consistent with a first-order operator-split approach.
"""
function implicit_diffusion_step!(field::Vector{Float64},
                                  grid::Grid,
                                  kappa::Vector{Float64},
                                  dt::Float64)
    nz = grid.nz
    a  = zeros(nz)   # sub-diagonal
    b  = zeros(nz)   # main diagonal
    c  = zeros(nz)   # super-diagonal

    # Interior levels
    for k in 2:nz-1
        dz_up = grid.z[k + 1] - grid.z[k]
        dz_dn = grid.z[k]     - grid.z[k - 1]
        K_up  = 0.5 * (kappa[k + 1] + kappa[k])
        K_dn  = 0.5 * (kappa[k]     + kappa[k - 1])
        L_up  = K_up / (dz_up * grid.dz[k])
        L_dn  = K_dn / (dz_dn * grid.dz[k])
        a[k]  = -dt * L_dn
        c[k]  = -dt * L_up
        b[k]  =  1.0 + dt * (L_up + L_dn)
    end

    # Bottom (k=1): zero-flux lower boundary; surface flux handled outside
    K_up  = 0.5 * (kappa[2] + kappa[1])
    L_up  = K_up / ((grid.z[2] - grid.z[1]) * grid.dz[1])
    a[1]  =  0.0
    c[1]  = -dt * L_up
    b[1]  =  1.0 + dt * L_up

    # Top (k=nz): zero-flux upper boundary
    K_dn  = 0.5 * (kappa[nz] + kappa[nz - 1])
    L_dn  = K_dn / ((grid.z[nz] - grid.z[nz - 1]) * grid.dz[nz])
    a[nz] = -dt * L_dn
    c[nz] =  0.0
    b[nz] =  1.0 + dt * L_dn

    A = Tridiagonal(a[2:nz], b, c[1:nz-1])
    field .= A \ field
    return nothing
end

"""
    compute_tke_budget(model, km, kh) → (Ps, B)

Compute column profiles of TKE shear production and buoyant destruction:

    P_s(z) = K_m(z) · S²(z)               (m²/s³)
    B(z)   = (g/θ₀) · K_h(z) · N²(z)      (m²/s³)

In near-equilibrium turbulence P_s ≈ B + ε (dissipation).  When a closure
excessively damps K, P_s drops while B stays elevated — this is the signature
of runaway cooling in stable polar boundary layers.

Monitoring P_s/B versus Ri gives a direct view of where the closure violates
the expected TKE budget balance.
"""
function compute_tke_budget(model::SCMModel, km::Vector{Float64}, kh::Vector{Float64})
    grid  = model.grid
    state = model.state
    nz    = grid.nz
    g     = 9.81
    theta_ref = sum(state.theta) / nz
    Ps = zeros(nz)
    B  = zeros(nz)
    for k in 2:nz-1
        dz_span = grid.z[k + 1] - grid.z[k - 1]
        d_theta = state.theta[k + 1] - state.theta[k - 1]
        du      = state.u[k + 1]     - state.u[k - 1]
        dv      = state.v[k + 1]     - state.v[k - 1]
        N2      = (g / theta_ref) * d_theta / dz_span
        S2      = (du^2 + dv^2) / dz_span^2
        Ps[k]   = km[k] * S2
        B[k]    = kh[k] * max(0.0, N2)   # only stable stratification destroys TKE
    end
    return Ps, B
end

"""
    apply_surface_fluxes!(model, forcing)

Apply prescribed surface heat and moisture fluxes to the first model layer and
update the diagnostic surface state.
"""
function apply_surface_fluxes!(model::SCMModel, forcing::Forcing)
    cfg = model.config
    dz1 = model.grid.dz[1]
    H, LE, Rn, G = resolve_surface_fluxes!(model, forcing)

    # Flux sign convention: positive H/LE = upward from surface into air.
    model.state.theta[1] += cfg.dt * H / (cfg.rho_air * cfg.cp_air * dz1)
    model.state.q[1] += cfg.dt * LE / max(cfg.rho_air * LV_WATER * dz1, 1e-12)

    # Slab prognostic update from residual surface energy.
    sp = model.slab_params
    slab_capacity = max(sp.rho * sp.cp * sp.depth, 1.0)
    dTs = cfg.dt * (Rn - H - LE + G) / slab_capacity
    model.slab.skin_temperature += dTs
    model.slab.deep_temperature += cfg.dt * (model.slab.skin_temperature - model.slab.deep_temperature) / (24.0 * 3600.0)

    if sp.material == :water
        model.slab.liquid_fraction = 1.0
    elseif sp.material == :ice
        model.slab.liquid_fraction = 0.0
    else
        # Smooth transition around freezing for mixed-phase surfaces.
        model.slab.liquid_fraction = clamp01((model.slab.skin_temperature - 271.15) / 4.0)
    end

    model.surface.temperature = model.slab.skin_temperature
end

"""
    step!(model, closure, forcing)

Advance the prognostic state by one timestep.  Solver sequence:

1.  Obtain diffusivity profiles from the closure.
2.  Add large-scale (non-stiff) external tendencies explicitly.
3a. Explicit path (`use_implicit=false`): forward-Euler diffusion.  Subject to
    the CFL constraint Δt ≤ Δz²/(2K).  Suitable for coarse grids and small K.
3b. Implicit path (`use_implicit=true`): backward-Euler tridiagonal solve.
    Unconditionally stable — required when dz < 5 m or K bursts intermittently.
4.  Apply surface fluxes to level 1 (operator-split, always explicit).
5.  Update surface diagnostics and model time.
"""
function step!(model::SCMModel, closure::AbstractClosure, forcing::Forcing)
    km, kh = diffusivities(closure, model)
    dt     = model.config.dt

    if model.config.use_implicit
        # Add non-stiff tendencies first, then solve diffusion implicitly
        model.state.theta .+= dt .* forcing.theta_tendency
        model.state.q     .+= dt .* forcing.q_tendency
        model.state.u     .+= dt .* forcing.u_tendency
        model.state.v     .+= dt .* forcing.v_tendency
        implicit_diffusion_step!(model.state.theta, model.grid, kh, dt)
        implicit_diffusion_step!(model.state.q,     model.grid, kh, dt)
        implicit_diffusion_step!(model.state.u,     model.grid, km, dt)
        implicit_diffusion_step!(model.state.v,     model.grid, km, dt)
    else
        theta_mix = diffusion_tendency(model.grid, model.state.theta, kh)
        q_mix     = diffusion_tendency(model.grid, model.state.q,     kh)
        u_mix     = diffusion_tendency(model.grid, model.state.u,     km)
        v_mix     = diffusion_tendency(model.grid, model.state.v,     km)
        model.state.theta .+= dt .* (theta_mix .+ forcing.theta_tendency)
        model.state.q     .+= dt .* (q_mix     .+ forcing.q_tendency)
        model.state.u     .+= dt .* (u_mix     .+ forcing.u_tendency)
        model.state.v     .+= dt .* (v_mix     .+ forcing.v_tendency)
    end

    apply_surface_fluxes!(model, forcing)
    model.surface.temperature = model.slab.skin_temperature
    model.time += dt
    return nothing
end

"""
    boundary_layer_proxy(model)

Return a crude boundary-layer depth estimate defined as the first level where
potential temperature exceeds the first-layer value by 0.5 K. This is only a
diagnostic placeholder and should be replaced for scientific use.
"""
function boundary_layer_proxy(model::SCMModel)
    theta0 = model.state.theta[1]
    for i in 2:model.grid.nz
        if model.state.theta[i] - theta0 > 0.5
            return model.grid.z[i]
        end
    end
    return model.grid.z[end]
end

"""
    initialize_history(n_steps)

Preallocate the full diagnostic time-series bundle, including TKE budget terms.
"""
function initialize_history(n_steps::Int)
    return SimulationHistory(
        zeros(n_steps), zeros(n_steps), zeros(n_steps),
        zeros(n_steps), zeros(n_steps),
        zeros(n_steps), zeros(n_steps))  # P_s_max, B_max
end

"""
    record!(history, index, model, km, kh)

Save selected model diagnostics including TKE shear production and buoyant
destruction for later plotting or notebook display.
"""
function record!(history::SimulationHistory, index::Int,
                 model::SCMModel, km::Vector{Float64}, kh::Vector{Float64})
    history.time[index]               = model.time
    history.theta_surface[index]      = model.surface.temperature
    history.theta_first_level[index]  = model.state.theta[1]
    history.q_first_level[index]      = model.state.q[1]
    history.boundary_layer_proxy[index] = boundary_layer_proxy(model)
    Ps, B = compute_tke_budget(model, km, kh)
    history.max_shear_production[index]    = maximum(Ps)
    history.max_buoyant_destruction[index] = maximum(B)
    return nothing
end

"""
    run_model([config]; closure=..., forcing_function=..., log_every=0)

High-level driver that initializes the model, advances it to `t_end`, and
returns both the final state and a lightweight history object. The forcing is
supplied through a callback so experiments can swap forcing scenarios without
touching the solver core.
"""
function run_model(config::ModelConfig=ModelConfig();
                   closure::AbstractClosure=ConstantDiffusivityClosure(0.8, 0.5),
                   forcing_function::Function=default_forcing,
                   log_every::Int=0)
    model   = initialize_model(config)
    n_steps = Int(floor(config.t_end / config.dt)) + 1
    history = initialize_history(n_steps)

    km0, kh0 = diffusivities(closure, model)
    record!(history, 1, model, km0, kh0)

    for step_index in 2:n_steps
        forcing = forcing_function(model)
        step!(model, closure, forcing)
        km, kh = diffusivities(closure, model)
        record!(history, step_index, model, km, kh)
        if log_every > 0 && (step_index - 1) % log_every == 0
            @printf("t = %6.2f h  θ_sfc = %6.2f K  h = %6.1f m  P_s = %.2e  B = %.2e\n",
                    model.time / 3600.0,
                    model.surface.temperature,
                    history.boundary_layer_proxy[step_index],
                    history.max_shear_production[step_index],
                    history.max_buoyant_destruction[step_index])
        end
    end

    return model, history
end

"""
    example_run()

Run a short demonstration case suitable for a first smoke test or Pluto/Jupyter
example.
"""
function example_run()
    config  = ModelConfig(z_top=500.0, nz=30, dt=20.0, t_end=6.0 * 3600.0,
                          use_implicit=false)
    closure = ConstantDiffusivityClosure(0.6, 0.4)
    return run_model(config; closure=closure, log_every=60)
end

end