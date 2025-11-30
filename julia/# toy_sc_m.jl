# toy_sc_m.jl
# Minimal Arctic single-column toy model (Julia)
# Prognostic: temperature (K), specific humidity (kg/kg), surface ice slab temp & thickness
# Processes: vertical diffusion (K-profile), simple shortwave/albedo, bulk surface fluxes,
#           simple longwave cooling (Stefan-Boltzmann), saturation adjustment for condensation.
#
# This is intentionally simple and explicit — reduce dt for stability or replace diffusion step
# with an implicit solver for production use.

using Printf
using LinearAlgebra
using Random

# ----- Physical constants -----
const g = 9.80665              # m/s^2
const cp = 1004.0              # J/(kg K)
const Rd = 287.05              # J/(kg K)
const Rv = 461.5               # J/(kg K)
const Lv = 2.5e6               # J/kg (latent heat of vaporization, approx)
const sigma = 5.670374419e-8   # Stefan-Boltzmann, W/m2/K4
const rho_air = 1.225          # reference air density kg/m3 (near surface)
const CH = 1.2e-3              # bulk transfer coeff for sensible heat (rough)
const CE = 1.0e-3              # bulk transfer coeff for latent heat

# ----- Grid & state -----
struct Grid
    z::Vector{Float64}    # center heights (m)
    dz::Vector{Float64}   # layer thickness (m)
    nz::Int
end

function create_stretched_grid(z_top::Float64, nz::Int; dz_min=5.0)
    # simple exponential stretching; fine near surface
    η = range(0.0, 1.0, length=nz)
    stretch = (exp.(3*η) .- 1.0) ./ (exp(3.0)-1.0)  # maps 0..1 to 0..1, clustered near 0
    z = z_top .* stretch
    # compute layer thickness as difference of cell centers (approx)
    dz = similar(z)
    dz[1] = max(dz_min, z[2]-z[1])
    for i in 2:nz-1
        dz[i] = z[i+1]-z[i-1]
    end
    dz[end] = z[end]-z[end-1]
    return Grid(z, dz, nz)
end

# Model state
mutable struct State
    T::Vector{Float64}      # potential temperature or T (K) (we'll use T for simplicity)
    q::Vector{Float64}      # specific humidity (kg/kg)
    Ts::Float64             # surface (ice) skin temperature (K)
    ice_h::Float64          # ice thickness (m)
end

# ----- Initial conditions -----
function initialize_state(grid::Grid)
    nz = grid.nz
    T = zeros(nz)
    q = zeros(nz)
    # simple Arctic-like profiles
    for i in 1:nz
        T[i] = 250.0 + 0.015 * grid.z[i]   # cold near-surface ~250 K, weak increase with height
        q[i] = 3e-4 * exp(-grid.z[i]/2000) # very dry
    end
    Ts = 258.0       # initial ice skin temp (K) ~ -15 C
    ice_h = 1.0      # 1 m ice
    return State(T, q, Ts, ice_h)
end

# ----- Albedo parameterization (very simple) -----
function albedo(ice_h::Float64, snow_mass::Float64=0.0, meltpond_frac::Float64=0.0)
    # base albedos:
    alb_ice = 0.6
    alb_snow = 0.85
    alb_pond = 0.2
    # weight by presence; snow dominates if any; meltpond reduces effective albedo
    if snow_mass > 0.0
        return alb_snow*(1.0-meltpond_frac) + alb_pond*meltpond_frac
    else
        # thin ice has slightly lower albedo
        frac = clamp(ice_h/0.3, 0.0, 1.0)  # if ice < 0.3m, darker
        return alb_pond*(1-frac) + alb_ice*frac
    end
end

# ----- Simple radiation: prescribe incoming SW, compute absorbed SW and LW up/down -----
function radiation_simple(Ts::Float64, alb::Float64, sw_down::Float64, lw_down::Float64)
    # sw_down: incoming shortwave at surface (W/m2)
    sw_absorbed = sw_down * (1.0 - alb)
    lw_up = sigma * Ts^4
    lw_net = lw_down - lw_up
    rnet = sw_absorbed + lw_net
    return rnet, sw_absorbed, lw_net
end

# ----- Bulk turbulent fluxes at surface -----
function surface_fluxes(Ts::Float64, T1::Float64, q1::Float64, U::Float64; rho=rho_air)
    # T1, q1 are air at first model level
    H = rho*cp*CH*U*(Ts - T1)      # positive when surface warmer than air (upward sensible flux)
    LE = rho*Lv*CE*U*(0.0 - q1)    # assume saturated surface q_s ~ 0 (ice surface evaporation negligible)
    # Note: for sea-ice sublimation/evap we could compute q_s from Ts, but keep simple
    return H, LE
end

# ----- Saturation mixing ratio (approx Clausius-Clapeyron, bolton formula variant) -----
function q_sat(T::Float64, p::Float64=900.0)  # pressure in hPa approx for Arctic surface low pressure
    # Tetens/Bolton approximation (kg/kg)
    es = 6.112 * exp((17.67*(T-273.15))/(T-29.65)) * 100.0  # Pa
    qsat = 0.622 * es / (p*100.0 - 0.378*es)
    return qsat
end

# ----- Simple K-profile (stability-aware) -----
function compute_K_profile(grid::Grid, state::State; K0=1.0, hbl=200.0, minK=1e-5)
    nz = grid.nz
    K = zeros(nz)
    for i in 1:nz
        z = grid.z[i]
        K[i] = K0 * exp(-z/hbl)  # decays with height
    end
    # adjust K for strong near-surface stability: if Ts < T1 (inversion) reduce K
    if state.Ts < state.T[1]
        # stronger inversion -> reduce K near surface
        inv_strength = state.T[1] - state.Ts
        factor = max(0.05, 1.0 - 0.2*inv_strength)  # clamp
        K .= K .* factor
    end
    K .= max.(K, minK)
    return K
end

# ----- Vertical diffusion tendency (explicit) for T and q -----
function vertical_diffusion_tend(grid::Grid, var::Vector{Float64}, K::Vector{Float64})
    nz = grid.nz
    tend = zeros(nz)
    # fluxes at interfaces using centered differences; need K at interfaces -> average
    # cell centers i=1..nz, interfaces 0..nz
    F = zeros(nz+1)
    # top/bottom boundary fluxes: F[1] at surface, F[end] at top
    for k in 2:nz
        Km = 0.5*(K[k] + K[k-1])
        dvar = (var[k] - var[k-1]) / (grid.z[k] - grid.z[k-1])
        F[k] = -Km * dvar   # downward positive? here flux upward is negative gradient * K
    end
    # boundaries: assume zero flux at top (F[nz+1]=0); surface flux handled separately -> set F[1] later
    F[1] = 0.0
    F[nz+1] = 0.0
    # divergences
    for i in 1:nz
        dz = grid.dz[i]
        tend[i] = -(F[i+1] - F[i]) / dz
    end
    return tend
end

# ----- Time stepping -----
function step!(grid::Grid, state::State, dt::Float64, forcings)
    # forcings: Dict with keys: sw_down, lw_down, U (10m wind), adv_T (vector tendencies), adv_q
    # Compute albedo
    alb = albedo(state.ice_h)
    # radiation
    rnet, sw_abs, lw_net = radiation_simple(state.Ts, alb, forcings[:sw_down], forcings[:lw_down])
    # surface fluxes
    H, LE = surface_fluxes(state.Ts, state.T[1], state.q[1], forcings[:U])
    # turbulence K
    K = compute_K_profile(grid, state; K0=0.5, hbl=100.0, minK=1e-6)
    # vertical diffusion tendencies for T and q
    tendT = vertical_diffusion_tend(grid, state.T, K)
    tendq = vertical_diffusion_tend(grid, state.q, K)
    # apply tendencies + large-scale advective tendencies if present
    for i in 1:grid.nz
        state.T[i] += dt*(tendT[i] + (forcings[:adv_T] === nothing ? 0.0 : forcings[:adv_T][i]))
        state.q[i] += dt*(tendq[i] + (forcings[:adv_q] === nothing ? 0.0 : forcings[:adv_q][i]))
    end
    # surface energy update: simple slab with heat capacity
    # G = ground conduction approximated by conductive flux into ice bottom ~ k_ice*(Tsurface - Tdeep)/h
    k_ice = 2.1      # W/m/K
    T_deep = 258.0   # deep ice temp (fixed)
    G = k_ice*(state.Ts - T_deep) / max(state.ice_h, 0.01)
    # net surface budget: rnet (SW+LW) - (H + LE) - G  -> heat into slab
    dTs = (rnet - (H + LE) - G) * dt / (2100.0 * state.ice_h)  # approximate heat capacity (ρ*c*h) with ρ*c ≈ 2100*J/m3/K
    state.Ts += dTs
    # surface-to-air exchange: modify first model level by H/(rho*cp*dz)
    # convert H (W/m2) to temperature tendency in first layer
    dz1 = grid.dz[1]
    state.T[1] += dt * (H / (rho_air * cp * dz1))
    # humidity tendency from LE (positive upward flux -> reduce near-surface q)
    state.q[1] += dt * (LE / (rho_air * Lv * dz1))  # sign/units simplified
    # simple saturation adjustment: if q > qsat -> condensate removed (latent cooling not fully implemented)
    for i in 1:grid.nz
        qsat = q_sat(state.T[i])
        if state.q[i] > qsat
            # remove excess to cloud (not tracked here), apply latent cooling to T
            dq = state.q[i] - qsat
            dTlat = - (Lv * dq) / cp
            state.T[i] += dTlat
            state.q[i] = qsat
        end
    end
    # optional: evaporation/sublimation from ice -> reduce ice mass if Ts > 273.15 (melt)
    if state.Ts > 273.15
        melt_rate = 1e-6*(state.Ts - 273.15)  # m/s arbitrary scaling
        state.ice_h = max(0.0, state.ice_h - melt_rate*dt)
    end
    return nothing
end

# ----- Driver / experiment -----
function run_experiment(; ztop=4000.0, nz=40, dt=60.0, tmax=24*3600.0)
    grid = create_stretched_grid(ztop, nz)
    state = initialize_state(grid)
    # forcings: simple diurnal-ish SW (peak during "day"), static LW_down
    sw_max = 200.0   # W/m2 at peak (very approximate for high lat in summer)
    lw_down = 250.0  # W/m2 typical downward LW
    U = 5.0          # m/s wind
    times = 0.0:dt:tmax
    nsteps = length(times)
    # recorders
    Ts_hist = zeros(nsteps)
    T1_hist = zeros(nsteps)
    ice_h_hist = zeros(nsteps)
    for (it, t) in enumerate(times)
        # simple daily cycle for sw_down, scaled by seasonality (here we just do a sine)
        dayfrac = (t / 86400.0) % 1.0
        sw_down = max(0.0, sw_max * sin(2π*dayfrac))  # simple day-night cycle
        forcings = Dict(:sw_down=>sw_down, :lw_down=>lw_down, :U=>U, :adv_T=>nothing, :adv_q=>nothing)
        step!(grid, state, dt, forcings)
        Ts_hist[it] = state.Ts
        T1_hist[it] = state.T[1]
        ice_h_hist[it] = state.ice_h
        if it % 240 == 1
            @printf("t=%.1f h  Ts=%.2f K  T1=%.2f K  ice_h=%.3f m  sw=%.1f W/m2\n", t/3600, state.Ts, state.T[1], state.ice_h, sw_down)
        end
    end
    return times, Ts_hist, T1_hist, ice_h_hist
end

# Run a short experiment (1 day)
if abspath(PROGRAM_FILE) == @__FILE__
    times, Ts_hist, T1_hist, ice_h_hist = run_experiment(ztop=4000.0, nz=40, dt=60.0, tmax=24*3600.0)
    println("Run complete.")
end


⸻

Notes, caveats & how to improve
   •   Simplicity vs. realism: this toy model intentionally simplifies microphysics, radiation, and surface processes. It’s meant for quick prototyping of feedbacks (albedo → SW absorption → Ts → melt) and for experimenting with ML closures later.
   •   Stability: the vertical diffusion uses an explicit scheme; for realistic K and dt you’ll need very small dt or an implicit solver (Crank–Nicolson) for diffusion. For production runs switch to an implicit tridiagonal solve for diffusion tendencies.
   •   Radiation: the 2-stream/longwave is replaced by a very simple Stefan–Boltzmann upwelling + prescribed downward LW. Replace with a two-stream solver (e.g., delta-Eddington) for realistic cloud radiative effects.
   •   Clouds & microphysics: saturation adjustment is very crude (instantaneous condensation). Add prognostic cloud liquid/ice and a microphysics scheme for mixed-phase clouds to capture Arctic cloud feedbacks.
   •   Surface fluxes: LE is approximated assuming a dry ice surface — for open water or melt ponds compute saturation q_s from Ts and use proper bulk formulas.
   •   K-profile: the simple exponential K can be replaced with an EDMF or a stability-dependent K from Monin–Obukhov theory for better representation of very stable Arctic BL.
   •   Units & pressure: here we used simplified constants and assumed near-surface pressure; for accuracy use a pressure profile and convert T→potential temperature if you want to conserve energy properly.

⸻

How to connect this to Clima.jl / TurbulenceConvection.jl
   •   The Clima stack would provide robust vertical discretization, an EDMF closure, radiation modules, and land/sea-ice components. Once you’re happy with the toy SCM behavior:
	1.	Replace the explicit diffusion vertical_diffusion_tend with the TurbulenceConvection EDMF tendencies (it will expose an API to compute tendencies given state & parameters).
	2.	Replace radiation_simple with a Clima radiation module (two-stream) and use ClimaLand for ice thermodynamics.
	3.	Use ClimaCore grids & fields instead of the lightweight Grid/State here — that lets you use existing solvers and time integrators.
   •   I can write a second script that shows an example glue layer calling TurbulenceConvection.jl APIs (if you confirm you have those packages available). I will not invent package APIs — I’ll read the actual package docs and adapt — if you want that, say “Yes — show Clima glue” and I’ll produce a safe adaptation.

⸻

Want me to:
   •   (A) Replace the explicit diffusion with an implicit Crank–Nicolson diffusion solver (safer dt), or
   •   (B) Expand the radiation to a two-stream shortwave + longwave scheme with cloud optical depth, or
   •   (C) Produce the Clima.jl/TurbulenceConvection.jl glue version (I’ll adapt to real package APIs if you confirm you have those libs installed), or
   •   (D) Translate this toy SCM to a Python script (xarray/numpy) so it plugs into your Python ML toolchain?

Pick one (A/B/C/D) and I’ll produce it now.