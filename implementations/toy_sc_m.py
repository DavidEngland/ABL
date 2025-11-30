# toy_sc_m.py
# Minimal Arctic single-column toy model (Python/Numpy)
# Direct translation of toy_sc_m.jl

import numpy as np

# ----- Physical constants -----
g = 9.80665
cp = 1004.0
Rd = 287.05
Rv = 461.5
Lv = 2.5e6
sigma = 5.670374419e-8
rho_air = 1.225
CH = 1.2e-3
CE = 1.0e-3


# ----- Grid -----
class Grid:
    def __init__(self, z, dz):
        self.z = z
        self.dz = dz
        self.nz = len(z)


def create_stretched_grid(z_top, nz, dz_min=5.0):
    η = np.linspace(0.0, 1.0, nz)
    stretch = (np.exp(3 * η) - 1.0) / (np.exp(3.0) - 1.0)
    z = z_top * stretch

    dz = np.zeros_like(z)
    dz[0] = max(dz_min, z[1] - z[0])

    for i in range(1, nz - 1):
        dz[i] = z[i + 1] - z[i - 1]

    dz[-1] = z[-1] - z[-2]
    return Grid(z, dz)


# ----- State -----
class State:
    def __init__(self, T, q, Ts, ice_h):
        self.T = T
        self.q = q
        self.Ts = Ts
        self.ice_h = ice_h


def initialize_state(grid):
    nz = grid.nz
    T = np.zeros(nz)
    q = np.zeros(nz)

    for i in range(nz):
        T[i] = 250.0 + 0.015 * grid.z[i]
        q[i] = 3e-4 * np.exp(-grid.z[i] / 2000.0)

    Ts = 258.0
    ice_h = 1.0
    return State(T, q, Ts, ice_h)


# ----- Albedo -----
def albedo(ice_h, snow_mass=0.0, meltpond_frac=0.0):
    alb_ice = 0.6
    alb_snow = 0.85
    alb_pond = 0.2

    if snow_mass > 0:
        return alb_snow * (1 - meltpond_frac) + alb_pond * meltpond_frac
    else:
        frac = np.clip(ice_h / 0.3, 0, 1)
        return alb_pond * (1 - frac) + alb_ice * frac


# ----- Radiation -----
def radiation_simple(Ts, alb, sw_down, lw_down):
    sw_abs = sw_down * (1.0 - alb)
    lw_up = sigma * Ts**4
    lw_net = lw_down - lw_up
    rnet = sw_abs + lw_net
    return rnet, sw_abs, lw_net


# ----- Surface fluxes -----
def surface_fluxes(Ts, T1, q1, U, rho=rho_air):
    H = rho * cp * CH * U * (Ts - T1)
    LE = rho * Lv * CE * U * (0.0 - q1)
    return H, LE


# ----- Saturation q -----
def q_sat(T, p=900.0):
    es = 6.112 * np.exp((17.67 * (T - 273.15)) / (T - 29.65)) * 100.0
    return 0.622 * es / (p * 100.0 - 0.378 * es)


# ----- K-profile -----
def compute_K_profile(grid, state, K0=1.0, hbl=200.0, minK=1e-5):
    K = K0 * np.exp(-grid.z / hbl)

    if state.Ts < state.T[0]:
        inv = state.T[0] - state.Ts
        factor = max(0.05, 1.0 - 0.2 * inv)
        K *= factor

    return np.maximum(K, minK)


# ----- Vertical diffusion -----
def vertical_diffusion_tend(grid, var, K):
    nz = grid.nz
    tend = np.zeros(nz)
    F = np.zeros(nz + 1)

    for k in range(1, nz):
        Km = 0.5 * (K[k] + K[k - 1])
        dvar = (var[k] - var[k - 1]) / (grid.z[k] - grid.z[k - 1])
        F[k] = -Km * dvar

    F[0] = 0.0
    F[-1] = 0.0

    for i in range(nz):
        tend[i] = -(F[i + 1] - F[i]) / grid.dz[i]

    return tend


# ----- Step -----
def step(grid, state, dt, forcings):
    alb = albedo(state.ice_h)
    rnet, sw_abs, lw_net = radiation_simple(
        state.Ts, alb, forcings["sw_down"], forcings["lw_down"]
    )

    H, LE = surface_fluxes(state.Ts, state.T[0], state.q[0], forcings["U"])

    K = compute_K_profile(grid, state, K0=0.5, hbl=100.0, minK=1e-6)

    tendT = vertical_diffusion_tend(grid, state.T, K)
    tendq = vertical_diffusion_tend(grid, state.q, K)

    adv_T = forcings.get("adv_T")
    adv_q = forcings.get("adv_q")

    for i in range(grid.nz):
        state.T[i] += dt * (tendT[i] + (0.0 if adv_T is None else adv_T[i]))
        state.q[i] += dt * (tendq[i] + (0.0 if adv_q is None else adv_q[i]))

    k_ice = 2.1
    T_deep = 258.0
    G = k_ice * (state.Ts - T_deep) / max(state.ice_h, 0.01)

    dTs = (rnet - (H + LE) - G) * dt / (2100.0 * state.ice_h)
    state.Ts += dTs

    dz1 = grid.dz[0]
    state.T[0] += dt * (H / (rho_air * cp * dz1))
    state.q[0] += dt * (LE / (rho_air * Lv * dz1))

    for i in range(grid.nz):
        qsat = q_sat(state.T[i])
        if state.q[i] > qsat:
            dq = state.q[i] - qsat
            dTlat = -(Lv * dq) / cp
            state.T[i] += dTlat
            state.q[i] = qsat

    if state.Ts > 273.15:
        melt_rate = 1e-6 * (state.Ts - 273.15)
        state.ice_h = max(0.0, state.ice_h - melt_rate * dt)


# ----- Driver -----
def run_experiment(ztop=4000.0, nz=40, dt=60.0, tmax=24 * 3600.0):
    grid = create_stretched_grid(ztop, nz)
    state = initialize_state(grid)

    sw_max = 200.0
    lw_down = 250.0
    U = 5.0

    times = np.arange(0.0, tmax + dt, dt)

    Ts_hist = np.zeros_like(times)
    T1_hist = np.zeros_like(times)
    ice_h_hist = np.zeros_like(times)

    for n, t in enumerate(times):
        dayfrac = (t / 86400.0) % 1.0
        sw_down = max(0.0, sw_max * np.sin(2 * np.pi * dayfrac))

        forc = dict(
            sw_down=sw_down,
            lw_down=lw_down,
            U=U,
            adv_T=None,
            adv_q=None,
        )

        step(grid, state, dt, forc)

        Ts_hist[n] = state.Ts
        T1_hist[n] = state.T[0]
        ice_h_hist[n] = state.ice_h

        if n % 240 == 0:
            print(f"t={t/3600:.1f} h  Ts={state.Ts:.2f}  T1={state.T[0]:.2f}  ice_h={state.ice_h:.3f}  sw={sw_down:.1f}")

    return times, Ts_hist, T1_hist, ice_h_hist


if __name__ == "__main__":
    times, Ts_hist, T1_hist, ice_h_hist = run_experiment()
    print("Run complete.")