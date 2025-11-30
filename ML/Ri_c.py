import numpy as np

EPS = 1e-6
RI_MIN, RI_MAX = 0.15, 2.0

def safe_div(x, y, eps=EPS):
    return x / np.clip(y, -eps, eps)

def clamp(x, lo=RI_MIN, hi=RI_MAX):
    return np.minimum(np.maximum(x, lo), hi)

def ri_c_star(state, coeff):
    # Example symbolic form (candidate)
    Gamma = np.maximum(state["dtheta_dz"], 0.0)           # K/m
    S = np.maximum(state["dU_dz_abs"], 0.0)               # 1/s
    zeta = state["z_over_L"]                              # unitless
    dz = state["delta_z"]                                 # m

    term1 = coeff.a0 + coeff.a1 * np.log1p(Gamma) + coeff.a2 * np.log1p(S)
    term2 = coeff.a3 * np.tanh(coeff.a4 * zeta) + coeff.a5 * safe_div(Gamma, S + EPS)
    raw = term1 + term2 + coeff.a6 * np.tanh(coeff.a7 * dz / (1.0 + dz))
    return clamp(raw)

def monotonicity_violations(state_grid, coeff):
    # Finite-difference checks on Gamma and S
    vio = {"Gamma_pos": 0, "S_neg": 0}
    for Gamma in state_grid["Gamma"]:
        for S in state_grid["S"]:
            s = dict(dtheta_dz=Gamma, dU_dz_abs=S, z_over_L=state_grid["zeta"], delta_z=state_grid["dz"])
            base = ri_c_star(s, coeff)
            s["dtheta_dz"] = Gamma * 1.05
            dGamma = ri_c_star(s, coeff) - base
            s["dtheta_dz"] = Gamma
            s["dU_dz_abs"] = S * 1.05
            dS = ri_c_star(s, coeff) - base
            vio["Gamma_pos"] += int(dGamma < -1e-4)   # should not decrease with Gamma (stable)
            vio["S_neg"]    += int(dS     >  1e-4)   # should not increase with shear
    return vio
