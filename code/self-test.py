# -----------------------------------------------------------------------------
# Optional self-test block (does not execute on import)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Simple diagnostic: evaluate BD_PL φ and Ri↔ζ transforms at small ζ
    pars = {"am": 0.25, "bm": 5.0, "ah": 0.5, "bh": 7.0}
    phi_m, phi_h = make_profile("BD_PL", pars)
    z = 0.05
    Ri = ri_from_zeta(z, phi_m, phi_h)
    print("ζ =", z, "Ri =", Ri)

    # Check Newton inversion
    z0 = zeta_from_ri_series(Ri, Delta=pars["ah"]*pars["bh"] - 2*pars["am"]*pars["bm"],
                             c1=pars["ah"]*(pars["bh"]**2) - 2*pars["am"]*(pars["bm"]**2))
    zN = zeta_from_ri_newton(Ri, phi_m, phi_h, z0)
    print("Newton ζ(Ri) =", zN)