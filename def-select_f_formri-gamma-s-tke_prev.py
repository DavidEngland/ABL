def select_f_form(Ri, Gamma, S, TKE_prev):
    # Compute dynamic Ri_c*
    Ric_star = dynamic_Ric(Gamma, S, TKE_prev)

    if Ri < 0.7 * Ric_star:
        # MOST regime: use series + Newton inversion
        return f_from_most_zeta(Ri, Delta, c1)
    elif Ri < 1.3 * Ric_star:
        # Blend zone: average Padé [1/1] and exponential
        f_pade = pade_11(Ri, a, b)
        f_exp = exp_decay(Ri, gamma, Ric_star)
        chi = blend_weight(Ri, Ric_star)
        return (1 - chi) * f_pade + chi * f_exp
    else:
        # Collapsed regime: exponential decay only
        return exp_decay(Ri, gamma, Ric_star)