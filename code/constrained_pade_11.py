def constrained_pade_11(Ri, b, am=6.5, bm=80, Delta=-6):
    # Fix a from neutral slope constraint
    a = am + b
    # Verify second-order match
    expected_c2 = bm - am * Delta
    actual_c2 = b**2 - a*b
    print(f"Second-order error: {abs(actual_c2 - expected_c2)}")
    return (1 + a*Ri) / (1 + b*Ri)

# Fit with constraint
popt_constrained, _ = curve_fit(
    lambda Ri, b: constrained_pade_11(Ri, b, am=6.5, bm=80, Delta=-6),
    Ri_data, fm_data
)