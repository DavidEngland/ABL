from spectral_scalar_closure import SpectralScalarClosure

# Initialize with momentum coefficients
model = SpectralScalarClosure(coeffs=my_momentum_coeffs, lambda_=0.25)

# Fit n_c and s0 for methane observations (ratio should be phi_CH4 / phi_m)
fit_results = model.fit_tracer_filter(zeta_obs, methane_ratio_obs)

print(f"Neutral ratio: {fit_results['s0']}")
print(f"Filter scale (n_c): {fit_results['n_c']}")

# Predict methane transfer ratio
methane_ratio_model = model.tracer_ratio(
    zeta_grid,
    s0=fit_results["s0"],
    n_c=fit_results["n_c"],
)