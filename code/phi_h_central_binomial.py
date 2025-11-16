from scipy.special import comb
import numpy as np

def phi_h_central_binomial(zeta, beta_h=9, N_exact=10, N_asymp=5):
    """
    Compute φ_h(ζ) = (1 - β_h ζ)^(-1/2) using:
    - First N_exact terms: exact central binomials
    - Next N_asymp terms: Stirling approximation
    """
    x = beta_h * zeta / 4
    phi = 0.0

    # Exact terms
    for n in range(N_exact):
        phi += comb(2*n, n, exact=True) * (x**n)

    # Asymptotic tail
    for n in range(N_exact, N_exact + N_asymp):
        stirling_approx = (4**n) / np.sqrt(np.pi * n)
        phi += stirling_approx * (x**n)

    return phi

# Test
zeta_test = -0.5  # Moderately unstable
phi_exact = (1 - 9*zeta_test)**(-0.5)
phi_series = phi_h_central_binomial(zeta_test, beta_h=9, N_exact=10, N_asymp=5)
error = abs(phi_series - phi_exact) / phi_exact
print(f"Exact: {phi_exact:.8f}, Series: {phi_series:.8f}, Error: {error:.2e}")