# Homework Problem: Central Binomial Series for Half-Integer MOST Profiles

**Course:** Advanced Boundary Layer Meteorology / Mathematical Methods in Atmospheric Science  
**Topic:** Analytic representations of Monin–Obukhov stability functions  
**Difficulty:** Graduate level  
**Prerequisites:** Taylor series, combinatorics, asymptotic analysis  
**Estimated Time:** 3–4 hours

---

## Problem Statement

The Businger–Dyer (1971) heat stability function for **unstable** conditions is
$$
\phi_h(\zeta) = (1 - \beta_h \zeta)^{-1/2}, \quad \zeta < 0
$$
with $\beta_h \approx 9$ (Kansas data; use 16 for cleaner arithmetic below).

Your task: **Derive the exact series expansion** using central binomial coefficients, analyze convergence, and compute the Richardson number curvature to specified precision.

---

## Part A: Series Derivation (40 points)

### A1. Binomial Theorem Review (10 points)

State the **generalized binomial theorem** for $(1 + x)^\alpha$ where $\alpha$ is not a positive integer:
$$
(1 + x)^\alpha = \sum_{n=0}^\infty \binom{\alpha}{n} x^n, \quad |x| < 1
$$

Define the **generalized binomial coefficient**:
$$
\binom{\alpha}{n} = \frac{\alpha(\alpha - 1)\cdots(\alpha - n + 1)}{n!}
$$

**Task:** Write out $\binom{-1/2}{n}$ for $n = 0, 1, 2, 3$ explicitly (reduce to fractions).

### A2. Central Binomial Connection (15 points)

The generating function for **central binomial coefficients** $\binom{2n}{n}$ is
$$
\frac{1}{\sqrt{1 - 4x}} = \sum_{n=0}^\infty \binom{2n}{n} x^n.
$$

**Prove** that for $\phi_h(\zeta) = (1 - \beta_h \zeta)^{-1/2}$:
$$
\phi_h(\zeta) = \sum_{n=0}^\infty \binom{2n}{n} \left(\frac{\beta_h \zeta}{4}\right)^n
$$

**Hint:** Substitute $x = \beta_h \zeta / 4$ into the generating function.

### A3. First Five Terms (15 points)

For $\beta_h = 16$ (convenient since $\beta_h / 4 = 4$), compute:
$$
\phi_h(\zeta) = \sum_{n=0}^4 \binom{2n}{n} (4\zeta)^n + O(\zeta^5)
$$

**Explicit values needed:**
- $\binom{0}{0} = 1$
- $\binom{2}{1} = 2$
- $\binom{4}{2} = 6$
- $\binom{6}{3} = 20$
- $\binom{8}{4} = 70$

Write the polynomial approximation $\phi_h^{(4)}(\zeta)$.

---

## Part B: Convergence and Domain (25 points)

### B1. Radius of Convergence (10 points)

The series $\sum \binom{2n}{n} x^n$ converges for $|x| < 1/4$.

**Show** that for $\phi_h(\zeta) = (1 - 16\zeta)^{-1/2}$, the series converges when
$$
|\zeta| < \frac{1}{16}.
$$

**Physical interpretation:** The power-law has a **pole** at $\zeta = 1/16 = 0.0625$. Beyond this, the profile is undefined. How does this compare to typical **unstable** $\zeta$ ranges in field data (e.g., SHEBA, ARM SGP)?

### B2. Asymptotic Growth Rate (15 points)

The central binomial coefficients grow as
$$
\binom{2n}{n} \sim \frac{4^n}{\sqrt{\pi n}}.
$$

**Verify numerically** for $n = 10, 20, 50$:
$$
R_n = \frac{\binom{2n}{n}}{\frac{4^n}{\sqrt{\pi n}}}
$$

Expected: $R_n \to 1$ as $n \to \infty$.

**Question:** At $\zeta = -0.05$ (moderately unstable), how many terms $N$ are needed for the series to converge within **1% relative error**? Use the truncation error estimate:
$$
\left|\phi_h - \phi_h^{(N)}\right| \lesssim \binom{2(N+1)}{N+1} \left(\frac{\beta_h |\zeta|}{4}\right)^{N+1}
$$

---

## Part C: Richardson Number Curvature (35 points)

### C1. Gradient Richardson Number (10 points)

Recall that in MOST:
$$
Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}
$$

For **momentum**, use (Businger et al. 1971):
$$
\phi_m(\zeta) = (1 - \beta_m \zeta)^{-1/4}, \quad \beta_m = 16
$$

**Derive** the central binomial series for $\phi_m(\zeta)$ using the generating function
$$
\frac{1}{(1 - 4x)^{1/4}} = \sum_{n=0}^\infty \binom{-1/4}{n} (-4x)^n
$$

**Simplify** the first three terms.

### C2. Series for $Ri_g$ (15 points)

Substitute the series for $\phi_m$ and $\phi_h$ into
$$
Ri_g(\zeta) = \zeta \frac{\phi_h}{\phi_m^2}
$$

**Compute** the **first four terms** of the series:
$$
Ri_g(\zeta) = a_1 \zeta + a_2 \zeta^2 + a_3 \zeta^3 + a_4 \zeta^4 + \cdots
$$

**Express** $a_1, a_2, a_3, a_4$ in terms of central binomial coefficients and generalized binomial coefficients.

### C3. Neutral Curvature (10 points)

The **neutral curvature invariant** is
$$
\left.\frac{d^2 Ri_g}{d\zeta^2}\right|_{\zeta=0} = 2\Delta = 2(\alpha_h \beta_h - 2\alpha_m \beta_m)
$$

For $\alpha_h = 1/2, \beta_h = 16, \alpha_m = 1/4, \beta_m = 16$:

**Verify** that
$$
2\Delta = 2\left(\frac{1}{2} \cdot 16 - 2 \cdot \frac{1}{4} \cdot 16\right) = 2(8 - 8) = 0
$$

**Physical interpretation:** The Businger et al. (1971) parameter set has **zero neutral curvature** in the unstable branch. What does this imply for the shape of $Ri_g(\zeta)$ near neutrality?

---

## Part D: Numerical Implementation (Bonus +10 points)

### D1. Python Series Evaluator

Write a function to compute $\phi_h(\zeta)$ using the central binomial series to order $N$:

```python
from scipy.special import comb
import numpy as np

def phi_h_series(zeta, beta_h=16, N=10):
    """
    Compute phi_h(zeta) = (1 - beta_h * zeta)^(-1/2) using central binomial series.
    
    Parameters:
    -----------
    zeta : float or array
        Stability parameter (negative for unstable)
    beta_h : float
        Coefficient (default 16 for clean arithmetic)
    N : int
        Number of terms in series
    
    Returns:
    --------
    phi_h : float or array
        Heat stability function
    """
    x = beta_h * zeta / 4
    phi = sum(comb(2*n, n, exact=True) * (x**n) for n in range(N+1))
    return phi

# Test
zeta_test = np.linspace(-0.06, 0, 100)
phi_exact = (1 - 16*zeta_test)**(-0.5)
phi_series = phi_h_series(zeta_test, N=10)

# Plot comparison
import matplotlib.pyplot as plt
plt.plot(zeta_test, phi_exact, 'k-', label='Exact (power-law)')
plt.plot(zeta_test, phi_series, 'r--', label=f'Series (N=10)')
plt.xlabel('ζ'); plt.ylabel('φ_h'); plt.legend(); plt.grid()
plt.title('Central Binomial Series vs Exact')
plt.show()
```

**Task:** Generate the plot for $N = 5, 10, 20$. Report the maximum relative error over $\zeta \in [-0.05, 0]$.

### D2. Curvature Computation

Modify the code to compute
$$
\frac{d^2 Ri_g}{d\zeta^2}\bigg|_{\zeta}
$$
using the series derivatives. Compare to the analytic formula from the lecture notes:
$$
\frac{d^2 Ri_g}{d\zeta^2} = F\left[2V_{\log} + \zeta(V_{\log}^2 - W_{\log})\right]
$$

**Plot** curvature vs $\zeta$ for both methods.

---

## Part E: Computational Efficiency (Bonus +15%)

### E1. Series vs Iterative Evaluation

**Task:** Compare computational cost and accuracy:

1. **Method A (Series + Stirling):**
   ```python
   phi_h = sum([comb(2*n, n) * (x**n) for n in range(10)])  # Exact
   phi_h += sum([(4**n/sqrt(pi*n)) * (x**n) for n in range(10,15)])  # Stirling tail
   ```

2. **Method B (Newton-Raphson ζ(Ri)):**
   ```python
   zeta = Ri - Delta * Ri**2  # Seed
   for _ in range(2):
       phi_h = (1 - beta_h * zeta)**(-0.5)
       # Derivative, update...
   ```

**Metrics:**
- FLOP count for single evaluation
- Vectorization efficiency (1000 points)
- Relative error vs exact $(1 - \beta\zeta)^{-1/2}$

**Expected Result:**
- Series: ~50 FLOPs, SIMD-friendly, error $< 10^{-8}$
- Newton: ~80 FLOPs, data-dependent convergence, error $< 10^{-10}$ (if converges)

### E2. Precomputed ζ(Ri) Table

**Task:** Generate lookup table via series reversion:
$$
\zeta = Ri - \Delta Ri^2 + \left(\frac{3}{2}\Delta^2 - \frac{1}{2}c_1\right)Ri^3 + O(Ri^4)
$$

1. Compute coefficients from central binomial Cauchy products.
2. Evaluate ζ at $Ri \in [0, 0.5]$ (1000 points).
3. Store as `zeta_ri_table.npy`.
4. Benchmark interpolation vs Newton solver.

**Deliverable:**
- Plot: interpolation error vs table resolution
- Timing comparison (1 million evaluations)

### E3. Machine Precision Limits

**Question:** At what order $N$ does the series saturate to machine precision ($\sim 10^{-16}$ for float64)?

**Method:**
- Compute φ_h(ζ) using $N = 5, 10, 15, 20, 30$ terms.
- Compare to quad-precision $(1 - \beta\zeta)^{-1/2}$.
- Plot error vs $N$.

**Expected:** Series converges exponentially until rounding error dominates at $N \sim 20$–25.

---

## Submission Guidelines

- **Format:** Typed (LaTeX/Markdown preferred); handwritten derivations acceptable if scanned clearly.
- **Code:** Commented Python/Julia/MATLAB with plots.
- **Length:** ~6–8 pages including derivations, plots, and discussion.

---

## Learning Objectives

By completing this problem, you will:

1. Master the connection between **special function generating functions** and atmospheric stability profiles.
2. Understand **domain restrictions** (poles) and convergence radii in physical applications.
3. Compute Richardson number curvature to **arbitrary precision** using series methods.
4. Develop intuition for when **exact series** outperform **finite-difference** approximations.

---

## Extension (Optional, +10%)

### E1. Hypergeometric Form

The power-law $\phi = (1 - \beta\zeta)^{-\alpha}$ can be written using the **Gaussian hypergeometric function**:
$$
\phi(\zeta) = {}_2F_1(\alpha, 1; 1; \beta\zeta)
$$

**Derive** the first five terms of the hypergeometric series and compare to the central binomial expansion for $\alpha = 1/2$.

### E2. Variable Exponent

Consider a **variable-exponent** form (VEXP):
$$
\phi_h(\zeta) = (1 - \beta_h \zeta)^{-\alpha_h(1 + \eta_h \zeta)}
$$

**Show** that for small $\eta_h$, the series begins:
$$
\phi_h(\zeta) \approx \sum_{n=0}^\infty \binom{2n}{n} \left(\frac{\beta_h \zeta}{4}\right)^n \left[1 + \eta_h \zeta \ln(1 - \beta_h\zeta) + \cdots\right]
$$

**Question:** Does this correction preserve the neutral curvature $2\Delta$?

---

## References

- Businger, J. A., et al., 1971: Flux–profile relationships in the atmospheric surface layer. *J. Atmos. Sci.*, **28**, 181–189.
- Abramowitz, M., and I. A. Stegun, 1964: *Handbook of Mathematical Functions*. National Bureau of Standards, Section 15.1 (hypergeometric functions).
- Graham, R. L., D. E. Knuth, and O. Patashnik, 1994: *Concrete Mathematics*, 2nd ed. Addison-Wesley, Chapter 5 (binomial coefficients).

---

**Instructor Note:** This problem can be adapted for **computational emphasis** (Part D as main focus) or **theoretical emphasis** (Parts A–C with analytic proofs). Solution key includes Mathematica/Maple symbolic verification.
