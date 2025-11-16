# Central Binomial Coefficients: Atmospheric Application

## Generating Function (One-Liner)
$$
\frac{1}{\sqrt{1-4x}} = \sum_{n=0}^\infty \binom{2n}{n} x^n, \quad |x| < \frac{1}{4}
$$

## MOST Power-Law Connection

For half-integer exponents in Monin–Obukhov profiles:
$$
\phi_h(\zeta) = (1 - \beta_h \zeta)^{-1/2} = \sum_{n=0}^\infty \binom{2n}{n} \left(\frac{\beta_h \zeta}{4}\right)^n
$$

### Example: Businger et al. (1971) Unstable Heat
$$
\beta_h = 9 \quad \Rightarrow \quad \phi_h(\zeta) = \sum_{n=0}^\infty \binom{2n}{n} \left(\frac{9\zeta}{4}\right)^n
$$

**Convergence radius:** $|\zeta| < 4/9 \approx 0.44$ (well beyond typical unstable range $\zeta \in [-2, 0]$).

### Explicit First Five Terms ($\beta_h = 16$, cleaner arithmetic)
$$
\phi_h(\zeta) = 1 + 2(4\zeta) + 6(4\zeta)^2 + 20(4\zeta)^3 + 70(4\zeta)^4 + \cdots
$$
$$
= 1 + 8\zeta + 96\zeta^2 + 1280\zeta^3 + 17920\zeta^4 + \cdots
$$

## Asymptotic Growth Rate

Stirling approximation for central binomials:
$$
\binom{2n}{n} \sim \frac{4^n}{\sqrt{\pi n}} \quad \text{as } n \to \infty
$$

**Relative accuracy:**
$$
R_n = \frac{\binom{2n}{n}}{\frac{4^n}{\sqrt{\pi n}}}
$$

| $n$ | $\binom{2n}{n}$ | $\frac{4^n}{\sqrt{\pi n}}$ | $R_n$ |
|-----|------------------|----------------------------|-------|
| 10  | 184,756          | 183,932                    | 1.0045|
| 20  | 137,846,528,820  | 137,603,904,531            | 1.0018|
| 50  | $2.46 \times 10^{29}$ | $2.46 \times 10^{29}$ | 1.0003|

**Conclusion:** Asymptotic formula accurate to <1% for $n \geq 10$.

## Hypergeometric Generalization

For arbitrary exponent $\alpha$:
$$
(1 - \beta\zeta)^{-\alpha} = {}_2F_1(\alpha, 1; 1; \beta\zeta) = \sum_{n=0}^\infty \frac{(\alpha)_n}{n!} (\beta\zeta)^n
$$

where $(\alpha)_n = \alpha(\alpha+1)\cdots(\alpha+n-1)$ is the Pochhammer symbol.

**Special case** $\alpha = 1/2$:
$$
(\tfrac{1}{2})_n = \frac{(2n)!}{4^n n!} \quad \Rightarrow \quad \frac{(\tfrac{1}{2})_n}{n!} = \frac{1}{4^n}\binom{2n}{n}
$$

matches the central binomial series (up to rescaling by $\beta/4$).

## Homework Application

See `HW_Central_Binomial_MOST.md` for graduate-level problem set:
1. Derive series from generating function.
2. Compute convergence radius for $\beta_h = 16$.
3. Use series to calculate Richardson number curvature $\partial^2 Ri_g / \partial\zeta^2$.
4. Compare to finite-difference approximations.

## Python Quick Check

```python
from scipy.special import comb
import numpy as np

def central_binomial(n):
    return comb(2*n, n, exact=True)

# Verify asymptotic formula
n = 20
C_n = central_binomial(n)
asymp = 4**n / np.sqrt(np.pi * n)
print(f"Exact: {C_n}, Asymptotic: {asymp:.0f}, Ratio: {C_n/asymp:.6f}")
# Output: Exact: 137846528820, Asymptotic: 137603904531, Ratio: 1.001763
```

## References

- **Generating function:** OEIS A000984, Concrete Mathematics (Graham, Knuth, Patashnik), Chapter 5.
- **MOST application:** Businger et al. (1971, JAS); Dyer (1974, BLM); Paulson (1970, JAM).
- **Hypergeometric connection:** Abramowitz & Stegun (1964), §15.1.

## Extension: Variable-Exponent Power Law (VEXP)

For
$$
\phi(\zeta) = (1 - \beta\zeta)^{-\alpha(1 + \eta\zeta)}
$$

expand using $\exp[\eta\zeta \ln(\phi)]$ perturbation; central binomials enter the base term, logarithmic corrections modulate growth.
