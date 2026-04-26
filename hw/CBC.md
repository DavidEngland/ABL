# Homework: Central Binomial Coefficients and the Structure of MOST Stability Functions

**Course:** Atmospheric Boundary Layer Physics  
**Topic:** CBC Expansions, Asymptotic Regimes, and Dynamic Critical Richardson Number  
**Prerequisites:** MOST_Derivations_Problems.md; basic power-series and binomial theorem; Python (scipy, numpy)

> **Why this matters.** The power-law MOST functions contain an exact structural identity:
> when \(b_m = b_h\) the gradient Richardson number equals the stability parameter,
> \(Ri_g = \zeta\), with no iterative solver. The proof goes through the *central binomial
> coefficient* (CBC) expansion term by term. This set of problems walks you from that
> expansion to asymptotics in the convective limit and finally to the dynamic critical
> Richardson number that unifies the unstable and stable regime limits.

---

## Background: Notation

Throughout, use the **code-consistent notation**:

| Symbol | Meaning | Typical values |
|---|---|---|
| \(b_m\) | Unstable momentum coefficient (literature: \(\gamma_m\)) | 15–19 |
| \(b_h\) | Unstable heat coefficient (literature: \(\gamma_h\)) | 9–16 |
| \(\beta_m\) | Stable momentum coefficient | 4.7–5.3 |
| \(\beta_h\) | Stable heat coefficient | 5.0–6.35 |
| \(\alpha_m\) | Unstable momentum exponent | 1/4 |
| \(\alpha_h\) | Unstable heat exponent | 1/2 |
| \(\zeta = z/L\) | Stability parameter | \(<0\) unstable, \(>0\) stable |
| \(\eta = -\zeta\) | Convective variable (\(\eta > 0\) unstable) | |

The MOST forms on the unstable branch (\(\zeta \le 0\)):
$$
\phi_m(\zeta) = (1 - b_m\,\zeta)^{-1/4}, \qquad
\phi_h(\zeta) = (1 - b_h\,\zeta)^{-1/2}.
$$

---

## Problem 1 — CBC Expansion of \(\phi_h\) (25 points)

### Part A: Derive the power series (10 pts)

The generalized binomial theorem states that for \(|x| < 1\):
$$
(1 - x)^{-1/2} = \sum_{n=0}^{\infty} \frac{(2n)!}{4^n\,(n!)^2}\,x^n
= \sum_{n=0}^{\infty} \binom{2n}{n} \frac{x^n}{4^n}.
$$

1. Substitute \(x = b_h\,\zeta\) to write \(\phi_h(\zeta)\) as a power series in \(\zeta\).
   Express each coefficient \(c_n^{(h)}\) in closed form.

2. Write out the first five coefficients \(c_0^{(h)}, \ldots, c_4^{(h)}\) for \(b_h = 16\)
   (Dyer/Brutsaert default).

3. What is the **radius of convergence** of this series? Express it in terms of \(b_h\).
   Explain physically what quantity sets this radius.

**Expected coefficients for \(b_h = 16\):** \(1,\;8,\;96,\;1280,\;17920,\;\ldots\)

### Part B: Two-term recurrence (5 pts)

Show that successive coefficients satisfy

$$
c_{n+1}^{(h)} = c_n^{(h)} \cdot \frac{b_h\,(2n+1)}{2(n+1)}, \qquad c_0^{(h)} = 1.
$$

Derive this from the ratio \(c_{n+1}/c_n\) using the explicit form in Part A.

### Part C: Numerical verification (10 pts)

Using `code/phi_h_central_binomial.py` as a starting point:

1. Compute \(\phi_h(-0.3)\) using the exact formula \((1 - b_h\zeta)^{-1/2}\) for \(b_h = 9\)
   and \(b_h = 16\).
2. Compute the same value using truncated CBC series with \(N = 5, 10, 20\) terms.
3. Plot relative error vs. \(N\) for \(\zeta = -0.3\) and \(\zeta = -0.05\).  
   What convergence rate do you observe?
4. What is the maximum \(|\zeta|\) for which the \(N=10\) truncation achieves
   relative error below \(10^{-4}\)? Relate your answer to the radius of convergence.

---

## Problem 2 — CBC Expansion of \(\phi_m\) (20 points)

### Part A: Generalized binomial coefficients (10 pts)

For exponent \(\alpha_m = 1/4\), the generalized binomial theorem gives

$$
(1-x)^{-1/4} = \sum_{n=0}^{\infty} \frac{\Gamma(n + 1/4)}{\Gamma(1/4)\,n!}\,x^n.
$$

1. Show that the recurrence for \(c_n^{(m)}\) (with \(x = b_m\zeta\)) is

$$
c_{n+1}^{(m)} = c_n^{(m)} \cdot \frac{b_m\,(4n+1)}{4(n+1)}, \qquad c_0^{(m)} = 1.
$$

2. Write out \(c_0^{(m)}, \ldots, c_4^{(m)}\) for \(b_m = 16\).

### Part B: Comparing coefficient growth (10 pts)

1. At large \(n\), the central-binomial coefficients behave as
   \(c_n^{(h)} \sim \bigl(\frac{b_h}{4}\bigr)^n \cdot \frac{4^n}{\sqrt{\pi n}}\) (Stirling).
   Derive the analogous large-\(n\) estimate for \(c_n^{(m)}\) using Stirling's approximation
   for \(\Gamma(n+1/4)/\Gamma(1/4)\,n!\).

2. For the same \(b_m = b_h = b\), show that \(c_n^{(h)} > c_n^{(m)}\) for all \(n \ge 1\),
   consistent with the slower singular response of \(\phi_m\) compared to \(\phi_h\).

3. Explain in one sentence the physical meaning of this inequality for
   momentum vs heat transport in convective conditions.

---

## Problem 3 — The Structural Identity \(Ri_g = \zeta\) (25 points)

### Part A: Algebraic proof (10 pts)

Assume \(b_m = b_h = b\).

1. Show that \(\phi_h(\zeta) = \phi_m(\zeta)^2\) as an identity of real functions on \(\zeta \le 0\).
2. Substitute into \(Ri_g = \zeta\,\phi_h / \phi_m^2\) to prove \(Ri_g = \zeta\).
3. Prove the same identity **term by term in the CBC series**: i.e., show that the CBC
   coefficients of \(\phi_h\) equal the convolution-square of the coefficients of \(\phi_m\).

   *Hint:* If \(\phi_m = \sum a_n \zeta^n\) then \(\phi_m^2 = \sum_n (\sum_{k=0}^n a_k a_{n-k})\zeta^n\).
   Check \(n=0,1,2\) explicitly.

### Part B: Degeneracy and what it means for modeling (5 pts)

1. In one paragraph: what is the **practical consequence** of \(Ri_g = \zeta\) for numerical
   models that compute the Obukhov length iteratively?
2. Under which observational dataset was this degeneracy *nearly* satisfied?  
   (Consult the parameter table in `notes/parameters.md`.)

### Part C: Breaking the identity (10 pts)

Now let \(b_m \neq b_h\) (general case). Let \(\delta = b_h - b_m\).

1. Expand \(Ri_g(\zeta)\) to first order in \(\delta\) and first order in \(\zeta\) to find:
$$
Ri_g \approx \zeta\!\left[1 + \frac{\delta}{4}\,\zeta + O(\zeta^2, \delta^2)\right].
$$

   *Hint:* Write \(\phi_h = \phi_m^2 \cdot (1-b_m\zeta)^{-\delta/4} \cdot \ldots\) and
   linearize in \(\delta\).

2. For Businger (1971) values \(b_m=15,\,b_h=9\):
   - Compute \(\delta = b_h - b_m\).
   - At \(\zeta = -0.5\) evaluate the first-order correction term.
   - Is the correction large or small relative to \(\zeta\)? What does this imply for the
     utility of the degenerate approximation at moderate instability?

---

## Problem 4 — Asymptotic Expansion in the Convective Limit (20 points)

**Background:** substitute \(\eta = -\zeta > 0\) so that the strongly unstable UBL
corresponds to \(\eta \to +\infty\).

### Part A: Leading-order asymptotics (10 pts)

1. Show that
$$
\phi_h(-\eta) = (1 + b_h\,\eta)^{-1/2} \sim (b_h\,\eta)^{-1/2}
\left[1 - \frac{1}{2b_h\,\eta} + \frac{3}{8(b_h\,\eta)^2} - \cdots\right]
\qquad \eta \to \infty.
$$
   Derive at least three terms by expanding \((1 + b_h\eta)^{-1/2} = (b_h\eta)^{-1/2}(1 + (b_h\eta)^{-1})^{-1/2}\)
   and then applying the CBC series to the inner factor.

2. Derive the analogous three-term asymptotic for \(\phi_m(-\eta)\).

3. Verify that both asymptotic series are **regular inverse-power** hierarchies (no
   oscillatory or exponential terms). Why is this physically reasonable?

### Part B: Decay rates and their implications (5 pts)

1. From the leading-order terms, confirm that \(\phi_h \propto \eta^{-1/2}\) and
   \(\phi_m \propto \eta^{-1/4}\).
2. The turbulent heat flux scales as \(w'\theta' \propto \phi_h^{-1}\); the friction
   velocity scales as \(u_* \propto \phi_m^{-1}\).  Based on the decay rates, which
   quantity collapses more slowly as convection intensifies?
3. Relate your answer to the phrase "scalar closures are intrinsically less sensitive
   to extreme instability than momentum closures."

### Part C: Python verification (5 pts)

Write a short Python function `phi_asymptotic(eta, b, alpha, n_terms)` that returns the
\(n\)-term asymptotic expansion for general exponent \(\alpha\) (use \(\alpha=1/2\) for
\(\phi_h\), \(\alpha=1/4\) for \(\phi_m\)).  

Plot exact vs asymptotic (3 terms) for \(\eta \in [1, 50]\) and both exponents.
At what \(\eta\) does the 3-term asymptotic achieve relative error below \(1\%\)?

---

## Problem 5 — Dynamic Critical Richardson Number (10 bonus points)

### Part A: Natural regime limits from the series

1. **Unstable limit:** The CBC series for \(\phi_m\) converges for \(|b_m\zeta| < 1\),
   i.e., the natural "UBL critical" in \(\zeta\)-space is \(\zeta_{c,UBL} = -1/b_m\)
   (beyond the radius of convergence in the negative direction).  In the degenerate case
   \(Ri_g = \zeta\), what is the corresponding \(Ri_c^{UBL}\)?  Express in terms of \(b_m\).

2. **Stable limit:** Starting from the linear stable inversion
   \(f_m(Ri) = (1 - \beta\,Ri)^2\), identify \(Ri_c^{SBL}\) as the zero of \(f_m\).
   Express in terms of \(\beta\).

3. Tabulate \(Ri_c^{UBL}\) and \(Ri_c^{SBL}\) for each canonical parameter set:

| Set | \(b_m\) | \(\beta_m\) | \(Ri_c^{UBL}\) | \(Ri_c^{SBL}\) |
|---|---|---|---|---|
| Businger 1971 | 15 | 4.7 | | |
| Dyer 1974 | 16 | 5.0 | | |
| Högström 1996 | 19 | 5.3 | | |

### Part B: Bridging the limits

The dynamic \(Ri_c^*\) framework proposes that the effective critical Richardson number
varies smoothly between \(Ri_c^{UBL}\) and \(Ri_c^{SBL}\) as a function of the local
flow state.

1. Propose a simple monotone interpolation formula for \(Ri_c^*(\zeta)\) that recovers
   \(Ri_c^{UBL}\) as \(\zeta \to -\infty\) and \(Ri_c^{SBL}\) as \(\zeta \to +\infty\).
   (No unique answer; justify your choice.)

2. How does your formula change \(\Delta = \alpha_h b_h - 2\alpha_m b_m\)?  Does it
   drive \(\Delta\) toward zero (the degenerate case)?

3. Describe in 2–3 sentences how this dynamic \(Ri_c^*\) would be implemented in the
   `URC` profile family in `code/profiles.py`.

---

## Summary Checklist

When you have completed this homework you should be able to:

- [ ] Derive and use the CBC recurrence for \(\phi_h\) (exponent 1/2) and \(\phi_m\) (exponent 1/4)
- [ ] State and prove the structural identity \(Ri_g = \zeta\) when \(b_m = b_h\)
- [ ] Perform term-by-term verification of the identity in the power series
- [ ] Write the 3-term asymptotic expansion in \(\eta = -\zeta\) for both functions
- [ ] Explain why \(\phi_h\) decays faster than \(\phi_m\) in the convective limit
- [ ] Relate the radius of convergence to a critical \(Ri\) in each regime
- [ ] Connect Businger/Dyer parameter differences to the perturbative correction \(\delta\)
- [ ] Sketch how a dynamic \(Ri_c^*\) bridges UBL and SBL limiting critical values

---

## References

- Businger et al. (1971) *J. Atmos. Sci.*
- Dyer (1974) *Bound.-Layer Meteor.*
- Brutsaert (1982) *Evaporation into the Atmosphere*
- Högström (1996) *Bound.-Layer Meteor.*
- `code/phi_h_central_binomial.py` — working CBC implementation
- `code/profiles.py` — profile registry including URC and Ri-based closures
- `notes/parameters.md` — canonical parameter table and CBC structural summary
- `Dynamic Critical Richardson Number Framework.md` — full dynamic \(Ri_c^*\) formulation
