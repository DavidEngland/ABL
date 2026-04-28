Here is a consolidated, research‑grade synthesis of the parameter sets, regimes, and specialized stable‑layer formulations you listed—organized to make the physics, numerical‑model relevance, and cross‑study consistency immediately clear.
(All citations refer to the sources you summarized; no new external claims are introduced.)

---

1. Canonical Parameter Sets for Businger–Dyer Type Similarity Functions

### Quick Reference Tables

| Set / Source | k | a_h^{-1} (or Pr0inv) | Unstable coeff. (m, h) | Stable coeff. (m, h) | Notes |
|---|---:|---:|---:|---:|---|
| Businger et al. (1971) | 0.35 | 1.35 | (15.0, 9.0) | (4.7, 6.35) | Kansas dataset baseline |
| Dyer (1974) / Brutsaert (1982) | 0.40 | 1.00 | (16.0, 16.0) | (5.0, 5.0) | Common WRF/LES default |
| Högström (1996) | 0.40 | 1.00 | (19.0, 11.6) | (5.3, 5.3) | Stronger unstable momentum curvature |

Notes:
Unstable coefficients are often denoted \(\gamma_m, \gamma_h\) in literature and \(b_m, b_h\) in code-oriented implementations.

| Regime Type | Range | Practical Meaning | Typical Modeling Guidance |
|---|---|---|---|
| Near-neutral unstable | \(-1 < \zeta < 0\) | Weak instability | Classical Businger-Dyer forms usually reliable |
| Strongly unstable | \(\zeta < -1\) | Convective dominance increases | Consider free-convection behavior checks |
| Formal unstable validity (classical) | \(-2 < \zeta < 0\) | Common cited calibration window | Avoid overextending without validation |
| Moderately stable | \(0 < \zeta \le 1\) | Stratification suppresses mixing | Linear stable forms often acceptable |
| Very/extremely stable | \(\zeta > 1\) | Intermittent, z-less turbulence regimes | Prefer specialized stable-layer formulations |
| Ri nearly neutral | \(0 < Ri < 0.02\) | Turbulent and close to neutral | Standard closures work well |
| Ri weakly stable | \(0.02 < Ri < 0.12\) | Increasing suppression of turbulence | Monitor stability-function sensitivity |
| Ri very stable | \(0.12 < Ri < 0.7\) | Intermittency more likely | Use robust/capped stable closures |
| Ri extremely stable | \(Ri > 0.7\) | Possible laminar tendency | Turbulence may collapse intermittently |
| Critical threshold | \(Ri_c \approx 0.21\text{–}0.25\) | Onset of laminarization risk | Useful switch point in higher-order closures |

The table you provided captures the core empirical constants used in MOST‑based surface‑layer schemes. These constants define:

• Unstable branch:\phi_m = (1 - \gamma_m \zeta)^{-1/4}, \qquad
\phi_h = a_h^{-1}(1 - \gamma_h \zeta)^{-1/2}

• Stable branch:\phi_m = 1 + \beta_m \zeta, \qquad
\phi_h = a_h^{-1} + \beta_h \zeta



Interpretation of the constants

• \(k\) (von Kármán constant): Controls the slope of the log‑law; values range 0.35–0.41, with modern consensus near 0.40.
• \(a_h^{-1}\): Neutral limit of the temperature gradient function; typically 0.95–1.35, depending on whether the source reports \(a_h\) or its inverse.
• \(\gamma_m, \gamma_h\): Strength of the unstable curvature; larger values steepen the gradient as \(\zeta \to -1/\gamma\).
• \(\beta_m, \beta_h\): Linear slope in stable conditions; classical values cluster around 5.


Cross‑study consistency

• Businger (1971) and Dyer (1974) remain the “Kansas standard,” widely used in NWP and LES surface‑layer schemes.
• Högström (1996) provides slightly larger unstable coefficients (\(\gamma_m=19\)), consistent with later field campaigns.
• Brutsaert (1982) and Paulson (1972) essentially reproduce the Dyer constants.


These constants are still embedded in WRF, ECMWF IFS, RUC, ARPS, and many land‑surface models.

---

2. Validity Domains and Regime Structure

Unstable regimes (`\(\zeta < 0\)`)

• Near‑neutral: \(-1 < \zeta < 0\).
Kansas‑type functions perform extremely well here.
• Strongly unstable: \(\zeta < -1\).
Classical forms begin to deviate; free‑convection scaling approaches\phi_h \sim (-\zeta)^{-1/3}.

• Formal validity: Many classical functions are cited only for \(-2 < \zeta < 0\).


Stable regimes (`\(\zeta > 0\)`)

• Moderately stable: \(0 < \zeta \le 1\).
Linear forms (\(\phi = 1 + 5\zeta\)) remain adequate.
• Very/extremely stable: \(\zeta > 1\).
Linear forms fail; turbulence becomes intermittent, non‑stationary, and “z‑less.”


Richardson‑number regimes

• Nearly neutral: \(0 < Ri < 0.02\)
• Weakly stable: \(0.02 < Ri < 0.12\)
• Very stable: \(0.12 < Ri < 0.7\)
• Extremely stable: \(Ri > 0.7\)
• Critical limit: \(Ri_c \approx 0.21–0.25\) (onset of laminarization)


These Ri thresholds are crucial for turbulence switches in higher‑order closures (e.g., MYNN, QNSE, EFB).

---

3. Specialized Stable‑Layer Formulations

These schemes were developed to avoid the decoupling problem (runaway cooling, vanishing fluxes) that plagues linear MOST functions at high stability.

Beljaars–Holtslag (1991)

• Uses exponential decay terms with constants \(a=1, b=0.667, c=5, d=0.35\).
• Smooths the transition into very stable regimes.
• Adopted in ECMWF and many land‑surface models.


Cheng & Brutsaert (2005)

• Momentum: \(a=6.1, b=2.5\)
• Heat: \(c=5.3, d=1.1\)
• Extends validity to \(\zeta \approx 20\).
• Adopted in MYNN and several Arctic‑focused schemes.


SHEBA (2007)

• Derived from Arctic sea‑ice data up to \(\zeta \approx 100\).
• Momentum: \(a_m=5, b_m=0.3\)
• Heat: \(a_h=5, b_h=0.4\)
• Essential for polar modeling where classical MOST fails.


---

4. How to Use These Parameters in Numerical Models

1. Choose a base similarity family (Businger–Dyer, Beljaars–Holtslag, Cheng–Brutsaert, SHEBA).
2. Select constants appropriate for the regime and surface type.
3. Ensure continuity across neutral, unstable, and stable branches.
4. Apply Ri‑based turbulence switches for very stable conditions.
5. Validate against site‑specific data, especially in Arctic or nocturnal SBL cases.

---

5. Central Binomial Coefficient Expansions — A Key Structural Result

**This is a big deal.** The MOST power-law functions

$$
\phi_m(\zeta) = (1 - b_m\,\zeta)^{-1/4}, \qquad \phi_h(\zeta) = (1 - b_h\,\zeta)^{-1/2}
$$

admit exact power-series representations near neutral via the generalized binomial theorem.
For \(\phi_h\) the series coefficients collapse to *central binomial coefficients*:

$$
\phi_h(\zeta) = \sum_{n=0}^{\infty} \binom{2n}{n}\!\left(\frac{b_h}{4}\right)^n \zeta^n
= 1 + \frac{b_h}{2}\,\zeta + \frac{3 b_h^2}{8}\,\zeta^2 + \cdots
$$

with the two-term recurrence
$$
c_{n+1} = c_n \cdot \frac{b_h\,(2n+1)}{2(n+1)}, \qquad c_0 = 1.
$$

For \(b_h = 16\) (Dyer default): coefficients are \(1,\,8,\,96,\,1280,\,\ldots\) — the central binomials scaled by powers of 4.

For \(\phi_m\) (exponent 1/4), the coefficients follow the analogous generalized-binomial recurrence:
$$
c_{n+1}^{(1/4)} = c_n^{(1/4)} \cdot \frac{b_m\,(4n+1)}{4(n+1)}, \qquad c_0^{(1/4)} = 1.
$$

**Radius of convergence:** both series converge for \(|\zeta| < 1/b_h\) (resp. \(1/b_m\)), set by the branch point.

### 5.1 The Structural Identity (degenerate case \(b_m = b_h = b\))

When the unstable coefficients are equal,

$$
\phi_h(\zeta) = (1-b\zeta)^{-1/2} = \bigl[(1-b\zeta)^{-1/4}\bigr]^2 = \phi_m(\zeta)^2.
$$

This is an **exact algebraic identity**, not an approximation. Substituting into the
MOST Richardson-number relation \(Ri_g = \zeta\,\phi_h / \phi_m^2\) gives immediately

$$
\boxed{Ri_g = \zeta} \qquad (b_m = b_h).
$$

The gradient Richardson number equals the stability parameter with no iterative solver required
on the unstable branch. This degeneracy is an exact structural consequence of the half-power
exponent relationship \(\alpha_h = 2\alpha_m\) and equal coefficients; it holds term by term
in the CBC expansion.

### 5.2 Asymptotic Expansion in the UBL (\(\eta = -\zeta\), \(\eta \to \infty\))

Substituting \(\eta = -\zeta > 0\) for the strongly unstable branch:

$$
\phi_m \sim (b_m\,\eta)^{-1/4}\!\left[1 - \tfrac{1}{4}(b_m\eta)^{-1}
+ \tfrac{5}{32}(b_m\eta)^{-2} - \cdots\right]
$$

$$
\phi_h \sim (b_h\,\eta)^{-1/2}\!\left[1 - \tfrac{1}{2}(b_h\eta)^{-1}
+ \tfrac{3}{8}(b_h\eta)^{-2} - \cdots\right]
$$

Both series are regular inverse-power hierarchies — no oscillatory terms.
\(\phi_h \propto \eta^{-1/2}\) decays faster than \(\phi_m \propto \eta^{-1/4}\):
scalar closures are intrinsically less sensitive to extreme instability than momentum closures.

### 5.3 Parameter Mismatch and the Dynamic \(Ri_c\) Resolution

In practice \(b_m \neq b_h\) (e.g., Businger: \(b_m=15,\,b_h=9\)), which **breaks** the
\(Ri_g=\zeta\) identity. The correction (first-order perturbation):

$$
Ri_g \approx \zeta\!\left[1 - (a_h - 2a_m)\,b\,\zeta + \cdots\right]
$$

where \(\Delta = \alpha_h b_h - 2\alpha_m b_m\) is the neutral curvature invariant.

The **dynamic critical Richardson number** framework reconciles the mismatch by allowing
\(Ri_c^*\) to vary continuously between regime limits:

| Regime | Natural limit | Origin |
|---|---|---|
| Strongly unstable (UBL) | \(Ri_c^{UBL} = -1/b_m\) | Convergence radius of \(\phi_m\) CBC series |
| Stable (SBL) | \(Ri_c^{SBL} = +1/\beta\) | Pole of stable-linear inversion |

A smooth transition \(Ri_c^* = Ri_c^*(S,\,\Gamma,\,\text{TKE})\) spanning
\([-1/b_m,\,+1/\beta]\) provides a parameter-consistent closure that avoids
the discontinuity at \(\zeta=0\) and eliminates the need for a fixed empirical
threshold. See `Dynamic Critical Richardson Number Framework.md` for the full formulation.

### 5.4 CBC vs Power Series vs Asymptotic: Comparison

| Representation | Valid range | Coefficients | Evaluation cost | Best use |
|---|---|---|---|---|
| Full CBC series (truncated to N terms) | \(|\zeta| < 1/b_h\) | Exact, 2-term recurrence | O(N) | Near-neutral regime |
| Stirling-approximated tail | \(n \gg 1\) | \(\approx 4^n/\sqrt{\pi n}\) | O(N_tail) | Accelerating convergence |
| Asymptotic in \(\eta=-\zeta\) | \(\eta \gg 1\) | Exact, inverse-power | O(few terms) | Strongly unstable UBL |
| Linear stable form | \(0 < \zeta \ll 1/\beta\) | Exact, trivial | O(1) | Near-neutral stable |

See `code/phi_h_central_binomial.py` for a working implementation of the CBC + Stirling hybrid.


---

If you want, I can assemble a unified parameter table, derive the corresponding flux–profile integrals, or produce a regime‑aware similarity function library (e.g., in Python, Fortran, or Julia) tailored to your Arctic inversion experiment.