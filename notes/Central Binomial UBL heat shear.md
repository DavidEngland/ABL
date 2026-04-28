Central Binomial and UBL heat-shear function.

Let
$$
\phi_h(\zeta)=\frac{1}{\sqrt{1-b_h\zeta}}=(1-b_h\zeta)^{-1/2},
$$
where \(b_h\) is the unstable heat coefficient. For Dyer/Brutsaert, \(b_h=16\).

**1. Core expansion (near neutral)**
For \(|b_h\zeta|<1\) (so \(|\zeta|<1/b_h\)):
$$
(1-x)^{-1/2}=\sum_{n=0}^{\infty}\frac{\binom{2n}{n}}{4^n}x^n,\qquad x=b_h\zeta,
$$
so
$$
\phi_h(\zeta)=\sum_{n=0}^{\infty}\binom{2n}{n}\left(\frac{b_h}{4}\right)^n\zeta^n.
$$

For \(b_h=16\):
$$
\phi_h(\zeta)=1+8\zeta+96\zeta^2+1280\zeta^3+17920\zeta^4+\cdots.
$$

This is the central-binomial series.

**2. Odd/even split**
Write
$$
\phi_h(\zeta)=F_{\text{even}}(\zeta)+F_{\text{odd}}(\zeta),
$$
with
$$
F_{\text{even}}=\sum_{k\ge0}c_{2k}\zeta^{2k},\qquad
F_{\text{odd}}=\sum_{k\ge0}c_{2k+1}\zeta^{2k+1},
$$
\(c_n=\binom{2n}{n}(b_h/4)^n\).

Why useful: parity-separated truncations can reduce bias when matching symmetric constraints around \(\zeta=0\), and are handy in closure fitting where you want separate slope-curvature control.

**3. Regime interpretation**

1. Unstable \((\zeta<0)\): real, bounded, decays for strong instability.
   Let \(\eta=-\zeta>0\):
   $$
   F=\frac{1}{\sqrt{1+b_h\eta}}
   \sim (b_h\eta)^{-1/2}\left(1-\frac{1}{2b_h\eta}+\frac{3}{8(b_h\eta)^2}-\cdots\right),\quad \eta\to\infty.
   $$
   So strong-convective limit is algebraic decay \(\propto |\zeta|^{-1/2}\).

2. Near neutral \((|\zeta|\ll1)\): use Taylor/binomial above.

3. Stable \((0<\zeta<1/b_h)\): function rises and blows up at \(\zeta=1/b_h\).
   That pole is usually not physically acceptable as a universal stable closure over all stable conditions.

4. “Laminar” very stable: if you need continuation past \(1/b_h\), you must replace/regularize the model (piecewise closure, capped/saturated form, or rational continuation). Raw formula becomes singular then complex for \(\zeta>1/b_h\).

**4. Continuation strategies (practical numerics/physics)**

1. Piecewise physics-consistent:
   - Use exact \(F\) for unstable and near-neutral.
   - Switch at small positive \(\zeta_s\) to a stable law without pole (e.g. linear, rational, or capped growth), enforcing \(C^1\) continuity.

2. Padé continuation:
   - Fit a \([m/n]\) rational approximant to neutral coefficients and impose no positive-real pole in target stable interval.
   - Better global behavior than finite Taylor truncation near \(\zeta\approx1/16\).

3. Matched asymptotics:
   - Neutral series for \(|\zeta|\ll1\).
   - Large-\(|\zeta|\) asymptotic for unstable tail.
   - Blend with smooth weight \(w(\zeta)\).

4. Chebyshev on bounded intervals:
   - For numerical robustness over \(\zeta\in[\zeta_{\min},\zeta_{\max}]\), Chebyshev minimax approximations reduce max error vs Taylor.

**5. Numerical techniques**

1. Stable evaluation near neutral:
   $$
   F=\exp\!\left(-\frac12\log(1-16\zeta)\right)
   $$
   using `log1p(-b_h\zeta)` to reduce cancellation.

2. Fast coefficient recurrence:
   If \(c_n=\binom{2n}{n}(b_h/4)^n\), then
   $$
   c_{n+1}=c_n\cdot \frac{b_h(2n+1)}{2(n+1)},\qquad c_0=1.
   $$
   Good for on-the-fly series generation without factorials.

3. Horner form for truncated polynomial to reduce roundoff and FLOPs.

**6. Advantages of central-binomial representation**

1. Exact analytic coefficients with simple recurrence.
2. Clear radius of convergence from branch point at \(\zeta=1/b_h\).
3. Easy to build higher-order neutral-consistent closures.
4. Coefficients encode curvature growth naturally (helpful in sensitivity analysis).
5. Good for symbolic manipulations (derivatives/integrals remain structured).

**Main caveat**
Central-binomial series is excellent locally near \(\zeta=0\), but alone is not a full-regime closure because the positive-\(\zeta\) pole limits physical and numerical utility in stable/laminar regimes.

**Parameter-dependent Richardson limits**

For recent dynamic-\(Ri_c\) work, the natural regime limits should be tied to the coefficients:

$$
Ri_{c,m}^{UBL} = -\frac{1}{b_m}, \qquad Ri_{c,h}^{UBL} = -\frac{1}{b_h}.
$$

These are the natural unstable momentum and heat branch radii implied by the CBC / Gegenbauer expansions. If a single unstable critical value is needed for a mixed closure, the conservative choice is the more restrictive one,
$$
Ri_c^{UBL} = -\frac{1}{\max(b_m,b_h)}.
$$

For the stable linear branch,
$$
Ri_c^{SBL} = \frac{1}{\beta},
$$
or component-wise \(1/\beta_m\), \(1/\beta_h\) if momentum and heat closures are treated separately.

Thus the old fixed \(1/16\) picture should be replaced by a parameter-aware one:
the unstable side is controlled by \(b_m\) and \(b_h\), while the stable side is controlled by \(\beta\).

If you want, I can draft a concrete 3-regime piecewise formula (unstable, neutral, stable-laminar) with \(C^1\) matching and tunable constants for WRF-style implementation.
