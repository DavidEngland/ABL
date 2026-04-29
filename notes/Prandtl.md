# Prandtl Number as Spectral Convolution

## Scope

This note reframes the turbulent Prandtl number $Pr_t(\zeta)=\phi_h/\phi_m$ as a spectral ratio, not a direct empirical curve fit.  
The intent is to explain the observed "hockey-stick" shape in stable conditions and provide a compact closure pathway using Gegenbauer/Legendre structure plus optional Pad\'e surrogates.

---

## 1. Spectral Formulation

Map stability to a bounded spectral coordinate,

$$
\xi(\zeta)=\tanh(\alpha\zeta), \qquad \xi\in[-1,1].
$$

Expand momentum in Gegenbauer modes:

$$
\phi_m(\zeta)=\sum_{n=0}^{N} a_n\,C_n^{(\lambda)}\!\bigl(\xi(\zeta)\bigr).
$$

Model scalar suppression by a modal filter $r_n$ ($0<r_n\le 1$, decreasing in $n$):

$$
\phi_h(\zeta)=Pr_0\sum_{n=0}^{N} r_n\,a_n\,C_n^{(\lambda)}\!\bigl(\xi(\zeta)\bigr),
$$

with a convenient choice

$$
r_n=\exp\!\left(-\frac{n}{n_c}\right).
$$

Then

$$
Pr_t(\zeta)=Pr_0\,
\frac{\sum_{n=0}^{N} r_n a_n C_n^{(\lambda)}(\xi)}{\sum_{n=0}^{N} a_n C_n^{(\lambda)}(\xi)}.
$$

Interpretation: higher-order modes (smaller scales) are damped more strongly for heat than momentum.

---

## 2. Why Hockey-Stick Curves Appear

The "hockey stick" has three zones:

1. Near-neutral: $Pr_t\approx Pr_0$ (flat segment).
2. Transition: rapid rise as heat modes are filtered faster than momentum modes.
3. Very stable: approach to a plateau (or slow growth), depending on asymptotic degree.

For a low-order truncation,

$$
\phi_m\approx a_0+a_1\xi+a_2\xi^2,
\qquad
\phi_h\approx Pr_0(a_0+r_1a_1\xi+r_2a_2\xi^2),
$$

with $r_2<r_1<1$, the quadratic and higher curvature terms are reduced more in $\phi_h$ than in $\phi_m$, causing a knee in $Pr_t(\zeta)$.

Useful controls:

- $Pr_0$: neutral intercept.
- $n_c$: filter aggressiveness (smaller $n_c$ gives earlier/steeper knee).
- $\alpha$: maps where in $\zeta$ the knee occurs.
- $\lambda$: redistributes spectral resolution across neutral-to-stable range.

---

## 3. Why Legendre/Gegenbauer Beats Direct Polynomial Fitting

Classical curve fitting in $\zeta$ often couples shape and asymptote poorly. The transform approach separates roles:

1. Coordinate map $\xi(\zeta)$ controls compactification.
2. Basis family $C_n^{(\lambda)}$ controls geometric weighting.
3. Filter $r_n$ controls physics of scalar suppression.

Legendre is the special case $\lambda=1/2$ with uniform weight.  
General Gegenbauer ($\lambda\ne 1/2$) provides a tunable weighting

$$
w_\lambda(\xi)=(1-\xi^2)^{\lambda-1/2},
$$

which shifts where curvature is represented most efficiently.

Practical reading:

- Larger $\lambda$: emphasizes interior structure; often sharpens the mid-stability knee.
- Smaller $\lambda$: increases sensitivity near endpoints; can delay sharp rise toward very stable conditions.

---

## 4. Convolution Link to Momentum-Heat Structure

In canonical power-law MOST with matched parameters ($a_h=1$, $b_m=b_h$):

$$
\phi_h=\phi_m^2.
$$

In ultraspherical space this is a product-to-convolution statement:

$$
C^{(1/2)} = C^{(1/4)}\otimes C^{(1/4)}.
$$

So heat is not independent of momentum; it is the convolved spectrum.  
The Prandtl ratio then measures departures from this exact squaring relation due to parameter mismatch and/or modal filtering.

---

## 5. Pad\'e Surrogate for Fast Models

For production codes, replace the spectral ratio with a monotone rational surrogate in $\zeta$:

$$
Pr_t^{[1/1]}(\zeta)=\frac{p_0+p_1\zeta}{1+q_1\zeta},
\qquad
Pr_t^{[2/2]}(\zeta)=\frac{p_0+p_1\zeta+p_2\zeta^2}{1+q_1\zeta+q_2\zeta^2}.
$$

Recommended constraints during fitting:

1. $Pr_t(0)=Pr_0$.
2. $Pr_t'(0)>0$ (stable increase).
3. No poles for $\zeta\ge 0$.
4. Match target plateau (if desired) from SHEBA/CASES-99 envelope.

Use spectral model output as training target, then fit Pad\'e coefficients once per parameter set.

---

## 6. Suggested Baseline Parameters (starting point)

These are a practical initial guess, not a final calibration:

| Parameter | Baseline | Role |
|---|---:|---|
| $Pr_0$ | 0.85 | Neutral value |
| $n_c$ | 1.2 | Heat-mode damping scale |
| $a_1/a_0$ | 2.0 | First-order response |
| $a_2/a_0$ | 0.5 | Curvature stiffness |
| $\alpha$ | 0.5 | Stability-to-spectral map width |
| $\lambda$ | 0.5 to 0.9 | Knee sharpness/location control |

Expected behavior target:

- $\zeta\lesssim 0.05$: $Pr_t\approx 0.8$ to 1.0.
- $0.1\lesssim\zeta\lesssim 1$: rapid rise.
- $\zeta\gtrsim 1$: weak-growth/plateau regime near 3 to 5.

---

## 7. Implementation Recipe

1. Fit $a_n$ for $\phi_m$ in chosen basis and map.
2. Choose $r_n$ form (exponential is robust).
3. Compute spectral $Pr_t(\zeta)$ on calibration grid.
4. If needed for speed, fit constrained Pad\'e surrogate.
5. Validate monotonicity, no-pole condition, and Richardson consistency.

This yields a physically interpretable closure with explicit controls for neutral limit, knee placement, and very-stable asymptotics.

---

## 8. Repeat the Same Procedure for $Ri_g$

Yes. The same transformed-spectral workflow can be applied directly to the Richardson mapping.

Start from

$$
Ri_g(\zeta)=\zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=\zeta\,\frac{Pr_t(\zeta)}{\phi_m(\zeta)}.
$$

Given the spectral forms in Sections 1 and 2,

$$
\phi_m(\zeta)=\sum_{n=0}^{N} a_n C_n^{(\lambda)}(\xi),
\qquad
\phi_h(\zeta)=Pr_0\sum_{n=0}^{N} r_n a_n C_n^{(\lambda)}(\xi),
$$

the Richardson mapping is obtained algebraically without introducing a new empirical curve.

### 8.1 Why this helps

1. Consistency: $Pr_t$, $Ri_g$, and $\zeta$ inversion all share the same basis and map.
2. Stability control: monotonicity of $Ri_g(\zeta)$ is checked at the same place as $Pr_t$ constraints.
3. Closure economy: one momentum expansion plus one filter gives both Prandtl and Richardson behavior.

### 8.2 Practical surrogate

For fast inversion, fit a constrained Pad\'e map directly for $Ri_g(\zeta)$:

$$
Ri_g^{[2/2]}(\zeta)=\frac{a_1\zeta+a_2\zeta^2}{1+b_1\zeta+b_2\zeta^2},
$$

with constraints:

1. $Ri_g(0)=0$,
2. $Ri_g'(0)>0$,
3. no poles on $\zeta\ge 0$,
4. enforce either finite asymptote ($Ri_c$) or no-asymptote growth, depending on chosen physics.

This mirrors the Prandtl surrogate workflow in Section 5.

---

## 9. Generalization to Other Scalars ($q$, water vapor, $CO_2$, tracers)

Use one shared momentum backbone and scalar-specific filters.

Let scalar $s$ denote temperature, humidity, water vapor, $CO_2$, etc. Define

$$
\phi_s(\zeta)=S_{0,s}\sum_{n=0}^{N} r_{n,s}\,a_n\,C_n^{(\lambda)}\!\bigl(\xi(\zeta)\bigr),
$$

where $S_{0,s}$ is the neutral scalar-transfer ratio and $r_{n,s}$ is scalar-specific damping.

Then scalar Schmidt/Prandtl analogs are

$$
Sc_{t,s}(\zeta)\equiv\frac{\phi_s}{\phi_m}
=S_{0,s}\,\frac{\sum r_{n,s}a_nC_n^{(\lambda)}(\xi)}{\sum a_nC_n^{(\lambda)}(\xi)}.
$$

### 9.1 Scalar hierarchy expectation

If small-scale scalar variance is increasingly suppressed from momentum to heat to passive tracers, one expects

$$
r_{n,m}(=1) \;\ge\; r_{n,h} \;\ge\; r_{n,q} \;\ge\; r_{n,CO_2}
$$

for $n\ge 1$ in strongly stable regimes (site and regime dependent).

Interpretation: the more aggressively filtered scalar gets an earlier/steeper hockey stick in $Sc_{t,s}(\zeta)$.

### 9.2 Minimal parameterization

A compact cross-scalar model is

$$
r_{n,s}=\exp\!\left(-\frac{n}{n_{c,s}}\right),
$$

with scalar-specific $n_{c,s}$ and shared $(a_n,\lambda,\alpha)$.

Recommended calibration order:

1. Fit momentum $(a_n,\lambda,\alpha)$ from wind/shear data.
2. Fit heat $(S_{0,h},n_{c,h})$.
3. Fit humidity/water vapor $(S_{0,q},n_{c,q})$.
4. Fit $CO_2$ or other tracers $(S_{0,s},n_{c,s})$ where data permit.

### 9.3 Why this is better than separate curve fits

1. Shared physics: all scalars inherit the same dynamical backbone.
2. Fewer ad hoc coefficients: scalar differences enter through interpretable filter strength.
3. Easier model integration: one transformation architecture, multiple scalar outputs.

This gives a general multi-scalar closure family with explicit control of neutral limits, transition knees, and very-stable asymptotics for each transported quantity.

---

## 10. Special Solvable ABL Case: UBL with $b=16$ and $\phi_h=\phi_m^2$

Yes. This is an important analytical anchor case, not just a historical parameter choice.

Take the unstable branch ($\zeta<0$) with matched canonical parameters

$$
b_m=b_h=b=16, \qquad a_h=1,
$$

so

$$
\phi_m(\zeta)=(1-b\zeta)^{-1/4},
\qquad
\phi_h(\zeta)=(1-b\zeta)^{-1/2}=\phi_m(\zeta)^2.
$$

### 10.1 Exact consequences

1. Exact Richardson identity:

$$
Ri_g(\zeta)=\zeta\frac{\phi_h}{\phi_m^2}=\zeta.
$$

No iterative inversion is needed: $\zeta=Ri_g$ exactly.

2. Exact Prandtl relation:

$$
Pr_t(\zeta)=\frac{\phi_h}{\phi_m}=\phi_m=(1-b\zeta)^{-1/4}.
$$

So $Pr_t$ is itself a generating-function object (Gegenbauer, $\lambda=1/4$).

3. Spectral closure collapses to convolution identity:

$$
C^{(1/2)}=C^{(1/4)}\otimes C^{(1/4)}.
$$

This means heat structure is not independently fitted; it is determined by momentum through convolution.

### 10.2 Why this is enlightening

This branch gives a "ground-truth" benchmark for model architecture:

1. Any closure aiming to represent matched-exponent UBL should recover $Ri_g=\zeta$ in this limit.
2. Any surrogate ($Pr_t$ Pad\'e, $Ri_g$ Pad\'e, neural map) can be constrained to reproduce this identity at calibration points.
3. Deviations from the identity become physically interpretable diagnostics (parameter mismatch, scalar filtering, intermittency effects).

### 10.3 Exact local expansions and asymptotic form

For $|b\zeta|<1$ (Taylor region):

$$
\phi_m=1+\frac{1}{4}b|\zeta|+\frac{5}{32}(b|\zeta|)^2+\cdots,
$$

$$
\phi_h=1+\frac{1}{2}b|\zeta|+\frac{3}{8}(b|\zeta|)^2+\cdots.
$$

For strong instability ($|\zeta|\gg 1/b$):

$$
\phi_m\sim (b|\zeta|)^{-1/4},
\qquad
\phi_h\sim (b|\zeta|)^{-1/2},
\qquad
Pr_t\sim (b|\zeta|)^{-1/4}.
$$

These give exact slope checks for numerics and clean asymptotic targets for reduced-order parameterizations.

### 10.4 Exact integral forms for profile diagnostics

Many layer diagnostics can be expressed in closed form under this branch because powers are rational.
For example, with $u=1-b\zeta$:

$$
\int (1-b\zeta)^{-1/4}\,d\zeta = -\frac{4}{3b}(1-b\zeta)^{3/4}+C,
$$

$$
\int (1-b\zeta)^{-1/2}\,d\zeta = -\frac{2}{b}(1-b\zeta)^{1/2}+C.
$$

So profile-level and bulk-interval checks can be done analytically in this special case before moving to generalized closures.

In short: yes, this UBL interval is a genuine special case where one gets exact identities, exact inversion, and useful closed-form integrals. It is a strong reference solution for evaluating any generalized Prandtl/Richardson/scalar framework.

---

## 11. What Can Be Determined Analytically, and What Needs Data?

The Gegenbauer / spherical analysis fixes more than a conventional curve fit does, but it does not determine everything.

### 11.1 Fixed or strongly constrained by analysis

1. **Backbone structure:** use one momentum expansion and derive heat/scalar behavior from filtered versions of the same modes.
2. **Anchor identities:** in the matched unstable case ($a_h=1$, $b_m=b_h=16$),

$$
\phi_h=\phi_m^2, \qquad Ri_g=\zeta, \qquad Pr_t=\phi_m.
$$

3. **Basis choice:** Legendre ($\lambda=1/2$) and Gegenbauer ($\lambda\ne 1/2$) are not arbitrary regressors; they encode how spectral resolution is distributed across the neutral-to-strongly-stable interval.
4. **Map shape:** the transformed coordinate $\xi(\zeta)=\tanh(\alpha\zeta)$ can be constrained by desired transition width and boundedness.

### 11.2 Usually needs data

1. Neutral scalar ratios ($Pr_0$, $S_{0,q}$, etc.).
2. Filter scales ($n_c$ for heat, $n_{c,q}$ for humidity, $n_{c,CO_2}$ for tracers).
3. Site- or regime-dependent departures in intermittency, wave activity, cloud influence, and roughness effects.

So the best workflow is not theory *or* data; it is theory first, then data for the remaining degrees of freedom.

### 11.3 Could ML help?

Yes, but in a narrow role.

Recommended use of ML:

1. Predict $n_c$ (or scalar-specific $n_{c,s}$) from metadata: roughness, season, cloud state, stratification class, intermittency indicators.
2. Learn residual corrections after fitting the constrained spectral model.
3. Classify regimes before applying different parameter subsets.

What ML should not do here:

1. Replace the anchor identities with an unconstrained black box.
2. Ignore monotonicity/no-pole constraints in $Pr_t(\zeta)$ or $Ri_g(\zeta)$.
3. Fit each scalar independently with no shared momentum backbone.

---

## 12. Practical Recommendation

For now, the best next step is:

1. Use the new prototype in [code/spectral_scalar_closure.py](../code/spectral_scalar_closure.py) as the calibration sandbox.
2. Fit $n_c$ for heat first using available field datasets.
3. Test whether humidity and other scalars can share $(a_n,\lambda,\alpha)$ with only scalar-specific $(S_{0,s},n_{c,s})$.
4. Only after that decide what belongs in production Fortran or WRF-facing code.

That keeps the mathematically exact anchor cases intact while letting observations determine the genuinely unknown transport constants.
