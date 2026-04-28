## 3.x Ultraspherical structure of the momentum similarity operator

### 3.x.1 Generating functions and polynomial families

The Businger–Dyer power-law similarity functions [@Businger_1971; @Dyer1974],

$$
\phi_m(\zeta) = (1 - b_m \zeta)^{-1/4}, \qquad
\phi_h(\zeta) = a_h^{-1}(1 - b_h \zeta)^{-1/2}, \qquad \zeta < 0,
$$

are generating functions of classical orthogonal polynomial sequences.
Expanding in a Taylor series about $\zeta = 0$ (absorbing $a_h^{-1}$ into the $n = 0$ term),

$$
\phi_m(\zeta) = \sum_{n=0}^{\infty} c_n^{(m)}\, \zeta^n, \qquad
c_n^{(m)} = \frac{b_m^n}{n!}\prod_{k=0}^{n-1}\!\left(\tfrac{1}{4}+k\right),
$$

identifies $\phi_m$ with the **Gegenbauer (ultraspherical) generating function** evaluated at the equator $x = 0$ with parameter $\lambda = 1/4$.
Likewise, the $a_h^{-1} = 1$ component of $\phi_h$ corresponds to $\lambda = 1/2$, the **Legendre** special case: $C_n^{(1/2)}(0) \equiv P_n(0)$.
Because $P_{2n}(0) = (-1)^n \binom{2n}{n}/4^n$, the scalar-flux Taylor coefficients are (signed) central binomial coefficients—a result that is exact and parameter-independent.

### 3.x.2 The Sturm–Liouville operator

The Gegenbauer polynomials $C_n^{(\lambda)}(x)$ are eigenfunctions of the self-adjoint **ultraspherical Sturm–Liouville operator**

$$
\mathcal{L}_\lambda[y] \;\equiv\;
-\frac{d}{dx}\!\left[(1-x^2)^{\lambda+1/2}\frac{dy}{dx}\right]
\;=\; n(n+2\lambda)\,(1-x^2)^{\lambda-1/2}\,y,
$$

or equivalently (in non-divergence form),

$$
(1-x^2)\,y'' - (2\lambda+1)\,x\,y' + n(n+2\lambda)\,y = 0,
$$

with orthogonality weight $w_\lambda(x) = (1-x^2)^{\lambda - 1/2}$ on $[-1,1]$.

The **Legendre case $\lambda = 1/2$** is the well-known special instance,

$$
(1-x^2)\,y'' - 2x\,y' + n(n+1)\,y = 0,
\qquad w_{1/2}(x) = 1,
$$

with uniform (Lebesgue) measure—the simplest possible self-adjoint singular-endpoint problem.
The **momentum case $\lambda = 1/4$** retains the same algebraic structure but with singular weight

$$
w_{1/4}(x) = (1-x^2)^{-1/4},
$$

which is integrable ($\in L^1[-1,1]$) yet concentrates measure near the endpoints.
Under the affine map from $\zeta$ to the spectral variable $x$, these endpoints correspond to the neutral limit ($\zeta \to 0$) and to the branch singularity ($\zeta \to b_m^{-1}$) of $\phi_m$.
The momentum operator therefore weights the near-neutral and near-critical stability regimes more heavily than the heat operator does—a spectral encoding of the physical observation that wind profiles respond more smoothly to stability changes than scalar profiles.

Other members of the Gegenbauer family arise at different exponents and serve as useful reference points.
The limiting case $\lambda \to 0$ yields **Chebyshev polynomials of the first kind** $T_n$, with the most singular weight $w_0 = (1-x^2)^{-1/2}$; this limit would correspond to a hypothetical similarity function $\phi \sim (1-b\zeta)^{0}$, i.e., a log profile with no buoyancy correction.
The case $\lambda = 1$ gives **Chebyshev polynomials of the second kind** $U_n$, weight $w_1 = (1-x^2)^{1/2}$, which could represent a profile with exponent $-1$.
No empirically supported MOST scheme falls at these extremes, but they bracket the physically relevant range $\lambda \in (0, 1/2]$.

**The Sturm–Liouville operator in stability space.**
The Legendre change of variables $x = b_h\zeta$ maps the standard interval $[-1,1]$ to the stability interval $\zeta \in [-1/b_h, 1/b_h]$, pulling the scalar Legendre operator back to (see also hw/Legendre.md, Problem 6)

$$
\mathcal{L}_h \;=\; \frac{d}{d\zeta}\!\left[(b_h^{-2} - \zeta^2)\frac{d}{d\zeta}\right],
$$

whose eigenfunctions are $P_n(b_h\zeta)$ with eigenvalues $-n(n+1)$, and whose domain endpoints $\zeta = \pm 1/b_h$ are the branch singularities of $\phi_h$.  The weight function for self-adjointness is $w(\zeta) = 1$—uniform in stability space, consistent with the Legendre ($\lambda = 1/2$) result above.  For momentum, the analogous map $x = b_m\zeta$ gives

$$
\mathcal{L}_m \;=\; \frac{d}{d\zeta}\!\left[(b_m^{-2} - \zeta^2)^{5/4}\frac{d}{d\zeta}\right],
$$

now with singular weight $w_m(\zeta) = (b_m^{-2}-\zeta^2)^{-1/4}$ (from $\lambda = 1/4$).  The operator $\mathcal{L}_h$ is exact for the heat similarity problem; $\mathcal{L}_m$ is the corresponding but non-standard operator for momentum.  Both operators become singular (in the Weyl limit-point sense) at the branch-point endpoints, which is the analytic origin of the instability cutoff in the similarity functions.  A coarse numerical grid that places no level near $\zeta = \pm 1/b$ effectively ignores the endpoint singularity and misrepresents the eigenvalue structure—one spectral-operator perspective on the grid-scale curvature bias documented in stable and very-unstable ABL simulations.

A practical consequence of the singular weight at $\lambda = 1/4$ is spectral accuracy: a truncated Gegenbauer series of degree $N$ approximates $\phi_m$ more faithfully than a Legendre series of the same degree approximates $\phi_h$.
This is consistent with the well-documented empirical finding that stable-boundary-layer models reproduce wind shear more accurately than scalar gradients [@Hogstrom1988; @Grachev_2012b].

### 3.x.3 The squaring identity and its spectral interpretation

When $a_h = 1$ and $b_m = b_h \equiv b$ [@Dyer1974], the two functions satisfy the exact identity

$$
\boxed{\phi_h(\zeta) = \phi_m(\zeta)^2.}
$$

In the Gegenbauer framework this is a **spectral convolution (Clebsch–Gordan) relation**: squaring the $\lambda = 1/4$ generating function yields the $\lambda = 1/2$ generating function,

$$
\left[(1-b\zeta)^{-1/4}\right]^2 = (1-b\zeta)^{-1/2},
$$

because the product of two generating functions with parameter $\lambda$ has parameter $2\lambda$.
Mode-by-mode this corresponds to the linearization formula for ultraspherical polynomials,

$$
C_k^{(1/4)}(0)\,C_{n-k}^{(1/4)}(0)
\;\longrightarrow\; C_n^{(1/2)}(0) = P_n(0),
$$

so the heat-flux spectrum is completely determined by convolution of the momentum spectrum with itself.
The complete algebraic chain is:

$$
\underbrace{b_m=b_h \;\Rightarrow\; \phi_h = \phi_m^2}_{\text{MOST identity} / \, Ri_g=\zeta}
\;\equiv\;
\underbrace{C^{(1/2)} = C^{(1/4)} \otimes C^{(1/4)}}_{\text{Gegenbauer product formula}}
$$

Any consistent closure that employs the canonical exponents $(1/4, 1/2)$ must respect this spectral coupling; it cannot adjust $\phi_m$ and $\phi_h$ independently without breaking the algebraic structure of the ultraspherical family.

### 3.x.4 General Richardson mapping: Prandtl number and the inequality $Ri_g \neq \zeta$

The gradient Richardson number is defined by

$$
Ri_g(\zeta) = \zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}.
$$

Substituting the general power-law forms gives

$$
Ri_g(\zeta) = a_h^{-1}\,\zeta
\left(\frac{1 - b_m\zeta}{1 - b_h\zeta}\right)^{1/2},
$$

which collapses to $Ri_g = \zeta$ only when simultaneously $a_h = 1$ and $b_m = b_h$.
Two independent mechanisms can therefore cause $Ri_g \neq \zeta$:

**Turbulent Prandtl number offset ($a_h \neq 1$).**
The neutral-limit turbulent Prandtl number is $Pr_t^0 = \phi_h(0)/\phi_m(0) = a_h^{-1}$.
Businger et al. [@Businger_1971] report $Pr_t^0 = 1/0.74 \approx 1.35$, so even at $\zeta = 0$ the Richardson mapping is scaled by a constant factor $a_h^{-1} \neq 1$.
This offset persists for all $\zeta$ and shifts the entire $Ri_g$–$\zeta$ curve uniformly in the near-neutral limit.
It is important to note, however, that Businger et al. [@Businger_1971] calibrated their functions using the von Kármán constant $\kappa = 0.35$, whereas the now-standard value is $\kappa = 0.40$.
Högström [@Hogstrom1988] reanalyzed the Kansas data with $\kappa = 0.40$ and obtained $b_m \approx 19.3$, $b_h \approx 11.6$, and $a_h \approx 0.95 \approx 1$, substantially reducing the apparent neutral Prandtl-number departure.
A further refinement advocated by Foken [@Foken2006] and others places $\kappa = 0.41$; at this value the parameter rescaling is minor ($\sim 2.5\%$ relative to $\kappa = 0.40$) but non-negligible for precision flux estimates, and the spread $\kappa \in \{0.40, 0.41\}$ contributes a corresponding uncertainty of roughly $\pm 2\%$ in the inferred $b$ values.
Consequently, the Businger et al. large-$Pr_t^0$ result should be understood as partly a $\kappa$-convention artifact rather than independent confirmation of a strongly non-unity neutral Prandtl number.

**Differential curvature ($b_m \neq b_h$).**
When $a_h = 1$ but $b_m \neq b_h$, the ratio $(1-b_m\zeta)/(1-b_h\zeta)$ departs from unity and introduces a $\zeta$-dependent distortion of the $Ri_g$ scale.
On the unstable branch ($\zeta < 0$, $\eta = -\zeta > 0$),

$$
Ri_g = -\eta\left(\frac{1+b_m\eta}{1+b_h\eta}\right)^{1/2}.
$$

If $b_m > b_h$ (momentum more strongly stratified than heat; e.g. Businger $b_m = 15, b_h = 9$; Högström $b_m = 19, b_h = 11.6$), then $(1+b_m\eta)/(1+b_h\eta) > 1$ for all $\eta > 0$, so

$$
|Ri_g| > |\zeta| \qquad (b_m > b_h,\; a_h = 1,\; \zeta < 0).
$$

The divergence grows with $|\zeta|$ and both effects compound when $a_h^{-1} > 1$, as in the original Businger et al.\ scheme.
Only the Dyer/Brutsaert choice $b_m = b_h = 16$, $a_h = 1$ restores exact equality and eliminates the need for iterative $\zeta \leftrightarrow Ri_g$ inversion [@Webb1970].

The full turbulent Prandtl number as a function of stability is

$$
Pr_t(\zeta) = \frac{\phi_h(\zeta)}{\phi_m(\zeta)} = a_h^{-1}
\frac{(1-b_h\zeta)^{-1/2}}{(1-b_m\zeta)^{-1/4}}
= a_h^{-1}(1-b_m\zeta)^{1/4}(1-b_h\zeta)^{-1/2}.
$$

Under the squaring condition ($a_h = 1$, $b_m = b_h$) this simplifies to $Pr_t(\zeta) = \phi_m(\zeta) = (1-b\zeta)^{-1/4}$, so the Prandtl number itself is a Gegenbauer generating function—momentum and scalar transport become spectrally identical up to a convolution.
When $b_m \neq b_h$ the two families decouple, $Pr_t$ acquires $\zeta$-dependence, and the $Ri_g$–$\zeta$ relationship requires iterative solution or a surrogate inversion scheme.
These conditions—$a_h = 1$, $b_m = b_h$—therefore represent more than a convenient parameterization choice; they are the precise algebraic conditions under which MOST admits the ultraspherical squaring structure described above.

### 3.x.5 Critical Richardson numbers from the similarity functions

Each similarity function carries its own natural stability scale set by its singularity or asymptotic behavior.  Three distinct critical Richardson numbers appear in practice.

**Unstable-branch convergence limits.**
The Taylor series for $\phi_m$ and $\phi_h$ (Section 3.x.1) converge absolutely only within the disk of radius equal to the distance from the origin to the nearest singularity.  For $(1-b_m\zeta)^{-1/4}$ the singularity is at $\zeta = 1/b_m$; for $(1-b_h\zeta)^{-1/2}$ it is at $\zeta = 1/b_h$.  These singularities lie on the stable side ($\zeta > 0$), but they define a natural scale that, mapped through the $Ri_g$--$\zeta$ relation, yields

$$
Ri_{c,m}^{\mathrm{UBL}} = -\frac{1}{b_m}, \qquad
Ri_{c,h}^{\mathrm{UBL}} = -\frac{1}{b_h}.
$$

Numerically: with Dyer parameters ($b_m = b_h = 16$), $Ri_{c,m}^{\mathrm{UBL}} = Ri_{c,h}^{\mathrm{UBL}} = -0.063$; with Businger parameters ($b_m = 15$, $b_h = 9$), $Ri_{c,m}^{\mathrm{UBL}} = -1/15 \approx -0.067$ and $Ri_{c,h}^{\mathrm{UBL}} = -1/9 \approx -0.111$; with Högström parameters ($b_m = 19.3$, $b_h = 11.6$), $Ri_{c,m}^{\mathrm{UBL}} \approx -0.052$ and $Ri_{c,h}^{\mathrm{UBL}} \approx -0.086$.  These values mark the onset of strong nonlinearity in the unstable similarity functions and represent the effective upper bound on the unstable branch for spectral-series representations.

**Stable-branch critical Richardson number.**
On the stable branch the linear forms $\phi_m = 1+\beta_m\zeta$ and $\phi_h = a_h^{-1}+\beta_h\zeta$ are unbounded as $\zeta\to\infty$, but the Richardson number approaches a finite asymptote.  For $a_h = 1$ and $\beta_m = \beta_h \equiv \beta$ (degenerate case),

$$
Ri_g(\zeta) = \zeta\,\frac{1+\beta\zeta}{(1+\beta\zeta)^2}
= \frac{\zeta}{1+\beta\zeta}
\;\longrightarrow\; \frac{1}{\beta} \equiv Ri_c^{\mathrm{SBL}}
\quad (\zeta\to\infty).
$$

This asymptotic value is a hard upper bound: the linear MOST closure predicts that $Ri_g$ cannot exceed $1/\beta$ no matter how large $\zeta$ becomes—a theoretical ceiling on turbulent stability in the linearly stratified regime.
Empirical choices for $\beta$ lead to different $Ri_c^{\mathrm{SBL}}$ values:

| Source | $\kappa$ | $\beta_m$ | $\beta_h$ | $Ri_c^{\mathrm{SBL}}$ | Notes |
|---|---:|---:|---:|---:|---|
| Businger et al. (1971) [@Businger_1971] | 0.35 | 4.7 | 6.35 | $1/4.7 \approx 0.213$ | $\kappa$-corrected $\to$ Högström |
| Dyer (1974) [@Dyer1974] | 0.40 | 5.0 | 5.0 | $1/5 = 0.200$ | Common WRF default |
| Högström (1996) [@Hogstrom1988] | 0.40 | 5.3 | 5.3 | $1/5.3 \approx 0.189$ | Kansas re-analysis |
| McNider et al. | 0.40 | 4.0 | 4.0 | $1/4 = 0.250$ | Used in stable-SBL studies |

The spread $0.19 \lesssim Ri_c^{\mathrm{SBL}} \lesssim 0.25$ reflects genuine empirical uncertainty amplified by the $\kappa$-convention issue: the Businger stable slope $\beta_m = 4.7$ measured at $\kappa = 0.35$ would shift somewhat under $\kappa = 0.40$ renormalization.  The theoretical "classical" value $Ri_c = 0.25$ from Kolmogorov/Richardson turbulence theory coincides with the McNider choice $\beta = 4$; the Dyer value $\beta = 5$ gives $Ri_c = 0.20$, consistent with the stricter observational estimates of Grachev et al.\ [@Grachev_2012b].  Neither the Gegenbauer structure nor the $Ri_g = \zeta$ identity constrain $\beta$ directly; the stable-branch critical Richardson number is an independent empirical input.  The relationship between the unstable-branch curvature coefficient $b$ and the stable-branch slope $\beta$ is discussed further in Section~3.x.8.

### 3.x.6 Asymptotic expansions of the similarity functions for large instability

For strongly unstable conditions, $\eta = -\zeta \gg 1/b$, the series representations of Section 3.x.1 diverge and must be replaced by the large-$\eta$ asymptotic expansions.  Writing $u = b_m\eta \gg 1$ and using the binomial series $(1+u)^{-\alpha} = u^{-\alpha}\sum_{k\geq 0}\binom{-\alpha}{k}u^{-k}$, the momentum function expands as

$$\phi_m(-\eta) = (b_m\eta)^{-1/4}\Bigl[1 - \frac{1}{4}(b_m\eta)^{-1} + \frac{5}{32}(b_m\eta)^{-2} - \frac{15}{128}(b_m\eta)^{-3} + \mathcal{O}\bigl((b_m\eta)^{-4}\bigr)\Bigr].$$

and the heat function as

$$\phi_h(-\eta) = (b_h\eta)^{-1/2}\Bigl[1 - \frac{1}{2}(b_h\eta)^{-1} + \frac{3}{8}(b_h\eta)^{-2} - \frac{5}{16}(b_h\eta)^{-3} + \mathcal{O}\bigl((b_h\eta)^{-4}\bigr)\Bigr].$$

Under the degenerate constraint $b_m = b_h = b$, these asymptotic series preserve the same product structure as the exact MOST identity. Squaring the momentum expansion gives

$$\phi_m(-\eta)^2 = (b\eta)^{-1/2}\Bigl[1 - \frac{1}{2}(b\eta)^{-1} + \frac{3}{8}(b\eta)^{-2} - \frac{5}{16}(b\eta)^{-3} + \mathcal{O}\bigl((b\eta)^{-4}\bigr)\Bigr] = \phi_h(-\eta).$$

so the first correction coefficient doubles from $-1/4$ to $-1/2$, exactly as required by the Gegenbauer product relation $C^{(1/2)} = C^{(1/4)} \otimes C^{(1/4)}$.

The leading terms $\phi_m \sim (b_m\eta)^{-1/4}$ and $\phi_h \sim (b_h\eta)^{-1/2}$ are the **free-convection scaling laws**: gradient Richardson-based turbulence transitions to buoyancy-dominated convection, and wind shear decays as $\eta^{-1/4}$ while scalar gradients decay more rapidly as $\eta^{-1/2}$.  These exponents are the direct spectral consequence of the Gegenbauer parameters $\lambda = 1/4$ and $\lambda = 1/2$, respectively.

Several implications follow from retaining the first correction term:

1. **Asymptotic Richardson mapping.** For $\eta \gg 1/b_m, 1/b_h$, the ratio $\phi_h/\phi_m^2$ approaches
$$
\frac{\phi_h}{\phi_m^2} \to \frac{(b_h\eta)^{-1/2}}{(b_m\eta)^{-1/2}} = \left(\frac{b_m}{b_h}\right)^{1/2},
$$
so $Ri_g \sim -\eta\,(b_m/b_h)^{1/2}$ as $\eta\to\infty$.  Only when $b_m = b_h$ does $Ri_g \to -\eta = \zeta$ in this limit, confirming the squaring identity holds asymptotically as well as locally.

2. **Onset of the asymptotic regime.** The first correction is $O((b\eta)^{-1})$ and is less than $5\%$ when $b\eta > 20$.  For Dyer $b = 16$ this requires $\eta > 1.25$ ($|\zeta| > 1.25$), a moderately strong instability.  At $\eta = 0.5$ the correction is $\sim 12\%$ for momentum and $\sim 6\%$ for heat.

3. **Accuracy crossover.** The Taylor series (Section 3.x.1) converges for $b\eta < 1$ and the asymptotic series is accurate for $b\eta > 5$.  In the intermediate range $1 < b\eta < 5$ neither representation is highly efficient and the exact power-law evaluation $(1-b\zeta)^{-\alpha}$ should be used directly—as all operational codes already do.  However, the asymptotic expansions provide a valuable check and a basis for analytical estimates of flux-profile integrals $\psi_m$, $\psi_h$ in the free-convection limit.

4. **Differential asymptotic decay.** The exponent $-1/2$ for $\phi_h$ is twice the magnitude of $-1/4$ for $\phi_m$, so $\phi_h$ decays to zero faster than $\phi_m$ in the free-convection limit.  This means scalar gradients vanish more rapidly than momentum gradients under strong convection—a well-documented observational result that is here seen as a direct consequence of $\lambda_h = 2\lambda_m$.

### 3.x.7 Numerical implications for simulation

The ultraspherical and asymptotic structure of the similarity functions has several practical consequences for numerical implementation.

**Evaluation strategy.** The exact power-law forms $(1-b\zeta)^{-1/4}$ and $(1-b\zeta)^{-1/2}$ are trivially cheap to evaluate and should always be used in production code. The spectral representations inform accuracy and error control:

- For $|b\zeta| < 0.5$: the $N$-term Taylor/CBC series converges with relative error $< (b|\zeta|)^N/(1-b|\zeta|)$. A 4-term series achieves $< 10^{-3}$ relative error throughout this range.
- For $b\eta > 5$: the 3-term asymptotic expansion is accurate to $< 0.1\%$. If only an approximate $\phi_h$ is needed (e.g., to seed a Newton iteration), the leading asymptote $(b_h\eta)^{-1/2}$ suffices.

**Richardson inversion.** The function $\zeta(Ri_g)$ requires iterative solution in general. The spectral structure suggests a two-branch strategy:

- *Degenerate case* ($a_h = 1$, $b_m = b_h$): use $\zeta = Ri_g$ directly — no iteration.
- *General case*: initialize Newton's method at $\zeta^{(0)} = Ri_g \bigl/ \bigl[a_h^{-1}(1+(b_m-b_h)|Ri_g|/2)^{1/2}\bigr]$ (first-order correction from the ratio formula in Section 3.x.4). This typically converges in 2--3 iterations.
- *Strongly unstable* ($Ri_g < -0.1$): use the asymptotic inverse $\eta \approx |Ri_g|(b_h/b_m)^{1/2}$ as the initial guess.

**Stable-branch singular behavior and the SL operator.** The operator $\mathcal{L}_h$ (Section 3.x.2) is singular at $\zeta = \pm 1/b_h$. In a discretized model, the numerical grid implicitly defines a finite-element approximation to this operator; grid points near the singularity are needed to represent the rapid changes in the eigenfunction structure near $|\zeta| \approx 1/b$. Grids that are geometrically stretched toward the surface (standard log-linear levels) do not align with the Gauss--Legendre nodes of $\mathcal{L}_h$, introducing a systematic bias that is largest in the moderately stable and moderately unstable regimes where $|b\zeta| \sim 0.3$--$0.7$.

**Parameter sensitivity.** Because $Ri_c^{\mathrm{UBL}} = -1/b$ and $Ri_c^{\mathrm{SBL}} = 1/\beta$, a $\pm 10\%$ uncertainty in $b$ or $\beta$ propagates directly to a $\pm 10\%$ uncertainty in the corresponding critical Richardson number, affecting the regime boundaries and any parameterization that conditions on $Ri_g$ crossing those thresholds.

### 3.x.8 Dynamic critical Richardson number: coupling between $b$ and $\beta$

The framework above reveals a natural algebraic constraint between the unstable-branch coefficient $b$ and the stable-branch slope $\beta$ that has not been widely exploited.

**Neutral slope matching.**  The derivative of the unstable power law at $\zeta = 0^-$ is

$$
\left.\frac{d\phi_m}{d\zeta}\right|_{\zeta=0^-} = \frac{b_m}{4},
$$

while the slope of the stable linear branch at $\zeta = 0^+$ is $\beta_m$.  Requiring **$C^1$ continuity** at neutral gives the matching condition

$$
\boxed{\beta = \frac{b}{4}.}
$$

Applied to standard parameter sets:

| Source | $b_m$ | $b_m/4$ | $\beta_m$ (empirical) | Match? |
|---|---:|---:|---:|---|
| Businger (1971, $\kappa=0.35$) | 15.0 | 3.75 | 4.7 | No ($\Delta \approx 25\%$) |
| Dyer (1974) | 16.0 | **4.0** | 5.0 | No ($\Delta = 25\%$) |
| Högström (1996) | 19.3 | **4.8** | 5.3 | Nearly ($\Delta \approx 10\%$) |
| McNider et al. | 16.0 | **4.0** | **4.0** | **Exact** |

The McNider choice $\beta = 4$ with $b = 16$ is the **only standard parameter set satisfying the neutral slope-matching condition**.  This implies $Ri_c^{\mathrm{SBL}} = 4/b$; for $b = 16$, $Ri_c^{\mathrm{SBL}} = 0.25$—the classical theoretical value.  The Dyer stable slope $\beta = 5$ requires $b = 20$ for matching, somewhat larger than the $b = 16$ Dyer unstable fit, explaining the persistent kink at $\zeta = 0$ in Dyer-family closures.

**Dynamic Ri$_c$ implications.**  If the similarity coefficients are allowed to vary with regime (e.g., as functions of local $TKE$, height, or stability itself in a dynamic-$Ri_c$ framework), the spectral structure imposes two consistency requirements:

1. The **ultraspherical squaring identity** $\phi_h = \phi_m^2$ is preserved if and only if $b_m = b_h$ is maintained at each dynamic state.  A dynamic scheme that updates $b_m$ independently of $b_h$ breaks the Gegenbauer convolution structure and necessarily introduces $Ri_g \neq \zeta$, requiring iterative inversion.

2. The **neutral slope-matching constraint** $\beta = b/4$ links the stable-branch critical Richardson number to the unstable-branch coefficient: $Ri_c^{\mathrm{SBL}} = 1/\beta = 4/b$.  If $b$ evolves dynamically (e.g., increasing toward the SHEBA/Arctic stable limit), then $\beta$ and $Ri_c^{\mathrm{SBL}}$ should co-evolve to avoid a spurious kink at neutral.  Concretely, $b \to b(\zeta)$ implies $\beta \to b(0)/4$ for continuity.

3. The **UBL critical value** $Ri_c^{\mathrm{UBL}} = -1/b$ moves with $b$: a larger $b$ (more sensitive momentum function) lowers $|Ri_c^{\mathrm{UBL}}|$, pushing the onset of the free-convection asymptotic regime toward neutral.  In practice this means models with large $b$ (e.g., Högström $b_m = 19.3$) enter the asymptotic regime at weaker instability, which affects the flux-profile integrals $\psi_m$, $\psi_h$ in the moderately unstable range.

Taken together, the ultraspherical framework suggests that a dynamically consistent MOST closure should parameterize a single master coefficient $b = b_m = b_h$ (enforcing the squaring identity), from which all critical Richardson numbers follow: $Ri_c^{\mathrm{UBL}} = -1/b$, $Ri_c^{\mathrm{SBL}} = 4/b$, and $Pr_t(\zeta) = (1-b\zeta)^{-1/4}$.  The remaining empirical freedom lies entirely in the value of $b$—one parameter rather than four.

### 3.x.9 Evaluation policy, safeguarded inversion, and WRF-style implementation

For implementation in production surface-layer and PBL code, we recommend an explicit policy that separates mathematically exact evaluation from acceleration aids:

1. **Small-instability branch** ($b\eta < 0.5$): evaluate $\phi_m,\phi_h$ with CBC/Gegenbauer recurrences.
2. **Transition branch** ($0.5 \le b\eta \le 5$): evaluate exact power laws directly, $(1+b\eta)^{-1/4}$ and $(1+b\eta)^{-1/2}$.
3. **Large-instability branch** ($b\eta > 5$): exact power laws remain primary; asymptotic expressions are used for diagnostics or initial guesses.
4. **Ri-inversion** ($Ri_g \mapsto \zeta$): use safeguarded Newton (Newton + bisection fallback) with branch-aware brackets.

The transition branch is intentionally direct-evaluation only: truncated Taylor and truncated asymptotic forms are both suboptimal there and can degrade monotonicity and inversion robustness.

#### Pseudocode (evaluation and inversion policy)

```text
function evaluate_phi_and_zeta(profile_id, Ri_target, b_m, b_h):
    # Step A: invert Ri_g(zeta) with safeguarded Newton
    zeta = invert_rig_safeguarded(profile_id, Ri_target)

    # Step B: choose evaluation branch by instability magnitude
    if zeta < 0:
        eta = -zeta
        s = max(b_m*eta, b_h*eta)
        if s < 0.5:
            phi_m = gegenbauer_recurrence_phi_m(zeta, b_m)
            phi_h = cbc_recurrence_phi_h(zeta, b_h)
        else:
            phi_m = (1 - b_m*zeta)^(-1/4)
            phi_h = (1 - b_h*zeta)^(-1/2)
    else:
        phi_m, phi_h = stable_profile(profile_id, zeta)

    # Step C: return closure functions
    f_m = clamp(1/(phi_m*phi_m), 0, 1)
    f_h = clamp(1/(phi_m*phi_h), 0, 1)
    return zeta, phi_m, phi_h, f_m, f_h
```

#### Safeguarded flowchart for $Ri_g \mapsto \zeta$

```mermaid
flowchart TD
    A[Input Ri_target, profile_id] --> B{Ri_target <= 0?}
    B -->|Yes| C[Set bracket z_lo=-20, z_hi=0]
    B -->|No| D[Set bracket z_lo=0, z_hi=50]
    C --> E[Evaluate f(z)=Ri_g(z)-Ri_target at bracket ends]
    D --> E
    E --> F{Sign change bracketed?}
    F -->|No| G[Fallback guarded Newton seed z=clamp(Ri_target)]
    F -->|Yes| H[Start hybrid loop]
    H --> I[Compute f(z) and finite-diff fprime(z)]
    I --> J{Newton step valid and in bracket?}
    J -->|Yes| K[Accept Newton candidate]
    J -->|No| L[Use bisection midpoint]
    K --> M[Evaluate f(candidate)]
    L --> M
    M --> N[Update bracket by sign]
    N --> O{abs(f) < tol or bracket < tol?}
    O -->|No| I
    O -->|Yes| P[Return converged zeta]
    G --> Q[Run limited iterations + clamp]
    Q --> R[Return zeta with warning flag]
```

#### WRF-style Fortran integration sketch

```fortran
! In PBL tendency loop, after Ri is computed:
REAL :: rig_local, zeta, phi_m, phi_h, fm_loc, fh_loc
LOGICAL :: ok
INTEGER :: it

rig_local = Ri(i,k)

CALL zeta_from_rig_safeguarded(PROFILE_BD71, rig_local, zeta, ok, it)
IF (.NOT. ok) THEN
   CALL zeta_from_rig_newton(PROFILE_BD71, rig_local, zeta)
END IF

CALL phi_from_zeta(PROFILE_BD71, zeta, phi_m, phi_h)

fm_loc = 1.0 / (phi_m * phi_m)
fh_loc = 1.0 / (phi_m * phi_h)
fm_loc = MAX(0.0, MIN(fm_loc, 1.0))
fh_loc = MAX(0.0, MIN(fh_loc, 1.0))

KM(i,k) = KM(i,k) * fm_loc
KH(i,k) = KH(i,k) * fh_loc
```

In this workflow, robustness is enforced at inversion time (safeguarded flow), while physical consistency is enforced at evaluation time (exact power-law core, recurrence only near neutral).  This separation is especially important for mesoscale operational grids where the transition band $1 < b\eta < 5$ is frequently occupied.
