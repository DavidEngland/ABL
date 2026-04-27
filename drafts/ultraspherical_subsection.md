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

**Differential curvature ($b_m \neq b_h$).**
When $a_h = 1$ but $b_m \neq b_h$, the ratio $(1-b_m\zeta)/(1-b_h\zeta)$ departs from unity and introduces a $\zeta$-dependent distortion of the $Ri_g$ scale.
On the unstable branch ($\zeta < 0$, $\eta = -\zeta > 0$),

$$
Ri_g = -\eta\left(\frac{1+b_m\eta}{1+b_h\eta}\right)^{1/2}.
$$

If $b_m > b_h$ (momentum more strongly stratified than heat; e.g.\ Businger $b_m = 15, b_h = 9$; Hogstrom $b_m = 19, b_h = 11.6$), then $(1+b_m\eta)/(1+b_h\eta) > 1$ for all $\eta > 0$, so

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
