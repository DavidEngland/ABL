
# Homework: Legendre and Gegenbauer Polynomials in MOST Stability Functions

**Course:** Atmospheric Boundary Layer Physics — Advanced Topics  
**Level:** Master's (Problems 1–3) and PhD (Problems 4–7)  
**Prerequisites:** `hw/CBC.md`; real analysis; orthogonal polynomials (Legendre, Gegenbauer); Python (scipy, numpy)

> **The big picture.**  
> The central binomial coefficient (CBC) expansion of the MOST heat-gradient function  
> $$\phi_h(-\eta) = (1+b_h\eta)^{-1/2} = \sum_{n=0}^{\infty} \binom{2n}{n}\frac{(-1)^n}{4^n}(-b_h\eta)^n$$  
> is **not merely a convenient Taylor series**.  Its coefficients are the values of the  
> Legendre polynomials at the *equatorial* angle $\theta=\pi/2$:  
> $$c_n^{(h)} = (-1)^n P_{2n}(0) \cdot b_h^n.$$  
> The function $\phi_h$ evaluated on the unstable branch is therefore a **generating function  
> for the even Legendre polynomials at the equator**, with the convective variable  
> $t = \sqrt{b_h\,\eta}$ as the parameter.  
>  
> Momentum $\phi_m$ (exponent $1/4$) does the same for **Gegenbauer (ultraspherical)  
> polynomials $C_n^{(1/4)}$** — a strict generalization that collapses to Legendre when  
> the exponent is $1/2$.  
>  
> The squaring identity $\phi_h = \phi_m^2$ is, at the polynomial level, a **Gegenbauer  
> product formula** (linearization identity). Every result in `hw/CBC.md` has an exact  
> counterpart in spherical harmonic theory.

---

## Notation (consistent with `hw/CBC.md` and `code/most.f`)

| Symbol | Meaning |
|---|---|
| $\eta = -\zeta > 0$ | Convective variable (unstable branch) |
| $b_m,\,b_h$ | Unstable Businger-Dyer coefficients |
| $\alpha_m = 1/4,\;\alpha_h = 1/2$ | Exponents |
| $P_n(x)$ | Legendre polynomial (Gegenbauer $C_n^{(1/2)}$) |
| $C_n^{(\lambda)}(x)$ | Gegenbauer (ultraspherical) polynomial, parameter $\lambda$ |
| $\hat{f}_n$ | Legendre spectral coefficient of $f$ |

Generating function identities (proved in any classical analysis text):
$$
\frac{1}{\sqrt{1-2xt+t^2}} = \sum_{n=0}^{\infty} P_n(x)\,t^n, \qquad |t|<1, \quad x\in[-1,1].
$$
$$
\frac{1}{(1-2xt+t^2)^{\lambda}} = \sum_{n=0}^{\infty} C_n^{(\lambda)}(x)\,t^n, \qquad \lambda > 0.
$$

---

## Problem 1 — Legendre Polynomials at the Equator and the CBC Coefficients (MS, 25 pts)

### Part A: Rodrigues' formula at $x = 0$ (10 pts)

Rodrigues' formula:
$$P_n(x) = \frac{1}{2^n n!}\frac{d^n}{dx^n}(x^2-1)^n.$$

1. Expand $(x^2-1)^{2m}$ by the binomial theorem and isolate the $x^{2m}$ term.
2. Differentiate $2m$ times and evaluate at $x=0$ to show that for odd $n$, $P_n(0)=0$,
	and for even $n=2m$:
	$$\boxed{P_{2m}(0) = (-1)^m \frac{\binom{2m}{m}}{4^m}.}$$
3. Verify numerically for $m=0,1,2,3$ using `scipy.special.legendre`.

### Part B: CBC coefficients as Legendre values (5 pts)

From `hw/CBC.md`, the CBC expansion of $\phi_h$ on the unstable branch ($\eta>0$) is
$$\phi_h(-\eta) = (1+b_h\eta)^{-1/2} = \sum_{n=0}^{\infty} c_n^{(h)}(-\eta)^n,
\quad c_n^{(h)} = \binom{2n}{n}\frac{b_h^n}{4^n}.$$

1. Show that the $n$-th coefficient satisfies $c_n^{(h)} = (-1)^n P_{2n}(0)\,b_h^n$.
2. Rewrite the series as
	$$\phi_h(-\eta) = \sum_{n=0}^{\infty} P_{2n}(0)\,(b_h\eta)^n.$$
	Note the signs cancel — the sum over even Legendre values has *positive* terms
	when $\eta>0$.

### Part C: Interpretation as the Legendre generating function (10 pts)

The Legendre generating function at $x=0$ with $t=s$ gives
$$\frac{1}{\sqrt{1+s^2}} = \sum_{n=0}^{\infty} P_{2n}(0)\,s^{2n}$$
(odd terms vanish by parity $P_{2k+1}(0)=0$).

1. Substitute $s = \sqrt{b_h\,\eta}$ to derive the exact identity
	$$\phi_h(-\eta) = G\!\left(\sqrt{b_h\eta},\;0\right),$$
	where $G(t,x)$ is the Legendre generating function.
2. What does it mean geometrically that $\phi_h$ evaluated on the unstable branch is a
	**slice of the Legendre generating function at the equator** ($\theta=\pi/2$,
	i.e., $x=\cos\theta = 0$)?
3. Explain in one paragraph why the *parity* of the relevant Legendre polynomials being
	**even** is connected to the unstable branch being a function of $\eta = |\zeta|$, not $\zeta$.

---

## Problem 2 — Gegenbauer Polynomials and $\phi_m$ (MS, 25 pts)

### Part A: The ultraspherical generalization (10 pts)

The Gegenbauer polynomial $C_n^{(\lambda)}(x)$ satisfies the generating function
$$\frac{1}{(1-2xt+t^2)^{\lambda}} = \sum_{n=0}^{\infty} C_n^{(\lambda)}(x)\,t^n.$$

At $x=0$ (equatorial slice):
$$\frac{1}{(1+t^2)^{\lambda}} = \sum_{n=0}^{\infty} C_{2n}^{(\lambda)}(0)\,t^{2n}, \qquad C_{2n+1}^{(\lambda)}(0) = 0.$$

1. Using the identity $C_{2n}^{(\lambda)}(0) = (-1)^n \dfrac{\Gamma(\lambda+n)}{\Gamma(\lambda)\,n!}$,  
	write the expansion of $(1+b_m\eta)^{-1/4}$ (i.e., $\phi_m(-\eta)$) in terms of
	$C_{2n}^{(1/4)}(0)$.

2. Show that the recurrence for $c_n^{(m)} = |C_{2n}^{(1/4)}(0)|\,b_m^n$ is exactly
	$$c_{n+1}^{(m)} = c_n^{(m)} \cdot \frac{b_m(4n+1)}{4(n+1)},\quad c_0^{(m)}=1,$$
	recovering the `hw/CBC.md` result from the Gamma-function ratio.

### Part B: Legendre as a special case (5 pts)

1. Show that $C_n^{(1/2)}(x) = P_n(x)$ (i.e., Legendre is Gegenbauer with $\lambda=1/2$).  
	*Hint:* Compare the generating functions.
2. Conclude that $\phi_h$ uses **Legendre** polynomials while $\phi_m$ uses
	**Gegenbauer $C_n^{(1/4)}$** — a "softer" ultraspherical family appropriate
	to the weaker ($-1/4$) singular response of momentum relative to heat.

### Part C: Comparing coefficient growth (10 pts)

1. For large $n$, use Stirling's approximation on the Gamma functions to show:
	$$C_{2n}^{(\lambda)}(0) \approx (-1)^n \frac{n^{\lambda-1}}{\Gamma(\lambda)},$$
	so that for $\lambda=1/2$, $|C_{2n}^{(1/2)}(0)| \sim 1/(\sqrt{\pi}\sqrt{n})$, and
	for $\lambda=1/4$, $|C_{2n}^{(1/4)}(0)| \sim 1/(\Gamma(1/4)\,n^{3/4})$.
2. Sketch both envelopes $|C_{2n}^{(\lambda)}(0)| b^n$ vs $n$ for $b=16$ near the
	radius of convergence.  Which family's terms grow faster to the singularity?
3. Connect your result to the physical statement: "scalar gradients respond more
	violently to strong instability than momentum gradients."

---

## Problem 3 — The Squaring Identity as a Gegenbauer Product Formula (MS/PhD, 20 pts)

### Part A: Algebraic proof (5 pts)

Recall the identity (proved in `hw/CBC.md` Problem 3):
when $b_m = b_h = b$,
$$\phi_h(\zeta) = \phi_m(\zeta)^2 \quad \Longleftrightarrow \quad (1-b\zeta)^{-1/2} = \bigl[(1-b\zeta)^{-1/4}\bigr]^2.$$

Verify this is an identity of generating functions: if
$\phi_m^{-1}(t) = C^{(1/4)}$-generating-function at equator, then  
$(\phi_m)^2$ corresponds to the convolution product of Gegenbauer series.

### Part B: Linearization of Gegenbauer products (10 pts)

The **linearization (Clebsch-Gordan) identity** for Gegenbauer polynomials states that
$$C_m^{(\lambda)}(x)\,C_n^{(\lambda)}(x) = \sum_{k=0}^{\min(m,n)} a_{mnk}^{(\lambda)}\,C_{m+n-2k}^{(2\lambda)}(x),$$
where the linearization coefficients $a_{mnk}^{(\lambda)}$ are known in closed form.

1. For the equatorial case $x=0$, simplify the identity using $C_n^{(\lambda)}(0) = 0$
	for odd $n$ to get a sum over even terms only.
2. Show that for $\lambda=1/4$ and the squaring $(m=n=2j)$:
	$$\left[C_{2j}^{(1/4)}(0)\right]^2 = a_{jj}^{(1/4)}\,C_{4j}^{(1/2)}(0) + \text{lower},$$
	i.e., the convolution maps $C^{(1/4)} \otimes C^{(1/4)} \to C^{(1/2)}$, which is the
	polynomial-level statement of $\phi_h = \phi_m^2$.
3. Verify numerically for $j=1$ (i.e., $n=2$): compute $[C_2^{(1/4)}(0)]^2$ and
	the right-hand side explicitly.

### Part C: Physical interpretation (5 pts)

1. Summarize in two sentences: what does the Gegenbauer product formula say about
	the **spectral relationship** between momentum and heat transport in the unstable ABL?
2. If we **perturb** $b_h \neq b_m$, the squaring identity breaks.  What changes in the
	generating-function language?  (A sentence or two is sufficient.)

---

## Problem 4 — Legendre Spectral Expansion of the Stability Profile (PhD, 30 pts)

**Setup:** Map the convergence domain of the CBC series to the standard Legendre interval.
Define the linear map
$$\mu = 1 - 2b_h|\zeta|, \quad |\zeta| = \frac{1-\mu}{2b_h}, \quad \mu \in [-1,1] \leftrightarrow |\zeta| \in [0,1/b_h].$$
Then $\theta$ defined by $\cos\theta = \mu$ ranges from $0$ (neutral) to $\pi$ (branch point).

### Part A: Profile in angular coordinates (5 pts)

1. Under this map, show that
	$$\phi_h\!\left(-\frac{1-\mu}{2b_h}\right) = \sqrt{\frac{2}{3-\mu}}.$$
2. Verify the limits: at $\mu=1$ ($\zeta=0$): $\phi_h = 1$.
	At $\mu = -1$ ($|\zeta| = 1/b_h$, branch point): $\phi_h = \sqrt{1/2}$. Wait — check:
	$3-(-1)=4$, so $\phi_h = \sqrt{2/4} = 1/\sqrt{2}$? Explain why the profile does **not**
	diverge at the branch point under this map (the full real singularity is at $\zeta=1/b_h$
	for the stable branch, and $1/(1+1) = 1/\sqrt{2}$ here).  
	*Note: the apparent mildness is because we are mapping $\zeta \in [-1/b_h,0]$, not crossing the stable-branch singularity.*
3. Repeat for $\phi_m$: show $\phi_m(-\frac{1-\mu}{2b_m}) = \left(\frac{2}{3-\mu}\right)^{1/4}$.

### Part B: Legendre spectral coefficients (15 pts)

Define the Legendre expansion
$$f_h(\mu) = \sqrt{\frac{2}{3-\mu}} = \sum_{n=0}^{\infty} \hat{f}_n^{(h)}\,P_n(\mu),$$
where $\hat{f}_n^{(h)} = \dfrac{2n+1}{2}\displaystyle\int_{-1}^{1} f_h(\mu)\,P_n(\mu)\,d\mu$.

1. Compute $\hat{f}_0^{(h)}$ analytically:
	$$\hat{f}_0^{(h)} = \frac{1}{2}\int_{-1}^{1}\sqrt{\frac{2}{3-\mu}}\,d\mu.$$
	*Hint:* Substitute $u = 3-\mu$, $du = -d\mu$, and evaluate $\int (u)^{-1/2}\,du$.
	Expected answer: $\hat{f}_0^{(h)} = 2\sqrt{2}-2$.

2. Compute $\hat{f}_1^{(h)}$ analytically using the Rodrigues formula to integrate by parts.
	*Hint:* $P_1(\mu) = \mu$; integrate $\mu/\sqrt{3-\mu}$ over $[-1,1]$.

3. Write a Python function `legendre_coeffs(f, N)` using `numpy.polynomial.legendre.leggauss`
	(Gauss-Legendre quadrature) to compute $\hat{f}_n$ numerically for $n=0,\ldots,N-1$.
	Plot the first 10 coefficients for both $f_h$ and $f_m$.

4. Estimate the **Legendre spectral truncation error** after $N$ terms as a function of $N$.
	What spectral decay rate do you observe?  Is it algebraic or exponential?  Why?
	(*Hint:* smoothness of $f_h$ on $[-1,1]$ governs the decay rate.)

### Part C: Physical interpretation of Legendre modes (10 pts)

Each Legendre coefficient $\hat{f}_n^{(h)}$ represents the projection of the stability profile
onto a mode of degree $n$ in the mapped stability space $\mu \in [-1,1]$.

1. The $n=0$ mode is the **mean** (spatially uniform part) of the stability function.
	The $n=1$ mode is the **linear tilt** (large-scale gradient in stability space).
	Interpret the higher modes physically in terms of the nonlinear curvature of the ABL profile.

2. In a numerical model that truncates at $N$ Legendre modes, how many modes are needed
	to represent $f_h$ with relative error below $1\%$ over the full convergence domain?
	What is the analogous requirement for $f_m$?  Which requires more modes?

3. Relate your answer to the **grid-scale bias** results in `MOST_Derivations_Problems.md`
	Problem 3: a coarse grid effectively projects the profile onto its lowest Legendre modes.

---

## Problem 5 — Gauss-Legendre Quadrature and the Bulk Richardson Number (PhD, 25 pts)

**Background.** The bulk Richardson number over a layer $[0,\zeta_*]$ is
$$Ri_b(\zeta_*) = \frac{1}{\zeta_*}\int_0^{\zeta_*} Ri_g(\zeta)\,d\zeta.$$
Under the degenerate case $Ri_g = \zeta$, this gives $Ri_b = \zeta_*/2$.
For the general case, numerical integration is required.

### Part A: GL quadrature setup (10 pts)

Map $\zeta \in [0,\zeta_*]$ to $\mu \in [-1,1]$ via $\zeta = \zeta_*(1+\mu)/2$.

1. Write the integral in the mapped variable:
	$$Ri_b = \frac{1}{2}\int_{-1}^{1} Ri_g\!\left(\frac{\zeta_*(1+\mu)}{2}\right)d\mu.$$

2. Using the $N$-point Gauss-Legendre rule with nodes $\mu_k$ and weights $w_k$
	(from `numpy.polynomial.legendre.leggauss(N)`), write the quadrature approximation.

3. Show that the **GL nodes $\mu_k$ are zeros of $P_N(\mu)$**.  Relate the *degree of
	exactness* of the $N$-point GL rule (exact for polynomials of degree $\leq 2N-1$) to
	the number of CBC terms needed for the same accuracy.

### Part B: Error analysis (10 pts)

1. Compute $Ri_b(\zeta_*)$ numerically for $\zeta_* = -0.5$ (unstable) using $N=2, 4, 8$
	GL nodes and compare to the exact integral (compute via `scipy.integrate.quad`).
2. Plot the relative error vs $N$ on a semilog axis.  Is convergence exponential?
3. For which value of $\zeta_*$ does the GL quadrature converge slowest?
	Relate your answer to the radius of convergence of the CBC series.

### Part C: Connection to Legendre zeros and CBC (5 pts)

The GL nodes $\mu_k = \cos\theta_k$ lie on the unit circle at angles that are roots of
$P_N(\cos\theta)=0$.  Show that for large $N$, the nodes cluster near $\theta = \pi/2$
(equatorial), the same point where the CBC generating function evaluation lives.

Speculate (one paragraph): does this imply that **Gauss-Legendre integration of $Ri_g$
is naturally aligned with the CBC structure of $\phi_h$**?  What would perfect alignment
require?

---

## Problem 6 — Rodrigues' Formula, Differential Equations, and Stability (PhD, 20 pts)

### Part A: The Legendre ODE for stability profiles (10 pts)

The Legendre polynomial $P_n(x)$ satisfies the self-adjoint ODE:
$$\frac{d}{dx}\!\left[(1-x^2)\frac{d P_n}{dx}\right] + n(n+1)P_n = 0.$$

1. Under the mapping $x = \mu = 1 - 2b_h|\zeta|$, write the ODE in terms of $\zeta$.
2. Identify the operator $\mathcal{L} = \frac{d}{d\zeta}[(b_h^{-2} - \zeta^2)\frac{d}{d\zeta}]$
	as a weighted Sturm-Liouville operator in stability space.
3. What are the **natural boundary conditions** at $\zeta=0$ (neutral) and $\zeta=-1/b_h$
	(branch point) that make the Legendre polynomials the appropriate eigenfunctions?

### Part B: Eigenfunction expansion of $Ri_g$ (10 pts)

In the degenerate case $b_m=b_h=b$, we have $Ri_g = \zeta$.  In the mapped variable,
$Ri_g = -\frac{1-\mu}{2b}$.  Expand this in Legendre polynomials.

1. Show that the only non-zero coefficients are $\hat{Ri}_0 = -\frac{1}{2b}$
	and $\hat{Ri}_1 = -\frac{1}{2b}$ (i.e., $Ri_g$ is a degree-1 polynomial in $\mu$,
	so only $P_0$ and $P_1$ contribute).
2. For the non-degenerate case $b_m \neq b_h$, how many Legendre modes does $Ri_g$ require?
	Is the expansion finite or infinite?
3. **The truncation interpretation:** A model that uses only $N$ Legendre modes effectively
	approximates $Ri_g(\zeta)$ with a degree-$N$ polynomial in $\mu$.  What is the leading
	error term in the Ri-inversion when the profile is truncated at $N=1$?

---

## Problem 7 — Research-Level: Spherical Harmonics and 3D ABL Structure (PhD Research, open-ended)

### Part A: The stability parameter as a polar angle

The Legendre generating function lives on the unit sphere in 3D.  The connection
$t = \sqrt{b_h\eta}$, $x = \cos\theta = 0$ suggests a specific geometric picture:

1. Draw a unit sphere.  Label $\theta = 0$ (north pole) as neutral ($\zeta=0$),  
	$\theta = \pi/2$ (equator) as the onset of CBC convergence failure ($|\zeta|=1/b$),  
	$\theta = \pi$ (south pole) as the strongly convective UBL limit.
2. In this picture, horizontal homogeneity (no $\phi$-dependence) means only **zonal**
	($m=0$) spherical harmonics $Y_n^0 \propto P_n(\cos\theta)$ contribute.
3. A **horizontally inhomogeneous** ABL (mesoscale variability, sea-breeze, urban heat island)
	breaks the zonal symmetry: all $(n,m)$ modes become active.  Speculate on what the
	first non-zonal correction ($m=1$ or $m=2$) would represent physically.

### Part B: Gegenbauer convolution as turbulence closure

The identity $\phi_h = \phi_m^2$ (when $b_m=b_h$) translates to:
the Gegenbauer spectral coefficients of $\phi_h$ are the *Cauchy convolution* of those
of $\phi_m$ with itself.  In turbulence language, the heat flux is the "square" of the
momentum flux in spectral stability space.

1. If we write $\phi_m = 1 + \sum_{n\ge1} a_n^{(m)} P_n(\cos\theta)$, express $\phi_h$
	as a Legendre series using the linearization formula from Problem 3B.
2. Show that the non-degenerate case ($b_m \neq b_h$, different Prandtl number) introduces
	**mixing of Legendre modes between $\phi_h$ and $\phi_m^2$** — a spectral Prandtl number
	effect.
3. Propose a criterion for "spectral closure quality" in terms of Legendre mode mixing:
	when is the mode coupling negligible, and when does it require explicit representation?

### Part C: Gauss-Legendre nodes as optimal stability sampling heights (open)

In a vertically discretized model with $N$ levels, the optimal placement of levels for
representing $Ri_g(\zeta)$ is given by the GL nodes $\zeta_k$ (zeros of $P_N$ in the
mapped domain).

1. Derive the optimal $N=4$ level heights $z_k$ (in meters) for a surface layer extending
	from $z_0 = 0.1$ m to $z_{top} = 200$ m, using the GL node mapping.
2. Compare to the **geometric mean height** $z_g = \sqrt{z_0 z_{top}}$ used in practice.
	Which better samples the non-linear curvature of $\phi_h$?
3. If the actual ABL profile is smooth in $\log z$, suggest an alternative mapping
	$\mu = f(\log z)$ that would make GL quadrature even more accurate.  Derive the
	Jacobian correction to the quadrature weights.

---

## Summary: The Complete Mathematical Chain

$$
\underbrace{\phi_h(-\eta) = (1+b_h\eta)^{-1/2}}_{\text{MOST power law}}
= \underbrace{G\!\left(\sqrt{b_h\eta},\,0\right)}_{\text{Legendre gen. fn. at equator}}
= \underbrace{\sum_{n=0}^{\infty} P_{2n}(0)\,(b_h\eta)^n}_{\text{even Legendre values}}
= \underbrace{\sum_{n=0}^{\infty} (-1)^n \binom{2n}{n}\frac{(b_h\eta)^n}{4^n}}_{\text{CBC series}}
$$

$$
\underbrace{\phi_m(-\eta) = (1+b_m\eta)^{-1/4}}_{\text{MOST power law}}
= \underbrace{G_{1/4}\!\left(\sqrt{b_m\eta},\,0\right)}_{\text{Gegenbauer }C^{(1/4)}\text{ at equator}}
= \underbrace{\sum_{n=0}^{\infty} |C_{2n}^{(1/4)}(0)|\,(b_m\eta)^n}_{\text{Gegenbauer values}}
$$

$$
\underbrace{b_m=b_h \;\Rightarrow\; \phi_h = \phi_m^2}_{\text{MOST identity / }Ri_g=\zeta}
\;\equiv\;
\underbrace{C^{(1/2)} = C^{(1/4)} \otimes C^{(1/4)}}_{\text{Gegenbauer product formula}}
$$

---

## References

- Abramowitz & Stegun (1972) §8: Legendre Functions; §22: Orthogonal Polynomials
- Szegő (1975) *Orthogonal Polynomials*, AMS — definitive reference on Gegenbauer
- DLMF §14 (Legendre), §18 (Gegenbauer/ultraspherical): https://dlmf.nist.gov
- Businger et al. (1971); Dyer (1974) — MOST parameter sets
- `hw/CBC.md` — prerequisite CBC expansion and $Ri_g=\zeta$ derivation
- `hw/MOST_Derivations_Problems.md` — MOST and grid-curvature background
- `code/phi_h_central_binomial.py` — numerical CBC implementation
- `notes/parameters.md` §5 — CBC structural summary and CBC vs asymptotic table
