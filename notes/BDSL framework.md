## Put Businger–Dyer into a Sturm–Liouville frame

Work in Monin–Obukhov form with \(\zeta = z/L\). For heat, the flux–gradient relation is

$$\frac{\kappa z}{\theta_*}\,\frac{d\overline{\theta}}{dz} = \phi_h(\zeta),$$


so

$$\frac{d\overline{\theta}}{d\zeta} = \frac{\theta_*}{\kappa}\,\frac{\phi_h(\zeta)}{\zeta}.$$


To see a Sturm–Liouville structure, introduce a “diffusivity” \(K_h(\zeta)\) and write a stationary balance

$$-\frac{d}{d\zeta}\Big(K_h(\zeta)\,\frac{d\overline{\theta}}{d\zeta}\Big) = \lambda\,w_*(\zeta),$$


or, more abstractly,

$$-\frac{d}{d\zeta}\Big(p(\zeta)\,\frac{d\theta}{d\zeta}\Big) + q(\zeta)\,\theta = \lambda\,w(\zeta)\,\theta,$$


with

• SL operator: $L[\theta] = -\frac{d}{d\zeta}(p(\zeta)\,\theta^\prime) + q(\zeta)\,\theta$,
• weight: $w(\zeta)$,
• eigenfunctions: $\{\varphi_n(\zeta)\}$ forming an orthogonal basis under $\langle f,g\rangle = \int w f g\,d\zeta$.


You can then:

• choose $p(\zeta)$ to encode a Businger–Dyer-like eddy diffusivity,
• or treat $\phi_h(\zeta)$ as a target function to be approximated in the SL eigenbasis.


The second route is cleaner analytically: treat $\phi_h(\zeta)$ (or its primitive) as a function on a finite interval, map to $[-1,1]$, and expand in Gegenbauer/Legendre eigenfunctions of a canonical SL problem.

---

## Central binomial expansion for the Businger–Dyer heat function

For stable conditions, a common Businger–Dyer form is

$$\phi_h(\zeta) = 1 + \beta_h\,\zeta,$$


or, in some variants,

$$\phi_h(\zeta) = (1 + \beta_h\,\zeta)^{1/2}.$$


The interesting central-binomial structure appears in the fractional power. Use the binomial series

$$(1+x)^{\alpha} = \sum_{n=0}^{\infty} \binom{\alpha}{n} x^n,
\quad
\binom{\alpha}{n} = \frac{\alpha(\alpha-1)\cdots(\alpha-n+1)}{n!}.$$


For \(\alpha = \frac{1}{2}\),

$$(1+x)^{1/2}
= \sum_{n=0}^{\infty} \binom{1/2}{n} x^n
= 1 + \frac{1}{2}x - \frac{1}{8}x^2 + \frac{1}{16}x^3 - \cdots,$$


with

$$\binom{1/2}{n}
= (-1)^{n-1}\,\frac{(2n-3)!!}{2^n\,n!}
= \frac{(-1)^{n-1}}{4^n}\,\frac{(2n)!}{(1-2n)\,(n!)^2}.$$


The central binomial coefficient

$$\binom{2n}{n} = \frac{(2n)!}{(n!)^2}$$


appears explicitly in the last form. So for

$$\phi_h(\zeta) = (1 + \beta_h \zeta)^{1/2},$$


you have

$$\phi_h(\zeta)
= \sum_{n=0}^{\infty} a_n\,\zeta^n,
\quad
a_n = \binom{1/2}{n}\,\beta_h^n
= \frac{(-1)^{n-1}}{4^n}\,\frac{\binom{2n}{n}}{1-2n}\,\beta_h^n.$$

For unstable forms like \(\phi_h(\zeta) = (1 - \gamma_h \zeta)^{-1/2}\), you get

$$(1-x)^{-1/2} = \sum_{n=0}^{\infty} \binom{-1/2}{n}(-x)^n
= \sum_{n=0}^{\infty} \frac{\binom{2n}{n}}{4^n}\,x^n,$$


so

$$\phi_h(\zeta)
= \sum_{n=0}^{\infty} \frac{\binom{2n}{n}}{4^n}\,\gamma_h^n\,\zeta^n.$$

That’s the clean central-binomial expansion you can carry into any orthogonal polynomial basis via a change of basis from monomials.

---

## Map to Gegenbauer/Legendre Sturm–Liouville problems

Gegenbauer polynomials \(C_n^{(\alpha)}(x)\) are eigenfunctions of

$$-\frac{d}{dx}\Big[(1-x^2)^{\alpha+\frac{1}{2}}\,\frac{d y}{dx}\Big]
= \lambda_n\,(1-x^2)^{\alpha-\frac{1}{2}}\,y,
\quad x\in(-1,1),$$


with weight

$$w_{\alpha}(x) = (1-x^2)^{\alpha-\frac{1}{2}}.$$


Legendre polynomials are the special case \(\alpha = \frac{1}{2}\):

$$-\frac{d}{dx}\Big[(1-x^2)\,\frac{dP_n}{dx}\Big] = n(n+1)\,P_n(x),$$
##


To use these for vertical profiles:

## Map the vertical coordinate.
Choose a finite \(\zeta\)-interval, say \([0,\zeta_{\max}]\), and map to $$[-1,1]$$ via $$x = 2\frac{\zeta}{\zeta_{\max}} - 1.$$

## Define the profile in \(x\).
For a temperature profile \(\overline{\theta}(\zeta)\), define $$\Theta(x) = \overline{\theta}(\zeta(x)).$$

## Expand in Gegenbauer/Legendre.
For Legendre: $$\Theta(x) = \sum_{n=0}^{\infty} a_n\,P_n(x),$$
with $$a_n = \frac{2n+1}{2}\int_{-1}^{1} \Theta(x)\,P_n(x)\,dx.$$
For Gegenbauer: $$\Theta(x) = \sum_{n=0}^{\infty} a_n\,C_n^{(\alpha)}(x),$$
with $$a_n = \frac{\int_{-1}^{1} \Theta(x)\,C_n^{(\alpha)}(x)\,w_{\alpha}(x)\,dx}
           {\int_{-1}^{1} [C_n^{(\alpha)}(x)]^2\,w_{\alpha}(x)\,dx}.$$



Because Gegenbauer/Legendre are SL eigenfunctions, this is exactly a Sturm–Liouville expansion of the profile.

---

## From central-binomial \(\phi_h\) to polynomial-profile expansions

You can connect the Businger–Dyer \(\phi_h\) expansion to a profile expansion in two steps:

## Integrate the series for \(\phi_h\).
From $$\frac{d\overline{\theta}}{d\zeta} = \frac{\theta_*}{\kappa}\,\frac{\phi_h(\zeta)}{\zeta},$$
and $$\phi_h(\zeta) = \sum_{n=0}^{\infty} a_n\,\zeta^n,$$
you get $$\frac{d\overline{\theta}}{d\zeta}
= \frac{\theta_*}{\kappa}\sum_{n=0}^{\infty} a_n\,\zeta^{n-1},$$
so (for \(n\ge 1\)) $$\overline{\theta}(\zeta)
= \overline{\theta}(0) + \frac{\theta_*}{\kappa}\Big(a_0\ln\zeta + \sum_{n=1}^{\infty} \frac{a_n}{n}\,\zeta^n\Big),$$
with the usual log term from the \(n=0\) contribution.
## Change basis from monomials to Gegenbauer/Legendre.
After mapping \(\zeta\mapsto x\), each \(\zeta^n\) becomes a polynomial in \(x\) of degree \(n\). You can write\zeta^n = \sum_{k=0}^{n} b_{n,k}\,P_k(x)
\quad\text{or}\quad
\zeta^n = \sum_{k=0}^{n} b_{n,k}^{(\alpha)}\,C_k^{(\alpha)}(x),
where the coefficients \(b_{n,k}\) are determined by orthogonality:b_{n,k} = \frac{2k+1}{2}\int_{-1}^{1} \zeta(x)^n\,P_k(x)\,dx,
and similarly for Gegenbauer with weight \(w_{\alpha}\).Then\overline{\theta}(x)
= \overline{\theta}(0) + \frac{\theta_*}{\kappa}\Big(a_0\ln\zeta(x)
  + \sum_{n=1}^{\infty} \frac{a_n}{n}\sum_{k=0}^{n} b_{n,k}\,P_k(x)\Big),
which you can reorganize as\overline{\theta}(x)
= \overline{\theta}(0) + \frac{\theta_*}{\kappa}\Big(a_0\ln\zeta(x)
  + \sum_{k=0}^{\infty} A_k\,P_k(x)\Big),
withA_k = \sum_{n=k}^{\infty} \frac{a_n}{n}\,b_{n,k}.



So the central-binomial structure of \(a_n\) feeds directly into closed-form (or at least systematically computable) expressions for the Legendre/Gegenbauer coefficients \(A_k\).

---

## A simpler Legendre-first route (if you want to prototype)

If you want a quick, clean prototype before going full Gegenbauer:

## Fix a finite \(\zeta\)-range (e.g. \([z_0/L, z_1/L]\)) and map to \([-1,1]\).
## Take a specific Businger–Dyer \(\phi_h\) (say \((1-\gamma_h\zeta)^{-1/2}\)), expand in central-binomial series.
## Integrate to get \(\overline{\theta}(\zeta)\) as a polynomial+log series.
## Numerically or symbolically project \(\overline{\theta}(x)\) onto Legendre polynomials using the standard inner product.


Once that pipeline feels right, you can:

• swap Legendre for Gegenbauer to tune the weight toward your preferred vertical emphasis,
• or invert the logic and design \(p(\zeta)\) in the SL operator so that its eigenfunctions approximate Businger–Dyer shapes directly.


---

If you tell me which exact Businger–Dyer form you want (stable vs unstable, exponent, coefficients), I can write out the first few explicit Legendre or Gegenbauer coefficients in closed form so you’ve got something concrete to drop into a manuscript.