SBL IBEx Expansions (Internal–Boundary–External)

A unified expansion framework for rational similarity forms in stable-boundary-layer (SBL) Richardson mappings

Scope

This note summarizes geometric-series–based IBEx expansions for rational kernels used in SBL Richardson mappings and stability‑function closures. The IBEx framework partitions the evaluation domain into:

• Internal region (\(|x|<1\)): OGF expansions about \(x=0\)
• External region (\(|x|>1\)): reciprocal expansions in powers of \(1/x\)
• Boundary region (\(|x|=1\)): nonconvergent zone requiring exact evaluation


Throughout, \(x\in\mathbb{R}\) unless otherwise stated.

---

1. Internal Region (`\(|x|<1\)`) — OGF Expansions

Internal expansions use ordinary generating functions centered at \(x=0\). These converge rapidly and admit stable recurrence evaluation.

1.1 Canonical geometric forms

\frac{1}{1-x}=\sum_{n=0}^\infty x^n,\qquad |x|<1.


\frac{1}{1+x}=\sum_{n=0}^\infty (-1)^n x^n,\qquad |x|<1.


1.2 Derived ratio

\frac{x}{1+x}
= x\sum_{n=0}^\infty (-1)^n x^n
= x - x^2 + x^3 - x^4 + \cdots,\qquad |x|<1.


These are the internal IBEx kernels.

---

2. External Region (`\(|x|>1\)`) — Reciprocal Expansions

External expansions re-express the rational form in powers of \(1/x\), preserving asymptotic behavior and avoiding cancellation.

2.1 Reciprocal-shift forms

\frac{1}{x-1}
= \frac{1}{x}\frac{1}{1-1/x}
= \sum_{n=1}^\infty x^{-n},\qquad |x|>1.


\frac{1}{x+1}
= \frac{1}{x}\frac{1}{1+1/x}
= \sum_{n=1}^\infty (-1)^{n-1}x^{-n},\qquad |x|>1.


2.2 Derived ratio

\frac{x}{x+1}
= \frac{1}{1+1/x}
= 1 - \frac{1}{x} + \frac{1}{x^2} - \frac{1}{x^3} + \cdots,\qquad |x|>1.


These are the external IBEx kernels.

---

3. SBL Form: `\(\displaystyle \frac{\zeta}{1+\beta\zeta}\)` with `\(\beta=1/Ri_c\)`

Assume SBL conditions: \(\zeta>0\), \(\beta>0\).

Define the nondimensional stability variable:

y=\beta\zeta=\frac{\zeta}{Ri_c},\qquad \beta=\frac{1}{Ri_c}.


Then

\frac{\zeta}{1+\beta\zeta}
= \frac{1}{\beta}\frac{y}{1+y}
= \frac{Ri_c\,\zeta}{Ri_c+\zeta}.


This is the canonical linear-branch Richardson mapping.

---

3.1 Internal region (`\(|\beta\zeta|<1\)`, i.e. `\(\zeta<Ri_c\)`)

\frac{\zeta}{1+\beta\zeta}
= \zeta\sum_{n=0}^\infty (-1)^n (\beta\zeta)^n
= \zeta - \beta\zeta^2 + \beta^2\zeta^3 - \beta^3\zeta^4 + \cdots.


Using \(\beta=1/Ri_c\):

\frac{\zeta}{1+\zeta/Ri_c}
= \zeta - \frac{\zeta^2}{Ri_c}
+ \frac{\zeta^3}{Ri_c^2}
- \frac{\zeta^4}{Ri_c^3} + \cdots,
\qquad \left|\frac{\zeta}{Ri_c}\right|<1.


This is the internal IBEx expansion for the SBL mapping.

---

3.2 External region (`\(|\beta\zeta|>1\)`, i.e. `\(\zeta>Ri_c\)`)

\frac{\zeta}{1+\beta\zeta}
= \frac{1}{\beta}\left(1 - \frac{1}{\beta\zeta}
+ \frac{1}{(\beta\zeta)^2}
- \frac{1}{(\beta\zeta)^3} + \cdots\right).


Using \(\beta=1/Ri_c\):

\frac{\zeta}{1+\zeta/Ri_c}
= Ri_c\left(1 - \frac{Ri_c}{\zeta}
+ \left(\frac{Ri_c}{\zeta}\right)^2
- \left(\frac{Ri_c}{\zeta}\right)^3 + \cdots\right),
\qquad \left|\frac{Ri_c}{\zeta}\right|<1.


Thus the SBL asymptote:

\frac{\zeta}{1+\beta\zeta}\to \frac{1}{\beta}=Ri_c
\quad\text{as}\quad \zeta\to +\infty.


This is the standard z‑less limit for linear Richardson closures.

---

4. Boundary Region (`\(|x|=1\)`) — Nonconvergent Zone

At \(|x|=1\), geometric series fail:

• \(x=+1\): divergence of \(1/(1-x)\)
• \(x=-1\): divergence of \(1/(1+x)\)


IBEx interpretation:

• The expansion domain is partitioned into interior (\(|x|<1\)), exterior (\(|x|>1\)), and boundary (\(|x|=1\)).
• Each rational kernel has a single relevant pole:• \(x=1\) for \(1/(1-x)\)
• \(x=-1\) for \(1/(1+x)\)

• In SBL applications with \(x=\beta\zeta>0\), the physical trajectory lies on the positive real axis and never approaches the pole at \(-1\).


Thus exact algebraic evaluation is used in the boundary region.

---

5. IBEx Policy for Production Code

A numerically stable, monotone evaluation strategy:

1. Internal branch
Use truncated OGF series when|x|\le x_{\mathrm{int}},\qquad x_{\mathrm{int}}\approx 0.5.

2. External branch
Use truncated reciprocal series when|x|\ge x_{\mathrm{ext}},\qquad x_{\mathrm{ext}}\approx 5.

3. Transition bridge
Forx_{\mathrm{int}}<|x|<x_{\mathrm{ext}},
evaluate the exact algebraic form to preserve monotonicity and avoid cancellation.
4. SBL ratio \(\zeta/(1+\beta\zeta)\)
In production SBL routines, always use the exact form in the transition region.
Internal/external expansions are reserved for asymptotics, analytic derivations, or reduced models.


---