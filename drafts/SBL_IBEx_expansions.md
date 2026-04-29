# SBL IBEx Expansions (Internal-Boundary-External)

## Scope

This note summarizes geometric-series-based IBEx expansions for rational forms used in stable-boundary-layer (SBL) Richardson mappings and closures.

Let x be real unless stated otherwise.

## 1. Internal Region ($|x|<1$)

Use ordinary generating-function (OGF) expansions about x = 0.

### 1.1 Canonical geometric forms

$$
\frac{1}{1-x} = \sum_{n=0}^{\infty} x^n = 1 + x + x^2 + x^3 + \cdots, \qquad |x|<1.
$$

$$
\frac{1}{1+x} = \sum_{n=0}^{\infty} (-1)^n x^n = 1 - x + x^2 - x^3 + \cdots, \qquad |x|<1.
$$

### 1.2 Derived ratio

$$
\frac{x}{1+x} = x\sum_{n=0}^{\infty}(-1)^n x^n = x - x^2 + x^3 - x^4 + \cdots, \qquad |x|<1.
$$

## 2. External Region ($|x|>1$)

Re-expand in powers of 1/x.

### 2.1 Reciprocal-shift forms

$$
\frac{1}{x-1} = \frac{1}{x}\frac{1}{1-1/x}
= \sum_{n=1}^{\infty} x^{-n}
= \frac{1}{x} + \frac{1}{x^2} + \frac{1}{x^3} + \cdots, \qquad |x|>1.
$$

$$
\frac{1}{x+1} = \frac{1}{x}\frac{1}{1+1/x}
= \sum_{n=1}^{\infty}(-1)^{n-1}x^{-n}
= \frac{1}{x} - \frac{1}{x^2} + \frac{1}{x^3} - \cdots, \qquad |x|>1.
$$

### 2.2 Derived ratio

$$
\frac{x}{x+1} = \frac{1}{1+1/x}
= 1 - \frac{1}{x} + \frac{1}{x^2} - \frac{1}{x^3} + \cdots, \qquad |x|>1.
$$

## 3. SBL Form: $\zeta/(1+\beta\zeta)$ with $\beta=1/Ri_c$

Assume SBL conditions throughout this section: $\zeta>0$ and $\beta>0$.

Define

$$
y = \beta\zeta = \frac{\zeta}{Ri_c}, \qquad \beta = \frac{1}{Ri_c}.
$$

$$
\frac{\zeta}{1+\beta\zeta} = \frac{1}{\beta}\frac{y}{1+y}.
$$

Equivalently,

$$
\frac{\zeta}{1+\beta\zeta}
= \frac{Ri_c\,\zeta}{Ri_c+\zeta}.
$$

### 3.1 Internal ($|\beta\zeta|<1$, i.e. $\zeta<Ri_c$)

$$
\frac{\zeta}{1+\beta\zeta}
= \zeta\sum_{n=0}^{\infty}(-1)^n(\beta\zeta)^n
= \zeta - \beta\zeta^2 + \beta^2\zeta^3 - \beta^3\zeta^4 + \cdots.
$$

Using $\beta=1/Ri_c$:

$$
\frac{\zeta}{1+\zeta/Ri_c}
= \zeta - \frac{\zeta^2}{Ri_c} + \frac{\zeta^3}{Ri_c^2} - \frac{\zeta^4}{Ri_c^3} + \cdots,
\qquad \left|\frac{\zeta}{Ri_c}\right|<1.
$$

### 3.2 External ($|\beta\zeta|>1$, i.e. $\zeta>Ri_c$)

$$
\frac{\zeta}{1+\beta\zeta}
= \frac{1}{\beta}\left(1 - \frac{1}{\beta\zeta} + \frac{1}{(\beta\zeta)^2} - \frac{1}{(\beta\zeta)^3} + \cdots\right).
$$

Using $\beta=1/Ri_c$:

$$
\frac{\zeta}{1+\zeta/Ri_c}
= Ri_c\left(1 - \frac{Ri_c}{\zeta} + \left(\frac{Ri_c}{\zeta}\right)^2 - \left(\frac{Ri_c}{\zeta}\right)^3 + \cdots\right),
\qquad \left|\frac{Ri_c}{\zeta}\right|<1.
$$

Hence,

$$
\frac{\zeta}{1+\beta\zeta} \to \frac{1}{\beta}\quad \text{as}\quad \zeta\to +\infty,
$$

equivalently,

$$
\frac{\zeta}{1+\zeta/Ri_c} \to Ri_c \quad \text{as}\quad \zeta\to +\infty.
$$

which is the standard SBL asymptote for the linear-branch Richardson mapping.

## 4. Boundary Region ($|x|=1$)

For real x, convergence fails at x = +1 for 1/(1-x) and at x = -1 for 1/(1+x).

Complex-unit-circle analysis (x = exp(i theta) = cos(theta) + i sin(theta)) can be handled with Abel/Cesaro methods and Bernoulli-polynomial machinery in periodic settings, but this is not required for current SBL implementation.

IBEx viewpoint: in complex (or higher-dimensional) language, the expansion domain is partitioned into inside the unit hypersphere ($|x|<1$), outside ($|x|>1$), and boundary ($|x|=1$). For each rational kernel there is one relevant pole (singularity): at $x=1$ for $1/(1-x)$, and at $x=-1$ for $1/(1+x)$. In SBL with $x=\beta\zeta>0$, the physical path stays on the positive real axis and does not cross the pole at $x=-1$.

## 5. IBEx Policy for Code

1. Internal branch: evaluate truncated OGF with recurrence when $|x| \le x_{\mathrm{int}}$ (for example $x_{\mathrm{int}}=0.5$).
2. External branch: evaluate truncated 1/x-series when $|x| \ge x_{\mathrm{ext}}$ (for example $x_{\mathrm{ext}}=5$).
3. Transition bridge: for $x_{\mathrm{int}}<|x|<x_{\mathrm{ext}}$, use exact algebraic form to avoid loss of monotonicity.
4. For SBL ratio $\zeta/(1+\beta\zeta)=Ri_c\zeta/(Ri_c+\zeta)$, exact evaluation in the transition bridge is preferred in production routines.

## 6. Ri_g-Centered View: Exact Quadratic Stability Function

For the canonical SBL linear shear/scalar form

$$
\phi_m=\phi_h=1+\beta\zeta, \qquad \beta>0,\ \zeta>0,
$$

the gradient Richardson number is

$$
Ri_g(\zeta)=\zeta\frac{\phi_h}{\phi_m^2}=\frac{\zeta}{1+\beta\zeta}
=\frac{Ri_c\,\zeta}{Ri_c+\zeta}, \qquad \beta=\frac{1}{Ri_c}.
$$

Inverting,

$$
\zeta=\frac{Ri_g}{1-\beta Ri_g}=\frac{Ri_c\,Ri_g}{Ri_c-Ri_g},
$$

and therefore

$$
f(R)\equiv \frac{1}{\phi_m^2}
=\frac{1}{(1+\beta\zeta)^2}
=(1-\beta R)^2
=\left(\frac{Ri_c-R}{Ri_c}\right)^2,
$$

where $R\equiv Ri_g$.

Key point: under linear SBL $\phi$, this $f(R)$ is already an exact quadratic polynomial.  No series truncation is needed once variables are written in $R$-space.

## 7. Two-Term vs Three-Term Truncation Implications

### 7.1 Truncating Ri_g as a function of zeta (internal branch)

Let $s=\beta\zeta=\zeta/Ri_c$ with $0<s<1$.

Exact:

$$
Ri_g=\zeta(1-s+s^2-s^3+\cdots).
$$

Two-term truncation:

$$
Ri_g^{(2)}=\zeta-\beta\zeta^2=\zeta(1-s).
$$

Three-term truncation:

$$
Ri_g^{(3)}=\zeta-\beta\zeta^2+\beta^2\zeta^3=\zeta(1-s+s^2).
$$

Because this is an alternating geometric remainder,

$$
|Ri_g-Ri_g^{(N)}| \le \zeta s^{N+1},
$$

and the first omitted term gives the dominant error for monotone SBL $s<1$.

Implication: moving from 2 terms to 3 terms reduces near-critical bias by one factor of $s=\zeta/Ri_c$; this is valuable only in the internal region, and rapidly loses value as $s\to 1$.

### 7.2 External branch truncation

For $q=Ri_c/\zeta<1$,

$$
Ri_g=Ri_c(1-q+q^2-q^3+\cdots).
$$

Two-term external approximation:

$$
Ri_g^{(2,ext)}=Ri_c\left(1-\frac{Ri_c}{\zeta}\right).
$$

Three-term external approximation:

$$
Ri_g^{(3,ext)}=Ri_c\left(1-\frac{Ri_c}{\zeta}+\frac{Ri_c^2}{\zeta^2}\right).
$$

Again, truncation error scales like the first omitted power in $q$.

Implication: external truncations are useful for asymptotic diagnostics and Newton seeds, not for production evaluation near the transition zone.

## 8. SBL Improvement Path Inspired by Ultraspherical Practice

The ultraspherical unstable work suggests the same robust pattern for SBL:

1. **Use exact algebraic forms in operational evaluation.**
	For linear SBL, evaluate $Ri_g=\zeta/(1+\beta\zeta)$ and $f(R)=((Ri_c-R)/Ri_c)^2$ exactly.

2. **Use truncated series only as branch-specific accelerators.**
	Internal 2-3 term forms are good for small $s$ preconditioners; external forms are good for large-$\zeta$ seeds.

3. **Keep exact bridge in transition region.**
	For $x_{int}<\beta\zeta<x_{ext}$, do not truncate; use exact rational form to preserve monotonicity and boundedness.

4. **Preserve critical-Ri structure explicitly.**
	Parameterize with $Ri_c$ first, then set $\beta=1/Ri_c$; this keeps asymptotic limits transparent and enforces the correct SBL ceiling.

5. **For non-linear SBL closures, upgrade around the exact quadratic core.**
	If higher-fidelity data require deviation from linear $\phi$, use a controlled correction on top of
	$f_0(R)=((Ri_c-R)/Ri_c)^2$ (for example, a bounded Pad\'e-style multiplicative factor), while retaining
	$f(0)=1$, $f(Ri_c)=0$, and monotone decay.

In short: the biggest SBL gain is usually not from adding one more truncated term, but from combining branch-aware asymptotic seeding with exact rational evaluation and safeguarded inversion.
