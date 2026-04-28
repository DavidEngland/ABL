# SBL IBEx Expansions (Internal-Boundary-External)

## Scope

This note summarizes geometric-series-based IBEx expansions for rational forms used in stable-boundary-layer (SBL) Richardson mappings and closures.

Let x be real unless stated otherwise.

## 1. Internal Region (|x| < 1)

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

## 2. External Region (|x| > 1)

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

## 3. SBL Form: zeta/(1 + beta zeta)

Define y = beta zeta.

$$
\frac{\zeta}{1+\beta\zeta} = \frac{1}{\beta}\frac{y}{1+y}.
$$

### 3.1 Internal (|beta zeta| < 1)

$$
\frac{\zeta}{1+\beta\zeta}
= \zeta\sum_{n=0}^{\infty}(-1)^n(\beta\zeta)^n
= \zeta - \beta\zeta^2 + \beta^2\zeta^3 - \beta^3\zeta^4 + \cdots.
$$

### 3.2 External (|beta zeta| > 1)

$$
\frac{\zeta}{1+\beta\zeta}
= \frac{1}{\beta}\left(1 - \frac{1}{\beta\zeta} + \frac{1}{(\beta\zeta)^2} - \frac{1}{(\beta\zeta)^3} + \cdots\right).
$$

Hence,

$$
\frac{\zeta}{1+\beta\zeta} \to \frac{1}{\beta}\quad \text{as}\quad \zeta\to +\infty,
$$

which is the standard SBL asymptote for the linear-branch Richardson mapping.

## 4. Boundary Region (|x| = 1)

For real x, convergence fails at x = +1 for 1/(1-x) and at x = -1 for 1/(1+x).

Complex-unit-circle analysis (x = exp(i theta) = cos(theta) + i sin(theta)) can be handled with Abel/Cesaro methods and Bernoulli-polynomial machinery in periodic settings, but this is not required for current SBL implementation.

## 5. IBEx Policy for Code

1. Internal branch: evaluate truncated OGF with recurrence when |x| <= x_int (for example x_int = 0.5).
2. External branch: evaluate truncated 1/x-series when |x| >= x_ext (for example x_ext = 5).
3. Transition bridge: for x_int < |x| < x_ext, use exact algebraic form to avoid loss of monotonicity.
4. For SBL ratio zeta/(1 + beta zeta), exact evaluation in the transition bridge is preferred in production routines.
