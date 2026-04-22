# Homework: From Central Binomial Expansions to SBL Stability and Shear Functions

## Audience and Level

- Bright undergraduates in STEM through graduate students in atmospheric science, applied math, physics, or engineering.
- Recommended prerequisites: calculus, series expansions, basic fluid mechanics, and introductory boundary-layer meteorology.

## Learning Goals

By the end of this assignment, you should be able to:

- Derive and use central binomial series for square-root singular forms.
- Explain why $(1-bx)^{-1/2}$ has a real-domain cutoff at $x=1/b$.
- Connect MOST variables $\zeta=z/L$, $\phi_m$, $\phi_h$ to Richardson-number closures.
- Build stable-branch $f_m(Ri)$ and $f_h(Ri)$ functions (linear and quadratic).
- Design a hybrid closure: MOST in unstable/near-neutral, Ri-based in stable conditions.
- Insert your own custom stable function in a clean software interface.

## Part 0: Setup and Notation

Use the following notation throughout:

- $\zeta = z/L$
- $Ri_g = \dfrac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2}$
- MOST relation: $Ri_g(\zeta)=\zeta\,\phi_h(\zeta)/\phi_m(\zeta)^2$
- Mixing-length closure form: $K_m=l^2 S f_m$, $K_h=l^2 S f_h$, where $S=|\partial U/\partial z|$

---

## Part 1: Central Binomial Expansion Warm-Up

### 1A. Derive the core generating function

Starting from the generalized binomial theorem,

$$
(1-z)^{-1/2} = \sum_{n=0}^{\infty}\binom{2n}{n}\frac{z^n}{4^n}, \quad |z|<1,
$$

show that

$$
\frac{1}{\sqrt{1-4x}} = \sum_{n=0}^{\infty} \binom{2n}{n}x^n, \quad |x|<\frac14.
$$

### 1B. Extend to $1/\sqrt{1-16x}$

Use $z=16x$ in the same theorem to prove

$$
\frac{1}{\sqrt{1-16x}} = \sum_{n=0}^{\infty} \binom{2n}{n}(4x)^n,
\quad |x|<\frac{1}{16}.
$$

Write the first 6 nonzero terms explicitly.

### 1C. Singularity and branch discussion

For

$$
f(x)=(1-bx)^{-1/2}, \quad b>0,
$$

answer:

1. What is the real-valued domain?
2. What happens at $x=1/b$?
3. Why is $x>1/b$ not real-valued (without analytic continuation)?
4. Give the principal complex continuation for $x>1/b$.

---

## Part 2: MOST and SBL Basics (Businger-Dyer style)

Assume the common stable linear forms:

$$
\phi_m(\zeta)=1+a_m\zeta, \qquad \phi_h(\zeta)=1+a_h\zeta, \qquad \zeta>0.
$$

Use typical values as one example: $a_m\approx 4.7$, $a_h\approx 7.8$.

### 2A. Derive $Ri_g(\zeta)$

Starting from

$$
Ri_g(\zeta)=\zeta\frac{\phi_h}{\phi_m^2},
$$

derive the explicit rational expression in $\zeta$.

### 2B. Near-neutral expansion and curvature invariant

Expand $Ri_g(\zeta)$ near $\zeta=0$ up to $O(\zeta^2)$ and identify

$$
\Delta = a_h-2a_m,
$$

then verify

$$
\left.\frac{d^2Ri_g}{d\zeta^2}\right|_{\zeta=0}=2\Delta.
$$

Interpret the sign of $\Delta$ for concavity and bulk-vs-point bias.

### 2C. Short interpretation

In 3-5 sentences, explain why concave-down $Ri_g(z)$ causes coarse-layer $Ri_b$ to underestimate point stability.

---

## Part 3: From MOST Functions to Ri-Based Stability Functions

In Ri-closure form, define

$$
f_m(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))^2},
\qquad
f_h(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\phi_h(\zeta(Ri_g))}.
$$

### 3A. Conceptual derivation

Explain why these formulas are consistent with the MOST diffusivities

$$
K_m=\frac{\kappa z u_*}{\phi_m},\qquad K_h=\frac{\kappa z u_*}{\phi_h},
$$

and with mixing-length form $K=l^2Sf$.

### 3B. Approximate inversion task

Use a near-neutral inversion (or symbolic expansion) to obtain an approximate $\zeta(Ri_g)$ through at least second order, then substitute into $f_m$ and $f_h$ to get low-order Ri-series forms.

---

## Part 4: Build Practical Stable Functions (Linear and Quadratic)

You want stable-branch functions that depend mainly on $Ri$:

### 4A. Linear model

$$
f_m^{(L)}(Ri)=\frac{1}{1+c_m Ri},
\qquad
f_h^{(L)}(Ri)=\frac{1}{1+c_h Ri}.
$$

Choose one coefficient set and justify it physically.

### 4B. Quadratic model

$$
f_m^{(Q)}(Ri)=\frac{1}{1+c_{m1}Ri+c_{m2}Ri^2},
\qquad
f_h^{(Q)}(Ri)=\frac{1}{1+c_{h1}Ri+c_{h2}Ri^2}.
$$

Pick coefficients that create a stronger high-stability tail damping than Part 4A.

### 4C. Comparison figure

Plot $f_m$ and $f_h$ for both models on $Ri\in[0,1.5]$. Briefly discuss:

- near-neutral slope behavior
- strong-stability behavior
- implications for $K_m$ and $K_h$

---

## Part 5: Hybrid Closure Design (Unstable via MOST, Stable via Ri)

Design a branch logic:

- If $Ri\le 0$: use MOST via $\zeta=z/L$.
- If $Ri>0$: use Ri-based stable functions (linear or quadratic).

### 5A. Pseudocode

Write pseudocode for

- input: $Ri,\zeta,\text{selector}$
- output: $f_m,f_h$
- options: `LINEAR`, `QUADRATIC`, and `CUSTOM`

### 5B. Custom insertion exercise

Add one custom stable model and document it in one paragraph.
Example template:

```text
CASE(CUSTOM)
  f_m = ...
  f_h = ...
```

Your model must satisfy all three:

1. $0\le f_m,f_h\le 1$
2. non-increasing with $Ri$ for $Ri>0$
3. finite for all tested $Ri\in[0,2]$

---

## Part 6: Short WRF-Oriented Reflection

Answer in 1 page max:

1. Why many operational schemes retain Ri-based stable functions even when rooted in MOST ideas?
2. What do you gain from hybrid logic in practice?
3. Which branch (MOST or Ri-based) is most sensitive to grid spacing and why?

---

## Deliverables

Submit one PDF or Markdown report containing:

- derivations for Parts 1-3
- chosen linear and quadratic forms with coefficients
- one plot panel for Part 4
- hybrid pseudocode and custom insertion for Part 5
- short reflection for Part 6

Optional code appendix (Python/Julia/Matlab/Fortran) is encouraged.

## Grading Rubric (100 pts)

- Part 1 (series + domain/branch): 20
- Part 2 (MOST derivation + curvature meaning): 20
- Part 3 (MOST-to-Ri mapping): 20
- Part 4 (function design + plot interpretation): 20
- Part 5 (hybrid logic + custom extension): 15
- Part 6 (WRF-oriented reflection clarity): 5

## Hints

- For central binomial terms, use

$$
\binom{2n}{n}=\frac{(2n)!}{(n!)^2}.
$$

- Radius of convergence follows nearest singularity in the complex plane.
- In stable SBL work, do not confuse:
- MOST similarity functions $\phi_m(\zeta),\phi_h(\zeta)$
- Ri-based closure functions $f_m(Ri),f_h(Ri)$

## Suggested Reading in This Repository

- Shear and stability relations:
- Shear and Stability Functions.md
- Ri-based candidates:
- Candidate Ri-Based Stability Functions.md
- Proposed_Ri_based_Stability_Functions.md
- Central binomial support:
- HW_Central_Binomial_MOST.md
- central binomial expansion.md

## Extension (Graduate Challenge)

Assume

$$
Ri_g(\zeta)=\zeta+\Delta\zeta^2+c_1\zeta^3+O(\zeta^4).
$$

1. Perform a formal inversion to get $\zeta(Ri_g)$ to $O(Ri_g^3)$.
2. Substitute into your chosen $f_m(\zeta)$ and $f_h(\zeta)$ expansions.
3. Show how coefficient choices in stable Ri functions can preserve a target neutral curvature invariant.
