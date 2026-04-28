# Homework: Reciprocal-Sum Mixing Length Model

## Topic

Construct and analyze a generalized mixing-length model of the form

$$
\frac{1}{l_m(z)} = \sum_{n=0}^{N} \frac{1}{\kappa_n z^n}.
$$

This assignment explores dimensional consistency, asymptotic behavior, stability implications, and practical closure design.

## Audience

Advanced undergraduate STEM students through graduate students in ABL/turbulence modeling.

## Learning Objectives

By the end of this homework, you should be able to:

- Check and enforce dimensional consistency in generalized closure formulas.
- Derive asymptotic limits of a reciprocal-sum mixing-length law.
- Compare this model with the classical neutral law $l_m=\kappa z$.
- Design coefficient choices that produce physically reasonable stable-boundary-layer behavior.
- Implement and test the model in a simple $K_m=l_m^2 S f_m(Ri)$ closure.

---

## Part 1. Dimensional Analysis (Required)

Given

$$
\frac{1}{l_m(z)} = \sum_{n=0}^{N} \frac{1}{\kappa_n z^n},
$$

do the following:

1. Determine the physical units of each $\kappa_n$ needed for each term to have units of inverse length.
2. Show explicitly that the classical neutral law is recovered if only one term remains.
3. Explain why this model is best interpreted as a rational approximation in $1/l_m$.

Hint:

- $[l_m]=L$ and therefore $[1/l_m]=L^{-1}$.
- For each term, require $[1/(\kappa_n z^n)] = L^{-1}$.

---

## Part 2. Minimal Two-Term Model

Study the two-term special case:

$$
\frac{1}{l_m(z)} = \frac{1}{\kappa z} + \frac{1}{\ell_\infty},
$$

where $\kappa$ is dimensionless and $\ell_\infty$ is a length scale.

1. Invert this expression to get $l_m(z)$ explicitly.
2. Prove:

- near-surface limit ($z \to 0$): $l_m \sim \kappa z$,
- outer limit ($z \to \infty$): $l_m \to \ell_\infty$.

1. Interpret physically why this gives near-wall linear growth with outer-layer saturation.

---

## Part 3. Three-Term Model and Curvature Control

Consider

$$
\frac{1}{l_m(z)} = \frac{1}{\kappa z} + \frac{1}{\ell_\infty} + \frac{z}{\ell_2^2}.
$$

1. Derive $l_m(z)$ explicitly.
2. Determine conditions on $\ell_\infty,\ell_2$ for $l_m(z)>0$ and bounded for $z>0$.
3. Compute $dl_m/dz$ and discuss monotonicity regions.
4. Show how the third term can force re-reduction of $l_m$ aloft (if applicable) and discuss when that is physically justified.

---

## Part 4. Embedding in a Turbulence Closure

Use

$$
K_m = l_m^2 S f_m(Ri), \qquad K_h = l_m^2 S f_h(Ri),
$$

with any stable $f_m, f_h$ model from your previous homework.

1. For fixed $S$ and $Ri$, compare $K_m(z)$ from:

- classical $l_m=\kappa z$,
- two-term reciprocal model,
- three-term reciprocal model.

1. Explain how reciprocal-sum $l_m$ can reduce over-mixing on coarse grids in stable regimes.
1. Discuss how this interacts with Ri-based suppression in $f_m(Ri)$ and $f_h(Ri)$.

---

## Part 5. Numerical Experiment

Create a small script (Python/Julia/Matlab) that:

1. Defines test heights $z\in[0.5,500]$ m.
2. Evaluates at least two parameter sets for each reciprocal model.
3. Plots:

- $l_m(z)$,
- $K_m(z)$ for a fixed $S$ and representative $Ri$ values.

1. Reports whether each parameter set satisfies:

- positivity,
- smoothness,
- near-surface linear scaling,
- realistic outer behavior.

---

## Part 6. Model Design Task (Open-Ended)

Propose your own reciprocal-sum form:

$$
\frac{1}{l_m(z)} = \sum_{n=0}^{N} \frac{1}{\kappa_n z^n},
$$

or an equivalent physically scaled variant.

Your proposal must include:

1. Formula and parameter definitions.
2. Dimensional consistency proof.
3. Asymptotic limits ($z\to 0$, $z\to\infty$).
4. Physical interpretation for SBL applications.
5. One paragraph on how you would calibrate it with LES/tower data.

---

## Deliverables

Submit one report containing:

- Derivations for Parts 1-3.
- Closure discussion for Part 4.
- Figures and script summary for Part 5.
- Your custom model from Part 6.

Optional appendix:

- Code listing.

---

## Grading Rubric (100 points)

- Part 1 (dimensional analysis): 20
- Part 2 (two-term derivation + interpretation): 20
- Part 3 (three-term behavior + constraints): 20
- Part 4 (closure integration discussion): 15
- Part 5 (numerical experiment quality): 15
- Part 6 (original model design): 10

---

## Hints and Checks

- If any term in $1/l_m$ is negative and large in magnitude, you can create nonphysical $l_m<0$.
- Reciprocal summation is often numerically robust because strong limiting mechanisms combine additively in resistance form.
- In stable conditions, keeping $l_m$ bounded can be as important as choosing $f_m(Ri)$.

---

## Graduate Extension

Assume a smooth target profile $l_m^{\text{target}}(z)$ from LES and fit the reciprocal model by minimizing

$$
J = \int_{z_1}^{z_2} \left(\frac{1}{l_m^{\text{model}}(z)}-\frac{1}{l_m^{\text{target}}(z)}\right)^2 dz.
$$

1. Derive normal equations for a linear least-squares fit in reciprocal space for fixed basis functions.
2. Compare this with fitting directly in $l_m$ space and discuss which is better conditioned.
