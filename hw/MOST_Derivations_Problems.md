# Homework: MOST Derivations and Stability Functions

**Course:** Atmospheric Boundary Layer Physics  
**Topic:** Monin–Obukhov Similarity Theory and Grid-Scale Bias  
**Prerequisites:** Vector calculus, atmospheric thermodynamics, turbulence fundamentals

---

## Problem 1: Deriving the Gradient Richardson Number (30 points)

**Background:** Monin–Obukhov Similarity Theory (MOST) provides the foundation for understanding turbulent fluxes in the surface layer through non-dimensional similarity functions.

### Part A (10 points)
Starting from the MOST definitions:
- Momentum gradient: $\phi_m(\zeta) = \frac{\kappa z}{u_*} \frac{\partial U}{\partial z}$
- Heat gradient: $\phi_h(\zeta) = \frac{\kappa z}{\theta_*} \frac{\partial \theta}{\partial z}$

And the standard definition of gradient Richardson number:
$$Ri_g = \frac{(g/\theta) \partial\theta/\partial z}{(\partial U/\partial z)^2}$$

**Derive** the relationship $Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}$.

**Show all algebraic steps** and clearly identify where you use the Obukhov length definition $L = \frac{u_*^2 \theta}{\kappa g \theta_*}$.

### Part B (10 points)
For **neutral conditions** ($\zeta = 0$), what is $Ri_g$? Explain physically why this makes sense.

### Part C (10 points)
Using the **linear stable profile** where $\phi_m = 1 + 5\zeta$ and $\phi_h = 1 + 5\zeta$:

1. Calculate $Ri_g$ at $\zeta = 0.1$
2. Calculate $Ri_g$ at $\zeta = 1.0$
3. Plot $Ri_g(\zeta)$ for $0 \leq \zeta \leq 2$ (use Python, MATLAB, or by hand)

---

## Problem 2: Stability Functions in Richardson Number Space (25 points)

**Background:** England and McNider (1995) reformulated stability functions to work directly in $Ri$-space, which is more convenient for numerical models.

### Part A (15 points)
Starting with the linear stable form $\phi_m = 1 + \beta \zeta$ where $\beta = 5$:

1. Use the relationship $Ri = \frac{\zeta \phi_h}{\phi_m^2}$ with $\phi_h = \phi_m$ (assume Pr = 1) to show that:
   $$\zeta = \frac{Ri}{1 - \beta Ri}$$

2. Substitute this into $f_m(Ri) = \frac{1}{\phi_m^2}$ to derive:
   $$f_m(Ri) = (1 - \beta Ri)^2$$

3. What is the **critical Richardson number** $Ri_{cr}$ where $f_m \to 0$? What does this mean physically?

### Part B (10 points)
The original 1995 paper had a sign error in the derivation. If someone incorrectly derived:
$$\zeta = \frac{Ri}{1 + \beta Ri}$$

What would the resulting (incorrect) $f_m(Ri)$ look like? Why would this be nonphysical for $Ri > 0$?

---

## Problem 3: Grid-Curvature Bias and Jensen's Inequality (30 points)

**Background:** Coarse vertical grids systematically underestimate stability because MOST functions are concave-down in the stable boundary layer.

### Part A (10 points)
The **neutral curvature invariant** is defined as $\Delta = a_h - 2a_m$, where $\phi_m = 1 + a_m \zeta$ and $\phi_h = 1 + a_h \zeta$.

For the **Businger-Dyer functions** with $a_m = 5$ and $a_h = 5$:
1. Calculate $\Delta$
2. Is the function concave-up or concave-down?
3. Will Jensen's Inequality predict that $\overline{Ri_g} > Ri_g(\bar{z})$ or $\overline{Ri_g} < Ri_g(\bar{z})$?

### Part B (15 points)
Consider a two-layer atmosphere:
- Layer 1: $z_0 = 10$ m to $z_{mid} = 50$ m, with $Ri_g(30\text{ m}) = 0.10$
- Layer 2: $z_{mid} = 50$ m to $z_1 = 100$ m, with $Ri_g(75\text{ m}) = 0.25$

If a **coarse model** computes a bulk Richardson number over the entire 10–100 m layer:

1. The geometric mean height is $z_g = \sqrt{10 \times 100} = 31.6$ m. Estimate $Ri_g(z_g)$ by interpolation.
2. If the bulk $Ri_b = 0.12$, calculate the **bias ratio** $B = Ri_g(z_g) / Ri_b$.
3. Does the model **over-mix or under-mix** turbulence? Explain.

### Part C (5 points)
Sketch a graph showing how the bias ratio $B$ changes as vertical grid spacing $\Delta z$ increases from 1 m to 100 m. Label key features (convergence, asymptote, etc.).

---

## Problem 4: Physics-Informed Neural Networks (15 points)

**Background:** PINNs embed physical laws into machine learning models via a composite loss function.

### Part A (10 points)
The PINN loss function is:
$$\mathcal{L} = \mathcal{L}_{data} + \lambda \mathcal{L}_{physics} + \mu \mathcal{L}_{BC}$$

Given:
- $N = 100$ observations with $\mathcal{L}_{data} = 0.05$ (MSE in K²)
- Physics residual at 50 collocation points: $\mathcal{L}_{physics} = 0.20$ (MSE in flux units)
- Boundary condition violations: $\mathcal{L}_{BC} = 0.01$

If $\lambda = 10$ and $\mu = 5$, calculate the **total loss** $\mathcal{L}$.

### Part B (5 points)
Why is it important that $\lambda > 1$? What would happen if $\lambda = 0$? (Answer in 2–3 sentences.)

---

## Problem 5: Real-World Application (Bonus: 10 points)

You are tasked with improving a weather forecast model that has a **2°C warm bias** at night in stable conditions.

1. Using concepts from this homework, explain **two physical mechanisms** that could cause this bias.
2. Propose **one solution** involving stability functions and **one solution** involving grid resolution.
3. How would you validate your solution using observational data?

(Answer in 1 page maximum, include a simple sketch if helpful.)

---

## Submission Guidelines

- **Format:** PDF or Jupyter Notebook (if using Python)
- **Due Date:** [Instructor fills in]
- **Collaboration:** You may discuss concepts, but all derivations and code must be your own work
- **Grading Rubric:**
  - Correct mathematical steps: 60%
  - Physical interpretation: 25%
  - Clarity and organization: 15%

---

## Helpful References

- Businger et al. (1971): Flux-Profile Relationships in the Atmospheric Surface Layer
- Grachev et al. (2012): The Critical Richardson Number and Limits of Applicability of Local Similarity Theory
- McCandless et al. (2022): Machine Learning for Improving Surface-Layer-Flux Estimates

*Hint: Many of these references are in your `foundational.md` file with DOIs for easy lookup!*
