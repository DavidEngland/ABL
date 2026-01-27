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

## Problem 6: Eddy Diffusivities and Stability Corrections (35 points)

**Background:** Turbulent mixing in the atmospheric boundary layer is parameterized using eddy diffusivities $K_m$ (momentum) and $K_h$ (heat). These depend critically on stability functions, which can be expressed either in MOST ($\phi_{m,h}$) or Richardson number ($f_{m,h}$) frameworks.

### Part A: Basic Diffusivity Formulation (10 points)

The **eddy diffusivity for momentum** in MOST is:
$$K_m = \frac{\kappa u_* z}{\phi_m(\zeta)}$$

And for **heat**:
$$K_h = \frac{\kappa u_* z}{\phi_h(\zeta)}$$

1. Starting from the definition $\tau = -K_m \frac{\partial U}{\partial z}$ and the MOST momentum gradient $\phi_m = \frac{\kappa z}{u_*} \frac{\partial U}{\partial z}$, **derive** the expression for $K_m$ above. (5 points)

2. If $u_* = 0.3$ m/s, $z = 50$ m, and $\phi_m = 2.0$ (stable conditions), calculate $K_m$. Use $\kappa = 0.4$. (2 points)

3. Explain physically why $K_m$ **decreases** as $\phi_m$ increases in stable conditions. (3 points)

### Part B: Richardson Number Formulation (12 points)

In the **Richardson number framework**, eddy diffusivities are expressed as:
$$K_m = l^2 S f_m(Ri_g)$$
$$K_h = l^2 S f_h(Ri_g)$$

where $l$ is the mixing length, $S = |\partial U/\partial z|$ is the wind shear magnitude, and $f_{m,h}$ are stability functions.

1. Show that for **neutral conditions** ($Ri_g = 0$), if $f_m(0) = 1$, then $K_m = l^2 S$. (3 points)

2. Using the derived form from Problem 2: $f_m(Ri) = (1 - 5Ri)^2$, calculate $K_m$ at:
   - $Ri_g = 0.05$ (weakly stable)
   - $Ri_g = 0.15$ (moderately stable)
   
   Assume $l = 40$ m and $S = 0.02$ s⁻¹. (4 points)

3. At what value of $Ri_g$ does $K_m \to 0$? What atmospheric regime does this represent? (5 points)

### Part C: Grid-Scale Corrections (13 points)

Coarse numerical models apply a **resolution-dependent correction factor** $G$ to diffusivities:
$$K_m^{corrected} = G(\zeta, \Delta z) \cdot K_m^{uncorrected}$$

where the correction factor has the form:
$$G(\zeta, \Delta z) = \exp\left[-D \left(\frac{\Delta z}{\Delta z_{ref}}\right)^p \left(\frac{\zeta}{\zeta_{ref}}\right)^q \right]$$

Given: $D = 0.5$, $p = 0.8$, $q = 1.2$, $\Delta z_{ref} = 20$ m, $\zeta_{ref} = 1.0$

1. Calculate $G$ for:
   - Fine grid: $\Delta z = 5$ m, $\zeta = 0.5$
   - Coarse grid: $\Delta z = 80$ m, $\zeta = 0.5$
   
   (4 points)

2. If the uncorrected $K_m = 10$ m²/s, what are the corrected values for each grid? (2 points)

3. Explain why $G \to 1$ as $\Delta z \to 0$. What physical principle does this preserve? (4 points)

4. **True or False:** "The correction factor $G$ increases mixing on coarse grids to compensate for unresolved turbulence." Justify your answer. (3 points)

---

## Problem 7: Comparing MOST and Richardson Formulations (Bonus: 15 points)

### Part A (10 points)

Create a **comparison table** showing how the same physical quantity is expressed in both frameworks:

| Quantity | MOST Framework $\phi(\zeta)$ | Richardson Framework $f(Ri)$ |
|----------|------------------------------|------------------------------|
| Momentum gradient | | |
| Heat gradient | | |
| Eddy diffusivity (momentum) | | |
| Critical stability limit | | |

Fill in all cells with correct mathematical expressions.

### Part B (5 points)

A forecaster provides you with observational data: $z = 10$ m, $u_* = 0.25$ m/s, $\theta_* = -0.05$ K, $\theta = 285$ K, and vertical gradients $\partial U/\partial z = 0.03$ s⁻¹ and $\partial \theta/\partial z = 0.02$ K/m.

1. Calculate $\zeta$ using the Obukhov length formula.
2. Calculate $Ri_g$ using its definition.
3. Verify that $Ri_g(\zeta) = \zeta \phi_h/\phi_m^2$ using Businger-Dyer functions ($\phi_m = \phi_h = 1 + 5\zeta$).

---

## Problem 8: General Transformation Procedure from $\phi(\zeta)$ to $f(Ri)$ (40 points)

**Background:** Numerical weather prediction models often prefer Richardson number formulations over traditional MOST because $Ri$ can be computed directly from model variables without knowing surface fluxes. This problem develops the systematic procedure to transform any given $\phi_{m,h}(\zeta)$ functions into their equivalent $f_{m,h}(Ri_g)$ forms.

### Part A: The Fundamental Relationships (8 points)

Given arbitrary stability functions $\phi_m(\zeta)$ and $\phi_h(\zeta)$, state the three fundamental relationships that connect the MOST and Richardson frameworks:

1. The **gradient Richardson number** in terms of $\zeta$ and $\phi$ functions: (3 points)
   $$Ri_g(\zeta) = ?$$

2. The **momentum stability function** transformation: (2 points)
   $$f_m(Ri) = ?$$

3. The **heat stability function** transformation: (3 points)
   $$f_h(Ri) = ?$$

Write each in the most general form (do not assume specific $\phi$ functions yet).

### Part B: Step-by-Step Algorithm (20 points)

Develop a **general algorithm** for transforming $\phi_m(\zeta)$ and $\phi_h(\zeta)$ into $f_m(Ri_g)$ and $f_h(Ri_g)$. Your algorithm should work for any monotonic stability functions.

**Fill in the following procedure:**

---

**ALGORITHM: MOST → Richardson Transformation**

**Given:** Stability functions $\phi_m(\zeta)$ and $\phi_h(\zeta)$ for $\zeta \geq 0$ (stable regime)

**Goal:** Derive $f_m(Ri_g)$ and $f_h(Ri_g)$

---

**Step 1:** Calculate the **intermediate function** $Ri_g(\zeta)$

Using the relationship from Part A, express:
$$Ri_g(\zeta) = \text{[student fills in]}$$

This gives you $Ri_g$ as an **explicit function of** $\zeta$. (3 points)

---

**Step 2:** **Invert** the relationship to obtain $\zeta(Ri_g)$

Solve the equation from Step 1 algebraically or numerically to get:
$$\zeta = \zeta(Ri_g)$$

*Note:* This is usually the hardest step. For simple $\phi$ (e.g., linear), you can solve analytically. For complex $\phi$, numerical inversion may be required. (5 points)

**Practical tip:** What condition must $Ri_g(\zeta)$ satisfy for the inverse to exist? (2 points)

---

**Step 3:** Derive $f_m(Ri_g)$

Substitute the inverted function $\zeta(Ri_g)$ from Step 2 into:
$$f_m(Ri_g) = \frac{1}{\phi_m(\zeta(Ri_g))^2}$$

Simplify to obtain $f_m$ as an **explicit function of $Ri_g$ only**. (3 points)

---

**Step 4:** Derive $f_h(Ri_g)$

Using the relationship between $f_m$, $f_h$, and the turbulent Prandtl number:
$$\frac{f_h}{f_m} = Pr_t(Ri_g) = \text{[student fills in]}$$

Therefore:
$$f_h(Ri_g) = f_m(Ri_g) \cdot Pr_t(Ri_g)$$

Express $Pr_t$ in terms of $\phi_m(\zeta)$ and $\phi_h(\zeta)$, then substitute $\zeta(Ri_g)$. (4 points)

---

**Step 5:** Verify limiting behavior

Check that your derived functions satisfy:

- **Neutral limit:** $f_m(0) = ?$ and $f_h(0) = ?$ (1 point)
- **Critical Richardson number:** Does $f_m(Ri_{cr}) = 0$ for some finite $Ri_{cr}$? If so, what is it? (2 points)

---

### Part C: Worked Example with Businger-Dyer Functions (12 points)

Apply your algorithm to the **Businger-Dyer stable functions**:
$$\phi_m(\zeta) = 1 + 5\zeta, \quad \phi_h(\zeta) = 1 + 5\zeta$$

**Execute each step:**

1. **Step 1:** Calculate $Ri_g(\zeta)$ explicitly. (2 points)

2. **Step 2:** Invert to find $\zeta(Ri_g)$. Show your algebra. (3 points)

3. **Step 3:** Derive $f_m(Ri_g) = (1 - 5Ri_g)^2$. (2 points)

4. **Step 4:** What is $Pr_t(Ri_g)$ for this case? Therefore, what is $f_h(Ri_g)$? (3 points)

5. **Step 5:** Verify that $f_m(0) = 1$ and find $Ri_{cr}$. (2 points)

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
- Dyer (1974): A Review of Flux-Profile Relationships
- Grachev et al. (2012): The Critical Richardson Number and Limits of Applicability of Local Similarity Theory
- McCandless et al. (2022): Machine Learning for Improving Surface-Layer-Flux Estimates
- Holtslag et al. (2013): Stable Atmospheric Boundary Layers and Simulating Them

*Hint: Many of these references are in your `foundational.md` file with DOIs for easy lookup!*
