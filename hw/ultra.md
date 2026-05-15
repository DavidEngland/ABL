### Part I: Physical Dimensionality and Scaling

**Problem 1: The Fractal Fingerprint of Momentum**
Using the relation $\lambda = (d - 2)/2$, where $d$ is the effective spectral dimension:
1. Calculate the required $\lambda$ for heat transport, assuming it is a space-filling 3D Euclidean process.
2. Calculate the required $\lambda$ for momentum transport, assuming it is an anisotropic shear layer living on a 2.5D fractal manifold.
3. **Physical Justification:** Explain why the factor $(1+\zeta)^{1/3}$ in the Grachev baseline is considered the "spectroscopic fingerprint" of a 2.5D energy cascade.

**Problem 2: The Reynolds Analogy Break**
In the Highly Stable Nocturnal Boundary Layer (HSNBL), momentum transport often undergoes a "dimensional collapse" where $\lambda \to 0$.
1. Identify the specific family of orthogonal polynomials associated with the limit $\lambda = 0$.
2. Explain physically why heat transport (scalar) typically maintains a dimension of $d=3.0$ even when the wind field flattens into "pancake turbulence" ($d=2.0$).

---

### Part II: Mathematical Basis and Polynomials

**Problem 3: Deriving the Basis**
Gegenbauer polynomials are defined by the three-term recurrence:
$C_0^{(\lambda)}(x) = 1$
$C_1^{(\lambda)}(x) = 2\lambda x$
$C_n^{(\lambda)}(x) = \frac{1}{n} [2x(n+\lambda-1)C_{n-1}^{(\lambda)}(x) - (n+2\lambda-2)C_{n-2}^{(\lambda)}(x)]$
*   **Task:** For the momentum case ($\lambda = 1/4$), derive the explicit form of $C_2^{(1/4)}(x)$ and $C_3^{(1/4)}(x)$. Compare your results to the monomial forms provided in the sources.

**Problem 4: Normalization at the Pole**
Using the identity $C_n^{(\alpha)}(1) = \frac{\Gamma(2\alpha + n)}{\Gamma(2\alpha) n!}$:
1. Calculate the value of the first four momentum basis functions ($\lambda = 1/4$) at the pole $x=1$.
2. **Physical Context:** In a log-tanh mapping where $\xi = \tanh(\alpha_\xi \ln(1+\zeta))$, what physical stability regime does the point $x=1$ represent?.

---

### Part III: Sturm-Liouville Theory and Operators

**Problem 5: Pullback Operators**
The Legendre operator for heat transport in stability space is defined as:
$\mathcal{L}_h = \frac{d}{d\zeta} \left[ (b_h^{-2} - \zeta^2) \frac{d}{d\zeta} \right]$
1. Show that the weight function $w(\zeta)$ for this operator is uniform ($w=1$).
2. Contrast this with the momentum operator $\mathcal{L}_m$, which has a singular weight $w_m(\zeta) = (b_m^{-2} - \zeta^2)^{-1/4}$. Explain why this singular weight is beneficial for approximating wind shear near the neutral limit.

**Problem 6: Generating Functions**
The Businger-Dyer similarity functions are identified as generating functions for ultraspherical polynomials.
*   **Task:** Expand the unstable heat function $\phi_h(\zeta) = (1 - b_h \zeta)^{-1/2}$ into a Taylor series about $\zeta=0$. Show that the coefficients are the central binomial coefficients.

---

### Part IV: Numerical Schemes and Inversion

**Problem 7: The Quadratic Seed for Inversion**
To invert the gradient Richardson number ($Ri_g$) to find the stability parameter ($\zeta$), models often use a Newton-Raphson scheme.
*   **Task:** Starting from the relationship $Ri_g = \frac{\zeta(1+\beta_h\zeta)}{(1+\beta\zeta)^2}$, derive the "Quadratic Seed" formula used for initialization:
$(\beta_h - Ri_g\beta^2)\zeta^2 + (1-2Ri_g\beta)\zeta - Ri_g = 0$.
Explain why this seed is more efficient than a trivial guess of $\zeta = Ri_g$.

**Problem 8: The Log-Tanh Mapping**
A linear mapping $\xi = \zeta / \zeta_{max}$ is often criticized for wasting resolution.
1. Sketch the Jacobian $d\xi/d\zeta$ for the compound map $\xi = \tanh(\alpha_\xi \ln(1+\zeta))$.
2. Demonstrate how this mapping provides roughly equal representation for each decade of stability (e.g., $\zeta = 0.1, 1, 10$).

---

### Part V: Statistical Inference and Heteroscedasticity

**Problem 9: The "Death of Efficiency"**
1. Define **heteroscedasticity** and explain why it is a diagnostic hallmark of the strongly stable boundary layer.
2. List three consequences of ignoring heteroscedasticity when performing standard OLS regression on ABL data.
3. **Mitigation:** Describe how **Stability-Weighted Ridge Regression** uses a weight pivot ($\zeta_{ref} \approx 2.0$) to protect the fit quality in the dynamically relevant range.

---

### Part VI: Subgrid Curvature (Advanced)

**Problem 10: The Curvature Operator $\chi$**
Given a profile projected onto a Legendre basis: $f(x) = A_0 P_0(x) + A_1 P_1(x) + A_2 P_2(x) + \dots$
1. At the midpoint $x=0$, show that the curvature is $f''(0) = 3A_2$.
2. Define the intrinsic curvature operator $\chi = \frac{\Theta_2}{\Theta_1} - 2\frac{U_2}{U_1}$.
3. **Application:** If a model cell finds a positive $\chi \Delta z$, explain how this functions as a "subgrid bifurcation detector" for identifying turbulence collapse in the upper half of a grid cell.