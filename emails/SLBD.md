I’ve been refining the mathematical framework for the Businger–Dyer functions and Sturm–Liouville structures, specifically focusing on the transition to ultraspherical boundary layer geometry. Following my analysis of the **SHEBA database** (where the Grachev profiles showed strong performance), I am now incorporating data from the **SMEAR I (Värriö)** and **SMEAR II (Hyytiälä)** stations in Finland to test for geographical and vegetative variance.

Below is a summary of the current methodological advancements.

### Diagnostic Findings: Heteroscedasticity in the SBL

Diagnostic runs in Julia on the SHEBA dataset confirm that the momentum similarity function $\phi_m$ exhibits explicit heteroscedasticity in the strongly stable regime ($Ri \ge 0.25$). To prevent the "Death of Efficiency" inherent in OLS, I am implementing **Stability-Weighted Ridge Regression**. We employ a weight pivot $\zeta_{ref} \approx 2.0$ to downweight noisy, wave-dominated data at extreme stabilities while preserving fit quality in the dynamically relevant range.

### Spectral Methodology: Hilbert Space and Convolution

In the attached deck, *“Ultraspherical Boundary Layer Geometry,”* I outline the transition from classical MOST to a spectral framework. We treat the boundary layer as a Hilbert space of Gegenbauer orthogonal polynomials.

We are interpreting the momentum-heat squaring identity, $\phi_h = \phi_m^2$, as a **Clebsch–Gordan spectral convolution**. This represents a physical "Dimensional Promotion," where intermittent momentum fluctuations on a $2.5\text{-D}$ fractal manifold ($\lambda = 1/4$) are rectified into a continuous $3\text{-D}$ Euclidean heat field ($\lambda = 1/2$, Legendre basis).

## Corollaries for Momentum, Heat, and Scalars

| Regime | Exponent $\alpha$ | Index $\lambda$ | Eigenfunctions $y_n(\zeta)$ | Weight $w_{\lambda}(\zeta)$ |
|---|---|---|---|---|
| **Momentum** | $\tfrac{1}{4}$ | $\tfrac{1}{4}$ | $C_n^{(1/4)}\!\left(\sqrt{-b_m\zeta}\right)$ | $(1+b_m\zeta)^{-1/4}\,\lvert\zeta\rvert^{-1/2}$ |
| **Heat** | $\tfrac{1}{2}$ | $\tfrac{1}{2}$ | $P_n\!\left(\sqrt{-b_h\zeta}\right)$ | $1$ (uniform/Lebesgue) |
| **Passive scalar** | $\alpha_s$ | $\alpha_s$ | $C_n^{(\alpha_s)}\!\left(\sqrt{-b_s\zeta}\right)$ | $(1+b_s\zeta)^{\alpha_s-1/2}\,\lvert\zeta\rvert^{-1/2}$ |

> **Note:** For heat ($\lambda=\tfrac{1}{2}$), the correct basis is Legendre: $C_n^{(1/2)} \equiv P_n$, with the associated 3D-Euclidean/uniform (Lebesgue) weight. By contrast, Chebyshev polynomials of the second kind satisfy $U_n = C_n^{(1)}$, i.e., the $\lambda=1$ ultraspherical case, which would correspond to a profile exponent of $-1$ and is not physically supported by canonical MOST heat scaling.

## Physical Interpretation

The ultraspherical index $\lambda$ controls endpoint singularity strength. Larger $\lambda$ corresponds to stronger buoyancy influence and sharper eigenmodes.

- **Momentum** ($\lambda=\tfrac{1}{4}$): shear-dominated, broader spectral footprint.
- **Heat** ($\lambda=\tfrac{1}{2}$): buoyancy-dominated, more localized modes.
- **Tracer** ($\lambda=\alpha_s$): interpolates between shear and buoyancy limits.

**Summary.** Businger-Dyer similarity functions can be interpreted as boundary traces of a unified ultraspherical spectral framework rather than purely empirical fits.

### Wigner 3-j Symbols and the Cross-Basis Tensor

The mathematical core is the **Transfer Tensor** $T_{m,n}^k$, which projects products of the Gegenbauer basis onto the Legendre basis. By employing Wigner 3-j symbols at zero magnetic quantum numbers, we transform the nonlinear relationship into a linear algebraic constraint in coefficient space:

$$b_k = \sum_{m,n} a_m a_n T_{m,n}^k$$

This ensures that the scalar gradient—the most numerically sensitive part of the inversion—is handled by the most stable, space-filling basis in the hierarchy.

### Progress on Inverse Transformations ($Ri_g \to \zeta$)

I am currently refining a safeguarded Newton scheme for the inverse transformation from the third order of infinite spaces:

* **Initialization:** Utilizing a closed-form quadratic seed for the stable branch to guarantee convergence within 1–2 iterations.
* **Jacobian Preconditioning:** By identifying the stability axis as a pullback of a compact hyperspherical manifold, the derivative chain "flattens" the singularity near neutrality ($Ri \to 0$), maintaining quadratic convergence.

The objective is a unified pipeline that enforces spectral consistency at the coefficient level before the solver executes.

**Resources:**

* **Repo:** [github.com/DavidEngland/ultra](https://github.com/DavidEngland/ultra)
* **SMEAR Data:** [atm.helsinki.fi/smear/](https://www.atm.helsinki.fi/smear/)

I intend to analyze more SMEAR data next week at different sites in Finland and hopefully reconcile the differences observed between the subarctic (Värriö) and southern boreal (Hyytiälä) results.  I feel that surface roughness and vegetation type may be key factors in the observed variance, and I am eager to test this hypothesis with the new datasets.

Cheers,

Dave