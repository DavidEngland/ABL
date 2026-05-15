# §X. Sturm–Liouville Structure of Businger–Dyer Functions

The classical **Businger-Dyer similarity functions** for momentum and heat exhibit power-law behavior near neutral stability. This section shows that these empirical exponents correspond to boundary behavior of eigenfunctions of a one-parameter family of **ultraspherical (Gegenbauer) Sturm-Liouville operators**.

## 1. Ultraspherical Operator on $(-1,1)$

For $\lambda > -\tfrac{1}{2}$, define the ultraspherical Sturm-Liouville operator

$$
\mathcal{L}_{\lambda}\big[y\big] = -\frac{d}{dx}\!\left[(1-x^2)^{\lambda+\frac{1}{2}}\,y'(x)\right] + \lambda(\lambda+1)(1-x^2)^{\lambda-\frac{1}{2}}\,y(x),
$$

with associated weight

$$
w_{\lambda}(x) = (1-x^2)^{\lambda-\frac{1}{2}}.
$$

The Gegenbauer polynomials $C_n^{(\lambda)}(x)$ satisfy

$$
\mathcal{L}_{\lambda}\big[C_n^{(\lambda)}\big] = n(n+2\lambda)\,C_n^{(\lambda)},
$$

and form a complete orthogonal basis in $L^2\big((-1,1),w_{\lambda}\big)$.

## 2. Pullback to MOST Stability Variable

Monin-Obukhov Similarity Theory (MOST) uses

$$
\zeta = \frac{z}{L}, \qquad \zeta < 0 \text{ in unstable conditions}.
$$

Introduce the canonical mapping

$$
x = \sqrt{-b\zeta}, \qquad b>0,
$$

where $b$ is the empirical Businger-Dyer constant. Under this transform, the ultraspherical operator induces

$$
\mathcal{L}_{\lambda}^{\zeta}[y] = -\frac{d}{d\zeta}\!\left[p_{\lambda}(\zeta)\,y'(\zeta)\right] + q_{\lambda}(\zeta)\,y(\zeta).
$$

### Coefficients and Weights

- **Stiffness:**

$$
p_{\lambda}(\zeta) = (1+b\zeta)^{\lambda+\frac{1}{2}}\,|\zeta|^{-1/2}
$$

- **Weight:**

$$
w_{\lambda}(\zeta) = (1+b\zeta)^{\lambda-\frac{1}{2}}\,|\zeta|^{-1/2}
$$

The singularity at $\zeta=0$ maps to the ultraspherical endpoint $x=0$.

## 3. Mapping Exponents to Ultraspherical Index $\lambda$

Near neutral stability, the Businger-Dyer form is

$$
\phi(\zeta) \sim (1-b\zeta)^{-\alpha}, \qquad \zeta \to 0^-.
$$

The ultraspherical generating behavior gives $(1-b\zeta)^{-\lambda}$, so the identification is

$$
\lambda = \alpha.
$$

### Theorem 1 (Ultraspherical Representation)

Let $\phi(\zeta)$ satisfy Businger-Dyer-type scaling and define $x=\sqrt{-b\zeta}$. Then:

1. The pullback $u(x)=\phi(\zeta(x))$ lies in the domain of $\mathcal{L}_{\lambda}$ with $\lambda=\alpha$.
2. $u(x)$ admits a convergent expansion

$$
u(x)=\sum_{n=0}^{\infty} a_n\,C_n^{(\lambda)}(x).
$$

3. The empirical exponent $\alpha$ uniquely determines the operator and its eigenbasis.

## 4. Corollaries for Momentum, Heat, and Scalars

| Regime | Exponent $\alpha$ | Index $\lambda$ | Eigenfunctions $y_n(\zeta)$ | Weight $w_{\lambda}(\zeta)$ |
|---|---|---|---|---|
| **Momentum** | $\tfrac{1}{4}$ | $\tfrac{1}{4}$ | $C_n^{(1/4)}\!\left(\sqrt{-b_m\zeta}\right)$ | $(1+b_m\zeta)^{-1/4}\,\lvert\zeta\rvert^{-1/2}$ |
| **Heat** | $\tfrac{1}{2}$ | $\tfrac{1}{2}$ | $P_n\!\left(\sqrt{-b_h\zeta}\right)$ | $1$ (uniform/Lebesgue) |
| **Passive scalar** | $\alpha_s$ | $\alpha_s$ | $C_n^{(\alpha_s)}\!\left(\sqrt{-b_s\zeta}\right)$ | $(1+b_s\zeta)^{\alpha_s-1/2}\,\lvert\zeta\rvert^{-1/2}$ |

> **Note:** For heat ($\lambda=\tfrac{1}{2}$), the correct basis is Legendre: $C_n^{(1/2)} \equiv P_n$, with the associated 3D-Euclidean/uniform (Lebesgue) weight. By contrast, Chebyshev polynomials of the second kind satisfy $U_n = C_n^{(1)}$, i.e., the $\lambda=1$ ultraspherical case, which would correspond to a profile exponent of $-1$ and is not physically supported by canonical MOST heat scaling.

## 5. Physical Interpretation

The ultraspherical index $\lambda$ controls endpoint singularity strength. Larger $\lambda$ corresponds to stronger buoyancy influence and sharper eigenmodes.

- **Momentum** ($\lambda=\tfrac{1}{4}$): shear-dominated, broader spectral footprint.
- **Heat** ($\lambda=\tfrac{1}{2}$): buoyancy-dominated, more localized modes.
- **Tracer** ($\lambda=\alpha_s$): interpolates between shear and buoyancy limits.

**Summary.** Businger-Dyer similarity functions can be interpreted as boundary traces of a unified ultraspherical spectral framework rather than purely empirical fits.
