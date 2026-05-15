# §X. Sturm–Liouville Structure of Businger–Dyer Functions
The classical **Businger–Dyer similarity functions** for momentum and heat exhibit power‑law behavior near neutral stability. This section demonstrates that these empirical exponents correspond precisely to the boundary behavior of eigenfunctions of a one‑parameter family of **ultraspherical (Gegenbauer) Sturm–Liouville operators**.
## 1. Ultraspherical Operator on (-1,1)
For \lambda > -\frac{1}{2}, define the ultraspherical Sturm–Liouville operator:
with the associated weight function:
The classical **Gegenbauer polynomials** C_n^{(\lambda)}(x) satisfy the eigenvalue problem:
forming a complete orthogonal basis in L^2((-1,1), w_\lambda).
## 2. Pullback to MOST Stability Variable
Monin-Obukhov Similarity Theory (MOST) employs the stability parameter \zeta = z/L, where \zeta < 0 in unstable conditions. We introduce the canonical mapping:
x = \sqrt{-b\zeta}, \qquad b>0,
where b > 0 is the empirical Businger–Dyer constant. Under this transformation, the ultraspherical operator induces a Sturm–Liouville operator on \zeta:
### Coefficients and Weights
 * **Stiffness:** p_\lambda(\zeta) = (1+b\zeta)^{\lambda+\frac{1}{2}} |\zeta|^{-1/2}
 * **Weight:** w_\lambda(\zeta) = (1+b\zeta)^{\lambda-\frac{1}{2}} |\zeta|^{-1/2}
The singularity at \zeta=0 corresponds to the ultraspherical endpoint x=0.
## 3. Mapping Exponents to Ultraspherical Index \lambda
The Businger–Dyer similarity functions have the near‑neutral form:

The ultraspherical generating function (1 - 2xt + t^2)^{-\lambda} reduces near x=0 to (1 - b\zeta)^{-\lambda}. Thus, the empirical exponent \alpha corresponds directly to the ultraspherical index:
### Theorem 1: Ultraspherical Representation
Let \phi(\zeta) satisfy a Businger–Dyer‑type scaling. Define x = \sqrt{-b\zeta}. Then:
 1. The pullback u(x) = \phi(\zeta(x)) lies in the domain of \mathcal{L}_\lambda with \lambda = \alpha.
 2. u(x) admits a convergent expansion: u(x) = \sum_{n=0}^\infty a_n\, C_n^{(\lambda)}(x).
 3. The empirical exponent \alpha uniquely determines the operator and its eigenbasis.
## 4. Corollaries for Momentum, Heat, and Scalars
| Regime | Exponent \alpha | Index \lambda | Eigenfunctions y_n(\zeta) | Weight w_\lambda(\zeta) |
|---|---|---|---|---|
| **Momentum** | 1/4 | 1/4 | C_n^{(1/4)}(\sqrt{-b_m\zeta}) | $(1+b_m\zeta)^{-1/4} |
| **Heat** | 1/2 | 1/2 | P_n(\sqrt{-b_h\zeta}) | 1 (uniform/Lebesgue) |
| **Passive Scalar** | \alpha_s | \alpha_s | C_n^{(\alpha_s)}(\sqrt{-b_s\zeta}) | $(1+b_s\zeta)^{\alpha_s-1/2} |
> **Note:** For Heat (\lambda=1/2), the correct basis is Legendre: C_n^{(1/2)} \equiv P_n, with uniform (Lebesgue) weight (3D-Euclidean transport regime). Chebyshev polynomials of the second kind satisfy U_n = C_n^{(1)}, i.e., the \lambda=1 case, which would imply a profile exponent of -1 and is not physically supported for canonical MOST heat scaling.
>
## 5. Physical Interpretation
The ultraspherical index \lambda controls the strength of the endpoint singularity. Larger \lambda values correspond to stronger buoyancy effects and sharper eigenmodes:
 * **Momentum (\lambda=1/4):** Shear‑dominated, broad spectral footprint.
 * **Heat (\lambda=1/2):** Buoyancy‑dominated, highly localized modes.
 * **Tracer (\lambda=\alpha_s):** Interpolates between shear and buoyancy regimes.
**Summary:** The Businger–Dyer similarity functions are not merely empirical fits; they are the boundary traces of a unified spectral framework.
Does this structure clarify the "b" constant's role as a scale parameter in the mapping, or would you like to further formalize the boundary condition at the free-convection limit?