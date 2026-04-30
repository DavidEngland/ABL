Spherical Analysis of Monin–Obukhov Similarity Theory

A Complete Manuscript Skeleton with Equations and Operator Derivations

---

1. Introduction

Monin–Obukhov Similarity Theory (MOST) provides the canonical scaling framework for the atmospheric surface layer (ASL). MOST expresses non‑dimensional gradients of mean wind and scalars as functions of the stability parameter

\zeta = \frac{z}{L},


where \(L\) is the Obukhov length.

The central thesis of this work is:

MOST similarity functions are generating functions of ultraspherical (Gegenbauer) polynomial families, and the physical validity of MOST corresponds to the degree of isotropy in the turbulence field.

This establishes a geometric and operator‑theoretic foundation for similarity theory.

---

2. Similarity Functions as Ultraspherical Generating Functions

2.1 Classical MOST forms

The canonical Businger–Dyer similarity functions are

\phi_m(\zeta) = (1 - b_m \zeta)^{-1/4}, \qquad
\phi_h(\zeta) = (1 - b_h \zeta)^{-1/2}.


2.2 Gegenbauer generating function

The generating function for Gegenbauer polynomials \(C_n^{(\lambda)}(x)\) is

(1 - 2xt + t^2)^{-\lambda}
= \sum_{n=0}^\infty C_n^{(\lambda)}(x)\, t^n.


MOST corresponds to the one‑sided specialization

(1 - b\zeta)^{-\lambda}
= \sum_{n=0}^\infty \binom{n+\lambda-1}{n} (b\zeta)^n,


which is the \(x=1\) limit of the ultraspherical generating function.

Thus:

• Momentum: \(\lambda = 1/4\)
• Heat/scalars: \(\lambda = 1/2\)


2.3 Squaring identity

When \(b_m = b_h = b\),

\phi_h(\zeta) = (1 - b\zeta)^{-1/2}
= \left[(1 - b\zeta)^{-1/4}\right]^2
= \phi_m(\zeta)^2.


This is the ultraspherical convolution identity

C_n^{(1/2)} = \sum_{k=0}^n C_k^{(1/4)} C_{n-k}^{(1/4)}.


2.4 Turbulent Prandtl number as a spectral ratio

Pr_t(\zeta) = \frac{\phi_h(\zeta)}{\phi_m(\zeta)}
= (1 - b\zeta)^{-1/4}.


Thus \(Pr_t\) is not a constant but a spectral transfer function.

---

3. Ultraspherical Sturm–Liouville Operator

3.1 Governing operator

Gegenbauer polynomials satisfy the Sturm–Liouville equation

\mathcal{L}_\lambda[y] = -n(n+2\lambda)\, y,


where

{
\mathcal{L}_\lambda[y]
= (1-x^2)\, y^\prime^\prime - (2\lambda+1)x\, y^\prime.
}


3.2 Mapping MOST to operator form

Let

x = b\zeta.


Then the MOST similarity function

\phi(\zeta) = (1 - b\zeta)^{-\lambda}
= (1-x)^{-\lambda}


satisfies

\mathcal{L}_\lambda[\phi] = \lambda(2\lambda+1)\,(1-x)^{-\lambda-1}.


3.3 Endpoint singularities

The operator becomes singular at

x = \pm 1 \quad \Longleftrightarrow \quad \zeta = \pm \frac{1}{b}.


These correspond to:

• Very unstable limit (\(\zeta \to -1/b\))
• Very stable limit (\(\zeta \to +1/b\))


The finite radius of convergence explains why MOST fails in extreme regimes.

---

4. Turbulence Anisotropy and the Sphere of Isotropy

4.1 Reynolds stress invariants

Let \(R_{ij} = \overline{u_i^\prime u_j^\prime}\).
Define the anisotropy tensor

a_{ij} = \frac{R_{ij}}{2k} - \frac{1}{3}\delta_{ij},
\qquad k = \frac{1}{2}R_{ii}.


The invariants are

II = -\frac{1}{2} a_{ij}a_{ij}, \qquad
III = \frac{1}{3} a_{ij}a_{jk}a_{ki}.


4.2 Barycentric map

The turbulence state is mapped to a point in the triangle spanned by:

• 1‑component (1C)
• 2‑component (2C)
• 3‑component isotropic (3C)


4.3 Isotropy as a sphere

The isotropic limit corresponds to equal energy in all directions:

R_{11} = R_{22} = R_{33}.


This is the center of the Barycentric triangle, and the natural basis for functions on this manifold is the ultraspherical basis.

4.4 MOST validity criterion

MOST implicitly assumes turbulence is near the 2C–3C edge, i.e., “near‑spherical” in energy distribution.

Breakdown occurs when turbulence approaches:

• 1C (line‑like)
• Axisymmetric contraction
• Highly stable stratification


because the ultraspherical symmetry collapses.

---

5. Hyperspherical Harmonics and Multi‑Point Scaling

5.1 Hyperspherical harmonics

Solutions to Laplace’s equation in \(d\)-dimensional space:

\nabla_d^2 Y_{n}^{(\lambda)}(\Omega) = -n(n+2\lambda) Y_{n}^{(\lambda)}(\Omega),


with \(\lambda = (d-2)/2\).

5.2 MOST as a 1D marginal

The similarity functions correspond to the zonal (m=0) ultraspherical harmonics evaluated at \(x=1\).

5.3 Multi‑point closures

Second‑order closures (EFB, MMO) produce kernels of the form

K(\zeta) = \sum_{n} A_n C_n^{(\lambda)}(b\zeta),


which reduce to MOST when truncated at \(n=0\).

Thus MOST is the lowest‑order truncation of a hyperspherical harmonic expansion.

---

6. Toward an Operator‑Consistent Similarity Theory

6.1 Closure economy

Enforcing the ultraspherical structure reduces empirical parameters:

• \(b\) sets the branch point
• \(\lambda\) sets the turbulence “dimension”
• Higher‑order coefficients follow from operator constraints


6.2 Spectral filtering interpretation

Scalar gradients correspond to the self‑convolution of momentum modes:

\phi_h = \phi_m * \phi_m.


6.3 Path forward

A unified MOST must:

1. Respect the ultraspherical Sturm–Liouville operator
2. Incorporate anisotropy via Reynolds‑stress invariants
3. Extend similarity theory to multi‑point, hyperspherical structures


---

7. Appendix: Key Derivations

A. Derivative of MOST similarity function

\phi(\zeta) = (1 - b\zeta)^{-\lambda}


\phi^\prime(\zeta) = \lambda b (1 - b\zeta)^{-\lambda-1}


\phi^\prime^\prime(\zeta) = \lambda(\lambda+1)b^2 (1 - b\zeta)^{-\lambda-2}.


B. Verification of Sturm–Liouville form

Insert into

(1-x^2)\phi^\prime^\prime - (2\lambda+1)x\phi^\prime


with \(x=b\zeta\) to obtain the operator identity in §3.

---

Next Steps:

Turn this skeleton into:

• a full LaTeX manuscript,
• a Methods section with full operator proofs, or
• a publication‑ready introduction with citations and narrative polish.
