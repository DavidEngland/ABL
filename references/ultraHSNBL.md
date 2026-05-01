# Ultraspherical / HSNBL Reference Guide

This note is a working reference map for the ultraspherical-harmonic interpretation of the highly stable nocturnal boundary layer (HSNBL). It is organized to keep three things separate:

1. established surface-layer and stable-boundary-layer literature,
2. mathematical references on Gegenbauer / ultraspherical structure,
3. the specifically proposed bridge from MOST exponents to ultraspherical parameters and effective turbulent dimension.

That separation matters. The HSNBL physics is well supported by the SHEBA and stable-boundary-layer literature. The mathematics of Gegenbauer generating functions and hyperspherical harmonics is also standard. The direct identification of MOST similarity exponents with ultraspherical spectral coordinates is the new synthesis and should be cited as a proposal unless and until it is published.

## 1. Core Atmospheric References

These are the papers to cite for the physical behavior of stable similarity functions, critical Richardson number limits, and the Arctic very-stable regime.

### 1.1 Foundational MOST and canonical similarity forms

```bibtex
@article{BusingerEtAl1971,
  author = {Businger, J. A. and Wyngaard, J. C. and Izumi, Y. and Bradley, E. F.},
  title = {Flux-Profile Relationships in the Atmospheric Surface Layer},
  journal = {Journal of the Atmospheric Sciences},
  volume = {28},
  number = {2},
  pages = {181--189},
  year = {1971},
  doi = {10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2}
}

@article{Dyer1974,
  author = {Dyer, A. J.},
  title = {A Review of Flux-Profile Relationships},
  journal = {Boundary-Layer Meteorology},
  volume = {7},
  number = {3},
  pages = {363--372},
  year = {1974},
  doi = {10.1007/BF00240838}
}

@article{Hogstrom1988,
  author = {Hogstrom, U.},
  title = {Non-dimensional Wind and Temperature Profiles in the Atmospheric Surface Layer: A Re-evaluation},
  journal = {Boundary-Layer Meteorology},
  volume = {42},
  number = {1-2},
  pages = {55--78},
  year = {1988},
  doi = {10.1007/BF00119875}
}
```

Use these when you need support for the classical forms

$$
\phi_m = 1 + \beta_m \zeta, \qquad \phi_h = \Pr_t + \beta_h \zeta
$$

or for the canonical unstable power-law families

$$
\phi_m = (1 - b_m \zeta)^{-1/4}, \qquad \phi_h = (1 - b_h \zeta)^{-1/2}.
$$

These papers support the empirical similarity framework. They do not, by themselves, support the ultraspherical interpretation.

### 1.2 Stable and very-stable SHEBA references

```bibtex
@article{Grachev2005,
  author = {Grachev, A. A. and Fairall, C. W. and Persson, P. O. G. and Andreas, E. L. and Guest, P. S.},
  title = {Stable Boundary-Layer Scaling Regimes: The SHEBA Data},
  journal = {Boundary-Layer Meteorology},
  volume = {116},
  number = {2},
  pages = {201--235},
  year = {2005},
  doi = {10.1007/s10546-004-2729-0}
}

@article{GrachevEtAl2007,
  author = {Grachev, A. A. and Andreas, E. L. and Fairall, C. W. and Guest, P. S. and Persson, P. O. G.},
  title = {{SHEBA} Flux-Profile Relationships in the Stable Atmospheric Boundary Layer},
  journal = {Boundary-Layer Meteorology},
  volume = {124},
  number = {3},
  pages = {315--333},
  year = {2007},
  doi = {10.1007/s10546-007-9177-6}
}

@article{Grachev2012,
  author = {Grachev, Andrey A. and Andreas, Edgar L. and Fairall, Christopher W. and Guest, Peter S. and Persson, P. Ola G.},
  title = {The Critical Richardson Number and Limits of Applicability of Local Similarity Theory in the Stable Boundary Layer},
  journal = {Boundary-Layer Meteorology},
  volume = {147},
  number = {1},
  pages = {51--82},
  year = {2012},
  doi = {10.1007/s10546-012-9771-0}
}

@article{Grachev2014,
  author = {Grachev, Andrey A. and Andreas, Edgar L. and Fairall, Christopher W. and Guest, Peter S. and Persson, P. Ola G.},
  title = {Similarity Theory Based on the Dougherty--Ozmidov Length Scale},
  journal = {Quarterly Journal of the Royal Meteorological Society},
  volume = {141},
  number = {690},
  pages = {1845--1856},
  year = {2014},
  doi = {10.1002/qj.2488}
}

@article{GryanikEtAl2020,
  author = {Gryanik, V. M. and Lupkes, C. and Grachev, A. and Sidorenko, D.},
  title = {New Modified and Extended Stability Functions for the Stable Boundary Layer Based on {SHEBA} and Parametrizations of Bulk Transfer Coefficients for Climate Models},
  journal = {Journal of the Atmospheric Sciences},
  volume = {77},
  number = {7},
  pages = {2687--2716},
  year = {2020},
  doi = {10.1175/JAS-D-19-0255.1}
}
```

These are the key references for the HSNBL side of the story.

They support:

- the breakdown of local similarity at large $\zeta$,
- multiple stable-scaling regimes,
- asymptotic or quasi-z-less behavior in the very-stable Arctic surface layer,
- the use of SHEBA as the main observational anchor,
- the need for non-classical closure behavior as $Ri_g$ approaches its effective limit.

If you want one minimal HSNBL citation stack, use [Grachev2005], [GrachevEtAl2007], [Grachev2012], and [GryanikEtAl2020].

## 2. Analytical / Inversion References from Your Existing ABL Line

These are important because they connect the empirical similarity functions to analytic inversion, Richardson-number structure, and closure design.

```bibtex
@article{EnglandMcNider1995,
  author = {England, D. E. and McNider, R. T.},
  title = {Stability Functions Based upon Shear Functions},
  journal = {Boundary-Layer Meteorology},
  volume = {74},
  number = {1-2},
  pages = {113--130},
  year = {1995},
  doi = {10.1007/BF00715711}
}
```

Use this when the argument shifts from fitting $\phi_m$ or $\phi_h$ directly to deriving or inverting functions of Richardson number, especially where neutral curvature, monotonicity, or analytical invertibility matters.

This paper supports the closure-analysis side of the project. It still does not directly establish the Gegenbauer interpretation, but it does justify treating the similarity forms as analytically structured rather than purely empirical fits.

## 3. Mathematical References for Ultraspherical / Gegenbauer Structure

These are the references that support the mathematical language: generating functions, recurrences, orthogonality, Sturm-Liouville structure, and hyperspherical interpretation.

### 3.1 Core orthogonal-polynomial references

```bibtex
@book{Szego1975,
  author = {Szego, Gabor},
  title = {Orthogonal Polynomials},
  publisher = {American Mathematical Society},
  address = {Providence, RI},
  edition = {4th},
  year = {1975}
}

@book{AbramowitzStegun1972,
  editor = {Abramowitz, Milton and Stegun, Irene A.},
  title = {Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables},
  publisher = {U.S. Government Printing Office},
  address = {Washington, DC},
  year = {1972}
}

@book{AndrewsAskeyRoy1999,
  author = {Andrews, George E. and Askey, Richard and Roy, Ranjan},
  title = {Special Functions},
  publisher = {Cambridge University Press},
  address = {Cambridge},
  year = {1999}
}
```

Use these for:

- the generating function

$$
\sum_{n=0}^{\infty} C_n^{(\lambda)}(x) t^n = (1 - 2xt + t^2)^{-\lambda},
$$

- the recurrence relations for $C_n^{(\lambda)}$,
- the orthogonality weight $(1-x^2)^{\lambda - 1/2}$,
- the identification of Legendre as the special case $\lambda = 1/2$.

These books are the right references for any equation-level mathematical statement in the ultraspherical section.

### 3.2 Generating-function and hyperspherical references

```bibtex
@article{Cohl2013,
  author = {Cohl, Howard S.},
  title = {On a Generalization of the Generating Function for {Gegenbauer} Polynomials},
  journal = {Integral Transforms and Special Functions},
  volume = {24},
  number = {11},
  pages = {916--927},
  year = {2013},
  doi = {10.1080/10652469.2012.761613}
}

@book{SteinWeiss1971,
  author = {Stein, Elias M. and Weiss, Guido},
  title = {Introduction to Fourier Analysis on Euclidean Spaces},
  publisher = {Princeton University Press},
  address = {Princeton, NJ},
  year = {1971}
}
```

Use these when you want to justify the language of hyperspherical harmonics, spectral expansions, and the geometric meaning of the Gegenbauer parameter.

One correction to the older draft: Stein and Weiss is a book, not an article, and it should not be cited as a direct atmospheric reference. It is a mathematical support citation only.

## 4. Related Modern Turbulence-Generalization Reference

```bibtex
@article{StiperskiCalaf2023,
  author = {Stiperski, I. and Calaf, M.},
  title = {Generalizing {Monin-Obukhov} Similarity Theory for Complex Atmospheric Turbulence},
  journal = {Physical Review Letters},
  volume = {130},
  pages = {124001},
  year = {2023},
  doi = {10.1103/PhysRevLett.130.124001}
}
```

This paper is useful as a conceptual bridge. It does not discuss ultraspherical polynomials, but it does strengthen the argument that classical MOST needs extra state variables or extra structure when the turbulence departs from equilibrium and isotropy. That makes it a good high-level support citation for why an expanded spectral state variable may be physically reasonable.

## 5. What the Literature Supports vs. What Is Currently Your Proposal

This distinction should stay explicit in papers, notes, and talks.

### 5.1 Supported directly by the literature

- MOST similarity functions are empirical but analytically structured.
- SHEBA shows that the very-stable Arctic boundary layer departs from classical local similarity.
- Stable and very-stable regimes exhibit intermittency, regime transitions, and effective breakdown of local closure assumptions.
- Gegenbauer polynomials have the exact generating-function and orthogonality properties needed for spectral representations on bounded intervals.

### 5.2 Proposed in your framework and should be cited as such

- identifying the unstable MOST power laws $(1-b\zeta)^{-1/4}$ and $(1-b\zeta)^{-1/2}$ as direct ultraspherical generating-function objects with atmospheric physical meaning,
- interpreting $\lambda = 1/4$ for momentum and $\lambda = 1/2$ for heat as effective turbulent spectral coordinates rather than only algebraic exponents,
- linking $\lambda$ to an effective turbulent dimension through

$$
\lambda = \frac{d-2}{2}, \qquad d = 2 + 2\lambda,
$$

- treating HSNBL residual structure as a low-order Gegenbauer correction about a SHEBA baseline,
- interpreting evolution of $\lambda^*$ as a regime diagnostic for anisotropy, intermittency, or dimensional collapse.

Those statements may be mathematically consistent and physically promising, but they are not established claims in the cited atmospheric literature. They should be presented as the new contribution.

## 6. Recommended Citation Stack by Use Case

### 6.1 If the sentence is about stable boundary-layer physics

Use:

- BusingerEtAl1971
- Dyer1974
- Grachev2005
- GrachevEtAl2007
- Grachev2012
- GryanikEtAl2020

### 6.2 If the sentence is about the mathematics of the spectral basis

Use:

- Szego1975
- AbramowitzStegun1972
- AndrewsAskeyRoy1999
- Cohl2013
- SteinWeiss1971

### 6.3 If the sentence is about analytical closure structure or inversion

Use:

- EnglandMcNider1995
- Grachev2014

### 6.4 If the sentence is about why MOST needs generalization

Use:

- Grachev2012
- GryanikEtAl2020
- StiperskiCalaf2023

## 7. References to Keep as Internal Working Notes, Not External Support

The following are valuable inside this repo but should not be cited as if they are peer-reviewed external evidence:

- [notes/bridge the gap.md](/Users/davidengland/Documents/GitHub/ABL/notes/bridge%20the%20gap.md)
- [drafts/ultraspherical_subsection.md](/Users/davidengland/Documents/GitHub/ABL/drafts/ultraspherical_subsection.md)
- [notes/Ultraspherical Approach for the HSNBL.md](/Users/davidengland/Documents/GitHub/ABL/notes/Ultraspherical%20Approach%20for%20the%20HSNBL.md)

These documents contain the project-specific synthesis and are the right place to develop the new claims before they are turned into a formal manuscript.

## 8. Minimal Shortlist for a Paper Introduction

If you need a concise reference list for a first manuscript draft, start with these ten:

```text
BusingerEtAl1971
Dyer1974
Hogstrom1988
Grachev2005
GrachevEtAl2007
Grachev2012
GryanikEtAl2020
EnglandMcNider1995
Szego1975
Cohl2013
```

That set is compact, defensible, and covers the atmospheric, closure, and mathematical sides without overclaiming that the ultraspherical interpretation is already standard literature.
