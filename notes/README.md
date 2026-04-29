# notes/ — Research Working Notes

Quick-reference sheets, parameter tables, and derivation notes for active research.  
These are working documents; polished versions migrate to [drafts/](../drafts/) or [manuscripts/](../manuscripts/).

---

## Index

| File | Description | Last updated |
|---|---|---|
| [parameters.md](parameters.md) | MOST parameter tables: Businger (1971), Dyer (1974), Högström (1988), McNider — b_m, b_h, a_h, κ values | Apr 2026 |
| [parameters2.md](parameters2.md) | Extended parameter reference including SHEBA, Cheng-Brutsaert, critical Ri_c values | Apr 2026 |
| [Legendre.md](Legendre.md) | Legendre and Gegenbauer (ultraspherical) orthogonality, generating functions, Rodrigues formulae, CBC equator identity | Apr 2026 |
| [Central Binomial UBL heat shear.md](Central%20Binomial%20UBL%20heat%20shear.md) | CBC series expansion for UBL heat-shear asymptotic coupling | Apr 2026 |
| [SBL-IBEx.md](SBL-IBEx.md) | Quick-reference sheet: IBEx (Internal/Boundary/External) series for SBL ζ/(1+βζ) mappings | Apr 2026 |
| [Rig.md](Rig.md) | Gradient Richardson number: definitions, exact forms, asymptotic limits | — |
| [Connections.csv](Connections.csv) | Cross-reference matrix of key results across files | — |

---

## Key Equations (at a glance)

**Power-law MOST functions (SBL):**
$$\phi_m = (1 - b_m \zeta)^{-1/4}, \quad \phi_h = a_h^{-1}(1 - b_h \zeta)^{-1/2}$$

**Gegenbauer generating function (λ = 1/4):**
$$(1 - 2xt + t^2)^{-1/4} = \sum_{n=0}^\infty C_n^{(1/4)}(x)\,t^n$$

**CBC equator identity:**
$$P_{2n}(0) = (-1)^n \frac{\binom{2n}{n}}{4^n}$$

**Safeguarded Newton derivative (for ζ inversion):**
$$\frac{dRi_g}{d\zeta} = a_h^{-1}\sqrt{\frac{1 - b_m\zeta}{1 - b_h\zeta}} \cdot \left[\frac{1}{1-b_h\zeta} - \frac{b_m}{2(1-b_m\zeta)}\right] \cdot$$

---

## Related

- [examples/](../examples/) — Fortran implementations of these identities
- [drafts/ultraspherical_subsection.md](../drafts/ultraspherical_subsection.md) — Full mathematical development
- [CANONICAL_GLOSSARY.md](../CANONICAL_GLOSSARY.md) — Symbol conventions
