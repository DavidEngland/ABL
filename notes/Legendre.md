Let me search for any existing Legendre content and review the CBC homework before building the advanced problem set.

Searched for text `Legendre|P_n|legendre`, 20 results

Read [](file:///Users/davidengland/Documents/GitHub/ABL/hw/Legendre.md)

No Legendre content in the repo — the file is empty. Let me work through the math carefully before writing, since the derivation links have multiple layers.

Read [](file:///Users/davidengland/Documents/GitHub/ABL/hw/CBC.md#1-1), lines 1 to 60

Read [](file:///Users/davidengland/Documents/GitHub/ABL/notes/parameters.md#80-80), lines 80 to 122

Good. Now I have everything needed. The math chain is:

1. **$\phi_h$:** CBC coefficients $= (-1)^n P_{2n}(0)$; unstable branch is the Legendre generating function at $x=0$ with $t=\sqrt{b_h\eta}$
2. **$\phi_m$:** same but Gegenbauer (ultraspherical) $C_n^{(1/4)}$ replaces Legendre $C_n^{(1/2)}$
3. **$\phi_h=\phi_m^2$** is a Gegenbauer convolution identity
4. **Legendre projection** of profiles via domain mapping $\mu = 1-2b|\zeta|$
5. **Gauss-Legendre quadrature** for bulk Ri integrals, nodes from Legendre zeros

Writing the file now:

Generating patch (378 lines) in [](file:///Users/davidengland/Documents/GitHub/ABL/hw/Legendre.md)


```

```
Done. Here's what Legendre.md contains — a progression from master's to PhD research level:

---

### Problem structure

| # | Level | Topic | Points |
|---|---|---|---|
| 1 | MS | $P_{2n}(0)$ via Rodrigues; CBC coefficients = equatorial Legendre values; $\phi_h$ as the Legendre generating function at $\theta=\pi/2$ | 25 |
| 2 | MS | Gegenbauer $C_n^{(1/4)}$ for $\phi_m$; Legendre as special case $\lambda=1/2$; coefficient growth comparison | 25 |
| 3 | MS/PhD | Squaring identity $\phi_h = \phi_m^2$ as Gegenbauer linearization (Clebsch-Gordan); term-by-term verification | 20 |
| 4 | PhD | Domain mapping $\mu = 1-2b|\zeta|$; analytic Legendre spectral coefficients; Python GL projection; physical interpretation of modes | 30 |
| 5 | PhD | Gauss-Legendre quadrature for $Ri_b$; GL nodes = zeros of $P_N$; connection back to CBC structure | 25 |
| 6 | PhD | Sturm-Liouville ODE in stability space; eigenfunction expansion of $Ri_g$; truncation error in Ri-inversion | 20 |
| 7 | Research | 3D spherical harmonics and horizontal ABL inhomogeneity; spectral Prandtl number as mode mixing; GL nodes as optimal vertical level placement | open |

**The key mathematical chain** that runs through all problems — the chain linking MOST physics → CBC → Legendre generating function → Gegenbauer product formula — is displayed as a summary equation block at the end.

Made changes.