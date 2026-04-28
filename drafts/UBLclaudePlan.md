Checking the current draft and related context before writing.

Read [](file:///Users/davidengland/Documents/GitHub/ABL/drafts/UBL%20heat%E2%80%93shear%20asymptotic.md)

Read [](file:///Users/davidengland/Documents/GitHub/ABL/notes/Central%20Binomial%20UBL%20heat%20shear.md#1-1), lines 1 to 80

Good context. The draft covers only $\phi_h$ and uses informal prose. I'll write a full paper-section that brings in $\phi_m$, tightens notation throughout, establishes the pair as UBL shear functions, derives the $Ri_g=\zeta$ result, and covers all regimes—then write it directly into the draft file.

Searched for regex `Businger|Dyer|Paulson|phi_m|phi_h|beta.*16|16.*zeta` (`**/*.bib`), 20 results

Read [](file:///Users/davidengland/Documents/GitHub/ABL/ABL_refs.bib#1-1), lines 1 to 20

Good — citation keys confirmed as `Businger1971`, `Dyer1974`, `Paulson1970`. Now I'll write the full tightened section.

Run in Terminal

The terminal is currently disabled, so I can't write directly to the file. Here is the complete, ready-to-paste notation-tightened section for your draft:

---

```markdown
<!-- UBL heat–shear similarity functions: analytic structure and asymptotics -->
<!-- For inclusion in turbulence-closure paper; LaTeX-compatible notation -->

## UBL Shear-Similarity Functions: Analytic Structure and Asymptotic Regimes

### 2.1 Definitions and the Businger–Dyer pair

Following the unstable surface-layer forms of \citet{Businger1971} and \citet{Dyer1974},
define the **UBL shear-similarity functions** for momentum and heat as

$$
\phi_m(\zeta) = (1 - \beta\,\zeta)^{-\alpha_m}, \qquad
\phi_h(\zeta) = (1 - \beta\,\zeta)^{-\alpha_h},
\qquad \zeta \le 0,
\tag{1}
$$

where $\zeta = z/L$ is the Obukhov stability parameter, $\beta = 16$ is the
Businger–Dyer empirical coefficient, and the exponents

$$
\alpha_m = \tfrac{1}{4}, \qquad \alpha_h = \tfrac{1}{2} = 2\alpha_m
\tag{2}
$$

are constrained by inertial-subrange similarity arguments.
The heat function $\phi_h = (1-\beta\zeta)^{-1/2}$ is the outer square root;
$\phi_m = (1-\beta\zeta)^{-1/4}$ is its double square root (smaller exponent,
weaker singular response).

### 2.2 Fundamental algebraic relationship

Because $\alpha_h = 2\alpha_m$ and both functions share the coefficient $\beta$,

$$
\boxed{\phi_h(\zeta) = \phi_m(\zeta)^2.}
\tag{3}
$$

Equation (3) is an exact identity, not an approximation; it holds for all
$\zeta \le 0$.

### 2.3 Richardson-number inversion

The gradient Richardson number satisfies the MOST relation

$$
Ri_g = \zeta\,\frac{\phi_h}{\phi_m^2}.
\tag{4}
$$

Substituting (3):

$$
Ri_g = \zeta\,\frac{\phi_h(\zeta)}{\phi_h(\zeta)} = \zeta.
\tag{5}
$$

**Equation (5) is the principal structural result**: for the Businger–Dyer exponent
pair the $Ri_g$–$\zeta$ conversion is trivially invertible with no iterative
solver on the unstable branch.  This degeneracy is an exact consequence of
$\alpha_h = 2\alpha_m$; perturbing either exponent breaks it.

To quantify sensitivity, let $\alpha_h = 2\alpha_m + \varepsilon$ with
$|\varepsilon| \ll 1$ and $\beta_m = \beta_h = \beta$.  Near neutrality,

$$
Ri_g = \zeta\bigl[1 - \varepsilon\,\beta\,\zeta + O(\zeta^2,\,\varepsilon^2)\bigr].
\tag{6}
$$

The leading correction is $O(\varepsilon\beta\zeta)$ and remains negligible in the
near-neutral regime for any modest $\varepsilon$.

### 2.4 Near-neutral expansion: central-binomial series

For $|\beta\zeta| < 1$ (convergence radius $|\zeta| < \beta^{-1} = 1/16$), the
generalized binomial identity gives

$$
\phi_x(\zeta) = \sum_{n=0}^{\infty} c_n^{(\alpha_x)}\,\zeta^n,
\qquad
c_n^{(\alpha_x)} = \frac{\Gamma(\alpha_x + n)}{\Gamma(\alpha_x)\,n!}\,\beta^n.
\tag{7}
$$

For $\alpha_h = \tfrac{1}{2}$ the Gamma factors collapse to central-binomial
coefficients:

$$
\phi_h(\zeta)
= \sum_{n=0}^{\infty} \binom{2n}{n} 4^n\,\zeta^n
= 1 + 8\zeta + 96\zeta^2 + 1280\zeta^3 + 17920\zeta^4 + \cdots,
\tag{8}
$$

with coefficients satisfying the two-term recurrence

$$
c_{n+1}^{(1/2)} = c_n^{(1/2)}\,\frac{8(2n+1)}{n+1}, \qquad c_0^{(1/2)} = 1.
\tag{9}
$$

For $\alpha_m = \tfrac{1}{4}$ the sequence $c_n^{(1/4)}$ grows more slowly,
consistent with the milder singularity of $\phi_m$, and satisfies the analogous
recurrence with $\alpha_m$ in place of $\alpha_h$.

**Advantages of the central-binomial representation:**
(i) exact closed-form coefficients via a simple two-term recurrence;
(ii) radius of convergence set transparently by the branch point at
$\zeta = \beta^{-1}$;
(iii) natural basis for $[m/n]$ Padé extensions across the neutral-to-stable
transition;
(iv) parity separation $\phi_h = \phi_h^{(\mathrm{even})} + \phi_h^{(\mathrm{odd})}$
decouples slope and curvature tuning in empirical calibration.

### 2.5 Asymptotic behavior: strongly unstable (UBL) limit

For $\zeta < 0$ set $\eta = -\zeta > 0$.  As $\eta \to \infty$:

$$
\phi_m(\zeta) = (\beta\eta)^{-1/4}
\!\left[1 - \tfrac{1}{4}(\beta\eta)^{-1}
+ \tfrac{5}{32}(\beta\eta)^{-2} - \cdots\right],
\tag{10a}
$$

$$
\phi_h(\zeta) = (\beta\eta)^{-1/2}
\!\left[1 - \tfrac{1}{2}(\beta\eta)^{-1}
+ \tfrac{3}{8}(\beta\eta)^{-2} - \cdots\right].
\tag{10b}
$$

Both series are regular inverse-power hierarchies; no oscillatory or exponential
terms appear.  The heat function decays as $|\zeta|^{-1/2}$ while $\phi_m$ decays
as $|\zeta|^{-1/4}$, so scalar-gradient closures built on $\phi_h$ are
intrinsically less sensitive to large instabilities than momentum-based closures.

Leading-order regime summary:

| Regime | $\zeta$ | $\phi_m$ | $\phi_h$ |
|---|---|---|---|
| Strong convection | $\zeta \to -\infty$ | $(\beta|\zeta|)^{-1/4}$ | $(\beta|\zeta|)^{-1/2}$ |
| Near neutral | $|\zeta| \ll \beta^{-1}$ | $1 + \tfrac{\beta}{4}\zeta + \cdots$ | $1 + \tfrac{\beta}{2}\zeta + \cdots$ |
| Stable onset | $0 < \zeta \ll \beta^{-1}$ | $\nearrow$, regular | $\nearrow$, faster |

### 2.6 Stable and laminar-limit continuation

For $\zeta > 0$ both functions increase monotonically toward the branch singularity
at $\zeta = \beta^{-1}$.  Beyond that point the expressions become complex; the
Businger–Dyer forms cannot serve as universal stable-regime closures.

Practical implementations require a piecewise or rational continuation matched at
a transition stability $\zeta_s \in (0,\,\beta^{-1})$ with at least $C^1$
continuity.  Candidate forms include the linear Beljaars–Holtslag expressions,
capped algebraic growth, or $[m/n]$ Padé approximants constructed from (8) with
no positive-real poles imposed.

### 2.7 Numerical implementation

1. **Near-neutral evaluation**: use

   $$
   \phi_x(\zeta) = \exp\!\bigl(-\alpha_x\,\mathrm{log1p}(-\beta\zeta)\bigr)
   \tag{11}
   $$

   to eliminate subtractive cancellation when $|\beta\zeta| \ll 1$.

2. **Series coefficients**: generate $c_n^{(\alpha_x)}$ via the recurrence (9);
   apply Horner's method for $O(N)$ flops in $N$-term truncations.

3. **Inversion**: on the unstable branch use $\zeta = Ri_g$ directly from (5).
   On a stable-continuation branch, solve $Ri_g = \zeta\,f(\zeta)$ by Newton
   iteration initialized from (6); convergence is rapid because the stable form
   is monotone and bounded.

4. **Exponent generalization**: if non-standard $(\alpha_m,\alpha_h)$ are used for
   sensitivity studies, verify $\alpha_h = 2\alpha_m$ before invoking (5);
   otherwise initialize Newton iteration from the perturbative estimate (6).
```

---

Key notation choices made:
- $\alpha_m, \alpha_h$ for exponents (separates coefficient from exponent cleanly)
- $\beta = 16$ explicitly named as Businger–Dyer coefficient (cite keys `Businger1971`, `Dyer1974`)
- $\eta = -\zeta$ for the unstable variable (avoids sign confusion in asymptotics)
- Equations (10a,b) give companion expansions for both functions explicitly
- The $\phi_h = \phi_m^2$ identity and $Ri_g = \zeta$ result are boxed as principal results
- Gamma-function form (7) unifies both exponent cases cleanly before specializing to central-binomial (8)