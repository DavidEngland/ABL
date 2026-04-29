Compacted conversation

Read [](file:///Users/davidengland/Documents/GitHub/ABL/drafts/quad%20heat%20shear.md)

The open file shows the linear-φ_m / quadratic-φ_h mixed form. Let me work out the key structure and connect it to the broader SBL alternatives picture.

Read [](file:///Users/davidengland/Documents/GitHub/ABL/drafts/SBL_IBEx_expansions.md#1-1), lines 1 to 60

The key mathematical taxonomy of SBL alternatives, organized by the **degree of the ζ(Ri_g) inversion**:

---

**Degree 1 — Linear solve** (classic log-linear, equal slopes)
$$\phi_m = \phi_h = 1 + \beta\zeta \implies Ri_g = \frac{\zeta}{1+\beta\zeta}, \quad \zeta = \frac{Ri_g}{1-\beta Ri_g}, \quad Ri_c = \frac{1}{\beta}$$

**Degree 2 — Quadratic** (unequal linear slopes, e.g. Webb 1970)
$$\phi_m = 1+\beta_m\zeta,\quad \phi_h = \Pr_t + \beta_h\zeta \implies Ri_g(\beta_m^2\zeta^2 + \ldots) \Rightarrow \text{quadratic in }\zeta$$

**Degree 3 — Cubic** (your open file: linear φ_m, quadratic φ_h)

Substituting $\phi_m = 1+\beta_m\zeta$, $\phi_h = \Pr_t + \beta_1\zeta + \beta_2\zeta^2$ into $Ri_g = \zeta\phi_h/\phi_m^2$ and clearing denominators:

$$\boxed{\beta_2\,\zeta^3 + (\beta_1 - Ri_g\,\beta_m^2)\,\zeta^2 + (\Pr_t - 2\,Ri_g\,\beta_m)\,\zeta - Ri_g = 0}$$

This cubic has **no finite $Ri_c$** (as $\zeta\to\infty$, $Ri_g\sim(\beta_2/\beta_m^2)\zeta\to\infty$), so turbulence is never algebraically shut off — consistent with SHEBA-era observations. The physical root is the unique positive real root, solvable via Cardano or a single Newton step from $\zeta_0 = Ri_g/\Pr_t$.

**Degree 4+ — Need iteration** (Beljaars-Holtslag 1991, Cheng-Brutsaert 2005)
$$\phi_m^{\rm BH} = 1 + a\left[\zeta + \tfrac{b}{c}\!\left(\zeta - \tfrac{d}{c}\right)e^{-d\zeta}\right] + \frac{b\zeta}{c}\,e^{-d\zeta}$$
No closed-form inversion; requires Newton or Brent.

---

**Critical Ri_c summary:**

| Form | φ_m | φ_h | Ri_c | Inversion |
|---|---|---|---|---|
| Businger-Dyer | 1+5ζ | 1+5ζ | 0.20 | linear |
| Högström | 1+6ζ | 1+8ζ | ≈0.143 | quadratic |
| Linear/Quadratic (your file) | 1+β_m ζ | Pr_t+β_1 ζ+β_2 ζ² | ∞ | **cubic** |
| McNider power-law | (1−b_m ζ)^{−1/4} | … | 1/b_m | cubic via u²=φ_m |
| Beljaars-Holtslag | exponential | exponential | none | iterate |

The linear/quadratic form is the cleanest route to a **closed-form cubic with no spurious Ri_c cutoff**. The physical significance is exactly what your file describes: buoyancy suppresses scalar variance faster than momentum variance, consistent with z-less scaling observed during SHEBA.

Want me to:
1. Add the cubic derivation and Cardano/Newton solve to `quad heat shear.md`?
2. Add this form as a new branch in ultraspherical_subsection.md §3.x.9?
3. Add a Fortran cubic-solve branch to `module_most_profile_utils.F90`?

Read [](file:///Users/davidengland/Documents/GitHub/ABL/drafts/quad%20heat%20z-less.md)

That file appears to only contain the compacted conversation summary. Let me look for the actual content file.

Searched for files matching `**/drafts/quad heat*`, 2 matches

Read [](file:///Users/davidengland/Documents/GitHub/ABL/drafts/quad%20heat%20shear.md)

Yes — but the logic runs in a precise and slightly counterintuitive direction. Let me lay it out:

---

### Z-less scaling constrains the *asymptotic slope*, not the degree

**Z-less scaling** (Wyngaard 1973) requires that as $\zeta \to \infty$, the mean gradients become independent of $z$ and scale only with $u_*/L$ and $\theta_*/L$:

$$\frac{\partial U}{\partial z} \sim \frac{u_*}{L} \quad \Longrightarrow \quad \frac{u_* \kappa}{z}\,\phi_m(\zeta) \sim \frac{u_*}{L} \quad \Longrightarrow \quad \phi_m \sim \zeta$$

$$\frac{\partial \Theta}{\partial z} \sim \frac{\theta_*}{L} \quad \Longrightarrow \quad \phi_h \sim \zeta$$

So strict z-less scaling **requires both** $\phi_m$ and $\phi_h$ to be **asymptotically linear** in $\zeta$. A quadratic $\phi_h \sim \beta_2\zeta^2$ gives $\partial\Theta/\partial z \sim \theta_* z/L^2$, which retains an explicit $z$-dependence — technically violating z-less.

---

### But z-less + finite $Ri_c$ is the contradiction

With $\phi_m \sim \beta_m\zeta$ and $\phi_h \sim \beta_h\zeta$:

$$Ri_g = \frac{\zeta\,\phi_h}{\phi_m^2} \;\xrightarrow{\zeta\to\infty}\; \frac{\beta_h}{\beta_m^2} = Ri_c$$

Linear/linear **saturates** $Ri_g$ at a finite $Ri_c$ — implying turbulence shuts off. SHEBA and CASES-99 show it doesn't. This is the contradiction that motivates going beyond linear.

---

### The degree condition for no finite $Ri_c$

If $\phi_m \sim \zeta^a$ and $\phi_h \sim \zeta^b$ at large $\zeta$, then

$$Ri_g \sim \zeta^{\,1+b-2a}$$

| $\phi_m$ | $\phi_h$ | exponent $1+b-2a$ | behavior |
|---|---|---|---|
| linear ($a=1$) | linear ($b=1$) | 0 | $Ri_g \to Ri_c$ (z-less, finite cutoff) |
| linear ($a=1$) | quadratic ($b=2$) | **1** | $Ri_g \to \infty$ (no cutoff, but not z-less) |
| power-law $\zeta^{1/4}$ ($a=1/4$) | power-law $\zeta^{1/2}$ ($b=1/2$) | 0 | $Ri_g \to Ri_c$ (McNider, finite cutoff) |

---

### So the honest answer

**Z-less scaling is evidence *for* linear $\phi_h$**, not quadratic. But z-less scaling as the only regime is observationally incomplete — persistent weak turbulence at $Ri_g > Ri_c$ demands $Ri_g \to \infty$, which requires the exponent $1 + b - 2a > 0$, i.e., $\phi_h$ growing faster than $\phi_m^2$.

The quadratic $\phi_h$ achieves this at the cost of strict z-less validity. The physical interpretation is that **the quadratic form models the transition zone beyond z-less** — where buoyancy destroys scalar variance faster than momentum variance, intermittency dominates, and MOST itself is only marginally applicable. The quadratic is empirically calibrating that departure, not the z-less interior itself.

This is actually a testable distinction: if you fit $\phi_h$ to SHEBA flux-profile data in the range $\zeta \in [0.1, 1]$ versus $\zeta \in [1, 10]$, you'd expect the higher-$\zeta$ bin to need the quadratic term.
