# Log Expansions and Logarithmic Profiles in the ABL
*Brainstorming notes — for discussion with McNider and Biazar*
*March 2026*

---

## 1. The Core Expansion

For a height $H \approx 10\,\text{m}$ (canopy height, displacement height, or reference level), the
key identity is

$$
\ln(z + H) = \ln z + \ln\!\left(1 + \frac{H}{z}\right)
$$

Using the standard Taylor series $\ln(1+x) = x - \tfrac{x^2}{2} + \tfrac{x^3}{3} - \cdots$ valid for $|x|<1$:

$$
\boxed{
\ln(z+H) = \ln z + \frac{H}{z} - \frac{H^2}{2z^2} + \frac{H^3}{3z^3} - \cdots
}
$$

**Note on signs:** The first correction is $+H/z$ (positive), not $-H/z$.  This makes sense
because $z+H > z \Rightarrow \ln(z+H) > \ln z$.  The series alternates in sign from the second
term onward.

Convergence requires $|H/z| < 1$, i.e., $z > H$.  The series diverges for $z < H$.

---

## 2. Three Height Regimes

### Regime 1 — Far-field: $z \gg H$ (e.g., $z > 50\,\text{m}$)

The asymptotic is excellent; retain only the leading correction:

$$
\ln(z+H) \approx \ln z + \frac{H}{z}
$$

Relative correction $\sim H/z \ll 1$.  For $z = 100\,\text{m}$, $H = 10\,\text{m}$: correction
is 10%.  For $z = 1000\,\text{m}$: 1%.  Higher-order terms ($H^2/2z^2$, ...) are negligible.

### Regime 2 — Transition: $z \sim H$ (e.g., $5 \lesssim z \lesssim 30\,\text{m}$)

**Neither asymptotic applies.**  Must use the exact $\ln(z+H)$.  The Taylor series is slowly
convergent at best; at $z = H$ each successive term is half the previous.  This is precisely
the roughness sublayer and the lowest few model levels — the most critical region for surface
flux parameterization.

### Regime 3 — Near-surface: $z \ll H$

Expand the *other* way: factor out $H$:

$$
\ln(z+H) = \ln H + \ln\!\left(1 + \frac{z}{H}\right) \approx \ln H + \frac{z}{H} - \frac{z^2}{2H^2} + \cdots
$$

Leading order: **the profile is approximately linear in $z$**, not logarithmic.

$$
\ln(z+H) \approx \ln H + \frac{z}{H}, \qquad z \ll H
$$

This is important: the log profile "flattens out" to a near-linear form very close to the surface.

---

## 3. Connection to Standard Logarithmic Profiles

The textbook neutral log profiles are

$$
u(z) = \frac{u_*}{\kappa} \ln\!\frac{z}{z_0}, \qquad
\theta(z) - \theta_s = \frac{\theta_*}{\kappa} \ln\!\frac{z}{z_0}
$$

where $z_0 \sim 0.01$–$0.1\,\text{m}$ is the roughness length.  The ratio $H/z_0 \gg 1$, so
the structure of $\ln(z/z_0)$ vs.\ $\ln((z+H)/z_0)$ matters.

$$
\ln\!\frac{z+H}{z_0} = \ln\!\frac{z}{z_0} + \frac{H}{z} - \frac{H^2}{2z^2} + \cdots
$$

- For $z \gg H$: the profiles with and without the offset $H$ agree to leading order.
- For $z \sim H$ (first model level?): the correction $H/z \sim \mathcal{O}(1)$ — a 100% error if
  the offset is ignored.
- For $z \ll H$: $\ln(z/z_0)$ may even go *negative* (if $z < z_0$, which is unphysical anyway),
  while $\ln((z+H)/z_0)$ stays well-behaved.

---

## 4. Displacement Height Formulation

The standard modified log profile shifts the origin by a displacement height $d$:

$$
u(z) = \frac{u_*}{\kappa} \ln\!\frac{z - d}{z_0}
$$

This is $\ln((z-d)/z_0)$, i.e., $H \to -d$.  Substituting $-d$ for $H$:

$$
\ln\!\frac{z-d}{z_0} = \ln\!\frac{z}{z_0} - \frac{d}{z} - \frac{d^2}{2z^2} - \cdots
$$

(All terms negative because the series now has $x = -d/z < 0$ and alternate signs cancel the
alternation.) The profile is **reduced** relative to the no-displacement case, as expected.

**Danger near $z = d$:** The log diverges to $-\infty$ as $z \to d^+$, which is unphysical.
The expansion is completely broken below $z = d$.  This is why MOST is only applied above
the roughness sublayer ($z > 2$–$5d$), exactly the near-$H$ regime flagged in Regime 2 above.

---

## 5. Gradient and Bulk Richardson Number Implications

The wind shear derived from the log profile is

$$
\frac{\partial u}{\partial z} = \frac{u_*}{\kappa z}
$$

If the "true" argument is $\ln((z+H)/z_0)$, the shear is

$$
\frac{\partial u}{\partial z} = \frac{u_*}{\kappa(z+H)}
$$

Taking the ratio:

$$
\frac{(\partial u/\partial z)_{\text{with }H}}{(\partial u/\partial z)_{\text{without }H}}
= \frac{z}{z+H} = \frac{1}{1 + H/z}
\approx 1 - \frac{H}{z} + \left(\frac{H}{z}\right)^2 - \cdots
$$

- At $z = 100\,\text{m}$: shear is underestimated by 10% if $H$ is ignored.
- At $z = 10\,\text{m}$ ($= H$): shear is underestimated by 50%.

The **bulk Richardson number** computed over a layer from $z_1$ to $z_2$ involves
$[\ln(z_2/z_0) - \ln(z_1/z_0)]^2 = \ln^2(z_2/z_1)$.  If either level is near $H$, the log ratio
is sensitive to whether $H$ is included.

---

## 6. Finite-Difference / Grid Discretization Link

This directly connects to the grid-curvature bias work.  A numerical model with levels
$z_k$, $z_{k+1} = z_k + \Delta z$ evaluates

$$
\frac{u_{k+1} - u_k}{\Delta z} = \frac{u_*}{\kappa \Delta z} \ln\!\frac{z_k + \Delta z}{z_k}
= \frac{u_*}{\kappa \Delta z} \ln\!\left(1 + \frac{\Delta z}{z_k}\right)
$$

Expanding (with $\Delta z$ playing the role of $H$):

$$
= \frac{u_*}{\kappa \Delta z}\left[\frac{\Delta z}{z_k} - \frac{(\Delta z)^2}{2z_k^2} + \cdots\right]
= \frac{u_*}{\kappa z_k}\left[1 - \frac{\Delta z}{2z_k} + \frac{(\Delta z)^2}{3z_k^2} - \cdots\right]
$$

The leading **grid-curvature error** is $-\Delta z/(2z_k)$.  Near the surface where $z_k \sim \Delta z$,
this is a $\sim 50\%$ bias in the computed shear — the same near-$H$ regime where the log
expansion breaks down.

> **Key takeaway:** The asymptotic failure of $\ln(z+H) \approx \ln z$ and the large
> discretization error of finite-difference shear estimates are *the same mathematical phenomenon*:
> both are failures of the first-order Taylor expansion of $\ln(1 + \epsilon)$ when $\epsilon$ is
> not small.

---

## 7. Practical Numbers ($H = 10\,\text{m}$)

| $z$ (m) | $H/z$ | $\ln(z+H)/\ln(z)$ error | Leading correction: $+H/z$ | 2nd order: $-H^2/(2z^2)$ |
|---------|--------|------------------------|---------------------------|--------------------------|
| 5       | 2.0    | series diverges        | —                         | —                        |
| 10      | 1.0    | converges slowly       | $+1.000$                  | $-0.500$                 |
| 20      | 0.5    | 11% error in $\ln z$   | $+0.500$                  | $-0.125$                 |
| 50      | 0.2    | 2.5% error             | $+0.200$                  | $-0.020$                 |
| 100     | 0.1    | 0.5% error             | $+0.100$                  | $-0.005$                 |
| 200     | 0.05   | 0.12% error            | $+0.050$                  | $-0.001$                 |

"Error in $\ln z$" = $|\ln(z+H) - \ln z| / |\ln z|$; uses $z_0 = 0.1\,\text{m}$ as reference to
give $\ln(z/z_0)$.

---

## 8. Open Questions / Things to Raise with McNider & Biazar

1. **Which specific profile forms produce $\ln(z+H)$?**
   - Displaced log profile $(z-d)$ vs.\ shifted origin $(z+H)$?
   - Is $H$ the blending height, canopy height, or first model level?

2. **Nocturnal stable BL:**  The nocturnal BL may be only 20–50 m deep.  Then $H \approx 10$ m is
   a substantial fraction of BL depth, and the asymptotic regime is never really reached.
   Does the log profile even apply, or does the $\ln(z+H) \approx \ln H + z/H$ (linear) form
   describe the shear better?

3. **Error in $Ri_g$ near the surface:**  If the shear is biased by $\Delta z/(2z_k)$ due to grid
   curvature, the Richardson number errors could trigger spurious regime transitions (laminar ↔
   turbulent).  Is this the same issue already identified in the curvature correction papers?

4. **Matched asymptotic approach:**  Could we construct a composite expansion valid across all
   three regimes by matching the far-field expansion and the near-surface expansion at some
   intermediate height $z \sim H$?  Something like:
   $$
   \ln(z+H) \approx \begin{cases}
     \ln H + z/H & z \ll H \\
     \ln\sqrt{z H} + \text{(cross terms)} & z \sim H \quad (\text{transition}) \\
     \ln z + H/z & z \gg H
   \end{cases}
   $$

5. **Temperature and humidity:**  Same expansion applies.  But the potential temperature profile
   under stable conditions has a much steeper gradient — does the near-$H$ breakdown affect
   flux-gradient relationships more strongly for heat than for momentum?

---

*End of brainstorming notes*
