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

### 4.1 Physical Meaning

The displacement height $d$ is the **effective level at which the canopy absorbs momentum**.
For a forest or urban canopy of geometric height $h_c$, common empirical estimates are:

$$
d \approx \tfrac{2}{3}\,h_c \quad \text{(vegetated canopies; Jackson 1981)}
\qquad
d \approx 0.5\text{–}0.8\,h_c \quad \text{(urban, dense canopies)}
$$

Below $z = d$, the log profile has no physical meaning — the wind is sheltered within the
canopy and turbulence structure is entirely different.  Above $z = d + z_0$, the adjusted
log profile applies:

$$
u(z) = \frac{u_*}{\kappa} \ln\!\frac{z - d}{z_0}
$$

Here $z_0$ is the **aerodynamic roughness length measured above the displacement level** —
typically $z_0 \approx h_c/10$ for vegetation.  A forest with $h_c = 15\,\text{m}$ might have
$d \approx 10\,\text{m}$ and $z_0 \approx 1.5\,\text{m}$.

### 4.2 Log Expansion with Displacement

Write $\ln((z-d)/z_0) = \ln(z/z_0) + \ln(1 - d/z)$.  Using $\ln(1+x) = x - x^2/2 + \cdots$ with
$x = -d/z$:

$$
\ln\!\frac{z-d}{z_0} = \ln\!\frac{z}{z_0} - \frac{d}{z} - \frac{d^2}{2z^2} - \frac{d^3}{3z^3} - \cdots
$$

**All correction terms are negative** (because $x = -d/z < 0$ and odd powers of a negative
number are negative, while the alternating sign of the series then makes all terms add
negatively).  This is opposite to Section 1 where the correction was positive.

**Convergence:** requires $d/z < 1$, i.e., the profile is only well-defined above $z = d$.
The series converges slowly near $z \sim d$ and diverges below.

**Danger near $z = d$:** The log diverges to $-\infty$ as $z \to d^+$, which is unphysical.
This is why MOST is applied only above the roughness sublayer; in practice the minimum valid
height is $z \gtrsim d + 2z_0$ to $d + 5z_0$.

### 4.3 Shear Error from Ignoring $d$

The correct local shear with displacement is:

$$
\frac{\partial u}{\partial z} = \frac{u_*}{\kappa(z - d)}
$$

If $d$ is ignored (set to zero):

$$
\frac{\partial u}{\partial z}\Big|_{\text{no-}d} = \frac{u_*}{\kappa z}
$$

The ratio of true to naive shear is:

$$
\frac{(\partial u/\partial z)_{\text{true}}}{(\partial u/\partial z)_{\text{naive}}} = \frac{z}{z - d}
= \frac{1}{1 - d/z} \approx 1 + \frac{d}{z} + \left(\frac{d}{z}\right)^2 + \cdots
$$

This is now $> 1$: **ignoring $d$ underestimates the shear** (the denominator $z - d < z$
means the true gradient is steeper than $u_*/\kappa z$).

Concrete example: forest canopy with $h_c = 15\,\text{m}$, $d = 10\,\text{m}$, $z_0 = 1.5\,\text{m}$:

| $z$ (m) | $z - d$ (m) | Error in shear: $z/(z-d) - 1$ |
|---------|------------|-------------------------------|
| 12      | 2          | 500% (catastrophic) |
| 15      | 5          | 200% |
| 20      | 10         | 100% |
| 30      | 20         | 50% |
| 50      | 40         | 25% |
| 100     | 90         | 11% |

Even at $z = 100\,\text{m}$ (well above the canopy), the error is 11%.  For urban canyons with
$d \sim 20$ m, an NWP first model level at $z = 30\,\text{m}$ has $z - d = 10\,\text{m}$ — a
200% shear error.

### 4.4 Effect on the Richardson Number

With displacement, the gradient Richardson number is:

$$
Ri_g(z) = \frac{(g/\overline{\Theta})\,\partial\Theta/\partial z}{(\partial u/\partial z)^2}
= \frac{(g/\overline{\Theta})\,\theta_*/\kappa(z-d)}{[u_*/\kappa(z-d)]^2}
= \frac{g\kappa\theta_*}{u_*^2\,\overline{\Theta}}
$$

Crucially, $Ri_g$ is **still height-independent** when displacement is included consistently
(the $(z-d)$ factors cancel just as the $z$ factors did in the no-displacement case).  The
danger arises when the **profile code uses $z$ but the flux code uses $z - d$**, or vice versa,
producing a spurious height dependence in $Ri_g$ and potentially a spurious regime transition.

### 4.5 Bulk Richardson Number with Displacement

Over a layer $[z_1, z_2]$, the bulk $Ri_B$ becomes:

$$
Ri_B = \frac{(g/\overline{\Theta})\,(\Theta_2 - \Theta_1)}{[u(z_2) - u(z_1)]^2/(z_2 - z_1)}
= \frac{g\kappa\theta_*\,(z_2 - z_1)}{u_*^2\,\overline{\Theta}\,\ln\!\frac{z_2-d}{z_1-d}}
$$

Ignoring $d$ gives $\ln(z_2/z_1)$ in the denominator.  The ratio:

$$
\frac{\ln(z_2/z_1)}{\ln\!\frac{z_2-d}{z_1-d}}
$$

is $> 1$ because $z/(z-d)$ decreases with height, so $\ln(z_2/z_1) < \ln((z_2-d)/(z_1-d))$.
Ignoring $d$ **over-estimates** $Ri_B$ — it biases toward stability.

Example: $z_1 = 2\,\text{m}$, $z_2 = 10\,\text{m}$, $d = 0$ vs.\ $d = 1\,\text{m}$
(low grassland with $h_c \approx 1.5\,\text{m}$):

$$
\ln(10/2) = 1.609, \qquad \ln(9/1) = 2.197
$$

Ratio = 1.37 — so $Ri_B$ is **overestimated by 37%** even for a modest 1 m displacement.

### 4.6 Three Scales and Their Ordering

In practice, there are three distinct length scales:

$$
z_0 \quad < \quad d \quad < \quad h_c
$$

with rough empirical ratios $z_0 \approx h_c/10$, $d \approx 2h_c/3$.  The log profile requires
$z \gg d \gg z_0$, i.e. the measurement must be well above the canopy displacement level.
The quantity $\ln(z/z_0)$ or $\ln((z-d)/z_0)$ is only large (making the log approximation
good) when this three-level ordering is satisfied by a comfortable margin.

For short grass ($h_c = 0.1\,\text{m}$, $d \approx 0.07\,\text{m}$, $z_0 \approx 0.01\,\text{m}$):
measurements at $z = 2\,\text{m}$ are easily in the valid regime because $z/d \approx 29$.

For urban terrain ($h_c = 15\,\text{m}$, $d \approx 10\,\text{m}$, $z_0 \approx 1.5\,\text{m}$):
measurements at $z = 20\,\text{m}$ give $z - d = 10\,\text{m}$ and $(z-d)/z_0 \approx 7$,
which is marginal.  Most standard-height meteorological masts are **inside or just above the
roughness sublayer** for dense urban or forest surfaces.

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

## 9. Pitfalls of Log-Profile Calculations Across ABL Regimes

Each stability regime introduces its own failure modes.  Many are directly traceable to the
log-expansion breakdown discussed above.

---

### 9.1 Neutral Regime

**Pitfall 1 — Valid height range violated.**
The constant-flux surface layer is typically the lowest 10–15% of the ABL.  In a well-mixed
daytime ABL of $h = 1\,\text{km}$, the surface layer extends to ~100–150 m.  In a shallow
nocturnal SBL of $h = 50\,\text{m}$, the surface layer ends at just ~5–7 m.  Applying the
neutral log profile above the surface layer imports an assumption that the flux is constant
that may be wrong by 50% or more at the top of the SBL.

**Pitfall 2 — Horizontal inhomogeneity.**
MOST requires horizontally homogeneous, stationary flow.  Internal boundary layers from
land-sea boundaries, urban edges, or orography can cause the flux to vary horizontally.
The neutral log profile then holds only in the equilibrium layer close to the surface
($z \lesssim 0.1 x$, where $x$ is the fetch distance from the surface change).

**Pitfall 3 — $z_0$ and $d$ retrieval is ill-conditioned.**
Extracting $z_0$ and $d$ from wind profiles requires fitting $\ln(z-d)$ to $u$, varying
three parameters ($u_*$, $z_0$, $d$).  The problem is poorly conditioned: small measurement
errors propagate to large errors in $d$, which then shift $z_0$ by an order of magnitude.
For a single pair of observations $(u_1, u_2)$ at $(z_1, z_2)$, only $u_*$ and the
combination $\ln((z-d)/z_0)$ are determined — there are infinitely many $(d, z_0)$ pairs
consistent with the data.

**Pitfall 4 — Roughness sublayer contamination.**
For $z < 2$–$5 h_c$, turbulence is directly influenced by individual roughness elements; $\phi_m \neq 1$
even under neutral conditions.  The log profile underestimates turbulent exchange here.
In practice this means the neutral log profile is often applied 5–15 m lower than it should be,
particularly in NWP models that do not resolve the roughness sublayer.

---

### 9.2 Unstable (Convective) Regime

**Pitfall 5 — Iterative inversion needed; poor convergence.**
Given observed wind shear and temperature gradient, finding $u_*$ and $\zeta$ requires solving
the implicit system:

$$
u_* = \frac{\kappa(U_2 - U_1)}{\ln(z_2/z_1) - \psi_m(\zeta)}, \qquad
L = -\frac{u_*^3 \overline{\Theta}}{\kappa g H_s}
$$

This must be iterated.  Starting from the neutral guess ($\psi_m = 0$) and iterating under
unstable conditions, convergence is typically fast.  But near the free-convection limit
($u_* \to 0$, $L \to 0^-$, $|\zeta| \to \infty$), $L$ becomes very sensitive to small
changes in $u_*$ and the iteration can oscillate.

**Pitfall 6 — Counter-gradient heat transport not captured.**
Deep in the convective mixed layer, large-scale updrafts carry warm air upward even where the
local temperature gradient is slightly positive.  K-theory (local gradient flux with $K_H > 0$)
gives the wrong sign of $H_s$.  The non-local transport by convective plumes is entirely missed
by MOST, which is a local surface-layer theory.  This is why bulk-transfer schemes must use
a separate entrainment parameterization for the mixed-layer top.

**Pitfall 7 — $\ln(z/z_0)$ vs.\ $[\ln(z/z_0) - \psi_m]$ confusion.**
Under strong instability, $\psi_m$ can be as large as 4–6 for $\zeta = -2$ to $-5$.
$\ln(z/z_0)$ might be $\sim 4$–$6$ for $z = 10\,\text{m}$, $z_0 = 0.05\,\text{m}$.
So **the correction term is as large as the log term itself** — entirely non-perturbative.
Using the neutral log profile in the unstable regime underestimates $u_*$ and overestimates $C_D$.

**Pitfall 8 — The free-convection singularity.**
As $|L| \to 0$ (e.g., calm winds with strong heating), $\zeta \to -\infty$ and the Businger–Dyer
functions go as $\phi_m \to (-16\zeta)^{-1/4} \to 0$.  Wind-based scaling breaks down entirely.
The correct scaling in this limit uses $w_*$ rather than $u_*$, and profile-based flux estimates
become unreliable when $z/|L| \gtrsim 2$.

---

### 9.3 Near-Neutral Regime

**Pitfall 9 — Ambiguity in the sign of stability.**
Near neutral, the temperature difference $\Delta\Theta = \Theta(z_2) - \Theta(z_1)$
may be within the measurement noise.  A thermocouple or aspirated psychrometer with
$\pm 0.1\,\text{K}$ resolution over a $\Delta z = 8\,\text{m}$ layer gives a gradient uncertainty
of $\pm 0.012\,\text{K m}^{-1}$ — which for $u_* = 0.4\,\text{m s}^{-1}$ corresponds to an
Obukhov length uncertainty of hundreds of meters.  The system could flip between stably and
unstably classified, changing $\psi_m$ from positive to negative.  In an NWP model running
with timesteps of minutes, this produces oscillations in $C_D$ and $H_s$.

**Pitfall 10 — $\psi_m$ corrections first-order but not negligible.**
For $|\zeta| = 0.1$ the correction to $C_D$ is $\approx 2\beta\zeta/\ln(z/z_0)$.  Over smooth
terrain ($z_0 = 0.01\,\text{m}$, $z = 10\,\text{m}$) $\ln(z/z_0) \approx 6.9$, giving a 14%
change in $C_D$.  Over rough terrain ($z_0 = 0.5\,\text{m}$) $\ln(z/z_0) \approx 3.0$, giving a
33% change.  Near-neutral should **not** be treated as exactly neutral for flux computation
over rough surfaces.

---

### 9.4 Weakly-to-Strongly Stable Regime

**Pitfall 11 — MOST breaks down for $\zeta \gtrsim 1$.**
Many observational campaigns (Kansas, Cabauw, CASES-99) show that $\phi_m$, $\phi_h$ depart
significantly from the Businger–Dyer linear forms for $\zeta > 0.5$–$1$.  Alternate forms
(Webb 1970, Beljaars–Holtslag 1991) attempt corrections, but no universal form is established.
Fluxes computed with the linear form at $\zeta = 2$ may be off by a factor of 2.

**Pitfall 12 — Log profile no longer describes the actual wind profile.**
In the strongly stable case the profile above the surface layer can be super-logarithmic,
sub-logarithmic, or even show a low-level jet maximum.  The surface layer itself may be only
2–5 m deep.  Applying the log profile from $z_0$ to $z = 10\,\text{m}$ then spans both the
true log layer and the stably stratified layer above it — the inferred $u_*$ is bogus.

**Pitfall 13 — Non-stationarity from internal gravity waves.**
Gravity waves generated by the stable density stratification cause periodic wind and temperature
oscillations at 10–30 minute periods.  Eddy-covariance flux measurements require 20–30 minute
averaging but the signals are non-stationary over these intervals — ogive curves fail to plateau.
The computed $u_*$ and $\theta_*$ carry large random errors ($\pm 50\%$ is common).

**Pitfall 14 — Iteration blows up near $Ri_c$.**
The formula $\zeta = Ri_g/(1 - \beta Ri_g)$ diverges as $Ri_g \to Ri_c = 1/\beta$.  In a model
that diagnoses $Ri_g$ from finite-difference shear (already biased toward low shear by
$-\Delta z/(2z_k)$), the spuriously reduced shear inflates $Ri_g$, which is then
magnified further by the $1/(1-\beta Ri_g)$ amplification near the singularity — a
compounding numerical disaster.  The grid-curvature correction and the MOST inversion are
coupled failure modes near $Ri_c$.

---

### 9.5 Laminar (Post-Collapse) Regime

**Pitfall 15 — Setting $K_M = 0$ exactly causes runaway cooling.**
If the model sets eddy diffusivity to zero (rather than to the molecular value $\nu$) when
$Ri \geq Ri_c$, the only heat exchange is by radiation, and the surface temperature can drop
unboundedly in a model timestep.  Real nocturnal minimum temperatures are limited by soil
heat conduction and intermittent turbulence; the pure laminar model has no such floor.

**Pitfall 16 — Intermittency is not captured by a mean parameterization.**
The real stable BL alternates between weakly turbulent events and nearly laminar quiescent
periods.  The ensemble mean flux is some fraction (perhaps 20–40%) of the fully turbulent
value.  A model that uses the time-mean $Ri$ to look up a smooth $f_m$ curve sees a single
point on the $Ri$ curve and misses the intermittent nature.  Two models with identical
time-mean forcings but different subgrid variability will produce different surface temperatures.

**Pitfall 17 — The flip-flop transition creates hysteresis; initial conditions matter.**
As documented by McNider et al., the stable BL can reside in either the turbulent or laminar
equilibrium under the same external forcing (radiation, wind).  Starting from different initial
conditions (e.g., evening transition vs.\ morning recovery) gives different nocturnal
temperatures.  This is not merely a parameterization deficiency — it is a real physical
bi-stability — but the model's ability to track the correct branch depends on having correctly
handled all the pitfalls above.

---

### 9.6 Summary Table of Pitfalls by Regime

| # | Regime | Pitfall | Root Cause |
|---|--------|---------|------------|
| 1 | Neutral | Surface layer top exceeded | Flux non-constant above $\sim 0.1h$ |
| 2 | Neutral | Inhomogeneity / limited fetch | MOST requires horizontal homogeneity |
| 3 | Neutral | Ill-conditioned $z_0$, $d$ retrieval | Log-profile fit underdetermined |
| 4 | Neutral | Roughness sublayer contamination | $\phi_m \neq 1$ for $z < 5h_c$ |
| 5 | Unstable | Iterative inversion diverges near $u_* \to 0$ | Free-convection singularity |
| 6 | Unstable | Counter-gradient fluxes missed | Non-local mixed-layer transport |
| 7 | Unstable | $\psi_m$ as large as $\ln(z/z_0)$ | Full correction, not perturbative |
| 8 | Unstable | Free-convection singularity | $|L| \to 0$, wind scaling fails |
| 9 | Near-neutral | Sign-of-stability ambiguity | Measurement noise in $\Delta\Theta$ |
| 10 | Near-neutral | $\psi_m$ non-negligible over rough terrain | Small $\ln(z/z_0)$ amplifies effect |
| 11 | Stable | Linear $\phi_m$ form incorrect for $\zeta > 1$ | No universal SBL functions |
| 12 | Stable | Profile not logarithmic in SBL | LLJ, super-log structure |
| 13 | Stable | Gravity wave non-stationarity | Periodic oscillations corrupt fluxes |
| 14 | Stable | Iteration amplifies shear bias near $Ri_c$ | Grid curvature × MOST singularity |
| 15 | Laminar | Runaway cooling if $K_M = 0$ exactly | No molecular diffusivity floor |
| 16 | Laminar | Intermittency not captured by mean $f_m$ | Subgrid variability in $Ri$ |
| 17 | Laminar | Path-dependence / bi-stability | Multiple equilibria (McNider et al.) |

---

*End of brainstorming notes*
