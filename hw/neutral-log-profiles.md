# Homework: Surface-Layer Profiles from Neutral to Laminar — MOST, Stability Functions, and Regime Transitions

**Course:** Atmospheric Boundary Layer Physics
**Topic:** Logarithmic and stability-corrected profiles; neutral, unstable, near-neutral, stable, and laminar regimes
**Prerequisites:** Turbulence fundamentals, MOST, dimensional analysis, K-theory

---

## Background

Under **neutral atmospheric conditions** ($Ri_g \approx 0$, buoyancy forces negligible relative to
shear production), Monin–Obukhov Similarity Theory (MOST) reduces to the classical
**logarithmic profiles** for mean wind speed and potential temperature:

$$
U(z) = \frac{u_*}{\kappa} \ln\!\left(\frac{z}{z_0}\right)
$$

$$
\Theta(z) - \Theta_s = \frac{\theta_*}{\kappa} \ln\!\left(\frac{z}{z_{0h}}\right)
$$

where $u_*$ is the friction velocity $[\text{m s}^{-1}]$, $\theta_*$ is the temperature scale
$[\text{K}]$, $\kappa \approx 0.40$ is the von Kármán constant, $z_0$ is the aerodynamic
roughness length for momentum, $z_{0h}$ is the scalar roughness length for heat, and $\Theta_s$
is the surface temperature. These profiles are valid in the **constant-flux surface layer**,
assumed to be the lowest 10–15% of the ABL.

Turbulent fluxes are defined as:

$$
\tau \equiv -\overline{u'w'} = u_*^2 \qquad \text{(kinematic momentum flux, m}^2\text{ s}^{-2}\text{)}
$$

$$
H_s \equiv -\overline{w'\theta'} = -u_* \theta_* \qquad \text{(kinematic sensible heat flux, K m s}^{-1}\text{)}
$$

**Regime Overview.** The neutral log profile is the baseline case ($\zeta = 0$, $|L| \to \infty$). As stratification develops, MOST extends the profiles through stability correction functions $\psi_m(\zeta)$ and $\psi_h(\zeta)$, where $\zeta \equiv z/L$ is the dimensionless stability parameter and $L$ is the Obukhov length. Four regimes span this problem set:

| Regime | $\zeta$ | Physics |
|---|---|---|
| Unstable (UBL) | $\zeta < 0$ | Buoyancy enhances mixing; $\phi_m < 1$; enhanced fluxes |
| Neutral | $\zeta \approx 0$ | Pure log profile; $\phi_m = 1$ |
| Weakly stable / near-neutral | $0 < \zeta \lesssim 0.1$ | Log-linear profile; $\phi_m > 1$; reduced fluxes |
| Strongly stable / laminar | $\zeta \gg 1$ or $Ri \geq Ri_c$ | Turbulence collapse; $f_m \to 0$; molecular diffusion |

Problems 1–5 develop the neutral baseline. Problems 6–8 extend to the full stability spectrum.

---

## Problem 1 — Vertical Wind Shear from the Log Profile (20 points)

### Part A (5 points)

Starting from $U(z) = \dfrac{u_*}{\kappa}\ln(z/z_0)$, analytically derive the vertical wind
shear:

$$
\frac{\partial U}{\partial z} = \; ?
$$

Show your steps. What is the functional dependence on $z$? Is the shear largest near the
surface or aloft?

### Part B (5 points)

The **non-dimensional shear** (MOST momentum function) is defined as:

$$
\phi_m(\zeta) \equiv \frac{\kappa z}{u_*}\frac{\partial U}{\partial z}
$$

Under neutral conditions $\zeta = z/L \to 0$.  Using your result from Part A, evaluate
$\phi_m(0)$ analytically and confirm it equals 1. Explain why $\phi_m(0) = 1$ is the
*definition* of the neutral limit.

### Part C (10 points)

Derive the **bulk wind shear** between two levels $z_1 < z_2$ in the surface layer:

$$
\Delta U = U(z_2) - U(z_1) = \frac{u_*}{\kappa}\ln\!\frac{z_2}{z_1}
$$

1. What happens to $\Delta U$ as $z_1 \to z_0^+$? Interpret physically.
2. A sonic anemometer gives $u_* = 0.35\,\text{m s}^{-1}$ and $z_0 = 0.05\,\text{m}$.
   Compute $\Delta U$ between $z_1 = 2\,\text{m}$ and $z_2 = 10\,\text{m}$.
3. Compute the local shear $\partial U/\partial z$ at those two heights and note how rapidly it
   changes.

---

## Problem 2 — Vertical Temperature Gradient and Sensible Heat Flux (20 points)

### Part A (5 points)

Starting from $\Theta(z) = \Theta_s + \dfrac{\theta_*}{\kappa}\ln(z/z_{0h})$, derive the
vertical temperature gradient $\partial \Theta / \partial z$.

### Part B (5 points)

The turbulent sensible heat flux is often expressed via K-theory:

$$
H_s = -K_H \frac{\partial \Theta}{\partial z}
$$

Under neutral MOST, $K_H = K_M = \kappa u_* z$.  Using your result from Part A and the
definition $H_s = -u_*\theta_*$, verify that K-theory is *self-consistent* with the log
profile. That is, show:

$$
-K_H \frac{\partial \Theta}{\partial z} = -u_* \theta_*
$$

### Part C (5 points)

**Lapse rate at the first model level.**  Typical NWP models place their lowest thermodynamic
level at $z_k = 10\,\text{m}$.  For a typical nocturnal case: $u_* = 0.20\,\text{m s}^{-1}$,
$\theta_* = +0.10\,\text{K}$ (stable, $\theta_*>0$ convention means downward heat flux).

1. Compute $\partial \Theta / \partial z$ at $z = 10\,\text{m}$.
2. Convert to a lapse rate in $\text{K}/100\,\text{m}$. How does it compare to the dry
   adiabatic lapse rate ($9.8\,\text{K}/\text{km}$)?

### Part D (5 points)

Derive a formula for the temperature difference between $z_1$ and $z_2$ in terms of $\theta_*$
and the roughness lengths.  Under what conditions is the temperature profile *independent* of
$z_{0h}$?

---

## Problem 3 — Richardson Number under Neutral Log Profiles (25 points)

### Part A — Gradient $Ri_g$ (10 points)

The gradient Richardson number is:

$$
Ri_g(z) = \frac{(g/\overline{\Theta}_v)(\partial \Theta_v/\partial z)}{(\partial U/\partial z)^2}
$$

Assume $\Theta_v \approx \Theta$ (dry air) and use your results from Problems 1 and 2.
Substituting the log-profile gradients:

1. Write $Ri_g(z)$ fully in terms of $u_*$, $\theta_*$, $\kappa$, $g$, $\overline{\Theta}$,
   and $z$.
2. Notice that the $z$-dependence cancels.  Show explicitly that $Ri_g$ is **height-independent
   throughout the log layer** for a given $(u_*, \theta_*)$.
3. Therefore $Ri_g$ is a property of the *flow state*, not the measurement height — why is this
   either expected or surprising from a MOST standpoint?

### Part B — Neutral limit (5 points)

Under strictly neutral conditions $\theta_* \to 0$.  What is $Ri_g$?  Using the MOST relation
$Ri_g = \zeta\phi_h/\phi_m^2$ (with $\phi_m = \phi_h = 1$ at neutral), verify this is consistent.

### Part C — Bulk Richardson Number (10 points)

The **bulk Richardson number** over a layer $[z_1, z_2]$ is:

$$
Ri_B = \frac{(g/\overline{\Theta})(\Theta(z_2) - \Theta(z_1))}{[U(z_2) - U(z_1)]^2 / (z_2 - z_1)}
$$

1. Substitute the neutral log profiles and derive $Ri_B$ in closed form.  Show:

$$
Ri_B = \frac{(g/\overline{\Theta})\,\theta_*\,\kappa^{-1} \ln(z_2/z_1)}{u_*^2\,\kappa^{-2}\,\ln^2(z_2/z_1)/(z_2-z_1)}
= \frac{g\,\kappa\,\theta_*\,(z_2-z_1)}{u_*^2\,\overline{\Theta}\,\ln(z_2/z_1)}
$$

2. Define $\overline{\Theta} = 290\,\text{K}$, $u_* = 0.30\,\text{m s}^{-1}$,
   $\theta_* = +0.05\,\text{K}$, $z_1 = 2\,\text{m}$, $z_2 = 10\,\text{m}$.
   Calculate $Ri_B$ numerically.  Is the layer stable, neutral, or unstable?

3. **Key difference:** $Ri_g$ is height-independent but $Ri_B$ depends on $z_1$, $z_2$ through
   the factor $(z_2-z_1)/\ln(z_2/z_1)$.  Show that as $z_2 \to z_1$ (thin layer),
   $Ri_B \to Ri_g$. *(Hint: use L'Hôpital's rule or the limit $\lim_{\epsilon\to0}\epsilon/\ln(1+\epsilon/z) = z$.)*

---

## Problem 4 — Momentum and Heat Fluxes: From Profiles Back to Fluxes (20 points)

### Part A — Inverse problem (10 points)

In practice, we *measure* the wind and temperature at two levels and must *infer* the fluxes.
Given only observations $(U_1, U_2)$ at $(z_1, z_2)$ and $(\Theta_1, \Theta_2)$ at the same
heights (assuming $z_{0h} = z_0$):

1. Derive formulas for $u_*$ and $\theta_*$ in terms of observables:

$$
u_* = \frac{\kappa(U_2 - U_1)}{\ln(z_2/z_1)}, \qquad
\theta_* = \frac{\kappa(\Theta_2 - \Theta_1)}{\ln(z_2/z_1)}
$$

2. Hence write the sensible heat flux $H_s = -u_*\theta_*$ as a function of observable bulk
   gradients only.  What class of parameterization is this?

3. A weather station gives $U_1 = 3.0\,\text{m s}^{-1}$ at 2 m, $U_2 = 5.2\,\text{m s}^{-1}$
   at 10 m, $\Theta_1 = 288.5\,\text{K}$ at 2 m, $\Theta_2 = 289.1\,\text{K}$ at 10 m.
   Calculate $u_*$, $\theta_*$, $H_s$, and $\tau$.

### Part B — Drag coefficient formulation (5 points)

Define a **bulk drag coefficient** $C_D$ such that $\tau = C_D U_r^2$, where $U_r = U(z_r)$ is
the wind at a reference height $z_r$.  Show:

$$
C_D(z_r) = \frac{\kappa^2}{\ln^2(z_r/z_0)}
$$

Evaluate $C_D$ for $z_r = 10\,\text{m}$ and $z_0 = 0.01\,\text{m}$ (smooth grassland).
Compare to the value for $z_0 = 0.5\,\text{m}$ (suburban canopy).  What does the ratio tell
you about surface drag contrast between these two surfaces?

### Part C — Heat transfer coefficient (5 points)

Analogously define a bulk heat transfer coefficient $C_H$ such that
$H_s = -C_H U_r (\Theta_r - \Theta_s)$.

1. Derive $C_H$ assuming $z_{0h} = z_0$.
2. The ratio $C_H/C_D$ is the **Stanton number** for neutral conditions.  Evaluate it; is it
   identically 1, and why (or why not)?
3. Now assume $z_{0h} = z_0/7$ (Reynolds-number-based Brutsaert formula).  How does $C_H$ change?

---

## Problem 5 — Spin-off: Log Expansion and Grid Discretization Error (15 points)

*This problem connects the log profile to the brainstorming memo on log expansions and to the
grid-curvature correction literature.*

### Part A (5 points)

Consider a model that estimates the local shear between adjacent levels $z_k$ and $z_{k+1} = z_k + \Delta z$ using a centered finite difference:

$$
\left.\frac{\partial U}{\partial z}\right|_{\text{FD}} = \frac{U(z_k + \Delta z) - U(z_k)}{\Delta z}
= \frac{u_*}{\kappa \Delta z}\ln\!\left(1 + \frac{\Delta z}{z_k}\right)
$$

Expand $\ln(1+\Delta z / z_k)$ in a Taylor series.  Show that the first-order bias in the finite
difference vs. the true derivative at $z_k$ is:

$$
\left.\frac{\partial U}{\partial z}\right|_{\text{FD}} \approx \frac{u_*}{\kappa z_k}\left(1 - \frac{\Delta z}{2z_k} + \frac{(\Delta z)^2}{3z_k^2} - \cdots\right)
$$

The leading error term is $-\Delta z / (2z_k)$.  Comment on the sign: does the FD method over- or under-estimate the shear?

### Part B (5 points)

Evaluate the relative error $|\text{bias}| / |\text{true shear}| \approx \Delta z / (2z_k)$ for:

| $z_k$ (m) | $\Delta z$ (m) | Relative shear error |
|-----------|---------------|----------------------|
| 5         | 10            |                      |
| 10        | 10            |                      |
| 20        | 10            |                      |
| 50        | 10            |                      |
| 100       | 10            |                      |

Fill in the table and sketch (roughly) how this error propagates into $Ri_g^2$ in the denominator
of $Ri_g$ (i.e., shear appears *squared*, so the Richardson number bias is approximately twice
the shear relative error).

### Part C (5 points)

Suppose you want the relative shear error to be $< 5\%$.  What constraint does this place on
the ratio $\Delta z / z_k$?  For a first model level at $z_k = 5\,\text{m}$, what maximum grid
spacing $\Delta z$ is required?  Is this feasible in current operational NWP models (typical
$\Delta z_1 \approx 10$–$30\,\text{m}$)?

---

---

## Problem 6 — Unstable Surface Layer: MOST $\psi$-Functions and Modified Log Profiles (25 points)

*Under unstable conditions ($\zeta < 0$, $L < 0$), buoyant convection enhances turbulent mixing.
MOST extends the log profile by integrating the gradient functions into correction functions
$\psi_m$ and $\psi_h$.*

The general stability-corrected profiles are:

$$
U(z) = \frac{u_*}{\kappa}\left[\ln\frac{z}{z_0} - \psi_m(\zeta)\right], \qquad
\Theta(z) - \Theta_s = \frac{\theta_*}{\kappa}\left[\ln\frac{z}{z_{0h}} - \psi_h(\zeta)\right]
$$

where:

$$
\psi_m(\zeta) = \int_0^{\zeta} \frac{1 - \phi_m(\zeta')}{\zeta'}\,d\zeta', \qquad
\psi_h(\zeta) = \int_0^{\zeta} \frac{1 - \phi_h(\zeta')}{\zeta'}\,d\zeta'
$$

For the **Businger–Dyer unstable forms** ($\zeta < 0$):

$$
\phi_m(\zeta) = (1 - 16\zeta)^{-1/4}, \qquad \phi_h(\zeta) = (1 - 16\zeta)^{-1/2}
$$

### Part A — Deriving $\psi_m$ and $\psi_h$ (8 points)

Let $x \equiv (1 - 16\zeta)^{1/4}$ so that $\phi_m = x^{-1}$, $\phi_h = x^{-2}$.

1. Change variable from $\zeta$ to $x$ in the $\psi_m$ integral and apply partial fraction
   decomposition to show the result is the **Paulson (1970)** form:
   $$
   \psi_m(\zeta) = 2\ln\frac{1+x}{2} + \ln\frac{1+x^2}{2} - 2\arctan(x) + \frac{\pi}{2}
   $$
   Verify $\psi_m(0) = 0$ (at $\zeta = 0$: $x = 1$).

2. Show that the heat correction evaluates to:
   $$
   \psi_h(\zeta) = 2\ln\frac{1+x^2}{2}
   $$
   Verify $\psi_h(0) = 0$.

### Part B — Numerical Evaluation at $\zeta = -1$ (7 points)

For $\zeta = -1$, $z_0 = 0.05\,\text{m}$, $u_* = 0.30\,\text{m s}^{-1}$, $z = 10\,\text{m}$:

1. Compute $x = 17^{1/4}$.  Evaluate $\phi_m$, $\phi_h$, $\psi_m$, $\psi_h$ numerically.
2. Compute $U(10\,\text{m})$ with and without the $\psi_m$ correction.  By what percentage does
   the unstable correction reduce the mean wind for the same $u_*$?  Is this physically
   expected — explain in terms of mixing efficiency.
3. Compute the local shear $\partial U/\partial z$ at $z = 10\,\text{m}$ via $\phi_m$ and
   compare to the neutral value.  Which regime has more wind shear near the surface?

### Part C — Richardson Number in the Unstable Case (5 points)

Evaluate $Ri_g = \zeta\,\phi_h(\zeta)/\phi_m(\zeta)^2$ for the Businger–Dyer forms.

1. Show that $\phi_h/\phi_m^2 = 1$ identically, so $Ri_g = \zeta$ for all $\zeta < 0$.
2. Is this a general MOST result or specific to the Businger–Dyer exponent choices?  Contrast
   with the stable linear forms you will derive in Problem 8A.
3. At $\zeta = -1$: what is $Ri_g$?  Is the flow turbulent by the Miles–Howard criterion
   ($Ri_g < 1/4$)?

### Part D — Breakdown of Wind-Based Scaling (5 points)

The **convective velocity scale** $w_* \equiv (g H_s h / \overline{\Theta})^{1/3}$ characterizes
mixed-layer convection.

1. Compute $w_*/u_*$ for: kinematic $H_s = 0.17\,\text{K m s}^{-1}$
   ($\approx 200\,\text{W m}^{-2}/\rho c_p$), $u_* = 0.30\,\text{m s}^{-1}$,
   $h = 1000\,\text{m}$, $\overline{\Theta} = 295\,\text{K}$.
2. When $w_* \gg u_*$, convective plumes dominate and the log profile scaling breaks down.
   Comment on whether log profiles are appropriate throughout a deep daytime convective
   boundary layer.

---

## Problem 7 — Near-Neutral Regime: Log-Linear Profiles and Regime Boundaries (20 points)

*Near neutral ($|\zeta| \ll 1$), the correction functions linearize and the profile transitions
from pure log to log-linear.  This regime governs the evening transition and early-morning
spin-up in NWP models.*

### Part A — The Log-Linear Profile (Weakly Stable) (8 points)

For weakly stable conditions use $\phi_m = \phi_h = 1 + \beta\zeta$ with $\beta = 5$.

1. Integrate $dU/dz = (u_*/\kappa z)(1 + \beta z/L)$ from $z_0$ to $z$ to show:
   $$
   U(z) = \frac{u_*}{\kappa}\left[\ln\frac{z}{z_0} + \beta\frac{z - z_0}{L}\right]
          \approx \frac{u_*}{\kappa}\left[\ln\frac{z}{z_0} + \beta\zeta\right]
   $$
   where the approximation drops the small $\beta z_0/L$ term.

2. Verify from the $\psi_m$ definition that $\psi_m(\zeta) = -\beta\zeta$ for the linear stable
   form, recovering the profile from the general expression $U = (u_*/\kappa)[\ln(z/z_0) - \psi_m]$.

3. At $z = 10\,\text{m}$, $z_0 = 0.10\,\text{m}$: at what $\zeta$ does the correction term
   $\beta\zeta$ equal 10% of the log term?  What is the corresponding Obukhov length $L$?

4. Compare the first-order stable correction to the first-order unstable correction from the
   linearized Businger–Dyer ($\phi_m \approx 1 + 4\zeta$ near $\zeta = 0^-$, noting $\zeta < 0$
   there).  Are the two-sided corrections symmetric about neutral?  Why might they differ
   physically?

### Part B — Obukhov Length and Regime Boundaries (5 points)

1. For the criterion $|\zeta| < 0.05$ at $z = 10\,\text{m}$, what is the minimum $|L|$?
2. Compute $L$ for $u_* = 0.40\,\text{m s}^{-1}$, $\theta_* = 0.02\,\text{K}$,
   $\overline{\Theta} = 290\,\text{K}$.  Is this near-neutral at 10 m?  At 50 m?
3. Sketch the regime boundaries in the $(u_*, \theta_*)$ plane and label the four regimes from
   the overview table in the Background.

### Part C — Bulk Drag Coefficient Correction (7 points)

With the stability-corrected profile the drag coefficient becomes:

$$
C_D(\zeta) = \frac{\kappa^2}{[\ln(z_r/z_0) - \psi_m(\zeta)]^2}
$$

1. For the stable near-neutral case ($\psi_m = -\beta\zeta, \zeta > 0$), show $C_D$
   **decreases** relative to neutral.  Explain the sign physically.

2. Expand $C_D(\zeta)$ to first order in $\beta\zeta/\ln(z_r/z_0) \ll 1$ to obtain:
   $$
   C_D(\zeta) \approx C_D^{\text{neutral}}\!\left(1 - \frac{2\beta\zeta}{\ln(z_r/z_0)}\right)
   $$

3. For $z_r = 10\,\text{m}$, $\zeta = +0.2$: evaluate the percentage reduction in $C_D$ for
   (a) $z_0 = 0.01\,\text{m}$ (smooth grassland) and (b) $z_0 = 0.5\,\text{m}$ (suburban
   canopy).  Does a rough or smooth surface show a larger **relative** correction, and why?

---

## Problem 8 — Stable BL: Critical Richardson Number, Turbulence Collapse, and the Laminar Regime (30 points)

*Under clear skies with light winds — the canonical nocturnal stable boundary layer — surface
radiative cooling drives an increasing Richardson number.  When $Ri_g$ exceeds the critical
value $Ri_c$, the MOST-based turbulence parameterization collapses, leaving only molecular
diffusion.  This regime transition is central to nocturnal surface temperature forecasting and
pollutant trapping.*

### Part A — Critical Richardson Number from MOST (10 points)

Given $\phi_m = \phi_h = 1 + \beta\zeta$ ($\beta = 5$, turbulent Prandtl number $= 1$):

1. Show:
   $$
   Ri_g(\zeta) = \frac{\zeta}{1 + \beta\zeta}
   $$

2. Show $Ri_g$ is monotonically increasing but **bounded above**:
   $$
   \lim_{\zeta \to \infty} Ri_g(\zeta) = \frac{1}{\beta} \equiv Ri_c = 0.2
   $$
   This asymptote is the **critical Richardson number**.  Explain why $Ri_g$ can never exceed
   $Ri_c$ within the MOST framework — i.e., what is the physical meaning of $\zeta \to \infty$?

3. Invert to find:
   $$
   \zeta = \frac{Ri_g}{1 - \beta\,Ri_g}, \qquad Ri_g < Ri_c
   $$
   What does $\zeta \to \infty$ as $Ri_g \to Ri_c^-$ imply about the Obukhov length scale
   relative to the measurement height?

4. The **Miles–Howard theorem** requires $Ri_g < 1/4$ for Kelvin–Helmholtz instability.
   Compare $Ri_c = 0.20$ to $1/4 = 0.25$.  Are these the same threshold?  Which controls
   sustained turbulence in the SBL, and why?

### Part B — Stability Function in $Ri$-Space: England–McNider (1995) Form (8 points)

The **stability reduction function** $f_m(Ri) \equiv 1/\phi_m^2$ enters the eddy diffusivity as
$K_M = \kappa u_* z\,f_m(Ri)$.

1. Substitute $\zeta = Ri/(1-\beta Ri)$ into $f_m = (1+\beta\zeta)^{-2}$ and simplify to:
   $$
   f_m(Ri) = \left(1 - \frac{Ri}{Ri_c}\right)^2, \quad 0 \leq Ri < Ri_c; \qquad f_m = 0, \quad Ri \geq Ri_c
   $$

2. Sketch $f_m(Ri)$.  Label the values at $Ri = 0$, $Ri_c/2$, and $Ri_c$.  Is the curve convex
   or concave?  What does this imply about the rate at which turbulence shuts down as $Ri$
   approaches $Ri_c$?

3. Show that the same derivation gives $f_h(Ri) = (1 - Ri/Ri_c)^2$ when $\phi_h = \phi_m$.
   How would $f_h$ differ if $\phi_h = 1 + \gamma\zeta$ with $\gamma \neq \beta$ (i.e., turbulent
   Prandtl number $\neq 1$)?

### Part C — The Laminar Regime (7 points)

When $Ri > Ri_c$, turbulence is fully suppressed and only molecular diffusion acts.

1. At $z = 5\,\text{m}$, $u_* = 0.30\,\text{m s}^{-1}$: compute turbulent $K_M = \kappa u_* z$.
   Compare to $\nu \approx 1.5 \times 10^{-5}\,\text{m}^2\,\text{s}^{-1}$.  By how many orders
   of magnitude does turbulence enhance transport relative to molecular diffusion?

2. Molecular thermal diffusivity $\alpha \approx 2.1 \times 10^{-5}\,\text{m}^2\,\text{s}^{-1}$.
   The depth of a diffusive layer growing over time $t$ scales as $\delta \sim \sqrt{\alpha t}$.
   Over a 6-hour night, how deep a layer does molecular conduction significantly affect?
   Compare to a typical stable BL depth of 50–200 m.

3. **Intermittent turbulence:** Real stable BLs exhibit bursts of turbulence even when the
   mean $\overline{Ri_g} > Ri_c$.  Name **two mechanisms** that can regenerate turbulence in
   this regime and describe each in one sentence.

### Part D — Turbulence Collapse: The Flip-Flop Transition (5 points)

1. Sketch $|H_s|$ vs.\ $Ri_B$ for the full stability range, showing:
   (a) the turbulent branch from $Ri = 0$ to $Ri_c$ (using $H_s \propto f_h(Ri)\,\Delta\Theta$),
   and (b) the laminar branch ($H_s \approx 0$) for $Ri > Ri_c$.

2. **Positive feedback (runaway cooling):** Starting from a state near $Ri \lesssim Ri_c$,
   explain in 4–5 sentences how a small perturbation in surface cooling can trigger an
   irreversible collapse to the laminar state.  Address: (a) how reduced turbulence changes $Ri$,
   (b) why the process is self-reinforcing, and (c) under what external conditions (wind speed,
   clouds, soil) the collapse is most likely.  This is the "flip-flop" transition of McNider
   et al.

3. Existence of two stable equilibria (turbulent and laminar) under some forcing implies
   path-dependence in the solution.  What does this mean for the ability of a deterministic NWP
   model to predict the nocturnal surface temperature?  Discuss requirements on model resolution
   and the sensitivity of the forecast to the value of $\beta$ (and hence $Ri_c$).

---

## Summary: Key Formulas

| Quantity | Neutral | Stability-modified |
|---|---|---|
| Wind profile | $U = \frac{u_*}{\kappa}\ln(z/z_0)$ | $U = \frac{u_*}{\kappa}[\ln(z/z_0) - \psi_m(\zeta)]$ |
| Temp profile | $\Theta - \Theta_s = \frac{\theta_*}{\kappa}\ln(z/z_{0h})$ | $\Theta - \Theta_s = \frac{\theta_*}{\kappa}[\ln(z/z_{0h}) - \psi_h(\zeta)]$ |
| Wind shear | $\partial U/\partial z = u_*/(\kappa z)$ | $\partial U/\partial z = (u_*/\kappa z)\,\phi_m(\zeta)$ |
| Temp gradient | $\partial \Theta/\partial z = \theta_*/(\kappa z)$ | $\partial \Theta/\partial z = (\theta_*/\kappa z)\,\phi_h(\zeta)$ |
| Eddy diffusivity | $K_M = \kappa u_* z$ | $K_M = \kappa u_* z\,f_m(Ri)$ |
| Gradient $Ri_g$ | $(g\kappa\theta_*)/(u_*^2 \overline{\Theta})$ — height independent | — |
| $Ri_g(\zeta)$ stable linear | — | $\zeta/(1+\beta\zeta)$; asymptote $Ri_c = 1/\beta = 0.2$ |
| $Ri_g(\zeta)$ unstable B-D | — | $Ri_g = \zeta$ (exact) |
| EM stability function | — | $f_m = (1-Ri/Ri_c)^2$, $Ri < Ri_c$ |
| Stable $\psi_m$ (log-linear) | — | $\psi_m = -\beta\zeta$ |
| Unstable $\psi_m$ (Paulson) | — | $2\ln\frac{1+x}{2}+\ln\frac{1+x^2}{2}-2\arctan x+\frac{\pi}{2}$, $x=(1{-}16\zeta)^{1/4}$ |
| Friction velocity | $u_* = \kappa\Delta U/\ln(z_2/z_1)$ | Iterative solution needed with $\psi_m$ |
| Drag coefficient | $C_D = \kappa^2/\ln^2(z_r/z_0)$ | $C_D = \kappa^2/[\ln(z_r/z_0)-\psi_m]^2$ |
| Momentum flux | $\tau = u_*^2$ | $\tau = u_*^2$ |
| Sensible heat flux | $H_s = -u_*\theta_*$ | $H_s = -u_*\theta_*$ |

---

## Grading

| Problem | Points |
|---------|--------|
| 1 — Wind shear | 20 |
| 2 — Temperature gradient and heat flux | 20 |
| 3 — Richardson number | 25 |
| 4 — Inverse problem and bulk coefficients | 20 |
| 5 — Grid discretization error | 15 |
| 6 — Unstable surface layer ($\psi$-functions, B-D forms) | 25 |
| 7 — Near-neutral regime and log-linear profiles | 20 |
| 8 — Stable BL, critical $Ri$, and laminar transition | 30 |
| **Total** | **175** |

---

*Prepared for ABL Physics — Graduate Level*
