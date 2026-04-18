I’m treating this as a model-design brainstorming task, starting by gathering context from the scaffold file and nearby project notes so I can frame a minimal but defensible slab model instead of inventing something disconnected from the repo.I’m checking the scaffold and a few relevant notes in parallel to anchor the brainstorming in whatever structure and terminology you already use.

Read [](file:///Users/davidengland/Documents/GitHub/ABL/param/SCAFFOLDING.md#1-1), lines 1 to 250

Searched for files matching `ABL` · `**/*slab*`, no matches

Searched for regex `snow|ice|slab|surface energy balance|heat flux|tower|machine learning|self calibrat` (`**/*.{md,tex,py}`), 80 results

I found the existing scaffold is surface-layer focused and already has the right abstractions for regime switching, diagnostics, and tower validation. Next I’m pulling just the most relevant notes on regimes and snow/ice extensions so I can map your slab idea onto that architecture without overbuilding it.

Read [](file:///Users/davidengland/Documents/GitHub/ABL/topics.md#60-60), lines 60 to 150

Read [](file:///Users/davidengland/Documents/GitHub/ABL/param/core/regimes.md)

Read [](file:///Users/davidengland/Documents/GitHub/ABL/CHANGELOG.md#30-30), lines 30 to 53

**Core Shape**

The simplest defensible path is to treat this as an extension of your existing surface-layer scaffold, not a separate model family. Keep the air side in SCAFFOLDING.md, and add a 1D substrate column underneath with a small number of prognostic slabs.

Use one common conservation form for any layer variable $\chi$:

$$
\frac{d}{dt}\left(\rho_i h_i \chi_i\right)=J_{i-\frac12}-J_{i+\frac12}+S_i
$$

For temperature, write it in enthalpy form so melting/freezing does not need special-case logic:

$$
\frac{d}{dt}\left(\rho_i h_i H_i\right)=F_{i-\frac12}-F_{i+\frac12}+Q_i
$$

with $H(T)$ carrying sensible plus latent heat. That gives you a clean way to handle snow, ice, wet snow, and bare ground with the same machinery.

The minimum useful state is:

1. Snow thickness $h_s \ge 0$.
2. One to three snow layers, each with $T_i$, density, conductivity, heat capacity.
3. One underlying substrate mode:
   land mode: soil slab or prescribed deep-soil boundary.
   ice-ocean mode: ice slab plus optional thin under-ice mixed layer.
4. Surface skin temperature $T_s$ diagnosed from the energy balance.

If $h_s = 0$, the snow block simply vanishes and the air couples directly to ice or land. Do not allow negative thickness; clamp with a prognostic source budget:
$\dot h_s = P_{snow} - M - S_{sub} - C_{comp}$.

**Flux Architecture**

Split the problem into two interfaces.

Air to surface:
Use your existing MOST-based framework for momentum, heat, and generic atmospheric scalar exchange. For heat:

$$
H = \rho_a c_p C_H U \left(\theta_a - T_s\right)
$$

For a generic atmospheric quantity $q$:

$$
F_q = \rho_a C_q U (q_a - q_s)
$$

where $q_s$ is either a diagnosed surface value or a resistance-limited value. This is the right place to use the regime logic already outlined in regimes.md, including stable, near-neutral, and laminar behavior.

Within and below the snow:
Use pure vertical diffusion/conduction first.

$$
F_T = -k \frac{\partial T}{\partial z}, \qquad
J_q = -K_q \frac{\partial q}{\partial z}
$$

That is enough for a first model. Do not try to make “generic $q$ in snow” too universal. Physically, different scalars behave differently:

1. Heat always diffuses through snow and ice.
2. Water vapor can diffuse in porous snow.
3. Salt matters below sea ice, not in dry land snow.
4. Most atmospheric chemistry tracers should be handled as surface deposition/resistance, not full snow-column transport, unless you specifically want snowpack chemistry.

For the “fresh layer below ice” idea, that is real enough to justify one optional thin sub-ice slab in ocean mode. Salt rejection during freezing can leave a sharp salinity structure, so the minimal extension is:

1. Under-ice layer temperature $T_w$.
2. Under-ice salinity $S_w$.
3. A mixing coefficient to the deeper ocean.

That gives you a way to represent a fresh or stratified under-ice layer without building a full ocean model.

**How To Keep It Simple But Realistic**

The biggest modeling mistake here would be to over-resolve fuzzy structure. Since real layers are not homogeneous, handle that with effective properties and occasional layer merge/split rules, not with many thin layers from the start.

A practical hierarchy is:

1. MVP:
   one snow slab + one substrate slab + diagnosed surface temperature.
2. Better:
   two or three snow slabs + enthalpy phase change + optional ice slab.
3. Arctic mode:
   add thin under-ice salinity/temperature slab.
4. Advanced:
   adaptive layer splitting/merging when gradients exceed a threshold.

For fuzzy interfaces, use blended material properties over a finite transition thickness instead of a hard discontinuity. In practice, that can be as simple as an interface conductivity based on harmonic means. You do not need a separate “fuzzy boundary model” at first.

I would also recommend a small state-machine for surface type:

1. Bare land
2. Snow over land
3. Bare ice
4. Snow over ice
5. Ice over ocean with optional under-ice stratified slab

That keeps the code organized and prevents impossible combinations.

**ML Role**

Do not let the ML module replace the physics. It should calibrate uncertain parameters online from tower data at a fixed site. That is much more defensible.

Let ML adjust things like:

1. $z_0$, $z_{0h}$, roughness Reynolds effects
2. Snow conductivity multiplier
3. Surface resistance for moisture or scalar exchange
4. Compaction and metamorphism rates
5. Under-ice exchange coefficient
6. Bias corrections to radiative forcing or near-surface stability classification

Keep hard physical constraints:

1. Energy conservation
2. Mass conservation
3. $h_s \ge 0$
4. Positive diffusivities and conductivities
5. Temperature and salinity bounds

A good first “ML” implementation is not a neural net. It is site-specific adaptive parameter estimation with memory, for example EnKF, recursive least squares, or Bayesian updating. Neural surrogates can come later if you want fast emulation.

**Recommended First Build**

If you want the fastest credible prototype, build this:

1. Existing MOST air-side fluxes from SCAFFOLDING.md.
2. One prognostic snow slab with enthalpy.
3. One underlying slab:
   soil for land mode, ice for sea-ice mode.
4. Optional fixed deep boundary:
   deep-soil temperature or ocean mixed-layer temperature.
5. Surface energy balance with shortwave, longwave, sensible, latent, and conductive flux.
6. Optional tower-based adaptive calibration of $z_{0h}$ and snow conductivity.

That is simple enough to implement and still realistic enough to capture the main behaviors you care about: insulation by snow, stable coupling to the atmosphere, zero-snow transitions, and different lower-boundary physics for land versus sea ice.

If useful, I can turn this into a concrete design note in the repo, for example a new param document with governing equations, state variables, and pseudocode.