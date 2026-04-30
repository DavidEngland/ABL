# Ultraspherical Harmonics for ABL: From Idea to Testable Science

## Why this is exciting

Ultraspherical (Gegenbauer) structure already appears in your MOST work through exponents that map to orthogonal families. That means you can treat stability functions not only as curve fits, but as a spectral basis with physically interpretable modes.

The short version:

- Earth crust reconstruction from many points is an existence proof that spectral basis methods can recover structure from sparse, irregular observations.
- ABL stratification can be treated similarly, but now with shear, heat, and tracer transfer coefficients across height and stability.
- This gives a path to local-to-global scaling, including Arctic amplification diagnostics and transfer to other planets/moons.

## Unifying concept

Represent closure behavior as a sum over orthogonal modes with stability-aware weights:

$$
F(\zeta, x, t) = \sum_{n=0}^{N} a_n(x,t)\,C_n^{(\lambda)}(\xi(\zeta))
$$

where:

- $F$ can be momentum, heat, or tracer transfer function,
- $C_n^{(\lambda)}$ are Gegenbauer polynomials,
- $\xi(\zeta)$ maps stability coordinate into the polynomial domain,
- coefficients $a_n$ carry space-time dynamics and data constraints.

This naturally links to your existing results on $\phi_m$, $\phi_h$, and Richardson mappings.

## Local to global: how to add more points

Think in three nested scales.

1. Local node (single met station)
- Inputs: wind, temperature, humidity, pressure, turbulence diagnostics if available.
- Output: local fitted modal coefficients and uncertainty.

2. Regional network (multiple stations, towers, buoys)
- Add data assimilation or regularized interpolation on coefficients, not raw variables.
- This smooths noise while preserving stability transitions.

3. Global/planetary shell
- Expand coefficients on the sphere with spherical/ultraspherical harmonics.
- Couple near-surface ABL shell to large-scale circulation constraints.

Practical note: start with station clusters in Alaska/Arctic to target amplification questions directly.

## Arctic amplification link (spelling corrected)

Arctic amplification can be reframed as a modal energy redistribution problem:

- Low-order modes: background stratification and seasonal forcing.
- Intermediate modes: synoptic variability and shear production shifts.
- High-order modes: intermittency, surface heterogeneity, and local decoupling events.

Testable hypothesis:

$$
\text{Amplification signal} \propto \Delta E_{\text{stable ABL modes}} + \Delta E_{\text{intermittent coupling modes}}
$$

If true, amplification is not only mean-state warming, but a change in modal partition across stability regimes.

## Fractional dimension and turbulence

Use effective/fractal dimension as a state variable for unresolved transport complexity.

Candidate definitions:

1. Structure-function scaling:
$$
S_q(r) = \langle |\Delta u(r)|^q \rangle \sim r^{\zeta_q},
$$
then infer intermittency and effective dimension from $\zeta_q$ curvature.

2. Spectral slope proxy:
$$
E(k) \sim k^{-p},
$$
with regime-dependent mapping between $p$ and effective dimension.

3. Information or correlation dimension from reconstructed phase-space attractors in high-frequency tower data.

For tracers, dimension may differ by species because chemistry, phase change, and deposition alter cascade pathways.

## Determining dimension for other tracers

Workflow per tracer (water vapor, CO2, methane, aerosols, reactive species):

1. Build synchronized high-frequency time series with quality control.
2. Compute multi-scale structure functions and wavelet spectra.
3. Estimate intermittency metrics and effective dimension with confidence intervals.
4. Regress dimension against stability class ($Ri_g$, $\zeta$), shear, and surface state.
5. Feed inferred dimension into tracer-specific closure terms.

Result: tracer-dependent shear/transport functions rather than one-size-fits-all similarity.

## Other planets and moons

The framework transfers cleanly by replacing forcing and fluid constants:

- gravity $g$
- atmospheric composition and molecular diffusivities
- surface roughness and thermal inertia
- radiative forcing cycle length

Targets with high payoff:

1. Mars: strong diurnal cycle, dust feedback, thin atmosphere.
2. Titan: thick atmosphere, methane cycle, exotic stable layers.
3. Venus upper cloud deck: super-rotation coupling questions.
4. Europa/Triton-style tenuous cases as edge tests for limit behavior.

## Shear functions roadmap (momentum, heat, tracers)

Generalize closure family to:

$$
\phi_{\chi}(\zeta) = \left(1 - b_{\chi}\zeta\right)^{-\alpha_{\chi}} G_{\chi}(D_{\chi}, I_{\chi})
$$

where $\chi \in \{m,h,q,c,...\}$ for momentum, heat, humidity, carbon, etc., and $D_{\chi}$ is effective dimension, $I_{\chi}$ intermittency.

This preserves classic MOST limits while adding scientifically meaningful degrees of freedom.

## Highest-yield direction forward

If your goal is fastest scientific return, do this first:

1. Build Arctic pilot with 2-3 station clusters and one reanalysis backbone.
2. Fit ultraspherical modal representation for momentum and heat.
3. Add one tracer (humidity) with tracer-specific dimension estimate.
4. Evaluate whether modal shifts explain Arctic amplification metrics better than baseline MOST.

Success criterion:

- improved prediction of stable ABL transitions,
- reduced bias in flux-gradient relations,
- interpretable modal diagnostics linked to amplification episodes.

## Homework-style deliverables

1. Derive the mapped Gegenbauer basis for your chosen $\zeta \to \xi$ transform.
2. Show identifiability of parameters $(\alpha, b, \lambda)$ under noisy station data.
3. Implement synthetic test where known coefficients are recovered from sparse observations.
4. Add real station case and compare against standard MOST.
5. Report uncertainty propagation into $Ri_g$ and flux estimates.

## Immediate next step this week

Create one baseline notebook that:

1. Ingests one Arctic station dataset,
2. Fits classic MOST and ultraspherical form side by side,
3. Computes one effective-dimension diagnostic,
4. Produces one figure: observed vs modeled transfer function across stability bins.

That single notebook will convert the idea into a falsifiable scientific program.

## Technical Q and A: tanh mapping, inversion, and SBL limits

### 1) Mapping to polynomial space with tanh

Yes, using a hyperbolic tangent map is a good default because it compresses an unbounded or semi-bounded stability coordinate into a bounded polynomial domain.

One practical choice is

$$
\xi = \tanh(\alpha_{\xi}\zeta), \qquad \xi \in (-1,1),
$$

with inverse

$$
\zeta = \frac{1}{\alpha_{\xi}}\operatorname{artanh}(\xi)
= \frac{1}{2\alpha_{\xi}}\ln\left(\frac{1+\xi}{1-\xi}\right).
$$

Notes:

1. $\alpha_{\xi}$ is a mapping-scale parameter controlling where neutral and strongly stable/unstable regimes fall in $\xi$.
2. The logarithmic inverse is useful for analytic asymptotics and Jacobian-aware fitting.
3. For branch handling near singular points, artanh forms are numerically safer than direct inversion of polynomial surrogates.

### 2) Can we choose $\alpha = 1$ without loss of generality?

Not in general. You can normalize amplitude constants without loss of generality, but exponent choices carry physical information.

For

$$
\phi_q(\zeta) = a_q\left(1-b_q\zeta\right)^{-1/\lambda_q},
$$

setting $\alpha=1$ (or equivalently fixing $1/\lambda_q=1$) changes curvature and asymptotic scaling unless your data are exactly consistent with that exponent. So:

1. WLOG for amplitude: usually yes.
2. WLOG for exponent: generally no.

Recommended practice is to fit or constrain exponent families by physics, then test identifiability.

### 3) artanh, coth, acoth and reciprocal behavior

Your note is well aligned with implementation practice.

1. artanh gives a log representation that is convenient for inversion and uncertainty propagation.
2. acoth can represent reciprocal-dominant regimes because
$$
\operatorname{acoth}(x)=\operatorname{artanh}(1/x), \quad |x|>1.
$$
3. coth/acoth formulations can stabilize numerical mappings when transformed variables naturally evolve in reciprocal coordinates.

For ABL closures, use these as coordinate transforms and inversion tools, not as replacements for physically constrained closure forms.

### 4) Single met station: implementing ultraspherical closure from curve-fit MOST

Given station time series:

1. Fit baseline MOST parameters $(a_q,b_q,\lambda_q)$ in your chosen regime split.
2. Map $\zeta \to \xi=\tanh(\alpha_{\xi}\zeta)$.
3. Expand residual or full transfer function in Gegenbauer basis:
$$
\phi_q(\zeta(t)) \approx \sum_{n=0}^{N} c_{n,q} C_n^{(\lambda_*)}(\xi(t)).
$$
4. Choose $N$ by out-of-sample validation and penalize high-order coefficients.
5. Report both physical parameters $(a,b,\lambda)$ and modal diagnostics $(c_n)$.

Interpretation:

- baseline parameters explain canonical scaling,
- modal coefficients capture local departures, intermittency, and regime transitions.

### 5) Cluster of met stations around the globe

Use a two-level expansion:

1. Vertical-stability basis at each site (ultraspherical in $\xi$).
2. Horizontal spherical basis for coefficient fields:
$$
c_{n,q}(\theta,\varphi,t)=\sum_{\ell,m} d_{n,q,\ell m}(t)Y_{\ell m}(\theta,\varphi).
$$

This is where ultraspherical methods help most:

1. Data assimilation on coefficient space (lower dimension, smoother priors).
2. Multi-scale simulations with interpretable mode energy transfer.
3. Visualization of global maps of low-order vs high-order mode power, highlighting Arctic transitions.

### 6) Domain of validity into SBL and beyond

Your statement is essentially right for the unstable-style power form.

For

$$
\phi_q(\zeta)=a_q(1-b_q\zeta)^{-1/\lambda_q},
$$

the branch point at $\zeta=1/b_q$ defines a mathematical boundary for that analytic branch. If that form is calibrated on the unstable side, it does not automatically specify behavior beyond that point on stable branches.

Implications:

1. It provides no guaranteed physical closure beyond the branch/singularity scale.
2. SBL usually requires separate stable formulations, blended transitions, or dynamic-$Ri_c$ constraints.
3. A practical approach is piecewise or smoothly matched closures with continuity of value and first derivative at regime boundaries.

### 7) Suggested student development path

1. Homework A: derive tanh mapping and Jacobian-weighted orthogonality.
2. Homework B: fit single-station MOST then ultraspherical residual model.
3. Homework C: compare $N=2,4,6$ truncations under cross-validation.
4. Project 1: Arctic station pair modal-comparison dashboard.
5. Project 2: global station cluster coefficient-map and mode-energy animation.
6. Project 3: tracer-dependent effective-dimension closure challenge.