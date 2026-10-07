# Biharmonic and Higher-Order Hyperdiffusion in Atmospheric Boundary-Layer Modeling

**To:** Dick McNider  
**From:** David England  
**Date:** 7 October 2026  
**Subject:** Scale-selective numerical dissipation, vertical-grid implementation, and similarity-coordinate diagnostics

## Executive Summary

Biharmonic hyperdiffusion is a fourth-order spatial filter that damps short wavelengths more strongly than long wavelengths. Triharmonic (sixth-order) and higher even-order filters sharpen this scale selectivity when compared at the same grid-scale damping time. They are numerical dissipation mechanisms, not substitutes for physically based turbulent transport.

For stable boundary-layer applications, the main concern is distinguishing numerical oscillations from resolved inversions, shallow shear layers, and low-level jets. A stronger or higher-order filter does not make that distinction automatically. Vertical filtering needs particular care because the smallest grid intervals and the sharpest physical gradients often occur together near the surface.

Three conclusions guide the discussion:

1. **Use an unambiguous dissipative sign convention.** Write the order-$2p$ tendency as $-\nu_{2p}(-\Delta)^p q$, where $\Delta=\nabla^2$ and $\nu_{2p}\geq0$. The sixth-order term is then $+\nu_6\Delta^3q$, whereas the fourth-order term is $-\nu_4\Delta^2q$.
2. **Construct the discrete filter to dissipate a defined norm.** On stretched grids, conservation, grid metrics, boundary conditions, and time integration matter at least as much as formal derivative order.
3. **Use similarity-coordinate derivatives diagnostically.** Bell-polynomial identities explain how closure derivatives and coordinate geometry contribute to physical-height derivatives. They do not, by themselves, establish that a filter is stabilizing locally or that a coordinate fold represents turbulence collapse.

The recommended starting experiment is a conservative, energy-dissipative fourth-order vertical filter, followed by a sixth-order comparison with matched shortest-wave damping. Higher orders should be considered only if those comparisons demonstrate a need.

## 1. Notation and Scope

| Symbol | Meaning | Units |
| --- | --- | --- |
| $q(\mathbf{x},t)$ | Model field to be filtered, such as a velocity component or potential temperature | Field dependent |
| $z$ | Physical height | m |
| $\Delta=\nabla^2$ | Laplacian in the specified physical directions | $\mathrm{m}^{-2}$ |
| $p$ | Positive integer; the spatial derivative order is $2p$ | Dimensionless |
| $\nu_{2p}$ | Constant hyperdiffusion coefficient in the basic analysis | $\mathrm{m}^{2p}\,\mathrm{s}^{-1}$ |
| $k$ | Physical angular wavenumber | $\mathrm{m}^{-1}$ |
| $\ell=2\pi/k$ | Wavelength | m |
| $\tau_g$ | Prescribed e-folding damping time for a chosen grid-scale mode | s |
| $Ri_g$ | Gradient Richardson number | Dimensionless |
| $\zeta=z/L$ | Similarity variable, where the Obukhov length $L$ is finite and nonzero | Dimensionless |
| $R(\zeta)$ | Similarity representation of $Ri_g$ | Dimensionless |
| $R^{(r)}$ | $\mathrm{d}^rR/\mathrm{d}\zeta^r$ | Dimensionless |
| $g_j$ | $\mathrm{d}^j\zeta/\mathrm{d}z^j$ | $\mathrm{m}^{-j}$ |

The notation $g_j$ avoids using $\partial^j$ for a coordinate derivative value: $\partial$ normally denotes an operator, not a coefficient. Physical-height derivatives are written explicitly as $\mathrm{d}^m Ri_g/\mathrm{d}z^m$.

For vertical filtering alone, replace $\Delta$ by $\partial_z^2$. Horizontal filtering instead uses $\Delta_h$. A three-dimensional operator is not interchangeable with independently applied horizontal and vertical filters: powers of $\Delta_h+\partial_z^2$ contain mixed derivatives. The analysis below assumes constant coefficients unless stated otherwise.

## 2. Why Biharmonic Diffusion Damps Short Waves

Consider the numerical-filter contribution to a prognostic equation:

$$
\left.\frac{\partial q}{\partial t}\right|_{\mathrm{filter}}
=\mathcal{H}_{2p}q,
\qquad
\mathcal{H}_{2p}=-\nu_{2p}(-\Delta)^p.
$$

For a Fourier mode with amplitude $\widehat q$ and wavenumber magnitude $k$, $\Delta$ has eigenvalue $-k^2$. Therefore,

$$
\frac{\mathrm{d}\widehat q}{\mathrm{d}t}
=-\nu_{2p}k^{2p}\widehat q,
\qquad
\widehat q(t)=\widehat q(0)\exp[-\nu_{2p}k^{2p}t].
$$

All nonzero modes decay, but shorter wavelengths decay faster. The constant mode is unchanged.

| Spatial order | Name | Dissipative tendency | Fourier decay rate |
| --- | --- | --- | --- |
| 2 | Laplacian diffusion | $+\nu_2\Delta q$ | $\nu_2k^2$ |
| 4 | Biharmonic hyperdiffusion | $-\nu_4\Delta^2q$ | $\nu_4k^4$ |
| 6 | Triharmonic hyperdiffusion | $+\nu_6\Delta^3q$ | $\nu_6k^6$ |
| 8 | Eighth-order hyperdiffusion | $-\nu_8\Delta^4q$ | $\nu_8k^8$ |
| $2p$ | General even-order hyperdiffusion | $(-1)^{p+1}\nu_{2p}\Delta^pq$ | $\nu_{2p}k^{2p}$ |

Thus, $-\nu_6\Delta^3q$ would amplify Fourier modes rather than damp them under this definition of $\Delta$. Expressions such as $-\nu_6\nabla^6q$ should be avoided unless the operator convention is specified.

### Comparing Orders Fairly

Let $k_g$ denote the chosen grid-scale reference wavenumber. Set

$$
\nu_{2p}=\frac{1}{\tau_g k_g^{2p}}.
$$

The decay rate at any other wavenumber is then

$$
\sigma_{2p}(k)=\frac{1}{\tau_g}\left(\frac{k}{k_g}\right)^{2p}.
$$

At half the reference wavenumber, the decay rates relative to the reference rate are $1/4$, $1/16$, $1/64$, and $1/256$ for orders 2, 4, 6, and 8. This is the advantage of higher order: less damping of longer waves for the same shortest-wave attenuation. It is not a perfect spectral cutoff, and all nonzero wavenumbers still experience some damping.

For a discrete Laplacian, calibrate using its actual eigenvalue rather than automatically substituting the continuum value $k^2$. On a uniform one-dimensional grid with the centered three-point Laplacian,

$$
\mu(k)=\frac{4}{h^2}\sin^2\!\left(\frac{kh}{2}\right),
\qquad
\sigma_{2p}(k)=\nu_{2p}\mu(k)^p.
$$

At the shortest represented wavelength $2h$, $\mu_g=4/h^2$ and the matched coefficient is $\nu_{2p}=1/(\tau_g\mu_g^p)$. This discrete calibration differs from using $k_g=\pi/h$ in the continuum formula.

## 3. Dissipation Is an Integral Property, Not a Pointwise Guarantee

For periodic boundaries, or boundary conditions that make $A=-\Delta$ a nonnegative self-adjoint operator, define the quadratic norm $E_q=\tfrac12\langle q,q\rangle$. Then

$$
\left.\frac{\mathrm{d}E_q}{\mathrm{d}t}\right|_{\mathrm{filter}}
=-\nu_{2p}\langle q,A^pq\rangle\leq0.
$$

In one dimension, with boundary terms eliminated, the same result is

$$
\left.\frac{\mathrm{d}E_q}{\mathrm{d}t}\right|_{\mathrm{filter}}
=-\nu_{2p}\int\left(\frac{\mathrm{d}^pq}{\mathrm{d}z^p}\right)^2\mathrm{d}z.
$$

For a velocity component this relates to kinetic-energy dissipation, with the appropriate density weighting in an atmospheric implementation. For potential temperature it is a quadratic scalar norm, not total thermodynamic energy.

The tendency can be positive at one height and negative at another while the domain-integrated norm decreases. Conversely, finding a nonzero fourth or sixth derivative at one location does not prove that the local tendency is restoring. Fourth- and higher-order filters generally lack the maximum principle of ordinary diffusion: new extrema and negative values can occur. Positivity of TKE, moisture, or other constrained variables requires separate attention, and a limiter must be assessed for its effect on conservation and dissipation.

## 4. Physical Transport Versus Numerical Filtering

Turbulent transport is commonly represented by a second-order flux-divergence term such as

$$
\left.\frac{\partial q}{\partial t}\right|_{\mathrm{turb}}
=\frac{\partial}{\partial z}\left(K_q\frac{\partial q}{\partial z}\right).
$$

Here $K_q$ is determined by a turbulence closure. A hyperdiffusion coefficient is instead selected to control numerical scales. The two contributions should be kept separate in budgets and sensitivity tests.

The gradient Richardson number is usually a diagnostic:

$$
Ri_g=\frac{N^2}{S^2},
\qquad
N^2=\frac{g}{\theta_{v,\mathrm{ref}}}\frac{\partial\theta_v}{\partial z},
\qquad
S^2=\left(\frac{\partial u}{\partial z}\right)^2
+\left(\frac{\partial v}{\partial z}\right)^2.
$$

Directly filtering $Ri_g$ is not equivalent to filtering $u$, $v$, and $\theta_v$ and then recomputing it. For a dynamical experiment, apply the chosen filter to appropriate prognostic fields and recompute diagnostics from the resulting state. If smoothing $Ri_g$ is used only to condition a closure, identify it as a separate closure modification and assess its effect on transport.

Near weak vector shear, $Ri_g$ can become very large even for a smooth state. A wind-speed maximum alone does not imply $S^2=0$, because wind direction may still change with height. Such features should not automatically be classified as grid noise.

## 5. Vertical Grids, Boundaries, and Time Integration

### Stretched-Grid Construction

A practical route is to start with a conservative discrete vertical Laplacian $L_h$, including the intended boundary treatment. Let $M$ contain positive cell-volume or mass weights and set $A_h=-L_h$. Require

$$
A_h\mathbf{1}=0,
\qquad
MA_h=A_h^{\mathsf T}M,
\qquad
\mathbf{q}^{\mathsf T}MA_h\mathbf{q}\geq0.
$$

Then $-\nu_{2p}A_h^p\mathbf q$ conserves the weighted mean and dissipates the weighted quadratic norm. These properties must be verified for the actual staggered grid and boundary closures, not inferred from the stencil order. Repeated application can implement the power without storing a wide stencil, but still involves a wider effective stencil and additional boundary or halo work.

Computational-grid stretching is distinct from the similarity mapping $\zeta(z)$. A derivative in a grid index is not a derivative in physical height; metric factors and their derivatives cannot be omitted.

If $\nu_{2p}$ varies vertically, simply multiplying $A_h^p$ by a diagonal coefficient matrix need not preserve dissipation or conservation. One alternative is a compatible weighted-adjoint construction,

$$
\mathcal H_h=-M^{-1}D_p^{\mathsf T}W\operatorname{diag}(\nu_{2p})D_p,
$$

where $D_p$ approximates the $p$th physical-height derivative, $W$ supplies positive quadrature weights, and boundary treatment is included consistently. Nonnegative coefficients give norm dissipation; $D_p\mathbf1=0$ additionally gives weighted-mean conservation. With constant coefficients, its interior continuum counterpart is the intended order-$2p$ operator, but its boundary realization must still be specified.

### Boundary Conditions

A direct order-$2p$ differential equation requires $2p$ scalar boundary conditions in one dimension. A powered discrete Laplacian instead inherits a particular boundary realization from its construction; that realization must be checked against the physical problem. Ordinary turbulent surface-flux conditions alone do not uniquely specify an added fourth- or sixth-order operator. Numerical-filter boundary fluxes must not silently change the imposed surface momentum or heat flux.

### Time Stepping

For forward Euler applied to the filter alone,

$$
\Delta t\,\nu_{2p}\lambda_{\max}(A_h)^p\leq2.
$$

On a uniform one-dimensional grid, $\lambda_{\max}\leq4/h^2$, giving $\Delta t\leq2h^{2p}/(4^p\nu_{2p})$. Other time integrators have different stability limits. The familiar $h^{2p}$ restriction assumes a fixed coefficient; if the coefficient is retuned with resolution to maintain a fixed grid-scale damping time, that scaling changes. Small cells on stretched grids still require checking the actual spectral bound.

Operator splitting does not remove this restriction if the filter step remains explicit. A backward-Euler filter step,

$$
\left(I+\Delta t\,\nu_{2p}A_h^p\right)\mathbf q^{n+1}=\mathbf q^*,
$$

is unconditionally stable for this linear dissipative subproblem, although accuracy and coupling to the remaining dynamics still constrain the step. It requires a linear solve, not necessarily a nonlinear solve. Crank-Nicolson is also linearly stable, but strongly damped modes can change sign between steps. GPU suitability depends on the actual stencil, communication, and solver costs; branch-free expressions alone do not establish an efficient implementation.

## 6. Similarity Coordinates and Coordinate Folds

Within a smooth local representation $Ri_g(z)=R(\zeta(z))$, the higher-order chain rule gives

$$
\frac{\mathrm{d}^mRi_g}{\mathrm{d}z^m}
=\sum_{r=1}^{m}R^{(r)}(\zeta)B_{m,r}(g_1,g_2,\ldots,g_{m-r+1}),
$$

where $B_{m,r}$ are partial exponential Bell polynomials. This representation is exact under the stated composition assumption; actual profiles with additional independent height-dependent closure parameters require their own chain-rule terms.

For fourth order,

$$
\frac{\mathrm{d}^4Ri_g}{\mathrm{d}z^4}
=R'g_4+R''(4g_1g_3+3g_2^2)
+6R'''g_1^2g_2+R^{(4)}g_1^4.
$$

If $L$ is constant, $g_1=1/L$ and all $g_j$ for $j\geq2$ vanish. The transformation is simply $\mathrm{d}^mRi_g/\mathrm{d}z^m=L^{-m}R^{(m)}$; there is no coordinate fold. A stretched numerical grid does not change this fact.

If a height-dependent local $L(z)$ is introduced, then

$$
g_1=\frac{1}{L}\left(1-\zeta\frac{\mathrm{d}L}{\mathrm{d}z}\right).
$$

A nondegenerate fold at $z=z_f$ has $g_1(z_f)=0$ and $g_2(z_f)\neq0$. This is a local stationary point of the similarity mapping, not a singularity of physical height. Extending classical constant-flux MOST using $L(z)$ requires physical justification; the algebra alone does not validate the extension.

### Fourth- and Sixth-Order Terms at a Fold

At the fold, the exact expressions reduce to

$$
\left.\frac{\mathrm{d}^4Ri_g}{\mathrm{d}z^4}\right|_{z_f}
=R'g_4+3R''g_2^2,
$$

$$
\left.\frac{\mathrm{d}^6Ri_g}{\mathrm{d}z^6}\right|_{z_f}
=R'g_6+R''(15g_2g_4+10g_3^2)+15R'''g_2^3.
$$

Every factor on the right is evaluated at the fold. Fourth-order differentiation can retain similarity curvature $R''$ even though the first physical derivative vanishes. Sixth-order differentiation additionally retains $R'''$; it does not bypass the $R''$ contribution.

More generally,

$$
B_{m,r}\big|_{g_1=0}=0\quad\text{for }m<2r,
\qquad
B_{2p,p}\big|_{g_1=0}
=\frac{(2p)!}{2^pp!}g_2^p=(2p-1)!!\,g_2^p.
$$

The highest similarity derivative that can appear in an order-$2p$ physical derivative at a fold is therefore $R^{(p)}$, with contribution $(2p-1)!!R^{(p)}g_2^p$. Lower similarity derivatives generally remain as well.

| Physical derivative order | Highest surviving similarity derivative | Its Bell coefficient at the fold |
| --- | --- | --- |
| 2 | $R'$ | $g_2$ |
| 4 | $R''$ | $3g_2^2$ |
| 6 | $R'''$ | $15g_2^3$ |
| 8 | $R^{(4)}$ | $105g_2^4$ |
| 10 | $R^{(5)}$ | $945g_2^5$ |

These coefficients describe composition, not relative filter strength. The change from 3 to 15 does not establish a fivefold increase in sensitivity: the derivative orders, coordinate powers, coefficient units, and prescribed damping rates all differ. Terms can vanish or cancel, so survival in the Bell formula does not guarantee a nonzero net tendency.

The finite Bell matrix is triangular with diagonal $g_1,g_1^2,\ldots,g_1^m$. It is singular at a fold. Forward evaluation remains meaningful, but recovering similarity derivatives by ordinary matrix inversion is not valid there and becomes poorly conditioned nearby.

## 7. Recommended Boundary-Layer Evaluation

The following sequence would give a defensible basis for deciding whether vertical hyperdiffusion improves a column model or an NWP configuration:

1. **Specify the target.** Identify the prognostic fields and numerical wavelengths requiring control. Keep horizontal and vertical filter choices separate and retain the physical turbulent transport formulation.
2. **Verify the operator in isolation.** Check constant-field preservation, weighted conservation, quadratic-norm decay, boundary behavior, and measured Fourier decay where periodic uniform-grid tests apply. Test fourth and sixth order explicitly for the correct signs.
3. **Match damping times.** Compare fourth and sixth order using the same reference discrete mode and $\tau_g$, rather than the same numerical coefficient. Document the coefficient units and the time-integration method.
4. **Exercise physical profiles.** Use smooth manufactured profiles, a sharp stable inversion, and a low-level jet, then a relevant stable-boundary-layer case such as GABLS. Inspect temperature, vector shear, TKE, turbulent fluxes, and diagnosed $Ri_g$ alongside the filter tendencies.
5. **Refine the vertical grid.** Check inversion thickness, jet height and speed, surface fluxes, and the depth of turbulent mixing. Separate operator-verification tests at fixed coefficients from production tests with grid-scaled coefficients; both answer useful but different questions.
6. **Audit numerical budgets.** Report filter-induced kinetic-energy loss and scalar-norm changes separately from turbulent tendencies. If kinetic-energy loss is not returned to internal energy, identify it as a numerical energy sink.
7. **Assess positivity and boundaries.** Verify that the filter does not create unacceptable negative TKE or moisture, alter imposed surface fluxes, or generate reflected oscillations at the model top.

For fold diagnostics, record $g_1$, $g_2$, and the individual Bell contributions alongside the physical profiles. A fold in $\zeta(z)$ and a jet maximum are distinct diagnostics. Neither proves turbulence extinction. Any proposed convergence of TKE minima, diffusivity reduction, and derivative extrema should be tested under grid refinement rather than assumed.

## 8. Recommendation to Dick

Biharmonic filtering is a useful first candidate because it provides greater scale selectivity than Laplacian diffusion without the additional derivative order and boundary complexity of a sixth-order operator. Triharmonic filtering is worth testing when fourth-order damping measurably erodes resolved inversion or jet structure at the damping rate needed to control grid noise. Eighth and higher orders offer still sharper spectral selectivity, but also wider effective stencils, more demanding boundary treatment, and potentially troublesome oscillatory or positivity behavior.

The Bell-polynomial results provide a sound explanation of how these operators sample a similarity-based profile near a coordinate fold. They should be presented as diagnostic identities, not as evidence that high-order numerical dissipation represents turbulence physics or automatically preserves real atmospheric features. The practical decision should rest on matched-damping experiments, discrete budget checks, and vertical-resolution convergence.

## Background Reading

- Durran, D. R. (2010). *Numerical Methods for Fluid Dynamics: With Applications to Geophysics*, second edition. Springer. General background on numerical dissipation and time integration.
- Boyd, J. P. (2001). *Chebyshev and Fourier Spectral Methods*, second edition. Dover. Background on spectral filtering and wavenumber selectivity.
- Beare, R. J., et al. (2006). An intercomparison of large-eddy simulations of the stable boundary layer. *Boundary-Layer Meteorology*, **118**, 247-272. Stable-boundary-layer benchmark context, not a validation of the proposed filters.
- Local working notes: [Scaling the Sky](SkyScaling.md) and [7 October 2026](7Oct2026.md). The present report corrects their sixth-order sign convention and qualifies their pointwise-damping and fold interpretations. No particular WRF, MPAS, or ICON implementation is asserted here; those comparisons require checking the relevant model version and configuration.
