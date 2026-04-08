Vertical Grid Design Note (for Dick)

Takeaway:
The SCM uses a cell-centered, exponentially stretched vertical grid. The mapping places the first prognostic level a few meters above the surface (not at z = 0), which improves surface-layer physics while preserving a compact and stable solver in physical z-space.

---

1. Grid definition in index form

Use a uniform computational coordinate eta in [0, 1], but place prognostic nodes at cell centers:

$$
\eta_k = \frac{k - \tfrac{1}{2}}{N} = \frac{2k - 1}{2N}, \qquad k = 1, \ldots, N
$$

Map eta to physical height with exponential stretching:

$$
z(\eta) = z_{\mathrm{top}}\,\frac{e^{s\eta} - 1}{e^s - 1}
$$

Equivalent closed form in k:

$$
z_k = z_{\mathrm{top}}\,\frac{\exp\!\left(s\,\frac{2k-1}{2N}\right)-1}{e^s-1}
$$

with stretch factor s > 0.  The boundaries z = 0 and z = z_top are interfaces used for boundary conditions, not prognostic nodes.

For N = 40, z_top = 1000 m, s = 3:

- k = 1 -> z_1 approx 2 m
- k = 2 -> z_2 approx 6 m
- k = 3 -> z_3 approx 11 m

This is the practical reason the grid is useful for stable ABL runs: the first model level sits where MOST assumptions are still credible.

---

2. Surface-layer advantage

Placing z_1 around 2 m is strategic. In stable nocturnal conditions, the constant-flux layer is shallow; putting the first prognostic level too high (for example 10-15 m on a coarse uniform grid) can bias diagnosed fluxes and over-smooth near-surface gradients. Exponential stretching concentrates degrees of freedom where shear and scalar curvature are strongest and where closure sensitivity is highest.

---

3. Jacobian and metric terms

The analytical Jacobian of the mapping is

$$
J(\eta) = \frac{dz}{d\eta} = z_{\mathrm{top}}\,\frac{s e^{s\eta}}{e^s-1}
$$

Storing J in the grid object is useful for a future terrain-following eta-formulation because transformed operators are metric-weighted:

$$
\frac{\partial}{\partial z} = \frac{1}{J}\frac{\partial}{\partial \eta},
\qquad
\nabla\!\cdot F = \frac{1}{J}\frac{\partial}{\partial \eta}(J F)
$$

Implementation nuance:

- For coordinate-transform algebra and diagnostics, analytical J is ideal.
- For strict finite-volume conservation in a discrete solver, a cell-wise geometric ratio Delta z / Delta eta can sometimes be safer than pointwise J at centers.

Both are consistent if used carefully; the key is to keep operator/discretization pairing consistent so 1/J factors do not drift during derivation.

---

4. Bernoulli generating-function connection

The factor e^{s\eta}/(e^s - 1) is the Bernoulli-polynomial generating structure, so stretching can be interpreted as a formal perturbation away from a uniform grid:

$$
\frac{J(\eta)}{z_{\mathrm{top}}}
= 1 + B_1(\eta)s + \frac{B_2(\eta)}{2!}s^2 + \frac{B_3(\eta)}{3!}s^3 + \cdots
$$

with B_1(eta) = eta - 1/2. Interpretation:

- Small s: near-uniform grid with a linear first correction, so second-order central differencing remains close to second order in physical space.
- Larger s: higher Bernoulli terms become active, and metric variation must be respected more carefully in truncation-error analysis.

This gives a clean bridge from analytic structure to numerical behavior.

---

5. Small-s implementation safeguard

Direct evaluation of

$$
\frac{e^{s\eta}-1}{e^s-1}
$$

has a 0/0 form as s -> 0 and can lose precision when |s| is very small. A robust implementation uses a branch:

- if |s| < s_tol: use Taylor form (uniform-grid limit)
- else: use exponential form

A numerically stable option is to use expm1:

$$
z(\eta) = z_{\mathrm{top}}\,\frac{\operatorname{expm1}(s\eta)}{\operatorname{expm1}(s)}
$$

which avoids catastrophic cancellation near s = 0.

---

6. Top boundary with stretched grids

Current SCM behavior is generally well-behaved at the top because the present dynamics are diffusion/relaxation-dominated rather than fully wave-resolving compressible dynamics. In that regime, large upper-level Delta z mostly reduces resolution aloft rather than creating strong wave reflections.

That said, as the framework evolves, reflections can appear if fast modes are introduced (for example, explicit gravity-wave-supporting dynamics) and the top boundary remains purely reflective or under-damped. If needed, mitigation is straightforward:

- add a top sponge/relaxation layer in the upper 10-20% of the column,
- enforce radiative/open-type top conditions for wave-supporting variables,
- and check sensitivity to z_top and top-layer damping timescale.

---

One-line summary:
Uniform eta, exponential eta->z mapping, cell-centered nodes, and stored Jacobian together provide a strong near-surface representation now and a clean path to terrain-following metrics later.
