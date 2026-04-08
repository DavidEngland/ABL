Here’s a clean, McNider‑ready explanation that keeps the physics front‑and‑center and avoids unnecessary numerics, while still showing exactly what your code is doing.

---

🌄 How to Explain the Vertical Grid to McNider

Takeaway:

You’re using a computational coordinate η that is uniform, and mapping it into a physical height z using an exponential stretching function. This gives fine resolution near the surface while keeping the solver in physical z‑space. The Jacobian \(J = \frac{dz}{d\eta}\) is stored so that, if you ever move to a terrain‑following coordinate system, the metric terms are already available.

---

1. What the grid is

The model does not solve equations in η‑coordinates.
It solves everything directly in physical height z.

But instead of specifying z‑levels manually, you generate them by stretching a uniform computational coordinate:

\eta \in [0,1],\quad \eta_i = \frac{i}{N-1}


Then you map η → z using an exponential function:

z(\eta) = z_{\text{top}} \frac{e^{s\eta} - 1}{e^s - 1}, \qquad s = \text{stretch}


• \(s > 0\) → finer near the ground
• \(s = 0\) → uniform grid
• \(s \to \infty\) → hyper‑fine near the surface


This is the same logic as many boundary‑layer LES and SCM grids: concentrate resolution where gradients and turbulence live.

---

2. Why this matters physically

Near the surface, the ABL has:

• sharp temperature inversions
• strong shear
• thin stable layers
• intermittent turbulence


So you want small Δz near the ground and larger Δz aloft.
The exponential map gives exactly that.

---

3. The Jacobian: the key geometric object

Even though the solver currently works in physical z, the mapping has a well‑defined derivative:

J(\eta) = \frac{dz}{d\eta}
       = z_{\text{top}} \frac{s e^{s\eta}}{e^s - 1}


This is stored as grid.jacobian.

Why store it?

Because if you ever switch to a terrain‑following coordinate (η‑coordinate) formulation, the operators transform as:

• [ \frac{\partial}{\partial z} = \frac{1}{J} \frac{\partial}{\partial \eta} ]
• [ \nabla \cdot F = \frac{1}{J} \frac{\partial}{\partial \eta}(J F) ]


Those are the standard metric terms in any stretched vertical coordinate system.

So the grid is “simple” in z, but the metric is already computed for future use.

---

4. What the code computes

Your function does three things:

(1) Generate η

eta = range(0, 1, length=nz)


(2) Map η → z using the exponential stretching

z = z_top * (exp.(s*eta) .- 1) / (exp(s) - 1)


(3) Compute the Jacobian analytically

jac = z_top * s * exp.(s*eta) / (exp(s) - 1)


(4) Compute Δz for finite‑difference operators

Centered differences in the interior, one‑sided at boundaries.

---

5. How to say it in one sentence to McNider

“We use a uniform computational coordinate η and map it to physical height z with an exponential stretching function. This gives fine resolution near the surface. The Jacobian \(dz/d\eta\) is stored so that, if we ever move to a terrain‑following coordinate system, the metric terms are already available. But the solver itself still works directly in physical z.”

---

6. If you want an even shorter version

“It’s a stretched vertical grid: uniform in η, exponentially mapped to z. We store the Jacobian so the grid can be used either in physical‑z mode (current) or η‑coordinate mode (future terrain‑following).”

---

prepare a one‑slide diagram, a whiteboard‑style derivation, or a McNider‑friendly analogy (e.g., “η is the ruler, z is the curved tape measure”).

Compacted conversation

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#1-1), lines 1 to 100

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#300-300), lines 300 to 400

Here is a draft email section you can drop straight in:

---

**Subject: Proposed Vertical Grid Strategy for the SCM — Exponential Stretching**

---

Dick,

Following our conversation about the vertical grid, I wanted to write up the approach I am proposing to try. Happy to talk through any of it.

---

### The Problem: Uniform Grids Waste Resolution

In a uniform $\Delta z$ grid with $N = 40$ levels and $z_\text{top} = 1000\,\text{m}$, the lowest level sits at $z_1 = 25\,\text{m}$. That is already above the nocturnal stable layer and well above the roughness sublayer where MOST is defined. We waste most of our vertical degrees of freedom in the free atmosphere, where gradients are weak, while under-resolving the surface layer, where sensible heat, moisture, and momentum fluxes are concentrated.

---

### The Solution: Exponential Stretching

Instead of specifying $z$-levels manually, I generate them from a _uniform_ computational coordinate $\eta \in [0,1]$:

$$
\eta_i = \frac{i-1}{N-1}, \qquad i = 1, \ldots, N
$$

and then map $\eta \to z$ with an exponential function:

$$
z(\eta) = z_\text{top} \cdot \frac{e^{s\,\eta} - 1}{e^s - 1}
$$

where $s > 0$ is a **stretch parameter** (I am currently using $s = 3$).

This single formula does three things:

1. $\eta = 0 \Rightarrow z = 0$ (surface anchored exactly)
2. $\eta = 1 \Rightarrow z = z_\text{top}$ (top anchored exactly)
3. The spacing between levels _increases exponentially with height_

**Concrete example** with $N = 40$, $z_\text{top} = 1000\,\text{m}$, $s = 3$:

| Level | $\eta_i$ | $z_i$ (m) |
|-------|----------|-----------|
| 1 | 0.000 | 0.0 |
| 2 | 0.026 | 1.0 |
| 3 | 0.051 | 2.3 |
| 5 | 0.103 | 6.4 |
| 10 | 0.231 | 28 |
| 20 | 0.487 | 190 |
| 40 | 1.000 | 1000 |

The bottom 10% of levels ($z < 30\,\text{m}$) covers the surface layer and roughness sublayer properly. The upper half of levels distributes across the entrainment zone and free atmosphere where exact spacing matters less.

---

### The Jacobian: Why I Store $dz/d\eta$

The exponential map has an analytically known derivative:

$$
J(\eta) = \frac{dz}{d\eta} = z_\text{top} \cdot \frac{s\, e^{s\eta}}{e^s - 1}
$$

This is the **Jacobian of the coordinate transformation** — the local ratio of physical to computational spacing. It is not needed by the current solver (which works entirely in physical $z$), but I store it in the grid struct for two reasons:

**Reason 1: Terrain-following upgrade.**
If we ever want to run over sloping terrain or with a variable surface height, we switch to an $\eta$-coordinate formulation. The differential operators transform as:

$$
\frac{\partial \phi}{\partial z} = \frac{1}{J}\frac{\partial \phi}{\partial \eta}, \qquad
\nabla \cdot \mathbf{F} = \frac{1}{J}\frac{\partial}{\partial \eta}(J\,F_\eta)
$$

Those metric terms are just $1/J$ and $J$. They are already in `grid.jacobian`. No grid recomputation is needed when we make that upgrade.

**Reason 2: Mathematical structure.**
For small stretch $s$, the Jacobian has a Taylor expansion in terms of Bernoulli polynomials $B_n(\eta)$:

$$
\frac{J(\eta)}{z_\text{top}} = 1 + B_1(\eta)\,s + \frac{B_2(\eta)}{2!}\,s^2 + \cdots
$$

This is actually the generating function identity $e^{s\eta}/(e^s - 1) = \sum_n B_n(\eta)\,s^{n-1}/n!$ — a known result from analytic combinatorics. In practice this means: for small $s$ the grid is near-uniform with a slight linear tilt ($B_1(\eta) = \eta - \frac{1}{2}$); for large $s$ the near-surface clustering is controlled by the exponential. There is a continuous family of grids parameterized by $s$ that we can tune without changing any other code.

---

### What Stays the Same

The finite-difference operators (centered differences in the interior, one-sided at the boundaries), the tridiagonal implicit solver, and the closure interface — all of these use the physical $z$ array and $\Delta z$ spacings directly. Changing $s$ from 0 to 3 to 5 does not require touching any physics code, only the `create_grid` call.

---

### Proposed Next Step

I would like to run a short sensitivity test:

- $s = 0$ (uniform)
- $s = 2$ (moderate stretching)
- $s = 4$ (strong near-surface clustering)

with the stable nocturnal case (surface cooling −1 K/hr, 8 m/s geostrophic) and compare BL height estimates and surface flux convergence. This should show whether the stretching actually improves the near-surface gradient representation at modest $N$.

Let me know if you would like to discuss or if you want me to run those tests first.

Best,
David

---

That gives McNider (1) a clear statement of the problem, (2) the math in one place, (3) a concrete resolution table, (4) the Jacobian rationale at two levels of depth, and (5) a specific next step to evaluate rather than just "try it." Let me know if you want any section expanded or shortened.
