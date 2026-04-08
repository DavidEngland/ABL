Dick & Arastoo,

Following our conversation about the vertical grid, I wanted to write up the approach I am proposing to try. Happy to talk through any of it.

### Exponential Stretching

Instead of specifying $z$-levels manually, I use a **cell-centred** layout: $N$ interface gaps span $\eta \in [0,1]$ uniformly, and prognostic levels sit at the cell centres

$$
\eta_i = \frac{i - \tfrac{1}{2}}{N}, \qquad i = 1, \ldots, N
$$

and then map $\eta \to z$ with an exponential function:

$$
z(\eta) = z_\text{top} \cdot \frac{e^{s\,\eta} - 1}{e^s - 1}
$$

where $s > 0$ is a **stretch parameter** (I am currently using $s = 3$).

Substituting the cell-centre $\eta_k$ directly gives the closed-form height of level $k$:

$$
\boxed{z_k = z_\text{top} \cdot \frac{\exp\!\left(\dfrac{s(2k-1)}{2N}\right) - 1}{e^s - 1}}, \qquad k = 1, \ldots, N
$$

For $k=1$ this gives $z_1 = z_\text{top}\,(e^{s/(2N)}-1)/(e^s-1)$, which for large $s$ scales as $z_\text{top}\,e^{s/(2N)}/e^s \approx z_\text{top}\,e^{-s(1-1/(2N))}$ — exponentially small relative to $z_\text{top}$, exactly the surface-layer clustering we want.

Key properties:

1. The surface $z = 0$ and the domain top $z = z_\text{top}$ are **interface** (boundary-condition) levels — not prognostic nodes.
2. The lowest prognostic level $z_1 > 0$ is above the ground.
3. The spacing between levels _increases exponentially with height_.

**Concrete example** with $N = 40$, $z_\text{top} = 1000\,\text{m}$, $s = 3$ (cell-centred, $\eta_i = (i-\tfrac{1}{2})/40$):

| Level $i$ | $\eta_i$ | $z_i$ (m) |
|-----------|----------|-----------|
| 1 | 0.0125 | 2.0 |
| 2 | 0.0375 | 6.2 |
| 3 | 0.0625 | 10.8 |
| 5 | 0.1125 | 21 |
| 10 | 0.2375 | 54 |
| 20 | 0.4875 | 174 |
| 40 | 0.9875 | 962 |

The first prognostic level is at $z_1 \approx 2\,\text{m}$ — well within the surface layer where MOST is valid. Compare this with the uniform cell-centred grid where $z_1 = 12.5\,\text{m}$. The bottom six levels (up to $\sim 20\,\text{m}$) would compress into the first two levels of the uniform grid.

---

### The Jacobian: Store $dz/d\eta$?

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


Best,
David

---