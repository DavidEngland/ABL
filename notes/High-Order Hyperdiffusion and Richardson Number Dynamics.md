# High-Order Hyperdiffusion and Richardson Number Dynamics in Atmospheric Boundary-Layer Modeling

---

## Pillar I: Physical Closure Mapping ($\zeta \longleftrightarrow Ri \longleftrightarrow (S_m, S_h)$)

### 1. Theoretical Foundations: Gradient Richardson Number Dynamics and Monin–Obukhov Linkage

Modern Numerical Weather Prediction (NWP) models and atmospheric boundary-layer parameterizations require computationally efficient methods to evaluate turbulent exchange and thermal stratification. Historically, non-dimensional Monin–Obukhov Similarity Theory (MOST) relied on iterative sub-loops to determine the dimensionless stability parameter $\zeta = z/L$ (where $z$ is physical height and $L$ is the Obukhov length scale). In high-performance computing environments—specifically vector processing units and GPUs—iterative loops induce severe performance bottlenecks due to conditional branching and warp divergence. Eliminating sub-iterations by mapping MOST functions directly into physical Gradient Richardson Number ($Ri$) space provides a direct, non-iterative, closed-form mapping. This mathematical bridge preserves exact algebraic fidelity with MOST while enabling unified, single-pass forward evaluations across all stability regimes.

### 1.1 Core Formulations & Exact $S_m(Ri) / S_h(Ri)$ Mapping

The parameterization of first-order turbulent exchange coefficients for momentum ($K_m$) and heat or scalars ($K_h$) links the Monin–Obukhov similarity framework to local gradient formulations:

* **Monin–Obukhov Formulation ($\zeta = z/L$):**

$$K_m = \frac{l_m u_*}{\phi_m(\zeta)}, \qquad K_h = \frac{l_m u_*}{\phi_h(\zeta)}$$

* **Gradient Richardson Number Formulation ($Ri$):**

$$K_m = S_m(Ri) \, l_m^2 s, \qquad K_h = S_h(Ri) \, l_m^2 s$$

where $l_m = \kappa (z - d)$ represents the surface-layer mixing length, $u_*$ is the friction velocity, $s = |\partial \mathbf{V} / \partial z|$ is the magnitude of vertical wind shear, and $\phi_m(\zeta), \phi_h(\zeta)$ are the non-dimensional MOST gradient functions. The dimensionless exchange stability functions $S_m(Ri)$ and $S_h(Ri)$ encapsulate the damping effects of thermal stratification on mechanical turbulence.

Equating friction velocity $u_* = S_m^{1/2} l_m s$ and scalar scale $\theta_* = (S_h / S_m^{1/2}) l_m (\partial \theta_v / \partial z)$ within the analytical definition of the Obukhov length $L = u_*^2 / (\kappa B \theta_*)$—where $B = g / \theta_{v,\mathrm{ref}}$ is the buoyancy parameter—yields the exact dual transformation identity:

$$\zeta(Ri) = Ri \frac{S_h(Ri)}{S_m^{3/2}(Ri)}$$

The exponent of $3/2$ on $S_m(Ri)$ in the denominator is a strict physical and mathematical requirement. When the turbulent Prandtl number $Pr_t = K_m / K_h = S_m(Ri) / S_h(Ri) = \phi_h(\zeta) / \phi_m(\zeta)$ deviates from unity ($Pr_t \neq 1$), any exponent other than $3/2$ violates dimensional consistency and creates algebraic discrepancies between friction velocity and local gradient definitions.

The table below provides a structural comparison between the Monin–Obukhov ($\zeta$) and Gradient Richardson Number ($Ri$) formulations:

| Parameter / Operational Concept | Monin–Obukhov Formulation ($\zeta = z/L$) | Gradient Richardson Number Formulation ($Ri$) |
| --- | --- | --- |
| **Primary Variable** | Non-dimensional height $\zeta = z/L$ | Local Gradient Richardson Number $Ri = \frac{N^2}{s^2}$ |
| **Exchange Coefficients ($K_m, K_h$)** | $K_m = \frac{l_m u_*}{\phi_m(\zeta)}, \quad K_h = \frac{l_m u_*}{\phi_h(\zeta)}$ | $K_m = S_m(Ri) l_m^2 s, \quad K_h = S_h(Ri) l_m^2 s$ |
| **Stability Functions** | Non-dimensional gradients $\phi_m(\zeta), \phi_h(\zeta)$ | Exchange functions $S_m(Ri) = \phi_m^{-2}(\zeta), \quad S_h(Ri) = \phi_m^{-1}(\zeta)\phi_h^{-1}(\zeta)$ |
| **Turbulent Prandtl Number ($Pr_t$)** | $Pr_t(\zeta) = \frac{\phi_h(\zeta)}{\phi_m(\zeta)}$ | $Pr_t(Ri) = \frac{S_m(Ri)}{S_h(Ri)}$ |
| **Inversion Mappings** | Classical iterative loops ($L$-convergence) | Closed-form identity $\zeta(Ri) = Ri \frac{S_h(Ri)}{S_m^{3/2}(Ri)}$ |

---

### 1.2 Regime-Specific Analytical Inversions

#### Stable Boundary Layer ($SBL, Ri > 0$)

For linear non-dimensional MOST gradients under stable stratification, $\phi_m(\zeta) = 1 + \beta_m \zeta$ and $\phi_h(\zeta) = \alpha_\theta + \beta_h \zeta$ (where $\alpha_\theta = Pr_t(0)$ is the neutral turbulent Prandtl number), inverting the stability relations yields an exact quadratic equation for $\zeta(Ri)$:

$$\zeta(Ri) = \frac{\alpha_\theta - 2 Ri \beta_m - \sqrt{\alpha_\theta^2 + 4(\beta_h - \alpha_\theta \beta_m)Ri}}{2(Ri \beta_m^2 - \beta_h)}$$

A crucial mathematical constraint in this quadratic solution is the selection of the negative root before the radical. Selecting the positive root yields $\lim_{Ri \to 0^+} \zeta_{+}(Ri) = \frac{\alpha_\theta}{-\beta_h} \neq 0$, creating an unphysical step discontinuity at neutral stability and violating the foundational physical boundary condition $\zeta(0) = 0$. Selecting the negative root strictly enforces $\zeta(0) = 0$ through L'Hôpital's rule:

$$\lim_{Ri \to 0^+} \zeta_{-}(Ri) = \lim_{Ri \to 0^+} \frac{\frac{d}{dRi}\left[\alpha_\theta - 2 Ri \beta_m - \sqrt{\alpha_\theta^2 + 4(\beta_h - \alpha_\theta \beta_m)Ri}\right]}{\frac{d}{dRi}\left[2(Ri \beta_m^2 - \beta_h)\right]} = \frac{1}{\alpha_\theta}$$

Evaluating the asymptotic limit of this closed-form mapping as $\zeta \to \infty$ reveals the closure-implied asymptotic Richardson-number limit:

$$\lim_{\zeta \to \infty} Ri(\zeta) = \frac{\beta_h}{\beta_m^2}$$

This parameter represents the mathematical asymptote of the specified rational function mapping, rather than an intrinsically universal dynamical critical Richardson number. It serves as an extinction threshold ($Ri_c$) if and only if the underlying turbulence closure explicitly enforces $S_m = S_h = 0$ at or above that bound.

When equal diffusivities at neutrality are assumed ($\alpha_\theta = 1$) alongside matching empirical coefficients ($\beta_m = \beta_h = \beta$), the asymptotic limit becomes $Ri_c = 1/\beta$. In this equal-diffusivity limit ($Pr_0 = 1$), the exchange stability functions take the piecewise form:

$$S_m(Ri) = S_h(Ri) = \begin{cases} \left(1 - \frac{Ri}{Ri_c}\right)^2, & 0 \le Ri < Ri_c \\ 0, & Ri \ge Ri_c \end{cases}$$

Evaluating derivatives across $Ri \to Ri_c^-$ yields $S(Ri_c^-) = 0$ and $S'(Ri_c^-) = 0$. Thus, $S(Ri)$ exhibits $C^1$ continuity, suppressing derivative jumps across the threshold. However, $S(Ri)$ is not $C^2$ continuous because $S''(Ri_c^-) = \frac{2}{Ri_c^2} \neq 0 = S''(Ri_c^+)$. Furthermore, $C^1$ derivative continuity suppresses gradient steps but does not guarantee Jacobian non-singularity or dynamical stability in implicit time integration: as exchange diffusivities vanish ($K_m, K_h \to 0$), the physical tangent operator can still become singular or ill-conditioned.

Modern Very Stable Boundary Layer (VSBL) research (e.g., Gryanik et al., 2020; Energy and Flux Budget [EFB] theory) addresses this by replacing zero-extinction limits with non-zero asymptotic tails ($\phi_m, \phi_h \sim \zeta^1$ as $\zeta \to \infty$), preventing nocturnal runaway surface cooling and physical land-surface decoupling.

#### Unstable Boundary Layer ($UBL, Ri < 0$)

Under unstable conditions, adopting standard Businger–Dyer empirical exponents ($\alpha_m = -1/4, \alpha_h = -1/2$) maps the relation between $\zeta$ and $Ri$ into a Master Cubic polynomial in $x \equiv \zeta$:

$$b_m \zeta^3 - \zeta^2 - b_h \left(\frac{Ri}{\alpha_\theta}\right)^2 \zeta + \left(\frac{Ri}{\alpha_\theta}\right)^2 = 0$$

When the momentum and scalar empirical exponents match identically ($b_m = b_h = B_u \approx 16$ and $\alpha_\theta = 1.0$), the stability functions evaluate to $S_m(Ri) = (1 - B_u \zeta)^{1/2}$ and $S_h(Ri) = (1 - B_u \zeta)^{3/4} = S_m^{3/2}(Ri)$. Substituting these into the dual transformation identity yields:

$$\zeta = Ri \frac{S_h(Ri)}{S_m^{3/2}(Ri)} = Ri \frac{(1 - B_u \zeta)^{3/4}}{(1 - B_u \zeta)^{3/4}} \equiv Ri$$

This "Dyer collapse" ($\zeta \equiv Ri$) is a closure-dependent coordinate collapse specific to matching exponents, rather than a universal physical equivalence.

For non-equal exponents ($b_m \neq b_h$), series inversion must be derived directly from the underlying rational mapping or master cubic to avoid precision loss near neutral stability ($\vert{}Ri\vert{} \ll 1$). For example, Taylor series expansion of the stable rational mapping $Ri(\zeta) = \frac{\zeta(\alpha_\theta + \beta_h \zeta)}{(1 + \beta_m \zeta)^2}$ yields the exact inversion:

$$\zeta(Ri) = \frac{Ri}{\alpha_\theta} + \frac{2\alpha_\theta\beta_m - \beta_h}{\alpha_\theta^3} Ri^2 + \frac{2(2\alpha_\theta\beta_m - \beta_h)^2 - \alpha_\theta^2(3\beta_m^2 - 2\beta_m\beta_h/\alpha_\theta)}{\alpha_\theta^5} Ri^3 + \mathcal{O}(Ri^4)$$

For unstable profiles with non-equal Businger–Dyer coefficients, the expansion coefficients are obtained by expanding the corresponding master cubic around $Ri = 0^-$.

#### Numerical Regularization & GPU Alignment

To maximize computational throughput on hardware architectures sensitive to thread divergence, conditional logic (`if`/`else` branching across $Ri = 0$) must be removed. Zero-Offset Hyperbolic Regularization (Z0HR) replaces piecewise logic by decomposing $Ri$ into $C^\infty$-smooth positive (stable) and negative (unstable) projections using an $\epsilon$-regularized absolute value:

$$\vert{}Ri\vert{}_\epsilon = \sqrt{Ri^2 + \epsilon^2}, \qquad Ri_+ = \frac{1}{2}\left(Ri + \vert{}Ri\vert{}_\epsilon - \epsilon\right), \qquad Ri_- = \frac{1}{2}\left(Ri - \vert{}Ri\vert{}_\epsilon + \epsilon\right)$$

Smooth weight functions $w_{\text{stable}} = \frac{1}{2}\left(1 + \frac{Ri}{\vert{}Ri\vert{}_\epsilon}\right)$ and continuous hyperbolic tangent embeddings:

$$S_m(Ri) = \frac{1}{2}\left[1 - \tanh(\sigma Ri)\right] S_{m,\text{uns}} + \frac{1}{2}\left[1 + \tanh(\sigma Ri)\right] S_{m,\text{sta}}$$

allow unified, branch-free execution. Removing conditional branching eliminates GPU warp divergence and CPU instruction pipeline stalls, accelerating execution across mass-parallel grid cells.

---

## Pillar II: Geometric Derivative Mapping ($\mathbf{R}_\zeta \xrightarrow{\mathbf{B}} \mathbf{R}_z$) & Kinematics

### 2. Higher-Order Spatial Derivatives and Partial Bell Matrix Mappings

Evaluating higher-order vertical spatial derivatives of the Gradient Richardson Number ($Ri_z, Ri_{zz}, Ri_{zzz}, Ri_{zzzz}$) is essential for high-order spatial discretizations, sub-grid scale diffusion models, and turbulence stability analyses. When Monin–Obukhov stability profiles $R(\zeta)$ are transformed into physical vertical coordinates $z$, the derivative chains generate non-linear couplings between intrinsic stability curvature and geometric coordinate variations. A mathematically rigorous mapping framework is required to transform non-dimensional derivative vectors into physical spatial derivative vectors without truncation errors.

### 2.1 Faà di Bruno’s Chain Rule and the Partial Bell Matrix ($\mathbf{B}$)

Differentiating a composite function $Ri(z) = R(\zeta(z))$ $m$ times with respect to physical height $z$ is governed by Faà di Bruno’s higher-order chain rule formula:

$$\frac{d^m}{dz^m} R(\zeta(z)) = \sum_{n=1}^{m} R^{(n)}(\zeta) \cdot B_{m,n}\left(\zeta_z, \zeta_{zz}, \dots, \zeta_{(m-n+1)z}\right)$$

where $B_{m,n}$ represents the partial exponential Bell polynomial (Bell polynomial of the second kind), which combinatorially counts the partitions of a set of $m$ elements into $n$ non-empty subsets.

Let $g_j \equiv \frac{d^j \zeta}{d z^j}$ denote spatial coordinate derivatives (where $g_1 = \zeta_z, g_2 = \zeta_{zz}$, etc.). Collecting derivative components up to fourth order constructs the lower-triangular $4 \times 4$ Partial Bell Transformation Matrix system $\mathbf{Ri}_z = \mathbf{B} \mathbf{R}_\zeta$:

$$\begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix} = \begin{bmatrix} g_1 & 0 & 0 & 0 \\ g_2 & g_1^2 & 0 & 0 \\ g_3 & 3g_1g_2 & g_1^3 & 0 \\ g_4 & 4g_1g_3 + 3g_2^2 & 6g_1^2g_2 & g_1^4 \end{bmatrix} \begin{bmatrix} R'(\zeta) \\ R''(\zeta) \\ R'''(\zeta) \\ R^{(4)}(\zeta) \end{bmatrix}$$

The individual matrix entries $B_{m,n}$ can be assembled using two equivalent recurrence relations:

* **Combinatorial Summation Recurrence (Faà di Bruno Identity):**

$$B_{m+1, n} = \sum_{j=n-1}^{m} \binom{m}{j} g_{m+1-j} B_{j, n-1} \quad \text{with } B_{0,0} = 1$$

* **Differential Recurrence:**

$$B_{m+1, n} = g_1 B_{m, n-1} + \frac{d}{dz} B_{m, n} \quad \text{with boundary conditions } B_{m,0} = 0 \text{ for } m \ge 1 \text{ and } B_{m,k} = 0 \text{ for } k > m$$

---

### 2.2 Inverse Bell Mapping ($\mathbf{B}^{-1}$) and Observational Parameter Calibration

Because $\mathbf{B}$ is lower-triangular with diagonal entries $B_{k,k} = g_1^k$, its determinant is:

$$\det(\mathbf{B}) = g_1^{1+2+3+4} = g_1^{10}$$

When $g_1 \neq 0$, inverting this matrix maps physical spatial derivatives back into non-dimensional similarity derivatives ($\mathbf{R}_\zeta = \mathbf{B}^{-1} \mathbf{Ri}_z$):

$$\begin{bmatrix} R'(\zeta) \\ R''(\zeta) \\ R'''(\zeta) \\ R^{(4)}(\zeta) \end{bmatrix} = \begin{bmatrix} \frac{1}{g_1} & 0 & 0 & 0 \\ -\frac{g_2}{g_1^3} & \frac{1}{g_1^2} & 0 & 0 \\ \frac{3g_2^2 - g_1 g_3}{g_1^5} & -\frac{3g_2}{g_1^4} & \frac{1}{g_1^3} & 0 \\ \frac{10g_1 g_2 g_3 - 15g_2^3 - g_1^2 g_4}{g_1^7} & \frac{10g_2^2 - 4g_1 g_3}{g_1^6} & -\frac{6g_2}{g_1^5} & \frac{1}{g_1^4} \end{bmatrix} \begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix}$$

Processing high-resolution field tower or Large-Eddy Simulation (LES) profile data through $\mathbf{B}^{-1}$ extracts non-dimensional similarity derivatives ($R', R'', R''', R^{(4)}$) *a priori* directly from measured physical profiles without prescribing functional forms. Near neutral stability ($\zeta \to 0$), empirical expansion slopes ($\beta_m, \beta_h$) and the neutral Prandtl number $\alpha_\theta = Pr_t(0)$ are diagnosed via:

$$\alpha_\theta = \left. \frac{Ri_z}{g_1} \right\vert{}_{\zeta \to 0}, \qquad 2(\beta_h - 2\alpha_\theta \beta_m) = \left. \frac{Ri_{zz} - \frac{g_2}{g_1} Ri_z}{g_1^2} \right\vert{}_{\zeta \to 0}$$

At or near a coordinate fold knee ($g_1 \to 0$), the transformation matrix becomes rank-deficient ($\det(\mathbf{B}) \to 0$). Consequently, observational parameter calibration routines cannot directly invert field profile data via $\mathbf{B}^{-1}$ at a knee without regularization (e.g., Tikhonov regularization $\mathbf{D}_\lambda^{(2)} = (\mathbf{A}^T \mathbf{A} + \lambda \mathbf{R})^{-1} \mathbf{A}^T$ or local chart transitions).

---

### 2.3 Asymptotic Limits and Profile Curvature

#### Near-Neutral Limit ($\zeta \to 0$)

Evaluating non-dimensional profile curvature $R''(0)$ near neutral stability reveals a key mathematical distinction between single-coefficient and dual-coefficient MOST closures:

* **Single-Coefficient Closure** ($R(\zeta) = \frac{\zeta}{1 + \beta \zeta}$, $\beta = 5.0$):

$$R''(0) = -2\beta = -10.0$$

* **Dual-Coefficient Closure** ($R(\zeta) = \frac{\zeta(\alpha_\theta + \beta_h \zeta)}{(1 + \beta_m \zeta)^2}$, $\alpha_\theta = 1.0, \beta_m = 6.0, \beta_h = 5.0$):

$$R''(0) = 2(\beta_h - 2\alpha_\theta \beta_m) = 2(5.0 - 12.0) = -14.0$$

The dual-coefficient model produces a ~40% steeper initial curvature spike near neutrality ($R''(0) = -14.0$ vs. $-10.0$). On discrete C-grid vertical meshes, this steeper initial curvature accelerates spatial gradient variations near the surface:

$$\lim_{\zeta \to 0} Ri_{zzzz} = g_4 + R''(0) (4g_1 g_3 + 3g_2^2) + R'''(0) (6g_1^2 g_2) + R^{(4)}(0) g_1^4$$

This provides the mathematical justification for mandatory softening noise floor parameters ($\epsilon_c Ri_c > 0$) in surface-layer modules to bound spatial derivative amplification near neutrality without corrupting physical turbulent fluxes.

#### High-Stability Limit ($\zeta \gg 1$)

Under strongly stable conditions ($\zeta \to \infty$), higher-order non-dimensional derivatives exhibit power-law decay:

$$R^{(n)}(\zeta) \approx (-1)^n n! \, \frac{\alpha_\theta \beta_m - 2\beta_h}{\beta_m^3} \zeta^{-(n+1)} \sim \mathcal{O}\left(\zeta^{-(n+1)}\right)$$

Assuming a locally constant Obukhov length $L$, physical vertical derivatives scale as $\frac{d^n Ri}{dz^n} \propto \frac{L}{z^{n+1}}$. For fourth-order derivatives, $Ri_{zzzz} \propto z^{-5}$. This rapid decay proves that standard MOST closures naturally suppress numerical noise aloft. Consequently, persistent grid-scale oscillations observed aloft in strongly stable layers during model execution indicate non-MOST physical processes (such as gravity wave breaking or shear instability) or boundary reflection, rather than inherent MOST model instability.

---

### 2.4 Kinematics: Jet "Nose" vs. Coordinate "Knee" and Fold Geometry

Generalized Similarity Profile Theory (GSPT) projects continuous state-space dynamics onto a 1D non-dimensional similarity coordinate $\zeta(z) = z / L(z)$. Distinguishing physical profile kinematics from mathematical coordinate artifacts is essential:

* **The Low-Level Jet (LLJ) "Nose":** A physical fluid-dynamical momentum maximum where total vertical vector shear magnitude $s = \vert{}\partial \mathbf{V} / \partial z\vert{} \to 0$ (or reaches a local minimum). Because shear appears in the denominator of $Ri = N^2/s^2$, vanishing vector shear causes point-gradient $Ri$ values to diverge toward $+\infty$. Note that a scalar wind speed maximum does not automatically imply $U_z = V_z = 0$ if directional shear is present.
* **The Coordinate "Knee" (Coordinate Fold Singularity):** A mathematical projection artifact where local flux divergence balances similarity scaling ($\zeta L' = 1$), causing the coordinate derivative to vanish ($g_1 = \zeta_z = 0$) and destroying local coordinate invertibility.

Total differentiation of $Ri(\zeta(z))$ yields the exact profile curvature decomposition:

$$\frac{d^2 Ri}{dz^2} = \underbrace{R''(\zeta) g_1^2}_{\text{Mapped Intrinsic Curvature Contribution}} + \underbrace{R'(\zeta) g_2}_{\text{Coordinate/Flux Geometry}} + \underbrace{\mathcal{E}_{\Delta z}}_{\text{Estimation Error}}$$

At a coordinate fold singularity ($g_1 = 0, g_2 \neq 0$), the Partial Bell Transformation Matrix evaluates to:

$$\mathbf{B}\Big\vert{}_{g_1 = 0} = \begin{bmatrix} 0 & 0 & 0 & 0 \\ g_2 & 0 & 0 & 0 \\ g_3 & 0 & 0 & 0 \\ g_4 & 3g_2^2 & 0 & 0 \end{bmatrix} \implies \begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix}_{g_1 = 0} = \begin{bmatrix} 0 \\ R'(\zeta) g_2 \\ R'(\zeta) g_3 \\ R'(\zeta) g_4 + 3 R''(\zeta) g_2^2 \end{bmatrix}$$

At $g_1 = 0$, physical gradient $Ri_z = 0$. The contribution of intrinsic similarity curvature $R''(\zeta)$ to second derivative $Ri_{zz}$ vanishes specifically because it is multiplied by $g_1^2 = 0$; the intrinsic similarity curvature $R''(\zeta)$ itself remains non-zero. Physical profile bending $Ri_{zz}$ at the fold reduces purely to coordinate geometry ($R'(\zeta) g_2$). In finite-neighborhood observational evaluations (e.g., CASES-99, GABLS3 field profiles), coordinate compression ($g_2 = \zeta_{zz}$) is the dominant driver of profile bending near coordinate knees.

At fourth order, the non-zero Bell entry $B_{4,2} = 3g_2^2$ retains the active intrinsic curvature term $3 R''(\zeta) g_2^2$. Higher-order combinatorial entries at $g_1 = 0$ up to tenth order follow the exact double-factorial relationship:

| Physical Derivative Order ($2p$) | Highest Surviving Similarity Derivative | Bell Polynomial Entry $B_{2p, p}\big\vert{}_{g_1=0}$ | Exact Combinatorial Entry |
| --- | --- | --- | --- |
| **2** | $R'(\zeta)$ | $g_2$ | $1 \cdot g_2$ |
| **4** | $R''(\zeta)$ | $3 g_2^2$ | $3 \cdot g_2^2 = (3!!) g_2^2$ |
| **6** | $R'''(\zeta)$ | $15 g_2^3$ | $15 \cdot g_2^3 = (5!!) g_2^3$ |
| **8** | $R^{(4)}(\zeta)$ | $105 g_2^4$ | $105 \cdot g_2^4 = (7!!) g_2^4$ |
| **10** | $R^{(5)}(\zeta)$ | $945 g_2^5$ | $945 \cdot g_2^5 = (9!!) g_2^5$ |

The entry $B_{2p, p}\big\vert{}_{g_1=0} = (2p-1)!! g_2^p$ reflects function composition geometry across coordinate transforms. It does not represent an intrinsic physical amplification of turbulence or numerical filtering strength.

---

### 2.5 The Triple-Point Convergence Hypothesis

To prevent misinterpreting coordinate fold artifacts or grid truncation errors as true physical turbulence collapse, GSPT establishes the Triple-Point Convergence Diagnostic. Physical turbulence extinction is evaluated by tracking the spatial convergence of three profile heights as vertical grid spacing vanishes ($\Delta z \to 0$):

1. **Diffusivity Extinction Height ($z_K$):** The level where exchange diffusivity reaches its background floor ($K_m \to K_{\min}$).
2. **TKE Floor Height ($z_e$):** The level where Turbulent Kinetic Energy reaches its minimum floor ($e \to e_{\min}$).
3. **Shear Stress / TKE Gradient Extremum Height ($z_{e_z}$):** The level where vertical stress divergence reaches a localized extremum ($\frac{\partial^2}{\partial z^2}(u_*^2) \to \text{max}$).

The collapse of spatial spread $\Delta z_{\text{TP}} = \max\vert{}z_i - z_j\vert{} \to 0$ under grid refinement ($\Delta z \to 0$) serves as a necessary numerical convergence diagnostic. To confirm true physical turbulence extinction, this spatial convergence must be paired with dynamical state criteria, such as background diffusivity floor convergence ($K_m \to K_{\min}$) or fast-mode Jacobian eigenvalue decay ($\lambda_f \to 0^+$).

---

## Pillar III: Numerical Filtering & Discretization ($q \xrightarrow{\mathcal{D}_{2p}} q^*$)

### 3. Mechanics of Biharmonic and Higher-Order Hyperdiffusion Operators

High-order hyperdiffusion operators provide scale-selective dissipation to attenuate high-wavenumber computational grid noise near the spatial truncation limit ($2\Delta z$) while preserving resolved boundary-layer structures.

### 3.1 Prognostic vs. Diagnostic Filtering Separation

A critical distinction must be enforced between prognostic time integration and diagnostic profile analysis:

$$\text{Prognostic State Variables } (u, v, \theta_v, q_v) \xrightarrow{\text{Numerical Hyperdiffusion Filter } \mathcal{D}_{2p}} (u^*, v^*, \theta_v^*, q_v^*) \xrightarrow{\text{Diagnostic Mapping}} Ri^*$$

Hyperdiffusion operators must be applied directly to prognostic state variables. Diagnostic profile operators (such as $\mathcal{D}_4(Ri)$) are strictly analytical profiling and diagnostic smoothing tools; they are not operational momentum or thermodynamic prognostic updates.

---

### 3.2 Scale-Selective Dissipation Mechanics & Sign Conventions

For a model prognostic field $q(\mathbf{x}, t)$, the order-$2p$ spatial filtering tendency is defined as:

$$\left.\frac{\partial q}{\partial t}\right\vert{}_{\mathrm{filter}} = -\nu_{2p}(-\Delta)^p q = (-1)^{p+1}\nu_{2p}\Delta^p q$$

where $\Delta = \nabla^2$ is the physical Laplacian operator, $p \ge 1$ is an integer defining spatial derivative order $2p$, and $\nu_{2p} \ge 0$ is the hyperdiffusion coefficient.

Applying this continuum operator to a spatial Fourier mode $\widehat{q} e^{i k z}$ yields the evolution equation $\frac{d\widehat{q}}{dt} = -\sigma_{2p}(k) \widehat{q}$, where the operator eigenvalue is $-\sigma_{2p}(k) = -\nu_{2p} k^{2p}$ and the positive damping rate magnitude is $\sigma_{2p}(k) = \nu_{2p} k^{2p}$.

| Spatial Derivative Order ($2p$) | Common Operator Name | Dissipative Tendency Expression | Positive Fourier Damping Rate ($\sigma_{2p}(k)$) |
| --- | --- | --- | --- |
| **2** | Laplacian Diffusion | $+\nu_2 \Delta q$ | $\nu_2 k^2$ |
| **4** | Biharmonic Hyperdiffusion | $-\nu_4 \Delta^2 q$ | $\nu_4 k^4$ |
| **6** | Triharmonic Hyperdiffusion | $+\nu_6 \Delta^3 q$ | $\nu_6 k^6$ |
| **8** | Eighth-Order Hyperdiffusion | $-\nu_8 \Delta^4 q$ | $\nu_8 k^8$ |
| **$2p$** | General Even-Order | $(-1)^{p+1}\nu_{2p} \Delta^p q$ | $\nu_{2p} k^{2p}$ |

Calibrating the hyperdiffusion coefficient to achieve a target e-folding damping time $\tau_g$ at reference grid-scale wavenumber $k_g = \pi / h$ gives $\nu_{2p} = \frac{1}{\tau_g k_g^{2p}}$, yielding $\sigma_{2p}(k) = \frac{1}{\tau_g} \left(\frac{k}{k_g}\right)^{2p}$. Increasing $2p$ sharpens scale selectivity, damping $2\Delta z$ noise while preserving longer resolved waves ($4\Delta z$ damping ratio scales as $(0.5)^{2p}$).

On a discrete uniform grid using a centered 3-point Laplacian operator, continuum wavenumber calibration $k^2$ is modified by the discrete Laplacian eigenvalue $\mu(k) = \frac{4}{h^2} \sin^2\left(\frac{kh}{2}\right)$:

$$\sigma_{2p}(k) = \nu_{2p} [\mu(k)]^p$$

At the Nyquist limit ($kh = \pi$), $\mu_g = 4/h^2$, yielding the discrete matched coefficient $\nu_{2p,\text{discrete}} = \frac{1}{\tau_g \mu_g^p} = \frac{h^{2p}}{4^p \tau_g}$. For fourth-order biharmonic hyperdiffusion ($p=2$):

$$\nu_{4,\mathrm{discrete}} = \frac{h^4}{16 \tau_g}$$

Comparing the continuum calibration formula ($\nu_{4,\mathrm{continuum}} = \frac{h^4}{\pi^4 \tau_g}$) to the discrete eigenvalue calibration yields:

$$\frac{\nu_{4,\mathrm{continuum}}}{\nu_{4,\mathrm{discrete}}} = \frac{\pi^4}{16} \approx 6.088$$

This factor of $6.09$ specifically quantifies the over-damping of the shortest grid scale ($2\Delta z$) when the continuum-calibrated coefficient is evaluated directly in a discrete 3-point factorized Laplacian operator without discrete eigenvalue recalibration.

---

### 3.3 Global Norm Dissipation vs. Local Pointwise Dynamics

For constant hyperdiffusivity $\mathcal{K}_z = \nu_4$, the domain-integrated quadratic norm $E_q = \frac{1}{2} \langle q, q \rangle$ satisfies:

$$\left.\frac{d E_q}{d t}\right\vert{}_{\mathrm{filter}} = -\nu_4 \int_{\Omega} \left(\frac{\partial^2 q}{\partial z^2}\right)^2 dz \le 0$$

For spatially varying hyperdiffusivity $\mathcal{K}_z(z)$, maintaining global energy dissipation requires formulating the operator in self-adjoint negative-semidefinite form:

$$\mathcal{D}_4(q) = -\frac{\partial^2}{\partial z^2} \left[ \mathcal{K}_z(z) \frac{\partial^2 q}{\partial z^2} \right] = -A^* \mathcal{K}_z A q \quad \text{where } A = \frac{\partial^2}{\partial z^2}$$

Integrating by parts under zero-flux boundary conditions yields:

$$\left.\frac{d E_q}{d t}\right\vert{}_{\mathrm{filter}} = \langle q, \mathcal{D}_4(q) \rangle = -\int_{\Omega} \mathcal{K}_z(z) \left(\frac{\partial^2 q}{\partial z^2}\right)^2 dz \le 0 \quad \text{for all } \mathcal{K}_z(z) \ge 0$$

Global norm dissipation guarantees overall variance decay but does not enforce a pointwise maximum principle. Biharmonic filters ($p \ge 2$) can generate localized unphysical overshoots (Gibbs ringing) near sharp gradients or inversions. Preserving positive-definiteness in scalar fields ($q_v, q_c, e$) requires pairing high-order filters with non-linear flux limiters (e.g., TVD or Zalesak multidimensional limiters).

---

### 3.4 Non-Uniform Staggered Discretization and Boundary Conditions

To enforce discrete conservation and accommodate variable grid spacing $\Delta z_k = z_{k+1} - z_k$, fourth-order hyperdiffusion is factorized into two sequential cell-centered second-order diffusion steps on a staggered grid:

1. **Intermediate Cell-Centered Curvature ($C_k$):**

$$C_k = \frac{2}{\Delta z_k + \Delta z_{k-1}} \left[ \frac{q_{k+1} - q_k}{\Delta z_k} - \frac{q_k - q_{k-1}}{\Delta z_{k-1}} \right]$$

1. **Secondary Fourth-Order Update ($\mathcal{D}_{4,k}$):**

$$\mathcal{D}_{4,k} = \frac{2}{\Delta z_k + \Delta z_{k-1}} \left[ \frac{\mathcal{K}_{z,k+1/2} \frac{C_{k+1} - C_k}{\Delta z_k} - \mathcal{K}_{z,k-1/2} \frac{C_k - C_{k-1}}{\Delta z_{k-1}}}{1} \right]$$

Hyperdiffusivity is scaled locally by cell thickness: $\nu_4(z_k) = \frac{(\Delta z_k)^4}{\alpha_d \Delta t}$.

Boundary conditions are defined in terms of curvature $C$ and physical hyperdiffusive curvature flux $F_4 = -\mathcal{K}_z \frac{\partial C}{\partial z}$:

> **Upper Boundary ($z = z_{\text{top}}$):**
>
> $$C_N = 0 \quad \text{and} \quad F_{4,N} = -\mathcal{K}_z \left.\frac{\partial C}{\partial z}\right\vert{}_{z_{\text{top}}} = 0 \quad (\text{Zero-Curvature \& Zero-Hyperdiffusive-Flux})$$
>
>
> $$\downarrow$$
>
>
>
> **Interior Domain:**
>
> $$\mathcal{D}_{4,k} = -\frac{\partial^2}{\partial z^2} \left[ \mathcal{K}_z(z) \frac{\partial^2 q}{\partial z^2} \right] \quad (\text{Biharmonic Damping})$$
>
>
> $$\uparrow$$
>
>
>
> **Lower Boundary ($z = z_0$, Surface Layer):**
>
> $$C_1 = 0 \quad \text{and} \quad \mathcal{D}_{4,1} = \mathcal{D}_{4,2} = 0 \quad (\text{Surface-Layer Clamping})$$
>
>

Surface-layer clamping prevents artificial numerical filtering within cells $k=1,2$, preserving surface fluxes ($H_s, u_*$) diagnosed via MOST. Upper boundary zero-hyperdiffusive-flux conditions prevent wave reflection at domain tops.

---

### 3.5 Time Integration Strategies and Stencil-Dependent CFL Stability Bounds

Explicit time integration imposes a stability bound determined by the discrete eigenvalue spectrum of the chosen spatial stencil and temporal scheme. For a factorized 3-point second-difference biharmonic operator with Forward Euler time stepping:

$$\lambda_{\max}(-\Delta_h) = \frac{4}{h^2} \implies \lambda_{\max}(\Delta_h^2) = \frac{16}{h^4}$$

Requiring $1 - \nu_4 \Delta t \lambda_{\max}(\Delta_h^2) \ge -1$ yields the explicit stability criterion:

$$\Delta t \le \frac{h^4}{8 \nu_4} = 0.125 \frac{h^4}{\nu_4}$$

Specific stability multipliers vary by stencil and time-stepping scheme (e.g., $\approx 0.11$ for RK3). On fine vertical grids ($\Delta z_{\min} < 10\text{ m}$), explicit restrictions require specialized time integration schemes:

1. **Super-Time-Stepping (STS):** Employs Chebyshev polynomials over $M$ internal sub-steps to relax explicit stability bounds by factors of $\mathcal{O}(M^2)$.
2. **Implicit Pentadiagonal Schemes:** Solves $[\mathbf{I} + \Delta t A^* \mathcal{K}_z A] q^{n+1} = q^*$ using pentadiagonal LU decomposition (generalized Thomas algorithm), achieving unconditional numerical stability.
3. **Exponential Time Differencing (ETD):** Evaluates $q^{n+1} = \exp(-\nu_4 \mathbf{A}^2 \Delta t) q^*$ via Krylov-subspace approximations, decoupling stiff linear hyperdiffusion from non-linear physical dynamics.

---

## Pillar Synthesis & Structural Framework

The operational and physical mechanisms in boundary-layer modeling operate across four distinct mathematical spaces:

$$\begin{array}{rcccl} \text{\bf Pillar I: Physical Closure} & \quad & \zeta \longleftrightarrow Ri \longleftrightarrow (S_m, S_h) & \quad & \text{Thermodynamic \& Turbulent Transport Physics} \\ \text{\bf Pillar II: Composition Geometry} & \quad & \mathbf{R}_\zeta \xrightarrow{\ \mathbf{B}\ } \mathbf{R}_z & \quad & \text{Coordinate Transformations \& Chain-Rule Maps} \\ \text{\bf Pillar III: Discrete Filtering} & \quad & q \xrightarrow{\ \mathcal{D}_{2p}\ } q^* & \quad & \text{Scale-Selective Grid Noise Attenuation} \\ \text{\bf Dynamical System State} & \quad & \mathbf{J} \mathbf{v} = \lambda \mathbf{v} & \quad & \text{Physical Tangent Operator \& Fast-Mode Transitions} \end{array}$$

The mathematical relationships across these four spaces are summarized below:

| Dimension / Aspect | I. Physical Closure Threshold | II. Geometric Coordinate Fold | III. Numerical Damping Scale | IV. Physical Fast-Mode Transition |
| --- | --- | --- | --- | --- |
| **Defining Condition** | $Ri \to Ri_c = \frac{\beta_h}{\beta_m^2}$ | $g_1 = \zeta_z = \frac{d\zeta}{dz} = 0$ | $k_g = \frac{\pi}{h} \implies 2\Delta z$ | $\text{Re}(\lambda_f) \to 0^+$ (Tangent Jacobian) |
| **Governing Mechanism** | Stratification damping of turbulence | Balance of coordinate scaling & flux divergence ($\zeta L' = 1$) | Discrete spatial derivative stencil operator | Dynamical system stability bifurcation |
| **Mathematical Representation** | Exchange diffusivity floor ($K_m \to K_{\min}$) | Determinant collapse ($\det(\mathbf{B}) = g_1^{10} \to 0$) | Filter decay rate ($\sigma_{2p} = \nu_{2p} k^{2p}$) | Tangent operator eigenvalue zero-crossing |
| **Invariance Class** | Closure-dependent physical parameterization | Coordinate-dependent transformation artifact | Discretization-dependent numerical filter | Coordinate-invariant dynamical property |
| **Numerical Impact** | Vanishing diffusivity ($K_m, K_h \to 0$) | Inverse map singularity ($\mathbf{B}^{-1}$ undefined) | Attenuates high-wavenumber grid modes | Onset of physical turbulence collapse / transition |

This separation establishes the central thesis unifying boundary-layer modeling physics:

$$\boxed{\text{Closure Threshold} \ \neq \ \text{Coordinate Fold} \ \neq \ \text{Numerical Damping Scale} \ \neq \ \text{Physical Fast-Mode Transition}}$$
