# High-Order Hyperdiffusion and Richardson Number Dynamics in Atmospheric Boundary-Layer Modeling

## 1. Theoretical Foundations: Gradient Richardson Number Dynamics and Monin–Obukhov Linkage

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

$$\zeta(Ri) = \frac{\alpha_\theta - 2 Ri \beta_m - \sqrt{\alpha_\theta^2 + 4(\beta_h - \alpha_\theta \beta_m)Ri}}{2(Ri \beta_m^2 - \beta_h)}, \qquad 0 \le Ri < Ri_c = \frac{\beta_h}{\beta_m^2}$$

A crucial mathematical constraint in this quadratic solution is the selection of the negative root before the radical. Selecting the positive root yields:

$$\lim_{Ri \to 0^+} \zeta_{+}(Ri) = \frac{2\alpha_\theta}{0} \to \frac{\alpha_\theta}{-\beta_h} \neq 0$$

Evaluating the positive root at $Ri = 0$ yields a non-zero constant, creating an unphysical step discontinuity at neutral stability and violating the foundational physical boundary condition $\zeta(0) = 0$. Selecting the negative root strictly enforces $\zeta(0) = 0$ through L'Hôpital's rule:

$$\lim_{Ri \to 0^+} \zeta_{-}(Ri) = \lim_{Ri \to 0^+} \frac{\frac{d}{dRi}\left[\alpha_\theta - 2 Ri \beta_m - \sqrt{\alpha_\theta^2 + 4(\beta_h - \alpha_\theta \beta_m)Ri}\right]}{\frac{d}{dRi}\left[2(Ri \beta_m^2 - \beta_h)\right]} = \frac{-2\beta_m - \frac{2(\beta_h - \alpha_\theta \beta_m)}{\alpha_\theta}}{-2\beta_h} = \frac{\beta_h/\alpha_\theta}{\beta_h} = \frac{1}{\alpha_\theta}$$

When equal diffusivities at neutrality are assumed ($\alpha_\theta = 1$) alongside matching empirical coefficients ($\beta_m = \beta_h = \beta$), the critical Richardson number becomes $Ri_c = 1/\beta$. In this equal-diffusivity limit ($Pr_0 = 1$), the exchange stability functions collapse into a smooth quadratic expression:

$$S_m(Ri) = S_h(Ri) = \begin{cases} \left(1 - \frac{Ri}{Ri_c}\right)^2, & 0 \le Ri < Ri_c \\ 0, & Ri \ge Ri_c \end{cases}$$

This smooth quadratic form guarantees continuous first derivatives across the extinction threshold ($Ri \to Ri_c^-$). Because $\partial S_m / \partial Ri \to 0$ as $Ri$ approaches $Ri_c$, the formulation preserves $C^1$ derivative continuity, preventing artificial limit-cycle oscillations, grid-scale numerical noise, and Jacobian singularities during implicit time integration.

While classical SBL parameterizations impose a hard critical cutoff ($S_m = 0$ for $Ri \ge Ri_c$), modern Very Stable Boundary Layer (VSBL) research (e.g., Gryanik et al., 2020; Energy and Flux Budget [EFB] theory, Casasanta et al., 2025) replaces zero-extinction limits with non-zero asymptotic tails ($\phi_m, \phi_h \sim \zeta^1$ as $\zeta \to \infty$). These non-zero tails account for wave-turbulence transport and intermittent mixing under strong thermal stratification, preventing catastrophic nocturnal runaway surface cooling and physical decoupling in land-surface models.

#### Unstable Boundary Layer ($UBL, Ri < 0$)

Under unstable conditions, adopting standard Businger–Dyer empirical exponents ($\alpha_m = -1/4, \alpha_h = -1/2$) maps the relation between $\zeta$ and $Ri$ into a Master Cubic polynomial in $x \equiv \zeta$:

$$b_m \zeta^3 - \zeta^2 - b_h \left(\frac{Ri}{\alpha_\theta}\right)^2 \zeta + \left(\frac{Ri}{\alpha_\theta}\right)^2 = 0$$

When the momentum and scalar empirical exponents match identically ($b_m = b_h = B_u \approx 16$ and $\alpha_\theta = 1.0$), the transformation identity experiences an exact algebraic collapse:

$$\phi_m(\zeta) = (1 - B_u \zeta)^{-1/4} \implies S_m(Ri) = (1 - B_u \zeta)^{1/2}$$

$$\phi_h(\zeta) = (1 - B_u \zeta)^{-1/2} \implies S_h(Ri) = (1 - B_u \zeta)^{3/4} = S_m^{3/2}(Ri)$$

Substituting $S_m(Ri)$ and $S_h(Ri)$ into the dual transformation identity yields:

$$\zeta = Ri \frac{S_h(Ri)}{S_m^{3/2}(Ri)} = Ri \frac{(1 - B_u \zeta)^{3/4}}{\left[(1 - B_u \zeta)^{1/2}\right]^{3/2}} = Ri \frac{(1 - B_u \zeta)^{3/4}}{(1 - B_u \zeta)^{3/4}} \equiv Ri$$

Under this "Dyer collapse," $\zeta$ equals $Ri$ identically across the entire unstable domain, eliminating mapping error and cubic solvers entirely.

For non-equal exponents ($b_m \neq b_h$), floating-point operations near neutral stability ($\vert{}Ri\vert{} \ll 1$) can suffer from numerical loss of precision if evaluated via trigonometric cubic roots. Numerical stability near $Ri \to 0^-$ is preserved by applying a series reversion via Lagrange inversion polynomials:

$$\zeta(Ri) \approx \frac{Ri}{\alpha_\theta} \left[ 1 + \lambda_1 \left(\frac{Ri}{\alpha_\theta}\right) + \lambda_2 \left(\frac{Ri}{\alpha_\theta}\right)^2 + \mathcal{O}(Ri^3) \right]$$

where $\lambda_1 = b_m - b_h$ and $\lambda_2 = (b_m - b_h)^2 + \frac{1}{2}b_m b_h$.

#### Numerical Regularization & GPU Alignment

To maximize computational throughput on hardware architectures sensitive to thread divergence, conditional logic (`if`/`else` branching across $Ri = 0$) must be removed. Zero-Offset Hyperbolic Regularization (Z0HR) replaces piecewise logic by decomposing $Ri$ into $C^\infty$-smooth positive (stable) and negative (unstable) projections using an $\epsilon$-regularized absolute value:

$$\vert{}Ri\vert{}_\epsilon = \sqrt{Ri^2 + \epsilon^2}, \qquad Ri_+ = \frac{1}{2}\left(Ri + \vert{}Ri\vert{}_\epsilon - \epsilon\right), \qquad Ri_- = \frac{1}{2}\left(Ri - \vert{}Ri\vert{}_\epsilon + \epsilon\right)$$

Smooth weight functions $w_{\text{stable}} = \frac{1}{2}\left(1 + \frac{Ri}{\vert{}Ri\vert{}_\epsilon}\right)$ and continuous hyperbolic tangent embeddings:

$$S_m(Ri) = \frac{1}{2}\left[1 - \tanh(\sigma Ri)\right] S_{m,\text{uns}} + \frac{1}{2}\left[1 + \tanh(\sigma Ri)\right] S_{m,\text{sta}}$$

allow unified, branch-free execution. Removing conditional branching eliminates GPU warp divergence and CPU instruction pipeline stalls, accelerating execution across mass-parallel grid cells.

---

## 2. Higher-Order Spatial Derivatives and Partial Bell Matrix Mappings

Evaluating higher-order vertical spatial derivatives of the Gradient Richardson Number ($Ri_z, Ri_{zz}, Ri_{zzz}, Ri_{zzzz}$) is essential for high-order spatial discretizations, sub-grid scale diffusion models, and turbulence stability analyses. When Monin–Obukhov stability profiles $R(\zeta)$ are transformed into physical vertical coordinates $z$, the derivative chains generate non-linear couplings between intrinsic stability curvature and geometric coordinate variations. A mathematically rigorous mapping framework is required to transform non-dimensional derivative vectors into physical spatial derivative vectors without truncation errors.

### 2.1 Faà di Bruno’s Chain Rule and the Partial Bell Matrix ($\mathbf{B}$)

Differentiating a composite function $Ri(z) = R(\zeta(z))$ $m$ times with respect to physical height $z$ is governed by Faà di Bruno’s higher-order chain rule formula:

$$\frac{d^m}{dz^m} R(\zeta(z)) = \sum_{n=1}^{m} R^{(n)}(\zeta) \cdot B_{m,n}\left(\zeta_z, \zeta_{zz}, \dots, \zeta_{(m-n+1)z}\right)$$

where $B_{m,n}$ represents the partial exponential Bell polynomial (Bell polynomial of the second kind), which combinatorially counts the partitions of a set of $m$ elements into $n$ non-empty subsets.

To maintain unambiguous notation between derivative operators and local coefficient values, let $g_j \equiv \frac{d^j \zeta}{d z^j}$ denote spatial coordinate derivatives (where $g_1 = \zeta_z, g_2 = \zeta_{zz}$, etc.). Collecting derivative components up to fourth order constructs the lower-triangular $4 \times 4$ Partial Bell Transformation Matrix system $\mathbf{Ri}_z = \mathbf{B} \mathbf{R}_\zeta$:

$$\begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix} = \begin{bmatrix} \zeta_z & 0 & 0 & 0 \\ \zeta_{zz} & \zeta_z^2 & 0 & 0 \\ \zeta_{zzz} & 3\zeta_z\zeta_{zz} & \zeta_z^3 & 0 \\ \zeta_{zzzz} & 4\zeta_z\zeta_{zzz} + 3\zeta_{zz}^2 & 6\zeta_z^2\zeta_{zz} & \zeta_z^4 \end{bmatrix} \begin{bmatrix} R'(\zeta) \\ R''(\zeta) \\ R'''(\zeta) \\ R^{(4)}(\zeta) \end{bmatrix}$$

The individual matrix entries $B_{m,n}$ can be assembled using two equivalent recurrence relations:

* **Combinatorial Summation Recurrence (Faà di Bruno Identity):**

$$B_{m+1, n} = \sum_{j=n-1}^{m} \binom{m}{j} g_{m+1-j} B_{j, n-1} \quad \text{with } B_{0,0} = 1$$

* **Differential Recurrence (Subroutine-Friendly Iteration):**

$$B_{m+1, n} = \zeta_z B_{m, n-1} + \frac{d}{dz} B_{m, n} \quad \text{with boundary conditions } B_{m,0} = 0 \text{ for } m \ge 1 \text{ and } B_{m,k} = 0 \text{ for } k > m$$

---

### 2.2 Inverse Bell Mapping ($\mathbf{B}^{-1}$) and Observational Parameter Calibration

Because $\mathbf{B}$ is lower-triangular with diagonal entries $B_{k,k} = \zeta_z^k$, its determinant is $\det(\mathbf{B}) = \zeta_z^{10}$. When $\zeta_z \neq 0$, inverting this matrix maps physical spatial derivatives back into non-dimensional similarity derivatives ($\mathbf{R}_\zeta = \mathbf{B}^{-1} \mathbf{Ri}_z$):

$$\begin{bmatrix} R'(\zeta) \\ R''(\zeta) \\ R'''(\zeta) \\ R^{(4)}(\zeta) \end{bmatrix} = \begin{bmatrix} \frac{1}{\zeta_z} & 0 & 0 & 0 \\ -\frac{\zeta_{zz}}{\zeta_z^3} & \frac{1}{\zeta_z^2} & 0 & 0 \\ \frac{3\zeta_{zz}^2 - \zeta_z \zeta_{zzz}}{\zeta_z^5} & -\frac{3\zeta_{zz}}{\zeta_z^4} & \frac{1}{\zeta_z^3} & 0 \\ \frac{10\zeta_z \zeta_{zz} \zeta_{zzz} - 15\zeta_{zz}^3 - \zeta_z^2 \zeta_{zzzz}}{\zeta_z^7} & \frac{10\zeta_{zz}^2 - 4\zeta_z \zeta_{zzz}}{\zeta_z^6} & -\frac{6\zeta_{zz}}{\zeta_z^5} & \frac{1}{\zeta_z^4} \end{bmatrix} \begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix}$$

This analytical inverse provides a diagnostic tool for observational dynamicists. Processing high-resolution field tower or Large-Eddy Simulation (LES) profile data through $\mathbf{B}^{-1}$ extracts non-dimensional similarity derivatives ($R', R'', R''', R^{(4)}$) *a priori* directly from measured physical profiles without prescribing functional forms. Near neutral stability ($\zeta \to 0$), empirical expansion slopes ($\beta_m, \beta_h$) and the neutral Prandtl number $\alpha_\theta = Pr_t(0)$ are diagnosed via:

$$\alpha_\theta = \left. \frac{Ri_z}{\zeta_z} \right\vert{}_{\zeta \to 0}, \qquad 2(\beta_h - 2\alpha_\theta \beta_m) = \left. \frac{Ri_{zz} - \frac{\zeta_{zz}}{\zeta_z} Ri_z}{\zeta_z^2} \right\vert{}_{\zeta \to 0}$$

Crucially, at or near a coordinate fold knee ($\zeta_z \to 0$), the transformation matrix becomes rank-deficient ($\det(\mathbf{B}) = 0$). Consequently, observational parameter calibration routines cannot directly invert field profile data via $\mathbf{B}^{-1}$ at a knee. Near coordinate folds, local chart transitions, singular-value decomposition (SVD), or Tikhonov regularization $\mathbf{D}_\lambda^{(2)} = (\mathbf{A}^T \mathbf{A} + \lambda \mathbf{R})^{-1} \mathbf{A}^T$ must be applied to prevent division-by-zero singularities.

---

### 2.3 Asymptotic Limits and NWP Solver Physics

#### Near-Neutral Limit ($\zeta \to 0$)

Evaluating non-dimensional profile curvature $R''(0)$ near neutral stability reveals a key mathematical distinction between single-coefficient and dual-coefficient MOST closures:

* **Single-Coefficient Closure** ($R(\zeta) = \frac{\zeta}{1 + \beta \zeta}$, $\beta = 5.0$):

$$R''(0) = -2\beta = -10.0$$

* **Dual-Coefficient Closure** ($R(\zeta) = \frac{\zeta(\alpha_\theta + \beta_h \zeta)}{(1 + \beta_m \zeta)^2}$, $\alpha_\theta = 1.0, \beta_m = 6.0, \beta_h = 5.0$):

$$R''(0) = 2(\beta_h - 2\alpha_\theta \beta_m) = 2(5.0 - 12.0) = -14.0$$

The dual-coefficient model produces a ~40% steeper initial curvature spike near neutrality ($R''(0) = -14.0$ vs. $-10.0$). On discrete C-grid vertical meshes, this steeper initial curvature accelerates spatial gradient variations near the surface:

$$\lim_{\zeta \to 0} Ri_{zzzz} = \zeta_{zzzz} + R''(0) (4\zeta_z \zeta_{zzz} + 3\zeta_{zz}^2) + R'''(0) (6\zeta_z^2 \zeta_{zz}) + R^{(4)}(0) \zeta_z^4$$

Because pre-factor damping terms approach unity near neutrality, grid-driven spatial derivatives excite numerical ringing on discrete vertical grids. This provides the mathematical justification for mandatory softening noise floor parameters ($\epsilon_c Ri_c > 0$) in surface-layer modules to bound spatial derivative amplification near neutrality without corrupting physical turbulent fluxes.

#### High-Stability Limit ($\zeta \gg 1$)

Under strongly stable conditions ($\zeta \to \infty$), higher-order non-dimensional derivatives exhibit power-law decay:

$$R^{(n)}(\zeta) \approx (-1)^n n! \, \frac{\alpha_\theta \beta_m - 2\beta_h}{\beta_m^3} \zeta^{-(n+1)} \sim \mathcal{O}\left(\zeta^{-(n+1)}\right)$$

Assuming a locally constant Obukhov length $L$, physical vertical derivatives scale as:

$$\frac{d^n Ri}{dz^n} \propto \frac{L}{z^{n+1}}$$

For fourth-order derivatives, $Ri_{zzzz} \propto z^{-5}$. This rapid decay proves that standard MOST closures naturally suppress numerical noise aloft. Consequently, if persistent grid-scale oscillations are observed aloft in strongly stable layers during model execution, they indicate non-MOST physical processes (such as gravity wave breaking or shear instability) or spurious wave reflection at domain boundaries, rather than inherent MOST model instability.

---

## 3. Mechanics of Biharmonic and Higher-Order Hyperdiffusion Operators

High-order hyperdiffusion operators are numerical dissipation mechanisms, not replacements for physically based turbulent transport closures. Their operational role in atmospheric models is to attenuate high-wavenumber grid noise near the spatial truncation limit ($2\Delta z$) while leaving larger, physically resolved boundary-layer structures unattenuated.

### 3.1 Scale-Selective Dissipation Mechanics

To prevent physical growth from being misinterpreted as damping, numerical filters must adhere to an unambiguous dissipative sign convention. For a model field $q(\mathbf{x}, t)$, the order-$2p$ spatial filtering tendency is defined as:

$$\left.\frac{\partial q}{\partial t}\right\vert{}_{\mathrm{filter}} = -\nu_{2p}(-\Delta)^p q = (-1)^{p+1}\nu_{2p}\Delta^p q$$

where $\Delta = \nabla^2$ is the physical Laplacian operator, $p \ge 1$ is an integer defining spatial derivative order $2p$, and $\nu_{2p} \ge 0$ is the constant hyperdiffusion coefficient.

The table below summarizes the operational properties across spatial orders:

| Spatial Derivative Order ($2p$) | Common Operator Name | Dissipative Tendency Expression | Fourier Decay Rate ($\sigma_{2p}(k)$) |
| --- | --- | --- | --- |
| **2** | Laplacian Diffusion | $+\nu_2 \Delta q$ | $\nu_2 k^2$ |
| **4** | Biharmonic Hyperdiffusion | $-\nu_4 \Delta^2 q$ | $\nu_4 k^4$ |
| **6** | Triharmonic Hyperdiffusion | $+\nu_6 \Delta^3 q$ | $\nu_6 k^6$ |
| **8** | Eighth-Order Hyperdiffusion | $-\nu_8 \Delta^4 q$ | $\nu_8 k^8$ |
| **$2p$** | General Even-Order | $(-1)^{p+1}\nu_{2p} \Delta^p q$ | $\nu_{2p} k^{2p}$ |

Applying the continuum operator to a spatial Fourier mode $\widehat{q} e^{i k z}$ yields the amplitude decay equation $\frac{d\widehat{q}}{dt} = -\nu_{2p} k^{2p} \widehat{q}$. Calibrating the hyperdiffusion coefficient to achieve a target e-folding damping time $\tau_g$ at reference grid-scale wavenumber $k_g = \pi / h$ gives:

$$\nu_{2p} = \frac{1}{\tau_g k_g^{2p}} \implies \sigma_{2p}(k) = \frac{1}{\tau_g} \left(\frac{k}{k_g}\right)^{2p}$$

Evaluating decay rates at half the reference wavenumber ($k / k_g = 0.5$, corresponding to $4\Delta z$ waves) demonstrates scale selectivity:

$$\frac{\sigma_2(0.5k_g)}{\sigma_2(k_g)} = \frac{1}{4}, \qquad \frac{\sigma_4(0.5k_g)}{\sigma_4(k_g)} = \frac{1}{16}, \qquad \frac{\sigma_6(0.5k_g)}{\sigma_6(k_g)} = \frac{1}{64}, \qquad \frac{\sigma_8(0.5k_g)}{\sigma_8(k_g)} = \frac{1}{256}$$

Increasing spatial derivative order sharpens scale selectivity, damping $2\Delta z$ noise while preserving longer resolved waves.

On a discrete uniform grid using a centered 3-point Laplacian, continuum wavenumber calibration $k^2$ is inaccurate. Calibration must use the discrete Laplacian eigenvalue $\mu(k) = \frac{4}{h^2} \sin^2\left(\frac{kh}{2}\right)$:

$$\sigma_{2p}(k) = \nu_{2p} [\mu(k)]^p$$

At the shortest grid wavelength $2h$ ($kh = \pi$), $\mu_g = 4/h^2$. The discrete matched coefficient is $\nu_{2p} = 1 / (\tau_g \mu_g^p) = h^{2p} / (4^p \tau_g)$. This discrete calibration differs substantially from using $k_g = \pi / h$ in the continuum formula ($\nu_{2p} = h^{2p} / (\pi^{2p} \tau_g)$). For fourth-order biharmonic hyperdiffusion ($p=2$), comparing the continuum coefficient to the discrete eigenvalue calibration yields:

$$\frac{\nu_{4,\mathrm{continuum}}}{\nu_{4,\mathrm{discrete}}} = \frac{h^4 / (\pi^4 \tau_g)}{h^4 / (16 \tau_g)} = \frac{\pi^4}{16} \approx \frac{97.41}{16} \approx 6.09$$

A naive application of the continuum calibration formula overdamps $2\Delta z$ grid modes by a factor of $\frac{\pi^4}{16} \approx 6.09$, severely eroding resolved mesoscale gradients and physically valid inversion structures.

---

### 3.2 Global Norm Dissipation vs. Local Pointwise Dynamics

For periodic domains or boundary conditions where $A = -\Delta$ is a non-negative self-adjoint operator, define the domain-integrated quadratic norm $E_q = \frac{1}{2} \langle q, q \rangle$. The domain-integrated filter tendency satisfies:

$$\left.\frac{d E_q}{d t}\right\vert{}_{\mathrm{filter}} = -\nu_{2p} \langle q, A^p q \rangle = -\nu_{2p} \int_{\Omega} \left(\frac{\partial^p q}{\partial z^p}\right)^2 dz \le 0$$

This proves that hyperdiffusion globally dissipates the quadratic norm (kinetic energy for velocity components, variance for scalar fields). However, global norm dissipation does not guarantee local pointwise monotonicity. Pointwise tendencies $\left.\frac{\partial q}{\partial t}\right\vert{}_{\mathrm{filter}}$ can be locally positive near sharp spatial gradients while domain-integrated energy decreases.

Unlike second-order diffusion ($p=1$), fourth- and higher-order filters ($p \ge 2$) lack a pointwise maximum principle. Consequently, biharmonic filters near sharp temperature inversions or Low-Level Jets generate localized unphysical overshoots (Gibbs ringing). In numerical models, these unphysical oscillations can drive positive-definite scalar fields—such as water vapor mixing ratio $q_v$, cloud water $q_c$, or Turbulent Kinetic Energy $e$—negative. Preserving physical non-negativity requires pairing high-order filters with non-linear flux limiters (such as TVD or Zalesak multidimensional limiters) evaluated specifically for mass conservation.

Physical turbulent transport ($\left.\frac{\partial q}{\partial t}\right\vert{}_{\mathrm{turb}} = \frac{\partial}{\partial z}\left(K_q \frac{\partial q}{\partial z}\right)$) must be kept strictly distinct from numerical filtering ($-\nu_4 \nabla^4 q$). Direct filtering of diagnostic fields such as $Ri$ is non-equivalent to filtering prognostic state variables ($u, v, \theta_v$) and recomputing diagnostics. In dynamical models, hyperdiffusion must be applied directly to prognostic state variables, with diagnostic fields recomputed from the filtered state.

---

## 4. Kinematics of Boundary-Layer Profiling: GSPT, Singularities, and Inflection Geometry

Generalized Similarity Profile Theory (GSPT) projects continuous state-space dynamics onto a 1D non-dimensional similarity coordinate $\zeta(z) = z / L(z)$. In real boundary layers with height-varying turbulent fluxes, coordinate stretching introduces kinematic features that must be distinguished from physical fluid dynamics.

### 4.1 Kinematic Features: The Jet "Nose" vs. The Coordinate "Knee"

Stable boundary layers under nocturnal cooling frequently generate Low-Level Jets (LLJs). Modern modeling requires distinguishing two distinct profile features:

* **The LLJ "Nose" (Kinematic Momentum Maximum):** A physical fluid-dynamical feature located at jet peak height $z_{\text{LLJ}}$, where mean vertical wind vector shear vanishes ($\frac{dU}{dz} = 0$). Because vertical shear appears squared in the denominator of the Richardson number ($Ri = \frac{g}{\theta_0} \frac{\partial \theta_v / \partial z}{(\partial U / \partial z)^2 + (\partial V / \partial z)^2}$), vanishing shear causes raw point-gradient $Ri$ values to diverge toward $+\infty$.
* **The Coordinate "Knee" (Coordinate Fold Singularity):** A mathematical projection artifact. Differentiating $\zeta(z) = z/L(z)$ with respect to physical height $z$ yields $g_1 = \zeta_z = \frac{d\zeta}{dz} = \frac{1 - \zeta L'}{L}$, where $L' = dL/dz$ represents vertical flux divergence. Where local flux divergence balances similarity scaling ($\zeta L' = 1$), the coordinate derivative vanishes ($\zeta_z = 0$), creating a 1D coordinate fold singularity (loss of local coordinate invertibility).

---

### 4.2 GSPT Three-Way Curvature Decomposition

Applying total differentiation to $Ri(\zeta(z))$ decomposes physical vertical profile curvature into three components:

$$\frac{d^2 Ri}{dz^2} = \underbrace{Ri_{\zeta\zeta} \zeta_z^2}_{\text{Intrinsic Stability Geometry}} + \underbrace{Ri_\zeta \zeta_{zz}}_{\text{Coordinate/Flux Geometry}} + \underbrace{\mathcal{E}_{\Delta z}}_{\text{Estimation Error}}$$

At a coordinate fold singularity ($\zeta_z = 0, \zeta_{zz} \neq 0$), the Partial Bell Transformation Matrix $\mathbf{B}$ undergoes structural collapse:

$$\mathbf{B}\Big\vert{}_{\zeta_z = 0} = \begin{bmatrix} 0 & 0 & 0 & 0 \\ \zeta_{zz} & 0 & 0 & 0 \\ \zeta_{zzz} & 0 & 0 & 0 \\ \zeta_{zzzz} & 3\zeta_{zz}^2 & 0 & 0 \end{bmatrix} \implies \begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix}_{\zeta_z = 0} = \begin{bmatrix} 0 \\ R'(\zeta) \zeta_{zz} \\ R'(\zeta) \zeta_{zzz} \\ R'(\zeta) \zeta_{zzzz} + 3 R''(\zeta) \zeta_{zz}^2 \end{bmatrix}$$

At the fold, physical gradients flatten ($Ri_z = 0$) and intrinsic thermodynamic stability curvature vanishes ($Ri_{\zeta\zeta} \zeta_z^2 = 0$). Physical profile bending reduces purely to coordinate geometry ($Ri_{zz} = R'(\zeta) \zeta_{zz}$). Tower evaluations (e.g., CASES-99, GABLS3) confirm that over 99% of observed profile bending at a coordinate knee is driven by coordinate compression ($\zeta_{zz}$) rather than intrinsic thermodynamic stability breakdown.

At fourth order, biharmonic hyperdiffusion reactivates at the fold via the non-zero Bell polynomial entry $B_{4,2} = 3\zeta_{zz}^2$. The resulting fourth derivative contains the active term $3 R''(\zeta) \zeta_{zz}^2$, allowing biharmonic filters ($-\nu_4 \nabla^4 Ri$) to smooth grid noise at the fold without misinterpreting coordinate stretching as physical turbulence generation.

Defining $g_j \equiv \frac{d^j \zeta}{dz^j}$, the table below lists Bell polynomial entries at a coordinate fold ($\zeta_z = 0$) up to tenth order, showing surviving non-dimensional derivatives $R^{(p)}$:

| Physical Derivative Order ($2p$) | Highest Surviving Similarity Derivative | Bell Polynomial Entry $B_{2p, p}\big\vert{}_{\zeta_z=0}$ | Exact Combinatorial Entry |
| --- | --- | --- | --- |
| **2** | $R'(\zeta)$ | $g_2$ | $1 \cdot g_2$ |
| **4** | $R''(\zeta)$ | $3 g_2^2$ | $3 \cdot g_2^2 = (3!!) g_2^2$ |
| **6** | $R'''(\zeta)$ | $15 g_2^3$ | $15 \cdot g_2^3 = (5!!) g_2^3$ |
| **8** | $R^{(4)}(\zeta)$ | $105 g_2^4$ | $105 \cdot g_2^4 = (7!!) g_2^4$ |
| **10** | $R^{(5)}(\zeta)$ | $945 g_2^5$ | $945 \cdot g_2^5 = (9!!) g_2^5$ |

The surviving coefficient for derivative order $2p$ at a fold follows the double-factorial formula $(2p-1)!! g_2^p$. Computational modelers must interpret these combinatorial entries as descriptions of mathematical function composition geometry across coordinate transforms. They do not represent a physical fivefold increase in numerical filter strength or localized physical turbulence generation.

---

### 4.3 The Triple-Point Invariant-Convergence Hypothesis

To prevent modelers from mistaking coordinate fold artifacts or grid truncation errors for true physical turbulence collapse, GSPT establishes the Triple-Point Invariant-Convergence Hypothesis. Physical turbulence extinction is verified by tracking the spatial convergence of three independent physical heights as vertical grid spacing vanishes ($\Delta z \to 0$):

1. **Diffusivity Extinction Height ($z_K$):** The level where exchange diffusivity reaches its background floor ($K_m \to K_{\min}$).
2. **TKE Floor Height ($z_e$):** The level where Turbulent Kinetic Energy reaches its minimum floor ($e \to e_{\min}$).
3. **Shear Stress / TKE Gradient Extremum Height ($z_{e_z}$):** The level where vertical stress divergence or TKE gradients reach a localized peak ($\frac{\partial^2}{\partial z^2}(u_*^2) \to \text{max}$).

If the spatial spread $\Delta z_{\text{TP}} = \max\vert{}z_i - z_j\vert{}$ collapses to zero as $\Delta z \to 0$, the feature represents a physical turbulence transition rather than a grid artifact.

---

## 5. Discretization, Boundary Layer Coupling, and Advanced Time-Stepping

Implementing fourth-order hyperdiffusion on operational NWP vertical grids presents numerical challenges. Vertical grids are non-uniform and stretched, with fine spacing near the surface ($\Delta z \sim 10\text{ m}$) that coarsens aloft ($\Delta z \sim 100\text{--}500\text{ m}$). Applying uniform finite-difference stencils to stretched grids introduces severe truncation errors and breaks discrete conservation laws.

### 5.1 Non-Uniform Vertical Grid Discretization

To enforce discrete conservation and accommodate variable grid spacing, the fourth-order hyperdiffusion operator is factorized into two sequential, cell-centered second-order diffusion steps on a staggered grid:

$$\mathcal{D}_4(Ri) = -\frac{\partial^2}{\partial z^2} \left[ \mathcal{K}_z(z) \frac{\partial^2 Ri}{\partial z^2} \right]$$

For a vertical grid with node heights $z_k$ and dual-grid intervals $\Delta z_k = z_{k+1} - z_k$:

1. **Intermediate Cell-Centered Curvature ($C_k$):** Compute intermediate physical profile curvature at cell centers (full levels) using non-uniform central differences:

$$C_k = \frac{2}{\Delta z_k + \Delta z_{k-1}} \left[ \frac{Ri_{k+1} - Ri_k}{\Delta z_k} - \frac{Ri_k - Ri_{k-1}}{\Delta z_{k-1}} \right]$$

1. **Secondary Fourth-Order Update ($\mathcal{D}_{4,k}$):** Apply the secondary curvature operator to $C_k$:

$$\mathcal{D}_{4,k} = \frac{2}{\Delta z_k + \Delta z_{k-1}} \left[ \frac{C_{k+1} - C_k}{\Delta z_k} - \frac{C_k - C_{k-1}}{\Delta z_{k-1}} \right]$$

On staggered C-grids, intermediate curvature $C_k$ is evaluated at cell centers, whereas hyperdiffusivities $\mathcal{K}_z(z)$ reside at cell interfaces. To ensure uniform attenuation of $2\Delta z$ grid noise across all vertical levels without over-damping coarse cells aloft, local hyperdiffusivity scaling uses local cell thickness $\Delta z_k$:

$$\nu_4(z_k) = \frac{(\Delta z_k)^4}{\tau_{\text{diff}}}$$

where $\tau_{\text{diff}} = \alpha_d \Delta t$ (with $\alpha_d \in [2.0, 5.0]$) represents the target grid-scale e-folding damping timescale.

---

### 5.2 Boundary Condition Engineering

Fourth-order differential operators require two boundary conditions at each domain interface:

> **Upper Boundary ($z = z_{\text{top}}$):**
>
> $$C_N = 0 \quad \text{and} \quad \frac{C_N - C_{N-1}}{\Delta z_{N-1}} = 0 \quad (\text{Zero-Curvature \& Zero-Flux})$$
>
>
> $$\downarrow$$
>
>
>
> **Interior Domain:**
>
> $$\mathcal{D}_{4,k} = -\frac{\partial^2}{\partial z^2} \left[ \mathcal{K}_z \frac{\partial^2 Ri}{\partial z^2} \right] \quad (\text{Biharmonic Damping})$$
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

#### Lower Boundary ($z = z_0$, Surface Layer)

Physics near the surface is governed by Monin–Obukhov Similarity Theory. Applying artificial numerical filters inside the lowest grid cells ($k=1, 2$) corrupts the diagnosed surface sensible heat flux $H_s$ and friction velocity $u_*$. Surface-layer clamping enforces $C_1 = 0$ and $\mathcal{D}_{4,1} = \mathcal{D}_{4,2} = 0$ for $z \le z_1$.

#### Upper Boundary ($z = z_{\text{top}}$, Free Atmosphere)

At the top model boundary, filters must prevent spurious reflection of high-wavenumber computational waves back down into the boundary layer. Enforce zero-curvature and zero-flux boundary conditions:

$$C_N = 0 \quad \text{and} \quad \frac{C_N - C_{N-1}}{\Delta z_{N-1}} = 0$$

---

### 5.3 Time Integration Strategies and Stability Bounds

Explicit time integration of fourth-order hyperdiffusion imposes a strict Courant–Friedrichs–Lewy (CFL) stability limit:

$$\Delta t \le \sigma_{\max} \frac{(\Delta z_{\min})^4}{\nu_{4,\max}}$$

where $\sigma_{\max} \approx 0.0625$ for Forward Euler and $\sigma_{\max} \approx 0.11$ for third-order Runge–Kutta (RK3). On fine vertical grids ($\Delta z_{\min} < 10\text{ m}$), this fourth-power scaling forces impractically small time steps.

To bypass explicit restrictions, three numerical integration solutions are utilized:

1. **Explicit Sub-Stepping & Super-Time-Stepping (STS):**
Divides the hyperdiffusion update into $M$ internal sub-steps $\Delta t_{\text{sub}} = \Delta t / M$ per dynamical time step, or employs Chebyshev-accelerated Super-Time-Stepping polynomials to relax the explicit stability limit by a factor of $\mathcal{O}(M^2)$. This maintains numerical stability across tight near-surface grid cells ($\Delta z_{\min} < 10\text{ m}$) without forcing the global dynamics solver to reduce its step size.
2. **Implicit & Semi-Implicit Pentadiagonal Schemes:**
Formulates the biharmonic filter implicitly at time level $n+1$:

$$\left[ \mathbf{I} + \Delta t \frac{\partial^2}{\partial z^2} \left( \mathcal{K}_z \frac{\partial^2}{\partial z^2} \right) \right] q^{n+1} = q^*$$

Because the factored spatial stencil spans five vertical points ($k-2$ to $k+2$), the resulting linear system forms a pentadiagonal matrix. Solving this system via specialized pentadiagonal LU decomposition (a generalized Thomas algorithm variant) provides unconditional numerical stability, removing the fourth-power CFL limit at minimal computational cost per vertical column.
3. **Operator Splitting & Exponential Time Differencing (ETD):**
Decouples stiff numerical hyperdiffusion from non-linear physical advection and turbulent mixing using fractional-step operator splitting. The linear biharmonic operator is integrated separately via precomputed matrix exponentials or Krylov-subspace approximations ($q^{n+1} = \exp(-\nu_4 \mathbf{A}^2 \Delta t) q^*$). This permits large dynamical time steps without incurring high-frequency numerical instability or unwanted damping of resolved physical modes.
