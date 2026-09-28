# Closed-Form Stability Derivative Mappings for Monin–Obukhov Similarity Theory via Partial Bell Matrices

**Document Type:** Technical Reference & Defense Note

**Target Audience:** Atmospheric Scientists, Boundary-Layer Dynamicists, and NWP Parameterization Developers

**Core Topic:** Mathematical Derivations, Matrix System Construction, and Asymptotic Regime Analysis for MOST Stability Closure Mappings

---

## 1. Executive Summary & Physical Motivation

In Numerical Weather Prediction (NWP) models and boundary-layer turbulence parameterizations, non-dimensional Monin–Obukhov Similarity Theory (MOST) stability functions $R(\zeta)$ must be mapped into physical vertical coordinates to evaluate spatial gradients of the Gradient Richardson Number ($Ri$) and flux profile functions:

$$\mathbf{Ri}_z = \mathbf{B} \, \mathbf{R}_\zeta \implies \begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix} = \begin{bmatrix} B_{1,1} & 0 & 0 & 0 \\ B_{2,1} & B_{2,2} & 0 & 0 \\ B_{3,1} & B_{3,2} & B_{3,3} & 0 \\ B_{4,1} & B_{4,2} & B_{4,3} & B_{4,4} \end{bmatrix} \begin{bmatrix} R'(\zeta) \\ R''(\zeta) \\ R'''(\zeta) \\ R^{(4)}(\zeta) \end{bmatrix}$$

Evaluating higher-order vertical derivatives ($Ri_{zz}, Ri_{zzz}, Ri_{zzzz}$) is required for:

1. **Higher-order spatial discretization** and grid-adaptive vertical diffusion.
2. **Biharmonic hyper-diffusion operators** ($-\nu_4 \nabla^4 Ri$) used to damp sub-grid numerical noise.
3. **Stability analysis of grid-scale turbulence transitions** near neutral boundaries and under strongly stable stratification.

This note establishes the exact mathematical identity connecting Faà di Bruno's higher-order chain rule to the lower-triangular **Partial Bell Matrix** $\mathbf{B}$, provides closed-form derivations for both single-coefficient and dual-coefficient closures, and proves the asymptotic properties governing near-neutral noise floors and high-stability noise suppression.

---

## 2. Mathematical Foundations: Faà di Bruno's Formula & Partial Bell Polynomials

### 2.1 Faà di Bruno's Higher-Order Chain Rule

When differentiating a composite function $Ri(z) = R(\zeta(z))$ $m$ times with respect to the physical coordinate $z$, standard chain-rule operations scale combinatorially according to **Faà di Bruno's formula**:

$$\frac{d^m}{dz^m} R(\zeta(z)) = \sum_{n=1}^m R^{(n)}(\zeta) \cdot B_{m,n}\left(\zeta_z, \zeta_{zz}, \dots, \zeta_{(m-n+1)z}\right)$$

where $B_{m,n}$ denotes the **partial exponential Bell polynomial** (or Bell polynomial of the second kind).

### 2.2 Algebraic Definition of Partial Bell Polynomials

Algebraically, $B_{m,n}(x_1, x_2, \dots, x_{m-n+1})$ is defined as:

$$B_{m,n}(x_1, x_2, \dots, x_{m-n+1}) = m! \sum \prod_{i=1}^{m-n+1} \frac{x_i^{j_i}}{(i!)^{j_i} j_i!}$$

where the summation is taken over all non-negative integer sequences $(j_1, j_2, \dots, j_{m-n+1})$ satisfying the dual combinatorial constraints:

$$\sum_{i=1}^{m-n+1} j_i = n \quad \text{and} \quad \sum_{i=1}^{m-n+1} i \cdot j_i = m$$

Combinatorially, $B_{m,n}$ counts the partitions of a set of $m$ elements into $n$ non-empty subsets, where $j_i$ specifies the number of blocks of size $i$.

---

## 3. Construction of the Partial Bell Matrix System ($\mathbf{B}$)

### 3.1 Recurrence Relations for Matrix Entries $B_{m,n}$

To construct the transformation matrix entries efficiently in code or symbolic algebra, two recurrence forms are available:

#### A. Combinatorial Summation Recurrence (Faà di Bruno Identity)

$$B_{m+1, n} = \sum_{j=n-1}^{m} \binom{m}{j} \left( \frac{d^{m+1-j}\zeta}{dz^{m+1-j}} \right) B_{j, n-1} \quad \text{with } B_{0,0} = 1$$

#### B. Differential Recurrence (Subroutine-Friendly Form)

Differentiating row $m$ directly with respect to $z$ yields an iterative differential rule:

$$B_{m+1, n} = \zeta_z B_{m, n-1} + \frac{d}{dz} B_{m, n}$$

subject to boundary conditions $B_{0,0} = 1$, $B_{m,0} = 0$ for $m \ge 1$, and $B_{m,n} = 0$ for $n > m$.

---

### 3.2 Step-by-Step Recursive Derivation of Matrix Entries & Derivatives

#### A. Base Row ($m = 1$): First Spatial Derivative $Ri_z$

For $m=1$ and $n=1$:

$$B_{1,1} = \zeta_z B_{0,0} + \frac{d}{dz}(B_{0,1}) = \zeta_z (1) + 0 = \mathbf{\zeta_z}$$

This yields $Ri_z = B_{1,1} R'(\zeta) = \mathbf{\zeta_z R'(\zeta)}$.

#### B. Recurrence Step to Row 2 ($m = 1 \to m+1 = 2$): Second Spatial Derivative $Ri_{zz}$

For the second spatial derivative, row $m=2$ evaluates as:

$$Ri_{zz} = B_{2,1} R'(\zeta) + B_{2,2} R''(\zeta)$$

1. **Calculate $B_{2,1}$ (coefficient for $R'$):** Set $m=1, n=1$ in the differential recurrence rule:

$$B_{2,1} = \zeta_z B_{1,0} + \frac{d}{dz}(B_{1,1})$$

Applying boundary condition $B_{1,0} = 0$ and base result $B_{1,1} = \zeta_z$:

$$B_{2,1} = \zeta_z (0) + \frac{d}{dz}(\zeta_z) = \mathbf{\zeta_{zz}}$$

1. **Calculate $B_{2,2}$ (coefficient for $R''$):** Set $m=1, n=2$ in the recurrence rule:

$$B_{2,2} = \zeta_z B_{1,1} + \frac{d}{dz}(B_{1,2})$$

Since $n > m$, $B_{1,2} = 0$:

$$B_{2,2} = \zeta_z (\zeta_z) + \frac{d}{dz}(0) = \mathbf{\zeta_z^2}$$

1. **Assemble $Ri_{zz}$:** Substituting $B_{2,1} = \zeta_{zz}$ and $B_{2,2} = \zeta_z^2$ into the matrix row equation yields:

$$Ri_{zz} = B_{2,1} R'(\zeta) + B_{2,2} R''(\zeta)$$

$$\mathbf{Ri_{zz} = R''(\zeta) \zeta_z^2 + R'(\zeta) \zeta_{zz}}$$

#### C. Higher-Order Rows ($m = 3$ and $m = 4$)

Applying the differential recurrence rule across higher orders generates:

* **Row 3 ($m = 3$):**
* $B_{3,1} = \zeta_z B_{2,0} + \frac{d}{dz}(B_{2,1}) = 0 + \frac{d}{dz}(\zeta_{zz}) = \mathbf{\zeta_{zzz}}$
* $B_{3,2} = \zeta_z B_{2,1} + \frac{d}{dz}(B_{2,2}) = \zeta_z \zeta_{zz} + \frac{d}{dz}(\zeta_z^2) = \mathbf{3\zeta_z \zeta_{zz}}$
* $B_{3,3} = \zeta_z B_{2,2} + \frac{d}{dz}(B_{2,3}) = \zeta_z(\zeta_z^2) + 0 = \mathbf{\zeta_z^3}$
* **Resulting Spatial Derivative:** $Ri_{zzz} = R'''(\zeta) \zeta_z^3 + 3 R''(\zeta) \zeta_z \zeta_{zz} + R'(\zeta) \zeta_{zzz}$

* **Row 4 ($m = 4$):**
* $B_{4,1} = \zeta_z B_{3,0} + \frac{d}{dz}(B_{3,1}) = \mathbf{\zeta_{zzzz}}$
* $B_{4,2} = \zeta_z B_{3,1} + \frac{d}{dz}(B_{3,2}) = \mathbf{4\zeta_z \zeta_{zzz} + 3\zeta_{zz}^2}$
* $B_{4,3} = \zeta_z B_{3,2} + \frac{d}{dz}(B_{3,3}) = \mathbf{6\zeta_z^2 \zeta_{zz}}$
* $B_{4,4} = \zeta_z B_{3,3} + \frac{d}{dz}(B_{3,4}) = \mathbf{\zeta_z^4}$
* **Resulting Spatial Derivative:** $Ri_{zzzz} = R^{(4)}(\zeta) \zeta_z^4 + 6 R'''(\zeta) \zeta_z^2 \zeta_{zz} + R''(\zeta) (4\zeta_z \zeta_{zzz} + 3\zeta_{zz}^2) + R'(\zeta) \zeta_{zzzz}$

---

### 3.3 Assembled Transformation Matrix $\mathbf{B}$

$$\mathbf{B} = \begin{bmatrix}  \zeta_z & 0 & 0 & 0 \\  \zeta_{zz} & \zeta_z^2 & 0 & 0 \\  \zeta_{zzz} & 3\zeta_z\zeta_{zz} & \zeta_z^3 & 0 \\  \zeta_{zzzz} & 4\zeta_z\zeta_{zzz} + 3\zeta_{zz}^2 & 6\zeta_z^2\zeta_{zz} & \zeta_z^4  \end{bmatrix}$$

---

## 4. Closure Vector Mappings ($\mathbf{Ri}_z = \mathbf{B}\mathbf{R}_\zeta$)

### 4.1 Single-Coefficient Businger–Dyer Closure

For the classical single-coefficient formulation:

$$R(\zeta) = \frac{\zeta}{1 + \beta \zeta}$$

The non-dimensional derivative vector $\mathbf{R}_\zeta$ is given by:

$$\mathbf{R}_\zeta = \begin{bmatrix} R'(\zeta) \\ R''(\zeta) \\ R'''(\zeta) \\ R^{(4)}(\zeta) \end{bmatrix} = \begin{bmatrix} (1 + \beta \zeta)^{-2} \\ -2\beta (1 + \beta \zeta)^{-3} \\ 6\beta^2 (1 + \beta \zeta)^{-4} \\ -24\beta^3 (1 + \beta \zeta)^{-5} \end{bmatrix}$$

Multiplying $\mathbf{B} \mathbf{R}_\zeta$ yields the closed-form physical spatial derivatives:

$$Ri_z = \frac{\zeta_z}{(1 + \beta \zeta)^2}$$

$$Ri_{zz} = \frac{1}{(1 + \beta \zeta)^2} \left[ \zeta_{zz} - \frac{2\beta \zeta_z^2}{1 + \beta \zeta} \right]$$

$$Ri_{zzz} = \frac{1}{(1 + \beta \zeta)^2} \left[ \zeta_{zzz} - \frac{6\beta \zeta_z \zeta_{zz}}{1 + \beta \zeta} + \frac{6\beta^2 \zeta_z^3}{(1 + \beta \zeta)^2} \right]$$

$$Ri_{zzzz} = \frac{1}{(1 + \beta \zeta)^2} \left[ \zeta_{zzzz} - \frac{2\beta (4\zeta_z \zeta_{zzz} + 3\zeta_{zz}^2)}{1 + \beta \zeta} + \frac{36\beta^2 \zeta_z^2 \zeta_{zz}}{(1 + \beta \zeta)^2} - \frac{24\beta^3 \zeta_z^4}{(1 + \beta \zeta)^3} \right]$$

---

### 4.2 Dual-Coefficient Monin–Obukhov Closure

For the generalized dual-coefficient formulation accounting for distinct momentum ($\beta_m$) and scalar ($\beta_h$) expansion slopes with neutral turbulent Prandtl number $\alpha = \text{Pr}_t(0)$:

$$R(\zeta) = \frac{\zeta(\alpha + \beta_h \zeta)}{(1 + \beta_m \zeta)^2}$$

The derivative components evaluated at $\zeta = 0$ yield:

$$R'(0) = \alpha$$

$$R''(0) = 2(\beta_h - 2\alpha \beta_m)$$

$$R'''(0) = 6\beta_m (3\alpha \beta_m - 2\beta_h)$$

$$R^{(4)}(0) = -24\beta_m^2 (4\alpha \beta_m - 3\beta_h)$$

---

## 5. Inverse Partial Bell Matrix & Parameter Estimation ($\mathbf{R}_\zeta = \mathbf{B}^{-1} \mathbf{Ri}_z$)

Because $\mathbf{B}$ is a lower-triangular matrix with non-zero diagonal entries $B_{k,k} = \zeta_z^k$ (provided $\zeta_z \neq 0$), $\mathbf{B}$ is strictly invertible with determinant $\det(\mathbf{B}) = \zeta_z^{10}$. Inverting this system maps physical spatial gradients back into non-dimensional similarity space:

$$\begin{bmatrix} R'(\zeta) \\ R''(\zeta) \\ R'''(\zeta) \\ R^{(4)}(\zeta) \end{bmatrix} = \begin{bmatrix} \frac{1}{\zeta_z} & 0 & 0 & 0 \\ -\frac{\zeta_{zz}}{\zeta_z^3} & \frac{1}{\zeta_z^2} & 0 & 0 \\ \frac{3\zeta_{zz}^2 - \zeta_z \zeta_{zzz}}{\zeta_z^5} & -\frac{3\zeta_{zz}}{\zeta_z^4} & \frac{1}{\zeta_z^3} & 0 \\ \frac{10\zeta_z \zeta_{zz} \zeta_{zzz} - 15\zeta_{zz}^3 - \zeta_z^2 \zeta_{zzzz}}{\zeta_z^7} & \frac{10\zeta_{zz}^2 - 4\zeta_z \zeta_{zzz}}{\zeta_z^6} & -\frac{6\zeta_{zz}}{\zeta_z^5} & \frac{1}{\zeta_z^4} \end{bmatrix} \begin{bmatrix} Ri_z \\ Ri_{zz} \\ Ri_{zzz} \\ Ri_{zzzz} \end{bmatrix}$$

### Observational Calibration Application

Given $R(\zeta) = \frac{\zeta \phi_h(\zeta)}{\phi_m^2(\zeta)}$, observational tower or Large-Eddy Simulation (LES) profile data can be processed via $\mathbf{B}^{-1}$ to diagnose local values of $\phi_m$ and $\phi_h$ without assuming specific functional profiles *a priori*. Near neutrality ($\zeta \to 0$):

$$\alpha = \left. \frac{Ri_z}{\zeta_z} \right\vert{}_{\zeta \to 0}, \qquad 2(\beta_h - 2\alpha \beta_m) = \left. \frac{Ri_{zz} - \frac{\zeta_{zz}}{\zeta_z} Ri_z}{\zeta_z^2} \right\vert{}_{\zeta \to 0}$$

---

## 6. Asymptotic Regime Analysis & NWP Solver Dynamics

### 6.1 Near-Neutral Limit ($\zeta \to 0$): Curvature Spike & Noise Floor

#### A. Curvature Spike Comparison at Neutrality

Evaluating $R''(0)$ for standard empirical parameters ($\alpha = 1.0, \beta_m = 6.0, \beta_h = 5.0$):

* **Dual-coefficient model:** $R''(0) = 2(5.0 - 2(1.0)(6.0)) = 2(5.0 - 12.0) = \mathbf{-14.0}$
* **Single-coefficient model** (with standard $\beta = 5.0$): $R''(0) = -2(5.0) = \mathbf{-10.0}$

> **Key Result:** Dual-coefficient closures exhibit a **~40% steeper initial curvature spike** near neutral stability ($R''(0) = -14$ vs. $-10$).

#### B. Justification for Softening Noise Floor Parameters ($\epsilon_c Ri_c$)

As $\zeta \to 0$, prefactor damping terms like $(1 + \beta \zeta)^{-n}$ approach $1.0$. Consequently, vertical spatial derivatives become unattenuated and are driven purely by geometric grid gradients:

$$\lim_{\zeta \to 0} Ri_{zzzz} = \zeta_{zzzz} + R''(0) (4\zeta_z \zeta_{zzz} + 3\zeta_{zz}^2) + R'''(0) (6\zeta_z^2 \zeta_{zz}) + R^{(4)}(0) \zeta_z^4$$

In discrete NWP solvers using stretched vertical grids, high spatial derivative values near the surface induce non-physical grid-scale oscillations ("ringing"). This provides the explicit mathematical proof for why boundary-layer schemes require a **noise floor parameter** $\epsilon_c Ri_c > 0$ to bound derivative growth near neutrality without distorting physical turbulent fluxes.

---

### 6.2 High-Stability Limit ($\zeta \gg 1$): Asymptotic Decay & Noise Damping

#### A. Asymptotic Derivative Scaling

For large stability parameter values ($\zeta \gg 1$), higher-order non-dimensional derivatives decay as:

$$R^{(n)}(\zeta) \approx (-1)^n n! \, \frac{\alpha \beta_m - 2\beta_h}{\beta_m^3} \zeta^{-(n+1)} \sim \mathcal{O}\left(\zeta^{-(n+1)}\right)$$

Assuming a constant local Obukhov length $L$, the physical vertical coordinate derivative scales as $\zeta_z = 1/L$, yielding:

$$\frac{d^n Ri}{dz^n} \propto \frac{L}{z^{n+1}}$$

#### B. Natural Suppression of Biharmonic Numerical Noise

In NWP models, fourth-order horizontal and vertical hyper-diffusion terms ($-\nu_4 \frac{\partial^4 Ri}{\partial z^4}$) are applied to suppress sub-grid computational noise.

Because $Ri_{zzzz} \propto z^{-5}$ aloft under strongly stable conditions, **MOST closures naturally damp higher-order numerical noise aloft**. Any persistent grid-scale oscillations observed aloft in stable layers are therefore diagnostic of:

1. Non-MOST physical processes (e.g., gravity wave breaking, shear instability, or intermittent turbulence).
2. Spurious numerical reflection at layer boundaries rather than inherent MOST model instability.

---

## 7. NWP Code Implementation Checklist

When embedding these transformations into atmospheric model subroutines (e.g., WRF, MPAS, ICON, or IFS boundary-layer packages):

| Step | Target Subroutine Component | Action Items & Algorithmic Rules |
| --- | --- | --- |
| **1** | **Grid Metric Pre-computation** | Pre-compute $4 \times 4$ matrix entries $\mathbf{B}(z)$ per column during metric initialization. Solve $\mathbf{R}_\zeta = \mathbf{B}^{-1} \mathbf{Ri}_z$ using fast forward substitution ($\mathcal{O}(m^2)$). |
| **2** | **Near-Surface Softening** | In surface cells ($z < 10\text{ m}$), apply the bound $\max(\zeta, \epsilon_c Ri_c)$ prior to evaluating $R''(\zeta)$ and $R'''(\zeta)$ to prevent curvature spike amplification. |
| **3** | **Staggered C-Grid Co-location** | Interpolate interface terms $\zeta_{z, i+1/2}$ to cell centers via thickness-weighted operators $(P \zeta_z)_i = \lambda_i^+ \zeta_{z, i+1/2} + \lambda_i^- \zeta_{z, i-1/2}$ before assembling non-linear products in Row 4. |
| **4** | **Biharmonic Time-Splitting** | Use biharmonic-modified semi-implicit time-stepping $\left(\mathbf{I} + \nu_{z4} \Delta t \boldsymbol{\delta}_{zzzz}\right) Ri^{n+1} = Ri^n + \dots$ to bypass explicit CFL restrictions ($\Delta t \propto \Delta z^4$). |
