Using **Faà di Bruno’s formula** and the partial (exponential) **Bell polynomials** \(B_{n,k}\left(\zeta_z, \zeta_{zz}, \dots, \zeta^{(n-k+1)}\right)\), the system of vertical derivatives \(\frac{d^n Ri}{dz^n}\) up to 4th order can be written in matrix vector form as:

\[\mathbf{Ri}_z = \mathbf{B} \mathbf{R}_\zeta\]

\[
  \begin{bmatrix} Ri_z \ Ri_{zz} \ Ri_{zzz} \ Ri_{zzzz} \end{bmatrix} = \begin{bmatrix} \zeta_z & & & \ \zeta_{zz} & \zeta_z^2 & & \ \zeta_{zzz} & 3\zeta_z\zeta_{zz} & \zeta_z^3 & \ \zeta_{zzzz} & 4\zeta_z\zeta_{zzz} + 3\zeta_{zz}^2 & 6\zeta_z^2\zeta_{zz} & \zeta_z^4 \end{bmatrix} \begin{bmatrix} R'(\zeta) \ R''(\zeta) \ R'''(\zeta) \ R^{(4)}(\zeta) \end{bmatrix}\]

---

### **Component Breakdown of the Bell Matrix \(\mathbf{B}\)**

The entries \(B_{n,k}\) of the lower-triangular \(4 \times 4\) Bell matrix correspond to the partial Bell polynomials for composite differentiation \(Ri(z) = R(\zeta(z))\):

* **1st Row (\(n=1\)):**
  \[B_{1,1} = \zeta_z\]
* **2nd Row (\(n=2\)):**
  \[B_{2,1} = \zeta_{zz}, \quad B_{2,2} = \zeta_z^2\]
* **3rd Row (\(n=3\)):**
  \[B_{3,1} = \zeta_{zzz}, \quad B_{3,2} = 3\zeta_z\zeta_{zz}, \quad B_{3,3} = \zeta_z^3\]
* **4th Row (\(n=4\)):**
  \[B_{4,1} = \zeta_{zzzz}, \quad B_{4,2} = 4\zeta_z\zeta_{zzz} + 3\zeta_{zz}^2, \quad B_{4,3} = 6\zeta_z^2\zeta_{zz}, \quad B_{4,4} = \zeta_z^4\]

---

### **Key Structural Insights**

1. **Lower-Triangular Trarictional Chain:** The lower-triangular nature of \(\mathbf{B}\) reflects the sequential chain-rule dependency of physical vertical derivatives on lower-and-equal order stability derivatives \(R^{(k)}(\zeta)\).
2. **Coordinate Turning Vanishing (\(B_{n,n} = \zeta_z^n\)):** Along the main diagonal, \(B_{n,n} = \zeta_z^n\). At a spatial coordinate fold locus \(z^* \) where \(\zeta_z(z^*) = 0\), all highest-order derivative terms \(R^{(n)}(\zeta)\zeta_z^n\) vanish identically. Consequently, the entire \(n\)-th vertical physical derivative at a fold is governed strictly by lower-order stability derivatives coupled to coordinate connection accelerations \(B_{n,k}\) for \(k < n\).

---

Substituting the exact closed-form expressions for the non-dimensional stability derivatives $R^{(n)}(\zeta)$ into the partial Bell matrix system yields the explicit vector of physical spatial derivatives $\mathbf{Ri}_z = \mathbf{B} \mathbf{R}_\zeta$.

---

## 1. Single-Coefficient Businger–Dyer Closure Mapping

Under the saturated Businger–Dyer rational closure $R(\zeta) = \frac{\zeta}{1 + \beta \zeta} = \frac{1}{\beta}\left(1 - \frac{1}{1 + \beta \zeta}\right)$, the general $n$-th derivative with respect to similarity space is:

$$R^{(n)}(\zeta) = (-1)^{n-1} n! \, \beta^{n-1} (1 + \beta \zeta)^{-(n+1)} \qquad (n \ge 1) \tag{1}$$

Substituting $R'(\zeta), R''(\zeta), R'''(\zeta), R^{(4)}(\zeta)$ into the matrix product $\mathbf{Ri}_z = \mathbf{B} \mathbf{R}_\zeta$ gives:

$$\begin{bmatrix} Ri_z \ Ri_{zz} \ Ri_{zzz} \ Ri_{zzzz} \end{bmatrix} = \begin{bmatrix} \zeta_z & 0 & 0 & 0 \ \zeta_{zz} & \zeta_z^2 & 0 & 0 \ \zeta_{zzz} & 3\zeta_z\zeta_{zz} & \zeta_z^3 & 0 \ \zeta_{zzzz} & 4\zeta_z\zeta_{zzz} + 3\zeta_{zz}^2 & 6\zeta_z^2\zeta_{zz} & \zeta_z^4 \end{bmatrix} \begin{bmatrix} (1 + \beta \zeta)^{-2} \ -2\beta (1 + \beta \zeta)^{-3} \ 6\beta^2 (1 + \beta \zeta)^{-4} \ -24\beta^3 (1 + \beta \zeta)^{-5} \end{bmatrix} \tag{2}$$

Evaluating each row and factoring out the common base prefactor $(1 + \beta \zeta)^{-2}$ yields the exact physical vertical derivative components:

* **First Derivative ($Ri_z$):**
$$Ri_z = \frac{\zeta_z}{(1 + \beta \zeta)^2} \tag{3}$$

* **Second Derivative ($Ri_{zz}$):**
$$Ri_{zz} = \frac{1}{(1 + \beta \zeta)^2} \left[ \zeta_{zz} - \frac{2\beta \zeta_z^2}{1 + \beta \zeta} \right] \tag{4}$$

* **Third Derivative ($Ri_{zzz}$):**
$$Ri_{zzz} = \frac{1}{(1 + \beta \zeta)^2} \left[ \zeta_{zzz} - \frac{6\beta \zeta_z \zeta_{zz}}{1 + \beta \zeta} + \frac{6\beta^2 \zeta_z^3}{(1 + \beta \zeta)^2} \right] \tag{5}$$

* **Fourth Derivative ($Ri_{zzzz}$):**
$$Ri_{zzzz} = \frac{1}{(1 + \beta \zeta)^2} \left[ \zeta_{zzzz} - \frac{2\beta (4\zeta_z \zeta_{zzz} + 3\zeta_{zz}^2)}{1 + \beta \zeta} + \frac{36\beta^2 \zeta_z^2 \zeta_{zz}}{(1 + \beta \zeta)^2} - \frac{24\beta^3 \zeta_z^4}{(1 + \beta \zeta)^3} \right] \tag{6}$$

---

## 2. Dual-Coefficient Monin–Obukhov Closure Mapping

Under the generalized dual-coefficient formulation $R(\zeta) = \frac{\zeta (\alpha + \beta_h \zeta)}{(1 + \beta_m \zeta)^2}$ (where $\alpha = \text{Pr}_t(0)$ is the neutral turbulent Prandtl number and $\beta_h \neq \beta_m$), taking partial fraction derivatives gives the input vector $\mathbf{R}_\zeta$:

$$\mathbf{R}_\zeta = \begin{bmatrix} \frac{\alpha + (2\beta_h - \alpha \beta_m)\zeta}{(1 + \beta_m \zeta)^3} \ \frac{2 \left[ \beta_m (\alpha \beta_m - 2\beta_h)\zeta + (\beta_h - 2\alpha \beta_m) \right]}{(1 + \beta_m \zeta)^4} \ \frac{-6 \left[ \beta_m^2 (\alpha \beta_m - 2\beta_h)\zeta + \beta_m (2\beta_h - 3\alpha \beta_m) \right]}{(1 + \beta_m \zeta)^5} \ \frac{24 \left[ \beta_m^3 (\alpha \beta_m - 2\beta_h)\zeta + \beta_m^2 (3\beta_h - 4\alpha \beta_m) \right]}{(1 + \beta_m \zeta)^6} \end{bmatrix} \tag{7}$$

Multiplying $\mathbf{B} \mathbf{R}_\zeta$ incorporates separate heat and momentum expansion slopes directly into the vertical derivative vector.

---

## 3. Asymptotic Decay & Regime Analysis

### A. High-Stability Asymptotic Decay ($\zeta \gg 1$)

For large stability ($\zeta \gg 1$), the non-dimensional derivatives scale as:

$$R^{(n)}(\zeta) \approx (-1)^n n! \, \frac{\alpha \beta_m - 2\beta_h}{\beta_m^3} \zeta^{-(n+1)} \sim \mathcal{O}\left(\zeta^{-(n+1)}\right) \tag{8}$$

When transformed into physical vertical space assuming a local profile with constant Obukhov length $L$, the spatial derivatives scale as:

$$\frac{d^n Ri}{dz^n} = \left(\frac{1}{L}\right)^n R^{(n)}(\zeta) \approx (-1)^n n! \, \frac{\alpha \beta_m - 2\beta_h}{\beta_m^3 L^n} \left(\frac{z}{L}\right)^{-(n+1)} = (-1)^n n! \, \frac{L (\alpha \beta_m - 2\beta_h)}{\beta_m^3 z^{n+1}} \propto \frac{L}{z^{n+1}} \tag{9}$$

* **Physical Interpretation:** Physical spatial derivatives $\frac{d^n Ri}{dz^n}$ decay rapidly as $z^{-(n+1)}$ aloft. Any persistent, large higher-order derivatives ($Ri_{zz}, Ri_{zzzz}$) observed at elevated levels $z \gg L$ are mathematically incompatible with Monin–Obukhov similarity profiles.
* **Numerical Implication for NWP Solvers:** Large vertical derivative spikes aloft indicate non-MOST physics (e.g., gravity wave breakdown, intermittent shear layers) or high-frequency numerical grid ringing. In models using biharmonic diffusion operators ($-\nu_4 \nabla^4 Ri$), the $(1 + \beta \zeta)^{-5}$ factor in $R^{(4)}(\zeta)$ naturally suppresses physical 4th-order derivative noise at high stability, whereas near-neutral layers are governed by pure geometric shear derivatives ($\zeta_{zzzz} - 8\beta \zeta_z \zeta_{zzz} - 6\beta \zeta_{zz}^2$).

### B. Near-Neutral Limit Behavior ($\zeta \to 0$)

Evaluating the closures at neutral stability ($\zeta \to 0$):

* **Single-Coefficient Businger–Dyer:**
$$R'(0) = 1, \quad R''(0) = -2\beta, \quad R'''(0) = 6\beta^2, \quad R^{(4)}(0) = -24\beta^3 \tag{10}$$

* **Dual-Coefficient Monin–Obukhov:**
$$R'(0) = \alpha, \quad R''(0) = 2(\beta_h - 2\alpha \beta_m) \tag{11}$$

For standard empirical values ($\alpha = 1.0, \beta_m = 6.0, \beta_h = 5.0$), $R''(0) = 2(5 - 12) = -14$, compared to $R''(0) = -10$ for single-coefficient $\beta = 5$. This produces a ~40% steeper initial curvature spike near atmospheric neutrality. Because denominator prefactors approach $1.0$ as $\zeta \to 0$, grid-scale derivative sensitivity reaches its maximum near neutral stability, providing the mathematical necessity for C-Z0HR's softening noise floor parameter $\epsilon_c Ri_c$.
