To clarify the complex interplay between atmospheric physics and numerical modeling, the following sections provide detailed derivations for the core diagnostics, stability functions, and grid-bias corrections discussed in the sources.

### **1. Derivation of the Gradient Richardson Number ($Ri_g$) from MOST**
The **Monin–Obukhov Similarity Theory (MOST)** relates non-dimensional vertical gradients of wind speed ($\phi_m$) and potential temperature ($\phi_h$) to the stability parameter $\zeta = z/L$.

*   **Momentum Gradient:** $\phi_m(\zeta) = \frac{\kappa z}{u_*} \frac{\partial U}{\partial z}$
*   **Heat Gradient:** $\phi_h(\zeta) = \frac{\kappa z}{\theta_*} \frac{\partial \theta}{\partial z}$

The **gradient Richardson number ($Ri_g$)** is defined as the ratio of buoyancy to shear:
$$Ri_g = \frac{(g/\theta) \partial\theta/\partial z}{(\partial U/\partial z)^2}$$

By substituting the expressions for $\partial U/\partial z$ and $\partial \theta/\partial z$ from the equations above, we derive $Ri_g$ as a function of $\zeta$:
$$Ri_g(\zeta) = \frac{g}{\theta} \cdot \frac{\phi_h \theta_* / \kappa z}{(\phi_m u_* / \kappa z)^2} = \frac{g \theta_* \kappa z}{\theta u_*^2} \cdot \frac{\phi_h}{\phi_m^2}$$

Recalling that the **Obukhov length** is $L = \frac{u_*^2 \theta}{\kappa g \theta_*}$, the term $\frac{g \theta_* \kappa z}{\theta u_*^2}$ simplifies to $z/L$, or $\zeta$. Thus:
$$\mathbf{Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}}$$.

---

### **2. Corrected Stability Function Formulation ($f_m, f_h$)**
The **Richardson number formulation** developed by **England and McNider (1995)** expresses stability functions directly in $Ri$-space. For a linear stable MOST profile where $\phi = 1 + \beta \zeta$, the sources derive the **momentum-limiting function** $f_m$:

1.  Start with the identity: $f_m(Ri) = \frac{1}{\phi_m(\zeta(Ri))^2}$.
2.  Substitute the linear form $\phi_m = 1 + \beta \zeta$: $f_m = \frac{1}{(1 + \beta \zeta)^2}$.
3.  Solve the relationship $Ri = \frac{\zeta}{1 + \beta \zeta}$ for $\zeta$: $\zeta = \frac{Ri}{1 - \beta Ri}$.
4.  Substitute this $\zeta$ back into the $f_m$ equation to reach the final **short-tail quadratic form**:
$$\boxed{\mathbf{f_m(Ri) = (1 - \beta Ri)^2 \quad \text{for } 0 < Ri < 1/\beta}}$$.

**Critical Corrections:** The sources note that the original 1995 publication required correction for a **"tipping error"** (a missing minus sign in the solution for the stability parameter $\zeta$) and a **wrongly derived relation** for the heat function ($f_h$) that otherwise produced nonphysical negative values.

---

### **3. Derivation of the Neutral Curvature Invariant ($\Delta$) and Bias Ratio ($B$)**
The **grid-curvature bias** occurs because the layer-averaged **bulk Richardson number ($Ri_b$)** underestimates the local $Ri_g$.

*   **Curvature Formula:** The second derivative of $Ri_g$ with respect to $\zeta$ is derived using logarithmic sensitivities:
$$\frac{d^2 Ri_g}{d\zeta^2} = F\left[2V_{\log} + \zeta(V_{\log}^2 - W_{\log})\right]$$.
*   **The Invariant ($\Delta$):** At the neutral limit ($\zeta \to 0$), the curvature reduces to **$2\Delta$**, where:
$$\mathbf{\Delta = a_h - 2a_m}$$ (for linear stable regimes).
*   **Bias Ratio ($B$):** For a layer defined by $[z_0, z_1]$, Jensen's Inequality dictates that for **concave-down** functions ($\Delta < 0$), the mean is less than the midpoint value. The bias ratio is defined as:
$$\mathbf{B = \frac{Ri_g(z_g)}{Ri_b} > 1}$$.
This mathematical divergence proves that **coarse grids overestimate mixing** (over-mixing).

---

### **4. Loss Function Derivation for Physics-Informed Neural Networks (PINNs)**
To rectify these biases, PINNs embed physical laws directly into the neural network's training via a **composite loss function ($\mathcal{L}$)**.

$$\mathbf{\mathcal{L} = \mathcal{L}_{data} + \lambda \mathcal{L}_{physics} + \mu \mathcal{L}_{BC}}$$.

*   **Data Loss ($\mathcal{L}_{data}$):** Calculated as the **Mean Squared Error (MSE)** between model predictions and observations: $\frac{1}{N}\sum(y_{pred} - y_{obs})^2$.
*   **Physics Loss ($\mathcal{L}_{physics}$):** Evaluates the **residual of the governing Partial Differential Equation (PDE)** at collocation points. For a system governed by $N[u] = f$, then $\mathcal{L}_{physics} = \text{MSE}(N[u_{pred}] - f)$.
*   **Regularization Terms ($\lambda, \mu$):** Weighting factors that prevent the network from violating **conservation laws** (mass, momentum, energy) or boundary conditions.

***

**Analogy**
Deriving these atmospheric diagnostics is like **calculating the fuel efficiency of a car**. $Ri_g$ is your **instantaneous fuel consumption** (stability at a single point), while $Ri_b$ is your **average fuel consumption over a long trip** (stability over a grid layer). If the road is mostly hills (curvature), simply averaging your starting and ending fuel levels will always give a different number than checking your gauge at every peak. The **Curvature Invariant ($\Delta$)** is the "hilly-ness" of the terrain that determines how much your average will deviate from the truth.