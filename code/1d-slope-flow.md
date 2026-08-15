For a 1-D surface-layer model, you should adopt a **Non-Uniform (Stretched) Grid**, specifically one that is **logarithmically or geometrically-spaced** near the bottom. This ensures you capture the sharp gradients and the high curvature you are interested in. 📉

---

## 📐 Grid Suggestions for Surface-Layer Fidelity

The key challenge is balancing resolution near the surface (where $z \to 0$ and turbulence is dominated by friction) with the need to reach a reasonable height (e.g., the top of the surface layer).

### 1. The Log-Linear Grid (Best Practice)

The most common and effective solution for the neutral and slightly stable/unstable surface layer is a grid that closely follows the expected $\ln(z)$ profile.

* **Definition:** The vertical spacing $\Delta z$ increases with height $z$.
* **Formula:** You define the cell center heights $z_k$ such that the ratio of adjacent cell heights is constant:
    $$
    \frac{z_{k+1}}{z_k} = C \quad (\text{where } C > 1)
    $$
* **Implementation:** You define your lowest cell center $z_1$ and your top height $H$.
    1. Calculate the constant ratio $C$: $C = (H/z_1)^{1/(N-1)}$, where $N$ is the number of grid points.
    2. Generate the grid: $z_k = z_1 \cdot C^{k-1}$.
* **Advantages:**
  * High resolution exactly where the $Ri_g$ curvature is largest ($\zeta \to 0$).
  * Efficient: Uses fewer total points than a uniform grid while maintaining accuracy near the surface.
  * Maintains a constant $\Delta z/z$ ratio, which is physically consistent with the $\ln(z)$ profile.

---

### 2. Geometric Mean Height for Fluxes

While not a full grid generation method, ensure that when you calculate the fluxes (which depend on $\phi_m, \phi_h$ and $Ri_g$) between two model levels, $z_k$ and $z_{k+1}$, you use the **geometric mean height** $\mathbf{z_g = \sqrt{z_k z_{k+1}}}$ to evaluate the stability functions.

* As discussed earlier, $z_g$ is dynamically more representative of the layer than the arithmetic mean because $\phi$ is logarithmic/power-law in $z$. This helps preserve the curvature fidelity, especially on coarser grids.

---

### 3. Stretching Function (Alternative)

If you need more control over the stretching, you can use a smooth function $f(\eta)$ to map a uniform $\eta \in [0, 1]$ coordinate to your physical $z \in [z_{min}, H]$ coordinate.

* **Example (Hyperbolic Tangent $\tanh$):**
    $$
    z_k = H \cdot \frac{\tanh(\alpha \cdot \eta_k)}{\tanh(\alpha)}
    $$
    where $\alpha$ is a stretching parameter ($\alpha > 1$ gives strong resolution near $z=0$). This gives a very smooth, continuous distribution that bunches points heavily near $z=0$.

**Recommendation:** Start with the **Log-Linear Grid (Method 1)**. It is simplest to implement, directly aligned with MOST theory, and proven to give accurate results for surface-layer diagnostics.

## 🛠️ Model Structure for Few Layers

Using only a few layers (e.g., $N=5$ to $N=10$) is possible because the MOST profiles are continuous analytic functions. However, it requires careful treatment of the surface layer:

### 1. Vertical Grid Setup (Log-Stretching)

You must use a **geometrically stretched grid** to ensure the few available points capture the critical near-surface gradients.

* **Lowest Level ($z_1$)**: Set the center of your first layer, $z_1$, very close to the surface, typically $z_1 \approx 0.5$ to $2$ meters, or even lower (e.g., $0.1$ m) if your model is very fine.
* **Stretching Constant ($C$)**: Define a constant ratio $C > 1$ between adjacent cell heights to ensure high resolution near $z_1$. For example, setting $C=1.2$ means $z_{k+1} = 1.2 \cdot z_k$.
* **Implementation**: Use the formula $z_k = z_1 \cdot C^{k-1}$ to define the center heights of your limited layers.

### 2. Flux-Based Monin-Obukhov Length ($L$)

The surface fluxes are incorporated by using them to define the **Monin-Obukhov length ($L$)**, which dictates the scaling of $\zeta = z/L$ and is required for evaluating the stability functions $\phi_{m,h}$.

The definition of $L$ is:
$$
L = -\frac{u_*^3 \bar{\theta}_{v}}{\kappa g \overline{w'\theta'_v}}
$$

| Variable | Description | Source (Input/Calculation) |
| :--- | :--- | :--- |
| $\overline{w'\theta'_v}$ | Surface kinematic heat flux | **Model Input** (or derived from $\Delta\theta$) |
| $u_*$ | Friction velocity ($\sqrt{\tau/\rho}$) | **Model Input** (or derived from $\Delta U$) |
| $\kappa$ | Von Kármán constant ($\approx 0.4$) | Constant |
| $g$ | Gravity | Constant |
| $\bar{\theta}_v$ | Mean virtual temperature | **Model Input** |

You must supply $u_*$ and $\overline{w'\theta'_v}$ (or quantities that derive them, like surface wind stress $\tau$ and heat flux $H$) as **inputs** to calculate a single, constant $L$ for the entire diagnostic profile.

### 3. First-Layer MOST Closure

The $\phi$ functions and $Ri_g$ are evaluated layer-by-layer based on the calculated $L$.

| Layer | Calculation | Notes |
| :--- | :--- | :--- |
| **All Layers ($k=1$ to $N$)** | $\zeta_k = z_k / L$ | Scales the height using the flux-derived $L$. |
| | $\phi_{m}(\zeta_k), \phi_{h}(\zeta_k)$ | Evaluate the chosen stability profiles (BD\_PL, QSBL, etc.) |
| | $Ri_g(\zeta_k) = \zeta_k \frac{\phi_h(\zeta_k)}{\phi_m(\zeta_k)^2}$ | Evaluate the Richardson number at the layer center. |

Since $L$ is a constant for a given surface flux state, the entire $Ri_g(\zeta)$ curve is fixed, and you can then proceed with the **Finite Difference curvature calculation** across your few discrete layers.

**Key point:** The accuracy of the $Ri_g$ curvature on a few layers depends entirely on having a high-quality (log-stretched) grid and an accurate **flux input** to anchor $L$.

### Curvature Hook
For each stretched grid level compute $v_m,v_h$ (power-law analytic), then
\[
V_{\log}=v_h-2v_m,\ W_{\log}=V_{\log}',\ \partial_\zeta^2 Ri_g=F[2V_{\log}+\zeta(V_{\log}^2-W_{\log})].
\]
Map to height with constant or variable L logic for slope-adjusted stratification.