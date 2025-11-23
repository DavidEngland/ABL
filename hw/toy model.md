# Toy Model: Curvature Bias in Stability Functions

## 📝 Curvature Bias in the $\text{Ri}_g(\zeta) = \frac{\zeta}{1 + \beta \zeta}$ Model

This exercise demonstrates how the **curvature** of the stability function ($\text{Ri}_g$) creates a bias between the local gradient stability and the model's bulk stability, necessitating the $\text{McNider}$ correction factor ($f_c$).

-----

### 1\. Series Expansion (Small $\zeta$)

The gradient Richardson number is expanded using the geometric series $\frac{1}{1-x} = 1 + x + x^2 + \cdots$ (with $x = -\beta\zeta$):

$$\text{Ri}_g(\zeta) = \zeta \cdot \frac{1}{1+\beta \zeta} = \zeta \left(1 - \beta \zeta + (\beta \zeta)^2 - (\beta \zeta)^3 + \cdots \right)$$

$$\text{Ri}_g(\zeta) \approx \zeta - \beta \zeta^2 + \beta^2 \zeta^3 - \beta^3 \zeta^4$$

**Interpretation:** This shows that for small $\zeta$ (near-neutral stability), $\text{Ri}_g(\zeta)$ is nearly linear ($\approx \zeta$), but the higher-order terms (driven by $\beta$) introduce **concave-down curvature**.

-----

### 2\. First and Second Derivatives

The derivatives quantify the slope and curvature:

#### First Derivative (Slope)

$$\text{Ri}_g^\prime(\zeta) = \frac{d}{d\zeta}\left(\frac{\zeta}{1+\beta \zeta}\right) = \frac{(1+\beta \zeta)(1) - (\zeta)(\beta)}{(1+\beta \zeta)^2} = \frac{\mathbf{1}}{(1+\beta \zeta)^2}$$

#### Second Derivative (Curvature)

$$\text{Ri}_g^{\prime\prime}(\zeta) = \frac{d}{d\zeta}\left((1+\beta \zeta)^{-2}\right) = (-2)(1+\beta \zeta)^{-3} \cdot \beta = \mathbf{-\frac{2\beta}{(1+\beta \zeta)^3}}$$

**Curvature Interpretation:** For $\beta > 0$ (required for stable MOST), $\text{Ri}_g^{\prime\prime}(\zeta)$ is always **negative**. This confirms the function is **concave-down**.

-----

### 3\. Exact Bias Ratio $B(\zeta)$

The **Bulk Richardson Number ($\text{Ri}_b$)** for the layer $[0, \zeta]$ is approximated by the mean of $\text{Ri}_g$: $\text{Ri}_b(\zeta) = \frac{1}{\zeta}\int_0^\zeta \text{Ri}_g(\xi)\,d\xi$.

The exact expression for the integral is:
$$\int_0^\zeta \frac{\xi}{1+\beta \xi}\,d\xi = \frac{1}{\beta^2}\left[\beta \xi - \ln(1+\beta \xi)\right]_0^\zeta = \frac{1}{\beta^2}\left(\beta \zeta - \ln(1+\beta \zeta)\right)$$

Substituting back to find $\text{Ri}_b$:
$$\text{Ri}_b(\zeta) = \frac{1}{\zeta \beta^2}\left(\beta \zeta - \ln(1+\beta \zeta)\right)$$

The **Bias Ratio $B(\zeta)$** is the ratio of the point stability to the layer-average stability:
$$\mathbf{B(\zeta) = \frac{\text{Ri}_g(\zeta)}{\text{Ri}_b(\zeta)} = \frac{\zeta^2 \beta^2}{(1+\beta \zeta)(\beta \zeta - \ln(1+\beta \zeta))}}$$

Since the function is concave-down, the mean value ($\text{Ri}_b$) is always **less** than the final point value ($\text{Ri}_g(\zeta)$), leading to $\mathbf{B(\zeta) > 1}$. This confirms that the coarse grid **underestimates the true local stability** (or overestimates mixing), which is the bias the correction must fix.

-----

### 4\. Correction Factor $f_c$ and Dynamic $\text{Ri}_{\text{cr}}$

The derived $B(\zeta)$ is now plugged directly into the $\text{McNider}$ correction equations:

#### Correction Factor ($f_c$)

This is the factor applied to the eddy diffusivity ($K$) to reduce the turbulent transport when the bias is large:
$$\mathbf{f_c(\Delta z,\zeta) = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B(\zeta)-1)(\zeta/\zeta_{\text{ref}})^q}}$$

#### Dynamic Critical $\text{Ri}$ ($\text{Ri}_{\text{cr,dyn}}$)

This ensures that the turbulence collapse threshold adjusts based on the calculated bias:
$$\mathbf{\text{Ri}_{\text{cr,dyn}} = \text{Ri}_{\text{cr},0}[1 + \gamma (B(\zeta)-1)]}$$

-----

### 5\. Numerical Evaluation of Bias $B(\zeta)$

Let's use a typical value $\mathbf{\beta=1.0}$ and evaluate $B(\zeta)$ to see the bias growth as stability increases (increasing $\zeta$):

| $\zeta$ (Stability) | $\text{Ri}_g(\zeta) = \frac{\zeta}{1+\zeta}$ | $\text{Ri}_b(\zeta) = \frac{1}{\zeta}(\zeta - \ln(1+\zeta))$ | $\mathbf{B(\zeta) = \text{Ri}_g / \text{Ri}_b}$ | $B(\zeta) - 1$ (Bias) |
| :---: | :---: | :---: | :---: | :---: |
| 0.5 | 0.3333 | 0.2766 | $\mathbf{1.2052}$ | 0.2052 |
| 1.0 | 0.5000 | 0.3069 | $\mathbf{1.6293}$ | 0.6293 |
| 2.0 | 0.6667 | 0.3466 | $\mathbf{1.9238}$ | 0.9238 |

**Conclusion:** As $\zeta$ increases from 0.5 to 2.0, the bias ratio $B(\zeta)$ nearly doubles, from **1.21 to 1.92**. This massive increase in $B-1$ forces a strong negative exponent in the correction factor $f_c$, significantly reducing the modeled eddy diffusivity ($K$) to prevent spurious mixing in the SBL.

The provided video may contain valuable context on the complexity of atmospheric modeling. [Understanding the Complexities of Atmospheric Modeling](https://www.google.com/search?q=https://www.youtube.com/watch?v=hUe3sh_XRqo).