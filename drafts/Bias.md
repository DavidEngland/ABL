## 📐 Bias ($B$) from $\text{Ri}_g$ as $\zeta \cdot F(\zeta)$

You are defining the gradient $\text{Richardson Number}$ ($\text{Ri}_g$) in terms of a general stability function $F(\zeta)$ (often related to $\phi_m$ and $\phi_h$ in $\text{MOST}$):

$$\text{Ri}_g(\zeta) = \zeta \cdot F(\zeta)$$

In the $\text{Pade}$ approximation example you provided earlier, the function was $F(\zeta) = \frac{1}{1 + \beta \zeta}$, so $\text{Ri}_g(\zeta) = \frac{\zeta}{1 + \beta \zeta}$.

### 1. The Bulk $\text{Richardson Number}$ ($\text{Ri}_b$)

The bulk $\text{Richardson Number}$ is the layer average of the gradient $\text{Ri}$ over the layer thickness $\Delta z$, which corresponds to a stability range of $[0, \zeta]$.

$$\text{Ri}_b(\zeta) \approx \frac{1}{\zeta} \int_0^\zeta \text{Ri}_g(\xi)\,d\xi = \frac{1}{\zeta} \int_0^\zeta \xi \cdot F(\xi)\,d\xi$$

### 2. The Bias Ratio ($B$)

The bias ratio $B$ is the ratio of the point (gradient) $\text{Ri}$ to the layer-averaged (bulk) $\text{Ri}$:

$$\mathbf{B(\zeta) = \frac{\text{Ri}_g(\zeta)}{\text{Ri}_b(\zeta)} = \frac{\zeta \cdot F(\zeta)}{\frac{1}{\zeta} \int_0^\zeta \xi \cdot F(\xi)\,d\xi}}$$

This simplifies the bias ratio to:

$$\mathbf{B(\zeta) = \frac{\zeta^2 \cdot F(\zeta)}{\int_0^\zeta \xi \cdot F(\xi)\,d\xi}}$$

---

## 3. The Pade Approximation Case

Let's apply this general formula to the specific **Pade Approximation** used in your homework example, where $F(\zeta) = \frac{1}{1 + \beta \zeta}$:

$$\text{Ri}_g(\zeta) = \frac{\zeta}{1 + \beta \zeta}$$

* **Numerator $(\text{Ri}_g(\zeta))$:**
    $$\text{Ri}_g(\zeta) = \frac{\zeta}{1 + \beta \zeta}$$

* **Denominator $(\text{Ri}_b(\zeta))$:**
    The integral is $\int_0^\zeta \frac{\xi}{1 + \beta \xi}\,d\xi = \frac{1}{\beta^2}\left(\beta \zeta - \ln(1+\beta \zeta)\right)$.
    $$\text{Ri}_b(\zeta) = \frac{1}{\zeta \beta^2}\left(\beta \zeta - \ln(1+\beta \zeta)\right)$$

* **Resulting Bias ($B$):**
    $$B(\zeta) = \frac{\frac{\zeta}{1 + \beta \zeta}}{\frac{1}{\zeta \beta^2}\left(\beta \zeta - \ln(1+\beta \zeta)\right)}$$

$$\mathbf{B(\zeta) = \frac{\zeta^2 \beta^2}{(1+\beta \zeta)(\beta \zeta - \ln(1+\beta \zeta))}}$$

**Conclusion:** This confirms that when $\text{Ri}_g(\zeta)$ has **concave-down curvature** (which occurs because $F(\zeta)$ decreases as $\zeta$ increases, making the $\text{Ri}_g$ curve flatten out), the layer integral ($\text{Ri}_b$) systematically **underestimates** the true point stability ($\text{Ri}_g$), thus $B>1$. The McNider correction factor is entirely dedicated to compensating for this numerical error.