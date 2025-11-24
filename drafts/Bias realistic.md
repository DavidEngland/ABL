## 🔬 Derivation with Realistic $\phi$ Functions

### 1. The Realistic $\text{Ri}_g(\zeta)$

Using the log-linear $\phi$ functions, the gradient $\text{Richardson Number}$ is:

$$\text{Ri}_g(\zeta) = \frac{\zeta \cdot \phi_h(\zeta)}{\phi_m(\zeta)^2} = \frac{\zeta (1 + a_h \zeta)}{(1 + a_m \zeta)^2}$$

This is the form that accurately captures stability in many operational models using $\text{MOST}$. Typical coefficients are $\mathbf{a_m \approx 5}$ and $\mathbf{a_h \approx 7.8}$.

---

### 2. The Bias Ratio $B(\zeta)$

The bias ratio $B(\zeta)$ is still defined as:

$$B(\zeta) = \frac{\text{Ri}_g(\zeta)}{\text{Ri}_b(\zeta)} = \frac{\zeta^2 \cdot F(\zeta)}{\int_0^\zeta \xi \cdot F(\xi)\,d\xi}$$

Now we must calculate the integral of $\text{Ri}_g(\xi)$ over the layer:

$$\int_0^\zeta \text{Ri}_g(\xi)\,d\xi = \int_0^\zeta \frac{\xi (1 + a_h \xi)}{(1 + a_m \xi)^2}\,d\xi$$

This integral is much more complex but can be solved analytically by decomposition.

---

### 3. Solving the Integral $\int \text{Ri}_g(\xi)\,d\xi$

Let $u = 1 + a_m \xi$, so $\xi = \frac{u-1}{a_m}$ and $d\xi = \frac{du}{a_m}$.
Also, $1 + a_h \xi = 1 + a_h \frac{u-1}{a_m} = \frac{a_m + a_h u - a_h}{a_m}$.

The integral transforms to:
$$\frac{1}{a_m^2} \int \frac{(u-1)}{u^2} \left( \frac{a_h}{a_m} u + \frac{a_m - a_h}{a_m} \right)\,du$$

This results in a closed-form solution (after a lot of algebra, using $a_h$ and $a_m$):

$$\int_0^\zeta \text{Ri}_g(\xi)\,d\xi = \frac{1}{a_m^2} \left[ (a_h - a_m)\left( \ln(1+a_m\zeta) + \frac{a_m\zeta}{1+a_m\zeta} \right) + a_m \zeta \left( 1 + \frac{a_h \zeta}{1+a_m\zeta} \right) \right]$$

---

### 4. Final Realistic Bias Ratio $B(\zeta)$

With the complexity of the integral, the final bias ratio $B(\zeta)$ is:

$$B(\zeta) = \frac{\frac{\zeta (1 + a_h \zeta)}{(1 + a_m \zeta)^2}}{\frac{1}{\zeta \cdot a_m^2} \left[ (a_h - a_m)\left( \ln(1+a_m\zeta) + \frac{a_m\zeta}{1+a_m\zeta} \right) + a_m \zeta \left( 1 + \frac{a_h \zeta}{1+a_m\zeta} \right) \right]}$$

**Takeaway:** While the algebra is complex, the physical interpretation remains the same:
* The bias ratio $B$ is a **highly non-linear function** of stability $\zeta$.
* For realistic values ($a_h > 2a_m$), the function $\text{Ri}_g(\zeta)$ is concave-down, ensuring $\mathbf{B(\zeta) > 1}$ (the bias exists).
* The $\text{McNider}$ correction must use this complex $B(\zeta)$ in the exponent to accurately compensate for the error introduced by the log-linear $\phi$ functions used in the model closure.