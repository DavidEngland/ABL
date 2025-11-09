Hello! I see you have a detailed set of notes and established procedures for working with the **Neutral Log Wind Profile**, the **Geometric Mean Height**, and their implications for estimating **Surface Roughness** ($z_0$), **Friction Velocity** ($u_*$), and the **Drag Coefficient** ($C_D$), extending into **Atmospheric Stability** corrections using a Polylogarithm/Power-law formulation.

Based on these comprehensive notes, I can answer specific questions, generate a practical workflow, or elaborate on the underlying principles.

To provide an accurate and authoritative response, especially concerning the underlying physics and the standard values of constants, I should perform a quick search to ensure the most current academic consensus aligns with the document's contents.

Here are the search queries I will use:
Yes, I can certainly elaborate on the implications of the **Geometric Mean Height** and the proposed **Polylogarithm Formulation** for atmospheric stability, based on the comprehensive framework you provided.

The key scientific principles outlined in your notes are consistent with established **Monin-Obukhov Similarity Theory (MOST)**, which is foundational in micrometeorology and boundary-layer modeling (Businger et al., 1971; Gryanik et al., 2020).

---

## 🌬️ Geometric Mean Height and Velocity Averaging

The notes correctly highlight a crucial point for multi-level measurements under **neutral stratification**:

* **The Arithmetic Mean of Velocities is a Log-Profile Velocity:** The arithmetic mean of mean horizontal wind speeds, $\bar{U}$, measured at heights $z_1, \dots, z_N$ is equal to the log-profile velocity evaluated at the **geometric mean height**, $z_g$ (Section 2).
    $$
    \bar{U} = \frac{1}{N}\sum_{i=1}^N U(z_i) = U(z_g), \quad \text{where} \quad z_g = \Big(\prod_{i=1}^N z_i\Big)^{1/N}
    $$
* **Bias in Drag Coefficient ($C_D$)**: Using the **arithmetic mean height** ($\bar{z}$) instead of the geometric mean height ($z_g$) when calculating the drag coefficient will systematically **underestimate** $C_D$, because $\bar{z} \ge z_g$ (Section 7). This is a significant source of structural bias, particularly when the measurement heights span a large range.

### Key Implication: Data Processing
To avoid bias, the best practice is to **not** pre-average the velocities and heights, but instead to:
1.  **Regress** individual $\mathbf{U}$ against $\mathbf{\ln(z-d)}$ to extract $u_*$ and $z_0$ (Sections 3, 4).
2.  **Calculate** the desired value (e.g., $C_D$) at a specific reference height ($z_{\mathrm{ref}}$), or at the calculated $z_g$ if a single, representative value for the averaged speed is required (Section 7).

---

## 🌡️ Polylogarithm Stability Function

Your proposed stability modifier, derived from the first-order polylogarithm ($\mathrm{Li}_1(x)$), provides an analytically flexible formulation that directly maps to the widely used power-law forms:

$$
\phi\left(\tfrac{z}{L}\right)=\exp\left(\alpha \,\mathrm{Li}_{1}\left(\beta\,\tfrac{z}{L}\right)\right)
$$

Since $\mathrm{Li}_1(x) = -\ln(1-x)$, the expression simplifies to the power-law form used in many models:
$$
\phi\left(\tfrac{z}{L}\right) = (1-\beta z/L)^{-\alpha} \quad \text{where } \zeta = z/L
$$

This formulation is consistent with the classic **Businger-Dyer** framework, where the stability correction function for momentum ($\phi_m$) is often written as a power law for unstable conditions ($\zeta < 0$): $\phi_m(\zeta) = (1-c\zeta)^{-p}$ (Li et al., 2012).

| Parameter | Interpretation in $\phi=(1-\beta z/L)^{-\alpha}$ |
| :--- | :--- |
| **$\alpha$** | Controls the **curvature** and strength of the stability correction. |
| **$\beta$** | Sets the effective **scaling** of the stability parameter $z/L$. |

### Connection to Richardson Number ($Ri$) Bias

The quadratic dependence of the Bulk Richardson number ($Ri$) on the inverse log term means that errors in the roughness parameters ($z_0, d$) or the omission of the stability correction ($\psi_m$) are significantly amplified (Section 10).

$$
Ri_{1/2} \propto \frac{1}{\left[\ln\frac{z_2-d}{z_0}-\psi_m(z_2/L)\right]^2} \propto \frac{1}{(C_D^{-1/2})^2}
$$

* **Bias Effect:** If stability corrections ($\psi_m$) are incorrectly assumed to be zero (i.e., neutral), the denominator term, $\left[\ln\frac{z_2-d}{z_0}\right]^2$, is typically *underestimated* for stable conditions, leading to an **overestimation of $Ri_{1/2}$** (Section 10).

---

## 📘 Von Kármán Constant ($\kappa$)

The value of the **von Kármán constant** ($\kappa$), which links the velocity profile slope to the friction velocity, is a critical parameter.

* **Consensus Value:** The notes use $\kappa \approx 0.40$, which is consistent with the values found in numerous micrometeorological studies, though values such as $\mathbf{\kappa = 0.41}$ are also commonly used and considered a consensus value in contemporary research (Gryanik et al., 2020; Kent et al., 2017).

The notes correctly state that $\kappa$ is an **empirical constant** whose value is still subject to small variations due to factors like sensor tilt, flow distortion, and the specific dataset used (Section 12). There is no accepted theoretical link between $\kappa$ and the Euler–Mascheroni constant ($\gamma$)—the proximity is purely incidental.

---

### References
Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971). Flux-profile relationships in the atmospheric surface layer. *Journal of the Atmospheric Sciences*, *28*(2), 181–189. [DOI](https://doi.org/10.1175/1520-0469(1971)028<0181:fprita>2.0.co;2)

Gryanik, V. M., Lüpkes, C., Grachev, A., & Sidorenko, D. (2020). New modified and extended stability functions for the stable boundary layer based on SHEBA and parametrizations of bulk transfer coefficients for climate models. *Journal of the Atmospheric Sciences*, *77*(7), 2687–2716. [DOI](https://doi.org/10.1175/jas-d-19-0255.1)

Kent, C. W., Grimmond, S., Barlow, J., Gatey, D., Kotthaus, S., Lindberg, F., & Halios, C. H. (2017). Evaluation of urban local-scale aerodynamic parameters: Implications for the vertical profile of wind speed and for source areas. *Boundary-Layer Meteorology*, *164*(1), 183–213. [DOI](https://doi.org/10.1007/s10546-017-0248-z)

Li, D., Katul, G. G., & Bou-Zeid, E. (2012). Mean velocity and temperature profiles in a sheared diabatic turbulent boundary layer. *Physics of Fluids*, *24*(10). [DOI](https://doi.org/10.1063/1.4757660)

---
This is a fascinating and fundamental problem, Dick. It sounds like you're aiming to improve the representation of the atmospheric surface layer (ASL) within a numerical model by using **Monin-Obukhov Similarity Theory (MOST)** to provide a more accurate boundary condition for the Bulk Richardson number ($Ri$), specifically focusing on its **curvature** near the surface.

To answer your question directly: **Yes, we absolutely can use the similarity forms to analytically derive the profiles of the gradient Richardson number ($Ri_g$) and, by extension, its curvature, $\partial^2 Ri_g / \partial z^2$.**

Using the analytical similarity functions is the most robust way to estimate these derivatives, especially when the coarse vertical resolution of the model grid (as hinted by your Figure 10 concern) fails to capture the sharp gradients near $z_0$.

---

## 📐 Deriving $Ri_g$ Curvature from MOST

The key to your approach is to use the known analytic forms for the derivatives of wind speed ($U$) and potential temperature ($\Theta$) from MOST.

### 1. The Gradient Richardson Number ($Ri_g$)

The gradient Richardson number is defined as:
$$
Ri_g(z) = \frac{g}{\Theta_v} \frac{\partial \Theta_v / \partial z}{(\partial U / \partial z)^2}
$$

In MOST, the vertical derivatives (gradients) are defined using the stability functions for momentum ($\phi_m$) and heat ($\phi_h$):

* **Momentum Gradient:** $\frac{\partial U}{\partial z} = \frac{u_*}{\kappa z} \phi_m(\zeta)$, where $\zeta = z/L$.
* **Heat Gradient:** $\frac{\partial \Theta_v}{\partial z} = \frac{\theta_*}{\kappa z} \phi_h(\zeta)$.

Substituting these into the $Ri_g$ definition, and noting that the mean virtual potential temperature $\Theta_v$ is approximately constant, we get the fundamental MOST closure:
$$
Ri_g(\zeta) = \frac{g \theta_* / \Theta_v}{u_*^2 / (\kappa z)} \frac{\phi_h(\zeta) / (\kappa z)}{\left( \frac{\phi_m(\zeta)}{\kappa z} \right)^2} \frac{1}{\phi_m^2(\zeta)} \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)}
$$
The dimensionless relationship is often simplified by the assumption that the ratio $Ri_g / \zeta" is directly related to the stability functions:
$$
Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)}
$$
This is the expression you need for the first derivative of the profile.

---

### 2. Calculating Curvature ($\partial^2 Ri_g / \partial z^2$)

To find the curvature, we need the second derivative of $Ri_g$ with respect to height $z$. It is mathematically cleaner to first find $\frac{\partial Ri_g}{\partial \zeta}$ and then use the chain rule:

$$\frac{\partial Ri_g}{\partial z} = \frac{\partial Ri_g}{\partial \zeta} \cdot \frac{\partial \zeta}{\partial z} = \frac{1}{L} \frac{\partial Ri_g}{\partial \zeta}$$

And the curvature is:
$$
\frac{\partial^2 Ri_g}{\partial z^2} = \frac{\partial}{\partial z} \left( \frac{1}{L} \frac{\partial Ri_g}{\partial \zeta} \right) = \frac{1}{L^2} \frac{\partial^2 Ri_g}{\partial \zeta^2}
$$

Since $L$ is constant in the ASL, the problem reduces to finding $\frac{\partial^2 Ri_g}{\partial \zeta^2}$.

Substituting $Ri_g(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)}$, let $F(\zeta) = \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)}$.
$$
Ri_g(\zeta) = \zeta F(\zeta)
$$

The derivatives are:
$$
\frac{\partial Ri_g}{\partial \zeta} = F(\zeta) + \zeta \frac{\partial F}{\partial \zeta}
$$
$$
\frac{\partial^2 Ri_g}{\partial \zeta^2} = 2 \frac{\partial F}{\partial \zeta} + \zeta \frac{\partial^2 F}{\partial \zeta^2}
$$

This requires two derivatives of the ratio $F(\zeta) = \phi_h(\zeta) / \phi_m^2(\zeta)$. Since $\phi_m$ and $\phi_h$ are standard similarity functions (e.g., Businger-Dyer or your Polylog forms), their derivatives are known.

### 3. Using Your Exponential/Polylog Approximation

You mentioned using an exponential to approximate the stability function. Let's assume you're approximating the $\phi$ functions themselves with a simple exponential form **near neutrality** ($\zeta \approx 0$):

$$\phi(\zeta) \approx e^{a \zeta} \quad (\text{Approximation for small } |\zeta|)$$

If you use the standard, more rigorous power-law forms you developed (which is better than a simple exponential):
$$\phi_m(\zeta) = (1-\beta_m \zeta)^{-\alpha_m}$$
$$\phi_h(\zeta) = (1-\beta_h \zeta)^{-\alpha_h}$$

These forms are highly convenient because they are analytic and easy to differentiate. You can perform the chain and quotient rule operations on $F(\zeta) = \phi_h(\zeta) / \phi_m^2(\zeta)$ to obtain the exact expression for the curvature $\partial^2 Ri_g / \partial \zeta^2$.

---

## 💡 Practical Recommendation for Modeling

Your approach is sound and provides a more physically grounded solution than a finite-difference approximation on a coarse grid.

1. **Analytic Curvature Calculation:** Select your best-fit **analytic MOST forms** for $\phi_m$ and $\phi_h$ (e.g., the $(\mathbf{1-\beta \zeta)^{-\alpha}}$ form).
2. **Derive the Curvature Formula:** Analytically compute $\frac{\partial^2 Ri_g}{\partial \zeta^2}$ using the derivatives of your chosen $\phi$ functions. This will yield a complex but exact algebraic expression in terms of $\zeta$, $\alpha_m$, $\beta_m$, $\alpha_h$, and $\beta_h$.
3. **Implement in Model:** In the model's ASL parameterization (below the first model level, $z_1$), use the derived $\frac{\partial^2 Ri_g}{\partial z^2} = \frac{1}{L^2} \frac{\partial^2 Ri_g}{\partial \zeta^2}$ as the estimate for the $Ri$ curvature term in your correction, effectively providing a **similarity-based boundary condition** for the curvature.

This approach will smoothly connect the surface layer dynamics (where the curvature is dominated by the MOST structure) to the lowest model level, mitigating the coarse-grid numerical error near $z_0$.

**Next Step:** Would you like me to help you set up the first derivative, $\partial F / \partial \zeta$, using the $\phi(\zeta) = (1-\beta \zeta)^{-\alpha}$ form so you can tackle the second derivative?

That's an excellent challenge! Let's tackle the second derivative, $\partial^2 Ri_g / \partial \zeta^2$, using your recommended power-law stability functions. This analytical approach provides the necessary **curvature** for your model's boundary condition near the surface.

We use the general power-law forms derived from your Polylog formulation:

* Momentum: $\phi_m(\zeta) = (1-\beta_m \zeta)^{-\alpha_m}$
* Heat: $\phi_h(\zeta) = (1-\beta_h \zeta)^{-\alpha_h}$
* Gradient Richardson Number: $Ri_g(\zeta) = \zeta F(\zeta)$, where $F(\zeta) = \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)}$.

---

## 1. Defining the Function $F(\zeta)$

Let's first define $F(\zeta)$ and its component parts:
$$
F(\zeta) = \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)} = (1-\beta_h \zeta)^{-\alpha_h} (1-\beta_m \zeta)^{2\alpha_m}
$$

To simplify the notation for the derivatives, let:

* $A = -\alpha_h$ and $B = \beta_h$
* $C = 2\alpha_m$ and $D = \beta_m$

Then:
$$
F(\zeta) = (1-B\zeta)^A (1-D\zeta)^C
$$

---

## 2. First Derivative of $F(\zeta)$ ($\partial F / \partial \zeta$)

Using the product and chain rules: $\frac{d}{d\zeta} (u v) = u'v + uv'$.

* $\frac{d}{d\zeta} (1-B\zeta)^A = A(1-B\zeta)^{A-1} (-B) = -AB (1-B\zeta)^{A-1}$
* $\frac{d}{d\zeta} (1-D\zeta)^C = C(1-D\zeta)^{C-1} (-D) = -CD (1-D\zeta)^{C-1}$

$$
\frac{\partial F}{\partial \zeta} = \underbrace{\left[-AB (1-B\zeta)^{A-1}\right] (1-D\zeta)^C}_{u'v} + \underbrace{(1-B\zeta)^A \left[-CD (1-D\zeta)^{C-1}\right]}_{uv'}
$$

Factor out $F(\zeta) = (1-B\zeta)^A (1-D\zeta)^C$:

$$
\frac{\partial F}{\partial \zeta} = F(\zeta) \left[ \frac{-AB}{1-B\zeta} - \frac{CD}{1-D\zeta} \right]
$$

---

## 3. Second Derivative of $F(\zeta)$ ($\partial^2 F / \partial \zeta^2$)

Now we take the derivative of the expression for $\partial F / \partial \zeta$. Again, using the product rule: $\frac{d}{d\zeta} (u v) = u'v + uv'$, where $u = F(\zeta)$ and $v = \left[ \dots \right]$.

$$
\frac{\partial^2 F}{\partial \zeta^2} = \underbrace{\left(\frac{\partial F}{\partial \zeta}\right) v}_{u'v} + \underbrace{F(\zeta) \left(\frac{\partial v}{\partial \zeta}\right)}_{uv'}
$$

### A. Calculate $\partial v / \partial \zeta$

Let $v(\zeta) = \left[ \frac{-AB}{1-B\zeta} - \frac{CD}{1-D\zeta} \right]$.
The derivative of $1/(1-B\zeta)$ is $-(-B)/(1-B\zeta)^2 = B/(1-B\zeta)^2$.

$$
\frac{\partial v}{\partial \zeta} = \frac{d}{d\zeta} \left[ -AB(1-B\zeta)^{-1} - CD(1-D\zeta)^{-1} \right]
$$

$$
\frac{\partial v}{\partial \zeta} = (-AB)\left[ -1(1-B\zeta)^{-2}(-B) \right] - (CD)\left[ -1(1-D\zeta)^{-2}(-D) \right]
$$

$$
\frac{\partial v}{\partial \zeta} = -A B^2 (1-B\zeta)^{-2} - C D^2 (1-D\zeta)^{-2}
$$

### B. Combine Terms

Substituting $u'v$ and $uv'$ back into the product rule:

$$
\frac{\partial^2 F}{\partial \zeta^2} = \underbrace{F(\zeta) \left[ \frac{-AB}{1-B\zeta} - \frac{CD}{1-D\zeta} \right]}_{u'} \underbrace{\left[ \frac{-AB}{1-B\zeta} - \frac{CD}{1-D\zeta} \right]}_{v} + \underbrace{F(\zeta)}_{u} \underbrace{\left[ \frac{-A B^2}{(1-B\zeta)^2} - \frac{C D^2}{(1-D\zeta)^2} \right]}_{v'}
$$

### C. Final Expression for $\partial^2 F / \partial \zeta^2$

The **final expression** for the second derivative is:

$$
\frac{\partial^2 F}{\partial \zeta^2} = F(\zeta) \left\{ \left[ \frac{-AB}{1-B\zeta} - \frac{CD}{1-D\zeta} \right]^2 - \left[ \frac{A B^2}{(1-B\zeta)^2} + \frac{C D^2}{(1-D\zeta)^2} \right] \right\}
$$

---

## 4. Curvature of $Ri_g$

Now we use the derivative relations for $Ri_g$:
$$
\frac{\partial^2 Ri_g}{\partial \zeta^2} = 2 \frac{\partial F}{\partial \zeta} + \zeta \frac{\partial^2 F}{\partial \zeta^2}
$$

Substitute the full expressions for $\frac{\partial F}{\partial \zeta}$ (Section 2) and $\frac{\partial^2 F}{\partial \zeta^2}$ (Section 3) into this equation, and remember that **the curvature of $Ri_g$ with respect to height $z$** is:
$$
\frac{\partial^2 Ri_g}{\partial z^2} = \frac{1}{L^2} \frac{\partial^2 Ri_g}{\partial \zeta^2}
$$

### Summary of the Curvature Term

The final analytical expression for the curvature is complex, but it is **exact** for your chosen $\phi$ functions and is the required term for your correction:

$$
\frac{\partial^2 Ri_g}{\partial \zeta^2} = F(\zeta) \left\{ 2 \cdot \left[ \frac{-AB}{1-B\zeta} - \frac{CD}{1-D\zeta} \right] + \zeta \left[ \left( \frac{-AB}{1-B\zeta} - \frac{CD}{1-D\zeta} \right)^2 - \left( \frac{A B^2}{(1-B\zeta)^2} + \frac{C D^2}{(1-D\zeta)^2} \right) \right] \right\}
$$

where:

* $A = -\alpha_h$
* $B = \beta_h$
* $C = 2\alpha_m$
* $D = \beta_m$
* $F(\zeta) = (1-\beta_h \zeta)^{-\alpha_h} (1-\beta_m \zeta)^{2\alpha_m}$
* $\zeta = z/L$

This derived analytical form should provide a far better estimate of $Ri$ curvature near the surface than finite differences on your coarse grid, as it accurately reflects the asymptotic behavior dictated by MOST.

**Next Step:** I can perform a check on the neutral limit ($\zeta \to 0$) for this expression to ensure the curvature correctly tends towards a simple value as you approach the surface. Would you like to review the neutral limit?

Absolutely. Using the definitions $A = -\alpha_h$, $B = \beta_h$, $C = 2\alpha_m$ and $D = \beta_m$, here is the analytical expression for the curvature of the gradient Richardson number, $\partial^2 Ri_g / \partial \zeta^2$, solely in terms of the $\alpha$ and $\beta$ stability parameters and $\zeta = z/L$.

---

## Curvature of $Ri_g$ in Terms of $\alpha$ and $\beta$

Recall the core relationship for $Ri_g$ and the general function $F(\zeta)$:
$$
Ri_g(\zeta) = \zeta F(\zeta), \quad \frac{\partial^2 Ri_g}{\partial \zeta^2} = 2 \frac{\partial F}{\partial \zeta} + \zeta \frac{\partial^2 F}{\partial \zeta^2}
$$

### 1. The Function $F(\zeta)$

$F(\zeta)$ is the ratio of stability functions:
$$
F(\zeta) = \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)} = (1-\beta_h \zeta)^{-\alpha_h} (1-\beta_m \zeta)^{2\alpha_m}
$$

### 2. First Derivative Term ($\partial F / \partial \zeta$)

$$
\frac{\partial F}{\partial \zeta} = F(\zeta) \left[ \frac{\alpha_h \beta_h}{1-\beta_h \zeta} - \frac{2\alpha_m \beta_m}{1-\beta_m \zeta} \right]
$$

### 3. Second Derivative Term ($\partial^2 F / \partial \zeta^2$)

$$
\frac{\partial^2 F}{\partial \zeta^2} = F(\zeta) \left\{ \left[ \frac{\alpha_h \beta_h}{1-\beta_h \zeta} - \frac{2\alpha_m \beta_m}{1-\beta_m \zeta} \right]^2 - \left[ \frac{\alpha_h \beta_h^2}{(1-\beta_h \zeta)^2} - \frac{2\alpha_m \beta_m^2}{(1-\beta_m \zeta)^2} \right] \right\}
$$

### 4. Final Curvature Expression ($\partial^2 Ri_g / \partial \zeta^2$)

Substituting the derivative terms back into the expression for $\partial^2 Ri_g / \partial \zeta^2$:

$$
\frac{\partial^2 Ri_g}{\partial \zeta^2} = F(\zeta) \left\{ 2 \cdot \left[ \frac{\alpha_h \beta_h}{1-\beta_h \zeta} - \frac{2\alpha_m \beta_m}{1-\beta_m \zeta} \right] + \zeta \left[ \left( \frac{\alpha_h \beta_h}{1-\beta_h \zeta} - \frac{2\alpha_m \beta_m}{1-\beta_m \zeta} \right)^2 - \left( \frac{\alpha_h \beta_h^2}{(1-\beta_h \zeta)^2} - \frac{2\alpha_m \beta_m^2}{(1-\beta_m \zeta)^2} \right) \right] \right\}
$$

where $F(\zeta)$ is defined above.

This formula provides the exact, analytical curvature necessary for your model's boundary condition, where the height curvature is related by $\frac{\partial^2 Ri_g}{\partial z^2} = \frac{1}{L^2} \frac{\partial^2 Ri_g}{\partial \zeta^2}$.

That's an excellent point for improving computational efficiency and clarity! You are absolutely right; distributing the $F(\zeta)$ terms and simplifying the resulting expression for the curvature $\partial^2 Ri_g / \partial \zeta^2$ will produce a **fully self-contained equation** that directly depends on the $\alpha$'s, $\beta$'s, and $\zeta$.

Here is the final, expanded, and simplified curvature expression.

---

## 🔬 Physical Interpretation of $Ri_g$ Curvature

The curvature term, $\frac{\partial^2 Ri_g}{\partial z^2}$, determines how rapidly the non-linear stability effects (i.e., the *change* in the stability correction) take hold as height increases. Near the surface, the curvature is dominated by the shape of the $\phi$ functions.

The general relationship is:
$$
\frac{\partial^2 Ri_g}{\partial z^2} = \frac{1}{L^2} \frac{\partial^2 Ri_g}{\partial \zeta^2}
$$
The dimensionless curvature $\frac{\partial^2 Ri_g}{\partial \zeta^2}$ is a function of $\zeta$ and the stability parameters $\alpha$ and $\beta$.

---

## 1. Structure Appealing to Atmospheric Scientists

We can express the result using the logarithmic derivatives of the stability functions, $\phi_m$ and $\phi_h$, which are more physically intuitive than the absolute derivatives of the product $F(\zeta)$.

Let:
$$
\Phi_{m}(\zeta) = \frac{\zeta}{\phi_m}\frac{d\phi_m}{d\zeta} \quad \text{and} \quad \Phi_{h}(\zeta) = \frac{\zeta}{\phi_h}\frac{d\phi_h}{d\zeta}
$$
(These are common ways to represent the stability corrections in terms of their sensitivity to $\zeta$).

For your power-law forms, $\phi(\zeta)=(1-\beta\zeta)^{-\alpha}$:
$$
\frac{d\phi}{d\zeta} = -\alpha(1-\beta\zeta)^{-\alpha-1}(-\beta) = \frac{\alpha\beta}{1-\beta\zeta} \phi(\zeta)
$$
Thus, the logarithmic derivative is:
$$
\frac{1}{\phi}\frac{d\phi}{d\zeta} = \frac{\alpha\beta}{1-\beta\zeta}
$$

The full structure of the curvature is then:
$$
\frac{\partial^2 Ri_g}{\partial \zeta^2} = \frac{1}{\zeta} Ri_g(\zeta) \cdot \left[ \mathbf{Term\,A} + \mathbf{Term\,B} + \mathbf{Term\,C} \right]
$$

### $\mathbf{Term\,A}$: Contribution from $\partial F / \partial \zeta$

This term relates to the initial change in the ratio $\phi_h/\phi_m^2$ and represents the immediate non-linearity:
$$
\mathbf{Term\,A} = 2 \left( \frac{1}{\phi_h} \frac{d\phi_h}{d\zeta} - \frac{2}{\phi_m} \frac{d\phi_m}{d\zeta} \right)
= 2 \left( \frac{\alpha_h \beta_h}{1-\beta_h \zeta} - \frac{2\alpha_m \beta_m}{1-\beta_m \zeta} \right)
$$

### $\mathbf{Term\,B}$: Contribution from $(\partial F / \partial \zeta)^2$

This is the term arising from the square of the logarithmic derivative of $F(\zeta)$, scaled by $\zeta$. It represents how the *rate of change* itself scales non-linearly.
$$
\mathbf{Term\,B} = \zeta \left( \frac{1}{\phi_h} \frac{d\phi_h}{d\zeta} - \frac{2}{\phi_m} \frac{d\phi_m}{d\zeta} \right)^2
= \zeta \left( \frac{\alpha_h \beta_h}{1-\beta_h \zeta} - \frac{2\alpha_m \beta_m}{1-\beta_m \zeta} \right)^2
$$

### $\mathbf{Term\,C}$: Contribution from $\partial^2 F / \partial \zeta^2$ (Intrinsic Curvature)

This term isolates the change in the logarithmic derivative itself (the intrinsic curvature of $\phi_h/\phi_m^2$).

Since $\frac{d}{d\zeta}\left(\frac{1}{\phi}\frac{d\phi}{d\zeta}\right) = \frac{\alpha\beta^2}{(1-\beta\zeta)^2}$, Term C is:
$$
\mathbf{Term\,C} = -\zeta \left( \frac{\alpha_h \beta_h^2}{(1-\beta_h \zeta)^2} - \frac{2\alpha_m \beta_m^2}{(1-\beta_m \zeta)^2} \right)
$$

### The Final, Physically Intuitive Form

The most appealing form presents the result as the ratio of rational functions, which highlights the singular points at $1/\beta_h$ and $1/\beta_m$.

Let:
$$
V_{\text{log}} = \frac{\alpha_h \beta_h}{1-\beta_h \zeta} - \frac{2\alpha_m \beta_m}{1-\beta_m \zeta}
$$
$$
W_{\text{log}} = \frac{\alpha_h \beta_h^2}{(1-\beta_h \zeta)^2} - \frac{2\alpha_m \beta_m^2}{(1-\beta_m \zeta)^2}
$$
The curvature is proportional to the overall stability function $F(\zeta)$ multiplied by a term that corrects the "linear" $\zeta$ growth:

$$
\frac{\partial^2 Ri_g}{\partial \zeta^2} = F(\zeta) \cdot \Biggl\{ 2 V_{\text{log}} + \zeta V_{\text{log}}^2 - \zeta W_{\text{log}} \Biggr\}
$$

---

## 2. The Neutral Limit ($\zeta \to 0$)

The behavior at neutrality is a great check, as $\phi_m(\zeta) \to 1$ and $\phi_h(\zeta) \to 1$.

As $\zeta \to 0$:

* $F(\zeta) \to (1-0)^{-\alpha_h}(1-0)^{2\alpha_m} = 1$.
* $V_{\text{log}} \to (\alpha_h \beta_h - 2\alpha_m \beta_m)$.
* $W_{\text{log}} \to (\alpha_h \beta_h^2 - 2\alpha_m \beta_m^2)$.

The dimensionless curvature at neutrality becomes:
$$
\left.\frac{\partial^2 Ri_g}{\partial \zeta^2}\right|_{\zeta=0} = 2 V_{\text{log}} + 0 \cdot [V_{\text{log}}^2 - W_{\text{log}}]
$$
$$
\left.\frac{\partial^2 Ri_g}{\partial \zeta^2}\right|_{\zeta=0} = 2 (\alpha_h \beta_h - 2\alpha_m \beta_m)
$$

This is a powerful result: The curvature at the surface ($z=0$) depends only on the **linear coefficients** of the stability function expansions (where $\phi \approx 1 + \alpha\beta \zeta + \dots$).

* If $h$ and $m$ are similar ($\alpha_h \beta_h \approx 2\alpha_m \beta_m$), the curvature is near zero (i.e., $Ri_g$ grows linearly with $\zeta$ near the surface).
* If **heat transfer is less corrected** than momentum ($\alpha_h \beta_h < 2\alpha_m \beta_m$), $\left.\frac{\partial^2 Ri_g}{\partial \zeta^2}\right|_{\zeta=0} < 0$, meaning the $Ri_g$ profile is **concave down** near the surface.
* If **heat transfer is more corrected** than momentum ($\alpha_h \beta_h > 2\alpha_m \beta_m$), $\left.\frac{\partial^2 Ri_g}{\partial \zeta^2}\right|_{\zeta=0} > 0$, meaning the $Ri_g$ profile is **concave up** near the surface.

This shows that the curvature near the surface is governed entirely by the relative departure from neutrality of the heat and momentum exchange coefficients. This confirms the physical interpretation of your analytical approach.

---

## Ri-based shear functions: f_m(Ri) and f_h(Ri)

Goal: work directly with shear functions depending on the gradient Richardson number Ri, instead of MOST functions depending on ζ=z/L.

MOST identity:
\[
Ri_g(\zeta)=\zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2},\qquad
\phi_m(\zeta)=(1-\beta_m\zeta)^{-\alpha_m},\ \ \phi_h(\zeta)=(1-\beta_h\zeta)^{-\alpha_h}.
\]
Define Ri-based shear functions by composition:
\[
f_m(Ri):=\phi_m(\zeta(Ri)),\qquad f_h(Ri):=\phi_h(\zeta(Ri)).
\]

Critical Richardson number Ri_c: pick ζ_c (onset threshold) and set
\[
Ri_c=\zeta_c\,\frac{\phi_h(\zeta_c)}{\phi_m(\zeta_c)^2}.
\]
Use s:=Ri/Ri_c∈[0,1) as the nondimensional stability ratio.

Equal-β closed forms (β_m=β_h=β):
Let r:=\frac{F(\zeta)}{F(\zeta_c)}=\frac{Ri/\zeta}{Ri_c/\zeta_c}, with
\(
F(\zeta)=\phi_h/\phi_m^2=(1-\beta\zeta)^{-(\alpha_h-2\alpha_m)}.
\)
Then
\[
\frac{F(\zeta)}{F(\zeta_c)}=\left(\frac{1-\beta\zeta}{1-\beta\zeta_c}\right)^{-(\alpha_h-2\alpha_m)}
\Rightarrow
1-\beta\zeta=(1-\beta\zeta_c)\,r^{1/(\,2\alpha_m-\alpha_h\,)}.
\]
Hence Ri-based shear functions are pure power laws in s (up to constants):
\[
f_m(Ri)=C_m\, s^{-e_m},\quad f_h(Ri)=C_h\, s^{-e_h},
\]
with exponents
\[
e_m=\frac{\alpha_m}{2\alpha_m-\alpha_h},\qquad e_h=\frac{\alpha_h}{2\alpha_m-\alpha_h},
\]
and constants
\[
C_m=(1-\beta\zeta_c)^{-\alpha_m},\qquad C_h=(1-\beta\zeta_c)^{-\alpha_h}.
\]
Optional normalization: use \(\tilde f_{m,h}(s)=s^{-e_{m,h}}\) and absorb C_{m,h} into calibration.

General case (β_m≠β_h): first-order linearized mapping
Let \(K=2\alpha_m\beta_m-\alpha_h\beta_h\). For small |ζ| (near-neutral),
\[
F(\zeta)\approx 1+K\zeta,\quad \frac{Ri}{Ri_c}\approx \frac{1+K\zeta}{1+K\zeta_c}
\Rightarrow
\zeta \approx \frac{Ri/Ri_c-1}{K}+\frac{\zeta_c(1-Ri/Ri_c)}{1+K\zeta_c}.
\]
Then
\[
f_m(Ri)\approx \left(1-\beta_m \zeta\right)^{-\alpha_m},\quad
f_h(Ri)\approx \left(1-\beta_h \zeta\right)^{-\alpha_h},
\]
reject if brackets ≤0 (out of domain).

Limits and masking
- Neutral: Ri→0 ⇒ s→0, f_{m,h}→1 (use constants C_{m,h} if not normalized).
- Critical: Ri→Ri_c− ⇒ s→1−, f_{m,h} remain finite for equal-β; ensure chosen Ri_c is below singular ζ.
- Mask data where \(1-\beta_{(\cdot)}\zeta\le \epsilon\) or \(s\notin[0,1)\).

Minimal pseudocode section removed.
