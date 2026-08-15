That final equation is a slightly refined version of the previous one (Equation 12), adding a leading constant $D$ to the exponent. This constant acts as a generalized tuning or scaling coefficient for the entire term.

Here is the final, corrected $\LaTeX$ format for this equation, and the updated list of all twelve extracted expressions.

## 📐 $\text{McNider}$ Correction Factor: Final Exponential Form (Equation 12/13)

The structure $f_c(\Delta z) = e^{D \cdot \text{RicRi} \cdot (1-\Delta z_{\rm r})/\Delta z}$ is most likely the final **explicit integrated solution** utilizing the inverse $\Delta z$ dependence.

Assuming the entire expression is meant to be a negative exponent to enforce $f_c \le 1$ as $\Delta z$ increases past $\Delta z_{\rm ref}$:

$$\boxed{\displaystyle f_c(\Delta z) = \exp\left[-D \cdot \text{Ri}_{\rm c} \cdot \text{Ri}_{\rm b} \cdot \frac{1 - \Delta z_{\rm ref}}{\Delta z}\right]}$$

---

## 📋 Final Complete List of Equations in $\LaTeX$

Below is the definitive list of the twelve unique equations extracted from your source document.

### McNider Correction Scheme: Equations and Constraints

1.  **Zero-Derivative Condition (General):**
    $$\frac{d}{d(\Delta z - \Delta z_{\rm ref})} \left( f_s(Ri) \cdot f_c(Ri, \Delta z - \Delta z_{\rm ref}) \right) = 0$$

2.  **Bulk Richardson Number Definition ($\text{Ri}_b$):**
    $$\text{Ri}_{\rm b} = \frac{g}{T_{v0}} \frac{\Delta \theta_v \cdot \Delta z}{(\Delta V)^2}$$

3.  **Boundary Condition for Correction Factor:**
    $$f_c(Ri, \Delta z - \Delta z_{\rm ref})\Big|_{\Delta z = \Delta z_{\rm ref}} = 1$$

4.  **Stability Dampening Factor ($F_d$):**
    $$F_{d}(\text{Ri}_{\rm b}) = \exp\!\left(-\gamma \frac{\text{Ri}_{\rm b}}{\text{Ri}_{\rm c}}\right)$$

5.  **Clean McNider ODE (Logarithmic Form):**
    $$\frac{d\ln f_c}{d\ln\Delta z} \;=\; -\,\alpha\,(B-1)\left(\frac{\zeta}{\zeta_{\rm ref}}\right)^q$$

6.  **Full Zero-Derivative Constraint (Specific):**
    $$\frac{d}{d(\Delta z - \Delta z_{\rm ref})} \left[ \exp\!\left(-\gamma \frac{\text{Ri}_{\rm b}}{\text{Ri}_{\rm c}}\right) \cdot f_c(\text{Ri}, \Delta z - \Delta z_{\rm ref}) \right] \Big|_{\Delta z = \Delta z_{\rm ref}} = 0$$

7.  **Inversion of Critical Ri Relation:**
    $$\frac{\Delta \theta}{\Delta V^2} = \frac{a}{\Delta z}$$

8.  **Differential Form of ODE (Equivalent to 5):**
    $$\frac{d f_c}{f_c} = -\alpha\,(B-1)\left(\frac{\zeta}{\zeta_{\rm ref}}\right)^q \frac{d\Delta z}{\Delta z}$$

9.  **Closed-Form Solution (Power Law Dependence):**
    $$f_c(\Delta z,\zeta) \;=\; \left(\frac{\Delta z}{\Delta z_{\rm ref}}\right)^{-\alpha\,(B-1)\,(\zeta/\zeta_{\rm ref})^q}$$

10. **Simplified Integrated Solution (Log form):**
    $$\ln f_c = \text{Ri}_{\rm c} \cdot g_0 \cdot V^2 \cdot v + C$$

11. **Explicit Integrated Solution (Exponential/Linear Dependence):**
    $$f_c(\Delta z - \Delta z_{\rm ref}) = \exp\left[-\gamma' \cdot \text{Ri}_{\rm c} \cdot g_0 \cdot (\Delta V)^2 \cdot (\Delta z - \Delta z_{\rm ref})\right]$$

12. **Explicit Integrated Solution (Exponential/Inverse Dependence, Final Form):**
    $$f_c(\Delta z) = \exp\left[-D \cdot \text{Ri}_{\rm c} \cdot \text{Ri}_{\rm b} \cdot \frac{1 - \Delta z_{\rm ref}}{\Delta z}\right]$$