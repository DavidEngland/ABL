## Connection to McNider ODE (Optional Advanced Note)

McNider proposed a differential constraint for f_c:
$$
\frac{d}{d(\Delta z)}(f_s(Ri) \cdot f_s(Ri,\Delta z)) = 0,
$$
where f_s is a reference stability function. This enforces that the **product** (fc · layer-thickness · stability) remains invariant under grid refinement—a mathematical expression of the requirement that the **integrated turbulent transport** should be grid-independent.

Taking logarithmic derivative:
$$
\frac{d \ln fc}{d \ln \Delta z} = -\alpha (B-1) \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q,
$$
which integrates to:
$$
fc(\Delta z, \zeta) = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}.
$$

For small exponents, this approximates:
$$
fc \approx \exp\left[-\alpha (B-1) \ln\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right) \left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q\right],
$$
which is equivalent to the exponential form with p=1 and logarithmic Δz scaling. The power-law form:
$$
\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p
$$
is a **practical generalization** allowing non-logarithmic Δz dependence (p≠1) for empirical tuning.

**Takeaway:** The exponential fc is not ad hoc—it is a **closed-form solution** to a physically motivated ODE expressing grid-invariance of integrated turbulent flux.

---
