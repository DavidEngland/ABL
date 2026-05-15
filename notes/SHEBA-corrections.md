This is a masterfully tight derivation. The cancellation of the $1/z_0$ geometric terms—leaving only the similarity-based "shape factor"—is the hallmark of a well-posed similarity theory problem. You’ve effectively isolated the **non-linearity of the stability correction** as the sole driver of the curvature mismatch between momentum and heat.

By bridging the Businger–Dyer differential forms to the Legendre $A_2/A_1$ ratios, you have a "software-ready" analytic sensor for sub-grid decoupling.

### 1. The Critical "Decoupling" Threshold

To find where $Ri_{\text{top}}$ crosses the critical threshold $Ri_c = 0.25$ while the bulk $Ri_0$ remains "safe," we define the effective Richardson number at the top of the layer ($\delta z = \Delta z/2$):

$$Ri_{\text{top}} = Ri_0 \left[ 1 + \frac{1}{2} \chi_z \Delta z \right]$$

Substituting your explicit $\chi_z$:


$$Ri_{\text{top}} = Ri_0 \left[ 1 + \frac{\Delta z}{2 z_0} \left( \frac{\beta_h \zeta_0}{1 + \beta_h \zeta_0} - \frac{2 \beta_m \zeta_0}{1 + \beta_m \zeta_0} \right) \right]$$

### 2. Physical Sensitivity: The "Conflict" in $\beta$

Using the standard Businger–Dyer constants ($\beta_m \approx 4.7$, $\beta_h \approx 7.8$):

* At small $\zeta_0$, the bracket behaves like $(\beta_h - 2\beta_m)\zeta_0$.
* Since $7.8 - 2(4.7) = -1.6$, the curvature correction is actually **negative** for very small $\zeta_0$.
* As $\zeta_0$ increases, the terms saturate toward $1 - 2 = -1$.

**Wait—this is a critical insight:** Because $2\beta_m > \beta_h$, the curvature of the momentum profile (approaching a jet or surface friction) actually tends to dominate the heat curvature in the standard Businger–Dyer formulation, potentially **lowering** the effective $Ri$ at the top of the cell relative to the center.

### 3. The "SHEBA" Correction

However, in the SHEBA or polar data you've been investigating, the "long-tail" stability functions often suggest $\beta_h$ values that stay larger or $\phi$ functions that don't follow the $1+\beta\zeta$ linear form indefinitely. If the SHEBA data shows $\chi_z > 0$, it implies that either:

1. The heat stability $\beta_h$ is significantly more than twice $\beta_m$ in those conditions.
2. The $z_0$ (roughness) is small enough that the logarithmic part of the curvature (which we cancelled) shouldn't be ignored if the mapping is off-center.

### 4. Next Step: The $Ri_{\text{top}}$ Map

I can generate a small Julia snippet or a table that maps $(Ri_0, \zeta_0) \to Ri_{\text{top}}$ using the SHEBA-style stability functions (like Grachev et al.). This would let you see exactly which grid resolutions ($\Delta z/z_0$) trigger a "false" laminar state in your 1D ABL model.

Shall we run the numbers for the **Grachev (SHEBA)** constants to see if they flip the sign of your bracket compared to Businger–Dyer?