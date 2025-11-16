...existing code...

## 7. Dynamic $Ri_c^*$ for Hybrid MOST/Ri Closures

**Concept:** Use $Ri_c^*$ as regime switch between MOST (below) and Ri closures (above).

**Formulation:**
$$
Ri_c^* = Ri_{c,0} \left[1 + \alpha_\Gamma \left(\frac{\Gamma}{\Gamma_{\text{ref}}} - 1\right) + \alpha_S \left(\frac{S}{S_{\text{ref}}} - 1\right) + \alpha_T \frac{\text{TKE}_{\text{prev}}}{\text{TKE}_{\text{ref}}}\right]
$$

**Regime Classification:**
- $Ri < 0.7 Ri_c^*$: MOST with φ functions + curvature correction
- $0.7 Ri_c^* \le Ri \le 1.3 Ri_c^*$: Blend zone
- $Ri > 1.3 Ri_c^*$: Ri-based $f_m, f_h$ (no L iteration)

**See:** `dynamic_Ric_strategy.md` for full implementation details.

...existing code...