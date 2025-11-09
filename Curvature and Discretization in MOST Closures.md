# Technical Note: Curvature and Vertical Discretization in MOST-Based Closures

## 0. Abstract
This note refines the interpretation of the gradient Richardson number curvature
\[
Ri_g(\zeta)=\zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2},\qquad \zeta=z/L,
\]
and its implications for stability-function formulation and vertical discretization in Monin–Obukhov Similarity Theory (MOST). We formalize near-neutral curvature control via the coefficient \(\Delta\), provide quantitative criteria for when a constant-\(L\) approximation is acceptable, and justify the geometric-mean height for layer evaluation. Implementation diagnostics and a ready figure caption are included.

## 1. Core Curvature Expression
Let
\[
F(\zeta)=\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2},\qquad
V_{\log}=(\phi_h'/\phi_h)-2(\phi_m'/\phi_m),\qquad
W_{\log}=dV_{\log}/d\zeta.
\]
Then
\[
\frac{d^2Ri_g}{d\zeta^2}=F(\zeta)\Big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\Big].
\]

Near neutrality (\(\zeta\to0\)):
\[
\Delta=V_{\log}(0)=\alpha_h\beta_h-2\alpha_m\beta_m,\quad
\left.\frac{d^2Ri_g}{d\zeta^2}\right|_{0}=2\Delta.
\]

Interpretation:
- \(\Delta>0\): initial concave-up (heat corrections dominate).
- \(\Delta<0\): initial concave-down (momentum corrections dominate).
- \(|\Delta|\) sets the strength of first departure from linear \(Ri_g(\zeta)\).

## 1A. Unified Log Terms
Component logs:
\[
v_m=\frac{\phi_m'}{\phi_m},\ v_h=\frac{\phi_h'}{\phi_h},\ V_{\log}=v_h-2v_m,\ W_{\log}=V_{\log}'.
\]
Curvature:
\[
\partial_\zeta^2 Ri_g = F[2V_{\log}+\zeta(V_{\log}^2-W_{\log})],\quad F=\phi_h/\phi_m^2.
\]
This replaces earlier generic G,G′ (G=V_log, G′=W_log).

## 2. Elevated Regime Behavior
Away from neutral (\(\zeta\gtrsim 0.1\) or approaching stable growth), the \(\zeta(V_{\log}^2-W_{\log})\) term drives rapid curvature escalation or decay. Parameter sets that keep \(V_{\log}^2\) and \(W_{\log}\) in partial balance delay curvature amplification, improving numerical stability in coarse vertical grids.

## 3. Discretization and Representative Height Selection
Given two model levels \(z_1<z_2\) with \(z_2/z_1\gg1\) possible under coarse grids:

1. Arithmetic mean \(z_a=(z_1+z_2)/2\) biases evaluation for MOST because profiles scale approximately with \(\ln z\).
2. Geometric mean \(z_g=\sqrt{z_1 z_2}\) satisfies
\[
\ln z_g = \frac{1}{z_2-z_1}\int_{z_1}^{z_2} \ln z\,dz
\]
minimizing the mean-square error of replacing the layer-average of any function linear in \(\ln z\) by a point evaluation.

Approximate curvature preservation error when using \(z_a\) instead of \(z_g\):
\[
\epsilon_{\text{geom}}\approx \frac{1}{2}\Big|\frac{d^2Ri_g}{dz^2}\Big|_{z_g}\left(z_a-z_g\right)^2
\]
while for \(z_g\) the first-order logarithmic bias term cancels.

Recommendation: use \(z_g\) for evaluating \(\phi_{m,h}\), \(Ri_g\), and curvature when \(\Delta z / z_1>0.2\).

## 4. Height Mapping with Variable \(L(z)\)
General chain rule:
\[
\frac{\partial^2Ri_g}{\partial z^2} = \left(\frac{d\zeta}{dz}\right)^2\frac{d^2Ri_g}{d\zeta^2} + \frac{d^2\zeta}{dz^2}\frac{dRi_g}{d\zeta},\quad \zeta=\frac{z}{L(z)}.
\]
Variation metrics:
\[
\varepsilon_1=\frac{z|L'|}{L},\qquad
\chi=\left|\frac{(d^2\zeta/dz^2)(dRi_g/d\zeta)}{(d\zeta/dz)^2(d^2Ri_g/d\zeta^2)}\right|.
\]
Use constant-\(L\) shortcut (\(\partial_z^2Ri_g\approx (1/L^2)\partial_\zeta^2Ri_g\)) if \(\varepsilon_1<0.05\) and \(\chi<0.05\); otherwise include full mapping.

## 5. Discretization Error Diagnostic
Layer reconstruction for a vertical interval \([z_i,z_{i+1}]\):
\[
E_i = \Big|\frac{Ri_g(z_{i+1}) - Ri_g(z_i)}{z_{i+1}-z_i} - \left.\frac{dRi_g}{dz}\right|_{z_g}\Big|
\]
Estimate using \(\frac{dRi_g}{dz} \approx (d\zeta/dz)dRi_g/d\zeta\). Flag layers with \(E_i / |dRi_g/dz|_{z_g} > \eta\) (e.g. \(\eta=0.2\)) for adaptive refinement or alternative φ-form (e.g. quadratic surrogate).

## 6. Practical Implementation Steps
1. Compute φ_m, φ_h at each level; evaluate \(F, V_{\log}, W_{\log}\).
2. Determine \(\Delta\) and neutral curvature; store for run metadata.
3. For each layer, compute \(z_g\); evaluate curvature at \(z_g\).
4. Assess \(L(z)\) variability; choose mapping logic (Section 4).
5. Compute \(E_i\) (Section 5); apply smoothing only if curvature spikes are noise (avoid suppressing true inflection).
6. Diagnostics output: \(\Delta, \partial_\zeta^2Ri_g|_{0},\) max \(|\partial_z^2Ri_g|\), first \(\zeta_{\text{inf}}\) if real and inside domain, fraction of layers with \(\chi>0.05\).

## 7. Blended / Piecewise Stability Functions
To enforce curvature continuity from near-neutral to elevated stable layers:
- Match neutral curvature: ensure surrogate (quadratic or regularized) has same \(2\Delta\).
- Continuity target: \(\left|\Delta_{\text{curv}}\right| = \left|(\partial_\zeta^2Ri_g^{\text{base}} - \partial_\zeta^2Ri_g^{\text{surrogate}})/\partial_\zeta^2Ri_g^{\text{base}}\right| < 0.05\) at blend height \(\zeta_b\).
- Use sigmoid blend in ζ (e.g. \(\sigma(\zeta)=1/(1+e^{-s(\zeta-\zeta_b)})\)) if a sharp switch creates curvature artifacts.

## 8. Recommended Diagnostics Summary
| Metric | Purpose | Threshold |
|--------|---------|-----------|
| \(2\Delta\) | Neutral curvature sign/magnitude | Report |
| \(\zeta_{\text{inf}}\) | Inflection location | Only if real & < domain limit |
| \(\max\lvert \partial_z^2 Ri_g \rvert\) | Vertical stability stress | Watch for spikes |
| \(\chi\) fraction | Variable-L impact | &lt;10% desirable |
| \(E_i\) distribution | Layer reconstruction error | &gt;0.2 flagged |
| Blend mismatch | Curvature continuity | &lt;5% |

## 9. Figure Caption (Curvature vs. \(\zeta\))
Caption candidate:
“Gradient Richardson number curvature \(\partial_\zeta^2Ri_g\) as a function of nondimensional height \(\zeta=z/L\). Neutral curvature \(2\Delta\) sets the initial concavity; rapid growth governed by \(\zeta(V_{\log}^2-W_{\log})\). Shaded band shows layers where \(\chi>0.05\), indicating significant \(L(z)\) variation. Dashed curve: quadratic stable surrogate matched on neutral curvature.”

## 10. Common Pitfalls
- Using arithmetic mean height with large \(\Delta z\): biases curvature and suppresses inflection detection.
- Over-smoothing derivatives: can erase legitimate curvature sign changes.
- Ignoring \(L(z)\) variability when \(\chi\) is large: misrepresents physical curvature scaling.
- Applying power-law φ beyond pole proximity without guard: can inflate curvature and collapse flux iterates.

## 11. Extension Notes
- Unstable branch: sign changes in \(\beta\zeta\) require separate evaluation; geometric mean justification still holds for log-scaling.
- Planetary adaptation: only \(L\) scaling changes; curvature expression identical in ζ.
- Ri-based closures: use ζ(Ri) inversion series + single Newton step to avoid iterative ζ loops per time step.

## 12. Summary
Curvature structure is controlled near neutrality by a single coefficient \(\Delta\); representative-height choice (geometric mean) reduces discretization bias for log-structured profiles. Quantitative metrics (\(\varepsilon_1,\chi,E_i\)) guide when constant-\(L\) assumptions and layer compression are valid. Implementing curvature-aware blending stabilizes SBL behavior and improves diagnostic fidelity.

## 13. (Optional) Quick Code Stub
```python
def layer_curvature(z1,z2,L1,L2,phi_m,phi_h):
    import math
    z_g = math.sqrt(z1*z2)
    # user supplies phi_m/h callables at ζ=z/L(z); approximate L(z_g) by geometric mean
    L_g = math.sqrt(L1*L2)
    zeta_g = z_g/L_g
    # central differences for derivatives (small h)
    h = 1e-6
    def F(zeta): return phi_h(zeta)/phi_m(zeta)**2
    def logF(zeta): return math.log(F(zeta))
    V = (logF(zeta_g+h)-logF(zeta_g-h))/(2*h)
    W = (logF(zeta_g+h)-2*logF(zeta_g)+logF(zeta_g-h))/(h*h)
    curv_zeta = F(zeta_g)*(2*V + zeta_g*(V*V - W))
    return curv_zeta
```

## 14. Critique of Original Draft (For Record)
- Mixed plain-text and LaTeX markers (e.g., `\partial^2 Ri_g/\partial \zeta^2`) without consistent math environments.
- Missing explicit definitions of \(V_{\log}, W_{\log}, \Delta\).
- No quantitative criteria for when constant-\(L\) approximation is valid.
- Geometric mean justification stated qualitatively only; now formalized via integral of \(\ln z\).
- Lacked discretization error metric; now \(E_i\) defined.
- No guidance for figure captions or diagnostics table.

All enhancements preserve original intent while adding formal rigor and actionable criteria.
