# Generalized Richardson Curvature via Log–Derivative (Bell Polynomial) Lemma

## 1. Purpose
Provide a compact, model‑agnostic formula for the curvature of the gradient Richardson number
\[
Ri_g(\zeta)=\zeta\,F(\zeta),\qquad F(\zeta)=\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2},
\]
valid for any differentiable stability functions \(\phi_m,\phi_h\) (MOST power law, quadratic truncation, regularized, capped, blended). This lets an atmospheric scientist (a) audit parameterizations, (b) quantify curvature sensitivity, (c) compare analytic vs. numerically differentiated model fields.

## 2. Log–Derivative Setup
Define the primary logarithmic derivative
\[
G(\zeta)\equiv \frac{F'(\zeta)}{F(\zeta)}=\frac{d}{d\zeta}\ln F(\zeta)
=\frac{1}{\phi_h}\frac{d\phi_h}{d\zeta}-2\frac{1}{\phi_m}\frac{d\phi_m}{d\zeta}.
\]
For the classic power‑law MOST forms \(\phi_{(\cdot)}=(1-\beta_{(\cdot)}\zeta)^{-\alpha_{(\cdot)}}\),
\[
G(\zeta)=\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta}=V_{\log}(\zeta),\qquad
G'(\zeta)=\frac{\alpha_h\beta_h^2}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^2}{(1-\beta_m\zeta)^2}=W_{\log}(\zeta).
\]
(Here \(W_{\log}\) equals \(G'\); earlier power‑law notation sometimes uses \(-W_{\log}\)—we adopt the sign above for consistency.)

## 2A. Mapping to MOST Notation
Let
\[
v_m=\frac{\phi_m'}{\phi_m},\quad v_h=\frac{\phi_h'}{\phi_h},\quad
V_{\log}=v_h-2v_m,\quad W_{\log}=V_{\log}'.
\]
Then the generic symbols of Section 2:
\[
G \equiv V_{\log},\qquad G' \equiv W_{\log},\qquad G'' \equiv \frac{dW_{\log}}{d\zeta}.
\]
Curvature formula equivalence:
\[
R''_g=F\big[2G+\zeta(G^{2}+G')\big]=F\big[2V_{\log}+\zeta(V_{\log}^{2}-W_{\log})\big].
\]

## 3. Bell Polynomial Lemma
Let \(x_1=G, x_2=G', x_3=G'',\dots\). The complete exponential Bell polynomials \(B_n\) give
\[
F^{(n)}(\zeta)=F(\zeta)\,B_n\big(G,G',\dots,G^{(n-1)}\big).
\]
First three:
\[
B_1=x_1,\quad B_2=x_1^2+x_2,\quad B_3=x_1^3+3x_1x_2+x_3.
\]

## 4. Curvature (Second Derivative of \(Ri_g\))
Starting with
\[
Ri_g=\zeta F,\quad R'_g=F+\zeta F',\quad R''_g=2F'+\zeta F''.
\]
Substitute \(F'=F G\), \(F''=F(G^2+G')\):
\[
\boxed{R''_g(\zeta)=F(\zeta)\Big[2G(\zeta)+\zeta\big(G(\zeta)^2+G'(\zeta)\big)\Big]}.
\]

## 5. Third Derivative (Change of Curvature)
Using \(F'''=F(G^3+3GG'+G'')\),
\[
R'''_g=3F''+\zeta F'''=F\Big[3\big(G^2+G'\big)+\zeta\big(G^3+3GG'+G''\big)\Big].
\]
Neutral (ζ→0) limits:
\[
R''_g(0)=2G(0),\qquad R'''_g(0)=3\big(G(0)^2+G'(0)\big).
\]

## 6. Atmospheric Interpretation
- Sign of \(G(0)\): initial concavity of \(Ri_g(\zeta)\). If \(G(0)>0\) ⇒ curvature positive (concave‑up); if \(G(0)<0\) ⇒ concave‑down.
- Magnitude of \(G'(0)\) modulates how quickly curvature strengthens or weakens with ζ.
- Inflection in \(Ri_g\) occurs when \(R''_g(\zeta)=0\); approximate near neutral by linearizing:
\[
R''_g(\zeta)\approx 2G(0)+\zeta\,(G(0)^2+G'(0))\Rightarrow 
\zeta_{\text{inf}}\approx -\frac{2G(0)}{G(0)^2+G'(0)}\quad (\text{valid if } | \zeta_{\text{inf}}|\ll \zeta_\text{domain}).
\]

## 6A. Climate / ABL Modeling Emphasis
- Unified \(V_{\log},W_{\log}\) let mixed profile families share identical curvature evaluation logic.
- Diagnostic: large \(|W_{\log}/V_{\log}^{2}|\) flags layers where stability correction curvature may dominate shear—candidate for vertical refinement.

## 7. Mapping to Height (Variable \(L(z)\))
With \(\zeta=z/L(z)\):
\[
\frac{d\zeta}{dz}=\frac{L-zL'}{L^2},\quad
\frac{d^2\zeta}{dz^2}=-\frac{2L'}{L^2}-\frac{zL''}{L^2}+\frac{2z(L')^2}{L^3}.
\]
Height‑space curvature:
\[
\frac{\partial^2 Ri_g}{\partial z^2}=\left(\frac{d\zeta}{dz}\right)^2 R''_g(\zeta)+\frac{d^2\zeta}{dz^2} R'_g(\zeta),
\quad R'_g=F(1+\zeta G).
\]
If vertical variation of \(L\) is weak, the simpler \(\partial_z^2 Ri_g \approx R''_g/L^2\) suffices; otherwise include full chain rule.

## 8. Practical Algorithm (Column Diagnostics)
1. Choose φ family (power law, quadratic, regularized, blended). Ensure domain validity (no pole or capped).
2. Evaluate φ_m, φ_h and their first / second ζ-derivatives (analytic preferred; else central difference with small step h).
3. Form \(G, G'\) (and \(G''\) if third derivative needed).
4. Compute curvature \(R''_g\) via boxed formula. Map to height if L(z) varies.
5. Smooth noisy derivatives (Savitzky–Golay or spline) before forming \(G'\) to avoid amplification.
6. Compare analytic curvature with numerically differentiated \(Ri_g(z)=N^2/S^2\) (bias, RMSE, sign agreement).
7. Store neutral metrics: \(G(0), G'(0)\) for parameter tuning.

## 9. Derivatives for Common φ Forms
Power law: \(\phi=(1-\beta\zeta)^{-\alpha}\)
\[
\phi'=\frac{\alpha\beta}{1-\beta\zeta}\phi,\quad
\phi''=\frac{\alpha\beta^2(1+\alpha-\beta\zeta)}{(1-\beta\zeta)^2}\phi.
\]

Quadratic (stable truncation): \(\phi=1+a\zeta+b\zeta^2\)
\[
\phi'=a+2b\zeta,\quad \phi''=2b.
\]

Regularized power law (example): \(\phi=(1+g)^\alpha,\ g=\frac{\beta\zeta}{1+\delta\beta\zeta}\)
Compute \(g', g''\) analytically, then
\[
\phi'=\alpha(1+g)^{\alpha-1}g',\quad
\phi''=\alpha(1+g)^{\alpha-2}\big[(\alpha-1)g'^2+(1+g)g''\big].
\]

Capped linear (stable tail): \(\phi=\min(\phi_{\text{base}}, 1+c\zeta)\).
Derivative piecewise; curvature only from active branch. Junction: treat \(G, G'\) using one‑sided derivatives or smooth blend.

## 10. Numerical Notes
- Step size h for finite differences: choose \(h\sim 10^{-5}\)–\(10^{-4}\) times ζ-range to balance truncation vs round‑off.
- Reject levels where \(|1-\beta\zeta|<\epsilon\) (pole proximity); use surrogate (quadratic or regularized) there.
- In variable L columns, large |L'| enhances second term of height curvature; apply mild vertical smoothing to L first.
- Provide uncertainty: propagate φ derivative errors through linearization:
\[
\delta R''_g \approx |F|\Big[2\,\delta G + \zeta(2|G|\,\delta G + \delta G')\Big].
\]

## 11. Applications for Atmospheric Scientist
- Parameter Tuning: Adjust (α,β) until neutral curvature \(2G(0)\) matches observed early‑ζ Richardson growth.
- Regime Classification: Use sign and magnitude of \(R''_g\) to tag transition toward intermittent turbulence.
- Closure Validation: Compare curvature shapes among competing φ families (power law vs quadratic vs regularized) under identical ζ spans.
- Data Assimilation: Include neutral curvature and inflection height as additional quality control measures on retrieved profiles.
- Stability Monitoring: Rapid analytic curvature prevents misclassification from noisy N²/S² finite differences at coarse vertical resolution.

## 12. Summary Formulas (Ready for Drop‑In)
\[
\boxed{
\begin{aligned}
Ri_g &= \zeta F,\quad F=\frac{\phi_h}{\phi_m^2},\\
G &= \frac{1}{\phi_h}\phi_h' - 2\frac{1}{\phi_m}\phi_m',\\
R''_g &= F\Big[2G+\zeta(G^2+G')\Big],\\
R'''_g &= F\Big[3(G^2+G')+\zeta(G^3+3GG'+G'')\Big].
\end{aligned}}
\]
Neutral:
\[
R''_g(0)=2G(0),\qquad R'''_g(0)=3\big(G(0)^2+G'(0)\big).
\]

## 13. Minimal Implementation Snippet
```python
def rig_curvature(zeta, phi_m, phi_h, dphi_m, dphi_h, d2phi_m=None, d2phi_h=None):
    F = phi_h/(phi_m**2)
    G  = (dphi_h/phi_h) - 2*(dphi_m/phi_m)
    # If analytic second derivatives available:
    if d2phi_m is not None and d2phi_h is not None:
        Gp = (d2phi_h/phi_h - (dphi_h/phi_h)**2) - 2*(d2phi_m/phi_m - (dphi_m/phi_m)**2)
    else:
        Gp = None  # supply via numerical differencing externally
    curv = F*(2*G + zeta*(G**2 + (Gp if Gp is not None else 0.0)))
    return curv, F, G, Gp
```

## 14. Checklist (Model Post‑Processing)
- [ ] Compute φ_m, φ_h and derivatives (analytic preferred).
- [ ] Evaluate \(F, G, G'\).
- [ ] Curvature \(R''_g\); map to height (include variable L if needed).
- [ ] Compare with centered-difference curvature from Ri_g(z).
- [ ] Record neutral metrics and inflection estimate.
- [ ] Flag levels with |curvature| spikes due to derivative noise.

## 15. Key Pitfalls
- Mixing sign convention for \(W_{\log}\) vs \(G'\) (resolved here by defining \(G'=W_{\log}\)).
- Using raw finite differences near ζ poles without surrogate (inflates curvature).
- Neglecting variable L contributions in strongly stratified shallow layers.
- Over‑smoothing derivatives (can erase genuine inflection signatures).

## 16. Reference Mapping
If you previously used power‑law formulas:
Replace every \(V_{\log}\) with \(G\), and every \((V_{\log}^2-W_{\log})\) with \((G^2+G')\) under the sign convention above.

## 17. Height Mapping Simplifications for \(L(z)\)

Core mapping:
\[
\partial_z^2 Ri_g=(\zeta'_z)^2 R''_g+\zeta''_z R'_g,\quad \zeta'_z=\frac{L-zL'}{L^2},\ \zeta''_z=-\frac{2L'}{L^{2}}-\frac{zL''}{L^{2}}+\frac{2z(L')^{2}}{L^{3}}.
\]

Cases:
- Constant \(L=L_0\): \(\zeta'_z=1/L_0,\ \zeta''_z=0 \Rightarrow \partial_z^2Ri_g=(1/L_0^2)R''_g.\)
- Affine \(L=L_0+\lambda z\): use analytic forms, add correction term.
- Power-law / exponential: same forms as curvature docs (reuse coefficients).

Neglect criterion:
\[
E_{\text{omit}}=\left|\frac{\zeta''_z R'_g}{(\zeta'_z)^2 R''_g}\right|<\epsilon\ \Rightarrow\ \text{drop } \zeta''_z\text{ term}.
\]

Quick selector:
```python
def map_curv(z, L, dL, d2L, R1, R2, eps=0.05):
    dzeta = (L - z*dL)/(L*L)
    d2zeta = -2*dL/(L*L) - z*d2L/(L*L) + 2*z*dL*dL/(L*L*L)
    omit_err = abs(d2zeta*R1)/(dzeta*dzeta*R2) if R2!=0 else 0
    if omit_err < eps:
        return dzeta*dzeta*R2
    return dzeta*dzeta*R2 + d2zeta*R1
```

Use \(R'_g=F(1+\zeta G)\), \(R''_g=F[2G+\zeta(G^2+G')]\) for substitution.

A single generalized curvature expression plus derivative recipes allows rigorous, transparent auditing of any MOST‑class or extended stable boundary layer profile within your modeling or observational analysis workflow.