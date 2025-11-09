# Curvature of the Gradient Richardson Number \(Ri_g\)

> Note on z vs. ζ: We ultimately need the curvature with respect to height z, ∂²Ri_g/∂z². We compute analytically in the nondimensional coordinate ζ=z/L, then convert at the end via the chain rule: ∂²/∂z² = (1/L²) ∂²/∂ζ² (treating L locally constant in the surface layer).

## 1. Purpose
Provide a clean, citation‑ready derivation and interpretation of \(\partial^{2}Ri_g/\partial \zeta^{2}\) for power‑law MOST stability functions to support discussion (McNider, Biazar).

## 2. Definitions
Dimensionless height: \(\zeta = z/L\).

Stability (momentum, heat):
\[
\phi_m(\zeta)=(1-\beta_m \zeta)^{-\alpha_m},\qquad
\phi_h(\zeta)=(1-\beta_h \zeta)^{-\alpha_h},
\quad (1-\beta_{(\cdot)}\zeta)>0.
\]

Gradient Richardson number (MOST form):
\[
Ri_g(\zeta)=\zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=\zeta\,F(\zeta),
\qquad
F(\zeta)= (1-\beta_h \zeta)^{-\alpha_h}(1-\beta_m \zeta)^{2\alpha_m}.
\]

Logarithmic derivative components (unified general form):
\[
v_m=\frac{1}{\phi_m}\frac{d\phi_m}{d\zeta},\qquad
v_h=\frac{1}{\phi_h}\frac{d\phi_h}{d\zeta},\qquad
V_{\log}=v_h-2v_m,\qquad
W_{\log}=\frac{dV_{\log}}{d\zeta}= \frac{dv_h}{d\zeta}-2\frac{dv_m}{d\zeta}.
\]
For the power-law case:
\[
v_m=\frac{\alpha_m\beta_m}{1-\beta_m\zeta},\quad
v_h=\frac{\alpha_h\beta_h}{1-\beta_h\zeta},\quad
V_{\log}=\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta},
\]
\[
W_{\log}=
\frac{\alpha_h\beta_h^{2}}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^{2}}{(1-\beta_m\zeta)^2}.
\]
Note: Earlier “G” and “G′” in generic analysis correspond to \(V_{\log}, W_{\log}\).

Gradient Richardson number (MOST form):
\[
Ri_g(\zeta)=\zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=\zeta\,F(\zeta),
\qquad
F(\zeta)= (1-\beta_h \zeta)^{-\alpha_h}(1-\beta_m \zeta)^{2\alpha_m}.
\]

Logarithmic derivative components:
\[
V_{\log}=\frac{1}{\phi_h}\frac{d\phi_h}{d\zeta}-\frac{2}{\phi_m}\frac{d\phi_m}{d\zeta}
=\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta},
\]
\[
W_{\log}=\frac{d}{d\zeta}V_{\log}
=\frac{\alpha_h\beta_h^{2}}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^{2}}{(1-\beta_m\zeta)^2}.
\]

Useful identity:
\[
\frac{1}{\phi}\frac{d\phi}{d\zeta}=\frac{\alpha\beta}{1-\beta\zeta}\quad\text{for}\quad \phi=(1-\beta\zeta)^{-\alpha}.
\]

## 3. Curvature Derivation
Starting from \(Ri_g=\zeta F\):
\[
\frac{dRi_g}{d\zeta}=F+\zeta F',\qquad
\frac{d^{2}Ri_g}{d\zeta^{2}}=2F'+\zeta F''.
\]

With
\[
\frac{F'}{F}=V_{\log},\qquad
\frac{F''}{F}=V_{\log}^2-W_{\log},
\]
we obtain the compact curvature form:
\[
\boxed{\frac{d^{2}Ri_g}{d\zeta^{2}}=F(\zeta)\big[2V_{\log}+\zeta(V_{\log}^{2}-W_{\log})\big]}.
\]

Expanded (self‑contained) form:
\[
\frac{d^{2}Ri_g}{d\zeta^{2}}=
(1-\beta_h \zeta)^{-\alpha_h}(1-\beta_m \zeta)^{2\alpha_m}
\left\{
\frac{2\alpha_h\beta_h}{1-\beta_h\zeta}
-\frac{4\alpha_m\beta_m}{1-\beta_m\zeta}
+\zeta\left[
\left(\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta}\right)^2
-\left(\frac{\alpha_h\beta_h^{2}}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^{2}}{(1-\beta_m\zeta)^2}\right)
\right]
\right\}.
\]

## 4. Neutral Limit (\(\zeta\to 0\))
\[
F(0)=1,\quad V_{\log}(0)=\alpha_h\beta_h-2\alpha_m\beta_m,\quad
W_{\log}(0)=\alpha_h\beta_h^{2}-2\alpha_m\beta_m^{2},
\]
\[
\left.\frac{d^{2}Ri_g}{d\zeta^{2}}\right|_{0}=2(\alpha_h\beta_h-2\alpha_m\beta_m).
\]

Interpretation:
- Positive: initial concave‑up \(Ri_g(\zeta)\) (heat “more corrected”).
- Negative: concave‑down (momentum “more corrected”).
- Near zero: quasi‑linear onset of \(Ri_g\).

## 5. Small-\(\zeta\) Series (optional)
Expand logs:
\[
\ln F(\zeta)= (\alpha_h\beta_h-2\alpha_m\beta_m)\zeta
+\tfrac{1}{2}(\alpha_h\beta_h^{2}-2\alpha_m\beta_m^{2})\zeta^{2}+O(\zeta^{3}).
\]
Then
\[
F(\zeta)=1+A_1\zeta+\tfrac{1}{2}(A_2+A_1^{2})\zeta^{2}+O(\zeta^{3}),
\quad A_1=\alpha_h\beta_h-2\alpha_m\beta_m,\ A_2=\alpha_h\beta_h^{2}-2\alpha_m\beta_m^{2}.
\]
So
\[
Ri_g(\zeta)=\zeta + A_1\zeta^{2} + O(\zeta^{3}),\qquad
\frac{d^{2}Ri_g}{d\zeta^{2}}=2A_1+O(\zeta).
\]

## 6. Singular / Domain Behavior
Poles at \(\zeta=1/\beta_h\) and \(\zeta=1/\beta_m\). As either is approached from below:
\[
V_{\log}\sim \frac{\alpha\beta}{1-\beta\zeta},\qquad
W_{\log}\sim \frac{\alpha\beta^{2}}{(1-\beta\zeta)^2},
\]
giving algebraic blow‑up; practical filtering typically restricts \(\zeta<0.7/\max(\beta_h,\beta_m)\).

## 7. Parameter Sensitivity (local)
\[
\frac{\partial}{\partial \alpha_h}\left(\frac{d^{2}Ri_g}{d\zeta^{2}}\right)=
\frac{d^{2}Ri_g}{d\zeta^{2}}\left[-\ln(1-\beta_h\zeta)+\frac{2}{2V_{\log}+\zeta(V_{\log}^2-W_{\log})}
\left(\frac{\beta_h}{1-\beta_h\zeta}+\frac{\zeta\beta_h}{1-\beta_h\zeta}V_{\log}-\frac{\zeta\beta_h^{2}}{(1-\beta_h\zeta)^2}\right)\right]
\]
(similar forms for \(\alpha_m,\beta_h,\beta_m\); retained symbolically to avoid clutter—can be expanded if needed).

Qualitatively: curvature most sensitive when denominators \(1-\beta_{(\cdot)}\zeta\) shrink (near singular approach) and when \(2\alpha_m\beta_m\) and \(\alpha_h\beta_h\) nearly cancel (magnifying relative uncertainty).

## 8. Recommended Reporting Set
For comparison across parameterizations:
- \(\left.d^{2}Ri_g/d\zeta^{2}\right|_{\zeta=0}\)
- Height of curvature sign change (if any).
- Fractional enhancement factor:
  \[
  \mathcal{C}(\zeta)=\frac{(d^{2}Ri_g/d\zeta^{2})}{2(\alpha_h\beta_h-2\alpha_m\beta_m)}.
  \]

## 9. Summary (for distribution)
A single, closed expression:
\[
\frac{d^{2}Ri_g}{d\zeta^{2}}=F(\zeta)\big[2V_{\log}+\zeta(V_{\log}^{2}-W_{\log})\big]
\]
encapsulates curvature; neutral limit depends only on the *linear* stability coefficients \(\alpha\beta\). Sign of \(\alpha_h\beta_h-2\alpha_m\beta_m\) determines initial concavity; proximity to singular factors drives rapid growth and parameter sensitivity. Series and sensitivity forms enable controlled inversion / fitting; suitable for manuscript discussion.

## 10. Action Items (next discussion)
- Confirm parameter ranges (site data) ⇒ evaluate neutral curvature spectrum.
- Decide acceptance band for \(\zeta\) before curvature saturation.
- Assess whether curvature sign change correlates with empirical \(Ri_c\).
- Prepare figure set: (i) curvature vs. \(\zeta\); (ii) normalized \(\mathcal{C}(\zeta)\); (iii) sensitivity heat map.

## 11. Physical Height Curvature (chain rule)
Using ζ=z/L and treating L locally constant over the layer of interest,
\[
\boxed{\frac{\partial^{2}Ri_g}{\partial z^{2}}=\frac{1}{L^{2}}\;\frac{\partial^{2}Ri_g}{\partial \zeta^{2}}.}
\]
If L varies with height and that variation is retained, additional terms involving dL/dz appear; in MOST applications we evaluate curvature at a level using the local L, so the 1/L² factor is appropriate.

## 11A. Variable-L(z) Mapping and Model Diagnostics
When \(L\) varies with height, let \(\zeta(z)=z/L(z)\). Then
\[
\frac{d\zeta}{dz}=\frac{1}{L}-\frac{zL'}{L^{2}}=\frac{L-zL'}{L^{2}},\qquad
\frac{d^{2}\zeta}{dz^{2}}=-\frac{2L'}{L^{2}}-\frac{zL''}{L^{2}}+\frac{2z\,L'^2}{L^{3}}.
\]
Chain rule:
\[
\frac{\partial Ri_g}{\partial z}=\frac{d\zeta}{dz}\,\frac{dRi_g}{d\zeta},\quad
\boxed{\frac{\partial^{2} Ri_g}{\partial z^{2}}=\Big(\frac{d\zeta}{dz}\Big)^{2}\frac{d^{2}Ri_g}{d\zeta^{2}}+\frac{d^{2}\zeta}{dz^{2}}\frac{dRi_g}{d\zeta}}.
\]
Here
\[
\frac{dRi_g}{d\zeta}=F+\zeta F' = F\big(1+\zeta V_{\log}\big),\qquad
\frac{d^{2}Ri_g}{d\zeta^{2}}=F\big[2V_{\log}+\zeta(V_{\log}^{2}-W_{\log})\big].
\]

Practical recipe (per model time/column)
- Inputs on model levels: height array z, θ (or θ_v), u,v, turbulent fluxes to compute L(z).
- Compute L(z) and finite-difference L'(z), L''(z) with light smoothing.
- Build ζ_i = z_i/L_i; evaluate theoretical dRi_g/dζ and d²Ri_g/dζ² from chosen φ-set at ζ_i.
- Map to height curvature using the boxed formula above.
- Independently estimate Ri_g(z)=N^2/S^2 (or MOST form from φ) and compute numeric ∂²Ri_g/∂z² via centered differences (optional Savitzky–Golay smoothing).
- Compare theoretical vs numeric curvature: bias, RMSE, correlation, and sign-agreement.

Minimal code sketch
```python
import numpy as np

def diffs(y, z):
    dz = np.gradient(z)
    y1 = np.gradient(y, z)
    y2 = np.gradient(y1, z)
    return y1, y2

def map_curvature_z(z, L, zeta, F, Vlog, Wlog):
    # theory in ζ
    dRi_dzeta  = F*(1 + zeta*Vlog)
    d2Ri_dzeta2= F*(2*Vlog + zeta*(Vlog*Vlog - Wlog))
    # L-derivatives
    L1, L2 = diffs(L, z)
    dzeta_dz  = 1.0/L - (z*L1)/(L*L)
    d2zeta_dz2= -2*L1/(L*L) - (z*L2)/(L*L) + 2*z*(L1*L1)/(L*L*L)
    # mapped curvature
    curv_z = (dzeta_dz*dzeta_dz)*d2Ri_dzeta2 + d2zeta_dz2*dRi_dzeta
    return curv_z

# numeric curvature from modeled Ri_g(z)
def numeric_curvature(Ri_g, z, window=None):
    # optional smoothing step could be inserted here
    _, d2 = diffs(Ri_g, z)
    return d2
```

Notes
- Use consistent staggering for z, θ, u,v when forming N² and S²; compute L at the same levels used for ζ.
- Regularize L'(z), L''(z) (e.g., weak SG filter) to reduce amplification of noise in d²ζ/dz².
- Report both mapped curvature and direct numeric curvature with uncertainty bands from differencing step.

## 11B. Special Cases for \(L(z)\): Constant vs Structured Variability

Notation:
\[
\zeta=\frac{z}{L(z)},\quad \frac{d\zeta}{dz}=\frac{L-zL'}{L^2},\quad \frac{d^{2}\zeta}{dz^{2}}=-\frac{2L'}{L^{2}}-\frac{zL''}{L^{2}}+\frac{2z(L')^{2}}{L^{3}}.
\]
General curvature mapping:
\[
\frac{\partial^{2}Ri_g}{\partial z^{2}}=\left(\frac{d\zeta}{dz}\right)^{2}\frac{d^{2}Ri_g}{d\zeta^{2}}+\frac{d^{2}\zeta}{dz^{2}}\frac{dRi_g}{d\zeta}.
\]

Define nondimensional variation measures (evaluated at height z):
\[
\varepsilon_1=\frac{z|L'|}{L},\qquad \varepsilon_2=\frac{z^{2}|L''|}{L},\qquad \chi=\frac{|d^{2}\zeta/dz^{2}|}{(d\zeta/dz)^2}.
\]
Rule of thumb: constant-L approximation (use \(\partial_z^{2}Ri_g\simeq (1/L^2)\partial_{\zeta}^{2}Ri_g\)) is acceptable if \(\varepsilon_1\ll1\) and \(\chi\ll1\).

Special profiles:

1. Constant \(L=L_0\):
\[
L'=0,\ L''=0 \Rightarrow \frac{d\zeta}{dz}=\frac{1}{L_0},\ \frac{d^{2}\zeta}{dz^{2}}=0,\ 
\boxed{\frac{\partial^{2}Ri_g}{\partial z^{2}}=\frac{1}{L_0^2}\frac{d^{2}Ri_g}{d\zeta^{2}}.}
\]

2. Affine \(L(z)=L_0+\lambda z\) (weak linear trend):
\[
L'= \lambda,\ L''=0,\ 
\frac{d\zeta}{dz}=\frac{L_0}{(L_0+\lambda z)^2},\quad
\frac{d^{2}\zeta}{dz^{2}}=-\frac{2\lambda L_0}{(L_0+\lambda z)^3}.
\]
Curvature correction term:
\[
\Delta_{\text{var}}=\frac{d^{2}\zeta}{dz^{2}}\frac{dRi_g}{d\zeta}.
\]

3. Power-law \(L(z)=L_0 (z/z_0)^p\) (e.g. evolving stratification layer):
\[
L'=\frac{pL}{z},\quad L''=\frac{p(p-1)L}{z^{2}},
\]
\[
\frac{d\zeta}{dz}=\frac{1-p}{L},\quad
\frac{d^{2}\zeta}{dz^{2}}=\frac{p( p-1)}{zL}-\frac{2p(1-p)}{zL}=\frac{p(2p-1)}{zL}.
\]
If \(p\to 0\) recover constant \(L\); if \(p<1\), sign of \(d\zeta/dz\) remains positive.

4. Exponential \(L(z)=L_0 e^{\lambda z}\) (rapid growth/decay region):
\[
L'=\lambda L,\quad L''=\lambda^{2}L,\quad
\frac{d\zeta}{dz}=\frac{e^{-\lambda z}}{L_0}-\frac{z\lambda e^{-\lambda z}}{L_0}
=\frac{e^{-\lambda z}}{L_0}(1-\lambda z),
\]
\[
\frac{d^{2}\zeta}{dz^{2}}=\frac{e^{-\lambda z}}{L_0}(\lambda^{2} z -2\lambda).
\]

5. Weakly varying (perturbative): write \(L(z)=L_0(1+\delta \ell(z))\), \(|\delta|\ll1\). Expand:
\[
\frac{d\zeta}{dz}=\frac{1}{L_0}\big[1-\delta(\ell+ z\ell')+O(\delta^2)\big],\quad
\frac{d^{2}\zeta}{dz^{2}}=\frac{-\delta}{L_0}\big[2\ell'+z\ell''-2z\ell'^2\big]+O(\delta^2).
\]
Leading relative curvature correction:
\[
\frac{\Delta_{\text{var}}}{(1/L_0^2)d^{2}Ri_g/d\zeta^{2}}
\approx -\delta L_0 \frac{2\ell'+z\ell''-2z\ell'^2}{d^{2}Ri_g/d\zeta^{2}}\,\frac{dRi_g}{d\zeta}.
\]

Error from ignoring variability:
\[
E_{\text{omit}}\equiv \left|\frac{d^{2}\zeta/dz^{2}\; dRi_g/d\zeta}{(d\zeta/dz)^2 d^{2}Ri_g/d\zeta^{2}}\right|=\chi.
\]
Set acceptance if \(E_{\text{omit}}<\epsilon\) (e.g. \(\epsilon=0.05\)).

Minimal helper:
```python
def curvature_height(z, L, dL, d2L, ri_zeta, ri_zeta1, ri_zeta2):
    # ri_zeta = Ri_g(ζ); ri_zeta1 = dRi_g/dζ; ri_zeta2 = d2Ri_g/dζ2
    dzeta_dz = (L - z*dL)/(L*L)
    d2zeta_dz2 = -2*dL/(L*L) - z*d2L/(L*L) + 2*z*dL*dL/(L*L*L)
    return (dzeta_dz*dzeta_dz)*ri_zeta2 + d2zeta_dz2*ri_zeta1

def curvature_height_constant(L0, ri_zeta2):
    return ri_zeta2/(L0*L0)
```
Use constant form if \(\varepsilon_1,\chi\) below thresholds.

## 12. Dimensionless Control Parameters and Inflection Structure
Define
\[
G=\frac{\alpha_h\beta_h}{2\alpha_m\beta_m},\qquad \Delta = \alpha_h\beta_h-2\alpha_m\beta_m=2\alpha_m\beta_m(G-1).
\]
Neutral curvature sign depends only on \(\Delta\). A (positive) interior inflection satisfies
\[
2V_{\log}+\zeta(V_{\log}^2-W_{\log})=0.
\]
Near neutrality (\(|\zeta|\ll1\)) expand \(V_{\log}= \Delta + c_1\zeta + O(\zeta^2)\) with
\[
c_1=\alpha_h\beta_h^2-2\alpha_m\beta_m^2 = W_{\log}(0).
\]
Retaining leading orders gives a quadratic approximation for the first root (if \(\Delta c_1<0\)):
\[
\zeta_{\text{inf}} \approx -\frac{2\Delta}{\Delta^2-c_1}.
\]
Validity requires \(|\zeta_{\text{inf}}|\ll \min(1/\beta_h,1/\beta_m)\); otherwise curvature is monotone up to singular approach.

## 13. Higher Derivatives
With \(Ri_g=\zeta F\),
\[
\frac{d^3Ri_g}{d\zeta^3}=3F''+\zeta F'''.
\]
Using logarithmic derivatives:
\[
\frac{F'''}{F}=V_{\log}^3-3V_{\log}W_{\log}-\frac{dW_{\log}}{d\zeta},
\quad
\frac{dW_{\log}}{d\zeta}= \frac{2\alpha_h\beta_h^3}{(1-\beta_h\zeta)^3}-\frac{4\alpha_m\beta_m^3}{(1-\beta_m\zeta)^3}.
\]
Thus
\[
\frac{d^3Ri_g}{d\zeta^3}=F\Big[3(V_{\log}^2-W_{\log}) + \zeta\big(V_{\log}^3-3V_{\log}W_{\log}-\tfrac{dW_{\log}}{d\zeta}\big)\Big].
\]
Neutral limit:
\[
\left.\frac{d^3Ri_g}{d\zeta^3}\right|_{0}=3(\Delta^2-c_1).
\]
Hence the sign of \(\Delta^2-c_1\) controls initial *change* of curvature magnitude.

## 14. Inversion: \(\zeta(Ri_g)\) Near Neutrality
Set \(Ri_g=\zeta(1+A_1\zeta+A_2'\zeta^2+O(\zeta^3))\) with
\[
A_1=A_1^{(0)}=\Delta,\qquad A_2'=\tfrac12\big(A_2+A_1^2\big)=\tfrac12\left(c_1+\Delta^2\right).
\]
Series inversion gives
\[
\zeta = Ri_g - A_1 Ri_g^2 + (2A_1^2-A_2')Ri_g^3 + O(Ri_g^4),
\]
i.e.
\[
\zeta = Ri_g - \Delta Ri_g^2 + \left(\tfrac32\Delta^2-\tfrac12 c_1\right)Ri_g^3+O(Ri_g^4).
\]
This expansion is needed to construct *Ri*-based shear functions without iterative ζ-solvers.

## 15. Turbulent Prandtl Number and Curvature
\[
Pr_t=\frac{\phi_h}{\phi_m} = (1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{\alpha_m}.
\]
Log derivative:
\[
\frac{1}{Pr_t}\frac{dPr_t}{d\zeta}= \frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{\alpha_m\beta_m}{1-\beta_m\zeta}.
\]
Near neutrality,
\[
Pr_t(\zeta)=1+(\alpha_h\beta_h-\alpha_m\beta_m)\zeta+O(\zeta^2).
\]
Note
\[
\Delta=(\alpha_h\beta_h-\alpha_m\beta_m)-(\alpha_m\beta_m).
\]
Thus deviation of \(Pr_t\) from unity partitions \(\Delta\) into a *Prandtl* component and a *momentum‑only* residual; curvature neutrality (\(\Delta=0\)) does not require \(Pr_t=1\).

## 16. Flux Richardson Number \(R_f\)
For MOST:
\[
R_f = \frac{-\zeta}{\phi_m^2/\phi_h}= \frac{-\zeta}{F(\zeta)}.
\]
Hence
\[
Ri_g = -R_f F(\zeta)^2.
\]
Near neutrality \(F\simeq 1+\Delta\zeta\); eliminating \(\zeta\) by inversion:
\[
R_f \simeq -Ri_g \left[1-2\Delta Ri_g + O(Ri_g^2)\right].
\]
Thus the *ratio* \(Ri_g/R_f \to -1\) with linear correction \(2\Delta Ri_g\).

## 17. Singular Asymptotics
As \(\zeta\to \zeta_h^-=1/\beta_h\) (heat singular first, if \(\beta_h>\beta_m\)):
\[
\phi_h\sim (\beta_h(\zeta_h-\zeta))^{-\alpha_h},\quad
F \sim (\beta_h(\zeta_h-\zeta))^{-\alpha_h}(1-\beta_m\zeta_h)^{2\alpha_m},
\]
\[
V_{\log}\sim \frac{\alpha_h}{\zeta_h-\zeta},\quad
W_{\log}\sim \frac{\alpha_h}{(\zeta_h-\zeta)^2}.
\]
Hence
\[
\frac{d^{2}Ri_g}{d\zeta^{2}} \sim (\beta_h(\zeta_h-\zeta))^{-\alpha_h}(1-\beta_m\zeta_h)^{2\alpha_m}
\left[
\frac{2\alpha_h}{\zeta_h-\zeta} + O(1)
\right]
\]
→ algebraic blow‑up exponent \(\alpha_h+1\). Provides *a priori* filter threshold for numerical stability.

## 18. Uniform Composite Approximation (Two-Term)
Define outer (regular) part \(Ri^{(o)}(\zeta)=\zeta(1+\Delta \zeta)\), inner stretched coordinate near singular heat pole: \(\eta=\frac{\zeta_h-\zeta}{\zeta_h}\). Inner form:
\[
Ri^{(i)}(\eta) \sim \zeta_h (1-\eta) (\beta_h\zeta_h \eta)^{-\alpha_h}C_m,
\quad C_m=(1-\beta_m\zeta_h)^{2\alpha_m}.
\]
Uniform (first order):
\[
Ri_g^{(u)}(\zeta)= Ri^{(o)} + Ri^{(i)} - \text{overlap},\quad \text{overlap} \sim \zeta_h\,( \beta_h\zeta_h)^{-\alpha_h} C_m.
\]
Improves analytic continuation toward large ζ without direct evaluation near singular point.

## 19. Error Bounds for Binomial / Curvature Series
Binomial remainder for \((1-\beta\zeta)^{-\alpha}=\sum_{n=0}^{N} \binom{\alpha+n-1}{n} (\beta\zeta)^n + R_{N+1}\).
With \(|\beta\zeta|<1\),
\[
|R_{N+1}|\le \binom{\alpha+N}{N+1} \frac{|\beta\zeta|^{N+1}}{(1-|\beta\zeta|)^{\alpha+N+1}}.
\]
Propagating to curvature (linear combination of products and rational factors) yields practical truncation target
\[
N \gtrsim \frac{\log(\varepsilon)}{\log|\beta \zeta|}-1,
\]
for tolerance \(\varepsilon\), provided \(|\beta\zeta|<\rho<1\). Use the *larger* of \(\beta_h,\beta_m\) for conservative choice.

## 20. Algorithmic Notes (Model Implementation)
1. Precompute \(\Delta, c_1\) for parameter set.
2. Use neutral series up to desired order if \(|\zeta|<\zeta_{th}\) (e.g. 0.05).
3. Else evaluate exact rational/log form guarding \(1-\beta_{(\cdot)}\zeta\).
4. Optional: If \(\Delta c_1<0\) and predicted \(\zeta_{\text{inf}}\) within layer, tabulate sign change height for diagnostic output.
5. For Ri → ζ mapping in Ri-based closures employ 3rd-order inversion above; refine by single Newton step on \(Ri_g(\zeta)-Ri_{target}=0\).

## 21. Summary Table (Key Neutral Coefficients)
\[
\begin{array}{l|l}
\text{Quantity} & \text{Neutral Value} \\
\hline
Ri_g/\zeta & 1 \\
\partial_\zeta^2 Ri_g & 2\Delta \\
\partial_\zeta^3 Ri_g & 3(\Delta^2-c_1) \\
Pr_t & 1 \\
\partial_\zeta Pr_t & (\alpha_h\beta_h-\alpha_m\beta_m) \\
Ri_g/R_f & -1 \\
\end{array}
\]

These coefficients give direct calibration handles for assessing parameter plausibility against high-resolution or LES diagnostics.

## 21A. Parameter Selection and Curvature Scaling
Typical MOST fits: \(\alpha_{m,h}\approx 0.45\text{–}0.55\), \(\beta_{m,h}\approx 14\text{–}16\).
Neutral curvature depends only on \(\Delta=\alpha_h\beta_h-2\alpha_m\beta_m\).
Example (symmetric choice): \(\alpha_m=\alpha_h=0.5,\ \beta_m=\beta_h=16\):
\[
\Delta=0.5\cdot 16 -2(0.5\cdot16)=8-16=-8,\qquad
\left.\partial_{\zeta}^2 Ri_g\right|_{0}=2\Delta=-16.
\]
Sign: negative ⇒ initial concave‑down (momentum corrections dominate).
Domain constraint: \(1-\beta\zeta>0\Rightarrow \zeta<1/\beta\approx 0.0625\); large \(\beta\) compresses usable ζ range while amplifying \(|\Delta|\).
Height curvature adds \(1/L^{2}\):
\[
\left.\partial_{z}^2 Ri_g\right|_{0}=\frac{-16}{L^{2}}.
\]
If \(\beta_h\neq\beta_m\), adjust: \(\Delta=\alpha_h\beta_h-2\alpha_m\beta_m\) (one large \(\beta\) can flip sign).
Series validity: with \(\rho=\max(\beta_m,\beta_h)\zeta\), choose ζ so \(\rho<0.9\) for tight error bounds; at \(\beta=16,\ \zeta=0.02\Rightarrow \rho=0.32\) small ⇒ rapid convergence.
Relative neutral curvature magnitude scales linearly with chosen \(\beta\); picking \(\beta=16\) yields integer-friendly coefficients in binomial expansion indexes.

Error bound (neutral series truncation after N terms, \(\rho<1\)):
\[
|R_{N+1}|\lesssim \frac{\rho^{N+1}}{(N+1)(1-\rho)}.
\]

## 22. Enhanced MOST / Ri Formulations (Model-Oriented)
Targets: remove singularities, enable Ri-based closure, incorporate nonlocal effects, stabilize curvature.

1. Regularized power law (RPL):
\[
\phi=(1+\gamma(\beta\zeta))^{\alpha},\quad \gamma(x)=\frac{x}{1+\delta x},\;\delta>0.
\]

2. Variable exponent (VEXP):
\[
\phi_m=(1-\beta_m\zeta)^{-\alpha_m(1+\eta_m\zeta)},\;
\phi_h=(1-\beta_h\zeta)^{-\alpha_h(1+\eta_h\zeta)}.
\]

3. Ri-conditioned blend (RB):
\[
\phi_m^{RB}=\phi_m(1-\chi)+\phi_m^{\infty}\chi,\quad
\chi(Ri)=\frac{Ri^p}{Ri^p+Ri_c^p}.
\]

4. Dynamic turbulent Prandtl:
\[
Pr_t=1+a_1 Ri + a_2 Ri^2,\;\phi_h=Pr_t\phi_m.
\]

5. Nonlocal augmentation:
\[
\phi_{m,h}^{NLM}=\phi_{m,h}\left(1+c_{m,h}\frac{z}{h_{mix}}\right).
\]

6. Curvature guard:
\[
\phi_{m,h}^{guard}=\frac{\phi_{m,h}}{1+\epsilon_{m,h}|\zeta^2 \partial_\zeta^2 Ri_g|}.
\]

7. Unified Ri closure:
\[
f_m(Ri)=\left(1+b_m\frac{Ri}{Ri_c}\right)^{-e_m},\; e_m=\frac{\alpha_m}{2\alpha_m-\alpha_h}.
\]

Pseudocode sketch:
```python
def stability_functions(zeta, p):
    if p.model=='RPL':
        g=(p.beta*zeta)/(1+p.delta*p.beta*zeta)
        return (1+g)**p.alpha
    elif p.model=='VEXP':
        return (1-p.beta*zeta)**(-p.alpha*(1+p.eta*zeta))
    # ...
```

Calibration steps:
- Fit neutral curvature and Pr_t slope.
- Enforce smooth $\partial_\zeta^2 Ri_g$ profile.
- Constrain high-Ri asymptote to LES/tower regime statistics.
- Validate inversion Ri ↔ ζ against power-law baseline.

Expected gains: fewer time-step restrictions (pole removed), reduced premature flux collapse in stable BL, more accurate Ri classification over coarse grids.

## 23. Appendix—Planetary Scaling and Polar Use
What changes off‑Earth
- Curvature form in ζ is unchanged:
  \[
  \partial_\zeta^2 Ri_g = F[2V_{\log}+\zeta(V_{\log}^2-W_{\log})].
  \]
- Physical curvature scales as
  \[
  \partial_z^2 Ri_g=\frac{1}{L^2}\partial_\zeta^2 Ri_g,\quad
  L=-\frac{u_*^3\,\theta_{\mathrm{ref}}}{\kappa\,g\,w'\theta'_v}.
  \]
  Planet‑specific g and θ_v (composition/humidity) only enter via L (and θ_v definition).
- When θ_v is ill‑defined (gas giants): use θ and N²; Ri_g and J=N²/S² provide parallel diagnostics; interpret ζ with an effective L from cloud‑top fluxes.

Polar regions
- Strong f, seasonal insolation gradients, katabatics: modify shear S and N², not the ζ‑curvature form.
- Analyze curvature sign/inflection vs. latitude/season to tag regime transitions (e.g., polar night jets).

Planet notes
- Mars (low g): larger |ζ| windows; dust radiative coupling alters θ′ budgets—treat L with radiative‑corrected heat flux.
- Venus (dense CO₂): small L near surface ⇒ amplified \(\partial_z^2 Ri_g\).
- Titan: redefine θ_v with methane; recompute Δ,c₁ and neutral curvature.

## 24. Quadratic SBL Truncation (Q‑SBL)
Goal: a pole‑free stable‑regime approximation for ζ>0 that preserves neutral slope/curvature and remains well‑behaved for large ζ.

Quadratic φ model (stable ζ≥0):
\[
\phi_m^{Q}(\zeta)=1+a_m\zeta+b_m\zeta^2,\qquad
\phi_h^{Q}(\zeta)=1+a_h\zeta+b_h\zeta^2,
\]
with coefficients mapped from power‑law MOST so that near‑neutral behavior is matched:
\[
a_m=\alpha_m\beta_m,\quad b_m=\tfrac12\alpha_m(\alpha_m+1)\beta_m^2,\qquad
a_h=\alpha_h\beta_h,\quad b_h=\tfrac12\alpha_h(\alpha_h+1)\beta_h^2.
\]

Resulting ratio and Ri_g to O(ζ^2/O(ζ^3)):
- Using the known series (already derived above),
\[
F(\zeta)=\frac{\phi_h}{\phi_m^2}\approx 1+\Delta\,\zeta+\tfrac12(\Delta^2+c_1)\zeta^2,
\]
\[
Ri_g^{Q}(\zeta)=\zeta\,F(\zeta)=\zeta+\Delta\zeta^2+\tfrac12(\Delta^2+c_1)\zeta^3,
\]
with
\[
\Delta=\alpha_h\beta_h-2\alpha_m\beta_m,\quad c_1=\alpha_h\beta_h^2-2\alpha_m\beta_m^2.
\]

Curvature (exact for the cubic Ri_g^Q):
\[
\frac{d^2 Ri_g^{Q}}{d\zeta^2}=2\Delta+3(\Delta^2-c_1)\,\zeta.
\]
Height‑coordinate curvature: ∂^2/∂z^2=(1/L^2)∂^2/∂ζ^2.

Domain and guards (stable SBL):
- Use on 0\le \zeta \le \zeta_{\max}, with \zeta_{\max}\in[0.2,\,0.5] (site dependent).
- Monotonic/convex stability: choose (a_{m,h}\ge0, b_{m,h}\ge0); optionally cap φ growth for extreme stability:
  \[
  \phi_{m,h}^{Q,\text{cap}}(\zeta)=\min\big(\phi_{m,h}^{Q}(\zeta),\,1+c_{m,h}\zeta\big)
  \]
  with c_{m}\approx c_{h}\approx 5 as a conservative linear cap (Businger–Dyer style).
- Blending (optional) for C^0 continuity:
  \[
  \phi_{(\cdot)}(\zeta)=
  \begin{cases}
  \phi_{(\cdot)}^{\text{power}}(\zeta), & \zeta\le \zeta_b,\\
  \phi_{(\cdot)}^{Q}(\zeta), & \zeta_b<\zeta\le \zeta_{\max},
  \end{cases}
  \]
  with small \(\zeta_b\sim 0.05\); both agree up to O(ζ^2) so mismatch is minimal.

Minimal code sketch (stable ζ):
```python
# Q-SBL phi and Ri_g/curvature (ζ≥0), with clamp
def qsbl_coeffs(alpha_m, beta_m, alpha_h, beta_h):
    am = alpha_m*beta_m
    bm = 0.5*alpha_m*(alpha_m+1)*beta_m**2
    ah = alpha_h*beta_h
    bh = 0.5*alpha_h*(alpha_h+1)*beta_h**2
    Delta = ah - 2*am
    c1 = beta_h**2*alpha_h - 2*beta_m**2*alpha_m
    return am,bm,ah,bh,Delta,c1

def phi_qsbl(zeta, a, b, zeta_max=0.5, cap_c=None):
    z = min(max(zeta, 0.0), zeta_max)
    val = 1.0 + a*z + b*z*z
    if cap_c is not None:
        val = min(val, 1.0 + cap_c*z)  # optional linear cap
    return val

def rig_qsbl(zeta, Delta, c1, zeta_max=0.5):
    z = min(max(zeta, 0.0), zeta_max)
    # cubic Ri_g with matched neutral slope/curvature
    Ri = z + Delta*z*z + 0.5*(Delta*Delta + c1)*z*z*z
    curv = 2.0*Delta + 3.0*(Delta*Delta - c1)*z
    return Ri, curv
```

When z/L » 0 (very stable), prefer the capped quadratic or switch to a prescribed asymptote to avoid excessive growth while maintaining positive φ and realistic Ri classification.

## 25. Multi‑Profile Extension: Recipe and Candidate φ Sets
The curvature expression
\[
\partial_\zeta^2 Ri_g = F\,[2V_{\log}+\zeta(V_{\log}^2-W_{\log})],\quad
F=\frac{\phi_h}{\phi_m^2}
\]
holds for any differentiable \(\phi_{m,h}(\zeta)\) in the MOST class. To repeat the analysis for other profiles, define:
\[
V_{\log}=\frac{1}{\phi_h}\phi_h' - 2\frac{1}{\phi_m}\phi_m',\qquad
W_{\log}=\frac{dV_{\log}}{d\zeta}.
\]
Then reuse all neutral limits, series, inversion, and sensitivity steps by substituting the new \(\phi_{m,h}\).

Candidate profiles (examples)
- Power‑law (Businger–Dyer; baseline, already derived)
  - Stable/unstable via signs in ζ:
    \(\phi_m=(1-\beta_m\zeta)^{-\alpha_m},\ \phi_h=(1-\beta_h\zeta)^{-\alpha_h}\).
- Beljaars–Holtslag (stable surrogate)
  - Typical polynomial/rational stable forms (site‑tuned):
    \(\phi_m=1+a_m\zeta+b_m\zeta^2,\ \phi_h=1+a_h\zeta+b_h\zeta^2\) (Q‑SBL already included); or capped linear variants.
- Cheng–Brutsaert (all‑stability, monotone, pole‑free)
  - Example pattern (parameters to be fit):
    \(\phi_m=(1+\gamma_m |\zeta|)^{p_m},\ \phi_h=(1+\gamma_h |\zeta|)^{p_h}\).
- Grachev–Fairall (very stable tail)
  - Use enhanced Prandtl: \(Pr_t=1+c_1 Ri + c_2 Ri^2\), then \(\phi_h=Pr_t\,\phi_m\) with \(\phi_m\) from chosen base.

Procedure (per profile)
1) Specify \(\phi_{m,h}(\zeta)\), compute \(\phi'_{m,h}\) and \(V_{\log}, W_{\log}\).
2) Evaluate neutral coefficients: \(\Delta=V_{\log}(0)\), \(c_1=W_{\log}(0)\) → neutral curvature \(2\Delta\), cubic term via \(\Delta^2\pm c_1\).
3) Reuse sections 5, 11, 12–16 for series, chain rule to z, inflection, inversion ζ(Ri), and diagnostics.
4) Validate against the same acceptance band for ζ and report the same minimal set (neutral curvature, first inflection, series inversion).

Pluggable evaluation sketch (extended profiles + ζ↔Ri)
```python
import math

def _central_diff(f, x, h=1e-6):  return (f(x+h)-f(x-h))/(2*h)
def _second_diff(f, x, h=1e-6):   return (f(x+h)-2*f(x)+f(x-h))/(h*h)

def make_profile(tag, pars):
    tag = tag.upper()
    if tag == 'BD_PL':
        am, bm, ah, bh = pars['am'], pars['bm'], pars['ah'], pars['bh']
        return (lambda z: (1 - bm*z)**(-am),
                lambda z: (1 - bh*z)**(-ah))
    if tag == 'BD_CLASSIC':
        a=pars.get('a',16.0); cm=pars.get('cm',5.0); ch=pars.get('ch',7.0)
        pm_e=pars.get('pm_exp',-0.25); ph_e=pars.get('ph_exp',-0.5)
        phi_m = lambda z: (1 - a*z)**(pm_e) if z<0 else (1 + cm*z)
        phi_h = lambda z: (1 - a*z)**(ph_e) if z<0 else (1 + ch*z)
        return phi_m, phi_h
    if tag == 'HOG88':
        cm=pars.get('cm',5.0); ch=pars.get('ch',7.8); c0h=pars.get('c0h',0.95)
        return (lambda z: 1 + cm*z, lambda z: c0h + ch*z)
    if tag == 'QSBL':
        am, bm, ah, bh = pars['am'], pars['bm'], pars['ah'], pars['bh']
        return (lambda z: 1 + am*z + bm*z*z,
                lambda z: 1 + ah*z + bh*z*z)
    if tag == 'CB':
        gm, pm, gh, ph = pars['gm'], pars['pm'], pars['gh'], pars['ph']
        return (lambda z: (1 + gm*abs(z))**pm,
                lambda z: (1 + gh*abs(z))**ph)
    if tag == 'RPL':
        am,bm,dm = pars['alpha_m'],pars['beta_m'],pars['delta_m']
        ah,bh,dh = pars['alpha_h'],pars['beta_h'],pars['delta_h']
        def g(b,d,z): return (b*z)/(1 + d*b*z)
        return (lambda z: (1 + g(bm,dm,z))**am,
                lambda z: (1 + g(bh,dh,z))**ah)
    if tag == 'VEXP':
        am,bm,em = pars['alpha_m'],pars['beta_m'],pars['eta_m']
        ah,bh,eh = pars['alpha_h'],pars['beta_h'],pars['eta_h']
        return (lambda z: (1 - bm*z)**(-am*(1 + em*z)),
                lambda z: (1 - bh*z)**(-ah*(1 + eh*z)))
    if tag == 'DTP':
        base_tag = pars['base_tag']; base_pars = pars['base_pars']
        return make_profile(base_tag, base_pars)
    if tag == 'URC':
        b_m, Ri_c, e_m = pars['b_m'], pars['Ri_c'], pars['e_m']
        fm = lambda Ri: (1 + b_m*Ri/Ri_c)**(-e_m)
        fh = None
        if 'b_h' in pars and 'e_h' in pars:
            fh = lambda Ri: (1 + pars['b_h']*Ri/pars.get('Ri_c_h', Ri_c))**(-pars['e_h'])
        return fm, fh
    raise ValueError(f'unknown tag {tag!r}')

def F_from(phi_m, phi_h): return lambda z: phi_h(z)/(phi_m(z)**2)
def ri_from_zeta(z, phi_m, phi_h): return z*F_from(phi_m, phi_h)(z)

def zeta_from_ri_series(Ri, Delta, c1):
    return Ri - Delta*Ri*Ri + (1.5*Delta*Delta - 0.5*c1)*(Ri**3)

def zeta_from_ri_newton(Ri_target, phi_m, phi_h, z0, tol=1e-10, maxit=20):
    F = F_from(phi_m, phi_h); z=z0
    for _ in range(maxit):
        Vlog = _central_diff(lambda zz: math.log(F(zz)), z)
        f  = z*F(z) - Ri_target
        fp = F(z) + z*F(z)*Vlog
        if fp == 0: break
        dz = f/fp; z -= dz
        if abs(dz) < tol: break
    return z

def ri_to_phi_wrappers(tag, pars, Delta=None, c1=None):
    if tag.upper() == 'URC':
        return make_profile(tag, pars)
    phi_m, phi_h = make_profile(tag, pars)
    def zeta_of_Ri(Ri):
        z0 = zeta_from_ri_series(Ri, Delta, c1) if (Delta is not None and c1 is not None) else Ri
        return zeta_from_ri_newton(Ri, phi_m, phi_h, z0)
    if tag.upper() == 'DTP':
        base_tag=pars['base_tag']; base_pars=pars['base_pars']
        a1=pars.get('a1',0.0); a2=pars.get('a2',0.0)
        base_m, _ = make_profile(base_tag, base_pars)
        def fm(Ri): z=zeta_of_Ri(Ri); return base_m(z)
        def fh(Ri): z=zeta_of_Ri(Ri); return (1+a1*Ri+a2*Ri*Ri)*base_m(z)
        return fm, fh
    def fm(Ri): z=zeta_of_Ri(Ri); return phi_m(z)
    def fh(Ri): z=zeta_of_Ri(Ri); return phi_h(z)
    return fm, fh

# Curvature in ζ for any φ-set (analytic or numeric derivatives)
def rig_curvature_generic(zeta, phi_m, phi_h, L):
    pm = phi_m(zeta); ph = phi_h(zeta)
    dpm = _central_diff(phi_m, zeta); dph = _central_diff(phi_h, zeta)
    Vlog = (dph/ph) - 2*(dpm/pm)
    # finite-diff W_log via V_log deriv
    Vp = ( ( _central_diff(phi_h, zeta+1e-6)/(phi_h(zeta+1e-6)) )
           - 2*( _central_diff(phi_m, zeta+1e-6)/(phi_m(zeta+1e-6)) ) )
    Vm = ( ( _central_diff(phi_h, zeta-1e-6)/(phi_h(zeta-1e-6)) )
           - 2*( _central_diff(phi_m, zeta-1e-6)/(phi_m(zeta-1e-6)) ) )
    Wlog = (Vp - Vm)/(2e-6)
    F = ph/(pm*pm)
    curv_zeta = F*(2*Vlog + zeta*(Vlog*Vlog - Wlog))
    return curv_zeta, curv_zeta/(L*L)
````