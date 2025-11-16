# Curvature of the Gradient Richardson Number $Ri_g$

> Note on $z$ vs. $\zeta$: We compute analytically in $\zeta=z/L$ (locally constant $L$) then map via $\partial^2/\partial z^2=(1/L^2)\,\partial^2/\partial \zeta^2$. For variable $L(z)$ use full chain rule (see Section 11).

## 1. Purpose
Provide a citation‑ready derivation and interpretation of $\partial^{2}Ri_g/\partial \zeta^{2}$ for power‑law MOST stability functions supporting curvature-aware correction discussions.

## 2. Definitions
$$
\zeta=z/L,\qquad
\phi_m=(1-\beta_m \zeta)^{-\alpha_m},\qquad
\phi_h=(1-\beta_h \zeta)^{-\alpha_h},\qquad
(1-\beta_{m,h}\zeta)>0.
$$
$$
Ri_g(\zeta)=\zeta\,\frac{\phi_h}{\phi_m^2}=\zeta F(\zeta),\qquad
F=(1-\beta_h \zeta)^{-\alpha_h}(1-\beta_m \zeta)^{2\alpha_m}.
$$

Historical note. The compact form $Ri_g(\zeta)=\zeta\,\phi_h(\zeta)/\phi_m(\zeta)^2$ is explicitly used in England & McNider (1995, BLM).

## 3. Curvature
$$
\frac{dRi_g}{d\zeta}=F\big(1+\zeta V_{\log}\big),\qquad
V_{\log}=\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta}.
$$
Evaluate variable-$L$ effects per the chain rule when needed.

## 4. Neutral Limit
\[
\Delta=V_{\log}(0)=\alpha_h\beta_h-2\alpha_m\beta_m,\quad
c_1=W_{\log}(0)=\alpha_h\beta_h^{2}-2\alpha_m\beta_m^{2},\quad
\left.\frac{d^{2}Ri_g}{d\zeta^{2}}\right|_{0}=2\Delta.
\]

## 5. Small-ζ Series (to \(O(\zeta^3)\))
\[
Ri_g(\zeta)=\zeta+\Delta\zeta^2+\tfrac12(\Delta^2+c_1)\zeta^3+O(\zeta^4).
\]
Inversion:
\[
\zeta = Ri_g - \Delta Ri_g^2 + \left(\tfrac32\Delta^2-\tfrac12 c_1\right)Ri_g^3 + O(Ri_g^4).
\]

## 6. Curvature Interpretation
- \(2\Delta<0\): initial concave‑down (momentum correction dominates heat).
- \(2\Delta>0\): initial concave‑up.
- Larger \(|\Delta|\): stronger early departure from linear \(Ri_g\approx \zeta\).

## 7. Inflection (if interior root exists)
Solve \(2V_{\log}+\zeta(V_{\log}^2-W_{\log})=0\). Small‑ζ approximation (if \(\Delta c_1<0\)):
\[
\zeta_{\text{inf}}\approx -\frac{2\Delta}{\Delta^2-c_1}.
\]
Discard if \(\zeta_{\text{inf}}\ge \min(1/\beta_m,1/\beta_h)\).

## 8. Domain Guard
Valid only for \(\zeta<1/\max(\beta_m,\beta_h)\). Practical safety: restrict to \(\zeta<0.7/\max(\beta)\). Near poles curvature blows up algebraically; mask before evaluation.

## 9. Singular Asymptotics (heat pole first, \(\beta_h>\beta_m\))
Let \(\zeta_h=1/\beta_h\), define \(\eta=\zeta_h-\zeta\):
\[
\phi_h \sim (\beta_h\eta)^{-\alpha_h},\quad
V_{\log}\sim \frac{\alpha_h}{\eta},\quad
\frac{d^{2}Ri_g}{d\zeta^{2}}\sim C\,\eta^{-(\alpha_h+1)},\quad C>0.
\]
Use filtering rather than evaluation as \(\eta\to0^+\).

## 10. Height Curvature (constant L)
\[
\frac{d^{2}Ri_g}{dz^{2}}=\frac{1}{L^{2}}\frac{d^{2}Ri_g}{d\zeta^{2}}.
\]

## 11. Variable \(L(z)\) Mapping
\[
\zeta(z)=\frac{z}{L(z)},\;
\frac{d\zeta}{dz}=\frac{L-zL'}{L^2},\;
\frac{d^{2}\zeta}{dz^{2}}=-\frac{2L'}{L^2}-\frac{zL''}{L^2}+\frac{2zL'^2}{L^3}.
\]
\[
\boxed{\frac{d^{2}Ri_g}{dz^{2}}=\left(\frac{d\zeta}{dz}\right)^2\frac{d^{2}Ri_g}{d\zeta^{2}}+\frac{d^{2}\zeta}{dz^{2}}\frac{dRi_g}{d\zeta}}.
\]
Omission metric (decide constant‑L shortcut):
\[
E_{\text{omit}}=\left|\frac{(d^{2}\zeta/dz^{2})(dRi_g/d\zeta)}{(d\zeta/dz)^2(d^{2}Ri_g/d\zeta^{2})}\right|.
\]
Use constant \(L\) form if \(E_{\text{omit}}<0.05\).

## 12. Reporting Set
- Neutral curvature \(2\Delta\)
- Presence/absence of interior inflection \(\zeta_{\text{inf}}\)
- Curvature amplification ratio \(A(\zeta)=| (d^{2}Ri_g/d\zeta^{2}) /(2\Delta) |\)
- Height curvature extrema \(\max|d^{2}Ri_g/dz^{2}|\)
- Omission metric fraction exceeding 0.05

## 13. Implementation Snippet
```python
def curvature_zeta(zeta, am,bm,ah,bh):
    # power-law φ
    dm = 1 - bm*zeta; dh = 1 - bh*zeta
    if dm <= 0 or dh <= 0: return float('nan')
    phi_m = dm**(-am); phi_h = dh**(-ah)
    F = phi_h / (phi_m*phi_m)
    V = (ah*bh)/dh - 2*(am*bm)/dm
    W = (ah*bh*bh)/dh**2 - 2*(am*bm*bm)/dm**2
    return F*(2*V + zeta*(V*V - W))
```

## 14. Neutral Example
\(\alpha_m=\alpha_h=0.5,\ \beta_m=\beta_h=16\) ⇒ \(\Delta=8-16=-8\), \(2\Delta=-16\) (concave‑down).

## 15. Inversion (Ri→ζ near neutral)
\[
\zeta \approx Ri_g - \Delta Ri_g^2 + (1.5\Delta^2 - 0.5 c_1) Ri_g^3.
\]
Refine with one Newton iteration on \(f(\zeta)=\zeta F(\zeta)-Ri_g\).

## 16. Summary
Curvature depends only on \(F,V_{\log},W_{\log}\); neutral behavior controlled by \(\Delta\). Proper domain guarding and optional variable‑L mapping ensure stable numerical use in coarse-grid SBL corrections.

## 17. Action Items
- Compute \(\Delta,c_1\) for candidate (α,β) sets.
- Tabulate \(A(\zeta)\) below 0.7/β.
- Evaluate \(E_{\text{omit}}\) for variable \(L(z)\) cases.

## Addendum: Pr_t Neutral Link
Since $Pr_t=\phi_h/\phi_m$ and $Ri_g=\zeta \phi_h/\phi_m^2$, neutral slopes ($a_h,a_m$) set both $Pr_t'(0)=a_h-a_m$ and curvature $d^2Ri_g/d\zeta^2|_0=2(a_h-2a_m)$; corrections must leave these unchanged.

Generic scalar mapping: $f_q(Ri_g)=1/(\phi_m\phi_q)$; neutral slopes $a_q$ enter via expansion $f_q\approx1+a_q Ri_g+(b_q-a_q\Delta+2a_m a_q)Ri_g^2$.