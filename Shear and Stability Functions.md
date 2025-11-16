# Turbulent Prandtl Number, Richardson Number, and MOST-Based Curvature Corrections

## 1. Turbulent Prandtl Number
Definition:
$$Pr_t = \frac{K_m}{K_h} = \frac{\phi_h}{\phi_m} = Pr_t(\zeta)=Pr_t(Ri_g),\qquad \zeta = z/L.$$
Near-neutral expansion (generic smooth φ):
$$\phi_{m,h}=1+a_{m,h}\zeta + b_{m,h}\zeta^2 + O(\zeta^3).$$
Hence
$$Pr_t=1+(a_h-a_m)\zeta + \left[b_h-b_m - a_m a_h + a_m^2\right]\zeta^2 + O(\zeta^3).$$
In power-law case (φ=(1-βζ)^(-α)):
$$a=\alpha\beta,\quad b=\tfrac12\alpha(\alpha+1)\beta^2.$$

Linear stable (common SBL surrogate): φ_m=1+a_m ζ, φ_h=Pr_0+a_h ζ ⇒
$$Pr_t(\zeta)=\frac{Pr_0+a_h\zeta}{1+a_m\zeta}.$$
Neutral slope:
$$\left.\frac{dPr_t}{d\zeta}\right|_0 = a_h - a_m.$$

Ri dependence through ζ(Ri):
$$Ri_g(\zeta)=\zeta\frac{\phi_h}{\phi_m^2}=\zeta F(\zeta).$$
Invert (near-neutral):
$$\zeta = Ri_g - \Delta Ri_g^2 + \Big(\tfrac32\Delta^2-\tfrac12 c_1\Big)Ri_g^3 + O(Ri_g^4),$$
$$\Delta = a_h - 2a_m,\qquad c_1 = b_h - 2b_m.$$
Thus
$$Pr_t(Ri_g)=Pr_t(\zeta(Ri_g))\approx 1+(a_h-a_m)Ri_g - (a_h-a_m)\Delta Ri_g^2 + O(Ri_g^3).$$

## 2. Exchange Coefficients (MOST vs Ri Closure)
MOST forms:
$$K_m = \frac{\kappa z u_*}{\phi_m},\qquad K_h = \frac{\kappa z u_*}{\phi_h},$$
with
$$u_*=\sqrt{-\overline{u'w'}},\qquad \theta_*=-\frac{\overline{w'\theta'}}{u_*},\qquad l_m=\frac{\kappa z}{\phi_m}.$$
Mixing-length (Ri-based):
$$K_m = l^2 S f_m(Ri_g),\quad K_h = l^2 S f_h(Ri_g),\quad S=\left|\frac{\partial U}{\partial z}\right|.$$
Consistency:
$$f_m(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))^2},\quad f_h(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\phi_h(\zeta(Ri_g))}.$$
Turbulent Prandtl:
$$Pr_t(Ri_g)=\frac{f_m}{f_h}=\frac{\phi_h}{\phi_m}.$$

## 3. Richardson Numbers
Gradient:
$$Ri_g=\frac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2}=\zeta\frac{\phi_h}{\phi_m^2}.$$
Bulk (layer z₀→z₁):
$$Ri_b=\frac{g}{\theta}\frac{\Delta\theta\,\Delta z}{(\Delta U)^2}.$$
Geometric mean height:
$$z_g=\sqrt{z_0 z_1},\quad Ri_g(z_g)\ge Ri_b \text{ if } \frac{d^2Ri_g}{dz^2}<0.$$

## 4. Curvature and Neutral Invariant
Log derivatives:
$$V_{\log}= \frac{\phi_h'}{\phi_h}-2\frac{\phi_m'}{\phi_m},\qquad W_{\log}=V_{\log}'.$$
Curvature:
$$\frac{d^2Ri_g}{d\zeta^2}=F\big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\big],\quad F=\frac{\phi_h}{\phi_m^2}.$$
Neutral curvature:
$$\left.\frac{d^2Ri_g}{d\zeta^2}\right|_{0}=2\Delta,\qquad \Delta=a_h-2a_m.$$

## 5. Basis for Curvature Correction
Concave-down (Δ<0) ⇒ underestimation:
$$B=\frac{Ri_g(z_g)}{Ri_b}>1.$$
Correction must:
1. Preserve $2\Delta$ (neutral physics).
2. Reduce bias for large ζ, coarse Δz.
3. Vanish as Δz→0.

Generic damping:
$$G(\zeta,\Delta z)=\exp\left[-D\left(\frac{\Delta z}{\Delta z_r}\right)^p\left(\frac{\zeta}{\zeta_r}\right)^q\right],\quad p\ge1,q\ge2,$$
$$K_{m,h}^*=K_{m,h}G.$$

## 6. Dynamic Critical Richardson Number
$$Ri_c^*=Ri_c\Big[1+\alpha_B\Big(\frac{Ri_{\text{bulk}}}{Ri_c}-1\Big)+\alpha_\Gamma\Big(\frac{\Gamma}{\Gamma_{\text{ref}}}-1\Big)\Big],$$
controls onset of collapse; curvature growth (large |Δ|, ζ(V_{\log}^2-W_{\log})) informs tuning of $\alpha$.

## 7. Linking Pr(Ri)=Pr(ζ)
Since $Ri_g$ and $\zeta$ monotonic near neutral:
$$\frac{dPr_t}{dRi_g} = \frac{dPr_t/d\zeta}{dRi_g/d\zeta}=\frac{(a_h-a_m)+O(\zeta)}{1+2\Delta\zeta+O(\zeta^2)}.$$
Negative Δ magnifies initial sensitivity; smoothing Pr_t(Ri) tails prevents excessive damping of K_h.

## 8. MOST ABLE
“MOST ABLE” (Monin–Obukhov Similarity Theory: Analytical Boundary Layer Extensions) — umbrella for pole-free φ forms, ζ↔Ri inversion, neutral curvature–preserving grid corrections.

## 9. Summary
- $Pr_t(Ri_g)$ derives directly from φ ratio evaluated at ζ(Ri).
- Curvature invariant $2\Delta$ anchors neutral regime; corrections must leave it unchanged.
- Bias $B>1$ from concave-down Ri_g motivates grid damping factor.
- Dynamic $Ri_c^*$ and Pr_t(Ri_g) co-tuned via curvature diagnostics.

## Generic Scalar Stability Function f_q(Ri_g)

Given scalar $q$ (water vapor, pollutant, trace gas) with similarity function $\phi_q(\zeta)=\kappa z\,(\partial q/\partial z)/q_*$:
\[
K_q=\frac{\kappa z u_*}{\phi_q},\qquad K_m=\frac{\kappa z u_*}{\phi_m}.
\]
Mixing-length form ($l_m=\kappa z/\phi_m$, $S=\partial U/\partial z$):
\[
K_q=l_m^2 S f_q(Ri_g),\quad f_q(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\,\phi_q(\zeta(Ri_g))}.
\]
Turbulent Schmidt/Prandtl:
\[
Sc_t^{(q)}=\frac{K_m}{K_q}=\frac{\phi_q}{\phi_m}.
\]
Near-neutral ($\phi_x=1+a_x\zeta+b_x\zeta^2$, $\zeta=Ri_g-\Delta Ri_g^2+\dots$):
\[
f_q(Ri_g)\approx 1+a_q Ri_g+\big(b_q-a_q\Delta+2a_m a_q\big)Ri_g^2.
\]
If $\phi_q=\phi_h$ ⇒ $f_q=f_h$ and $Sc_t^{(q)}=Pr_t$.
Reactive attenuation (optional):
\[
K_q^{\text{eff}}=K_q/(1+\tau_r/\tau_t),\quad \tau_t=(\kappa z/u_*)/\phi_m.
\]