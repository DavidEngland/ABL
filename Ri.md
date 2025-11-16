# Gradient and Bulk Richardson Numbers: ζ Mapping and Curvature Basis

## 1. Definitions
Gradient:
$$Ri_g=\frac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2}.$$
Bulk (layer z₀→z₁):
$$Ri_b=\frac{g}{\theta}\frac{\Delta\theta\,\Delta z}{(\Delta U)^2}.$$

## 2. MOST Relation
$$Ri_g(\zeta)=\zeta\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=\zeta F(\zeta),\qquad \zeta=z/L.$$

## 3. Near-Neutral Series
Let
$$\phi_{m,h}=1+a_{m,h}\zeta+b_{m,h}\zeta^2+O(\zeta^3).$$
Then
$$Ri_g=\zeta + \Delta\zeta^2 + \tfrac12(\Delta^2+c_1)\zeta^3+O(\zeta^4),$$
$$\Delta=a_h-2a_m,\quad c_1=b_h-2b_m.$$

## 4. Inversion ζ(Ri)
$$\zeta=Ri_g - \Delta Ri_g^2 + \Big(\tfrac32\Delta^2 - \tfrac12 c_1\Big)Ri_g^3+O(Ri_g^4).$$
Seed for Newton refinement when evaluating φ at given Ri.

## 5. Curvature
Log derivatives:
$$V_{\log}=\frac{\phi_h'}{\phi_h}-2\frac{\phi_m'}{\phi_m},\quad W_{\log}=V_{\log}'.$$
Curvature:
$$\frac{d^2Ri_g}{d\zeta^2}=F\big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\big],\quad F=\frac{\phi_h}{\phi_m^2}.$$
Neutral:
$$\left.\frac{d^2Ri_g}{d\zeta^2}\right|_0=2\Delta.$$

## 6. Bulk vs Point Bias
Concave-down ($\Delta<0$) ⇒
$$Ri_b < Ri_g(z_g),\qquad z_g=\sqrt{z_0 z_1},\quad B=\frac{Ri_g(z_g)}{Ri_b}>1.$$

## 7. Correction Principle
Grid damping $G(\zeta,\Delta z)$ with $G(0,\Delta z)=1$, $\partial_\zeta G|_0=0$ preserves 2Δ while reducing $B$ at coarse Δz.

## 8. Critical Richardson Number
Fixed $Ri_c$ vs dynamic $Ri_c^*$ informed by curvature growth (magnitude of ζ(V_{\log}^2-W_{\log})) and inversion strength.

## 9. Key Identities
$$Pr_t=\frac{\phi_h}{\phi_m},\quad f_m(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))^2},\quad f_h(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\phi_h(\zeta(Ri_g))}.$$

## 10. Generic Scalar Closure

For any scalar $q$:
\[
f_q(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\,\phi_q(\zeta(Ri_g))},\quad K_q=l_m^2 S f_q.
\]
Schmidt number:
\[
Sc_t^{(q)}=\frac{\phi_q}{\phi_m}.
\]
Series (near-neutral):
\[
f_q\approx 1+a_q Ri_g+(b_q-a_q\Delta+2a_m a_q)Ri_g^2.
\]

Concise algorithm (near-neutral):
1. Compute Ri_g.
2. ζ ≈ Ri_g − Δ Ri_g² (cubic if needed).
3. Evaluate φ_m, φ_h; obtain K_m,K_h.

## 11. Numerical Estimation of Ri_g and Ri_b

Given discrete $z_k, U_k, \theta_k$:

**Point gradient (centered):**
$$
\partial U/\partial z \Big|_{z_k} \approx (U_{k+1} - U_{k-1})/(z_{k+1} - z_{k-1}).
$$

**Bulk Ri_b (layer $[z_0,z_1]$):**
- Definition: $Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz$.
- Trapezoid: $Ri_b \approx \frac{1}{2}[Ri_g(z_0) + Ri_g(z_1)]$.
- Simpson (3-pt): $Ri_b \approx \frac{1}{6}[Ri_g(z_0) + 4Ri_g(z_g) + Ri_g(z_1)]$.

**Representative heights:**
- $z_g = \sqrt{z_0 z_1}$ (geometric mean, midpoint in $\ln z$).
- $z_L = (z_1 - z_0)/\ln(z_1/z_0)$ (logarithmic mean, exact for $\Delta U$ in log wind).
- $z_a = (z_0 + z_1)/2$ (arithmetic mean, biases high for log profiles).

Use $z_g$ for point $Ri_g$ evaluation; use $z_L$ for exact layer-averaged gradient matching.

## 12. Summary
Series + curvature supply analytic bridge for Ri-based closures, bias diagnostics, and neutral-invariant correction design.
