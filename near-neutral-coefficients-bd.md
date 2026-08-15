## Near-neutral coefficients from Businger–Dyer

For the classic Businger–Dyer (BD) stable formulations near \(\zeta\to 0^+\):
\[
\phi_m = 1 + a_m \zeta,\quad \phi_h = 1 + a_h \zeta,
\]
with commonly used values \(a_m\approx 4.7\) and \(a_h\approx 7.8\). If you retain only the linear terms (i.e., \(b_m=b_h=0\)), then
\[
\Delta = a_h - 2a_m = 7.8 - 9.4 = -1.6,\qquad c_1 = b_h - 2b_m = 0.
\]
Implications:

- **Neutral curvature:** \(\left.\frac{d^2 Ri_g}{d\zeta^2}\right|_0 = 2\Delta = -3.2\) → concave-down.
- **Bias direction:** concave-down implies \(Ri_b < Ri_g(z_g)\) and a bulk-versus-point bias factor \(B>1\).
- **ζ inversion (cubic):**
  \[
  \zeta \approx Ri_g - \Delta Ri_g^2 + \Big(\tfrac32\Delta^2 - \tfrac12 c_1\Big)Ri_g^3
  = Ri_g + 1.6\,Ri_g^2 + 3.84\,Ri_g^3.
  \]

For unstable BD, \(\phi\) functions are non-linear (e.g., \(\phi_m=(1-16\zeta)^{-1/4}\), \(\phi_h=(1-16\zeta)^{-1/2}\)). A Taylor expansion about \(\zeta=0^-\) yields finite linear coefficients as well:
\[
\phi_m \approx 1 + 4\zeta + 10\zeta^2 + \dots,\quad
\phi_h \approx 1 + 8\zeta + 48\zeta^2 + \dots,
\]
so near-neutral on the unstable side you’d have \(a_m\approx 4\), \(a_h\approx 8\), giving \(\Delta\approx 0\) (specifically \(\Delta = 8 - 2\cdot 4 = 0\)), i.e., weak curvature in the immediate neutral limit and rapidly increasing nonlinearity at larger \(|\zeta|\). If you prefer other empirical constants (e.g., 5 and 5 for modified BD), update \(a_m,a_h\) and recompute \(\Delta,c_1\) accordingly.

---

## Curvature mapping: ζ to z

- **If L is uniform in the layer:**
  - **Label:** Derivative scaling
  - \[
    \frac{d^2 Ri_g}{dz^2} = \frac{1}{L^2}\frac{d^2 Ri_g}{d\zeta^2}
    \]
  - **Result:** Same concavity and bias logic; only magnitude rescales by \(1/L^2\).

- **If L varies with z:**
  - **Label:** Extra terms
  - \[
    \frac{d^2 Ri_g}{dz^2} = \frac{1}{L^2}\frac{d^2 Ri_g}{d\zeta^2}
    - 2\frac{1}{L^3}\frac{dL}{dz}\frac{d Ri_g}{d\zeta}
    - \text{terms with } \frac{d^2 L}{dz^2}
    \]
  - **Result:** Curvature in \(z\) combines profile shape and stability variation; using \(\zeta\) isolates the MOST shape. Practically, use \(\zeta\)-curvature for diagnostics and treat \(L(z)\) variability via layer splitting or effective \(L\).

---

## Practical guidance

- **Near-neutral diagnostics:** Use BD-derived \(a_m,a_h\) to compute \(\Delta\) and neutral curvature \(2\Delta\). This sets the expected sign of the bulk vs point bias.
- **Bias correction:** Apply damping \(G(\zeta,\Delta z)\) designed so \(G(0)=1\) and \(G'(0)=0\), preserving neutral curvature while reducing coarse-grid bias.
- **Representative height:** Evaluate point \(Ri_g\) at \(z_g=\sqrt{z_0z_1}\) and use Simpson/trapezoid integration for \(Ri_b\) when curvature is non-negligible.
- **ζ inversion for closures:** Use \(\zeta(Ri_g)\) from the series as a seed for Newton to evaluate \(\phi_m,\phi_h\) robustly when Ri is the control variable.

---

## Quick check with BD-stable

- **Coefficients:** \(a_m=4.7,\ a_h=7.8 \Rightarrow \Delta=-1.6\).
- **Neutral curvature:** \(2\Delta=-3.2\) → concave-down; expect \(B>1\).
- **ζ seed:** \(\zeta \approx Ri_g + 1.6\,Ri_g^2 + 3.84\,Ri_g^3\).
- **Action:** Prefer Simpson for \(Ri_b\); if \(B>1.05\), apply mild damping; if \(B>1.3\) with strong inversion, reduce mixing length and consider raising \(Ri_c^*\).

If you’re using a different BD constant set, share them and I’ll plug them in to give you the updated \(\Delta\), curvature, and ζ-inversion coefficients.
