# Key Equations (MOST & Ri Curvature)

Gradient Richardson (MOST):
Ri_g(\zeta)=\zeta\,\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=\zeta F(\zeta),\quad F=\phi_h/\phi_m^2.

Near-neutral expansions (φ_x=1+a_x ζ + b_x ζ² + …):
Ri_g=\zeta+\Delta\zeta^2+\tfrac12(\Delta^2+c_1)\zeta^3+O(\zeta^4),
Δ=a_h-2a_m,\ c_1=b_h-2b_m.

Inversion:
\zeta=Ri_g - Δ Ri_g^2 + (1.5Δ^2 - 0.5 c_1)Ri_g^3 + O(Ri_g^4).

Curvature (log-derivative form):
V_log=(\phi_h'/\phi_h) - 2(\phi_m'/\phi_m),\ W_log=V_log',
d²Ri_g/d\zeta² = F [ 2 V_log + \zeta (V_log² - W_log) ].
Neutral curvature: (d²Ri_g/d\zeta²)|₀ = 2Δ.

Turbulent Prandtl:
Pr_t=\phi_h/\phi_m≈1+(a_h-a_m)\zeta ⇒ Pr_t(Ri_g)≈1+(a_h-a_m)Ri_g - (a_h-a_m)Δ Ri_g².

Bias ratio (Jensen concave-down):
B=Ri_g(z_g)/Ri_b,\ z_g=√(z_0 z_1),\ Ri_b<Ri_g(z_g) if Δ<0.

Log-mean height (shear match):
z_L=(z_1 - z_0)/ln(z_1/z_0),\ z_g ≤ z_L ≤ z_a=(z_0+z_1)/2.

Generic damping (neutral-preserving):
G(\zeta,Δz)=exp[-D (Δz/Δz_ref)^p (ζ/ζ_ref)^q ], q≥2 ⇒ G(0)=1,G'_ζ(0)=0.

Power-law φ (stable):
φ_x=(1-β_x ζ)^(-α_x), domain ζ<1/β_x.

Variable-L mapping (constant L shortcut):
d²Ri_g/dz² ≈ (1/L²) d²Ri_g/d\zeta² if |(d²ζ/dz²)(dRi_g/d\zeta)| ≪ |(d\zeta/dz)²(d²Ri_g/d\zeta²)|.

Critical Ri (dynamic concept):
Ri_c^* = Ri_c [1 + α_B (Ri_b/Ri_c -1) + α_Γ (Γ/Γ_ref -1)].

Mixing length:
l_m = κ z / φ_m,\ K_m = l_m² S f_m(Ri_g), f_m=1/φ_m(ζ(Ri_g))².

Scalar closure:
f_q(Ri_g)=1/[φ_m(ζ(Ri_g)) φ_q(ζ(Ri_g))], K_q = l_m² S f_q.

Neutral invariants to preserve in corrections:
1. 2Δ (curvature)
2. (a_h - a_m) (Pr_t slope)
3. Signs of Δ, c_1.

Tuning ranges (literature norms):
a_m ≈ 4–5, a_h ≈ 7–8 (BD stable); Δ ≈ -1 to -2; Ri_c ≈ 0.2–0.25; Pr_t neutral ≈ 0.9–1.0.


