# Equations from Pub2 Document titled "Stability Functions based upon Shear Functions" 1995 BLM

Authors:  "Me and Dick McNider"

Just equations, which will have is, will be a

#) equation number sometimes, which might have an a or b

It has come to light that may be more errors than know.
Have all the Text, but equations were not done in TeX,
most will have to be created.

Might be followed by some text with an update.
Will somehow need a database.  Thinking currently that equations need to be indexed and cross-referenced for accuracy and ease of use.  Somehow equations are inserted much like BibTeX and *.bib files somehow.

## MOST

1a) $$K_m = l_m^2 s / \phi_m^2,$$
1b) $$K_h = l_h^2 s / \phi_h^2,$$

E1–E7: Core MOST scales and diffusivities.
E8–E10: Ri_g linkage.
E11–E15: Stability families.
E16–E19: Grid & curvature corrections.
E20–E23: Tail damping, dynamic Ric*, final Km, Kh, defaults.
E24–E25: Generic scalar similarity and MOST diffusivity.
E26: Ri-based scalar stability function f_q(Ri_g)=1/(φ_m φ_q).
E27: Turbulent Schmidt number Sc_t^{(q)}=φ_q/φ_m.
E28: Near-neutral series for f_q(Ri_g).
E29: Effective diffusivity with reaction/deposition modifier χ_q.

E-Pr1: Turbulent Prandtl number
$$Pr_t=\frac{K_m}{K_h}=\frac{\phi_h}{\phi_m}.$$

E-Ri1: Gradient Richardson (MOST)
$$Ri_g(\zeta)=\zeta\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=\zeta F(\zeta).$$

E-Ri2: Near-neutral expansion
$$Ri_g=\zeta+\Delta\zeta^2+\tfrac12(\Delta^2+c_1)\zeta^3+O(\zeta^4),\quad \Delta=a_h-2a_m,\ c_1=b_h-2b_m.$$

E-Ri3: ζ inversion (series)
$$\zeta=Ri_g - \Delta Ri_g^2 + \Big(\tfrac32\Delta^2-\tfrac12 c_1\Big)Ri_g^3+O(Ri_g^4).$$

E-K1: Exchange coefficients (MOST)
$$K_m=\frac{\kappa z u_*}{\phi_m},\quad K_h=\frac{\kappa z u_*}{\phi_h}.$$

E-f1: Ri-based closure mappings
$$f_m(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))^2},\quad f_h(Ri_g)=\frac{1}{\phi_m(\zeta(Ri_g))\phi_h(\zeta(Ri_g))}.$$

E-Curv1: Curvature
$$\frac{d^2Ri_g}{d\zeta^2}=F\big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\big],\quad F=\frac{\phi_h}{\phi_m^2}.$$

E-Curv2: Neutral curvature invariant
$$\left.\frac{d^2Ri_g}{d\zeta^2}\right|_0=2\Delta.$$

E-Bias1: Bulk vs point bias ratio
$$B=\frac{Ri_g(z_g)}{Ri_b},\quad z_g=\sqrt{z_0 z_1},\quad B>1\text{ if }\Delta<0.$$

E-G1: Grid damping (neutral-preserving)
$$G(\zeta,\Delta z)=\exp\left[-D\left(\frac{\Delta z}{\Delta z_r}\right)^p\left(\frac{\zeta}{\zeta_r}\right)^q\right],\ p\ge1,q\ge2.$$

E-Rc1: Dynamic critical Richardson
$$Ri_c^*=Ri_c\left[1+\alpha_B\left(\frac{Ri_{\text{bulk}}}{Ri_c}-1\right)+\alpha_\Gamma\left(\frac{\Gamma}{\Gamma_{\text{ref}}}-1\right)\right].$$

E-Pr2: Pr_t near-neutral in Ri
$$Pr_t(Ri_g)\approx 1+(a_h-a_m)Ri_g - (a_h-a_m)\Delta Ri_g^2+O(Ri_g^3).$$


