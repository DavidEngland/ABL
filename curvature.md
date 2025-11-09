# Gradient Richardson Number Curvature: Formal Derivation and Analytical Toolkit

> Note on z vs. ζ: We seek ∂²Ri_g/∂z², but derive in ζ=z/L and convert at the end using ∂²/∂z²=(1/L²)∂²/∂ζ² (L treated locally constant).

## 1. Definitions
Let
\[
\phi_m(\zeta)=(1-\beta_m\zeta)^{-\alpha_m},\qquad
\phi_h(\zeta)=(1-\beta_h\zeta)^{-\alpha_h},\qquad \zeta=\frac{z}{L},
\]
with domains \(1-\beta_{m,h}\zeta>0\). Define
\[
Ri_g(\zeta)=\zeta\,\frac{\phi_h}{\phi_m^2}=\zeta F(\zeta),\quad
F=(1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{2\alpha_m}.
\]

## 2. Logarithmic Structure
\[
\frac{F'}{F}=V_{\log}=\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta},\quad
W_{\log}=\frac{dV_{\log}}{d\zeta}= \frac{\alpha_h\beta_h^2}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^2}{(1-\beta_m\zeta)^2}.
\]

## 2A. Unified Log-Derivative Notation
\[
v_m=\frac{\phi_m'}{\phi_m},\quad v_h=\frac{\phi_h'}{\phi_h},\quad
V_{\log}=v_h-2v_m,\quad W_{\log}=V_{\log}'=v_h'-2v_m'.
\]
Generic curvature:
\[
\frac{d^{2}Ri_g}{d\zeta^{2}}=F\left[2V_{\log}+\zeta(V_{\log}^{2}-W_{\log})\right].
\]
This matches earlier “G” (\(=V_{\log}\)) and “G′” (\(=W_{\log}\)).

## 3. Curvature
\[
\frac{d^2 Ri_g}{d\zeta^2}=F\Big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\Big].
\]
(This matches and subsumes earlier expanded forms.)

## 4. Neutral Expansions
Coefficients:
\[
\Delta = \alpha_h\beta_h-2\alpha_m\beta_m,\quad c_1=\alpha_h\beta_h^2-2\alpha_m\beta_m^2.
\]
Series:
\[
Ri_g=\zeta + \Delta \zeta^2 + \tfrac12(\Delta^2+c_1)\zeta^3 + O(\zeta^4).
\]
Curvature:
\[
\partial_{\zeta}^2 Ri_g=2\Delta +3(\Delta^2-c_1)\zeta + O(\zeta^2).
\]

## 5. Inversion
\[
\zeta = Ri_g - \Delta Ri_g^2 + \Big(\tfrac{3}{2}\Delta^2-\tfrac12 c_1\Big)Ri_g^3 + O(Ri_g^4).
\]

## 6. Inflection Criterion
Inflection solves \(2V_{\log}+\zeta(V_{\log}^2-W_{\log})=0\).
Approximate interior root (if exists):
\[
\zeta_{\text{inf}}\approx -\frac{2\Delta}{\Delta^2-c_1}.
\]

## 7. Singular Asymptotics
Heat singular (\(\beta_h>\beta_m\)):
\[
Ri_g \sim (\beta_h(\zeta_h-\zeta))^{-\alpha_h},\quad \zeta_h=1/\beta_h;
\quad \frac{d^2Ri_g}{d\zeta^2}\sim C (\zeta_h-\zeta)^{-(\alpha_h+1)}.
\]

## 8. Flux Richardson Relation
\[
R_f=-\frac{\zeta}{F},\quad Ri_g = -R_f F^2.
\]
Near neutrality: \(Ri_g=-R_f+2\Delta R_f^2+O(R_f^3)\).

## 9. Turbulent Prandtl Number
\[
Pr_t=\frac{\phi_h}{\phi_m}=(1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{\alpha_m},\quad
Pr_t=1+(\alpha_h\beta_h-\alpha_m\beta_m)\zeta+O(\zeta^2).
\]

## 10. Error Control
Binomial remainder for φ functions ensures absolute error bound for truncated curvature series; choose truncation \(N\) s.t.
\[
\max(\beta_m,\beta_h)^N |\zeta|^N < \varepsilon (1-\max \beta|\zeta|)^{\alpha_h+2\alpha_m+N}.
\]

## 11. Implementation Outline
1. If \(|\zeta|<\zeta_{th}\), use series up to target order.
2. Else evaluate rational form with guards.
3. Optional Newton refinement for ζ(Ri) starting from series inversion.
4. Report diagnostics: \(\Delta, c_1, \zeta_{\text{inf}}\) (if admissible), curvature ratio \(\mathcal{C}(\zeta)\).

## 11A. Operational / Climate Model Notes
- Provide \(V_{\log},W_{\log}\) directly for multi-profile testing (power-law, quadratic, regularized).
- Use \(W_{\log}\) sign changes to anticipate intermittent regime onset.
- Arctic amplification: elevated static stability ⇒ larger \(v_h\) gradients; curvature aids vertical diffusion constraint tuning.

## 12. Summary Identities (Ready for Coding)
\[
\boxed{
\begin{aligned}
Ri_g &= \zeta(1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{2\alpha_m},\\
\partial_{\zeta}^2 Ri_g &= (1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{2\alpha_m}
\Big( \frac{2\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{4\alpha_m\beta_m}{1-\beta_m\zeta}\\
&\quad + \zeta \big[
(\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta})^2
- (\frac{\alpha_h\beta_h^2}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^2}{(1-\beta_m\zeta)^2})
\big]\Big).
\end{aligned}}
\]

## 13. Minimal Data Products
- Neutral curvature \(2\Delta\).
- (Optional) Inflection height (if inside domain).
- Series inversion coefficients \(\Delta, \tfrac32\Delta^2-\tfrac12 c_1\).

## 14. Notes
All previous explanatory prose retained in prior document; this file is the condensed mathematical reference.

## 15. Height-Coordinate Curvature (chain rule)
Given $\zeta=z/L$,
\[
\boxed{\frac{\partial^{2}Ri_g}{\partial z^{2}}=\frac{1}{L^{2}}\;\frac{\partial^{2}Ri_g}{\partial \zeta^{2}}.}
\]
If desired, assess sensitivity to vertical variation of $L$ separately; for MOST surface-layer diagnostics we evaluate with the local $L$.

## 15A. Parameter Choice: α≈0.5, β≈14–16
Neutral curvature fixed by \(\Delta=\alpha_h\beta_h-2\alpha_m\beta_m\).
Example symmetric choice: \(\alpha_m=\alpha_h=0.5,\ \beta_m=\beta_h=16\):
\[
\Delta=-8,\quad \partial_{\zeta}^2 Ri_g|_{0}=-16,\quad \partial_{z}^2 Ri_g|_{0}=-16/L^2.
\]
Large \(\beta\) increases \(|\Delta|\) but restricts ζ domain (\(\zeta<1/\beta\)).
Series regime: with \(\beta=16\), practical ζ seldom exceeds 0.05 ⇒ \(\rho=\beta\zeta\le 0.8\) ensures fast convergence.
If \(\beta_h\neq\beta_m\) (e.g. \(\beta_h=16,\beta_m=14\)):
\[
\Delta=0.5\cdot16-2(0.5\cdot14)=8-14=-6,\quad \partial_{\zeta}^2 Ri_g|_{0}=-12.
\]
Thus curvature is immediately determined once (\(\alpha,\beta\)) pair is fixed; only vertical coordinate change introduces the \(1/L^2\) scaling.

## 15B. Quadratic SBL Truncation (Q‑SBL)
A pole‑free stable‑regime surrogate that preserves neutral coefficients:
- φ model:
  \[
  \phi_m^{Q}=1+a_m\zeta+b_m\zeta^2,\ \ \phi_h^{Q}=1+a_h\zeta+b_h\zeta^2,
  \]
  with
  \[
  a_{m,h}=\alpha_{m,h}\beta_{m,h},\quad b_{m,h}=\tfrac12\alpha_{m,h}(\alpha_{m,h}+1)\beta_{m,h}^2.
  \]
- Ri and curvature (cubic form):
  \[
  Ri_g^{Q}=\zeta+\Delta\zeta^2+\tfrac12(\Delta^2+c_1)\zeta^3,\quad
  \partial_\zeta^2 Ri_g^{Q}=2\Delta+3(\Delta^2-c_1)\zeta,
  \]
  where \(\Delta=\alpha_h\beta_h-2\alpha_m\beta_m,\ c_1=\alpha_h\beta_h^2-2\alpha_m\beta_m^2\).

Use on ζ∈[0,ζ_max] with ζ_max≈0.2–0.5; enforce monotone/convex behavior via a_{m,h},b_{m,h}≥0 and (optionally) a linear cap \(\phi\le 1+c\,\zeta\) (e.g., c≈5). Blend with the exact power law for ζ≤ζ_b to maintain continuity.

Implementation note
- Height‑coordinate curvature: ∂^2/∂z^2=(1/L^2)∂^2/∂ζ^2.
- This Q‑SBL matches neutral slope/curvature and avoids the finite‑height pole, stabilizing SBL integrations.

## 16. Enhanced MOST / Ri Function Candidates (Implementation Focus)
| Tag | Form | Purpose |
|-----|------|---------|
| RPL | $(1+\frac{\beta\zeta}{1+\delta\beta\zeta})^{\alpha}$ | Remove singularity |
| VEXP | $(1-\beta\zeta)^{-\alpha(1+\eta\zeta)}$ | Tune curvature slope |
| RB | Blend MOST with shear-asymptote via $\chi(Ri)$ | Smooth regime transition |
| DTP | $Pr_t=1+a_1 Ri+a_2 Ri^2$ | Dynamic heat/momentum ratio |
| NLM | $\phi(1 + c z/h_{mix})$ | Nonlocal mixing depth |
| GUARD | $\phi/(1+\epsilon|\zeta^2 \partial_\zeta^2 Ri_g|)$ | Dampen spikes |
| URC | $f_m(Ri)=(1+b_m Ri/Ri_c)^{-e_m}$ | Direct Ri-based closure |

Neutral calibration targets:
- Curvature: $2(\alpha_h\beta_h-2\alpha_m\beta_m)$
- Prandtl slope: $(\alpha_h\beta_h-\alpha_m\beta_m)$
- Optional inflection: $\zeta_{inf}$

Minimal integration pseudocode (extensible):
```python
def phi_m(zeta, prm):
    if prm.tag=='RPL':
        g=(prm.beta_m*zeta)/(1+prm.delta*prm.beta_m*zeta)
        return (1+g)**prm.alpha_m
    elif prm.tag=='URC':
        # Ri provided externally
        return (1+prm.b_m*prm.Ri/prm.Ri_c)**(-prm.e_m)
    # ...
```

Validation metrics:
- Bias reduction in modeled Ri vs observed (tower/lidar).
- Stability of curvature under Δz coarsening.
- Maintained neutral flux accuracy.

## 17. Planetary Notes and Polar Focus
- Replace (g,R,c_p,θ_v) by planet values; recompute \(L=-u_*^3\theta_{\mathrm{ref}}/(\kappa g w'\theta'_v)\).
- ζ‑space formulas unchanged; physical curvature gains 1/L² scaling.
- Gas giants: use θ, N², S from cloud‑level retrievals; treat L as effective flux scale; Ri_g pairs with J=N²/S².
- Polar diagnostics: map curvature sign and \(\zeta_{\mathrm{inf}}\) vs. latitude/season; relate to jet and vortex edges.

Quick constants set (examples)
- Mars: low g, thin CO₂ → larger ζ span; strong diurnal cycles.
- Venus: dense CO₂ → small L, stable near‑surface.
- Titan: methane humidity → modify θ_v and L.

## 18. Beyond ABL—Rotating Core/MHD (Sketch)
- Dimensionless controls: Ro, E, Λ (Elsasser), Rm; define shear–stratification J=N²/S².
- Use HS tables to precompute closure/sensitivity across (Ro,E,Λ); adopt curvature‑style classifiers for regime tagging.
- Observational anchors: gravity (mass redistribution), magnetic secular variation (flow constraints).

Actionables
- Build planet packs (g, R, c_p, θ_v composition) and re‑run curvature spectra.
- Prototype MHD toy layers to test “curvature” classifiers vs. (Ro,E,Λ).

## 19. Multi‑Profile Pack (Reusable Curvature Recipe)
General recipe (any φ)
\[
F=\phi_h/\phi_m^2,\quad V_{\log}=(\phi_h'/\phi_h)-2(\phi_m'/\phi_m),\quad
W_{\log}=dV_{\log}/d\zeta,\quad
\partial_\zeta^2 Ri_g = F[2V_{\log}+\zeta(V_{\log}^2-W_{\log})].
\]

Profiles to test (examples; fit coefficients to site/LES)
- Power‑law (Businger–Dyer baseline).
- Quadratic stable (Q‑SBL) with optional cap (Section 15B).
- Cheng–Brutsaert‑style: \(\phi_{m,h}=(1+\gamma_{m,h}|\zeta|)^{p_{m,h}}\).
- Dynamic‑Prandtl: pick \(\phi_m\) baseline, then \(\phi_h=Pr_t\phi_m\) with \(Pr_t=1+a_1 Ri+a_2 Ri^2\).

Pluggable skeleton
```python
import math

def _central_diff(f, x, h=1e-6):
    return (f(x+h)-f(x-h))/(2*h)

def _second_diff(f, x, h=1e-6):
    return (f(x+h)-2*f(x)+f(x-h))/(h*h)

def make_profile(tag, pars):
    """
    Returns (phi_m(zeta), phi_h(zeta)) callables.
    Tags and expected pars (example defaults shown as comments):

    BD_PL   : power-law (Businger–Dyer power-law form)
        pars: am, bm, ah, bh   # α_m, β_m, α_h, β_h
    BD_CLASSIC : classic BD (unstable power, stable linear)
        pars: a=16.0, cm=5.0, ch=7.0, pm_exp=-0.25, ph_exp=-0.5
    HOG88  : Högström-like stable linear
        pars: cm=5.0, ch=7.8, c0h=0.95
    QSBL   : quadratic stable (Beljaars–Holtslag style)
        pars: am, bm, ah, bh   # linear & quadratic coeffs (a,b) for each
    CB     : Cheng–Brutsaert-style all-stability
        pars: gm, pm, gh, ph   # φ=(1+γ|ζ|)^p (even in sign), use sign if desired
    RPL    : regularized power law
        pars: alpha_m, beta_m, delta_m; alpha_h, beta_h, delta_h
    VEXP   : variable exponent
        pars: alpha_m, beta_m, eta_m; alpha_h, beta_h, eta_h
    DTP    : dynamic turbulent Prandtl (φ_h = Pr_t(Ri)*φ_m)
        pars: base_tag, base_pars, a1, a2
    URC    : Ri-based closure (direct f_m(Ri); φ_m only)
        pars: b_m, Ri_c, e_m; optionally b_h, e_h for φ_h
    """
    tag = tag.upper()

    if tag == 'BD_PL':
        am, bm, ah, bh = pars['am'], pars['bm'], pars['ah'], pars['bh']
        return (lambda z: (1 - bm*z)**(-am),
                lambda z: (1 - bh*z)**(-ah))

    if tag == 'BD_CLASSIC':
        a     = pars.get('a', 16.0)     # unstable |ζ|
        cm    = pars.get('cm', 5.0)     # stable slope momentum
        ch    = pars.get('ch', 7.0)     # stable slope heat
        pm_e  = pars.get('pm_exp', -0.25)  # unstable exponent momentum
        ph_e  = pars.get('ph_exp', -0.5)   # unstable exponent heat
        def phi_m(z):
            return (1 - a*z)**(pm_e) if z < 0 else (1 + cm*z)
        def phi_h(z):
            return (1 - a*z)**(ph_e) if z < 0 else (1 + ch*z)
        return phi_m, phi_h

    if tag == 'HOG88':
        # Högström (1988) – typical stable fits; tune per site
        cm = pars.get('cm', 5.0)
        ch = pars.get('ch', 7.8)
        c0h= pars.get('c0h', 0.95)  # intercept slightly below 1
        return (lambda z: 1 + cm*z,
                lambda z: c0h + ch*z)

    if tag == 'QSBL':
        # Beljaars–Holtslag-like stable polynomial
        am, bm, ah, bh = pars['am'], pars['bm'], pars['ah'], pars['bh']
        return (lambda z: 1 + am*z + bm*z*z,
                lambda z: 1 + ah*z + bh*z*z)

    if tag == 'CB':
        gm, pm, gh, ph = pars['gm'], pars['pm'], pars['gh'], pars['ph']
        return (lambda z: (1 + gm*abs(z))**pm,
                lambda z: (1 + gh*abs(z))**ph)

    if tag == 'RPL':
        am, bm, dm = pars['alpha_m'], pars['beta_m'], pars['delta_m']
        ah, bh, dh = pars['alpha_h'], pars['beta_h'], pars['delta_h']
        def g(b, d, z): return (b*z)/(1 + d*b*z)
        return (lambda z: (1 + g(bm, dm, z))**am,
                lambda z: (1 + g(bh, dh, z))**ah)

    if tag == 'VEXP':
        am, bm, em = pars['alpha_m'], pars['beta_m'], pars['eta_m']
        ah, bh, eh = pars['alpha_h'], pars['beta_h'], pars['eta_h']
        return (lambda z: (1 - bm*z)**(-am*(1 + em*z)),
                lambda z: (1 - bh*z)**(-ah*(1 + eh*z)))

    if tag == 'DTP':
        # φ_h = Pr_t(Ri)*φ_m; need Ri mapping => wrapper provided below
        # Here we return base φ assuming Pr_t will be applied externally.
        base_tag = pars['base_tag']; base_prs = pars['base_pars']
        return make_profile(base_tag, base_prs)

    if tag == 'URC':
        # φ_m(Ri) directly; φ_h optional
        b_m, Ri_c, e_m = pars['b_m'], pars['Ri_c'], pars['e_m']
        fm = lambda Ri: (1 + b_m*Ri/Ri_c)**(-e_m)
        if 'b_h' in pars and 'e_h' in pars:
            b_h, e_h = pars['b_h'], pars['e_h']
            fh = lambda Ri: (1 + b_h*Ri/pars.get('Ri_c_h', Ri_c))**(-e_h)
        else:
            fh = None
        # Return ζ-agnostic wrappers; use Ri→ζ transformers below when needed
        return fm, fh

    raise ValueError(f'unknown profile tag {tag!r}')

# === ζ ↔ Ri transformers and utilities ===

def F_from(phi_m, phi_h):
    return lambda z: phi_h(z)/(phi_m(z)**2)

def ri_from_zeta(zeta, phi_m, phi_h):
    F = F_from(phi_m, phi_h)
    return zeta * F(zeta)

def zeta_from_ri_series(Ri, Delta, c1):
    # ζ = Ri - Δ Ri^2 + (1.5Δ^2 - 0.5 c1) Ri^3
    return Ri - Delta*Ri*Ri + (1.5*Delta*Delta - 0.5*c1)*(Ri**3)

def zeta_from_ri_newton(Ri_target, phi_m, phi_h, z0, tol=1e-10, maxit=20):
    F = F_from(phi_m, phi_h)
    z = z0
    for _ in range(maxit):
        # numerical V_log, W_log
        Vlog = _central_diff(lambda zz: math.log(F(zz)), z)
        Wlog = _second_diff(lambda zz: math.log(F(zz)), z)
        f  = z*F(z) - Ri_target
        fp = F(z) + z*F(z)*Vlog
        if fp == 0: break
        dz = f/fp
        z -= dz
        if abs(dz) < tol:
            return z
    return z  # return last iterate

def ri_to_phi_wrappers(tag, pars, Delta=None, c1=None):
    """
    Returns (f_m(Ri), f_h(Ri)) built from a ζ-profile and ζ(Ri).
    If DTP is requested, applies Pr_t(Ri)=1+a1 Ri+a2 Ri^2 to φ_h.
    """
    if tag.upper() == 'URC':
        return make_profile(tag, pars)

    phi_m, phi_h = make_profile(tag, pars)
    # choose series seed if Delta,c1 provided else small-Ri seed
    def zeta_of_Ri(Ri):
        z0 = zeta_from_ri_series(Ri, Delta, c1) if (Delta is not None and c1 is not None) else Ri
        return zeta_from_ri_newton(Ri, phi_m, phi_h, z0)

    if tag.upper() == 'DTP':
        base_tag = pars['base_tag']; base_pars = pars['base_pars']
        a1, a2   = pars.get('a1', 0.0), pars.get('a2', 0.0)
        base_m, base_h = make_profile(base_tag, base_pars)
        def fm(Ri): 
            z = zeta_of_Ri(Ri); return base_m(z)
        def fh(Ri):
            z = zeta_of_Ri(Ri)
            Prt = 1 + a1*Ri + a2*Ri*Ri
            return Prt*base_m(z)
        return fm, fh

    def fm(Ri):
        z = zeta_of_Ri(Ri); return phi_m(z)
    def fh(Ri):
        z = zeta_of_Ri(Ri); return phi_h(z)
    return fm, fh

## 15D. Constant vs Structured \(L(z)\) — Simplifications

General mapping (from 15C):
\[
\partial_z^2 Ri_g=(\zeta'_z)^2\,\partial_\zeta^2 Ri_g+\zeta''_z\,\partial_\zeta Ri_g,\quad
\zeta'_z=\frac{L-zL'}{L^2},\ \zeta''_z=-\frac{2L'}{L^{2}}-\frac{zL''}{L^{2}}+\frac{2z(L')^{2}}{L^{3}}.
\]

Special cases:

1. Constant \(L=L_0\):
\[
\partial_z^2 Ri_g=\frac{1}{L_0^2}\partial_\zeta^2 Ri_g.
\]

2. Affine \(L=L_0+\lambda z\):
\[
\zeta'_z=\frac{L_0}{(L_0+\lambda z)^2},\quad
\zeta''_z=-\frac{2\lambda L_0}{(L_0+\lambda z)^3}.
\]

3. Power-law \(L=L_0(z/z_0)^p\):
\[
L'=\frac{pL}{z},\ L''=\frac{p(p-1)L}{z^2},\quad
\zeta'_z=\frac{1-p}{L},\ \zeta''_z=\frac{p(2p-1)}{zL}.
\]

4. Exponential \(L=L_0 e^{\lambda z}\):
\[
\zeta'_z=\frac{e^{-\lambda z}}{L_0}(1-\lambda z),\quad
\zeta''_z=\frac{e^{-\lambda z}}{L_0}(\lambda^{2} z -2\lambda).
\]

Neglect test:
\[
E_{\text{omit}}=\left|\frac{\zeta''_z\,\partial_\zeta Ri_g}{(\zeta'_z)^2\,\partial_\zeta^2 Ri_g}\right|.
\]
Use constant-\(L\) form if \(E_{\text{omit}}<\epsilon\) (e.g. 5%).

Perturbative small variation: \(L=L_0(1+\delta\ell(z)), |\delta|\ll1\):
\[
\zeta'_z=\frac{1}{L_0}\big[1-\delta(\ell+z\ell')\big]+O(\delta^2),\quad
\zeta''_z=\frac{-\delta}{L_0}\big[2\ell'+z\ell''-2z\ell'^2\big]+O(\delta^2).
\]

Implementation hook (switch):
```python
def height_curvature(z, L, dL, d2L, ri1, ri2):
    dzeta = (L - z*dL)/(L*L)
    d2zeta = -2*dL/(L*L) - z*d2L/(L*L) + 2*z*dL*dL/(L*L*L)
    return dzeta*dzeta*ri2 + d2zeta*ri1

def height_curvature_if_constant(L0, ri2):
    return ri2/(L0*L0)
```

Apply constant shortcut where \(|z dL/L|\) and \(E_{\text{omit}}\) both small.
