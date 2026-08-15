# Physical Interpretation and Implementation Plan: Curvature of \(Ri_g\)

> Note on z vs. ζ: Target quantity is ∂²Ri_g/∂z². We derive in ζ=z/L for analytic clarity, then apply the chain rule at the end: ∂²/∂z² = (1/L²) ∂²/∂ζ² (with L treated locally constant).

## 0. Purpose
Provide a clean analytic form for \(\partial^{2}Ri_g/\partial \zeta^{2}\), its neutral limit, physical meaning, and a concrete work plan for implementation, testing, and manuscript preparation.

## 1. Definitions
Let \(\zeta = z/L\). For momentum (m) and heat (h) stability functions:
\[
\phi_m(\zeta)=(1-\beta_m \zeta)^{-\alpha_m},\qquad \phi_h(\zeta)=(1-\beta_h \zeta)^{-\alpha_h},
\]
domain: \(1-\beta_{(\cdot)}\zeta>0\) (avoid singularities at \(\zeta=1/\beta_m,1/\beta_h\)).

Gradient Richardson number (scaled form) has factor
\[
F(\zeta)=\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}=(1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{2\alpha_m}.
\]

## 1A. Unified Log-Derivatives
Previous notation used \(V_{\log}\) directly; we now define component logs:
\[
v_m=\frac{\phi_m'}{\phi_m},\quad v_h=\frac{\phi_h'}{\phi_h},\quad
V_{\log}=v_h-2v_m,\quad W_{\log}=V_{\log}'=v_h'-2v_m'.
\]
Curvature:
\[
\partial_{\zeta}^2 Ri_g=F\big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\big].
\]
Implementation: compute \(v_{m,h}\) analytically (power-law) or numerically, then form \(V_{\log},W_{\log}\) for robustness across profile families.

Logarithmic derivative components:
\[
V_{\log}=\frac{1}{\phi_h}\frac{d\phi_h}{d\zeta}-\frac{2}{\phi_m}\frac{d\phi_m}{d\zeta}
= \frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta},
\]
\[
W_{\log}=\frac{d}{d\zeta}\left(\frac{1}{\phi_h}\frac{d\phi_h}{d\zeta}-\frac{2}{\phi_m}\frac{d\phi_m}{d\zeta}\right)
= \frac{\alpha_h\beta_h^2}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^2}{(1-\beta_m\zeta)^2}.
\]

## 2. Curvature Expression
\[
\frac{\partial^2 Ri_g}{\partial \zeta^2}=F(\zeta)\Big[2V_{\log}+\zeta V_{\log}^2-\zeta W_{\log}\Big].
\]
Physical height curvature:
\[
\frac{\partial^2 Ri_g}{\partial z^2}=\frac{1}{L^2}\frac{\partial^2 Ri_g}{\partial \zeta^2}.
\]

## 3. Neutral Limit (\(\zeta\to 0\))
\[
F(0)=1,\quad V_{\log}(0)=\alpha_h\beta_h-2\alpha_m\beta_m,\quad W_{\log}(0)=\alpha_h\beta_h^2-2\alpha_m\beta_m^2,
\]
\[
\left.\frac{\partial^2 Ri_g}{\partial \zeta^2}\right|_{\zeta=0}=2(\alpha_h\beta_h-2\alpha_m\beta_m).
\]
Sign determines initial concavity:
- Positive: heat correction dominates momentum ⇒ concave up.
- Negative: momentum correction dominates ⇒ concave down.
- Near zero: nearly linear start of \(Ri_g(\zeta)\).

## 4. Singular Behavior
As \(\zeta\to 1/\beta_{(\cdot)}^{-}\), \(V_{\log},W_{\log}\) blow up like first and second order poles; truncation schemes must flag data approaching a chosen fraction (e.g. \(\zeta>0.7/\beta_{(\cdot)}\)).

## 5. Implementation Outline (Code)
1. Precompute arrays \(\zeta_i\).
2. For each \(\zeta_i\):
   - \(P_m=1-\beta_m\zeta_i,\ P_h=1-\beta_h\zeta_i\).
   - Reject if \(P_{(\cdot)}\le \epsilon\) (domain guard).
   - \(F_i=P_h^{-\alpha_h}P_m^{2\alpha_m}\).
   - \(V_i=\alpha_h\beta_h/P_h-2\alpha_m\beta_m/P_m\).
   - \(W_i=\alpha_h\beta_h^2/P_h^2-2\alpha_m\beta_m^2/P_m^2\).
   - Curvature: \(C_i=F_i(2V_i+\zeta_i V_i^2-\zeta_i W_i)\).
3. Height curvature: \(C^{(z)}_i=C_i/L^2\).
4. Store diagnostic tuple \((\zeta_i,F_i,V_i,W_i,C_i)\) for QA.

Pseudo-function sketch:
```python
def curvature_Ri_g(zeta, a_m, b_m, a_h, b_h, L, eps=1e-8):
    Pm = 1 - b_m*zeta
    Ph = 1 - b_h*zeta
    if Pm <= eps or Ph <= eps:
        return float('nan')
    F  = Ph**(-a_h) * Pm**(2*a_m)
    V  = a_h*b_h/Ph - 2*a_m*b_m/Pm
    W  = a_h*b_h**2/Ph**2 - 2*a_m*b_m**2/Pm**2
    Cζ = F*(2*V + zeta*V**2 - zeta*W)
    return Cζ / L**2
```

## 6. Grad Student Task List
1. Literature scan: summarize reported ranges of \(\alpha_m,\beta_m,\alpha_h,\beta_h\) (stable vs. unstable) + canonical forms (Businger–Dyer, Beljaars–Holtslag).
2. Data selection: produce QC filter for neutral vs. weakly stable/unstable segments; document thresholds for \(|\zeta|\).
3. Parameter estimation: fit \(\phi_m,\phi_h\) using nonlinear least squares; quantify covariance matrix for (\(\alpha,\beta\)).
4. Curvature computation: implement function above; benchmark vectorized vs. loop performance.
5. Sensitivity analysis: partial derivatives \(\partial C/\partial \alpha_{(\cdot)},\partial C/\partial \beta_{(\cdot)}\) using symbolic or automatic differentiation; rank contributions.
6. Neutral verification: confirm numerical curvature at smallest \(|\zeta|\) matches analytic \(2(\alpha_h\beta_h-2\alpha_m\beta_m)\) within tolerance.
7. Singularity proximity study: plot curvature growth as \(\zeta\) approaches \(0.6/\beta_{(\cdot)}\); define exclusion policy.
8. Error propagation: derive linearized uncertainty of curvature using parameter covariance; produce table of median, IQR.
9. Documentation: prepare figures (curvature vs. \(\zeta\); neutrality zoom; sensitivity bars) and a short methods section draft.
10. Manuscript integration: write subsection “Curvature Diagnostics” referencing role in early stability classification and model tuning.

## 7. Validation & QA
- Cross-check against numerical second differences of \(Ri_g(\zeta)\) computed directly from synthetic profiles.
- Use randomized perturbations of parameters to test robustness (Monte Carlo).
- Compare curvature sign classification with observed regime tagging (stable vs. unstable onset).

## 8. Extensions (Optional)
- Replace power-law \(\phi\) with HS-generated polylog expansions for non-neutral regimes; re-derive \(V_{\log},W_{\log}\).
- Multi-parameter regularization: penalize excessive curvature in inversion of (\(\alpha,\beta\)).

## 9. Summary
A compact curvature formula enables rapid stability diagnostics; neutral limit provides immediate physical interpretation; outlined workflow supports reproducible estimation, validation, and publication.

# Physical Interpretation of Ri_g Curvature

What curvature tells you
- Ri_g(ζ)=ζ φ_h/φ_m² rises roughly linearly near neutrality; the curvature term d²Ri_g/dζ² controls how quickly nonlinearity “turns on”.
- Neutral limit: 2(α_hβ_h − 2α_mβ_m) sets initial concavity.
  - >0: concave-up (heat corrections dominate).
  - <0: concave-down (momentum corrections dominate).
  - ≈0: quasi-linear onset (highest grid sensitivity).

Reading the terms
- F(ζ)=(1−β_hζ)^−α_h(1−β_mζ)^{2α_m} scales magnitude.
- V_log=α_hβ_h/(1−β_hζ) − 2α_mβ_m/(1−β_mζ) captures differential correction between heat and momentum.
- W_log removes intrinsic curvature double counting; singular factors 1−β_{•}ζ indicate approach to model-domain limits.

Inflection and domains
- Interior inflection (if any) solves 2V_log+ζ(V_log² − W_log)=0 and marks a curvature sign change.
- Practical masks: enforce 1−β_{•}ζ>0 and avoid ζ beyond ~0.7/max(β_m,β_h).

Operational guidance
- Use curvature-aware boundary condition at the first model level to anchor second-order behavior.
- Report: neutral curvature value; any curvature sign-change height; normalized curvature ratio
  C(ζ) = (d²Ri_g/dζ²)/[2(α_hβ_h − 2α_mβ_m)].
- For Ri-based closures, use near-neutral ζ(Ri) series + single Newton step; fall back to equal-β power laws if 2α_m≈α_h.

Validation cues
- LES: compare analytic vs diagnosed curvature; RMSE/bias across stability bins.
- Field: bootstrap CIs for curvature metrics; check neutral-limit intercepts and masked regions near 1/β.

Implications
- Correct curvature reduces grid-dependent Ri bias and stabilizes bulk/flux diagnostics.
- Provides interpretable diagnostics for parameter sets (α,β) and for QC of stability functions.

## References
- Monin, A. S., & Obukhov, A. M. (1954). Basic laws of turbulent mixing in the ground layer. Trudy Geofiz. Inst. AN SSSR, 151, 163–187.
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971). Flux–profile relationships in the ASL. J. Atmos. Sci., 28, 181–189.
- Paulson, C. A. (1970). Mathematical representation of wind and temperature profiles in unstable ASL. J. Appl. Meteor., 9, 857–861.
- Dyer, A. J. (1974). A review of flux–profile relationships. Boundary-Layer Meteorol., 7, 363–372.
- Högström, U. (1988). Non-dimensional wind and temperature profiles—A re-evaluation. Boundary-Layer Meteorol., 42, 55–78.
- Garratt, J. R. (1992). The Atmospheric Boundary Layer. Cambridge Univ. Press.
- Kaimal, J. C., & Finnigan, J. J. (1994). Atmospheric Boundary Layer Flows. Oxford Univ. Press.
- Stull, R. B. (1988). An Introduction to Boundary Layer Meteorology. Kluwer Academic.
- Li, D., Katul, G. G., & Bou-Zeid, E. (2012). Mean velocity and temperature profiles in a sheared diabatic TBL. Phys. Fluids, 24, 105105.
- Gryanik, V. M., Lüpkes, C., Grachev, A., & Sidorenko, D. (2020). Modified stable SBL functions (SHEBA). J. Atmos. Sci., 77, 2687–2716.
- Mahrt, L. (2014). Stably stratified atmospheric boundary layers. Annu. Rev. Fluid Mech., 46, 23–45.
- England, D. E., & McNider, R. T. (1995). Stability functions based upon shear functions. Boundary-Layer Meteorol., 72, 115–140.