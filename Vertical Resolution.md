## Role of Model Vertical Resolution to Improve Understanding of Modeled AA
As discussed above, the SBL is sensitively dependent on model resolution. Under the present project, we will attempt to develop improved SBL parameterizations within the context of equation (1) by creating a correction to the short-tailed stability function fs(Ri) dependent on model resolution. Within the SBL community, the role of model resolution has been discussed as a fundamental issue with general SBL parameterizations going back to Shir and Bornstein (1977), Blackadar (1979), and McNider and Pielke (1981). More recently, Esau and Byrkjedal (2007) provided an excellent review of the resolution issue, providing several examples of the sensitivity error to model resolution. Byrkjedal et al. (2008) showed GCM performance in the Arctic deteriorated with coarser resolution.  Thus, resolution is a critical issue in global climate modeling, in that the coarse resolution in GCMs may not capture the correct strength/depth of the Arctic inversion, the entrainment warming discussed above, or the correct energy budgets (MCN12). In particular, Vignon et al. (2017) showed that while the shorter-tailed (EM) form performs best in higher resolution models, they couldn’t recommend these for GCM resolutions since they produced runaway cooling over the Antarctic Plateau.
Development of Grid-Dependent Correction Function
The improvement of parameterization within the context of (1) will be carried out by developing a grid resolution correction to the stability function fs(Ri). Within the SBL community, the role of model resolution has long been discussed as a fundamental issue with general SBL parameterizations going back to Shir and Bornstein (1977) and Blackadar (1979). As noted above, the most fundamental theoretical forms for the stability functions in the stable boundary layer are the short–tailed forms in Fig. 4 (such as England and McNider and Duynkerke), which are derived from the Monin-Obukhov surface profiles (Businger 1973).
We believe that examining the discretized form of  (1), which is how it is used in numerical models, may be a somewhat fresh approach to this long-running SBL issue.  For analytical simplicity, all the stability function forms in Fig. 4 can be well approximated by an exponential function
$$
f_s(Ri) = \exp\left(- \frac{\gamma Ri}{Ri_c}\right)
$$
where large γ yields a shorter-tailed form.  A value of =3.2 approximates the shorter tail England McNider form seen in Figure 4 and was used in the present study. While research models with high vertical resolution perform well using the more fundamental short-tail forms, coarse resolution weather forecasting and climate models using the same forms appeared not to have enough mixing. As noted above, mixing was added, often in an ad hoc manner, to improve model performance. Ideally, the results of parameterizations should not depend on the grid spacing. Thus, to preserve the integrity of the simulation in coarse grid models, fs(Ri) may need to be a function of vertical grid spacing ∆z.  This was discussed early on by Blackadar (1979), Shir and Bornstein (1977), McNider and Pielke (1981), and more recently by Walters et al (2007), though no general forms were proposed. Thus, we desire to find a correction function fcRi,∆z such that
$$
f_{c}(\Delta z) \;=\; e^{D\left(\frac{\gamma}{Ri_c}\right)\;Ri\left(1-\frac{\Delta z_r}{\Delta z}\right)}
$$

$$
f_s*f_c = e^{-\frac{\gamma}{Ri_c}Ri} * e^{D\left(\frac{\gamma}{Ri_c}\right)Ri\left(1 - \frac{\Delta z_r}{\Delta z}\right)}
$$

$$
f_{is}=f_s*f_c = e^{-\frac{\gamma}{Ri_c}Ri\left[\frac{(1-D)\Delta z + D \Delta z_r}{\Delta z} \right]}
$$
Which conveys the weighted averaging of ∆z and  ∆zr, being normalized by∆z (where
 ∆z≥∆zr). This form also conveys the dependence on D. For D=0, the fraction is 1 (no correction to the short-tail stability function), and for D=1, the exponent is adjusted by ∆zr∆z, yielding a long-tailed stability function for larger grid spacing. Thus, this form represents a weighted combination of short- and long-tail functions, with D controlling the relative contribution of each.
Figure 6 shows the behavior of fs*fc for three different values of D. While the assumption of V2 being a constant is a simplification and not likely to occur in general, it is interesting that the form of the analytical recovery reflects the need for longer-tailed forms at larger grid spacing. This is consistent with added mixing discussed by Savijarvi (2009) and Beljaars and Holtslag (1991), such as the Beljaars-Holtslag and Louis forms in Figure 4 that have been used in operational and coarse models. Also, it is consistent with the approach to increase Ric to add more mixing (McNider and Pielke 1981) in coarse grid settings in models or bulk Ri in observations (Stull 1988).

Preliminary Tests of the Analytically Derived Correction Function
In order to test the correction function in (8) we carry out a series of numerical simulations for the Arctic using the UAH single column model (see Figure 5) discussed in McN12. Note that the UAH model uses an explicit diffusion scheme to avoid numerical diffusion (see MCN12) rather than a semi-implicit scheme used in many operational models. We employ the GABLS1 single-column model inter-comparison framework (see Cuxart et al., 2006), which is based on an idealized Arctic stable boundary layer (Kosovic and Curry, 2000). The initial components of the wind are set equal to those of the geostrophic wind (u=ug, v=vg). The potential temperature (θ) equals 265 K up to 100 m, then it increases at a rate of 0.01 K/m until the domain top (400 m), where a value of 268K is reached. The surface temperature is set to 265 K initially, decreasing at a constant rate of 0.25 K/h. The value of the aerodynamic roughness length (z0) is set to 0.1 m. Note, this is a very idealized case. Actual Arctic boundary layers may have fine-scale features with more complex and deeper stable atmospheric structure.
Here we explore the sensitivity of the solutions to grid spacing and then examine whether we can use the correction function in (8) to make the solution more independent of grid resolution.  Figure 7 shows the model simulation of potential temperature after 10 hours of simulation for different grid spacings with no correction. (Case 1). Note that the mixing in the coarser grids is less than the values for the 2 m grid. This is consistent with our original formulation for Ri (see (4) above), in that Ri would be larger, indicating a reduction in fs(Ri) (since Ri would be shifted to the right as in Figure 4).  The coarser grid also shows the mixing increasing with height, but near the 250 m level, it has slightly less mixing than the finer grids. While the differences in mixing are small, the behavior of the temperature profile is quite different. Showing the sensitive non-linear behavior of the SBL (McNider et al. 1995 and Van de Wiel et al. 2002a,b).

Preliminary Tests of the Analytically Derived Correction Function
In order to test the correction function in (8) we carry out a series of numerical simulations for the Arctic using the UAH single column model (see Figure 5) discussed in McN12. Note that the UAH model uses an explicit diffusion scheme to avoid numerical diffusion (see MCN12) rather than a semi-implicit scheme used in many operational models. We employ the GABLS1 single-column model inter-comparison framework (see Cuxart et al., 2006), which is based on an idealized Arctic stable boundary layer (Kosovic and Curry, 2000). The initial components of the wind are set equal to those of the geostrophic wind (u=ug, v=vg). The potential temperature (θ) equals 265 K up to 100 m, then it increases at a rate of 0.01 K/m until the domain top (400 m), where a value of 268K is reached. The surface temperature is set to 265 K initially, decreasing at a constant rate of 0.25 K/h. The value of the aerodynamic roughness length (z0) is set to 0.1 m. Note, this is a very idealized case. Actual Arctic boundary layers may have fine-scale features with more complex and deeper stable atmospheric structure.
Here we explore the sensitivity of the solutions to grid spacing and then examine whether we can use the correction function in (8) to make the solution more independent of grid resolution.  Figure 7 shows the model simulation of potential temperature after 10 hours of simulation for different grid spacings with no correction. (Case 1). Note that the mixing in the coarser grids is less than the values for the 2 m grid. This is consistent with our original formulation for Ri (see (4) above), in that Ri would be larger, indicating a reduction in fs(Ri) (since Ri would be shifted to the right as in Figure 4).  The coarser grid also shows the mixing increasing with height, but near the 250 m level, it has slightly less mixing than the finer grids. While the differences in mixing are small, the behavior of the temperature profile is quite different. Showing the sensitive non-linear behavior of the SBL (McNider et al. 1995 and Van de Wiel et al. 2002a,b).

Figure 9 shows the result with D = 0.36 (Case 3). It shows that the different grids can be corrected to have a solution mostly like the 2m reference case. However, the lower part of the profile correction is not as improved as the D = 0.7 case in figure 11.  Note also in figure 9 (right) that the mixing coefficients have been corrected to be more in line with the 2m reference case.

mproved Specifications to the Parameter “D” in the Analytical Solution ( 8)
The analytical correction in equation (8) is appealing because of its structure and computational efficiency. It only takes two additional lines of code - one line to compute fc (8) and one to compute the product fs*fc.  The correction, which effectively yields longer tailed stability functions for larger grids is based on V2 being a constant, which in effect is for Ri being constant.  This is not in general true. The examples in figures 8 and 9 show that the correction functions can make model predictions less grid dependent. However, it also shows that different values of D can impact the correction. Given the fact that the ODE is ill-posed for linear gradients, it seems that one path might be to make D a function of the curvature in V2 or equivalently Ri (Beljaars 2020). Here we explore D as a function of curvature in Ri, i.e. D (d2 (Ri)/dz2).  Under this proposed framework, D becomes smaller as gradients in Ri becomes linear (i.e. no need for correction).

Figure 10 shows a plot of Ri for the case of “No Correction” given in Figure 7.  Note this is the average value of Ri over the length of the integration (10 hours), consistent with the average K values in figures 7-9.  Visually, the average curvature (concave upward) is largest below about 40 m, and between 40 m and 90 m, there is less curvature. Near 75 m, the curvature changes sign (concave down). Above 100 m, the curvature becomes less but remains concave downward. Note also that the curvature is smaller in the larger grids.

As noted above D=0.36 proves a correction which makes the coarser grids look more like the finer grid solution. In adding curvature, we explored D being a linear function of curvature.  We used a linear function
D = D0 + M*|(d2 (Ri)/dz2|
where D0 is the intercept (minimum value) and M is the slope.  For curvature, d2 (Ri)/dz2, we employ the absolute value of the curvature since the sign of the curvature should not be relevant. Since a value of D=0.36 as seen in Figure 9, provides an improvement, but a value of D = 0.7, as seen in Figure 8, provides too much mixing, we tried linear functions that bracket these values. Specifically, we chose through trials that D0 = 0.3 and M= 300 so that
 D = 0.3 +300* |(d2 (Ri)/dz2|  with D < 0.7 (8)
provided the best results.
Figure 11 shows the results of this curvature-dependent correction. It shows that the curvature-dependent correction for D (8), is an improvement over D=0.36. The curvature correction improves the upper part of the temperature profile, bringing it in line with the 2m grid case. Near the bottom of the profile, there seems to be a need for further correction. However, in looking at the curvature, as seen in Figure 10, it can be seen that the coarser grids do not capture the amount of curvature as seen in the smaller grids. In these simulations, we used similarity theory (England and McNider 1995) to specify heat and momentum from the surface up to the first model layer. As discussed below, a refinement in the technique in the future could utilize the curvature from surface similarity in terms of Ri rather than the grid calculated curvature employed here.

Future Work and Summary
The background on Arctic Amplification (AA) and the modeling of the SBL in coarse grid climate models suggests that treatments of the SBL and deeper SL in GCMs may be partly responsible for the spread seen in GCM predictions in high latitudes. This may also be a part of the reason models have underestimated the strength of AA. As discussed in section 2.5, if models have too much mixing in the base-state, they may not have the sensitivity to enhance GHG energy as models with less mixing. Also, given the difficulties with handling the SBL and SL in coarse grid settings, GCMs may have other fixes that add mixing and underestimate the sensitivity to GHG energy.
We believe the development of a grid correction and the initial tests presented here show promise in improving the problem of grid dependence in modeling the SBL and SL.  This may provide a path to allow coarse grid models to make improved climate predictions in the Arctic. However, the initial results show that further improvements are needed in making the analytical adjustment perform better. Also, the simple Arctic simulation carried out here with such a limited vertical domain needs to be expanded to more realistic settings where grid spacing in models may be much greater. The present paper and its grid correction function provide a potential path for making model performance better in the future. As a final caveat, we note that the approach presented here only addresses the explicit dependence of Ri on ∆z  (4). It does not address issues of unresolved structure (Blackadar, 1979) or mesoscale energy (Mahrt 1996; Steeneveld et al 2008b), which may add mixing.
We hope to carry out some future model experiments to look at specifying better correction functions and more realistic settings.  However, because of our team’s lack of experience with full climate models, we feel publishing these preliminary results can perhaps allow other investigators to make progress in making the modeling SL in Arctic environments less grid-dependent.  We outline some paths for improvement below.
Metrics for evaluating model performance
In testing the derived fc to determine the best value for variable D, we largely used a visual test for model improvement. For future, more complex recovery of D or numerical recovery of fc (see below), we believe two basic metrics should be considered.

(1) Quantification of improvements in the wind and temperature profile compared to the high-resolution model (i.e., model run at the ∆zR spacing) and to data sets (LES and observational). This can be accomplished by developing statistical estimates of the functional form that provides the best fit (RMSE, BIAS) in correcting the coarser grid profiles to agree with the 2m or (∆zR) reference grid. This is similar to the approach throughout the history of boundary layer similarity and modeling.

(2) Assess whether the model with the grid corrected stability function can preserve the transition from weak/strong stability. Vignon et al. (2017) and Acevedo et al (2021) have transition postulated as a critical element in SBL parametrizations. McNider et al. (1995) and Van de Wiel et al. (2002a,b) showed that the SBL is a non-linear dynamical system with the transition from/to strong/weak stability sensitively dependent on external parameters. Vignon et al. (2017) showed that the transition over the Antarctic was also highly dependent on the form of the stability function and was at the heart of capturing the SBL in the Antarctic. Vignon et al. used scatter plots of low-level temperature gradients versus 9 m wind speed to capture the transition. It showed that short-tailed forms best captured the transition between strongly/weakly stable conditions as found in the observations. One could use similar measures to see whether the grid corrected stability function maintains the transition for larger grids. One could also follow Acevedo et al. (2021) in calculating the windspeed, VR, where the SBL transitions from weakly to strongly stable occurs. The VR could be calculated in both the fine grid and in grid corrected runs to see if they are invariant. One might use the GABLS1 and possibly Dome C framework with varying geostrophic winds to calculate the VR for different grid spacings to see whether the VR in the fine (2m) resolution is maintained in the corrected coarser grid runs.

Numerical Solutions to Find Correction (fc)
The analytical solution derived and tested here  (8) is appealing due to its framework for providing corrected stability functions, fs, which are longer-tailed (see figure 6). Previous investigators have shown that longer-tailed forms are necessary to provide additional mixing in coarser grids (Beljaars and Holtslag 1991), Louis 1979), Savijarvi (2009). While this paper lays out a framework for correction, its assumption of a constant V2, or equivalently Ri, does not generally hold. As one moves to more complex and deeper test cases, a different approach may be required (if finding a variable D as described above is not sufficient). Another approach, would be to carry out an “on- the–fly” numerical solution of the ODE given in (5). Here, an ODE solver could be used to solve for  fcRi,∆z.  Because of the suspected dependence on curvature in  V2 ,  the stencil in the ODE solver will have to capture the curvature. The ODE might initially be solved at each time step in the model simulation. However, it might be solved less often, such as for how radiation is not called every time step in some models. The correction function, as calculated numerically, would be used to calculate the corrected stability function fs*fc.  While solving the ODE will add computational cost, it should be far less than running a model at a very high vertical resolution.  Similar metrics as described above would be used to test the numerical fc.


### Reconciliation: McNider–Biazar D versus Bias B

Original bulk-driven ODE (paper draft form):
\[
\frac{d\ln f_c}{d\ln \Delta z} = - D \left(\frac{Ri_b}{Ri_{\text{ref}}}\right)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q
\]
with layer bulk Richardson number \(Ri_b\), reference \(Ri_{\text{ref}}\) (e.g. 0.25), exponent \(q\ge1\).

Bias-based form (current framework):
\[
\frac{d\ln f_c}{d\ln \Delta z} = - \alpha (B-1)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q,\qquad
B=\frac{Ri_g(z_g)}{Ri_b},\ z_g=\sqrt{z_0 z_1}.
\]

Parameter mapping (equate RHS):
\[
\alpha (B-1) = D \frac{Ri_b}{Ri_{\text{ref}}} \;\Rightarrow\;
D = \alpha (B-1)\frac{Ri_{\text{ref}}}{Ri_b}, \qquad
B-1 = \frac{D}{\alpha}\frac{Ri_b}{Ri_{\text{ref}}}.
\]

Implication:
Since typically \(Ri_b < Ri_g(z_g)\) in concave-down stable layers (\(B>1\)), using \(Ri_b\) directly underestimates needed correction amplitude by roughly a factor \(1/B\). Driving the ODE with \((B-1)\) captures curvature-induced bulk vs point deficit explicitly.

Numeric example:
Layer: \(Ri_g(z_g)=0.60,\ Ri_b=0.30 \Rightarrow B=2.0,\ B-1=1.0\).
Take \(\alpha=1.0,\ Ri_{\text{ref}}=0.25\).
Then \(D = 1.0 \cdot 1.0 \cdot (0.25/0.30) \approx 0.83\).
If manuscript reports \(D=1.0\), equivalent \(\alpha \approx 1.2\) by inversion.

Recommended operational ODE:
\[
\frac{d\ln f_c}{d\ln \Delta z} = - \alpha (B-1)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q,\quad
f_c(\Delta z,\zeta)=\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}.
\]

Practical workflow:
1. Compute \(Ri_b\) (bulk) and \(Ri_g(z_g)\) (gradient at geometric mean).
2. Form \(B\). If \(B \le B_{\text{thr}}\) (≈1.05) set \(f_c=1\).
3. Else apply exponential or power-law form with driver \((B-1)\).

Code snippet (bulk-to-bias conversion):
```python
# Given D (bulk form) from paper and measured Ri_b, Ri_ref, alpha target:
# Recover alpha or translate D to bias-based alpha
alpha_equiv = D * (Ri_b / Ri_ref) / (B - 1.0)

# Preferred direct usage:
exponent = -alpha * (B - 1.0) * (dz/dz_ref)**p * (zeta/zeta_ref)**q
f_c = math.exp(exponent)
```

Summary:
- Bulk-driven D hides curvature strength when \(Ri_b\) is suppressed.
- Bias form with \(B-1\) is dimensionless, directly tied to Jensen concavity effect, and transfers across regimes.
- Report both \(B\) and chosen \(\alpha\) for reproducibility; provide conversion via \(D = \alpha (B-1)(Ri_{\text{ref}}/Ri_b)\).

## Improved Mathematical Reconstruction and Critique of the McNider–Biazar D Formulation

### 1. Original handwritten form (as entered)
\[
f_c(\Delta z) = e^{D\left(\frac{\gamma}{Ri_c}\right) Ri \left(1 - \frac{\Delta z_r}{\Delta z}\right)}
\]
Issues:
- Sign: positive exponent increases mixing for larger Ri; for stability damping a negative sign is normally required.
- Units: \((\gamma/Ri_c)Ri\) is dimensionless only if \(\gamma\) itself is dimensionless (OK), but \(D\) must be dimensionless too; confirm.
- Grid factor: \(1 - \Delta z_r/\Delta z\) → 0 at fine resolution (good) but saturates (→1) at coarse resolution; lacks ζ (stability-depth coupling).
- Driver: uses Ri (unspecified point or bulk); if Ri_b is used, curvature bias is embedded and weakens response relative to Ri_g.

### 2. Corrected “linear tail-extension” form (keeping their structure)
Damping (physically consistent) version:
\[
f_c^{(lin)}(\Delta z) = \exp\left[-\,D \left(\frac{\gamma}{Ri_c}\right) Ri \left(1 - \frac{\Delta z_r}{\Delta z}\right)\right].
\]
Properties:
- \(f_c^{(lin)} \to 1\) as \(\Delta z \to \Delta z_r\).
- Monotone decrease with \(\Delta z\) for \(D,\gamma,Ri>0\).
Limitations:
- No neutral-preserving slope guarantee (\(\partial_\zeta f_c|_{\zeta=0} \neq 0\) unless Ri→0).
- No curvature awareness (Ri'' absent).
- No explicit bias correction; purely grid-scalar.

### 3. Bias / ODE-driven derivation (recommended)
Starting local invariance condition (grid independence of effective stability):
\[
\frac{d\ln f_c}{d\ln \Delta z} = - \alpha (B - 1)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q.
\]
Integrates to (power law):
\[
f_c^{(pl)}(\Delta z,\zeta) = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}.
\]
Exponential surrogate (smoother tail):
\[
f_c^{(exp)}(\Delta z,\zeta) = \exp\left[-\alpha (B-1)\left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^p\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q\right].
\]
Neutral preservation: choose \(q\ge 2\) so \(\partial_\zeta f_c|_0=0\). Grid convergence: \(f_c\to1\) as \(\Delta z\to0\).

### 4. Curvature-aware extension for D
Observed: D constant fails when curvature varies strongly with height. Define:
\[
D_{\text{eff}}(z) = D_0 + M \frac{|Ri''(z)|}{|Ri''(z)| + C},
\]
with:
- \(Ri''(z)\): second derivative (height or ζ-based, consistent scaling).
- \(C\): curvature soft-threshold (prevents runaway at small |Ri''|).
- Bound: \(D_{\min} \le D_{\text{eff}} \le D_{\max}\) (e.g. 0.2 ≤ D_eff ≤ 0.7).
Use \(Ri''\) from analytic MOST curvature if φ known; else numerical second difference with smoothing.

### 5. Reconciliation (verification)
Earlier mapping:
\[
\alpha (B-1) = D \frac{Ri_b}{Ri_{\text{ref}}} \Rightarrow D = \alpha (B-1)\frac{Ri_{\text{ref}}}{Ri_b}.
\]
Still valid if their ODE used \(Ri_b/Ri_{\text{ref}}\) as driver. If they reported D from bulk-Ri fits, recover \(\alpha\) via:
\[
\alpha = D \frac{Ri_b}{Ri_{\text{ref}}(B-1)}.
\]
Caution: use time/height-averaged Ri_b and Ri_g(z_g) over representative stable intervals; instant noisy values distort mapping.

### 6. Critique Summary
| Aspect | Handwritten form | Recommended |
|--------|------------------|-------------|
| Sign | Growth (mixing increase) | Damping (negative exponent) |
| Driver | Ri (ambiguous; likely bulk) | (B - 1) explicit bias |
| Neutral invariance | Not guaranteed | Enforced (q ≥ 2) |
| Curvature use | None | Optional D_eff(Ri'') |
| Grid/ζ scaling | Linear (1 - Δz_r/Δz) only | Power/exponential with ζ coupling |
| Calibration clarity | Weak (single D) | Multi-parameter (α, p, q, bounds) |

### 7. Path Forward (detailed steps)
1. Archive current formulations (keep original fc text for provenance).
2. Implement bias-based exponential form \(f_c^{(exp)}\) with initial parameters: \(\alpha=1.0, p=1, q=2, \Delta z_{\text{ref}}=10\text{ m}, \zeta_{\text{ref}}=0.5\).
3. Compute per-layer B = Ri_g(z_g)/Ri_b; apply fc only if B > B_thresh (e.g. 1.05).
4. Add curvature diagnostic: compute Ri'' analytically (if φ known) or numerically; derive D_eff; replace α by α·(D_eff/D_ref) if needed.
5. Bound fc: fc = max(fc, fc_min) with fc_min ≈ 0.2.
6. Validation loop:
   - Metrics: bias reduction (ΔB), neutral curvature error (<5%), temperature profile RMSE vs fine grid, inversion height bias.
   - Sensitivity: vary α, p, q, and curvature parameters (D_0,M,C).
7. Reconcile historical D values: produce table converting published D to α for key cases.
8. Reporting: log (Ri_b, Ri_g, B, Ri'', D_eff, fc) each timestep for evaluation.
9. Optional refinement: switch from Ri_b to Ri_g(z_g) in curvature zone to avoid bulk underestimation bias feeding back into D.
10. Transition plan: after stable performance, test dynamic Ri_c*; ensure fc and Ri_c* changes do not conflict (log both).

### 8. Practical Implementation Snippet (concise)
```python
# inputs: z0,z1,Ri_b,Ri_g_zg,zeta,alpha,p,q,dz,dz_ref,zeta_ref,B_thresh,fc_min
B = Ri_g_zg / Ri_b if Ri_b > 1e-6 else 1.0
if B <= B_thresh:
    fc = 1.0
else:
    # optional curvature weighting
    D_eff = D0 + M * (abs(Ri_ddz)/(abs(Ri_ddz)+C))  # Ri_ddz precomputed
    alpha_eff = alpha * (D_eff / D_ref)
    exponent = -alpha_eff * (B - 1) * (dz/dz_ref)**p * (zeta/zeta_ref)**q
    fc = max(math.exp(exponent), fc_min)
K_new = K_old * fc
```

### 9. Where We Can Help (targeted)
- Supply analytic Ri'' for chosen φ set to avoid noisy numerical second differences.
- Provide automated blending of bias-based fc with curvature-weighted D_eff.
- Build validation scripts: bias statistics, curvature preservation, profile RMSE.
- Convert legacy D values to α and populate a reconciliation table for manuscript appendix.
- Advise on dynamic Ri_c* integration (sequence: stabilize fc first, then Ri_c*).

### 10. Key Risks and Mitigations
| Risk | Mitigation |
|------|------------|
| Noisy Ri'' inflates D_eff | Apply small centered smoothing window; cap derivative |
| Over-damping (fc too small) | Enforce fc_min; monitor flux sensitivity |
| Misuse of bulk Ri for driver | Replace with (B - 1) bias; retain mapping for comparison only |
| Neutral distortion | Verify \(\partial_\zeta fc|_0 \approx 0\) numerically each run |

### 11. Final Summary (one paragraph)
The original D-based correction as entered lacks a damping sign, explicit curvature coupling, and a bias-normalized driver; adopting a bias-triggered, neutral-preserving exponential (or power-law) form with optional curvature-weighted D_eff provides a controlled, dimensionless, and verifiable path that reduces grid-induced stability underestimation while safeguarding near-neutral invariants (2Δ). Mapping between historical D and current α allows continuity; curvature diagnostics refine spatial targeting; a structured validation loop consolidates performance before dynamic Ri_c* extension.
