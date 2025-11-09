# Surface Roughness, Log Wind Profile, and Geometric Mean Heights

## 1. Neutral Log Profile
Under neutral stratification (ignore stability corrections):
\[
U(z) = \frac{u_*}{\kappa}\,\ln\!\left(\frac{z}{z_0}\right),
\]
where:
- \(u_*\): friction velocity
- \(\kappa \approx 0.40\): von Kármán constant
- \(z_0\): aerodynamic roughness length
- \(z\): measurement height (above displacement height if applicable)

## 2. Averaging Velocities Across Multiple Heights
Let heights \(z_1,\dots,z_N\) and corresponding mean horizontal speeds \(U_i = U(z_i)\). Then:
\[
\bar{U} = \frac{1}{N}\sum_{i=1}^N U(z_i) = \frac{u_*}{\kappa}\left[\frac{1}{N}\sum_{i=1}^N \ln\!\left(\frac{z_i}{z_0}\right)\right]
= \frac{u_*}{\kappa}\ln\!\left(\frac{( \prod_{i=1}^N z_i )^{1/N}}{z_0}\right).
\]
Define the geometric mean height:
\[
z_g = \Big(\prod_{i=1}^N z_i\Big)^{1/N}.
\]
Then:
\[
\bar{U} = U(z_g).
\]
Conclusion: The arithmetic mean of the velocities equals the log-profile velocity evaluated at the geometric mean of heights, not at the arithmetic mean height. Any interpretation of “representative” height for the averaged velocity should use \(z_g\).

## 3. Implication for Roughness and Friction Velocity Estimation
If one regress \(U_i\) against \(\ln z_i\):
\[
U_i = A \ln z_i + B,\quad A = \frac{u_*}{\kappa},\quad B = -\frac{u_*}{\kappa}\ln z_0.
\]
Thus:
- \(u_* = \kappa A\).
- \(z_0 = \exp\!\left(-B/A\right)\).

Averaging velocities first and pairing with an arithmetic mean height biases the inferred \(z_0\) (systematically misrepresents the intercept), especially when heights span orders of magnitude.

## 4. Recommended Procedure
1. Use individual height measurements (do not pre-average across heights).
2. Restrict to neutral periods (e.g., small |Obukhov length|; apply stability corrections otherwise: add \(\psi_m(z/L)\) term).
3. Perform linear regression of \(U\) vs \(\ln z\).
4. Extract \(u_*\) and \(z_0\) from slope and intercept.
5. If reporting a single “representative height” for an averaged velocity, specify \(z_g\), not \(\bar{z}\).

## 5. Practical Notes
- Displacement height \(d\): For tall vegetation, replace \(z\) with \(z-d\) everywhere; geometric mean uses \(z_i-d\).
- Quality control: Exclude low-level sensor data if within roughness sublayer (~first 2–3 \(z_0\)).
- Uncertainty: Propagate regression standard errors to \(u_*\) and \(z_0\); log-space linearity simplifies error estimates.

## 6. Summary
Neutral-layer logarithmic dependence means velocity averaging across levels maps to geometric mean height; roughness and friction velocity estimates should honor the log-linear structure, avoiding arithmetic-height substitutions.

## 7. Implications for Drag Coefficient (C_D)
Definition:
\[
C_D(z) \equiv \frac{\tau}{\tfrac{1}{2}\rho\,U(z)^2} = \left(\frac{u_*}{U(z)}\right)^2.
\]
Under neutral conditions (outside the roughness sublayer) with displacement height \(d\):
\[
U(z) = \frac{u_*}{\kappa}\ln\!\frac{z-d}{z_0}
\quad\Rightarrow\quad
C_D(z) = \left[\frac{\kappa}{\ln\!\frac{z-d}{z_0}}\right]^2.
\]
With stability correction \(\psi_m\):
\[
C_D(z,L) = \left[\frac{\kappa}{\ln\!\frac{z-d}{z_0}-\psi_m(z/L)}\right]^2,
\]
and \(\psi_m\equiv 0\) for neutral.

- Averaging across heights: since \(\bar{U}=U(z_g)\), the drag coefficient implied by the averaged speed is
  \[
  C_D(\bar{U}) = \left[\frac{\kappa}{\ln\!\frac{z_g-d}{z_0}}\right]^2
  \quad(\text{neutral}).
  \]
- Using the arithmetic mean height \( \bar{z} \) instead yields
  \[
  C_D(\bar{z})=\left[\frac{\kappa}{\ln\!\frac{\bar{z}-d}{z_0}}\right]^2
  \le C_D(z_g),
  \]
  because \( \bar{z}\ge z_g \Rightarrow \ln\!\frac{\bar{z}-d}{z_0}\ge \ln\!\frac{z_g-d}{z_0}\).
  Hence pairing \(\bar{U}\) with \(\bar{z}\) systematically underestimates \(C_D\).

Bias (neutral):
\[
\frac{C_D(\bar{z})}{C_D(z_g)}=\left[\frac{\ln\!\frac{z_g-d}{z_0}}{\ln\!\frac{\bar{z}-d}{z_0}}\right]^2 \le 1.
\]
For small spread with \(\delta=\ln\!\frac{\bar{z}}{z_g}\ll 1\) and \(L_g=\ln\!\frac{z_g-d}{z_0}\):
\[
\frac{C_D(\bar{z})}{C_D(z_g)} \approx \left(1-\frac{\delta}{L_g}\right)^2 \approx 1-\frac{2\,\delta}{L_g},
\]
so relative underestimation is approximately \( \approx 2\,\ln(\bar{z}/z_g)/\ln((z_g-d)/z_0)\).

Recommended workflow:
1. Estimate \(u_*\) and \(z_0\) by regressing \(U\) on \(\ln(z-d)\) (or include \(\psi_m\) when non-neutral).
2. Choose a reference height \(z_{\mathrm{ref}}\) and compute
   \[
   C_D(z_{\mathrm{ref}})=\left[\frac{\kappa}{\ln\!\frac{z_{\mathrm{ref}}-d}{z_0}-\psi_m(z_{\mathrm{ref}}/L)}\right]^2.
   \]
3. If reporting a single drag coefficient for averaged multi-level speeds, report it at \(z_g\) (or explicitly at \(z_{\mathrm{ref}}\)), not at \(\bar{z}\).
4. Avoid using data within the roughness sublayer (~first 2–3 \(z_0\)); ensure stationarity and neutral filtering when invoking the neutral form.

## 8. Stability Corrections (ψ) and Roughness for Heat
Stability-corrected integrated forms (with displacement d):
\[
U(z)=\frac{u_*}{\kappa}\Big[\ln\frac{z-d}{z_0}-\psi_m(z/L)\Big],\quad
\theta(z)-\theta_s=\frac{\theta_*}{\kappa}\Big[\ln\frac{z-d}{z_{0h}}-\psi_h(z/L)\Big].
\]
Canonical choices:
- Unstable (\(\zeta<0\)):
  \[
  \phi_m=(1-16\zeta)^{-1/4},\ \phi_h=(1-16\zeta)^{-1/2},\ x=(1-16\zeta)^{1/4},
  \]
  \[
  \psi_m=2\ln\frac{1+x}{2}+\ln\frac{1+x^2}{2}-2\arctan x+\frac{\pi}{2},\quad
  \psi_h=2\ln\frac{1+x^2}{2}.
  \]
- Stable (\(\zeta>0\)): \(\phi_{m,h}=1+5\zeta,\ \psi_{m,h}=-5\zeta\) (variants with mild curvature are used in very stable conditions).

Heat roughness \(z_{0h}\):
- Typically \(z_{0h}\le z_0\); over land often parameterized via roughness Reynolds or fixed \(z_{0h}/z_0\).
- Use separate \(\zeta_{0h}=z_{0h}/L\) in integrated temperature.

## 9. Bulk Transfer Coefficients with Stability
At reference height \(z_r\) (outside roughness sublayer):
\[
C_D(z_r)=\Bigg[\frac{\kappa}{\ln\frac{z_r-d}{z_0}-\psi_m(\zeta_r)}\Bigg]^2,\quad
C_H(z_r)=\frac{\kappa^2}{\big[\ln\frac{z_r-d}{z_0}-\psi_m(\zeta_r)\big]\big[\ln\frac{z_r-d}{z_{0h}}-\psi_h(\zeta_r)\big]}.
\]
Workflow: estimate \(u_*,z_0,z_{0h}\) by stability-aware regression; evaluate \(C_{D,H}\) at a specified \(z_r\) (or at \(z_g\) for averaged winds).

## 10. Link to Ri_g Curvature Diagnostics
Stability adjustments (non-neutral): add ψ_m, ψ_h. Their ζ-gradients feed \(v_m,v_h\) ⇒ \(V_{\log},W_{\log}\) for curvature:
\[
V_{\log}=v_h-2v_m,\quad W_{\log}=V_{\log}'.
\]
Including curvature at the first level reduces drag coefficient bias under stable stratification (see curvature toolkit).

## 11. Heat vs Momentum Roughness (z0h ≠ z0) — Why and When
Key points
- Temperature is a passive (or weakly active) scalar near the surface and must be continuous at the interface (soil/water/skin) and below it, while velocity goes to zero at the roughness elements. Consequently, the effective scalar (heat) roughness length z0h is often much smaller than the momentum roughness z0.
- Use kB−1 ≡ ln(z0/z0h) as the diagnostic. Typical ranges:
  - Smooth/water: kB−1 ≈ 10–30 (cool-skin/warm-layer physics; z0h ≪ z0).
  - Short crops/grass: kB−1 ≈ 2–7 (z0h < z0).
  - Shrub/forest/urban: kB−1 ≈ 5–15 (displacement height d critical; within-canopy sources break simple logs).
- Setting z0h = z0 is a neutral, smooth-surface convenience; it overestimates surface heat transfer when the molecular sublayer controls the scalar exchange (common in stable nights, calm wind, water).

Practical choice
- Prefer separate z0h (and z0q for moisture) when:
  - Stable/very stable stratification (ζ > 0, weak turbulence).
  - Over water/ice (cool-skin/warm-layer, surface renewal).
  - Sparse canopies/rough urban fabrics (different source heights for momentum vs scalars).
- z0h ≈ z0 is acceptable when:
  - Neutral to mildly unstable over aerodynamically rough, homogeneous terrain with strong turbulence and shallow molecular sublayers.

kB−1 parameterizations (for reference)
- Empirical relations exist as functions of roughness Reynolds number Re∗ and stability. A simple diagnostic is to tune kB−1 so that observed sensible heat flux and profile match via the bulk formula in Section 9.

## 12. Where the Temperature Log-Profile Holds, and Where It Breaks
Log-layer validity (MOST)
- Height range: Roughness sublayer top to 10–20% of boundary layer depth. Practically, (z − d)/z0 ≳ 30 and z/hmix ≲ 0.1.
- Stratification: MOST with φh(ζ) holds for stationary, horizontally homogeneous flow with continuous turbulence.

Breakdown regimes
- Molecular sublayer (very near surface): Scalar transfer is diffusion-limited; the log law for θ is not valid within the viscous/conductive sublayer. For heat, the scalar sublayer thickness δh is O(ν/u∗)Pr−2/3; if measurement z ≲ O(δh), use surface renewal or skin models (especially over water/ice).
- Very stable nights: Intermittent/laminar patches, suppressed turbulence; φh grows rapidly; log similarity degrades. Use capped/regularized φh or multi-layer resistance schemes.
- Within canopies/urban arrays: Displacement height d required; below/within canopy the profile is not logarithmic; use canopy-resistance models (source/sink distributions).
- Convective mixed layer aloft: The log segment is confined near the surface; above, nonlocal transport dominates (−w′θ′ terms not captured by local φh alone).

Continuity at the interface
- Temperature is continuous across the surface and into the substrate; the correct boundary condition is the surface (skin) temperature with conductive heat flux into soil/water. Thus, z0h is an effective parameter of the turbulent-convective layer, not a physical “zero-crossing” height like z0 for momentum.

## 13. Which “Temperature” to Use (θ, θv, q, others?)
- θ (potential temperature): use for scalar temperature profiles and sensible heat flux H over land; this is the standard in MOST for φh and ψh in the absence of strong moisture effects.
- θv (virtual potential temperature): use for buoyancy, Obukhov length L, Ri and N2 (stability metrics) because buoyancy depends on moisture; compute L from θv (and heat flux w′θ′v).
- q (specific humidity) or r (mixing ratio): treat as a separate scalar with its own roughness z0q (typically z0q ≈ z0h but not guaranteed). Latent heat flux uses q with φq analogous to φh.
- Sonic temperature Ts: diagnostic from sonic anemometers; relates to θ and q (Ts ≈ θ(1 + 0.51q) near surface); use with care in moist conditions.
- Over the ocean/ice: use skin SST/IST and cool-skin/warm-layer corrections; z0h and z0q often differ substantially from z0 due to surface renewal and molecular sublayers.

Rule of thumb
- Use θ for scalar profiles and H, θv for L, Ri, and buoyancy-related terms; treat humidity as a separate scalar with its own roughness.

## 14. Implications for Modeling and Retrievals
- Bulk transfer (Section 9) already supports z0h ≠ z0:
  \[
  C_H(z_r)=\frac{\kappa^2}{\big[\ln\frac{z_r-d}{z_0}-\psi_m(\zeta_r)\big]\big[\ln\frac{z_r-d}{z_{0h}}-\psi_h(\zeta_r)\big]}.
  \]
  If z0h is set equal to z0 in regimes with large kB−1, CH is biased low, underestimating H.
- Retrieval consistency: When fitting z0 and z0h, enforce that θ is continuous at the surface and use observed surface temperature (skin) as the lower boundary; avoid extrapolating θ to z = z0h as if it were a physical temperature at that height.
- Stability coupling: Since L uses θv, a mistaken choice of z0h can indirectly bias φh(ζ), L, and curvature diagnostics; calibrate kB−1 jointly with φh parameters (αh, βh) using flux and profile data.

## 15. Quick Guidance (Rules of Thumb)
- If in doubt, estimate kB−1 = ln(z0/z0h) from flux–profile pairs; expect kB−1 > 0 in most cases; water/ice often much larger than land.
- Keep z0h = z0 only for rough, neutral to mildly unstable flow with strong turbulence and shallow molecular sublayers.
- Use θ for temperature profiles and H; use θv for L, Ri, and curvature; treat humidity (q) with its own roughness z0q.
- Over canopies/urban, include displacement d and consider multi-layer resistance; within-canopy temperatures do not follow a log law.

## 16. Polar / Arctic Implications of Setting z0h = z0
Problem context
- Polar surfaces (sea ice, snow, thin melt ponds) often have molecular/renewal processes making z0h ≪ z0 (large kB−1). Forcing z0h = z0 removes this contrast.
Key impacts
- Sensible heat flux H overestimated in very stable cases if model compensates via tuned φh (or underestimated if using bulk transfer directly); either misrepresents surface cooling rate → inversion strength error.
- Obukhov length L error: θv flux term biased → ζ=z/L mis-scaled → shifts stability function curvature (Ri_g growth). Neutral-curvature invariant (2Δ) applied at wrong ζ → grid-dependent mixing biases.
- Ri_g curvature misclassification: inflated (or deflated) ∂²Ri_g/∂z² alters triggering of SBL turbulence shutdown; affects nocturnal and polar night inversion persistence.
- Feedback to Arctic Amplification: incorrect partition of turbulent vs radiative cooling modifies surface temperature trend; biases amplify through ice–albedo seasonal transitions (earlier melt onset or delayed freeze).
- Large-scale oscillations (AO/NAO): surface stress and heat flux anomalies modify lower-tropospheric temperature gradients, influencing baroclinicity; cumulative bias in stable-season fluxes can shift mean AO phase statistics.
- Vertical diffusion tuning compensation: models may lower critical Ri or lengthen φ tails to counter z0h=z0 simplification → introduces artificial grid dependence.
- Energy budget closure: latent vs sensible flux partition error (if z0q also set equal) leads to spurious moisture tendencies, altering cloud formation and longwave feedbacks.
Diagnostic signals
- kB−1 near zero in outputs over ice/snow.
- Systematic L too small (or large) compared with tower/ship observations.
- Curvature spectrum (∂²Ri_g/∂z²) flatter than observational retrievals in strong stability windows.
Recommended actions
- Implement distinct z0h (and z0q) parameterization tied to roughness Reynolds number or ice/snow surface state.
- Use curvature-aware correction (neutral invariance + tail modifier) with proper ζ scaling from corrected L.
- Validate kB−1 climatology and Ri_g curvature against Arctic sites (e.g. SHEBA, Barrow) before AO/NAO attribution studies.
- Include HS-assisted ζ↔Ri inversion to reduce numerical damping artifacts in stable polar layers.

Net implication
Using z0h = z0 in polar modeling distorts heat flux and stability scaling, propagating through Ri_g curvature to boundary-layer mixing, amplifying or muting Arctic Amplification signals and modulating large-scale oscillation statistics; separation of scalar and momentum roughness is therefore required for physically consistent SBL feedbacks.

