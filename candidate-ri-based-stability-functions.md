Candidate Ri-Based Stability Functions for ABL Mixing Models

1. Functional Families

1. Exponential–Rational

$$
f_s(Ri) \;=\; \exp\!\left[-\,\frac{a\,Ri}{1 + b\,Ri}\right]
$$
- Short-tail suppression with nonlinear tapering.
- Reduces to classic $\exp(-a\,Ri)$ for small $Ri$.

⸻

2. Logistic–Exponential

$$
f_s(Ri) \;=\; \frac{\exp(-\gamma\,Ri)}{1 + (Ri/Ri_c)^p}
$$
- Logistic denominator lengthens the tail.
- Good for coarse grids where long-tail behavior delays collapse of mixing.

⸻

3. Rational Polynomial

$$
f_s(Ri) \;=\; \frac{1}{1 + c\,Ri + d\,Ri^2}
$$
- Simple, explicit curvature control via $(c,d)$.
- Often underestimates neutral slope unless tuned.

⸻

4. Curvature-Aware Hybrid

$$
A \;=\; 1 - \exp\!\left(-\frac{\kappa}{\kappa_0}\right),\qquad
f_s(Ri) \;=\; \exp\!\big[-a\,Ri(1-A)\big]\;\exp\!\left[-\,\frac{a\,Ri}{1+b\,Ri}\,A\right]
$$
- Interpolates between exponential and exp–rational depending on curvature amplitude.
- $A$ acts as a curvature gate:
  - $A\to0$: exponential regime
  - $A\to1$: rational–exponential regime

⸻

5. Grid-Corrected ($\Delta z$-Aware)

$$
f_{s,\mathrm{eff}}(Ri,\Delta z) \;=\; f_s(Ri)\;
\exp\!\left[D \left(\frac{\gamma}{Ri_c}\right) Ri \left(1 - \frac{\Delta z_r}{\Delta z}\right)\right]
$$
- Compensates for coarse-grid under-suppression.
- Vanishes ($f_c\to1$) when grid spacing approaches the reference resolution.

⸻

6. Dynamic Critical Richardson Number

Apply to any of the above:
$$
Ri_c \;\rightarrow\; Ri_c^*(z,t)
$$
- Enables flow-regime feedback (e.g., inversion depth, shear intermittency).
- Useful for maintaining inversion strength in warming-sensitive climate runs.

⸻

2. Parameter Roles

Parameter — Role  
$a,\gamma$: neutral-limit slope; initial suppression rate.  
$b,p$: tail-length / curvature shaping.  
$Ri_c, Ri_c^*$: transition scale from weak to strong stability.  
$c,d$: explicit polynomial curvature control.  
$D$: grid-coarsening correction amplitude.  
$\kappa_0$: curvature sensitivity threshold for hybrid formulation.

⸻

3. Physical & Numerical Consistency Checks

Boundary and Neutral Behavior
- $f_s(0)=1$
- $f_s'(0) = -a_{\text{target}}$ within tolerance (neutral calibration).

Monotonicity & Tail Behavior
- $df_s/dRi < 0$, except controlled plateaus for intermittent mixing regimes.
- No artificial shoulder near $Ri \approx Ri_c$.

Grid-Resolution Invariance
- $\Delta f_{s,\mathrm{eff}} / \Delta z$ small across typical $\Delta z$ changes.

Diffusivity & Prandtl Constraints
- $K_m,K_h > 0$
- $f_m/f_h = 1/Pr$ under chosen parameterization.
- Avoid unphysical Prandtl spikes at high $Ri$.

Transition Wind Invariance
- Maintain $\big|V_R(\Delta z_i)-V_R(\Delta z_{\rm ref})\big| < \varepsilon_V$ to ensure consistent regime transition.

Energy Budget Closure
- Vertical integral of flux divergences remains within target tolerance.

⸻

4. Default Parameter Sets

High-Resolution Research Runs
- $a = 3.2$, $b = 0.4$, $D = 0$

Coarse Climate Models
- $a = 3.0$, $b = 0.8$, $D = 0.4$, $Ri_c^* = 0.28$

Strong Surface Inversions
- Choose $\kappa_0$ so that $A \approx 0.8$ at surface curvature peak.

⸻

5. Tail-Control Strategy

To reduce excessive mixing for $Ri>1$:
$$
f_s \;\leftarrow\; f_s \,\exp\!\big[-\lambda\,(Ri-1)\big],\qquad \lambda \ge 0
$$
Shape-preserving for low $Ri$, damping only in the high-$Ri$ tail.

⸻

6. Minimal Implementation Skeleton

```text
// Inputs: Ri, Δz, κ
fs_base  = exp(-a * Ri / (1 + b * Ri))
A        = 1 - exp(-κ/κ0)
fs_curv  = exp(-a*Ri*(1-A)) * exp(-a*Ri/(1+b*Ri)*A)
fc       = exp(D*(γ/Ric)*Ri*(1 - Δzr/Δz))
fs_eff   = fs_curv * fc

Pr       = Pr0 + Pr1*Ri
fh       = fs_eff
fm       = fh / Pr

Km = fm * l*l * s
Kh = fh * l*l * s
```

⸻

7. Recommended Adoption Path
1) Calibrate $(a,b,D)$ using single-column benchmarks at multiple $\Delta z$.  
2) Enable curvature gating (hybrid or $\kappa$-based) only when $\Delta f_s / \Delta z$ exceeds tolerance.  
3) Test dynamic $Ri_c^*$ for inversion preservation vs sensitivity.  
4) Adopt formulation that passes both invariance and transition-wind criteria.