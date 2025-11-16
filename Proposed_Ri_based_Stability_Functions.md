# Candidate Ri-Based Stability Functions (fs)

## Catalog
1. Exponential–Rational: fs = exp( - a Ri / (1 + b Ri) )
2. Logistic–Exponential: fs = exp(-γ Ri) / (1 + (Ri/Ric)^p )
3. Rational: fs = 1 / (1 + c Ri + d Ri^2 )
4. Curvature-Aware Hybrid: fs = exp( - a Ri [1 - A] ) * exp( - a Ri/(1 + b Ri) A ), A=1-exp(-κ/κ0)
5. Grid-Corrected: fs_eff = fs * exp( D (γ/Ric) Ri (1 - Δzr/Δz) )
6. Dynamic Ric*: replace Ric by Ric* (in any above).

## Parameter Roles
- a,γ: neutral slope / initial suppression.
- b,p: tail length shaping.
- Ric,Ric*: transition scaling.
- c,d: polynomial damping curvature.
- D: grid-tail amplification.
- κ0: curvature sensitivity threshold.

## Physical / Numerical Checks
- fs(0)=1, df/dRi|0 = -a_target within tolerance.
- fs monotone, no plateaus unless intended for intermittent regime.
- fs_eff invariance: Δ(fs_eff)/Δz small.
- Km,Kh positive; fm/fh realistic Prandtl behavior.
- Transition wind V_R preserved across Δz.
- Energy budget closure (integral flux differences < threshold).

## Suggested Default Sets
- Research high-res: a=3.2, b=0.4, D=0.
- Coarse climate: a=3.0, b=0.8, D=0.4, Ric*=0.28.
- Strong inversion w/ curvature: κ0 tuned so A≈0.8 at surface curvature peak.

## Pros/Cons Snapshot
- Short-tail (exp-rational): physical, may under-mix coarse grids without fc.
- Long-tail (logistic mix): robust in coarse grids, risks overwarming.
- Rational: simple, can misfit neutral gradient if not tuned.
- Hybrid curvature: adaptive, needs curvature estimate quality.
- Grid-corrected: minimal code cost, depends on Δz and Ric accuracy.

## Minimal Code Skeleton (language-neutral)
```
// Given Ri, Δz, κ
fs_base = exp(-a*Ri/(1+b*Ri))
A = 1 - exp(-κ/κ0)
fs_curv = exp(-a*Ri*(1-A)) * exp(-a*Ri/(1+b*Ri)*A)
fc = exp(D*(γ/Ric)*Ri*(1 - Δzr/Δz))
fs_eff = fs_curv * fc
Pr = Pr0 + Pr1*Ri
fh = fs_eff
fm = fh / Pr
Km = fm * l*l*s
Kh = fh * l*l*s
```

## Adoption Path
1. Calibrate a,b,D on single-column benchmark (multi-Δz runs).
2. Introduce curvature gating only if Δ(fs)/Δz > tolerance.
3. Evaluate dynamic Ric* effect on inversion retention vs warming sensitivity.
4. Promote formulation passing invariance and transition tests.

## Transition Metric
Maintain |V_R(Δz_i) - V_R(Δz_ref)| < ε_V for regime shift (weak→strong stability).

## Tail Control Strategy
If overshoot in mixing at Ri>1: apply multiplicative damping exp(-λ (Ri-1)) for Ri>1 (λ≥0) preserving lower Ri slope.

