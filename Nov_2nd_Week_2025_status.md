# Weekly Status — Week of Nov 10, 2025
Recipients: R. T. McNider, A. P. Biazar
Owner: D. E. England
Repo: https://github.com/DavidEngland/ABL

## 1) Executive Summary (1–2 lines)
- Curvature framework stabilized: compact d²Ri_g/dζ² expression, neutral invariant (2Δ) preserved, variable-L mapping added. Draft figures and code stubs ready; need decisions on stable φ form and calibration dataset to proceed to paper 1A draft.

## 2) Highlights (technical/milestones)
- Theory
  - Finalized compact curvature: d²Ri_g/dζ² = F[2V_log + ζ(V_log² − W_log)] with neutral limit 2Δ.
  - Derived ζ(Ri) near-neutral inversion + Newton refinement (1–2 iters to machine precision).
  - Variable L(z) mapping and omission metric E_omit implemented in notes.
  - Q‑SBL surrogate (pole‑free) matched to {Δ,c1} for ζ ≤ 0.2–0.3.
- Implementation
  - Profiles/curvature snippets added; φ-pack scaffolding and Ri-based closures outlined.
  - Geometric-mean height justification and diagnostics (bias B = Ri_g(z_g)/Ri_b).
- Publications/Docs
  - Curvature notes consolidated (curvature.md; Richardson number curvature.md).
  - SBL corrections playbook drafted (SBL corrections.md).
  - Publication roadmap and roles (topics.md) completed.
  - Urban/remote-sensing case skeleton (High-Resolution Ri_g in Megacity Boundary Layer.md).

## 3) Decisions Needed (by next check-in)
- D1: Stable φ baseline for ζ>0
  - Option A: Linear-stable (Högström/BH-style).
  - Option B: Q‑SBL surrogate (preferred for coarse Δz safety).
  - Option C: BH91 hybrid (slightly more complex).
- D2: Calibration dataset priority
  - Pick two: ARM NSA (Alaska), SHEBA (ice), GABLS1 LES, Urban (325 m tower).
- D3: Grid-damping factor G(ζ,Δz) defaults
  - Confirm initial (D=1.0, p=1.5, q=2, ζ_r=0.3, Δz_r=10 m) for first tests.
- D4: Target journal for Paper 1A
  - Boundary-Layer Meteorology (primary) vs JAS.

## 4) Next 2 Weeks (deliverables)
- Figs: curvature vs ζ (baseline vs Q‑SBL), bias ratio B vs Δz (10–100 m), constant‑L vs variable‑L (E_omit bands).
- Code: minimal Python module (curvature, ζ↔Ri, Q‑SBL, G factor) + notebook reproducing figures.
- Data: one site preprocessed profile set (ARM NSA or GABLS1).
- Draft: Paper 1A outline with equations, methods, and fig placeholders.

## 5) Risks / Blockers
- Data access/QA (tower/fluxes) may delay calibration; mitigation: start with GABLS1 LES.
- Parameter transfer from unstable fits can overstate stable curvature; mitigation: use SBL-only fits or Q‑SBL.
- Variable L(z) significance unknown per case; use E_omit < 0.05 rule to shortcut else map fully.

## 6) Metrics to Report
- Neutral curvature preservation: |(2Δ* − 2Δ)/(2Δ)| < 5%.
- Curvature error reduction at Δz=100 m: ≥ 40% vs baseline.
- ζ(Ri) inversion residual: < 1e−10 after Newton (2 steps).
- Fraction of layers with E_omit > 0.05.

## 7) Links (current)
- Theory: curvature.md; Richardson number curvature.md
- Implementation notes: SBL corrections.md; curvature code snippets inside docs
- Strategy: topics.md
- Urban case: High-Resolution Ri_g in Megacity Boundary Layer.md

## 8) Proposed Meeting Agenda (20–30 min)
1) Confirm φ baseline (ζ>0) and dataset sequence.
2) Approve G(ζ,Δz) defaults and diagnostics.
3) Paper 1A scope/timeline and author roles.
4) Student involvement (if applicable).

# Pros and Cons — McNider Similarity-Based Ri Curvature (Quick Brief)

Context
- Idea: use MOST stability functions near the surface (φ_m, φ_h; e.g., linear or exponential approximations) to compute Ri_g(ζ) = ζ φ_h/φ_m² and its curvature d²Ri_g/dζ² analytically; map to height via 1/L² (constant L) or full chain rule for L(z).

Pros
- Physics-grounded: directly ties curvature to φ_m, φ_h; preserves neutral invariant 2Δ if φ near-neutral slopes are respected.
- Closed-form and cheap: no ζ-iteration; easy to evaluate at representative height (prefer z_g).
- Drop-in: compatible with existing MOST codes; can be used to set curvature-aware lower-boundary conditions.
- Transparent diagnostics: 2Δ, c1, V_log, W_log available for QC and tuning.

Cons
- Sensitivity to φ choice: linear vs exponential vs BH/SHEBA fits change Δ and c1; misfit propagates to curvature.
- Domain/guards: power-law poles (ζ<1/β) or exponential tails need guards; variable L(z) can bias height mapping if ignored (use E_omit).
- Near-surface grid bias remains if evaluated at arithmetic means; must use geometric/log means for layers.
- Limited aloft: similarity-only curvature may miss secondary inflections or nonlocal effects above the surface layer.
- Pr_t dependence: assuming Pr_t≈1 simplifies algebra but can bias heat curvature in strong SBL unless calibrated.

Recommendations
- Use similarity-based curvature at z_g for the first layer; preserve 2Δ; apply constant-L map only if E_omit<0.05.
- Prefer Q‑SBL (quadratic) surrogate for ζ≤0.2–0.3 in very stable cases to avoid pole artifacts.
- If residual bias at coarse Δz, add neutral-preserving grid damping G(ζ,Δz); keep G(0)=1 and ∂G/∂ζ|0=0.

Decision asks
- Confirm φ baseline (linear/Q‑SBL/BH) for experiments.
- Approve z_g usage and E_omit threshold (0.05) for height mapping choice.

