Title: Curvature‑Aware, Grid‑Dependent Corrections to MOST: Preserving Neutral Invariance and Reducing Coarse‑Grid Bias in Stable Boundary Layers

Authors:
- David E. England (draft)
- [Coauthor 2]
- [Coauthor 3]

Affiliations:
- (1) Department / Institution
- (2) ...

Date: YYYY‑MM‑DD

Abstract
--------
We present a concise, implementable formulation for grid‑aware corrections to Monin–Obukhov similarity functions that preserves the neutral curvature invariant (2Δ) while reducing coarse‑grid bias in estimates of the gradient Richardson number Ri_g. The correction is a multiplicative damping factor G(ζ,Δz) applied to eddy diffusivities (or directly to φ functions). We summarize analytic motivation, provide operational constraints on G (G(0)=1, ∂G/∂ζ|₀=0, monotonic tail damping), describe validation with LES/tower cases, and outline a minimal ML surrogate workflow for fast deployment. Draft is for internal iteration; figures and full numerical tables pending.

Keywords: Monin–Obukhov, Richardson number, curvature, stable boundary layer, grid correction, ML surrogate

Table of contents
-----------------
- 1 Introduction
- 2 Theory: Ri_g curvature and neutral invariant
- 3 Grid bias from layer averaging
- 4 Correction design and constraints
- 5 Implementation options
  - 5.1 Multiplicative G on K
  - 5.2 φ‑modifier alternative
  - 5.3 Q‑SBL surrogate
- 6 Validation & results (summary)
- 7 ML surrogate recipe (optional)
- 8 Discussion
- 9 Future work and outside help
- 10 Acknowledgements
- 11 Data & code availability
- 12 References / Bibtex

1 Introduction
--------------
- Motivation: coarse Δz causes systematic underestimation of near‑surface stability in the SBL → overmixing, warm biases.
- Short literature context: MOST, Businger/Beljaars, Jensen effects from concave Ri_g.
- Aim: provide a neutral‑preserving correction G(ζ,Δz), practical calibration & deployment notes.

2 Theory: Ri_g curvature and neutral invariant
----------------------------------------------
- Define ζ = z/L; Ri_g(ζ)=ζ φ_h/φ_m^2.
- Introduce logarithmic sensitivities V_log, W_log; present compact curvature formula:
  d²Ri_g/dζ² = F[2V_log + ζ(V_log² − W_log)]  (refer derivation in repo).
- Define neutral curvature invariant Δ and note d²Ri_g/dζ²|₀ = 2Δ.
- Physical interpretation: Δ<0 ⇒ concave‑down ⇒ Jensen bias (Ri_b < Ri_g(z_g)).

3 Grid bias from layer averaging
-------------------------------
- Explain geometric mean z_g = √(z0 z1) as representative height.
- Define bulk Ri_b and bias ratio B = Ri_g(z_g) / Ri_b.
- Summarize observed ranges (from experiments): typical B for coarse Δz (e.g., 50–100 m).
- Diagnostic procedure (short pseudocode or reference to notebook).

4 Correction design and constraints
----------------------------------
- Design principles:
  - Preserve neutral curvature 2Δ (no change at ζ→0).
  - Dampen tail curvature for ζ>0 proportionally to Δz (grid‑convergent).
  - Monotonic, pole‑free, computationally cheap.
- Proposed template:
  G(ζ,Δz) = exp[ − D (Δz/Δz_ref)^p (ζ/ζ_ref)^q ], with p≥1, q≥2, calibrated D.
- Operational constraints and guardrails (clipping, fallback to background K when laminar classifier triggers).

5 Implementation options
-----------------------
5.1 Multiplicative G on K
- K* = K · G applied in existing diffusion update; minimal code changes.
- Preserve K behavior at ζ→0 (G→1) and Δz→0.

5.2 φ‑modifier alternative
- φ*_m = φ_m · f_c(ζ,Δz) etc.; pros/cons (unified but needs chain‑rule care).

5.3 Quadratic SBL surrogate (Q‑SBL)
- Use quadratic φ for ζ∈[0,ζ_max]; blends to standard forms aloft; preserves neutral coefficients.

6 Validation & results (summary)
--------------------------------
- Short results summary (populate numeric values from experiments):
  - LES/GABLS1: bias reduction X% (Δz=100 m → B reduced from A→B).
  - Tower cases (ARM/SHEBA): RMSE of surface T reduced by Y K when using G.
  - Computational cost: ~Z% extra per timestep.
- Table and figure placeholders below.

Figure placeholders and captions
-------------------------------
- Fig 1 (placeholder): Ri_g(ζ) curves for selected φ sets and Δ values. Caption: "Analytic Ri_g (ζ) demonstrating neutral curvature 2Δ and sample concave‑down SBL profiles."
- Fig 2 (placeholder): Bias ratio B vs Δz (before/after correction). Caption: "Bias amplification B as a function of layer thickness Δz; G reduces B across tested range."
- Fig 3 (placeholder): Time series from tower/LES showing surface Ts with and without correction. Caption: "Surface temperature evolution: baseline vs curvature correction."
- Fig 4 (placeholder): ML surrogate parity plot (G_true vs G_pred) and AUC for laminar classifier. Caption: "ML surrogate performance metrics."

7 ML surrogate recipe (optional)
--------------------------------
- Use case: fast inference of Ĝ(ζ,Δz,Ri_b,features).
- Data: synthetic analytic sweep + LES/tower labeled targets.
- Models: small LightGBM / ONNX export or tiny MLP; enforce neutral constraint via loss augmentation or data augmentation at ζ≈0 (target G=1).
- Classifier: P_laminar with high‑precision threshold (e.g., 0.999) to gate K→K_background.
- Minimal training snippet and evaluation metrics: RMSE, MAE, AUC‑ROC, AUC‑PR, F1.

8 Discussion
-----------
- Strengths: physics‑anchored, simple deployability, robust in tested regimes.
- Limitations: dependence on φ choice and calibration; inflection handling when curvature sign changes within a layer; variable L(z) effects (use omission metric E_omit).
- Integration notes: patch points for WRF/other models; unit test suggestions.

9 Future work and outside help
------------------------------
- Immediate next steps: populate figures, full numerical tables, and code snippets for model insertion.
- External expertise recommended:
  - Field data / LES: request raw LES cases (GABLS, SHEBA) — contact LES maintainers or ARM data team.
  - Co‑author review: invite R.T. McNider and A. Biazar for model‑integration and evaluation.
  - ML calibration: small collaboration with ML practitioner (experience with LightGBM/ONNX) for surrogate export.
- Where to seek help:
  - ARM data access and QA — ARM Data Center.
  - GABLS LES runs — the GABLS community (contact leads listed in GABLS papers).
  - Model integration testing — colleagues managing WRF/CMAQ instances at our institution.

10 Acknowledgements
-------------------
- Funding placeholder.
- Data providers (ARM, LES authors).
- Colleagues for early feedback.

11 Data & code availability
--------------------------
- Link to repository: GitHub: /Users/davidengland/Documents/GitHub/ABL (add release/DOI when ready).
- Notebook references: notebooks/Curvature_Demo.ipynb and toy_sc_m results.

12 References / Bibtex
----------------------
- Include BibTeX entries; minimal skeleton:
```bibtex
@article{businger1971,
  author = {Businger, J. A. and et al.},
  title = {Flux-Profile Relationships in the Atmospheric Surface Layer},
  journal = {J. Atmos. Sci.},
  year = {1971},
  volume = {28},
  pages = {181--189}
}
@article{beljaars1991,
  author = {Beljaars, A. C. M. and Holtslag, A. A. M.},
  title = {Flux-Profile Relationships},
  journal = {Quart. J. Roy. Meteorol. Soc.},
  year = {1991}
}
```
- Extend with England & McNider (1995), Högström (1988), SHEBA/GABLS citations.

Notes for authors / TODOs
------------------------
- Replace placeholders with experiment figures and exact numbers.
- Decide author list and order; add affiliations and funding text.
- Add figure files to visuals/ and cite them in draft using relative paths.
- After internal review, convert to LaTeX using the existing TeX templates in /drafts when ready.

Outside help (summary)
----------------------
- LES & observational data access: ARM Data Center, GABLS community.
- Co‑author technical review: R. T. McNider (physics/interpretation), A. Biazar (operational validation).
- ML surrogate packaging: LightGBM/ONNX expert for low‑overhead runtime export.

