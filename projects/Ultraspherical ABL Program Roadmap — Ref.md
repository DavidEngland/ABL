Ultraspherical ABL Program Roadmap — Refined, Execution‑Ready

Program Objective

Develop a next‑generation boundary‑layer closure using ultraspherical spectral structure, tracer‑specific transport scaling, and fractional‑dimension diagnostics to reduce stable‑regime errors, explain Arctic amplification pathways, and enable cross‑planet transferability.

---

Strategic Prioritization

Highest expected gain per unit effort:

1. Arctic stable ABL pilot with existing tower/flux data
(largest current-model failure + your strongest analytic leverage)
2. Tracer‑aware extension
(humidity → CO₂ → aerosols)
3. Planetary transfer experiments
(after Earth validation)


Rationale: Stable-regime misrepresentation drives both forecast error and Arctic amplification bias; ultraspherical structure directly targets the modal redistribution responsible for these failures.

---

Three Project Tracks

Track A — Core Science (Immediate)

A1. Define basis and mapping

• Select \(\xi(\zeta)\) transform (likely monotone, bounded, neutral‑limit preserving).
• Fix truncation order \(N\) with low‑order bias control.
• Lock a baseline MOST comparator (BD74/Högström + one stable branch).


A2. Estimate modal coefficients

• Fit \(a_n\) by stability class (ζ‑bins or Ri‑bins).
• Use bootstrap or Bayesian posterior intervals for uncertainty.
• Add cross‑site validation to prevent regional overfitting.


A3. Validate against physical outcomes

• Flux–gradient residuals.
• Stable collapse/recovery timing.
• Forecast skill during Arctic events (e.g., radiative clear nights, cold pools, LLJ suppression).


---

Track B — Tracer & Fractional Dimension (Next)

B1. Effective dimension estimator

Implement three independent estimators:

• Structure‑function scaling.
• Spectral‑slope method.
• Attractor/correlation dimension (optional but high insight).


B2. Tracer‑specific closures

• Humidity first (largest data availability + strong coupling to stability).
• CO₂ next (slower response, useful for modal persistence).
• Aerosols/reactive tracers when data allow.


B3. Coupled evaluation

Test whether tracer‑dependent effective dimension reduces:

• Flux‑gradient bias,
• Intermittency misclassification,
• Transition timing errors.


---

Track C — Planetary Transfer (Later, High Upside)

C1. Mars benchmark

• Low density, strong radiative forcing, shallow SBL.


C2. Titan benchmark

• High density, methane humidity, slow dynamics.


C3. Cross‑world nondimensional collapse

• Compare ultraspherical modal redistribution across atmospheres.


---

90‑Day Execution Plan — Refined

Days 1–30 — Foundation

1. Build baseline notebook + reproducible data pipeline.
2. Reproduce classic MOST metrics for one Arctic region.
3. Fit ultraspherical model for momentum & heat.


Exit criteria:

• Fully reproducible run.
• One publication‑quality figure (modal structure vs MOST).
• Documented uncertainty for \(a_n\).


---

Days 31–60 — Tracer & Dimension

1. Add humidity tracer + effective dimension metric.
2. Run ablation studies (with/without dimension term).
3. Identify event classes with largest improvement.


Exit criteria:

• Statistically significant improvement in ≥1 key metric.
• Clear diagnostic linking modal redistribution to humidity‑dependent dimension.


---

Days 61–90 — Synthesis & Externalization

1. Draft manuscript‑style methods/results section.
2. Build external AI‑ready notebook baseline.
3. Prepare proposal for expanded station network/collaboration.


Exit criteria:

• Submission‑ready preprint skeleton.
• External notebook package for collaborators.


---

Risk Controls (Strengthened)

1. Overfitting• Enforce cross‑site validation.
• Start with low \(N\); increase only when modal curvature justifies it.

2. Data sparsity• Regional assimilation of coefficients before increasing model order.
• Use hierarchical priors for stability‑class pooling.

3. Interpretability• Every added parameter must map to a physical mechanism (e.g., modal steepening ↔ enhanced shear production).



---

Deliverables for Collaborators

1. Concept graphic
Earth shell + ABL shell + ultraspherical modal decomposition.
2. One‑page overview
“Why ultraspherical now” — emphasize stability, curvature, and modal redistribution.
3. Baseline notebook kit
References, starter tests, reproducible pipeline.


---

Collaboration Targets

1. Arctic observational teams (tower/flux networks).
2. Planetary boundary‑layer modelers.
3. Data assimilation & inverse‑problem specialists.


---

Success Metrics

1. Reduced stable‑regime flux‑gradient error.
2. Improved representation of transition/intermittency.
3. Mechanistic link between modal redistribution & Arctic amplification.
4. Transferability to at least one non‑Earth atmosphere.


---

One actionable next step

Would you like me to turn this into a manuscript‑ready “Project Overview” section, a proposal‑ready two‑page summary, or a technical architecture for the notebook pipeline?