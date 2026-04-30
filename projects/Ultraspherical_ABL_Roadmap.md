# Ultraspherical ABL Program Roadmap

## Program objective

Develop a next-generation boundary-layer closure framework using ultraspherical spectral structure, tracer-specific transport scaling, and fractional-dimension diagnostics to improve stable ABL prediction and explain Arctic amplification pathways.

## Decision rule: which direction yields most fruit

Prioritize by expected gain per unit effort:

1. Arctic stable ABL pilot with existing station data (highest).
2. Tracer-aware extension (humidity first, then CO2).
3. Planetary transfer experiments (Mars/Titan) after Earth validation.

Reason: Arctic stable-regime error reduction is where current closures fail most and where your existing work already provides leverage.

## Three project tracks

## Track A: Core science (now)

1. Define basis and mapping
- Choose $\xi(\zeta)$ transform and truncation order $N$.
- Lock baseline MOST comparator.

2. Estimate modes from observations
- Fit $a_n$ coefficients by stability class.
- Add uncertainty via bootstrap/Bayesian intervals.

3. Validate against outcomes
- Flux-gradient residuals,
- transition timing for stable collapse/recovery,
- forecast skill during Arctic events.

## Track B: Tracer and fractional dimension (next)

1. Effective dimension estimator package
- Structure functions,
- spectral slope method,
- attractor/correlation dimension option.

2. Tracer-specific closures
- humidity first,
- then CO2,
- then aerosols/reactive tracers as data allow.

3. Coupled evaluation
- Test whether tracer-dependent dimension lowers bias vs shared closure.

## Track C: Planetary transfer (later, high upside)

1. Mars benchmark case.
2. Titan benchmark case.
3. Compare nondimensional parameter collapse across worlds.

## 90-day execution plan

## Days 1-30

1. Build baseline notebook and data pipeline.
2. Reproduce classic MOST metrics for one Arctic region.
3. Fit ultraspherical model for momentum and heat.

Exit criteria:
- Reproducible run,
- one publishable-quality figure,
- documented uncertainty.

## Days 31-60

1. Add humidity tracer and effective dimension metric.
2. Perform ablation studies (with/without dimension term).
3. Identify event classes where improvement is largest.

Exit criteria:
- Statistically significant improvement in at least one key metric.

## Days 61-90

1. Draft manuscript-style methods/results section.
2. Build external AI notebook baseline reference set.
3. Prepare proposal for expanded station network or collaboration.

Exit criteria:
- Submission-ready preprint skeleton and external notebook package.

## Risk controls

1. Overfitting risk: enforce cross-site validation and low-order truncation first.
2. Data sparsity risk: assimilate coefficients regionally before increasing model order.
3. Interpretability risk: require each added parameter to map to a physical mechanism.

## Deliverables for collaborators

1. Concept graphic: Earth shell plus ABL shell modal decomposition.
2. One-page overview: why ultraspherical now.
3. Baseline notebook kit with references and starter tests.

## Collaboration targets

1. Arctic observational teams (tower and flux networks).
2. Planetary boundary-layer modelers.
3. Data assimilation and inverse-problem specialists.

## Success metrics

1. Error reduction in stable ABL flux-gradient prediction.
2. Better representation of transition/intermittency events.
3. Clear mechanistic link between modal redistribution and Arctic amplification diagnostics.
4. Transferability score to at least one non-Earth atmosphere.