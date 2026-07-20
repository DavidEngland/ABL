# Paper-Series Roadmap (Not Paper 1 Draft)

Status: program-level planning document
Last updated: 2026-07-01

## 1. Organizing Narrative for Paper 1

To keep the Gegenbauer basis as the organizing principle, the Paper 1 narrative sequence is fixed as:

1. Observational Evidence: establish scale-dependent SBL spectral roughness using high-resolution eddy-covariance data from SHEBA and SMEAR.
2. Methodological Innovation: introduce adaptive estimation of the Gegenbauer parameter $\lambda$ and show how tunable $\lambda$ absorbs anomalous scaling without ad hoc adjustments.
3. Spectral Representation: demonstrate improved capture of SBL spectral structure versus fixed-basis alternatives.
4. Geometric Interpretation: frame regime transitions as a consequence of state-space attractor geometry.
5. Validation and Limitations: validate on out-of-sample intervals and define the model's operational envelope.

---

## 2. Table 2: Conservative State-Space EWS and Falsification

The Table 2 framework is intentionally standalone and publishable as a conservative baseline testing whether classical indicators distinguish true fold approach from background noise.

### Table 2: Early Warning Signal (EWS) Baseline and Falsification Criteria

| Metric | Expected Behavior Approaching Fold | Failure Mode | Falsification Criterion |
| --- | --- | --- | --- |
| Variance | As the dominant stable eigenvalue approaches zero near the folded slow manifold, perturbations decay more slowly, increasing state variance under stochastic forcing. | Intermittent Bursts: high-amplitude localized turbulent excursions mimic variance inflation without regime shift. | Variance increase is not accompanied by corresponding deceleration of the mean state trajectory toward the fold line. |
| Lag-1 Autocorrelation | Near-zero eigenvalue reduces memory decay, causing lag-1 autocorrelation to approach 1. | Stationary Persistence: naturally long-lived stable turbulent structures skew short-window estimates. | Autocorrelation remains high without geometric approach to the attractor fold boundary in state space. |
| Recovery Rate | Return-to-equilibrium rate decreases (critical slowing down). | Measurement Gaps: downtime or filtering artificially flattens the relaxation curve. | Recovery-rate confidence intervals overlap substantially with baseline stationary-regime intervals. |
| Composite Indicator | Integrated metric (variance + lag-1 + recovery) flags high-probability transitions. | False Positives: coincident but decoupled sensor-noise spikes trigger alarms. | Discrimination falls below the ROC/AUC threshold set on the training set. |

---

## 3. Grand Arc: Papers 1 Through 4

- Paper 1 (Current): Adaptive Gegenbauer spectral closure and observational geometry. Establishes adaptive basis, SHEBA/SMEAR validation, and classical state-space EWS diagnostics.
- Paper 2 (Theoretical Progression): Dynamical systems diagnostics of folded slow-manifold transitions. Tests whether finite-time geometric diagnostics detect loss of normal hyperbolicity earlier than classical EWS.
- Paper 3 (Computational Realization): Computational algorithms and synthetic verification. Introduces reduced-state FTLE, finite-time growth rates, and synthetic-system benchmarking.
- Paper 4 (Physical Scaling): Extension to physical-space transport and LES applications. Connects reduced-order state-space geometry to physical transport in $\mathbb{R}^3$ using Haller LCS on LES fields.

---

## 4. Immediate Next-Writes (Paper 1)

Choose one and execute first:

1. Mathematical Framework subsection: stochastic forcing interacting with decaying eigenvalues near the fold.
2. Introduction/Narrative hook: SHEBA/SMEAR thread and motivation framing.

Recommended order: (1) then (2), so introduction claims are anchored by finalized equations.
