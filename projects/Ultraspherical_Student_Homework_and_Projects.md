# Ultraspherical ABL: Homework and Project Track for STEM Students

## Audience

Upper-level undergraduate, early graduate, and interdisciplinary STEM students (atmospheric science, applied math, data science, physics).

## Learning goals

1. Connect MOST physics with orthogonal polynomial methods.
2. Implement data-driven spectral closures from real meteorological data.
3. Interpret stability-regime transitions using modal diagnostics.
4. Communicate results through reproducible notebooks and visualizations.

## Homework sequence

## HW1: Coordinate transforms and inversion

Tasks:

1. Derive $\xi=\tanh(\alpha_{\xi}\zeta)$ and $\zeta=\alpha_{\xi}^{-1}\operatorname{artanh}(\xi)$.
2. Express artanh in logarithmic form and discuss numerical conditioning.
3. Compare tanh vs rational mapping for stability compression.

Deliverable:
- short derivation note and one figure of mapping behavior.

## HW2: Baseline MOST fit at one station

Tasks:

1. Fit $(a_q,b_q,\lambda_q)$ for momentum and heat transfer functions.
2. Quantify fit uncertainty and regime-dependent residuals.
3. Test whether fixed exponents are supported by data.

Deliverable:
- parameter table and residual diagnostics.

## HW3: Ultraspherical augmentation

Tasks:

1. Expand fitted function or residual in Gegenbauer basis.
2. Select truncation order by cross-validation.
3. Compare skill vs baseline MOST.

Deliverable:
- notebook with side-by-side model comparison.

## HW4: Fractional-dimension diagnostic

Tasks:

1. Estimate one effective-dimension metric from turbulence/tracer data.
2. Regress dimension vs stability class and shear.
3. Evaluate whether adding dimension improves closure skill.

Deliverable:
- short report and reproducible code cell outputs.

## Project sequence

## Project A: Arctic pilot

1. Use one month from one Arctic site.
2. Build baseline-plus-ultraspherical model stack.
3. Produce one publication-quality diagnostic figure.

## Project B: Multi-station regional cluster

1. Fit site-level modal coefficients.
2. Interpolate/assimilate coefficients across region.
3. Visualize evolving mode fields during selected events.

## Project C: Global network and planetary transfer

1. Build global coefficient maps from station clusters.
2. Compare mode-energy distributions by latitude band.
3. Test transfer protocol on Mars/Titan parameter sets.

## Visualization ideas

1. Stability-to-mode heatmaps for each station.
2. Animated globe with mode-energy overlays.
3. Event timeline showing shifts between low and high spectral modes.
4. Tracer-specific closure comparison plots.

## Assessment rubric

1. Physical consistency with MOST limits.
2. Statistical validity (cross-validation and uncertainty).
3. Reproducibility (clean notebook and metadata).
4. Interpretability of modal diagnostics.
5. Communication quality (figures and narrative clarity).

## Capstone prompt

Can a low-order ultraspherical closure with tracer-dependent effective dimension outperform conventional MOST during stable Arctic transition events while preserving interpretability?