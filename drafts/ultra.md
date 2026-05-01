# Ultraspherical Closure Strategy for the Highly Stable Nocturnal Boundary Layer Prior to SHEBA Calibration

The primary obstacle in extending Monin-Obukhov similarity theory (MOST) into the highly stable nocturnal boundary layer (HSNBL) is not a lack of functional forms, but a mismatch between the assumptions of local equilibrium theory and the actual dynamics of very stable flow. The canonical Kansas and Businger-Dyer relations were derived and calibrated in regimes where turbulent transport remains sufficiently continuous for a universal local scaling to be meaningful [@Businger_1971; @Dyer1974; @Hogstrom1988]. In the HSNBL, however, turbulence becomes intermittent, the separation between turbulent and submeso motions becomes blurred, and gravity-wave or drainage-flow contributions can contaminate the flux-gradient relationship. Observational analyses of stable boundary-layer data, especially the SHEBA studies and later reassessments of local similarity limits, show that the classical linear stable forms lose explanatory power as $\zeta$ increases beyond order unity [@Grachev2005; @GrachevEtAl2007; @Grachev_2012b; @GryanikEtAl2020].

The ultraspherical approach is motivated by the observation that the unstable MOST power laws,

$$
\phi_m(\zeta) = (1-b_m\zeta)^{-1/4},
\qquad
\phi_h(\zeta) = a_h^{-1}(1-b_h\zeta)^{-1/2},
$$

are not merely empirical fits but algebraically structured generating functions. In particular, the momentum exponent $1/4$ and the heat exponent $1/2$ place the two similarity functions in the Gegenbauer family, with heat occupying the Legendre special case. This observation suggests that similarity functions should be treated as members of a spectral hierarchy rather than as unrelated curve fits. The immediate practical consequence is that closure error can be represented as modal structure: once a physically constrained baseline is chosen, the remaining discrepancy can be expanded in a low-order ultraspherical basis rather than absorbed into ad hoc changes of empirical constants.

We therefore propose to write the observed similarity function in the form

$$
\phi_q^{\mathrm{obs}}(\zeta) =
\phi_q^{\mathrm{base}}(\zeta)
+
\sum_{n=0}^{N} c_n\,C_n^{(\lambda^*)}(\xi(\zeta)),
$$

where $\phi_q^{\mathrm{base}}$ carries the known neutral-limit and asymptotic behavior, $\xi(\zeta)$ is a bounded stability coordinate, and the coefficients $c_n$ represent the residual structure not captured by the baseline. This decomposition is attractive for two reasons. First, it respects the fact that some parts of the closure problem are already constrained by similarity theory and should not be relearned from noisy data. Second, it gives a compact mechanism for representing regime-dependent curvature, intermittency, and non-local transport effects without sacrificing interpretability. The role of the ultraspherical correction is therefore not to replace MOST wholesale, but to repair the part of the closure that MOST leaves unresolved.

The parameter $\lambda^*$ may also be interpreted as a spectral-dimension coordinate through the formal relation

$$
\lambda = \frac{d-2}{2},
\qquad
d = 2 + 2\lambda,
$$

so that the momentum-like case $\lambda = 1/4$ corresponds to an effective dimension $d = 5/2$, while the Legendre case $\lambda = 1/2$ corresponds to $d = 3$. We do not regard this dimensional interpretation as established atmospheric doctrine; rather, it is a useful organizing hypothesis for the present program. In that interpretation, a fitted value $\lambda^* < 1/4$ in very stable conditions would indicate a reduction in effective transport dimensionality as turbulence collapses toward strongly anisotropic, intermittent structures. This is precisely the kind of physics that is difficult to encode in a one-parameter linear stable branch but is naturally expressed in a spectral framework.

A practical consequence of this viewpoint is that the HSNBL problem should be approached in two stages. Before SHEBA data are in hand, the correct task is not to claim an Arctic-calibrated similarity law, but to establish that the ultraspherical closure architecture is statistically and physically viable on generic stable-capable tower data. For that phase, a simple baseline of the form

$$
\phi(\zeta) = a\,(1-b\zeta)^{-1/\lambda_q}
$$

is sufficient, provided the fit is validated out of sample and the ultraspherical correction is kept low order. Once strong-stability data become available, the baseline can be upgraded to a SHEBA-type form,

$$
\phi_q(\zeta) = 1 + \frac{a_q\zeta(1+\zeta)^{1/3}}{1+b_q\zeta},
$$

which is specifically designed to represent the flattening and partial saturation observed in the very stable Arctic boundary layer [@GrachevEtAl2007; @GryanikEtAl2020]. In other words, the absence of SHEBA today does not block the theory; it simply limits the strength of the asymptotic claims we can responsibly make.

The choice of stability coordinate is central to making the spectral correction numerically useful. For generic station data spanning moderate stability, the linear compactification

$$
\xi = \tanh(\alpha_\xi \zeta)
$$

is adequate. For HSNBL data, where $\zeta$ may span multiple decades, a logarithmic compactification is preferable,

$$
\xi = \tanh\!\bigl(\alpha_\xi\ln(1+\zeta)\bigr),
$$

because it allocates spectral resolution more evenly across weakly stable, transitional, and very stable regimes. This distinction is already reflected in the software design: the generic Julia prototype uses the simpler tanh map as a first-pass educational and diagnostic tool, whereas the HSNBL-specific driver adopts the log-tanh map together with stability weighting and blocked validation. That division of labor is scientifically appropriate. The first script asks whether low-order spectral residuals improve a baseline closure at all; the second asks whether a strong-stability baseline plus spectral correction can explain the HSNBL tail without overfitting it.

This staging also clarifies how the present work should be described to students and collaborators. The generic workflow should be presented as a proof-of-concept residual-closure experiment on any stable-capable site, not as a final Arctic closure. The SHEBA-oriented workflow should be reserved for data sets that genuinely probe the very stable regime. In practical terms, the criterion for success at the current stage is straightforward: the ultraspherical correction should improve held-out error with a small number of modes, stable hyperparameters, and a physically smooth residual. If the improvement appears only at large mode counts or under unstable parameter estimates, then the model is fitting noise rather than transport structure.

The value of the present framework is therefore not that it already solves the HSNBL closure problem, but that it organizes the next steps coherently. It preserves the analytical structure emphasized in the similarity-theory tradition [@EnglandMcNider1995], admits a direct path to HSNBL-specific asymptotics once Arctic data are available [@GrachevEtAl2007; @Grachev_2012b], and aligns with the broader need to generalize MOST when turbulence departs from local equilibrium [@StiperskiCalaf2023]. For the McNider-Biazar program, the defensible scientific claim at this stage is that ultraspherical closure provides a mathematically structured and testable extension of MOST, one that can be calibrated incrementally now and specialized to the Arctic HSNBL when the appropriate benchmark data arrive.

