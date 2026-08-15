Mathematical Framework: Stochastic Forcing and Eigenvalue Decay Near the Fold Boundary
To rigorously isolate the geometric mechanism behind Critical Slowing Down (CSD) within our adaptive spectral closure, we analyze the linearized state-space dynamics of the reduced system as it approaches the fold line of the slow manifold. We assume the system is driven by high-frequency, small-amplitude turbulent fluctuations, which we model as a stochastic perturbation.

1. Local Normal Form and Time-Scale Separation
Let $Z = (\eta_1, \eta_2, \eta_3)^T \in \mathbb{R}^3$ represent the vector of Gegenbauer spectral coefficients governing the resolved planetary boundary layer (PBL) state. Near a fold parameter threshold, the local dynamics of $Z$ can be decomposed into a fast-slow canonical system via Geometric Singular Perturbation Theory (GSPT).
By shifting the coordinates to the fold line origin, the local deterministic dynamics are governed by the generic saddle-node bifurcation (or fold) normal form coupled to a slowly evolving parameterizing variable:
$$\begin{aligned} d\eta_1 &= \left( \alpha - \eta_1^2 \right) dt \\ d\eta_2 &= -\gamma_2 \eta_2 dt \\ d\alpha &= \epsilon \, g(Z) dt \end{aligned}$$
where:

* $\eta_1$ is the fast dynamic variable representing the primary stability state (e.g., related to the local bulk Richardson number or turbulent kinetic energy closure equilibrium).
* $\eta_2$ captures the remaining orthogonal fast degrees of freedom, which decay rapidly with a stable, large deterministic rate $\gamma_2 > 0$.
* $\alpha$ is the slowly drifting parameter representing macro-scale surface cooling, modulated by the small scale-separation parameter $0 < \epsilon \ll 1$.
The critical manifold $M_0$ is defined in the limit $\epsilon \to 0$ by setting the fast dynamics to zero:
$$M_0 = \left\{ Z \mid \alpha = \eta_1^2, \, \eta_2 = 0 \right\}$$
The manifold consists of two branches meeting at the fold line $\eta_1 = 0$: a stable branch $M_0^+$ where $\eta_1 > 0$, and an unstable branch $M_0^-$ where $\eta_1 < 0$. For $\epsilon > 0$, Fenichel's theorem guarantees the persistence of a perturbed slow manifold $M_\epsilon$ within the domain where $M_0$ remains normally hyperbolic—specifically, away from the immediate neighborhood of the fold line $\eta_1 = 0$.

2. Linearized Jacobian and Eigenvalue Decay
We now consider the system operating on the stable branch $M_\epsilon^+$ under the influence of stochastic environmental forcing. Let $Z^*(t) = (\eta_1^*(t), \eta_2^*(t), \alpha^*(t))^T$ be a reference trajectory tracking the slow evolution along the manifold toward the fold. Linearizing the fast subsystem about $Z^*(t)$ yields the time-dependent Jacobian matrix $J(Z^*)$:
$$J(Z^*) = \begin{pmatrix} \frac{\partial \dot{\eta}_1}{\partial \eta_1} & \frac{\partial \dot{\eta}_1}{\partial \eta_2} \\ \frac{\partial \dot{\eta}_2}{\partial \eta_1} & \frac{\partial \dot{\eta}_2}{\partial \eta_2} \end{pmatrix}_{Z^*} = \begin{pmatrix} -2\eta_1^* & 0 \\ 0 & -\gamma_2 \end{pmatrix}$$
The dominant eigenvalue governing the relaxation rate of perturbations back to the slow manifold is:
$$\lambda_1(t) = -2\eta_1^*(t) = -2\sqrt{\alpha^*(t)}$$
As the macro-scale surface cooling drives the system toward the catastrophic transition point ($\alpha^* \to 0^+$), the dominant stable eigenvalue scales according to a classic power law:
$$\lambda_1(t) \propto \left( \alpha^*(t) \right)^{1/2} \longrightarrow 0^-$$
This mathematical decay of $\lambda_1$ signifies the structural breakdown of normal hyperbolicity. The normal attraction velocity perpendicular to $M_\epsilon$ diminishes, forcing the characteristic relaxation time $\tau \sim |\lambda_1|^{-1}$ to diverge toward infinity.
2. Stochastic Variational Response
To model the high-frequency background noise observed in the SHEBA and SMEAR eddy-covariance data, we introduce a stochastic perturbation to the fast subsystem. The localized state space dynamics are written as a system of Stochastic Differential Equations (SDEs):
$$dZ = f(Z)dt + \sigma dW(t)$$
where $W(t)$ is a standard multi-dimensional Wiener process and $\sigma = \text{diag}(\sigma_1, \sigma_2)$represents the vector of noise intensities.
Let $x(t) = \eta_1(t) - \eta_1^*(t)$ denote a small deviation from the slow manifold trajectory along the primary fast direction. Near the fold, the localized Ornstein-Uhlenbeck (OU) approximation for the deviation $x(t)$ reduces to:
$$dx(t) = \lambda_1(t) x(t) dt + \sigma_1 dW_1(t)$$
Assuming the parameter drift $\dot{\alpha}$ is slow relative to the local relaxation time (justified when $\epsilon \ll |\lambda_1|$), we apply a quasi-stationary approximation. The analytical solution for the state deviation covariance is given by:
$$\langle x(t)^2 \rangle = \int_{-\infty}^{t} \sigma_1^2 e^{2\lambda_1(t)(t-s)} ds = -\frac{\sigma_1^2}{2\lambda_1(t)}$$
Substituting the eigenvalue scaling relation $\lambda_1 = -2\sqrt{\alpha}$ yields the explicit state variance as a function of the distance to the geometric fold:
$$\mathbb{E}[x^2] = \frac{\sigma_1^2}{4\sqrt{\alpha}}$$
Similarly, the lag-$\Delta t$ autocorrelation function $R(\Delta t)$ under this stochastically forced framework is derived as:
$$R(\Delta t) = \frac{\langle x(t)x(t+\Delta t)\rangle}{\langle x(t)^2 \rangle} = e^{\lambda_1(t)\Delta t} = e^{-2\sqrt{\alpha(t)}\Delta t}$$
3. Geometric Grounding of Critical Slowing Down (CSD)
Equations (8) and (9) establish the precise mathematical mapping between the state-space geometry and the classical early warning signals formalized in Table 2:
4. Variance Inflation: As $\alpha \to 0^+$, the denominator in Equation (8) vanishes, driving $\mathbb{E}[x^2] \to \infty$. This variance explosion is purely a geometric consequence of the weakening restorative vector field near the fold boundary, allowing small-scale stochastic fluctuations $\sigma_1$ to accumulate power.
5. Autocorrelation Memory Retention: As $\alpha \to 0^+$, the exponent in Equation (9) approaches zero, forcing $R(\Delta t) \to 1$. The system loses its ability to forget perturbations, creating the statistical illusion of persistent, long-memory stationary structures.
This derivation isolates the foundational diagnostic trap of Paper 1: a localized turbulent burst (an inflation of $\sigma_1$ via structural shear) can mimic the variance inflation of Equation (8) without any reduction in $\alpha$. Therefore, to robustly confirm a upcoming regime transition, the variance inflation must be validated against the explicit trajectory deceleration ($|dZ/dt| \to 0$) derived from the underlying deterministic field.
Ready for Next Phase
With the mathematical equations locked down and the scaling laws ($\mathbb{E}[x^2] \propto \alpha^{-1/2}$) formally established, the theoretical anchor for your Paper 1 results is set.
Shall we proceed immediately to Item 4.2: The Introduction / Narrative Hook, using this exact mathematical landscape to frame the SHEBA/SMEAR scale-dependent roughness motivation?

—-

This is a strong and coherent framework. The progression

1. fold normal form,
2. linearization,
3. stochastic forcing,
4. observable early-warning statistics,

is exactly the order that readers from applied mathematics, fluid dynamics, and atmospheric science expect.

There are, however, several places where tightening the mathematics will make the paper substantially more rigorous.

1. State-space dimensions

There is one inconsistency.

You define

Z=(\eta_1,\eta_2,\eta_3)^T\in\mathbb R^3,

but later write

d\alpha=\epsilon g(Z)\,dt.

Throughout the derivation, the third state variable is actually \alpha, not \eta_3.

I would simply write

Z=(\eta_1,\eta_2,\alpha)^T\in\mathbb R^3.

That eliminates possible confusion when defining the Jacobian.

⸻

1. Fast Jacobian

You correctly linearize only the fast subsystem.

I would explicitly state

f_f(\eta_1,\eta_2;\alpha)=
\begin{pmatrix}
\alpha-\eta_1^2\\
-\gamma_2\eta_2
\end{pmatrix},

with

J_f
=

D_{(\eta_1,\eta_2)}f_f.

That tells the reader immediately why \alpha does not appear in the Jacobian.

Without this clarification, some reviewers may wonder why the full 3\times3 Jacobian wasn’t used.

⸻

1. Normal hyperbolicity

This sentence is excellent:

Fenichel’s theorem guarantees persistence…

I would strengthen it slightly.

Rather than saying

“…within the domain where M_0 remains normally hyperbolic…”

say

“…on every compact subset of M_0 bounded away from the fold where the fast Jacobian possesses no eigenvalues with zero real part.”

This is essentially the precise hypothesis of Neil Fenichel’s theorem.

⸻

1. Eigenvalue scaling

The derivation

\lambda_1=-2\eta_1^*

then

=-2\sqrt{\alpha}

is correct only because you’re restricting to

M_\epsilon^+.

I would remind the reader:

Since the trajectory remains on the attracting branch M_\epsilon^+, we choose the positive root

\eta_1=\sqrt{\alpha}.

That prevents reviewers from asking where the negative branch went.

⸻

1. Breakdown of normal hyperbolicity

I like

“structural breakdown of normal hyperbolicity.”

I would make it slightly more mathematically precise.

Instead of

structural breakdown

say

loss of normal hyperbolicity at the fold.

That is the standard dynamical systems terminology.

⸻

1. OU approximation

This section is particularly good.

One suggestion:

Rather than

localized OU approximation

say

linear stochastic variational equation

followed by

which is an Ornstein–Uhlenbeck process.

That language aligns better with stochastic dynamical systems literature.

⸻

1. Covariance integral

Your calculation

\int_{-\infty}^{t}
\sigma^2
e^{2\lambda(t-s)}
ds
=

-\frac{\sigma^2}{2\lambda}

is correct only because

\lambda<0.

I would state

Since \lambda_1<0 on the attracting branch, the integral converges.

That removes one hidden assumption.

⸻

1. Quasi-stationary approximation

This is perhaps the single most important mathematical assumption.

Instead of

justified when

\epsilon\ll |\lambda_1|

I’d write

\epsilon
=

O(\dot\alpha)
\ll
|\lambda_1|,

meaning

the slow parameter evolves on a timescale much longer than the local relaxation time.

That wording is standard in stochastic slow-fast systems.

⸻

1. Variance scaling

Your result

\mathbb E[x^2]
=

\frac{\sigma^2}{4\sqrt\alpha}

is exactly the classic saddle-node scaling.

One sentence I’d add:

Thus the variance diverges algebraically as the inverse square root of the distance to the fold.

That immediately tells readers this is not an exponential divergence.

⸻

1. Autocorrelation

Excellent.

I might additionally point out

R(\Delta t)
=

e^{-2\sqrt{\alpha}\Delta t}

implies

R(\Delta t)
=

1
-

2\sqrt{\alpha}\,\Delta t
+
O(\alpha),

near the fold.

Atmospheric scientists often think in terms of “autocorrelation approaching one,” and this asymptotic expansion makes that explicit.

⸻

1. Diagnostic trap

This is my favorite part.

The statement

localized turbulent burst can mimic variance inflation

is, in my opinion, the central contribution of Paper 1.

I would sharpen it slightly.

Rather than

statistical illusion

I’d say

statistical degeneracy.

Specifically,

Variance inflation alone is not uniquely identifiable as evidence of an approaching fold because identical variance growth can arise from transient increases in stochastic forcing amplitude.

That language is stronger mathematically.

⸻

1. Connection to adaptive spectral closure

One thing I’d explicitly connect is the geometry to your adaptive basis.

Right now the derivation applies to any reduced-order model.

One paragraph should remind readers:

Because the reduced coordinates are the leading Gegenbauer spectral coefficients, the fold geometry is represented directly in spectral coefficient space rather than in physical space. Consequently, the early-warning indicators emerge as intrinsic properties of the reduced spectral dynamics rather than empirical statistics imposed after projection.

That paragraph ties this entire section back to what makes your framework distinct.

Overall assessment

This section is already at a level suitable for a mathematically oriented applied journal. It combines ideas from Geometric Singular Perturbation Theory, stochastic differential equations, and reduced-order spectral modeling in a coherent way, and the derivation of the \alpha^{-1/2} variance scaling provides a clear theoretical bridge between fold geometry and observable critical slowing down.

The next logical section is indeed the Introduction/Narrative Hook. The mathematics you’ve established gives you a strong foundation to motivate why observations from the Surface Heat Budget of the Arctic Ocean project and SMEAR campaigns exhibit apparent early-warning signals, while also explaining why variance alone can be misleading unless interpreted in the context of the underlying deterministic trajectory through the reduced spectral state space.

---

# 1. Organizing Narrative for Paper 1

To keep the Gegenbauer basis as the organizing principle, the Paper 1 narrative sequence is fixed as:

1. Observational Evidence: establish scale-dependent SBL spectral roughness using high-resolution eddy-covariance data from SHEBA and SMEAR.
2. Methodological Innovation: introduce adaptive estimation of the Gegenbauer parameter $\lambda$ and show how tunable $\lambda$ absorbs anomalous scaling without ad hoc adjustments.
3. Spectral Representation: demonstrate improved capture of SBL spectral structure versus fixed-basis alternatives.
4. Geometric Interpretation: frame regime transitions as a consequence of state-space attractor geometry.
5. Validation and Limitations: validate on out-of-sample intervals and define the model's operational envelope.

---

# 2. Table 2: Conservative State-Space EWS and Falsification

The Table 2 framework is intentionally standalone and publishable as a conservative baseline testing whether classical indicators distinguish true fold approach from background noise.

### Table 2: Early Warning Signal (EWS) Baseline and Falsification Criteria

| Metric | Expected Behavior Approaching Fold | Failure Mode | Falsification Criterion |
| --- | --- | --- | --- |
| Variance | As the dominant stable eigenvalue approaches zero near the folded slow manifold, perturbations decay more slowly, increasing state variance under stochastic forcing. | Intermittent Bursts: high-amplitude localized turbulent excursions mimic variance inflation without regime shift. | Variance increase is not accompanied by corresponding deceleration of the mean state trajectory toward the fold line. |
| Lag-1 Autocorrelation | Near-zero eigenvalue reduces memory decay, causing lag-1 autocorrelation to approach 1. | Stationary Persistence: naturally long-lived stable turbulent structures skew short-window estimates. | Autocorrelation remains high without geometric approach to the attractor fold boundary in state space. |
| Recovery Rate | Return-to-equilibrium rate decreases (critical slowing down). | Measurement Gaps: downtime or filtering artificially flattens the relaxation curve. | Recovery-rate confidence intervals overlap substantially with baseline stationary-regime intervals. |
| Composite Indicator | Integrated metric (variance + lag-1 + recovery) flags high-probability transitions. | False Positives: coincident but decoupled sensor-noise spikes trigger alarms. | Discrimination falls below the ROC/AUC threshold set on the training set. |

---

# 3. Grand Arc: Papers 1 Through 4

* Paper 1 (Current): Adaptive Gegenbauer spectral closure and observational geometry. Establishes adaptive basis, SHEBA/SMEAR validation, and classical state-space EWS diagnostics.
* Paper 2 (Theoretical Progression): Dynamical systems diagnostics of folded slow-manifold transitions. Tests whether finite-time geometric diagnostics detect loss of normal hyperbolicity earlier than classical EWS.
* Paper 3 (Computational Realization): Computational algorithms and synthetic verification. Introduces reduced-state FTLE, finite-time growth rates, and synthetic-system benchmarking.
* Paper 4 (Physical Scaling): Extension to physical-space transport and LES applications. Connects reduced-order state-space geometry to physical transport in $\mathbb{R}^3$ using Haller LCS on LES fields.

---

# 4. Finalized Section 3 Text (Paper 1)

## Mathematical Framework: Stochastic Forcing and Eigenvalue Decay Near the Fold Boundary

To rigorously isolate the geometric mechanism behind Critical Slowing Down (CSD) within our adaptive spectral closure, we analyze the linearized state-space dynamics of the reduced system as it approaches the fold line of the slow manifold. We assume the system is driven by high-frequency, small-amplitude turbulent fluctuations, which we model as a stochastic perturbation.

### 1. Local Normal Form and Time-Scale Separation

Let $Z = (\eta_1, \eta_2, \alpha)^T \in \mathbb{R}^3$ represent the vector of state variables governing the resolved planetary boundary layer (PBL) state. Near a fold parameter threshold, the local dynamics of $Z$ can be decomposed into a fast-slow canonical system via Geometric Singular Perturbation Theory (GSPT). By shifting the coordinates to the fold line origin, the local deterministic dynamics are governed by the generic saddle-node bifurcation normal form coupled to a slowly evolving parameterizing variable:

$$\begin{aligned}
d\eta_1 &= \left( \alpha - \eta_1^2 \right) dt \\
d\eta_2 &= -\gamma_2 \eta_2 dt \\
d\alpha &= \epsilon \, g(Z) dt
\end{aligned}$$

where $\eta_1$ is the fast dynamic variable representing the primary stability state (e.g., related to the local bulk Richardson number or turbulent kinetic energy closure equilibrium), $\eta_2$ captures the remaining orthogonal fast degrees of freedom, which decay rapidly with a stable, large deterministic rate $\gamma_2 > 0$, and $\alpha$ is the slowly drifting parameter representing macro-scale surface cooling, modulated by the small scale-separation parameter $0 < \epsilon \ll 1$.

The critical manifold $M_0$ is defined in the limit $\epsilon \to 0$ by setting the fast dynamics to zero:

$$M_0 = \left\{ Z \in \mathbb{R}^3 \mid \alpha = \eta_1^2, \, \eta_2 = 0 \right\}$$

The manifold consists of two branches meeting at the fold line $\eta_1 = 0$: an attracting branch $M_0^+$ where $\eta_1 > 0$, and an unstable branch $M_0^-$ where $\eta_1 < 0$. For $\epsilon > 0$, Fenichel's theorem guarantees the persistence of a perturbed slow manifold $M_\epsilon$ on every compact subset of $M_0$ bounded away from the fold where the fast Jacobian possesses no eigenvalues with zero real part.

### 2. Fast Jacobian and Eigenvalue Scaling

We isolate the fast subsystem by defining the vector field $f_f(\eta_1, \eta_2; \alpha) = (\alpha - \eta_1^2, -\gamma_2 \eta_2)^T$. We now consider the system operating on the stable branch $M_\epsilon^+$ under the influence of stochastic environmental forcing. Let $Z^*(t) = (\eta_1^*(t), \eta_2^*(t), \alpha^*(t))^T$ be a reference trajectory tracking the slow evolution along the manifold toward the fold. Linearizing only the fast subsystem about $Z^*(t)$ yields the time-dependent fast Jacobian matrix $J_f = D_{(\eta_1, \eta_2)}f_f$:

$$J_f(Z^*) = \begin{pmatrix}
\frac{\partial \dot{\eta}_1}{\partial \eta_1} & \frac{\partial \dot{\eta}_1}{\partial \eta_2} \\
\frac{\partial \dot{\eta}_2}{\partial \eta_1} & \frac{\partial \dot{\eta}_2}{\partial \eta_2}
\end{pmatrix}_{Z^*} = \begin{pmatrix} -2\eta_1^* & 0 \\ 0 & -\gamma_2 \end{pmatrix}$$

The dominant eigenvalue governing the relaxation rate of perturbations back to the slow manifold is given by $\lambda_1(t) = -2\eta_1^*(t)$. Since the reference trajectory remains on the attracting branch $M_\epsilon^+$, we choose the positive root $\eta_1^* = \sqrt{\alpha^*(t)}$, yielding:

$$\lambda_1(t) = -2\sqrt{\alpha^*(t)}$$

As macro-scale surface cooling drives the system toward the catastrophic transition point ($\alpha^* \to 0^+$), the dominant stable eigenvalue scales according to a classic power law:

$$\lambda_1(t) \propto \left( \alpha^*(t) \right)^{1/2} \longrightarrow 0^-$$

This mathematical decay of $\lambda_1$ signifies the loss of normal hyperbolicity at the fold. The normal attraction velocity perpendicular to $M_\epsilon$ diminishes, forcing the characteristic relaxation time $\tau \sim |\lambda_1|^{-1}$ to diverge toward infinity.

### 3. Stochastic Variational Response

To model the high-frequency background noise observed in the SHEBA and SMEAR eddy-covariance data, we introduce a stochastic perturbation to the fast subsystem. Let $x(t) = \eta_1(t) - \eta_1^*(t)$ denote a small deviation from the slow manifold trajectory along the primary fast direction. Near the fold, the linear stochastic variational equation for the deviation $x(t)$ reduces to:

$$dx(t) = \lambda_1(t) x(t) dt + \sigma_1 dW_1(t)$$

where $W_1(t)$ is a standard Wiener process and $\sigma_1$ represents the noise intensity. This formulation describes an Ornstein–Uhlenbeck process. Assuming the slow parameter evolves on a timescale much longer than the local relaxation time ($\epsilon = O(\dot{\alpha}) \ll |\lambda_1|$), we apply a quasi-stationary approximation. Since $\lambda_1 < 0$ on the attracting branch, the covariance integral converges, yielding the analytical solution:

$$\langle x(t)^2 \rangle = \int_{-\infty}^{t} \sigma_1^2 e^{2\lambda_1(t)(t-s)} ds = -\frac{\sigma_1^2}{2\lambda_1(t)}$$

Substituting the eigenvalue scaling relation $\lambda_1 = -2\sqrt{\alpha}$ yields the explicit state variance as a function of the distance to the geometric fold:

$$\mathbb{E}[x^2] = \frac{\sigma_1^2}{4\sqrt{\alpha}}$$

Thus, the variance diverges algebraically as the inverse square root of the distance to the fold. Similarly, the lag-$\Delta t$ autocorrelation function $R(\Delta t)$ is derived as:

$$R(\Delta t) = \frac{\langle x(t)x(t+\Delta t)\rangle}{\langle x(t)^2 \rangle} = e^{\lambda_1(t)\Delta t} = e^{-2\sqrt{\alpha(t)}\Delta t}$$

Expanding this expression asymptotically near the fold ($\alpha \to 0^+$) yields:

$$R(\Delta t) = 1 - 2\sqrt{\alpha}\,\Delta t + O(\alpha)$$

### 4. Connection to Adaptive Spectral Closure and the Diagnostic Trap

Because the reduced coordinates $(\eta_1, \eta_2)$ are the leading Gegenbauer spectral coefficients, the fold geometry is represented directly in spectral coefficient space rather than in physical space. Consequently, the early-warning indicators emerge as intrinsic properties of the reduced spectral dynamics rather than empirical statistics imposed after projection. Equations (8) and (9) establish the precise mathematical mapping between this spectral state-space geometry and the classical indicators formalized in Table 2: variance inflation ($\mathbb{E}[x^2] \to \infty$) and memory retention ($R(\Delta t) \to 1$).

However, this mathematical clarity reveals a profound statistical degeneracy that serves as the central diagnostic trap of boundary layer forecasting. Variance inflation alone is not uniquely identifiable as evidence of an approaching fold because identical variance growth can arise from transient increases in stochastic forcing amplitude ($\sigma_1$), such as a localized turbulent burst driven by local wind shear. A true geometric approach to the fold requires a concurrent structural change in the underlying deterministic field. Therefore, within this adaptive spectral framework, variance inflation cannot be interpreted in isolation; it must be cross-validated against the explicit trajectory deceleration ($|dZ/dt| \to 0$) along the mean state trajectory toward the fold line.

---

# 5. Immediate Next-Write: Introduction / Narrative Hook (Draft)

Below is the draft for the **Introduction/Narrative Hook** section of Paper 1, using the specific SHEBA/SMEAR observational constraints to motivate the adaptive Gegenbauer framework.

## 1. Introduction and Observational Motivation

The parameterization of the Stable Boundary Layer (SBL) remains one of the most persistent sources of uncertainty in modern numerical weather prediction and climate modeling. When solar radiation diminishes and the surface undergoes net longwave cooling, the boundary layer frequently transitions from a well-mixed convective regime into a highly stratified state. Under strongly stable conditions, turbulence becomes weak, intermittent, and localized, often decoupling from the surface entirely. Classical boundary layer theories founded on Reynolds averaging and isotropic Monin-Obukhov Similarity Theory (MOST) break down in this regime, as they fail to capture the non-local transport dynamics, gravity wave interactions, and abrupt regime transitions that characterize the nocturnal SBL.

This modeling failure is intimately tied to the complex spectral structure of SBL turbulence. High-resolution eddy-covariance datasets from long-term field campaigns—specifically the Surface Heat Budget of the Arctic Ocean (SHEBA) project and the Station for Measuring Ecosystem-Atmosphere Relations (SMEAR) campaigns—reveal a highly variable and scale-dependent spectral roughness. Unlike the predictable $-5/3$ Kolmogorov scaling observed in convective conditions, SBL velocity and temperature spectra exhibit anomalous scaling laws that continuously shift in response to the ambient stratification. Fixed-basis spectral projections or traditional grid-based closures fail to adapt to these shifting multi-scale features, generating numerical grid-dependencies or requiring highly tuned, ad-hoc stability functions to prevent the spurious runaway cooling of the modeled surface.

To overcome these limitations, this paper introduces a low-dimensional spectral closure rooted in an adaptive **Gegenbauer polynomial basis**. Rather than forcing the atmospheric state onto fixed spatial grids or rigid orthogonal modes, we project the governing momentum, thermodynamic, and tracer equations into a unified Hilbert space $\mathcal{H}$, where the basis functions are modulated by a tunable parameter $\lambda$. By continuously estimating $\lambda$ directly from high-frequency eddy-covariance metrics, the underlying mathematical basis adaptively shifts its weighting to match the exact spectral roughness and scale-dependent anisotropy observed in the SHEBA and SMEAR data.

Crucially, this adaptive spectral projection transforms the problem of SBL regime transitions from an empirical parameterization challenge into an unyielding problem of state-space geometry. By treating macro-scale surface cooling as a slow parameter driving a reduced-order system of spectral coefficients, the sudden collapse or recovery of turbulence emerges naturally as a transition across a folded slow manifold. While classical statistical indicators such as variance inflation and critical slowing down have been floated as early warning signals for these boundary layer transitions, we leverage our geometric framework to expose a significant statistical degeneracy: localized, transient turbulent bursts can mimic traditional early-warning metrics without a true geometric approach to a tipping point.

By anchoring our mathematical framework in the deterministic trajectories of spectral coefficient space, this paper provides both an optimal spectral representation of the stratified boundary layer and a mathematically rigorous baseline for distinguishing genuine regime transitions from transient background noise.

---

### Strategy Review

The Introduction is now cleanly aligned with Section 3, explicitly preparing the reader for the structural properties of the state-space geometry while keeping the focus on the **Gegenbauer basis** as the primary methodological contribution.

Would you like to expand the **Observational Evidence section** next to detail exactly how the $\lambda$ parameter is estimated from the SHEBA and SMEAR datasets, or should we refine the **Validation and Limitations section**?
