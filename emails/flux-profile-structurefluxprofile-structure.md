Unified MOST parameter framework and flux–profile structure

Draft note for McNider & Biazar

“The table you provided captures the core empirical constants used in MOST‑based surface‑layer schemes. These constants define:
\(\phi_m = (1 - \gamma_m \zeta)^{-1/4}, \quad \phi_h = a_h^{-1}(1 - \gamma_h \zeta)^{-1/2}\) (unstable) and
\(\phi_m = 1 + \beta_m \zeta, \quad \phi_h = a_h^{-1} + \beta_h \zeta\) (stable).”

---

1. Master parameter and regime table

1.1 Canonical Businger–Dyer type sets

Family / Source	Branch	\(k\)	\(a_h^{-1}\)	\(b_m=\gamma_m\)	\(b_h=\gamma_h\)	\(\beta_m\)	\(\beta_h\)	Nominal \(\zeta\) range	Notes
Businger et al. (1971)	Unstable	0.35	1.35	15.0	9.0	–	–	\(-2 \lesssim \zeta < 0\)	Kansas baseline; strong scalar–momentum mismatch
Businger et al. (1971)	Stable	0.35	1.35	–	–	4.7	6.35	\(0 < \zeta \lesssim 1\)	Slightly asymmetric stable slopes
Dyer (1974) / Brutsaert (1982)	Unstable	0.40	1.00	16.0	16.0	–	–	\(-2 \lesssim \zeta < 0\)	Degenerate case \(b_m=b_h\Rightarrow \phi_h=\phi_m^2\)
Dyer (1974) / Brutsaert (1982)	Stable	0.40	1.00	–	–	5.0	5.0	\(0 < \zeta \lesssim 1\)	Symmetric, widely used in WRF/LES
Högström (1996)	Unstable	0.40	1.00	19.0	11.6	–	–	\(-2 \lesssim \zeta < 0\)	Stronger momentum curvature
Högström (1996)	Stable	0.40	1.00	–	–	5.3	5.3	\(0 < \zeta \lesssim 1\)	Slightly steeper stable branch


Regime structure (for reference)

Regime	\(\zeta\) range	Ri range	Physical character
Near‑neutral unstable	\(-1 < \zeta < 0\)	\(Ri \lesssim 0\)	Weak buoyant enhancement; log‑law dominant
Strongly unstable	\(\zeta < -1\)	\(Ri \ll 0\)	Convective plumes; free‑convection tendencies
Moderately stable	\(0 < \zeta \le 1\)	\(0.02 < Ri < 0.12\)	Suppressed mixing, still continuous turbulence
Very/extremely stable	\(\zeta > 1\)	\(0.12 < Ri < 0.7\)	Intermittent, z‑less turbulence; decoupling risk
Critical	–	\(Ri_c \approx 0.21\text{–}0.25\)	Onset of laminarization / flux collapse


These ranges are the natural “operating window” for Kansas‑type similarity; beyond them, the physics becomes regime‑dependent (intermittency, wave–turbulence coupling, drainage).

---

1.2 Specialized stable‑layer formulations

“These schemes were developed to avoid the decoupling problem (runaway cooling, vanishing fluxes) that plagues linear MOST functions at high stability.”

We summarize them in a common parameter language. Here \(\phi = 1 + a\,\zeta + b\,\zeta^2 + \ldots\) is schematic; each scheme has its own exact form, but the constants below control curvature and asymptotic behavior.

Scheme	Variable	Core constants	Valid \(\zeta\)	Physics / usage
Beljaars–Holtslag (1991)	m, h	\(a=1,\; b=0.667,\; c=5,\; d=0.35\)	\(0 < \zeta \lesssim 10\)	Exponential tail; smooths transition into very stable SBL; avoids flux shutdown in ECMWF‑type schemes
Cheng & Brutsaert (2005)	m	\(a=6.1,\; b=2.5\)	\(0 < \zeta \lesssim 20\)	Stronger curvature; tuned to extended stable datasets
Cheng & Brutsaert (2005)	h	\(c=5.3,\; d=1.1\)	\(0 < \zeta \lesssim 20\)	Scalar more strongly damped than momentum; consistent with nocturnal inversions
SHEBA (2007)	m	\(a_m=5,\; b_m=0.3\)	\(0 < \zeta \lesssim 100\)	Arctic sea‑ice SBL; captures z‑less, intermittent turbulence
SHEBA (2007)	h	\(a_h=5,\; b_h=0.4\)	\(0 < \zeta \lesssim 100\)	Strong scalar suppression; essential for polar inversions


For your Arctic inversion experiment, these specialized forms are the natural candidates for the “very/extremely stable” branch, with Kansas‑type linear forms reserved for \(0 < \zeta \lesssim 1\).

---

2. Flux–profile relationships and physics

2.1 MOST backbone

Define the usual similarity variables

\zeta = \frac{z}{L}, \qquad
\phi_m(\zeta) = \frac{\kappa z}{u_*}\frac{\partial U}{\partial z}, \qquad
\phi_h(\zeta) = \frac{\kappa z}{\theta_*}\frac{\partial \Theta}{\partial z},


with \(u_*\) the friction velocity, \(\theta_* = -\overline{w^\prime\theta^\prime}/u_*\), and \(L\) the Obukhov length.

The flux–profile integrals are

\psi_m(\zeta) = \int_0^\zeta \left(1 - \frac{1}{\phi_m(\zeta^\prime)}\right)\,\frac{d\zeta^\prime}{\zeta^\prime},\qquad
\psi_h(\zeta) = \int_0^\zeta \left(1 - \frac{1}{\phi_h(\zeta^\prime)}\right)\,\frac{d\zeta^\prime}{\zeta^\prime}.


Then, for a reference height \(z\) and roughness lengths \(z_{0m}, z_{0h}\),

U(z) - U(z_{0m}) = \frac{u_*}{\kappa}\Big[\ln\!\frac{z}{z_{0m}} - \psi_m(\zeta) + \psi_m(\zeta_0)\Big],


\Theta(z) - \Theta(z_{0h}) = \frac{\theta_*}{\kappa}\Big[\ln\!\frac{z}{z_{0h}} - \psi_h(\zeta) + \psi_h(\zeta_0)\Big].


The physics is encoded entirely in \(\phi\) (local mixing efficiency) and thus in \(\psi\) (non‑local integrated effect).

---

2.2 Unstable Businger–Dyer family

For the canonical power‑law forms

\phi_m(\zeta) = (1 - b_m \zeta)^{-1/4}, \qquad
\phi_h(\zeta) = a_h^{-1}(1 - b_h \zeta)^{-1/2}, \qquad \zeta < 0,


define

x = (1 - b_m \zeta)^{1/4}, \qquad
y = (1 - b_h \zeta)^{1/2}.


The standard integrals (Kansas forms) are:

• Momentum\psi_m(\zeta) = 2\ln\!\frac{1+x}{2} + \ln\!\frac{1+x^2}{2} - 2\arctan x + \frac{\pi}{2}.
This is the familiar “Kansas \(\psi_m\)” used in WRF and many LES codes when \(b_m=16\).
• Heat\psi_h(\zeta) = 2\ln\!\frac{1+y}{2}.
When \(a_h^{-1}=1\) and \(b_h=16\), this reduces to the classical Dyer scalar integral.


Physics injection

• The exponents \(-1/4\) and \(-1/2\) encode how quickly the mixing length “inflates” with instability. Momentum responds more weakly (\(\phi_m \sim (-\zeta)^{-1/4}\)) than scalars (\(\phi_h \sim (-\zeta)^{-1/2}\)), so velocity profiles remain closer to log‑law than temperature profiles in strong convection.
• The coefficients \(b_m, b_h\) set the curvature scale: \(|\zeta| \sim 1/b\) is where the neutral‑regime CBC series breaks down and the UBL asymptotics take over. This is a natural, parameter‑consistent definition of a “dynamic” convective limit.


In the degenerate Dyer/Brutsaert case \(b_m=b_h=b\) and \(a_h^{-1}=1\),

\phi_h = \phi_m^2,\qquad
Ri_g = \zeta,


so the gradient Richardson number equals the stability parameter exactly on the unstable branch—no iterative inversion. This is a structurally clean limit case that you can exploit analytically in testing inversion algorithms.

---

2.3 Stable linear branch (Kansas‑type)

For the classical linear forms

\phi_m(\zeta) = 1 + \beta_m \zeta,\qquad
\phi_h(\zeta) = a_h^{-1} + \beta_h \zeta,\qquad \zeta > 0,


the integrals are elementary.

• Momentum\psi_m(\zeta) = -\beta_m \zeta.
This is the usual linear stable correction: the mean wind profile steepens linearly with \(\zeta\), reflecting reduced turbulent mixing.
• HeatIt is convenient to factor out the neutral limit:\phi_h(\zeta) = a_h^{-1}(1 + a_h \beta_h \zeta),
so that\psi_h(\zeta) = -a_h \beta_h \zeta.



Physics injection

• Linear \(\phi\) implies a constant reduction of eddy diffusivity with height‑scaled stability: \(K_m \propto u_* z / \phi_m\) decays roughly as \(1/(1+\beta_m \zeta)\). This is adequate for \(0<\zeta\lesssim 1\), where turbulence is still continuous and wave activity modest.
• Beyond \(\zeta\sim 1\), linear forms over‑suppress mixing, leading to runaway surface cooling and decoupling—exactly the nocturnal SBL pathology that motivated BH91, CB05, and SHEBA.


---

2.4 Beljaars–Holtslag, Cheng–Brutsaert, SHEBA: stable flux profiles

Each of these schemes modifies \(\phi(\zeta)\) so that:

• Near neutral (\(\zeta\to 0^+\)), they reduce to Kansas‑type linear behavior with \(\phi \approx 1 + \mathcal{O}(\zeta)\).
• At large \(\zeta\), they saturate or grow sub‑linearly, preventing \(K\to 0\) and maintaining finite fluxes in very stable conditions.


A generic BH91‑type structure can be written as

\phi_m(\zeta) = 1 + a\frac{\zeta}{1 + b\zeta},\qquad
\phi_h(\zeta) = 1 + c\frac{\zeta}{1 + d\zeta},


with \((a,b,c,d)\) as in the table. The corresponding integrals are

\psi_m(\zeta) = -a\ln(1 + b\zeta)/b,\qquad
\psi_h(\zeta) = -c\ln(1 + d\zeta)/d,


up to neutral‑matching constants.

For CB05 and SHEBA, the exact analytic forms are more cumbersome, but the same pattern holds:

• Near neutral: \(\psi \sim -\beta \zeta\), matching Kansas.
• Very stable: \(\psi\) grows only logarithmically or sub‑linearly with \(\zeta\), so \(K\) asymptotes to a finite fraction of its neutral value instead of collapsing.


Physics injection

• These forms effectively impose a “floor” on turbulent mixing controlled by the constants \(b,d\): as \(\zeta\to\infty\), \(\phi\to 1+a/b\) (or analogous), so \(K\) remains finite. This is a parameterized representation of intermittent, patchy turbulence that persists even in very strong stratification.
• In Arctic inversions and nocturnal SBLs, this prevents the model from entering an unrealistically laminar state where radiative cooling runs away and the surface layer decouples from the air aloft—exactly the failure mode McNider‑style drainage and column models are sensitive to.


---

2.5 Linking to dynamic `\(Ri_c\)` and inversion experiments

The structural picture that emerges:

• Unstable side: The CBC expansion and UBL asymptotics define a natural lower bound \(Ri_c^{UBL} \sim -1/b_m\), tied to the branch point of \(\phi_m\).
• Stable side: Linear or BH91/CB05/SHEBA forms define an upper bound \(Ri_c^{SBL} \sim 1/\beta\) (or its generalized analogue), tied to the pole or saturation scale of \(\phi\).


A dynamic critical Richardson number \(Ri_c^*\) that varies smoothly between these limits,

Ri_c^* \in [-1/b_m,\; 1/\beta_{\text{eff}}],


and depends on SBL state variables (shear, stratification, TKE) is then not an ad‑hoc switch but a direct consequence of the chosen similarity family and its analytic structure.

For your Arctic stability‑function inversion:

• The unified table above gives a clean parameter space to explore (Kansas vs BH91 vs CB05 vs SHEBA).
• The flux–profile integrals provide the forward operator from \((u_*,\theta_*,L)\) to observed profiles.
• The physics‑aware bounds on \(Ri_c^*\) give a principled way to regularize the inversion in regimes where classical MOST is formally invalid but turbulence persists.


If you’d like, the next step could be: (i) a compact LaTeX table and appendix‑ready derivation of the \(\psi\) functions, or (ii) a small Python library that implements these families with regime‑aware selection and exposes \(Ri_c^*\) as a diagnostic.