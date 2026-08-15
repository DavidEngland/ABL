Curvature of the Gradient Richardson Number (citation-ready note)

\zeta\equiv\frac{z}{L},\qquad
\phi_m(\zeta)=(1-\beta_m\zeta)^{-\alpha_m},\qquad
\phi_h(\zeta)=(1-\beta_h\zeta)^{-\alpha_h},
with 1-\beta_{(\cdot)}\zeta>0 in the domain of interest.
MOST gradient Richardson number (MOST form)
Ri_g(\zeta)=\zeta\frac{\phi_h(\zeta)}{\phi_m(\zeta)^2}\equiv \zeta\,F(\zeta).

Logarithmic derivative shorthand

Define
V_{\log}(\zeta)=\frac{1}{F}\frac{dF}{d\zeta}
=\frac{1}{\phi_h}\frac{d\phi_h}{d\zeta}-\frac{2}{\phi_m}\frac{d\phi_m}{d\zeta}
=\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta},
W_{\log}(\zeta)=\frac{dV_{\log}}{d\zeta}
=\frac{\alpha_h\beta_h^{2}}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^{2}}{(1-\beta_m\zeta)^2}.
(Useful identity: for \phi=(1-\beta\zeta)^{-\alpha}, \phi’/\phi=\alpha\beta/(1-\beta\zeta).)

Curvature (compact form)

From Ri_g=\zeta F:
\frac{dRi_g}{d\zeta}=F+\zeta F’,\qquad
\frac{d^{2}Ri_g}{d\zeta^{2}}=2F’+\zeta F’’,
and using F’/F=V_{\log}, F’’/F=V_{\log}^2-W_{\log}, we obtain the fundamental boxed expression
\boxed{\frac{d^{2}Ri_g}{d\zeta^{2}}=F(\zeta)\Big[\,2V_{\log}(\zeta)+\zeta\big(V_{\log}(\zeta)^2-W_{\log}(\zeta)\big)\Big].}

Equivalently, expanded in the power-law parameters:
\begin{aligned}
\frac{d^{2}Ri_g}{d\zeta^{2}}&=
(1-\beta_h \zeta)^{-\alpha_h}(1-\beta_m \zeta)^{2\alpha_m}\times\\
&\quad\left\{\frac{2\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{4\alpha_m\beta_m}{1-\beta_m\zeta}
+\zeta\Big[\Big(\frac{\alpha_h\beta_h}{1-\beta_h\zeta}-\frac{2\alpha_m\beta_m}{1-\beta_m\zeta}\Big)^2\right.\\
&\qquad\left.-\Big(\frac{\alpha_h\beta_h^{2}}{(1-\beta_h\zeta)^2}-\frac{2\alpha_m\beta_m^{2}}{(1-\beta_m\zeta)^2}\Big)\Big]\right\}.
\end{aligned}

Neutral limit (\zeta\to0)

Let
\Delta\equiv\alpha_h\beta_h-2\alpha_m\beta_m,\qquad
c_1\equiv\alpha_h\beta_h^{2}-2\alpha_m\beta_m^{2}.
Then
F(0)=1,\quad V_{\log}(0)=\Delta,\quad W_{\log}(0)=c_1,
and
\boxed{\left.\frac{d^{2}Ri_g}{d\zeta^{2}}\right|_{\zeta=0}=2\Delta.}
Interpretation:
	•	\Delta>0: initial concave-up Ri_g(\zeta) — heat stabilization corrections dominate.
	•	\Delta<0: initial concave-down — momentum corrections dominate.
	•	\Delta\approx0: quasi-linear near onset.

Series (small \zeta):
\ln F(\zeta)=\Delta\zeta+\tfrac12 c_1\zeta^2+O(\zeta^3),
F(\zeta)=1+\Delta\zeta+\tfrac12(\Delta^2+c_1)\zeta^2+O(\zeta^3),
so
Ri_g(\zeta)=\zeta+\Delta\zeta^2+O(\zeta^3),\qquad
\frac{d^{2}Ri_g}{d\zeta^{2}}=2\Delta+O(\zeta).

Singular behavior and domain

Poles at \zeta_h=1/\beta_h and \zeta_m=1/\beta_m. As \zeta\to\zeta_h^- (assuming \beta_h>\beta_m),
\frac{d^{2}Ri_g}{d\zeta^{2}}\sim (\beta_h(\zeta_h-\zeta))^{-\alpha_h}\,C\Big[\frac{2\alpha_h}{\zeta_h-\zeta}+O(1)\Big],
so curvature exhibits algebraic blow-up with exponent \alpha_h+1. In practice restrict \zeta away from poles (e.g. \zeta\lesssim 0.7/\max(\beta)) or use pole-free surrogate forms (Q-SBL, RPL, etc.) for large \zeta.

Physical height conversion (chain rule)

Treating local L constant,
\boxed{\frac{\partial^2 Ri_g}{\partial z^2}=\frac{1}{L^2}\frac{d^{2}Ri_g}{d\zeta^{2}}.}
If L=L(z) is allowed to vary, extra terms involving dL/dz (and d^2L/dz^2) appear; for MOST diagnostics the local-L approximation above is standard.

Inflection (interior root)

Inflection occurs when
2V_{\log}(\zeta)+\zeta\big(V_{\log}^2(\zeta)-W_{\log}(\zeta)\big)=0.
Near neutrality (|\zeta|\ll1) expand V_{\log}=\Delta+c_1\zeta+O(\zeta^2). Retaining leading orders gives the quadratic approximation for the first interior root (if it exists and is small):
\zeta_{\mathrm{inf}}\approx -\frac{2\Delta}{\Delta^2-c_1},
valid only when |\zeta_{\mathrm{inf}}|\ll\min(1/\beta_h,1/\beta_m).

Third derivative (initial change of curvature)

With the same notation,
\frac{d^3Ri_g}{d\zeta^3}=F\Big[3(V_{\log}^2-W_{\log})+\zeta\big(V_{\log}^3-3V_{\log}W_{\log}-W_{\log}’\big)\Big],
and at neutrality
\left.\frac{d^3Ri_g}{d\zeta^3}\right|_{0}=3(\Delta^2-c_1),
so the sign of \Delta^2-c_1 controls the initial evolution of curvature with \zeta.

Turbulent Prandtl and relations

Pr_t(\zeta)=\frac{\phi_h}{\phi_m}=(1-\beta_h\zeta)^{-\alpha_h}(1-\beta_m\zeta)^{\alpha_m},
Pr_t(\zeta)=1+(\alpha_h\beta_h-\alpha_m\beta_m)\zeta+O(\zeta^2).
Note \Delta=(\alpha_h\beta_h-\alpha_m\beta_m)-(\alpha_m\beta_m); thus \Delta=0 (curvature neutral) does not imply Pr_t=1.

# Turbulent Prandtl Link (Addendum)

Given $Pr_t=\phi_h/\phi_m$ and $Ri_g=\zeta \phi_h/\phi_m^2$:
Near-neutral series:
$$Pr_t=1+(a_h-a_m)\zeta + O(\zeta^2),\quad Ri_g=\zeta+\Delta\zeta^2+O(\zeta^3).$$
Eliminate ζ:
$$Pr_t(Ri_g)\approx 1+(a_h-a_m)Ri_g - (a_h-a_m)\Delta Ri_g^2.$$
Second derivative:
$$\frac{d^2Pr_t}{dRi_g^2}\Big|_0=-2(a_h-a_m)\Delta.$$
Curvature correction leaving $\Delta$ unchanged preserves this intrinsic Pr_t response.

Scalar closure addendum: $K_q=\kappa z u_*/\phi_q$, $f_q(Ri_g)=1/(\phi_m\phi_q)$ shares curvature-driven bias via $\Delta$ in series inversion.

Diagnostic & reporting suggestions

Report as a minimal set for comparing parameterizations:
	1.	Neutral curvature: \left.\dfrac{d^{2}Ri_g}{d\zeta^{2}}\right|_{\zeta=0}=2\Delta.
	2.	First inflection \zeta_{\mathrm{inf}} (if interior and within acceptance band).
	3.	Fractional enhancement factor
\mathcal{C}(\zeta)=\frac{(d^{2}Ri_g/d\zeta^{2})}{2\Delta}.
	4.	Map to physical height via \partial_z^2=(1/L^2)\partial_\zeta^2 and report using representative local L values.

Figure suggestions (manuscript)
	•	(i) Ri_g(\zeta) and its curvature vs. \zeta for several parameter sets (show poles and Q-SBL surrogate).
	•	(ii) \mathcal{C}(\zeta) normalized to the neutral curvature.
	•	(iii) Sensitivity heat map: \partial (d^{2}Ri_g/d\zeta^{2})/\partial\alpha_{h,m} and /\partial\beta_{h,m} on a reasonable grid.
	•	(iv) Inversion example \zeta(Ri_g) near neutrality comparing Taylor inversion vs. Newton refine.

⸻

Short summary (one-line)

\boxed{\displaystyle \frac{d^{2}Ri_g}{d\zeta^{2}}=F(\zeta)\Big[2V_{\log}+\zeta(V_{\log}^2-W_{\log})\Big]} — neutral curvature equals 2(\alpha_h\beta_h-2\alpha_m\beta_m); proximity to parameter poles controls rapid curvature growth and parameter sensitivity.

⸻

If you want, I can immediately:
	•	produce a ready-to-include LaTeX snippet (with \begin{equation} / \label{} and citation placeholders) — formatted and minimal; or
	•	generate publication-quality plots (curvature vs \zeta, \mathcal{C}(\zeta), sensitivity maps) for a set of parameter choices you prefer (give me a few parameter sets), or
	•	expand the parameter-sensitivity derivatives into explicit symbolic forms for \partial/\partial\alpha_h,\partial/\partial\beta_m etc.

Which of those would you like first?