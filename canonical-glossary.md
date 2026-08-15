ABL Closure Project — Canonical Glossary

Term	Symbol	Definition & units	Project context / notes
Obukhov length	L	Scale where shear-driven turbulence and buoyancy balance (m).	Fundamental length scale defining SBL stability; used to form \zeta=z/L.
Dimensionless height	\zeta	\zeta = z/L (unitless).	Non-dimensional height; \zeta>0 indicates stable conditions; used in MOST functions.
Geometric mean height	z_g	z_g=\sqrt{z_0 z_1} (m).	Representative height for log-like profiles; often where point-value of \mathrm{Ri}_g is evaluated.
Logarithmic mean height	z_L	z_L=(z_1-z_0)/\ln(z_1/z_0) (m).	Exact representative height for velocity difference \Delta U in a log wind law.
Gradient Richardson number (local)	\mathrm{Ri}_g	\displaystyle \mathrm{Ri}_g(z)=\frac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2} (unitless).	Local stability diagnostic; used to identify collapse point Ri_c.
Bulk Richardson number (layer)	\mathrm{Ri}_b	\displaystyle \mathrm{Ri}_b=\frac{(g/\bar\theta)\,\Delta\theta\;\Delta z}{(\Delta U)^2} (unitless).	Layer-averaged stability diagnostic; typical input for NWP closures.
Bias ratio	B	B=\mathrm{Ri}_g(z_g)/\mathrm{Ri}_b (unitless).	Primary project metric. B>1 indicates \mathrm{Ri}_b underestimates local stability.
Grid damping factor	G	G(\zeta,\Delta z) (multiplier, unitless).	ML target: multiplicative modifier applied to eddy diffusivity (K_{\text{new}}=K_{\text{old}}\cdot G).
Eddy diffusivities	K_m, K_h	Momentum and heat diffusivities (m^2 s^{-1}).	First-order closure diffusivities; apply G directly.
Critical Richardson number	\mathrm{Ri}_c,\ \mathrm{Ri}_c^*	Unitless threshold.	\mathrm{Ri}_c^* denotes dynamic, state-dependent threshold (symbolic regression target).
Neutral curvature invariant	\Delta	\Delta=a_h-2a_m (unitless).	Determines sign of \mathrm{Ri}_g curvature near neutrality; sets initial bias direction in B.
Turbulent kinetic energy	TKE	0.5\left(\overline{u'^2}+\overline{v'^2}+\overline{w'^2}\right) (m^2 s^{-2}).	Used in higher-order closures for turbulence memory and intermittency; often base for K calculations.


⸻

Minor suggestions / clarifications
	•	Make symbols consistent: you use both g and g_0 in different places — pick one (I recommend g for gravity) and state its value (e.g., g=9.81\ \mathrm{m\,s^{-2}}).
	•	For \mathrm{Ri}_g and \mathrm{Ri}_b include whether \theta is potential temperature or absolute temperature; typically these use potential temperature \theta or \bar\theta (mean).
	•	Note the assumed Prandtl number (you used Pr\approx1 in derivations); put that assumption in the glossary or a separate “assumptions” section.
	•	Consider adding a short “units & constants” block (values of k, g, Pr, typical z_0 ranges) so anyone reusing the glossary has numerical context.

