Here is a consolidated, research‑grade synthesis of the parameter sets, regimes, and specialized stable‑layer formulations you listed—organized to make the physics, numerical‑model relevance, and cross‑study consistency immediately clear.
(All citations refer to the sources you summarized; no new external claims are introduced.)

---

1. Canonical Parameter Sets for Businger–Dyer Type Similarity Functions

The table you provided captures the core empirical constants used in MOST‑based surface‑layer schemes. These constants define:

• Unstable branch:\phi_m = (1 - \gamma_m \zeta)^{-1/4}, \qquad
\phi_h = a_h^{-1}(1 - \gamma_h \zeta)^{-1/2}

• Stable branch:\phi_m = 1 + \beta_m \zeta, \qquad
\phi_h = a_h^{-1} + \beta_h \zeta



Interpretation of the constants

• \(k\) (von Kármán constant): Controls the slope of the log‑law; values range 0.35–0.41, with modern consensus near 0.40.
• \(a_h^{-1}\): Neutral limit of the temperature gradient function; typically 0.95–1.35, depending on whether the source reports \(a_h\) or its inverse.
• \(\gamma_m, \gamma_h\): Strength of the unstable curvature; larger values steepen the gradient as \(\zeta \to -1/\gamma\).
• \(\beta_m, \beta_h\): Linear slope in stable conditions; classical values cluster around 5.


Cross‑study consistency

• Businger (1971) and Dyer (1974) remain the “Kansas standard,” widely used in NWP and LES surface‑layer schemes.
• Högström (1996) provides slightly larger unstable coefficients (\(\gamma_m=19\)), consistent with later field campaigns.
• Brutsaert (1982) and Paulson (1972) essentially reproduce the Dyer constants.


These constants are still embedded in WRF, ECMWF IFS, RUC, ARPS, and many land‑surface models.

---

2. Validity Domains and Regime Structure

Unstable regimes (`\(\zeta < 0\)`)

• Near‑neutral: \(-1 < \zeta < 0\).
Kansas‑type functions perform extremely well here.
• Strongly unstable: \(\zeta < -1\).
Classical forms begin to deviate; free‑convection scaling approaches\phi_h \sim (-\zeta)^{-1/3}.

• Formal validity: Many classical functions are cited only for \(-2 < \zeta < 0\).


Stable regimes (`\(\zeta > 0\)`)

• Moderately stable: \(0 < \zeta \le 1\).
Linear forms (\(\phi = 1 + 5\zeta\)) remain adequate.
• Very/extremely stable: \(\zeta > 1\).
Linear forms fail; turbulence becomes intermittent, non‑stationary, and “z‑less.”


Richardson‑number regimes

• Nearly neutral: \(0 < Ri < 0.02\)
• Weakly stable: \(0.02 < Ri < 0.12\)
• Very stable: \(0.12 < Ri < 0.7\)
• Extremely stable: \(Ri > 0.7\)
• Critical limit: \(Ri_c \approx 0.21–0.25\) (onset of laminarization)


These Ri thresholds are crucial for turbulence switches in higher‑order closures (e.g., MYNN, QNSE, EFB).

---

3. Specialized Stable‑Layer Formulations

These schemes were developed to avoid the decoupling problem (runaway cooling, vanishing fluxes) that plagues linear MOST functions at high stability.

Beljaars–Holtslag (1991)

• Uses exponential decay terms with constants \(a=1, b=0.667, c=5, d=0.35\).
• Smooths the transition into very stable regimes.
• Adopted in ECMWF and many land‑surface models.


Cheng & Brutsaert (2005)

• Momentum: \(a=6.1, b=2.5\)
• Heat: \(c=5.3, d=1.1\)
• Extends validity to \(\zeta \approx 20\).
• Adopted in MYNN and several Arctic‑focused schemes.


SHEBA (2007)

• Derived from Arctic sea‑ice data up to \(\zeta \approx 100\).
• Momentum: \(a_m=5, b_m=0.3\)
• Heat: \(a_h=5, b_h=0.4\)
• Essential for polar modeling where classical MOST fails.


---

4. How to Use These Parameters in Numerical Models

1. Choose a base similarity family (Businger–Dyer, Beljaars–Holtslag, Cheng–Brutsaert, SHEBA).
2. Select constants appropriate for the regime and surface type.
3. Ensure continuity across neutral, unstable, and stable branches.
4. Apply Ri‑based turbulence switches for very stable conditions.
5. Validate against site‑specific data, especially in Arctic or nocturnal SBL cases.


---

If you want, I can assemble a unified parameter table, derive the corresponding flux–profile integrals, or produce a regime‑aware similarity function library (e.g., in Python, Fortran, or Julia) tailored to your Arctic inversion experiment.