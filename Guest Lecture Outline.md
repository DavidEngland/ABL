Guest Lecture Outline — MOST, Stable BL, and Curvature-Aware Corrections

Audience: Graduate students & scientists (STEM)
Goal: Intuition → Equations → Implications → Open Problems
Duration: 60–90 min

---

🕐 Time Plan

• 0–5 min: Motivation & objectives
• 5–20 min: MOST primer (intuitive → formal)
• 20–35 min: Gradient Richardson number & curvature
• 35–50 min: Grid sensitivity & correction
• 50–65 min: Results & validation
• 65–80 min: Open questions & discussion
• 80–90 min: Live demo (optional) + Q&A


---

🎯 Learning Objectives

• Understand similarity variables & Obukhov length \(L\).
• See how \(Ri_g(\zeta)\) curvature explains coarse‑grid bias.
• Learn neutral‑curvature invariance (2Δ) and why it matters.
• Review correction strategies & validation gaps.


---

1. Motivation

Hook:
“Forecasting Arctic temperature or urban air quality depends on the first 100 m above ground. That’s where models fail hardest — and where we can fix them.”

Problem:

• Stable BLs are grid‑sensitive → coarse Δz misclassifies stability.
• Models over‑mix at night → warm bias, wrong inversion height, pollutant errors.


Goal:

• Resolution‑aware fix that preserves neutral physics, reduces bias 40%+, and scales across applications.


---

2. MOST Primer

Anchor:
“Near the surface, turbulence collapses into one number: \(\zeta = z/L\).”

Key Definitions:

• Obukhov length: \(L=-\frac{u_*^3\theta}{\kappa g\,\overline{w^\prime\theta^\prime}}\)
• Dimensionless height: \(\zeta=z/L\)
• Stability functions: \(\phi_m, \phi_h\)


Stability Trinity Table:

Quantity	Formula	Meaning
Obukhov Length	\(L\)	Shear ≈ buoyancy height
Dimensionless Height	\(\zeta=z/L\)	“How many L up?”
Richardson Number	\(Ri_g=\zeta \phi_h/\phi_m^2\)	Local stability


---

3. Curvature of `\(Ri_g\)`

Question: How fast does stability grow with height?

Invariant:

\Delta = \alpha_h \beta_h - 2\alpha_m \beta_m


• Δ < 0 → concave down (typical SBL)
• Δ = 0 → linear (neutral)
• Δ > 0 → concave up (rare)


Impact:

• Coarse grids average over curvature → underestimate stability → over‑mixing.


---

4. Grid Sensitivity & Correction

Villain: Layer averaging.

• Jensen’s inequality: bulk Ri < local Ri → bias.
• Geometric mean height \(z_g=\sqrt{z_0 z_1}\) is natural center.


Fix Strategy:

1. Preserve 2Δ (neutral anchor).
2. Dampen tail for ζ > 0.
3. Converge to 1 as Δz → 0.


---

5. Results Snapshot

• ✅ Neutral curvature preserved within 2%.
• ✅ Bias reduction 40–55% at Δz = 60–100 m.
• ✅ Cost <3% overhead.
• ✅ Robust across 100+ cases.


Validation Example:
ARM NSA tower → corrected scheme reduces RMSE by 28%, inversion error −15 m → −5 m.

---

6. Dynamic Critical Richardson Number

Problem: Fixed \(Ri_c=0.25\) misses intermittent turbulence.

Solution:

Ri_c^* = 0.25 + \alpha \frac{\Gamma}{\Gamma_{ref}} + \beta(1 - TKE_{rel})


Paths:

• Mixing length (McNider).
• Diffusivity (Biazar).


---

7. Open Questions

• Terrain curvature vs Ri_g curvature.
• Urban BL: variable L(x,y,z).
• Polar BL: sea‑ice leads, roughness.
• Planetary BLs: Mars, Titan, exoplanets.


---

8. Closing

Three Takeaways:

1. Curvature matters (Δ sets departure from neutrality).
2. Use geometric mean height \(z_g\).
3. Preserve 2Δ in corrections.


Final Thought:
“The atmosphere doesn’t care about our grid spacing — but our forecasts do. Respect curvature, honor the physics.”

---
