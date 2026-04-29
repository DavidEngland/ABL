# drafts/ — Working Manuscript Drafts

In-progress technical writing for the McNider-Biazar collaboration and related papers.  
For near-submission manuscripts, see [manuscripts/](../manuscripts/).

---

## Active Drafts

### Ultraspherical / Gegenbauer Theory

| File | Description | Status |
|---|---|---|
| [ultraspherical_subsection.md](ultraspherical_subsection.md) | §3.x subsection for McNider-Biazar paper: Gegenbauer generating functions, SL operator in stability space, Clebsch-Gordan squaring identity, Ri_g mapping, asymptotic expansions, dynamic Ri_c, algebra-first inversion policy | Active — 9 subsections |
| [Momentum as a Gegenbauer ultraspherical problem.md](Momentum%20as%20a%20Gegenbauer%20ultraspherical%20problem.md) | Source derivation notes | Reference |
| [Gegebauer mods.md](Gegebauer%20mods.md) | Modification notes for Gegenbauer approach | Working |
| [SBL_IBEx_expansions.md](SBL_IBEx_expansions.md) | IBEx (Internal/Boundary/External) series framework for SBL ζ/(1+βζ) Richardson mappings | Active |

### SBL Alternatives

| File | Description | Status |
|---|---|---|
| [quad heat shear.md](quad%20heat%20shear.md) | Linear φ_m / quadratic φ_h mixed form: profiles, Pr_t variation, Ri_g behavior, decoupling | Active |
| [quad heat z-less.md](quad%20heat%20z-less.md) | Z-less scaling and the degree constraint on φ_h | Active |

### UBL / Unstable Regime

| File | Description | Status |
|---|---|---|
| [UBL Shear-Similarity Functions.md](UBL%20Shear-Similarity%20Functions.md) | Unstable boundary layer shear functions | Active |
| [UBL heat–shear asymptotic.md](UBL%20heat%E2%80%93shear%20asymptotic.md) | Asymptotic analysis for UBL heat-shear coupling | Active |
| [UBLclaudePlan.md](UBLclaudePlan.md) | UBL extension plan | Planning |

### Grid Curvature Corrections

| File | Description |
|---|---|
| [Grid-Curvature Correction SBL.md](Grid-Curvature%20Correction%20SBL.md) | Core grid-curvature correction derivation |
| [Grid-Curvature_Correction_MASTER.tex](Grid-Curvature_Correction_MASTER.tex) | Master LaTeX version |
| [Curvature_Correction_Draft.md](Curvature_Correction_Draft.md) | Draft curvature correction |
| [Grid_Curvature_Bias_Summary.md](Grid_Curvature_Bias_Summary.md) | Summary of grid-curvature bias results |
| [Vertical grid spacing ($\Delta z$) Bias.md](Vertical%20grid%20spacing%20%28%24%5CDelta%20z%24%29%20Bias.md) | Grid spacing bias analysis |

### Other Working Documents

| File | Description |
|---|---|
| [summary.md](summary.md) | High-level project summary |
| [numerics.md](numerics.md) | Numerical methods notes |
| [asymptotics.txt](asymptotics.txt) | Asymptotic expansion scratch work |
| [log-expansions.md](log-expansions.md) | Log-profile series expansions |
| [gradient Richardson number.md](gradient%20Richardson%20number.md) | Ri_g reference derivation |
| [Adaptive_Regime_Transitions_Draft.md](Adaptive_Regime_Transitions_Draft.md) | Adaptive regime transition framework |
| [Stability Parameter Tipping Error.md](Stability%20Parameter%20Tipping%20Error.md) | Tipping-point error analysis |
| [Generalized Exponential Form.md](Generalized%20Exponential%20Form.md) | Exponential stability function forms (BH-type) |

---

## Key Mathematical Results (quick reference)

**Ultraspherical identity:**
φ_m = (1 − b_m ζ)^{−1/4} generates Gegenbauer polynomials C_n^{(1/4)};  
φ_h = (1 − b_h ζ)^{−1/2} generates Legendre polynomials P_n = C_n^{(1/2)}.  
When a_h = 1, b_m = b_h: φ_h = φ_m² exactly (Clebsch-Gordan convolution in polynomial space).

**Ri_g inversion degrees:**
- Linear φ_m, φ_h → linear solve, finite Ri_c = 1/β
- Linear φ_m, linear φ_h (unequal slopes) → quadratic in ζ
- Linear φ_m, quadratic φ_h → cubic in ζ, no finite Ri_c
- Power-law φ_m → cubic via u = √φ_m substitution

**IBEx series (SBL canonical form):**  
ζ/(1 + βζ) = Σ (−β)^{n−1} ζ^n for |βζ| < 1 (internal);  
= 1/β · Σ (−1)^{n−1} (βζ)^{−n} for |βζ| > 1 (external, z-less asymptote ζ → Ri_c).
