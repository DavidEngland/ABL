# McNider--Biazar Collaboration Action Items

## 🔵 McNider-Led Sections (Theoretical Physics & Historical Context)

### 1. Grid-Invariance Principle Expansion (Section 3.1--3.2)

**Current Status:** Basic ODE derivation present, lacks historical motivation.

**Required Additions:**
- [ ] **Historical narrative:** Why did the 1995 Kansas/Alabama urban simulations first reveal grid sensitivity? What specific forecast failures motivated the investigation?
- [ ] **Mathematical rigor:** Full chain-rule derivation from Eq. 13 to Eq. 16, showing how $\partial \ln \Rig / \partial \ln \Delta z$ emerges from the MOST structure.
- [ ] **Approximation justification:** Under what conditions is the power-law form (Eq. 17) superior to alternative formulations (polynomial, rational)? Show convergence analysis.
- [ ] **Variable $L(z)$ treatment:** Current draft assumes constant $L$ over $\Delta z$. Derive the omission metric $E_{\text{omit}}$ rigorously and specify when layer-splitting is mandatory.

**Deliverables:**
- Expanded Section 3 (target: 3 pages → 5 pages)
- New Figure: "Evolution of grid-invariance concept 1995--2025"
- Supplementary Material: Derivation of $f_c$ for non-log-linear $\phi$ functions

---

### 2. Curvature Analysis Deep Dive (Section 2.2--2.3)

**Current Status:** Neutral invariant $\Delta$ introduced, but logarithmic sensitivity formalism ($V_{\log}$, $W_{\log}$) not fully explained.

**Required Additions:**
- [ ] **Logarithmic sensitivity derivation:** Show full expansion of Eq. 11 from first principles, connecting to the "compact curvature formula" from earlier McNider work.
- [ ] **Inflection point handling:** Derive conditions under which $\Rig(\zeta)$ can exhibit inflection (sign change in curvature). Propose adaptive layer-splitting algorithm.
- [ ] **Generalization to non-log-linear forms:** Extend $\Delta$ invariant to Beljaars--Holtslag, QNSE, Gryanik--Lüpkes formulations. Create comparison table.
- [ ] **Physical interpretation:** Why is $\Delta < 0$ universal in SBL parameterizations? Connect to turbulent Prandtl number evolution.

**Deliverables:**
- Expanded Section 2 (target: 2 pages → 4 pages)
- New Table: "$\Delta$ values for 10 canonical MOST formulations"
- Appendix A: "Full curvature derivation for arbitrary $\phi_m, \phi_h$"

---

### 3. Q-SBL Surrogate Scheme (Section 5.3)

**Current Status:** Mentioned as Option 3, no mathematical detail.

**Required Additions:**
- [ ] **Quadratic surrogate construction:** Derive $\phi_m^{\text{Q}}(\zeta) = 1 + a_m \zeta + b_m \zeta^2$ coefficients that achieve $\Delta^{\text{Q}} = 0$ while matching observations.
- [ ] **Blending strategy:** How to transition from Q-SBL (near-surface) to log-linear (aloft) without introducing discontinuities? Propose smooth transition function.
- [ ] **Validation plan:** What observational datasets are needed to tune $a_m, b_m$? Suggest SHEBA, ARM, or FLUXNET-specific campaigns.

**Deliverables:**
- New Section 5.4: "Zero-Curvature Surrogate Design" (2 pages)
- Implementation pseudocode for WRF/CESM2
- Jupyter notebook demonstrating Q-SBL vs. log-linear comparison

---

### 4. Unstable Branch Investigation (Section 8.2)

**Current Status:** Open question posed, no preliminary analysis.

**Required Additions:**
- [ ] **Businger--Dyer curvature analysis:** Compute $\Delta_{\text{USL}}$ for canonical unstable parameters. Is it truly zero, or slightly positive?
- [ ] **Entrainment zone bias:** Does the CBL top (capping inversion) exhibit analogous grid sensitivity? Preliminary LES analysis.
- [ ] **Hypothesis:** If $\Delta_{\text{USL}} \approx 0$, this explains why daytime forecasts are less grid-sensitive than nocturnal ones.

**Deliverables:**
- New Section 8.2: "Convective BL Curvature: Preliminary Analysis" (1.5 pages)
- Figure: "$\Rig(\zeta)$ curvature for unstable vs. stable regimes"

---

## 🟠 Biazar-Led Sections (Numerical Methods & Closure Equations)

### 5. ADM Mathematical Foundation (Section 4.2--4.3)

**Current Status:** ADM introduced conceptually, but detailed mathematical proofs absent.

**Required Additions:**
- [ ] **Convergence theorem:** Prove that the ADM series $e = \sum_{n=0}^\infty e_n$ converges for typical SBL parameter ranges ($S, N, \ell, C_\epsilon$).
- [ ] **Adomian polynomial derivation:** Show explicit calculation of $A_0, A_1, A_2$ for the $e^{3/2}$ dissipation term. Extend to $e^2$ (Kolmogorov-type) if needed.
- [ ] **Comparison with Newton--Raphson:** Demonstrate failure modes of iterative methods in near-laminar limit ($e \to 0$). Show ADM stability.
- [ ] **Computational cost analysis:** FLOPS count for ADM vs. iteration. Show break-even point.

**Deliverables:**
- Expanded Section 4 (target: 2 pages → 4 pages)
- Appendix B: "Adomian Polynomial Recursion Formulas"
- Algorithm Box: "ADM TKE Solver (Fortran/Python pseudocode)"

---

### 6. TKE Closure Integration (Section 4.4)

**Current Status:** Conceptual link to McNider correction established, but implementation pathway unclear.

**Required Additions:**
- [ ] **End-to-end workflow:** Step-by-step procedure for computing $K_{\text{final}}$ in a prognostic TKE model:
  1. ADM solves for $e(z)$
  2. Compute $K_{\text{base}} = C_K \ell \sqrt{e}$
  3. Diagnose $\Rig, B$ from ADM-derived gradients
  4. Apply McNider $f_c$ correction
- [ ] **Consistency checks:** How to ensure $K_{\text{base}}$ from ADM is consistent with the MOST-diagnosed $\Rib$? Propose iterative refinement if needed.
- [ ] **Boundary condition handling:** Surface flux BC vs. Dirichlet BC for TKE. Which is more stable with ADM?

**Deliverables:**
- New Section 4.5: "ADM--McNider Integration Protocol" (2 pages)
- Flowchart: "TKE closure with ADM + curvature correction"
- Test case: GABLS1 with full TKE budget (compare vs. algebraic closure)

---

### 7. Parameter Sensitivity and Inverse Estimation (Section 4.6, new)

**Current Status:** Not addressed.

**Required Additions:**
- [ ] **Variational approach:** Use ADM to compute $\partial e / \partial C_\epsilon$, $\partial e / \partial \ell$ analytically. Feed into data assimilation framework to estimate optimal closure coefficients from tower observations.
- [ ] **Machine learning integration:** Train PINN to predict $C_\epsilon(\zeta, \Delta z)$ rather than using fixed constants. Does this improve generalization across sites?
- [ ] **Uncertainty quantification:** Use ADM series truncation error as a proxy for closure uncertainty. Propagate through ensemble forecasting.

**Deliverables:**
- New Section 4.6: "Data-Driven Parameter Estimation via ADM" (2 pages)
- Jupyter notebook: "Inverse estimation of $C_\epsilon$ from ARM NSA data"
- Comparison: Fixed vs. adaptive $C_\epsilon$ in WRF runs

---

### 8. Advanced ADM Applications (Section 8.3, future work)

**Current Status:** Mentioned as open question.

**Required Additions:**
- [ ] **Time-dependent TKE:** Can ADM handle the transient term $\partial e / \partial t$ in Eq. 20? Propose semi-implicit scheme.
- [ ] **Mellor--Yamada Level 2.5:** ADM for algebraic second-moment closures. Show how to decompose $\overline{u'w'}$, $\overline{w'\theta'}$ nonlinearities.
- [ ] **Coupled land--atmosphere:** Use ADM to solve joint soil heat diffusion + atmospheric TKE system. Application: permafrost modeling.

**Deliverables:**
- New Section 8.3: "ADM for Higher-Order Closures" (1.5 pages)
- Proof-of-concept: Time-dependent TKE with ADM in idealized 1D column

---

## 🟢 Joint Sections (Collaborative Synthesis)

### 9. Validation Strategy Expansion (Section 6)

**Current Lead:** England (LES benchmarks), but needs McNider (field campaign interpretation) and Biazar (closure diagnostics).

**Required Additions:**
- [ ] **McNider:** Interpret MOSAiC results in context of Arctic climate biases. Why does correction fail at sea ice leads? Propose heterogeneity parameterization.
- [ ] **Biazar:** Add diagnostic plots showing ADM-computed $e(z)$ vs. iterative solution for GABLS1. Demonstrate convergence advantage.
- [ ] **Joint:** Taylor diagrams showing skill improvement across multiple metrics (T, H, mixing ratio, TKE) simultaneously.

**Deliverables:**
- Enhanced Section 6 (target: 3 pages → 6 pages)
- 5 new figures (MOSAiC time series, ADM convergence, Taylor diagrams)
- Supplementary Material: Full GABLS suite results (GABLS2--4)

---

### 10. WRF/CESM2 Implementation Guide (Section 5.1--5.2)

**Current Lead:** England (code snippets), but needs McNider (MOST module modifications) and Biazar (TKE scheme integration).

**Required Additions:**
- [ ] **McNider:** Detailed WRF `module_sf_mynn.F` modifications. Show exact line numbers where $f_c$ is inserted.
- [ ] **Biazar:** CESM2 CAM7 TKE closure (`eddy_diff.F90`) integration. How to replace iterative solver with ADM?
- [ ] **Joint:** Namelist parameter recommendations. Provide default values for $D, p, q, \Delta z_{\text{ref}}$ based on validation campaign.

**Deliverables:**
- New Appendix C: "WRF Implementation Guide" (4 pages)
- New Appendix D: "CESM2 Implementation Guide" (4 pages)
- GitHub pull requests to WRF/CESM2 repositories with reference implementation

---

## 📊 Figure/Table Assignments

| Figure # | Lead | Description | Deadline |
|----------|------|-------------|----------|
| 1 | McNider | $\Rig(\zeta)$ for $\Delta = -2.2, 0, +2.2$ with Jensen gap | Week 1 |
| 2 | McNider | Cascade diagram (bias → diffusivity → flux → temperature) | Week 1 |
| 3 | McNider | Contour plot of $f_c(\Delta z, \zeta)$ | Week 2 |
| 4 | England | MOSAiC validation scatter plots (before/after) | Week 2 |
| 5 | Biazar | ADM convergence for GABLS1 TKE profile | Week 3 |
| 6 | Biazar | SHAP feature importance (PINN analysis) | Week 3 |
| 7 | Joint | Taylor diagram (WRF skill across 10 cases) | Week 4 |

| Table # | Lead | Description | Deadline |
|---------|------|-------------|----------|
| 1 | England | GABLS1 performance metrics | Week 1 |
| 2 | McNider | $\Delta$ values for canonical MOST formulations | Week 2 |
| 3 | Biazar | ADM vs. Newton--Raphson computational cost | Week 3 |
| 4 | Joint | WRF/CESM2 namelist parameter recommendations | Week 4 |

---

## 🎯 Publication Timeline

| Milestone | Lead | Deadline |
|-----------|------|----------|
| **Draft v2.0** (with expanded Sections 2--4) | McNider + Biazar | Feb 15, 2025 |
| **Full figures/tables** | All | Feb 28, 2025 |
| **Internal review** (UAH colleagues) | England | Mar 10, 2025 |
| **Submit to JAMC** | England | Mar 20, 2025 |
| **Expected reviews** | — | May 2025 |
| **Resubmit** | All | Jun 2025 |
| **Publication** | — | Aug 2025 |

---

## 💡 Key Discussion Points for Next Meeting

1. **Q-SBL vs. Correction:** Should we pursue the zero-curvature surrogate (Section 5.3) as the primary recommendation, or position it as a "future alternative" to the $f_c$ correction? McNider's call.

2. **ADM in Operational Models:** Is it realistic to replace Newton--Raphson with ADM in WRF 4.6/CESM3? Or should ADM remain a "research tool" for offline diagnosis? Biazar's assessment.

3. **Authorship Order:** Proposal: **England (first), McNider (second), Biazar (third)**. Rationale: England did initial synthesis and validation, McNider provides theoretical foundation, Biazar contributes numerical methods. Discuss.

4. **Follow-Up Papers:**
   - **Paper 2 (McNider lead):** Q-SBL surrogate scheme with global validation
   - **Paper 3 (Biazar lead):** ADM for Mellor--Yamada and second-moment closures
   - **Paper 4 (England lead):** PINN surrogate deployment in NOAA/ECMWF operational systems

---

## 📝 Notes from Previous Discussions

- McNider: "The neutral invariant $\Delta$ is the key insight—it's what makes the correction physics-based rather than empirical."
- Biazar: "ADM gives us a way to avoid the iterative hell that plagues very stable cases. It's not just faster—it's more robust."
- England: "We need to emphasize that this is **not** a new turbulence scheme. It's a resolution-aware wrapper for existing schemes."

---

**Last Updated:** January 20, 2025  
**Next Meeting:** January 27, 2025, 2:00 PM CST (Zoom)
