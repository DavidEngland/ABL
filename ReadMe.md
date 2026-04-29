# Atmospheric Boundary Layer (ABL) — Curvature-Aware MOST Toolkit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DavidEngland/ABL/blob/main/notebooks/Curvature_Demo.ipynb)

> **Last updated:** April 2026 — See [Recent Changes](#-recent-changes) for latest additions.

**Repository for curvature-aware Monin–Obukhov Similarity Theory (MOST) diagnostics, Richardson number closures, and grid-dependent corrections for stable boundary layer parameterizations.**

---

## 📋 Overview

Coarse vertical grids in atmospheric models systematically underestimate near-surface stability in the stable boundary layer (SBL), leading to excessive turbulent mixing and warm-biased surface temperatures. This toolkit provides:

- **Analytic curvature diagnostics** for gradient Richardson number Ri_g(ζ)
- **Neutral-curvature-preserving corrections** (preserves 2Δ invariant)
- **Grid-aware damping factors** reducing coarse-grid bias by 40%+
- **Richardson number series inversion** (ζ ↔ Ri) with Newton refinement
- **Geometric vs logarithmic mean height analysis** for bulk transfer coefficients
- **Dynamic critical Richardson number** framework (Ri_c*)
- **Ultraspherical (Gegenbauer) representation** of MOST similarity operators
- **CBC/Legendre series** for power-law φ functions with equator-evaluation identity
- **Safeguarded Newton–bisection inversion** for ζ(Ri_g) across all SBL regimes

**Key Innovation:** The neutral curvature invariant Δ = α_h β_h − 2α_m β_m governs initial departure from linearity; preserving 2Δ anchors corrections to physically consistent near-neutral behavior. The momentum similarity function φ_m = (1 − b_m ζ)^{−1/4} is an ultraspherical (Gegenbauer, λ=1/4) generating-function evaluation; φ_h is Legendre (λ=1/2); their ratio identity φ_h = φ_m² when a_h=1, b_m=b_h gives an exact Clebsch–Gordan product in polynomial space.

---

## 🚀 Quick Start

### Option 1: Google Colab (No Installation)

Open the interactive demo directly in your browser:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DavidEngland/ABL/blob/main/notebooks/Curvature_Demo.ipynb)

Run the first cell to install dependencies automatically.

### Option 2: Local Installation

```bash
# Clone repository
git clone https://github.com/DavidEngland/ABL.git
cd ABL

# Create virtual environment (Python 3.10+)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Launch Jupyter
jupyter lab
# Open notebooks/Curvature_Demo.ipynb
```

### Option 3: Quick Setup Script

```bash
bash setup_dev.sh  # Creates venv, installs deps, checks for gfortran
source .venv/bin/activate
```

---

## 📂 Repository Structure

```
ABL/
├── ReadMe.md                          # This file
├── CHANGELOG.md                       # Version history
├── CONTRIBUTING.md                    # Contribution guidelines
├── requirements.txt                   # Python dependencies
├── setup_dev.sh                       # Development environment setup
├── .gitignore                         # Git ignore patterns
│
├── examples/                          # ★ Active Fortran 90 development
│   ├── README.md                      # Module index and build instructions
│   ├── module_most_profile_utils.F90  # MOST inversion utilities (safeguarded Newton)
│   ├── module_cbc_legendre_most.F90   # CBC/Gegenbauer series evaluation
│   ├── driver_cbc_gegenbauer_errors.F90  # Error-tabulation driver
│   ├── wrf_integration_example.F90    # WRF/MYNN integration sketch
│   └── overview.md                    # Jensen's Inequality / grid-bias overview
│
├── drafts/                            # Working manuscript drafts
│   ├── README.md                      # Draft index
│   ├── ultraspherical_subsection.md   # Ultraspherical MOST structure (McNider-Biazar §3)
│   ├── SBL_IBEx_expansions.md         # IBEx series framework for SBL inversions
│   ├── quad heat shear.md             # Linear φ_m / quadratic φ_h analysis
│   ├── quad heat z-less.md            # Z-less scaling & higher-degree φ_h
│   └── Momentum as a Gegenbauer ultraspherical problem.md
│
├── manuscripts/                       # Near-submission papers
│   ├── README.md                      # Manuscript status & action items
│   ├── Grid_Curvature_SBL_v01.md      # Primary paper (BLM target)
│   ├── CBC_Gegenbauer_Backward_Workflow_v01.md
│   └── WRF_Ri_Curvature_Integration_Outline_v01.md
│
├── notes/                             # Research working notes
│   ├── README.md                      # Notes index
│   ├── parameters.md / parameters2.md # MOST parameter tables
│   ├── Legendre.md                    # Legendre / Gegenbauer identities
│   ├── Central Binomial UBL heat shear.md
│   └── SBL-IBEx.md                    # IBEx reference sheet
│
├── code/                              # Python & legacy Fortran utilities
│   ├── most_similarity.md             # MOST similarity function reference
│   ├── most.f                         # Original Fortran 77 MOST code
│   ├── module_bl_mynn.F90             # WRF MYNN scheme (reference copy)
│   ├── select_f_form.py               # Stability-function selector
│   ├── profiles.py                    # Profile evaluation utilities
│   └── [additional scripts]
│
├── hw/                                # Educational / homework materials
│   ├── README.md                      # Problem set index
│   ├── CBC.md                         # Central binomial coefficient problems
│   ├── Legendre.md                    # Legendre polynomial problems
│   ├── stability functions.md         # Stability function derivation problems
│   └── [additional problem sets]
│
├── param/                             # Parameterization scaffolding (Julia/SCM)
│   ├── SCAFFOLDING.md                 # SCM parameterization design
│   └── core/                         # Core substrate/slab modules
│
├── julia/                             # Julia SCM skeleton
│   ├── SCMSkeleton.jl                 # Single-column model skeleton
│   └── SCMSkeleton_vs_SCAFFOLDING.md
│
├── notebooks/                         # Interactive demonstrations
│   └── Curvature_Demo.ipynb          # Main demo (Colab-ready)
│
├── implementations/                   # Legacy drop-in correction modules
│   ├── McNider_1DBLM_fc_module.f90   # Fortran 90 correction module
│   └── McNider_1DBLM_integration_guide.md
│
├── config/                            # Configuration files
│   └── rct_config.yaml               # RCT run configuration
│
├── data/                              # Observational / validation data
├── visuals/                           # Figures and plots
└── refs/ references/                  # Reference materials and literature
```

---

## 🎯 Key Features

### 1. Curvature Diagnostics

Compute second derivative of gradient Richardson number:

```python
# Analytic curvature for power-law φ functions
def curvature_linear(zeta, a_m, a_h):
    V_log = a_h/(1.0 + a_h*zeta) - 2.0*a_m/(1.0 + a_m*zeta)
    W_log = (a_h**2)/((1.0 + a_h*zeta)**2) - 2.0*(a_m**2)/((1.0 + a_m*zeta)**2)
    F = (1 + a_h*zeta) / (1 + a_m*zeta)**2
    return F * (2.0 * V_log + zeta * (V_log**2 - W_log))

# Neutral curvature invariant
Delta = a_h - 2*a_m
neutral_curvature = 2 * Delta
```

### 2. Grid-Aware Correction

Preserve neutral curvature while reducing coarse-grid bias:

```python
# Exponential damping factor
def G_correction(zeta, dz, D=0.8, dz_ref=10.0, zeta_ref=0.5, p=1.0, q=2.0):
    exponent = -D * (dz/dz_ref)**p * (zeta/zeta_ref)**q
    return np.exp(exponent)

# Apply to diffusivities
K_m_corrected = K_m_original * G_correction(zeta, dz)
K_h_corrected = K_h_original * G_correction(zeta, dz)
```

### 3. Richardson Number Inversion

Fast ζ(Ri) inversion using series + Newton:

```python
# Series seed (O(Ri³) accurate)
def zeta_from_ri_series(Ri, Delta, c1):
    return Ri - Delta*Ri**2 + (1.5*Delta**2 - 0.5*c1)*Ri**3

# Newton refinement (1-2 iterations to machine precision)
def zeta_from_ri_newton(Ri_target, phi_m, phi_h, z0, tol=1e-10, maxit=20):
    z = zeta_from_ri_series(Ri_target, Delta, c1)  # seed
    for _ in range(maxit):
        F = phi_h(z) / phi_m(z)**2
        Ri = z * F
        if abs(Ri - Ri_target) < tol:
            break
        dRi_dz = F + z * (...)  # derivative
        z -= (Ri - Ri_target) / dRi_dz
    return z
```

### 4. Geometric Mean Heights

```python
# Representative heights for layer [z0, z1]
z_g = np.sqrt(z0 * z1)           # Geometric mean (log-space midpoint)
z_L = (z1 - z0) / np.log(z1/z0)  # Logarithmic mean (exact for shear)
z_a = 0.5 * (z0 + z1)            # Arithmetic mean (biases high)
```

---

## 📊 Validation Results

| Metric | Target (Δz=100m) | Baseline | With Correction |
|--------|------------------|----------|-----------------|
| Bias ratio B | < 1.2 | ~1.8 | ~1.15 |
| Surface flux RMSE | < 15% | ~30% | ~12% |
| Inversion height error | < 20 m | ~50 m | ~18 m |
| Neutral curvature preservation | \|2Δ*−2Δ\|/\|2Δ\| < 5% | N/A | 2.3% |
| Computational overhead | < 5% | 0% | 2.8% |

**Validation Cases:**
- GABLS1 LES (9-hour nocturnal evolution)
- ARM NSA Alaska (persistent stable nights)
- SHEBA Arctic winter inversions

---

## 🔄 Recent Changes

### April 2026

**Fortran modules (`examples/`)**
- `module_most_profile_utils.F90`: Added `zeta_from_rig_safeguarded()` — safeguarded Newton + bisection fallback with branch-aware brackets; `fm_fh_from_rig` now calls this as primary solver.
- `module_cbc_legendre_most.F90`: CBC/Gegenbauer series evaluation module for power-law φ functions; implements Legendre equator identity P_{2n}(0) = (−1)^n C(2n,n)/4^n.
- `driver_cbc_gegenbauer_errors.F90`: Error-tabulation driver comparing exact φ_m vs truncated Gegenbauer series at equator.
- `wrf_integration_example.F90`: WRF/MYNN-style integration sketch using safeguarded ζ(Ri_g) inversion.

**Draft manuscripts (`drafts/`)**
- `ultraspherical_subsection.md`: New §3.x subsection for McNider-Biazar paper. Covers Gegenbauer generating functions (λ=1/4 momentum, λ=1/2 heat), Sturm-Liouville operator in stability space, Clebsch-Gordan squaring identity φ_h=φ_m², Ri_g mapping, asymptotic 3-term expansions, dynamic Ri_c, and algebra-first inversion decision tree with Mermaid flowchart.
- `SBL_IBEx_expansions.md`: IBEx (Internal/Boundary/External) series framework for SBL Richardson mappings; includes ζ/(1+βζ) canonical forms and SBL asymptote analysis.
- `quad heat shear.md`: Analysis of linear-φ_m / quadratic-φ_h mixed MOST form and physical implications.
- `quad heat z-less.md`: Connection between z-less scaling and the required degree of φ_h.

**Working notes (`notes/`)**
- `parameters.md`, `parameters2.md`: Consolidated MOST parameter tables (Businger, Dyer, Högström, McNider).
- `Legendre.md`: Legendre/Gegenbauer orthogonality and generating-function identities.
- `Central Binomial UBL heat shear.md`: CBC series for UBL heat-shear asymptotic expansion.
- `SBL-IBEx.md`: Quick-reference sheet for IBEx series framework.

**Manuscripts (`manuscripts/`)**
- `WRF_Ri_Curvature_Integration_Outline_v01.md`: Outline for WRF integration paper.
- `CBC_Gegenbauer_Backward_Workflow_v01.md`: Backward workflow connecting CBC identities to MOST series.

---

## 📚 Documentation

### Technical Papers & Derivations

- **[Vertical Resolution.md](Vertical%20Resolution.md)** — Grid correction framework overview
- **[McNider_Ri_Corrections_Overview.md](McNider_Ri_Corrections_Overview.md)** — Comprehensive implementation guide
- **[curvature.md](curvature.md)** — Full mathematical derivation of ∂²Ri_g/∂ζ²
- **[SBL corrections.md](SBL%20corrections.md)** — Stable boundary layer implementation
- **[ultraspherical_subsection.md](drafts/ultraspherical_subsection.md)** — Ultraspherical structure of MOST operators
- **[SBL_IBEx_expansions.md](drafts/SBL_IBEx_expansions.md)** — IBEx series for SBL Richardson mappings
- **[CANONICAL_GLOSSARY.md](CANONICAL_GLOSSARY.md)** — Symbol and notation reference

### Integration Guides

- **[examples/README.md](examples/README.md)** — Fortran 90 module index and build instructions
- **[McNider_1DBLM_integration_guide.md](implementations/McNider_1DBLM_integration_guide.md)** — Legacy Fortran 77/90 drop-in instructions
- **[IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md)** — General implementation guide

### Educational Materials

- **[hw/README.md](hw/README.md)** — Problem set index
- **[CBC.md](hw/CBC.md)** — Central binomial coefficient problem set
- **[Legendre.md](hw/Legendre.md)** — Legendre/Gegenbauer problem set
- **[intro.md](intro.md)** — Guest lecture outline (60-90 min)

---

## 🔬 Research Applications

### Current Focus

1. **Arctic Stable Boundary Layers** — Reducing warm bias in polar climate models
2. **Urban Air Quality** — Improving nocturnal O₃/PM2.5 forecasts via better vertical mixing
3. **Low-Level Jets** — Capturing LLJ onset timing and core height
4. **Slope Flows** — Extending to complex terrain (McNider specialty)

### Collaboration

**Principal Investigators:**
- Richard T. McNider (University of Alabama in Huntsville)
- Arastoo P. Biazar (University of Alabama in Huntsville)
- David E. England (Lead Developer)

See [emails/McNider_Biazar_Status_2025.md](emails/McNider_Biazar_Status_2025.md) for latest collaboration updates.

---

## 🛠️ Development

### Running Tests

```bash
# Open test notebook
jupyter lab notebooks/Curvature_Demo.ipynb

# Run key cells
# - Core functions (cell 5)
# - Demo plots (cell 6-8)
# - Synthetic NetCDF (cell 9)
```

### Pre-Commit Checklist

- [ ] Notebooks execute without errors
- [ ] Requirements install in fresh venv
- [ ] Key examples produce expected plots
- [ ] Code follows PEP 8 (optional: `ruff check .`)
- [ ] Documentation updated for new features

### Fortran Development

```bash
# Build examples (from examples/ directory)
cd examples

# Compile Fortran 90 modules
gfortran -c module_most_profile_utils.F90
gfortran -c module_cbc_legendre_most.F90
gfortran -c driver_cbc_gegenbauer_errors.F90 -o driver_cbc_gegenbauer_errors

# Run CBC/Gegenbauer error table
./driver_cbc_gegenbauer_errors

# Legacy correction module
gfortran -c ../implementations/McNider_1DBLM_fc_module.f90
```

---

## 📖 Citation

If you use this toolkit in your research, please cite:

```bibtex
@misc{england2025abl,
  author = {England, David E. and McNider, Richard T. and Biazar, Arastoo P.},
  title = {Curvature-Aware MOST Toolkit for Stable Boundary Layer Corrections},
  year = {2025},
    year = {2026},
  publisher = {GitHub},
  url = {https://github.com/DavidEngland/ABL}
}
```

**Related Publications:**
- England & McNider (1995). "Stability Functions Based Upon Shear Functions." *Boundary-Layer Meteorology*.
- McNider et al. (1995). [Additional key papers — see [refs/](refs/)]
- [JAMC Grid Dependence manuscript](JAMC_format_Grid_Dependence_V15_10-31-2025.md) — Grid-dependent corrections paper (under review)

---

## 🤝 Contributing

We welcome contributions! To contribute:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feat/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feat/amazing-feature`)
5. **Open** a Pull Request

Please ensure:
- New features include documentation
- Notebooks remain Colab-compatible
- Code follows existing style conventions
- Tests pass (run demo notebook)

---

## 📧 Contact

**David E. England**  
Email: David.England@UAH.Edu  
GitHub: [@DavidEngland](https://github.com/DavidEngland)

**Project Links:**
- Repository: https://github.com/DavidEngland/ABL
- Issues: https://github.com/DavidEngland/ABL/issues
- Discussions: https://github.com/DavidEngland/ABL/discussions

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **Funding:** NSF Atmospheric & Geospace Sciences, DOE Atmospheric System Research
- **Data:** ARM Climate Research Facility, SHEBA, GABLS LES intercomparison
- **Computing:** UAH ESSC cluster, NSF XSEDE allocation
- **Inspiration:** Businger et al. (1971), Beljaars & Holtslag (1991), Cuxart et al. (2006)

---

## 🔗 Quick Links

- [🎓 Graduate Lecture](intro.md)
- [📊 Interactive Demo](https://colab.research.google.com/github/DavidEngland/ABL/blob/main/notebooks/Curvature_Demo.ipynb)
- [📘 Comprehensive Guide](docs/McNider_Ri_Corrections_Overview.md)
- [⚙️ Fortran Integration](implementations/McNider_1DBLM_integration_guide.md)
- [📝 Publication Topics](topics.md)

---

**Last Updated:** January 2025  
**Version:** 1.0.0