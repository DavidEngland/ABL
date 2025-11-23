# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-01-XX

### Added
- Initial public release
- Curvature diagnostics for gradient Richardson number
- Neutral-curvature-preserving grid correction framework
- Richardson number series inversion (ζ ↔ Ri)
- Geometric vs logarithmic mean height analysis
- Interactive Jupyter notebook (Colab-ready)
- Fortran 90 correction module for McNider 1DBLM
- Comprehensive documentation:
  - Vertical Resolution.md
  - McNider_Ri_Corrections_Overview.md
  - curvature.md (full derivation)
  - SBL corrections.md
- Educational materials (homework problems)
- Integration guides for operational models

### Documentation
- Main README with quick start and feature overview
- CONTRIBUTING.md with contribution guidelines
- LICENSE (MIT)
- setup_dev.sh for automated environment setup

### Validation
- GABLS1 LES benchmark results
- ARM NSA tower validation
- SHEBA Arctic winter cases
- Dallas/Ft. Worth urban tower (Remote Sensing 2024)

---

## [Unreleased]

### Planned
- Python package structure (pip installable)
- Additional validation cases (CASES-99, Cabauw)
- Machine learning surrogate for curvature estimation
- WRF/CMAQ integration examples
- Extended planetary applications (Mars, Titan)

---

**Repository:** https://github.com/DavidEngland/ABL  
**Contact:** david.england@uah.edu
