# Contributing to ABL Toolkit

Thank you for your interest in contributing to the Atmospheric Boundary Layer Curvature-Aware MOST Toolkit!

## 🚀 Quick Start for Contributors

1. **Fork** the repository on GitHub
2. **Clone** your fork locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/ABL.git
   cd ABL
   ```
3. **Set up development environment**:
   ```bash
   bash setup_dev.sh
   source .venv/bin/activate
   ```
4. **Create a feature branch**:
   ```bash
   git checkout -b feat/your-feature-name
   ```

## 📝 Contribution Types

### Bug Reports
- Use GitHub Issues with the "bug" label
- Include: OS, Python version, minimal reproducible example
- Expected vs actual behavior

### Feature Requests
- Use GitHub Issues with the "enhancement" label
- Describe use case and proposed API

### Code Contributions
- Follow existing code style (PEP 8 for Python, standard Fortran conventions)
- Add docstrings for new functions
- Ensure notebooks remain Colab-compatible
- Run smoke tests (demo notebook should execute without errors)

### Documentation
- Markdown files for new features
- Inline comments for complex algorithms
- Jupyter notebook examples encouraged

## 🔍 Code Review Process

1. **Submit Pull Request** against `main` branch
2. **CI checks** must pass:
   - Notebooks execute without error
   - Requirements install cleanly
   - No merge conflicts
3. **Review** by maintainers (McNider, Biazar, England)
4. **Iterate** on feedback
5. **Merge** when approved

## 📋 Style Guidelines

### Python
- PEP 8 compliant (use `ruff` or `flake8`)
- Type hints encouraged for public APIs
- Docstrings: NumPy style

### Fortran
- Standard Fortran 90+ conventions
- UPPERCASE for module/subroutine names
- Inline comments for physics/math

### Notebooks
- Clear markdown headers
- Reproducible (no hardcoded paths)
- Outputs cleared before commit (optional)

## 🧪 Testing

Run the main demo notebook:
```bash
jupyter lab notebooks/Curvature_Demo.ipynb
# Execute all cells; verify plots appear correctly
```

For Fortran:
```bash
cd implementations
gfortran -c McNider_1DBLM_fc_module.f90
gfortran -o test_driver mcnb_fc_driver.f90 McNider_1DBLM_fc_module.o
./test_driver
```

## 📧 Communication

- **GitHub Issues**: Bug reports, feature requests
- **GitHub Discussions**: General questions, ideas
- **Email**: Direct collaboration inquiries to david.england@uah.edu

## 🙏 Recognition

Contributors will be acknowledged in:
- README.md contributors section
- Relevant publication acknowledgments (if applicable)
- Release notes

## 📜 License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

**Questions?** Open a GitHub Discussion or contact the maintainers.
