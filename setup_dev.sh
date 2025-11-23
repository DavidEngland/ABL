#!/usr/bin/env bash
set -euo pipefail

VENV_DIR="${1:-.venv}"
REQ_FILE="$(dirname "$0")/requirements.txt"

echo "Creating venv at ${VENV_DIR} (Python from PATH)..."
python3 -m venv "${VENV_DIR}"
# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"

echo "Upgrading pip and installing Python requirements..."
python -m pip install --upgrade pip setuptools wheel
if [ -f "${REQ_FILE}" ]; then
    python -m pip install -r "${REQ_FILE}"
else
    echo "Warning: requirements.txt not found at ${REQ_FILE}"
fi

echo
echo "Checking for Fortran compiler (gfortran)..."
if command -v gfortran >/dev/null 2>&1; then
    echo "gfortran found: $(gfortran --version | head -n1)"
else
    echo "gfortran not found on PATH."
    echo
    echo "Install instructions (choose one):"
    echo "  Debian/Ubuntu: sudo apt update && sudo apt install -y gfortran"
    echo "  macOS (Homebrew): brew install gcc   # provides gfortran as gfortran-<ver> or gfortran"
    echo "  RHEL/CentOS/Fedora: sudo dnf install -y gcc-gfortran"
    echo "  WSL: use Debian/Ubuntu instructions inside WSL"
    echo
    echo "After installing gfortran, re-run this script or simply activate the venv and continue."
fi

echo
echo "Done. To activate the venv in this shell run:"
echo "  source ${VENV_DIR}/bin/activate"
echo "Then compile Fortran code (example):"
echo "  gfortran -O2 -std=f2008 -o mcnb_driver implementations/McNider_1DBLM_fc_module.f90 implementations/mcnb_fc_driver.f90"
