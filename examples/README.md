# examples/ — Fortran 90 MOST Modules

Active Fortran 90 development for MOST inversion, CBC/Gegenbauer series, and WRF integration.

**Compiler:** gfortran (free-form F90, double precision)  
**Last updated:** April 2026

---

## Modules

### `module_most_profile_utils.F90`
Core MOST profile utilities. Public interface:

| Subroutine | Purpose |
|---|---|
| `zeta_from_rig_safeguarded(profile_id, rig_target, zeta, converged, n_iter)` | Branch-aware safeguarded Newton + bisection fallback for ζ(Ri_g) |
| `fm_fh_from_rig(profile_id, rig, fm, fh)` | Compute φ_m, φ_h from Ri_g (calls safeguarded solver first) |
| `phi_m(zeta, profile_id)` / `phi_h(zeta, profile_id)` | Evaluate similarity functions at ζ |

**Profile IDs:** 1 = Businger-Dyer linear, 2 = power-law McNider, 3 = Högström, 4 = Cheng-Brutsaert  
**Solver strategy:** closed-form seed → safeguarded Newton (analytic dRi/dζ) → bisection fallback within branch-aware bracket.

---

### `module_cbc_legendre_most.F90`
Central Binomial Coefficient (CBC) / Gegenbauer (ultraspherical) series for power-law φ functions.

Key identity used:
$$P_{2n}(0) = (-1)^n \frac{\binom{2n}{n}}{4^n}$$

| Function | Purpose |
|---|---|
| `phi_m_gegenbauer(zeta, N)` | φ_m = (1−b_m ζ)^{−1/4} via Gegenbauer C_n^{(1/4)} series, N terms |
| `phi_h_legendre(zeta, N)` | φ_h = (1−b_h ζ)^{−1/2} via Legendre P_n series, N terms |
| `cbc_coeff(n)` | Central binomial coefficient C(2n,n) |
| `gegenbauer_equator(lambda, N)` | Series sum at ζ_0 = b_m^{-1} (stability-space equator) |

---

### `driver_cbc_gegenbauer_errors.F90`
Driver for convergence and error tabulation. Produces two tables:

- **Table A:** φ_h (Legendre) exact vs N-term series across ζ ∈ [0, 0.9/b_h]
- **Table B:** φ_m (Gegenbauer, λ=1/4) exact vs series at the equator ζ_0 = 1/b_m

Run with:
```bash
./driver_cbc_gegenbauer_errors
# or filter to Table B momentum rows:
./driver_cbc_gegenbauer_errors | awk '/Table B/{f=1;next} f && /^[[:space:]]*-[0-9]/{print}'
```

---

### `wrf_integration_example.F90`
WRF/MYNN-style sketch demonstrating how to plug `zeta_from_rig_safeguarded` into a model's stability-function call.  
See [manuscripts/WRF_Ri_Curvature_Integration_Outline_v01.md](../manuscripts/WRF_Ri_Curvature_Integration_Outline_v01.md) for context.

---

## Build

```bash
cd examples/

# Compile modules (order matters: utils before cbc, cbc before driver)
gfortran -c module_most_profile_utils.F90
gfortran -c module_cbc_legendre_most.F90
gfortran -c driver_cbc_gegenbauer_errors.F90

# Link driver
gfortran -o driver_cbc_gegenbauer_errors \
    driver_cbc_gegenbauer_errors.o \
    module_cbc_legendre_most.o \
    module_most_profile_utils.o

# Run
./driver_cbc_gegenbauer_errors
```

---

## Related

- [drafts/ultraspherical_subsection.md](../drafts/ultraspherical_subsection.md) — Mathematical derivation of ultraspherical structure
- [drafts/SBL_IBEx_expansions.md](../drafts/SBL_IBEx_expansions.md) — IBEx series framework
- [notes/Legendre.md](../notes/Legendre.md) — Legendre/Gegenbauer identities reference
- [overview.md](overview.md) — Jensen's Inequality and grid-bias motivation
