# Integration Guide: fc Correction into McNider 1DBLM

## Quick Start (3 Steps)

### Step 1: Add Module to Compilation
Add `McNider_1DBLM_fc_module.f90` to your Makefile or compile command:
```bash
gfortran -c McNider_1DBLM_fc_module.f90
gfortran -o model model.f90 McNider_1DBLM_fc_module.o
```

### Step 2: Modify PROFILEK Subroutine
In your existing `model.f90`, find the `PROFILEK` subroutine and add the module call at the **end** (after K-profile is computed but before returning):

```fortran
SUBROUTINE PROFILEK(U, V, T, Q, Z, DELZ, KZ, STABFUNC, ZI, LAT, T0, &
                    TSTAR, USTAR, L, DELT, VONKAR, KM, KH)
    USE FC_CORRECTION  ! Add this line at top of subroutine
    
    ! ...existing code that computes KM, KH...
    
    ! Add this block just before END SUBROUTINE
    LOGICAL, SAVE :: FIRST_FC_CALL = .TRUE.
    CALL APPLY_FC_TO_PROFILE(Z, KM, KH, L, T, U, V, KZ, G, FIRST_FC_CALL)
    
END SUBROUTINE PROFILEK
```

### Step 3: Review Diagnostics
After first run, check `FC_DIAGNOSTICS.dat`:
- B values should be ~1.0–1.5 in near-neutral, ~1.5–2.5 in strong SBL
- FC values should be ~0.8–1.0 near-neutral, ~0.3–0.7 in strong SBL
- KM_NEW should be less than KM_OLD when B > 1.05

---

## Detailed Implementation Notes

### Where NOT to Apply fc
**Do not apply to:**
1. Surface layer (first level) — fluxes prescribed by USTAR, TSTAR
2. Top boundary (level KZ) — often prescribed or open boundary
3. Unstable conditions (L < 0) — correction designed for stable only

**Current implementation:** Module applies to levels 2 through KZ-1 automatically.

### Tuning Parameters (Advanced)

If default fc is too strong (excessive damping):
```fortran
! In your main program or PROFILEK, pass custom alpha
REAL :: MY_ALPHA = 0.7  ! Reduce from default 1.0
CALL COMPUTE_FC(B, ZETA, DELTA_Z, FC, ALPHA=MY_ALPHA)
```

If default fc is too weak (insufficient bias reduction):
```fortran
REAL :: MY_ALPHA = 1.3  ! Increase from default 1.0
```

**Recommended tuning sequence:**
1. Run GABLS1 case with default α=1.0
2. Compare surface temperature bias: target 0.5–1.0 K improvement
3. Adjust α in 0.1 increments; rerun
4. Validate on 2–3 additional stable nights before finalizing

### Connecting to Your Existing Variables

**Your variable → Module expectation:**

| Your Variable | Module Needs | Notes |
|---------------|--------------|-------|
| `L` | Monin-Obukhov length | Already computed in SFCFLUX |
| `Z(KZ)` | Heights | From INITIALIZE |
| `T(KZ)` | Potential temperature | θ (not T_v) |
| `U(KZ), V(KZ)` | Winds | Components, not speed |
| `KM(KZ), KH(KZ)` | Diffusivities | Modified in-place |
| `G` | Gravity constant | = 9.81 m/s² |

**If your model uses different names:**
```fortran
! Example: if you have TH instead of T
CALL APPLY_FC_TO_PROFILE(Z, KM, KH, L, TH, U, V, KZ, G, FIRST_FC_CALL)
```

### Validating Neutral Preservation

**Unit test (add to your test suite):**
```fortran
! Set up neutral conditions
L_TEST = 1.0E5  ! Large positive (near-neutral)
ZETA_TEST = 0.01
B_TEST = 1.02   ! Minimal bias
CALL COMPUTE_FC(B_TEST, ZETA_TEST, 50.0, FC_TEST)
IF (ABS(FC_TEST - 1.0) .GT. 0.05) THEN
    WRITE(*,*) 'WARNING: Neutral preservation violated, fc=', FC_TEST
END IF
```

**Expected:** fc ≈ 0.98–1.00 when ζ < 0.05 and B < 1.1

---

## Reconciling with Your D Parameter

If you have existing D values from empirical fits (e.g., VB spreadsheet or paper draft):

**Conversion formula:**
```
alpha_equivalent = D * (Ri_b / 0.25) / (B - 1.0)
```

**Example:**
- Your fitted D = 0.3
- Measured Ri_b = 0.20, Ri_g = 0.50 → B = 2.5
- Then α = 0.3 * (0.20/0.25) / 1.5 ≈ 0.16

**Action:** Use this α value as starting point instead of default 1.0.

---

## Diagnostic Outputs Explained

**FC_DIAGNOSTICS.dat columns:**

1. **LEVEL**: Vertical level index (2 to KZ-1)
2. **Z_G(m)**: Geometric mean height of layer
3. **B**: Bias ratio Ri_g(z_g)/Ri_b
4. **FC**: Applied correction factor
5. **KM_OLD**: Original momentum diffusivity
6. **KM_NEW**: Corrected momentum diffusivity
7. **KH_OLD**: Original heat diffusivity
8. **KH_NEW**: Corrected heat diffusivity

**What to look for:**
- Stable nights: B should increase with height (up to jet nose)
- Strong inversions: fc should drop to ~0.3–0.5 in first few levels
- Neutral periods: B ≈ 1.0, fc ≈ 1.0 (no correction applied)

**Red flags:**
- B < 1.0 persistently → Check Ri calculation signs
- fc = FC_MIN often → α too large or C too small (over-damping)
- No variation in fc → Check that L is being updated each timestep

---

## Troubleshooting Common Issues

### Issue: "Module not found" compilation error
**Fix:** Ensure module file is compiled **before** main program:
```bash
gfortran -c McNider_1DBLM_fc_module.f90
gfortran -c model.f90  # Now it can find the module
gfortran -o model model.o McNider_1DBLM_fc_module.o
```

### Issue: Negative KM or KH after correction
**Fix:** Floor is already applied (fc_min=0.2), but if you see negative values:
```fortran
! Add explicit guard in APPLY_FC_TO_PROFILE after correction
IF (KM(I) .LT. 0.0) KM(I) = 1.0E-6
IF (KH(I) .LT. 0.0) KH(I) = 1.0E-6
```

### Issue: Surface temperature still too warm
**Possible causes:**
1. α too small → Increase to 1.2–1.5
2. Correction not reaching first interior layer → Check loop bounds
3. Surface flux issue unrelated to K-profile → Validate USTAR, TSTAR separately

### Issue: Model crashes with "floating point exception"
**Check:**
1. L = 0 case handled? (module guards against this)
2. Division by zero in shear? (module uses TINY safeguard)
3. Exponential overflow? (should not occur with fc_min floor, but add trap):
```fortran
IF (EXPONENT .LT. -20.0) EXPONENT = -20.0  ! Before FC = EXP(EXPONENT)
```

---

## Performance Impact

**Expected overhead:** <3% additional CPU time
- Bias ratio computation: negligible (simple arithmetic)
- Exponential evaluation: ~10 FLOPS per level per timestep
- Diagnostic I/O: minimal (writes every 10 levels only)

**To disable diagnostics** (production runs):
```fortran
! Comment out WRITE statement in APPLY_FC_TO_PROFILE
! WRITE(99,'(I5,7F8.3)') I, Z_G, B, FC, KM_OLD, KM(I), KH_OLD, KH(I)
```

---

## Next Steps After Integration

1. **Baseline run:** Run unchanged model on GABLS1, save outputs
2. **Corrected run:** Activate fc module, rerun same case
3. **Compare:** Surface T bias, inversion height, LLJ timing
4. **Tune:** Adjust α if needed (see tuning sequence above)
5. **Validate:** Test on 3–5 additional stable nights
6. **Document:** Record final α, fc_min choices in BL.DAT comments

---

## Contact for Support

David E. England  
david.england@uah.edu  

**Please provide when requesting help:**
1. Sample of FC_DIAGNOSTICS.dat (first 20 lines)
2. Typical L value range for your cases
3. STABFUNC choice (1, 2, or 3) from your BL.DAT
4. Whether using England-McNider (FLUXTYPE=2) or Businger (FLUXTYPE=1)
