!==============================================================================
! Example: Integrating Ri curvature corrections into WRF
!
! This shows how to add corrections to module_bl_ysu.F (YSU scheme)
! Similar approach works for MYNN, ACM2, or other PBL schemes.
!
! STEPS:
!   1. Add correction module to WRF/phys/
!   2. Call from existing Ri computation
!   3. Tune parameters based on GABLS1 validation
!==============================================================================

MODULE module_ri_correction
   !---------------------------------------------------------------------------
   ! Richardson number curvature corrections for WRF PBL schemes
   ! Based on: England et al. (2025) "Grid-Dependent Stability Function
   !           Corrections for Arctic Boundary Layers"
   !---------------------------------------------------------------------------
   
   IMPLICIT NONE
   PRIVATE
   
   ! Public interfaces
   PUBLIC :: apply_simple_ri_correction
   PUBLIC :: apply_stability_function_correction
   PUBLIC :: estimate_curvature_proxy
   
   ! Parameters (tunable, see rct_config.yaml for defaults)
   REAL, PARAMETER :: D_default = 0.5      ! Damping amplitude
   REAL, PARAMETER :: p_exponent = 0.5     ! Grid scaling exponent
   REAL, PARAMETER :: q_exponent = 1.0     ! Ri scaling exponent
   REAL, PARAMETER :: dz_ref = 10.0        ! Reference grid spacing (m)
   REAL, PARAMETER :: Ri_ref = 0.2         ! Reference Richardson number
   
CONTAINS

   !---------------------------------------------------------------------------
   ! Simple multiplicative correction (STRATEGY 1)
   !---------------------------------------------------------------------------
   SUBROUTINE apply_simple_ri_correction(Ri_bulk, dz, Ri_corrected, D_param)
      !------------------------------------------------------------------------
      ! Apply grid-dependent correction to bulk Richardson number.
      !
      ! Inputs:
      !   Ri_bulk     - Original bulk Richardson number
      !   dz          - Grid spacing for this layer (m)
      !   D_param     - Damping amplitude (optional, default=0.5)
      !
      ! Output:
      !   Ri_corrected - Corrected Richardson number
      !------------------------------------------------------------------------
      
      REAL, INTENT(IN)  :: Ri_bulk, dz
      REAL, INTENT(IN), OPTIONAL :: D_param
      REAL, INTENT(OUT) :: Ri_corrected
      
      REAL :: D, C, grid_ratio, ri_ratio
      
      ! Use provided D or default
      IF (PRESENT(D_param)) THEN
         D = D_param
      ELSE
         D = D_default
      END IF
      
      ! Compute correction factor C
      grid_ratio = (dz / dz_ref)**p_exponent
      ri_ratio = (ABS(Ri_bulk) / Ri_ref)**q_exponent
      
      C = 1.0 + D * grid_ratio * ri_ratio
      
      ! Clamp C to reasonable range (prevents extreme values)
      C = MAX(0.5, MIN(C, 2.0))
      
      ! Apply correction
      Ri_corrected = Ri_bulk * C
      
   END SUBROUTINE apply_simple_ri_correction
   
   
   !---------------------------------------------------------------------------
   ! Stability function correction (STRATEGY 2)
   !---------------------------------------------------------------------------
   SUBROUTINE apply_stability_function_correction(Ri, dz, f_base, f_corrected, D_param)
      !------------------------------------------------------------------------
      ! Modify stability function to account for grid-dependent curvature.
      !
      ! Inputs:
      !   Ri          - Richardson number
      !   dz          - Grid spacing (m)
      !   f_base      - Base stability function value f(Ri)
      !   D_param     - Damping amplitude (optional)
      !
      ! Output:
      !   f_corrected - Grid-aware stability function f*(Ri, Δz)
      !------------------------------------------------------------------------
      
      REAL, INTENT(IN)  :: Ri, dz, f_base
      REAL, INTENT(IN), OPTIONAL :: D_param
      REAL, INTENT(OUT) :: f_corrected
      
      REAL :: D, C, grid_ratio, ri_ratio
      
      ! Use provided D or default
      IF (PRESENT(D_param)) THEN
         D = D_param
      ELSE
         D = D_default
      END IF
      
      ! Compute correction factor (same as simple correction)
      grid_ratio = (dz / dz_ref)**p_exponent
      ri_ratio = (ABS(Ri) / Ri_ref)**q_exponent
      
      C = 1.0 + D * grid_ratio * ri_ratio
      C = MAX(0.5, MIN(C, 2.0))
      
      ! Modify stability function
      ! This extends the tail: f* = f / C
      ! Near Ri_c, this allows more mixing on coarse grids
      f_corrected = f_base / C
      
      ! Ensure f* stays in physical range [0, 1]
      f_corrected = MAX(0.0, MIN(f_corrected, 1.0))
      
   END SUBROUTINE apply_stability_function_correction
   
   
   !---------------------------------------------------------------------------
   ! Curvature proxy estimation
   !---------------------------------------------------------------------------
   SUBROUTINE estimate_curvature_proxy(theta, U, V, z, kappa)
      !------------------------------------------------------------------------
      ! Estimate curvature proxy from profile shape.
      !
      ! Inputs:
      !   theta(1:3)  - Potential temperature at 3 consecutive levels (K)
      !   U(1:3)      - U wind at 3 levels (m/s)
      !   V(1:3)      - V wind at 3 levels (m/s)
      !   z(1:3)      - Heights of the 3 levels (m)
      !
      ! Output:
      !   kappa       - Dimensionless curvature proxy
      !------------------------------------------------------------------------
      
      REAL, INTENT(IN)  :: theta(3), U(3), V(3), z(3)
      REAL, INTENT(OUT) :: kappa
      
      REAL :: dz1, dz2, dtheta_dz1, dtheta_dz2, d2theta_dz2
      REAL :: dU_dz1, dU_dz2, d2U_dz2
      REAL :: dV_dz1, dV_dz2, d2V_dz2
      REAL :: kappa_theta, kappa_U, kappa_V
      
      ! Grid spacings
      dz1 = z(2) - z(1)
      dz2 = z(3) - z(2)
      
      ! First derivatives
      dtheta_dz1 = (theta(2) - theta(1)) / dz1
      dtheta_dz2 = (theta(3) - theta(2)) / dz2
      dU_dz1 = (U(2) - U(1)) / dz1
      dU_dz2 = (U(3) - U(2)) / dz2
      dV_dz1 = (V(2) - V(1)) / dz1
      dV_dz2 = (V(3) - V(2)) / dz2
      
      ! Second derivatives (central difference)
      d2theta_dz2 = (dtheta_dz2 - dtheta_dz1) / (0.5 * (dz1 + dz2))
      d2U_dz2 = (dU_dz2 - dU_dz1) / (0.5 * (dz1 + dz2))
      d2V_dz2 = (dV_dz2 - dV_dz1) / (0.5 * (dz1 + dz2))
      
      ! Curvature contributions (dimensionless)
      ! kappa ~ |d²f/dz²| / (|df/dz| / Δz)
      IF (ABS(dtheta_dz1 + dtheta_dz2) > 1e-6) THEN
         kappa_theta = ABS(d2theta_dz2) / (ABS(dtheta_dz1 + dtheta_dz2) / (dz1 + dz2))
      ELSE
         kappa_theta = 0.0
      END IF
      
      IF (ABS(dU_dz1 + dU_dz2) > 1e-6) THEN
         kappa_U = ABS(d2U_dz2) / (ABS(dU_dz1 + dU_dz2) / (dz1 + dz2))
      ELSE
         kappa_U = 0.0
      END IF
      
      IF (ABS(dV_dz1 + dV_dz2) > 1e-6) THEN
         kappa_V = ABS(d2V_dz2) / (ABS(dV_dz1 + dV_dz2) / (dz1 + dz2))
      ELSE
         kappa_V = 0.0
      END IF
      
      ! Combined curvature (simple average)
      kappa = (kappa_theta + kappa_U + kappa_V) / 3.0
      
      ! Clamp to reasonable range
      kappa = MAX(0.0, MIN(kappa, 5.0))
      
   END SUBROUTINE estimate_curvature_proxy
   
END MODULE module_ri_correction


!==============================================================================
! INTEGRATION INTO module_bl_ysu.F
!==============================================================================
!
! In subroutine ysu (around line 200-300 where Ri is computed):
!
! BEFORE (original code):
! -------------------------
!     DO k = kts, kte-1
!        dz = z(k+1) - z(k)
!        dtheta = th(k+1) - th(k)
!        du = u(k+1) - u(k)
!        dv = v(k+1) - v(k)
!        
!        ! Bulk Richardson number
!        shear2 = du**2 + dv**2
!        IF (shear2 > 1.e-10) THEN
!           Rib(k) = g * dtheta * dz / (th(k) * shear2)
!        ELSE
!           Rib(k) = 0.0
!        END IF
!        
!        ! Apply cutoff at Ri_c = 0.25
!        IF (Rib(k) < 0.25) THEN
!           ! ... compute mixing ...
!        ELSE
!           ! ... suppress mixing ...
!        END IF
!     END DO
!
!
! AFTER (with corrections):
! -------------------------
!     USE module_ri_correction, ONLY: apply_simple_ri_correction
!     
!     DO k = kts, kte-1
!        dz = z(k+1) - z(k)
!        dtheta = th(k+1) - th(k)
!        du = u(k+1) - u(k)
!        dv = v(k+1) - v(k)
!        
!        ! Bulk Richardson number (original)
!        shear2 = du**2 + dv**2
!        IF (shear2 > 1.e-10) THEN
!           Rib(k) = g * dtheta * dz / (th(k) * shear2)
!        ELSE
!           Rib(k) = 0.0
!        END IF
!        
!        ! *** NEW: Apply curvature correction ***
!        CALL apply_simple_ri_correction(Rib(k), dz, Rib_corrected, D_param=0.5)
!        
!        ! Use corrected Ri for rest of calculation
!        IF (Rib_corrected < 0.25) THEN
!           ! ... compute mixing (now uses Rib_corrected) ...
!        ELSE
!           ! ... suppress mixing ...
!        END IF
!     END DO
!
!
! NOTES:
!   - Only 1 line added (the CALL statement)
!   - Can toggle on/off with namelist flag: apply_ri_correction = .true.
!   - D_param can be tuned via namelist (default 0.5)
!   - Validation: Run GABLS1, check if h_BL and LLJ height improve
!
!==============================================================================


!==============================================================================
! ALTERNATIVE: Stability function correction (more accurate)
!==============================================================================
!
! If your PBL scheme uses explicit stability functions (like MYNN),
! modify the function evaluation instead:
!
! BEFORE:
! -------
!     f_m = stability_function(Rib(k))  ! e.g., f_m = 1/(1+5*Ri)
!     K_m(k) = l_mix**2 * shear * f_m
!
!
! AFTER:
! ------
!     USE module_ri_correction, ONLY: apply_stability_function_correction
!     
!     f_m_base = stability_function(Rib(k))
!     CALL apply_stability_function_correction(Rib(k), dz, f_m_base, f_m_corrected)
!     K_m(k) = l_mix**2 * shear * f_m_corrected
!
!
! This extends the f(Ri) tail for coarse grids, allowing more mixing
! near the critical Richardson number.
!
!==============================================================================


!==============================================================================
! TESTING CHECKLIST
!==============================================================================
!
! 1. Neutral profile test:
!    - Run with constant theta, linear U profile
!    - Expect: No change from corrections (C → 1)
!
! 2. GABLS1 case:
!    - Single-column mode, 9-hour simulation
!    - Check: Boundary layer height h_BL, low-level jet height z_LLJ
!    - With corrections: Should match LES better (h_BL ~ 200m, z_LLJ ~ 250m)
!
! 3. Parameter sensitivity:
!    - Try D = [0.3, 0.5, 0.7]
!    - Plot h_BL vs. time for each
!    - Choose D that best matches observations
!
! 4. 3D case:
!    - Run full 3D WRF simulation
!    - Check stability: timestep crashes, CFL violations
!    - Validate: Surface fluxes, 2m temperature, 10m wind
!
!==============================================================================
