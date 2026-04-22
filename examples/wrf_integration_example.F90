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
   PUBLIC :: compute_dynamic_ric
   PUBLIC :: update_turbulence_state
   PUBLIC :: compute_stability_functions_hybrid
   PUBLIC :: STAB_FORM_LINEAR, STAB_FORM_QUADRATIC
   
   ! Parameters (tunable, see rct_config.yaml for defaults)
   REAL, PARAMETER :: D_default = 0.8            ! Damping amplitude
   REAL, PARAMETER :: p_exponent = 1.0           ! Grid scaling exponent
   REAL, PARAMETER :: q_exponent = 2.0           ! Stability exponent; q>=2 preserves neutral slope
   REAL, PARAMETER :: dz_ref = 10.0              ! Reference grid spacing (m)
   REAL, PARAMETER :: zeta_ref = 0.5             ! Reference z/L for stable damping
   REAL, PARAMETER :: Ri_ref = 0.2               ! Reference Richardson number
   REAL, PARAMETER :: fc_min_default = 0.2       ! Lower bound for damping factor
   REAL, PARAMETER :: B_threshold_default = 1.05 ! Curvature-bias activation threshold

   ! Dynamic critical Richardson number parameters
   REAL, PARAMETER :: Ri_c0 = 0.25
   REAL, PARAMETER :: Ri_c_min = 0.20
   REAL, PARAMETER :: Ri_c_max = 1.00
   REAL, PARAMETER :: gamma_ref = 0.010          ! K m^-1
   REAL, PARAMETER :: shear_ref = 0.050          ! s^-1
   REAL, PARAMETER :: tke_ref = 0.10             ! m^2 s^-2
   REAL, PARAMETER :: alpha_gamma = 0.35
   REAL, PARAMETER :: alpha_shear = 0.25
   REAL, PARAMETER :: alpha_tke = 0.60
   REAL, PARAMETER :: alpha_kappa = 0.15

   ! Stable-function selector IDs (for compute_stability_functions_hybrid)
   INTEGER, PARAMETER :: STAB_FORM_LINEAR = 1
   INTEGER, PARAMETER :: STAB_FORM_QUADRATIC = 2
   
CONTAINS

   !---------------------------------------------------------------------------
   ! Simple multiplicative correction (STRATEGY 1)
   !---------------------------------------------------------------------------
   SUBROUTINE apply_simple_ri_correction(Ri_bulk, dz, Ri_corrected, zeta, B_ratio, D_param, fc_min, B_threshold)
      !------------------------------------------------------------------------
      ! Apply grid-dependent correction to bulk Richardson number.
      !
      ! Inputs:
      !   Ri_bulk     - Original bulk Richardson number
      !   dz          - Grid spacing for this layer (m)
      !   zeta        - Stability parameter z/L (optional)
      !   B_ratio     - Curvature bias ratio Ri_g(z_g)/Ri_b (optional)
      !   D_param     - Damping amplitude (optional, default=0.5)
      !   fc_min      - Minimum damping factor floor (optional)
      !   B_threshold - Activation threshold for B_ratio (optional)
      !
      ! Output:
      !   Ri_corrected - Corrected Richardson number
      !------------------------------------------------------------------------
      
      REAL, INTENT(IN)  :: Ri_bulk, dz
      REAL, INTENT(IN), OPTIONAL :: zeta, B_ratio
      REAL, INTENT(IN), OPTIONAL :: D_param
      REAL, INTENT(IN), OPTIONAL :: fc_min, B_threshold
      REAL, INTENT(OUT) :: Ri_corrected
      
      REAL :: fc
      
      CALL compute_grid_damping_factor(Ri_bulk, dz, fc, zeta, B_ratio, D_param, fc_min, B_threshold)
      
      ! Apply damping as an Ri increase: Ri* = Ri / fc, with 0 < fc <= 1
      Ri_corrected = Ri_bulk / fc
      
   END SUBROUTINE apply_simple_ri_correction
   
   
   !---------------------------------------------------------------------------
   ! Stability function correction (STRATEGY 2)
   !---------------------------------------------------------------------------
   SUBROUTINE apply_stability_function_correction(Ri, dz, f_base, f_corrected, zeta, B_ratio, D_param, fc_min, B_threshold)
      !------------------------------------------------------------------------
      ! Modify stability function to account for grid-dependent curvature.
      !
      ! Inputs:
      !   Ri          - Richardson number
      !   dz          - Grid spacing (m)
      !   f_base      - Base stability function value f(Ri)
      !   zeta        - Stability parameter z/L (optional)
      !   B_ratio     - Curvature bias ratio Ri_g(z_g)/Ri_b (optional)
      !   D_param     - Damping amplitude (optional)
      !   fc_min      - Minimum damping factor floor (optional)
      !   B_threshold - Activation threshold for B_ratio (optional)
      !
      ! Output:
      !   f_corrected - Grid-aware stability function f*(Ri, Δz)
      !------------------------------------------------------------------------
      
      REAL, INTENT(IN)  :: Ri, dz, f_base
      REAL, INTENT(IN), OPTIONAL :: zeta, B_ratio
      REAL, INTENT(IN), OPTIONAL :: D_param
      REAL, INTENT(IN), OPTIONAL :: fc_min, B_threshold
      REAL, INTENT(OUT) :: f_corrected
      
      REAL :: fc
      
      CALL compute_grid_damping_factor(Ri, dz, fc, zeta, B_ratio, D_param, fc_min, B_threshold)
      
      ! Damped tail correction for coarse stable layers: f* = f * fc
      f_corrected = f_base * fc
      
      ! Ensure f* stays in physical range [0, 1]
      f_corrected = MAX(0.0, MIN(f_corrected, 1.0))
      
   END SUBROUTINE apply_stability_function_correction


   !---------------------------------------------------------------------------
   ! Dynamic critical Richardson number (recent framework)
   !---------------------------------------------------------------------------
   SUBROUTINE compute_dynamic_ric(Gamma, S, TKE_prev, kappa, Ri_c_star)
      !------------------------------------------------------------------------
      ! Compute local dynamic critical Richardson number Ri_c*.
      !
      ! Inputs:
      !   Gamma       - dtheta/dz (K/m)
      !   S           - Shear magnitude (1/s)
      !   TKE_prev    - Previous-step TKE (m^2/s^2)
      !   kappa       - Curvature proxy from estimate_curvature_proxy
      !
      ! Output:
      !   Ri_c_star   - Dynamic critical Richardson number in [0.20, 1.00]
      !------------------------------------------------------------------------

      REAL, INTENT(IN)  :: Gamma, S, TKE_prev, kappa
      REAL, INTENT(OUT) :: Ri_c_star

      REAL :: gamma_term, shear_term, tke_term, kappa_term

      gamma_term = alpha_gamma * ((Gamma / gamma_ref) - 1.0)
      shear_term = alpha_shear * ((S / shear_ref) - 1.0)
      tke_term = alpha_tke * (MAX(TKE_prev, 0.0) / tke_ref)
      kappa_term = alpha_kappa * MAX(kappa, 0.0)

      Ri_c_star = Ri_c0 * (1.0 + gamma_term + shear_term + tke_term + kappa_term)
      Ri_c_star = MAX(Ri_c_min, MIN(Ri_c_star, Ri_c_max))

   END SUBROUTINE compute_dynamic_ric


   !---------------------------------------------------------------------------
   ! Hysteresis transition helper for intermittent turbulence
   !---------------------------------------------------------------------------
   SUBROUTINE update_turbulence_state(Ri_local, Ri_c_star, was_turbulent, is_turbulent)
      !------------------------------------------------------------------------
      ! Uses two-threshold hysteresis:
      !   suppress if Ri > 1.5*Ri_c*
      !   restart  if Ri < 0.5*Ri_c*
      !------------------------------------------------------------------------

      REAL, INTENT(IN)    :: Ri_local, Ri_c_star
      LOGICAL, INTENT(IN) :: was_turbulent
      LOGICAL, INTENT(OUT):: is_turbulent

      IF (was_turbulent) THEN
         is_turbulent = (Ri_local <= 1.5 * Ri_c_star)
      ELSE
         is_turbulent = (Ri_local < 0.5 * Ri_c_star)
      END IF

   END SUBROUTINE update_turbulence_state


   !---------------------------------------------------------------------------
   ! Hybrid stability functions: MOST in unstable, Ri-based in stable
   !---------------------------------------------------------------------------
   SUBROUTINE compute_stability_functions_hybrid(Ri, zeta, f_m, f_h, stable_form)
      !------------------------------------------------------------------------
      ! Unstable/near-neutral branch (Ri <= 0):
      !   Use MOST functions parameterized by zeta = z/L.
      ! Stable branch (Ri > 0):
      !   Use Ri-based stability functions with selectable forms.
      !
      ! Inputs:
      !   Ri          - Richardson number
      !   zeta        - z/L (used on unstable branch)
      !   stable_form - Selector (optional):
      !                 STAB_FORM_LINEAR (default)
      !                 STAB_FORM_QUADRATIC
      !
      ! Outputs:
      !   f_m, f_h    - Stability functions used in K = l_mix^2 * S * f
      !------------------------------------------------------------------------

      REAL, INTENT(IN) :: Ri, zeta
      INTEGER, INTENT(IN), OPTIONAL :: stable_form
      REAL, INTENT(OUT) :: f_m, f_h

      INTEGER :: form_sel
      REAL :: phi_m, phi_h, ri_loc

      form_sel = STAB_FORM_LINEAR
      IF (PRESENT(stable_form)) form_sel = stable_form

      IF (Ri <= 0.0) THEN
         ! Businger-Dyer style unstable MOST (zeta < 0).
         phi_m = (1.0 - 16.0 * MIN(zeta, 0.0))**(-0.25)
         phi_h = (1.0 - 16.0 * MIN(zeta, 0.0))**(-0.50)
         f_m = 1.0 / (phi_m * phi_m)
         f_h = 1.0 / (phi_m * phi_h)
      ELSE
         ri_loc = MAX(Ri, 0.0)
         SELECT CASE (form_sel)
         CASE (STAB_FORM_LINEAR)
            ! Linear-in-Ri denominator: tune coefficients per host scheme.
            f_m = 1.0 / (1.0 + 5.0 * ri_loc)
            f_h = 1.0 / (1.0 + 7.5 * ri_loc)
         CASE (STAB_FORM_QUADRATIC)
            ! Quadratic tail for stronger damping in very stable conditions.
            f_m = 1.0 / (1.0 + 4.7 * ri_loc + 12.0 * ri_loc * ri_loc)
            f_h = 1.0 / (1.0 + 7.8 * ri_loc + 20.0 * ri_loc * ri_loc)
         CASE DEFAULT
            ! Safe fallback
            f_m = 1.0 / (1.0 + 5.0 * ri_loc)
            f_h = 1.0 / (1.0 + 7.5 * ri_loc)
         END SELECT

         ! Hook for user-defined stable functions:
         ! Add your CASE(...) branch above and assign f_m, f_h.
      END IF

      ! Keep in physical bounds
      f_m = MAX(0.0, MIN(f_m, 1.0))
      f_h = MAX(0.0, MIN(f_h, 1.0))

   END SUBROUTINE compute_stability_functions_hybrid
   
   
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


   !---------------------------------------------------------------------------
   ! Internal helper: neutral-preserving, stable-only exponential damping
   !---------------------------------------------------------------------------
   SUBROUTINE compute_grid_damping_factor(Ri, dz, fc, zeta, B_ratio, D_param, fc_min, B_threshold)
      REAL, INTENT(IN) :: Ri, dz
      REAL, INTENT(IN), OPTIONAL :: zeta, B_ratio, D_param, fc_min, B_threshold
      REAL, INTENT(OUT) :: fc

      REAL :: D, zeta_loc, b_loc, b_thresh, b_excess, grid_ratio, zeta_ratio, fc_floor

      D = D_default
      IF (PRESENT(D_param)) D = D_param

      fc_floor = fc_min_default
      IF (PRESENT(fc_min)) fc_floor = fc_min
      fc_floor = MAX(0.05, MIN(fc_floor, 1.0))

      b_thresh = B_threshold_default
      IF (PRESENT(B_threshold)) b_thresh = B_threshold

      ! Stable-only gating: do not alter neutral/unstable layers.
      IF (Ri <= 0.0 .OR. dz <= 0.0) THEN
         fc = 1.0
         RETURN
      END IF

      zeta_loc = Ri
      IF (PRESENT(zeta)) zeta_loc = MAX(zeta, 0.0)

      ! If B is unavailable, fall back to Ri-based activation.
      b_loc = Ri / Ri_ref
      IF (PRESENT(B_ratio)) b_loc = MAX(B_ratio, 0.0)

      b_excess = MAX(0.0, b_loc - b_thresh)
      IF (b_excess <= 0.0) THEN
         fc = 1.0
         RETURN
      END IF

      grid_ratio = dz / dz_ref
      zeta_ratio = zeta_loc / zeta_ref

      fc = EXP(-D * (grid_ratio**p_exponent) * (zeta_ratio**q_exponent) * b_excess)
      fc = MAX(fc_floor, MIN(fc, 1.0))

   END SUBROUTINE compute_grid_damping_factor
   
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
!           S = SQRT(shear2) / MAX(dz, 1.e-3)
!           l_mix = MIN(l_max, karman * z(k))
!           ! Use your scheme's Ri-based stability functions here
!           f_m = stability_function_m(Rib(k))
!           f_h = stability_function_h(Rib(k))
!           K_m(k) = l_mix**2 * S * f_m
!           K_h(k) = l_mix**2 * S * f_h
!        ELSE
!           ! Suppressed but not zero: keep small background diffusivity
!           K_m(k) = K_min
!           K_h(k) = K_min
!        END IF
!     END DO
!
!
! AFTER (with corrections):
! -------------------------
!     USE module_ri_correction, ONLY: apply_simple_ri_correction, &
!                                      compute_dynamic_ric, update_turbulence_state, &
!                                      compute_stability_functions_hybrid, STAB_FORM_LINEAR
!     
!     DO k = kts, kte-1
!        dz = z(k+1) - z(k)
!        dtheta = th(k+1) - th(k)
!        du = u(k+1) - u(k)
!        dv = v(k+1) - v(k)
!        S = SQRT(MAX(1.e-10, (du/dz)**2 + (dv/dz)**2))
!        
!        ! Bulk Richardson number (original)
!        shear2 = du**2 + dv**2
!        IF (shear2 > 1.e-10) THEN
!           Rib(k) = g * dtheta * dz / (th(k) * shear2)
!        ELSE
!           Rib(k) = 0.0
!        END IF
!        
!        ! Optional: estimate zeta and curvature-bias ratio B from your host scheme
!        ! zeta_k = z(k) / L(k)
!        ! B_k = Ri_g(z_g) / MAX(Rib(k), 1.e-6)
!
!        ! *** NEW: neutral-preserving, stable-only Ri correction ***
!        CALL apply_simple_ri_correction(Rib(k), dz, Rib_corrected, &
!                                        zeta=zeta_k, B_ratio=B_k, D_param=0.8)
!
!        ! *** NEW: dynamic critical Richardson number ***
!        CALL compute_dynamic_ric(Gamma=dtheta/dz, S=S, TKE_prev=tke_old(k), &
!                                 kappa=kappa_prof(k), Ri_c_star=Ric_star)
!
!        ! *** NEW: hysteresis for intermittent turbulence ***
!        CALL update_turbulence_state(Rib_corrected, Ric_star, turb_on_prev(k), turb_on(k))
!        
!        ! Use dynamic Ric* instead of fixed 0.25
!        IF (turb_on(k)) THEN
!           S = SQRT(shear2) / MAX(dz, 1.e-3)
!           l_mix = MIN(l_max, karman * z(k))
!           ! Hybrid closure: MOST(zeta) for unstable, Ri-form for stable
!           CALL compute_stability_functions_hybrid(Rib_corrected, zeta_k, f_m, f_h, &
!                                                   stable_form=STAB_FORM_LINEAR)
!           K_m(k) = l_mix**2 * S * f_m
!           K_h(k) = l_mix**2 * S * f_h
!        ELSE
!           ! In suppressed regime, only background mixing remains
!           K_m(k) = K_min
!           K_h(k) = K_min
!        END IF
!     END DO
!
!
! NOTES:
!   - Stable-only damping avoids changing neutral/unstable layers
!   - Dynamic Ric* replaces fixed 0.25 cutoff and improves regime transitions
!   - Keep feature flags in namelist: apply_ri_correction, use_dynamic_ric
!   - D_param can be tuned via namelist (start with 0.8)
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
!     USE module_ri_correction, ONLY: apply_stability_function_correction, &
!                                      compute_stability_functions_hybrid, STAB_FORM_QUADRATIC
!     
!     ! Base function can be linear or quadratic on stable branch
!     CALL compute_stability_functions_hybrid(Rib(k), zeta_k, f_m_base, f_h_base, &
!                                             stable_form=STAB_FORM_QUADRATIC)
!     CALL apply_stability_function_correction(Rib(k), dz, f_m_base, f_m_corrected, &
!                                              zeta=zeta_k, B_ratio=B_k, D_param=0.8)
!     K_m(k) = l_mix**2 * shear * f_m_corrected
!
!
! This damps over-mixing for coarse stable layers while preserving
! the neutral limit.
!
! NOTE:
!   - In this module, unstable/near-neutral uses MOST functions of zeta=z/L.
!   - Stable branch uses Ri-based functions (linear or quadratic by selector).
!   - To insert your own stable function, add a CASE(...) branch in
!     compute_stability_functions_hybrid and set f_m, f_h there.
!
!==============================================================================


!==============================================================================
! TESTING CHECKLIST
!==============================================================================
!
! 1. Neutral profile test:
!    - Run with constant theta, linear U profile
!    - Expect: No change from corrections (fc → 1)
!
! 2. GABLS1 case:
!    - Single-column mode, 9-hour simulation
!    - Check: Boundary layer height h_BL, low-level jet height z_LLJ
!    - With corrections: Should match LES better (h_BL ~ 200m, z_LLJ ~ 250m)
!
! 3. Parameter sensitivity:
!    - Try D = [0.5, 0.8, 1.2]
!    - Plot h_BL vs. time for each
!    - Choose D that best matches observations
!
! 4. Dynamic Ric* sensitivity:
!    - Compare fixed Ri_c=0.25 vs dynamic Ri_c*
!    - Check intermittent turbulence and LLJ timing
!
! 5. 3D case:
!    - Run full 3D WRF simulation
!    - Check stability: timestep crashes, CFL violations
!    - Validate: Surface fluxes, 2m temperature, 10m wind
!
!==============================================================================
