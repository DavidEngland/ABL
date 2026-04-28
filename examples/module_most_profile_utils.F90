!==============================================================================
! Utility module: MOST profile and Ri_g<->zeta helpers for WRF integration
!
! Purpose:
!   - Centralize common phi_m/phi_h profile evaluation.
!   - Provide Ri_g(zeta) and zeta(Ri_g) utilities for stable closures.
!   - Streamline construction of f_m(Ri_g), f_h(Ri_g) from MOST-consistent forms.
!
! Notes:
!   - This module is kept lightweight and self-contained for easy drop-in use.
!   - The default coefficients are illustrative and should be calibrated per scheme.
!==============================================================================

MODULE module_most_profile_utils
   USE module_cbc_legendre_most, ONLY: phi_h_cbc_unstable, phi_m_gegenbauer_unstable, &
                                       compute_cbc_critical_ri_limits
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: PROFILE_BD71, PROFILE_HOG88, PROFILE_BH91
   PUBLIC :: phi_from_zeta
   PUBLIC :: rig_from_zeta
   PUBLIC :: zeta_from_rig_newton
   PUBLIC :: zeta_from_rig_safeguarded
   PUBLIC :: fm_fh_from_rig
   PUBLIC :: cbc_critical_ri_limits

   INTEGER, PARAMETER :: PROFILE_BD71 = 1
   INTEGER, PARAMETER :: PROFILE_HOG88 = 2
   INTEGER, PARAMETER :: PROFILE_BH91 = 3

CONTAINS

   SUBROUTINE phi_from_zeta(profile_id, zeta, phi_m, phi_h)
      !------------------------------------------------------------------------
      ! Return phi_m(zeta), phi_h(zeta) for selected profile.
      ! For zeta < 0, BD71 unstable branch is used for all profile IDs.
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: profile_id
      REAL, INTENT(IN)    :: zeta
      REAL, INTENT(OUT)   :: phi_m, phi_h

      REAL, PARAMETER :: b_m_unst = 16.0
      REAL, PARAMETER :: b_h_unst = 16.0
      INTEGER, PARAMETER :: n_cbc = 12
      REAL :: eta_mag

      IF (zeta < 0.0) THEN
         ! Near neutrality, use the CBC / Gegenbauer recurrences directly.
         eta_mag = -zeta
         IF (MAX(b_m_unst * eta_mag, b_h_unst * eta_mag) <= 0.5) THEN
            phi_m = phi_m_gegenbauer_unstable(zeta, b_m_unst, n_cbc)
            phi_h = phi_h_cbc_unstable(zeta, b_h_unst, n_cbc)
         ELSE
            phi_m = (1.0 - b_m_unst * zeta)**(-0.25)
            phi_h = (1.0 - b_h_unst * zeta)**(-0.50)
         END IF
         RETURN
      END IF

      SELECT CASE (profile_id)
      CASE (PROFILE_BD71)
         ! Classical stable linear coefficients
         phi_m = 1.0 + 5.0 * zeta
         phi_h = 1.0 + 7.8 * zeta
      CASE (PROFILE_HOG88)
         ! HOG88-style linear stable variant
         phi_m = 1.0 + 5.0 * zeta
         phi_h = 0.95 + 7.8 * zeta
      CASE (PROFILE_BH91)
         ! Beljaars-Holtslag-inspired hybrid stable form
         phi_m = 1.0 + 5.0 * zeta + 4.0 * zeta * (1.0 + 0.5 * zeta)**(1.0 / 3.0)
         phi_h = 1.0 + 5.0 * zeta + 4.0 * zeta * (1.0 + 0.5 * zeta)**(1.0 / 3.0)
      CASE DEFAULT
         phi_m = 1.0 + 5.0 * zeta
         phi_h = 1.0 + 7.8 * zeta
      END SELECT
   END SUBROUTINE phi_from_zeta


   SUBROUTINE rig_from_zeta(profile_id, zeta, rig)
      !------------------------------------------------------------------------
      ! Compute gradient Richardson number:
      !   Ri_g = zeta * phi_h / phi_m^2
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: profile_id
      REAL, INTENT(IN)    :: zeta
      REAL, INTENT(OUT)   :: rig

      REAL :: phi_m, phi_h

      CALL phi_from_zeta(profile_id, zeta, phi_m, phi_h)
      rig = zeta * phi_h / (phi_m * phi_m)
   END SUBROUTINE rig_from_zeta


   SUBROUTINE zeta_from_rig_newton(profile_id, rig_target, zeta)
      !------------------------------------------------------------------------
      ! Solve Ri_g(zeta)=rig_target using Newton iteration.
      ! Uses finite-difference derivative for generality across profiles.
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: profile_id
      REAL, INTENT(IN)    :: rig_target
      REAL, INTENT(OUT)   :: zeta

      INTEGER, PARAMETER :: maxit = 30
      REAL, PARAMETER    :: tol = 1.e-8
      REAL, PARAMETER    :: h = 1.e-5

      INTEGER :: it
      REAL :: z, f0, fp, rigp, rigm

      ! Near-neutral seed (works well for small/moderate stability)
      z = rig_target
      z = MAX(z, -5.0)
      z = MIN(z, 10.0)

      DO it = 1, maxit
         CALL rig_from_zeta(profile_id, z, f0)
         f0 = f0 - rig_target

         CALL rig_from_zeta(profile_id, z + h, rigp)
         CALL rig_from_zeta(profile_id, z - h, rigm)
         fp = (rigp - rigm) / (2.0 * h)

         IF (ABS(fp) < 1.e-12) EXIT

         z = z - f0 / fp

         ! Guard against runaway iterates
         z = MAX(z, -20.0)
         z = MIN(z,  50.0)

         IF (ABS(f0) < tol) EXIT
      END DO

      zeta = z
   END SUBROUTINE zeta_from_rig_newton


   SUBROUTINE zeta_from_rig_safeguarded(profile_id, rig_target, zeta, converged, n_iter)
      !------------------------------------------------------------------------
      ! Solve Ri_g(zeta)=rig_target using safeguarded Newton with bisection.
      !
      ! Strategy:
      !   1) Build a monotone bracket [z_lo, z_hi] by branch.
      !   2) Attempt Newton update when derivative is well-behaved and step
      !      remains inside bracket.
      !   3) Otherwise fall back to bisection (guaranteed contraction).
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: profile_id
      REAL, INTENT(IN)    :: rig_target
      REAL, INTENT(OUT)   :: zeta
      LOGICAL, INTENT(OUT), OPTIONAL :: converged
      INTEGER, INTENT(OUT), OPTIONAL :: n_iter

      INTEGER, PARAMETER :: maxit = 50
      REAL, PARAMETER    :: tol_f = 1.e-8
      REAL, PARAMETER    :: tol_z = 1.e-6
      REAL, PARAMETER    :: h = 1.e-5

      INTEGER :: it
      REAL :: z_lo, z_hi, z
      REAL :: f_lo, f_hi, f, fp, rigp, rigm
      REAL :: z_new, f_new
      LOGICAL :: use_newton, done

      ! Branch-aware initial bracket for typical surface-layer usage
      IF (rig_target <= 0.0) THEN
         z_lo = -20.0
         z_hi =  0.0
      ELSE
         z_lo =  0.0
         z_hi = 50.0
      END IF

      CALL rig_from_zeta(profile_id, z_lo, f_lo)
      CALL rig_from_zeta(profile_id, z_hi, f_hi)
      f_lo = f_lo - rig_target
      f_hi = f_hi - rig_target

      ! If the target is not bracketed, fall back to guarded Newton path.
      IF (f_lo * f_hi > 0.0) THEN
         CALL zeta_from_rig_newton(profile_id, rig_target, z)
         z = MAX(z, z_lo)
         z = MIN(z, z_hi)
         zeta = z
         IF (PRESENT(converged)) converged = .FALSE.
         IF (PRESENT(n_iter)) n_iter = 0
         RETURN
      END IF

      ! Near-neutral seed, then clamp into bracket.
      z = rig_target
      z = MAX(z, z_lo)
      z = MIN(z, z_hi)

      done = .FALSE.
      DO it = 1, maxit
         CALL rig_from_zeta(profile_id, z, f)
         f = f - rig_target

         IF (ABS(f) < tol_f) THEN
            done = .TRUE.
            EXIT
         END IF
         IF (ABS(z_hi - z_lo) < tol_z) THEN
            done = .TRUE.
            EXIT
         END IF

         CALL rig_from_zeta(profile_id, z + h, rigp)
         CALL rig_from_zeta(profile_id, z - h, rigm)
         fp = (rigp - rigm) / (2.0 * h)

         use_newton = (ABS(fp) > 1.e-12)
         IF (use_newton) THEN
            z_new = z - f / fp
            IF (z_new <= z_lo .OR. z_new >= z_hi) use_newton = .FALSE.
         END IF

         IF (.NOT. use_newton) THEN
            z_new = 0.5 * (z_lo + z_hi)
         END IF

         CALL rig_from_zeta(profile_id, z_new, f_new)
         f_new = f_new - rig_target

         ! Bracket update preserves sign change over [z_lo, z_hi].
         IF (f_lo * f_new <= 0.0) THEN
            z_hi = z_new
            f_hi = f_new
         ELSE
            z_lo = z_new
            f_lo = f_new
         END IF

         z = z_new
      END DO

      zeta = z
      IF (PRESENT(converged)) converged = done
      IF (PRESENT(n_iter)) n_iter = it
   END SUBROUTINE zeta_from_rig_safeguarded


   SUBROUTINE fm_fh_from_rig(profile_id, rig, f_m, f_h)
      !------------------------------------------------------------------------
      ! Build Ri-based closure functions from MOST profile by inversion:
      !   zeta <- Ri_g^{-1}(rig)
      !   f_m = 1 / phi_m^2
      !   f_h = 1 / (phi_m * phi_h)
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: profile_id
      REAL, INTENT(IN)    :: rig
      REAL, INTENT(OUT)   :: f_m, f_h

      REAL :: zeta, phi_m, phi_h
      LOGICAL :: ok
      INTEGER :: it

      CALL zeta_from_rig_safeguarded(profile_id, rig, zeta, ok, it)
      IF (.NOT. ok) THEN
         ! Safety fallback for unusual profiles/targets.
         CALL zeta_from_rig_newton(profile_id, rig, zeta)
      END IF
      CALL phi_from_zeta(profile_id, zeta, phi_m, phi_h)

      f_m = 1.0 / (phi_m * phi_m)
      f_h = 1.0 / (phi_m * phi_h)

      f_m = MAX(0.0, MIN(f_m, 1.0))
      f_h = MAX(0.0, MIN(f_h, 1.0))
   END SUBROUTINE fm_fh_from_rig


   SUBROUTINE cbc_critical_ri_limits(b_m, b_h, beta_s, ri_c_ubl_m, ri_c_ubl_h, ri_c_sbl)
      !------------------------------------------------------------------------
      ! Wrapper exposing the parameter-dependent critical Ri limits implied by
      ! the CBC / stable-linear representations.
      !------------------------------------------------------------------------
      REAL, INTENT(IN)  :: b_m, b_h, beta_s
      REAL, INTENT(OUT) :: ri_c_ubl_m, ri_c_ubl_h, ri_c_sbl

      CALL compute_cbc_critical_ri_limits(b_m, b_h, beta_s, ri_c_ubl_m, ri_c_ubl_h, ri_c_sbl)
   END SUBROUTINE cbc_critical_ri_limits

END MODULE module_most_profile_utils
