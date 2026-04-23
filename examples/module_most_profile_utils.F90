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
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: PROFILE_BD71, PROFILE_HOG88, PROFILE_BH91
   PUBLIC :: phi_from_zeta
   PUBLIC :: rig_from_zeta
   PUBLIC :: zeta_from_rig_newton
   PUBLIC :: fm_fh_from_rig

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

      REAL, PARAMETER :: a_unst = 16.0

      IF (zeta < 0.0) THEN
         ! Businger-Dyer unstable branch
         phi_m = (1.0 - a_unst * zeta)**(-0.25)
         phi_h = (1.0 - a_unst * zeta)**(-0.50)
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

      CALL zeta_from_rig_newton(profile_id, rig, zeta)
      CALL phi_from_zeta(profile_id, zeta, phi_m, phi_h)

      f_m = 1.0 / (phi_m * phi_m)
      f_h = 1.0 / (phi_m * phi_h)

      f_m = MAX(0.0, MIN(f_m, 1.0))
      f_h = MAX(0.0, MIN(f_h, 1.0))
   END SUBROUTINE fm_fh_from_rig

END MODULE module_most_profile_utils
