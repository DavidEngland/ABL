!==============================================================================
! Driver: exact vs CBC/Gegenbauer approximation errors for unstable MOST branch
!
! Build example:
!   gfortran -ffree-form -o driver_cbc_gegenbauer_errors \
!     module_cbc_legendre_most.F90 driver_cbc_gegenbauer_errors.F90
!
! Run:
!   ./driver_cbc_gegenbauer_errors
!==============================================================================

PROGRAM driver_cbc_gegenbauer_errors
   USE module_cbc_legendre_most, ONLY: phi_h_cbc_unstable, phi_m_gegenbauer_unstable, &
                                        central_binomial_scaled, legendre_p_even_zero,  &
                                        compute_cbc_critical_ri_limits
   IMPLICIT NONE

   INTEGER, PARAMETER :: nz = 8, nn = 5
   REAL, PARAMETER :: b_h = 16.0, b_m = 16.0, beta_s = 5.0
   REAL, DIMENSION(nz) :: zeta_vals
   INTEGER, DIMENSION(nn) :: n_terms_list

   REAL :: zeta, phi_h_exact, phi_m_exact, phi_h_approx, phi_m_approx
   REAL :: err_h, err_m
   REAL :: ri_c_ubl_m, ri_c_ubl_h, ri_c_sbl
   INTEGER :: i, j, n

   zeta_vals = (/ -0.001, -0.005, -0.010, -0.020, -0.050, -0.100, -0.200, -0.400 /)
   n_terms_list = (/ 4, 6, 8, 12, 16 /)

   WRITE(*,'(A)') '============================================================='
   WRITE(*,'(A)') 'CBC/Gegenbauer Error Driver (Unstable MOST branch)'
   WRITE(*,'(A)') '============================================================='
   WRITE(*,'(A,F8.3,A,F8.3,A,F8.3)') 'Parameters: b_m=', b_m, ', b_h=', b_h, ', beta=', beta_s

   CALL compute_cbc_critical_ri_limits(b_m, b_h, beta_s, ri_c_ubl_m, ri_c_ubl_h, ri_c_sbl)
   WRITE(*,'(A,F10.6)') 'Ri_c^UBL,m = ', ri_c_ubl_m
   WRITE(*,'(A,F10.6)') 'Ri_c^UBL,h = ', ri_c_ubl_h
   WRITE(*,'(A,F10.6)') 'Ri_c^SBL   = ', ri_c_sbl
   WRITE(*,'(A)') '-------------------------------------------------------------'

   WRITE(*,'(A)') 'Legendre-equator check: P_{2n}(0) vs (-1)^n * binom(2n,n)/4^n'
   WRITE(*,'(A)') '   n        P_{2n}(0)              recurrence-value'
   DO n = 0, 8
      WRITE(*,'(I4,2X,ES20.10,2X,ES20.10)') n, legendre_p_even_zero(n),  &
           (-1.0)**REAL(n) * central_binomial_scaled(n)
   END DO
   WRITE(*,'(A)') '-------------------------------------------------------------'

   WRITE(*,'(A)') 'Table A: phi_h exact vs CBC series'
   WRITE(*,'(A)') ' zeta      N    phi_h_exact      phi_h_approx     rel_err'
   DO i = 1, nz
      zeta = zeta_vals(i)
      phi_h_exact = (1.0 - b_h * zeta)**(-0.5)
      DO j = 1, nn
         phi_h_approx = phi_h_cbc_unstable(zeta, b_h, n_terms_list(j))
         err_h = rel_err(phi_h_approx, phi_h_exact)
         WRITE(*,'(F8.3,1X,I4,2X,ES14.6,2X,ES14.6,2X,ES12.4)') zeta, n_terms_list(j), &
              phi_h_exact, phi_h_approx, err_h
      END DO
      WRITE(*,'(A)') ' '
   END DO

   WRITE(*,'(A)') 'Table B: phi_m exact vs Gegenbauer-equator series'
   WRITE(*,'(A)') ' zeta      N    phi_m_exact      phi_m_approx     rel_err'
   DO i = 1, nz
      zeta = zeta_vals(i)
      phi_m_exact = (1.0 - b_m * zeta)**(-0.25)
      DO j = 1, nn
         phi_m_approx = phi_m_gegenbauer_unstable(zeta, b_m, n_terms_list(j))
         err_m = rel_err(phi_m_approx, phi_m_exact)
         WRITE(*,'(F8.3,1X,I4,2X,ES14.6,2X,ES14.6,2X,ES12.4)') zeta, n_terms_list(j), &
              phi_m_exact, phi_m_approx, err_m
      END DO
      WRITE(*,'(A)') ' '
   END DO

   WRITE(*,'(A)') 'Done.'

CONTAINS

   PURE REAL FUNCTION rel_err(v, vref)
      REAL, INTENT(IN) :: v, vref
      REAL, PARAMETER :: tiny = 1.0e-20
      rel_err = ABS(v - vref) / MAX(ABS(vref), tiny)
   END FUNCTION rel_err

END PROGRAM driver_cbc_gegenbauer_errors
