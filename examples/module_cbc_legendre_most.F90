!==============================================================================
! CBC / Legendre / Gegenbauer utilities for MOST stability functions
!
! Purpose:
!   - Fast recurrence-based evaluation of central-binomial coefficients.
!   - Equatorial Legendre/Gegenbauer values linked to MOST unstable branches.
!   - Parameter-dependent critical Ri limits for unstable and stable branches.
!
! Notes:
!   - phi_h on the unstable branch uses the CBC / Legendre-equator expansion.
!   - phi_m on the unstable branch uses the Gegenbauer (lambda=1/4) expansion.
!   - The series are most useful near neutrality; direct power laws remain robust
!     deeper into the unstable branch.
!==============================================================================

MODULE module_cbc_legendre_most
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: central_binomial_scaled
   PUBLIC :: legendre_p_even_zero
   PUBLIC :: gegenbauer_even_zero
   PUBLIC :: legendre_p_value
   PUBLIC :: phi_h_cbc_unstable
   PUBLIC :: phi_m_gegenbauer_unstable
   PUBLIC :: compute_cbc_critical_ri_limits

CONTAINS

   PURE REAL FUNCTION central_binomial_scaled(n)
      !------------------------------------------------------------------------
      ! Return binom(2n,n) / 4^n using a stable recurrence.
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: n
      INTEGER :: k

      central_binomial_scaled = 1.0
      DO k = 1, n
         central_binomial_scaled = central_binomial_scaled * REAL(2 * k - 1) / REAL(2 * k)
      END DO
   END FUNCTION central_binomial_scaled


   PURE REAL FUNCTION legendre_p_even_zero(n)
      !------------------------------------------------------------------------
      ! Return P_{2n}(0) = (-1)^n * binom(2n,n) / 4^n.
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: n

      legendre_p_even_zero = central_binomial_scaled(n)
      IF (MOD(n, 2) == 1) legendre_p_even_zero = -legendre_p_even_zero
   END FUNCTION legendre_p_even_zero


   PURE REAL FUNCTION gegenbauer_even_zero(n, lambda)
      !------------------------------------------------------------------------
      ! Return C_{2n}^{(lambda)}(0) = (-1)^n Gamma(n+lambda)/(Gamma(lambda)n!).
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN)    :: lambda
      INTEGER :: k

      gegenbauer_even_zero = 1.0
      DO k = 1, n
         gegenbauer_even_zero = gegenbauer_even_zero * (lambda + REAL(k - 1)) / REAL(k)
      END DO
      IF (MOD(n, 2) == 1) gegenbauer_even_zero = -gegenbauer_even_zero
   END FUNCTION gegenbauer_even_zero


   PURE REAL FUNCTION legendre_p_value(n, x)
      !------------------------------------------------------------------------
      ! Three-term recurrence for P_n(x).
      !------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN)    :: x
      INTEGER :: k
      REAL :: p_nm2, p_nm1, p_n

      IF (n == 0) THEN
         legendre_p_value = 1.0
         RETURN
      END IF

      IF (n == 1) THEN
         legendre_p_value = x
         RETURN
      END IF

      p_nm2 = 1.0
      p_nm1 = x
      DO k = 2, n
         p_n = (REAL(2 * k - 1) * x * p_nm1 - REAL(k - 1) * p_nm2) / REAL(k)
         p_nm2 = p_nm1
         p_nm1 = p_n
      END DO

      legendre_p_value = p_nm1
   END FUNCTION legendre_p_value


   PURE REAL FUNCTION phi_h_cbc_unstable(zeta, b_h, n_terms)
      !------------------------------------------------------------------------
      ! CBC / Legendre-equator expansion for phi_h = (1 - b_h zeta)^(-1/2)
      ! on the unstable branch zeta <= 0.
      !------------------------------------------------------------------------
      REAL, INTENT(IN)    :: zeta, b_h
      INTEGER, INTENT(IN) :: n_terms
      INTEGER :: n
      REAL :: eta, term

      IF (zeta >= 0.0) THEN
         phi_h_cbc_unstable = (1.0 - b_h * zeta)**(-0.5)
         RETURN
      END IF

      eta = -zeta
      phi_h_cbc_unstable = 1.0
      term = 1.0

      DO n = 0, n_terms - 2
         term = term * (-(REAL(2 * n + 1)) / REAL(2 * n + 2)) * (b_h * eta)
         phi_h_cbc_unstable = phi_h_cbc_unstable + term
      END DO
   END FUNCTION phi_h_cbc_unstable


   PURE REAL FUNCTION phi_m_gegenbauer_unstable(zeta, b_m, n_terms)
      !------------------------------------------------------------------------
      ! Gegenbauer-equator expansion for phi_m = (1 - b_m zeta)^(-1/4)
      ! on the unstable branch zeta <= 0.
      !------------------------------------------------------------------------
      REAL, INTENT(IN)    :: zeta, b_m
      INTEGER, INTENT(IN) :: n_terms
      INTEGER :: n
      REAL :: eta, term

      IF (zeta >= 0.0) THEN
         phi_m_gegenbauer_unstable = (1.0 - b_m * zeta)**(-0.25)
         RETURN
      END IF

      eta = -zeta
      phi_m_gegenbauer_unstable = 1.0
      term = 1.0

      DO n = 0, n_terms - 2
         term = term * (-(REAL(4 * n + 1)) / REAL(4 * n + 4)) * (b_m * eta)
         phi_m_gegenbauer_unstable = phi_m_gegenbauer_unstable + term
      END DO
   END FUNCTION phi_m_gegenbauer_unstable


   SUBROUTINE compute_cbc_critical_ri_limits(b_m, b_h, beta_s, ri_c_ubl_m, ri_c_ubl_h, ri_c_sbl)
      !------------------------------------------------------------------------
      ! Natural critical Ri limits implied by the parameterization:
      !   unstable momentum  : Ri_c = -1 / b_m
      !   unstable heat      : Ri_c = -1 / b_h
      !   stable linear branch: Ri_c = 1 / beta_s
      !------------------------------------------------------------------------
      REAL, INTENT(IN)  :: b_m, b_h, beta_s
      REAL, INTENT(OUT) :: ri_c_ubl_m, ri_c_ubl_h, ri_c_sbl

      ri_c_ubl_m = -1.0 / b_m
      ri_c_ubl_h = -1.0 / b_h
      ri_c_sbl   =  1.0 / beta_s
   END SUBROUTINE compute_cbc_critical_ri_limits

END MODULE module_cbc_legendre_most
