<!-- filepath: /Users/davidengland/Documents/GitHub/ABL/implementations/McNider_1DBLM_fc_module.f90 -->
MODULE FC_CORRECTION
!-------------------------------------------------------------------------------
! Grid-aware correction module for McNider 1DBLM
! Computes bias ratio B = Ri_g(z_g)/Ri_b and applies neutral-preserving 
! multiplicative damping fc to eddy diffusivities KM, KH
!-------------------------------------------------------------------------------
IMPLICIT NONE

! Default tuning parameters (adjust after validation)
REAL, PARAMETER :: ALPHA_DEFAULT = 1.0      ! Correction strength
REAL, PARAMETER :: P_DEFAULT     = 1.0      ! Dz scaling exponent
REAL, PARAMETER :: Q_DEFAULT     = 2.0      ! Zeta scaling exponent (>=2 for neutral preservation)
REAL, PARAMETER :: DZ_REF        = 10.0     ! Reference layer thickness (m)
REAL, PARAMETER :: ZETA_REF      = 0.5      ! Reference stability
REAL, PARAMETER :: FC_MIN        = 0.2      ! Floor to prevent over-damping
REAL, PARAMETER :: B_THRESH      = 1.05     ! Bias threshold for triggering correction

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE COMPUTE_BIAS_RATIO(Z0, Z1, T0, T1, U0, U1, V0, V1, L, G, B, Z_G)
!-------------------------------------------------------------------------------
! Computes bias ratio B = Ri_g(z_g) / Ri_b for a single layer
! Uses geometric mean height z_g = sqrt(z0*z1) for point Ri evaluation
!
! INPUTS:
!   Z0, Z1     : Layer bounds (m)
!   T0, T1     : Potential temperature at bounds (K)
!   U0, U1     : U-wind at bounds (m/s)
!   V0, V1     : V-wind at bounds (m/s)
!   L          : Monin-Obukhov length (m)
!   G          : Gravitational acceleration (m/s^2)
! OUTPUTS:
!   B          : Bias ratio (dimensionless)
!   Z_G        : Geometric mean height (m)
!-------------------------------------------------------------------------------
    REAL, INTENT(IN)  :: Z0, Z1, T0, T1, U0, U1, V0, V1, L, G
    REAL, INTENT(OUT) :: B, Z_G
    
    REAL :: THETA_REF, DELTA_THETA, DELTA_Z, DELTA_U, DELTA_V
    REAL :: SPEED_DIFF_SQ, RI_B, RI_G_ZG, ZETA_G, PHI_M, PHI_H
    REAL, PARAMETER :: TINY = 1.0E-10
    
    ! Geometric mean height (log-space midpoint)
    Z_G = SQRT(Z0 * Z1)
    
    ! Bulk Richardson number
    THETA_REF   = 0.5 * (T0 + T1)
    DELTA_THETA = T1 - T0
    DELTA_Z     = Z1 - Z0
    DELTA_U     = U1 - U0
    DELTA_V     = V1 - V0
    SPEED_DIFF_SQ = DELTA_U**2 + DELTA_V**2
    
    IF (SPEED_DIFF_SQ .LT. TINY) SPEED_DIFF_SQ = TINY  ! Guard against zero shear
    
    RI_B = (G / THETA_REF) * (DELTA_THETA * DELTA_Z) / SPEED_DIFF_SQ
    
    ! Point Richardson at geometric mean using MOST
    IF (ABS(L) .GT. TINY) THEN
        ZETA_G = Z_G / L
        ! Use existing phi functions from your model (Businger or England-McNider)
        CALL GET_PHI_M(ZETA_G, PHI_M)  ! Replace with your actual function
        CALL GET_PHI_H(ZETA_G, PHI_H)
        RI_G_ZG = ZETA_G * PHI_H / (PHI_M**2)
    ELSE
        ! Neutral or near-neutral fallback
        RI_G_ZG = RI_B
    END IF
    
    ! Bias ratio (guard against division by zero)
    IF (ABS(RI_B) .GT. TINY) THEN
        B = RI_G_ZG / RI_B
    ELSE
        B = 1.0  ! Neutral default
    END IF
    
    ! Ensure B >= 1 (concave-down assumption)
    IF (B .LT. 1.0) B = 1.0

END SUBROUTINE COMPUTE_BIAS_RATIO

!-------------------------------------------------------------------------------
SUBROUTINE GET_PHI_M(ZETA, PHI_M)
!-------------------------------------------------------------------------------
! Placeholder for momentum stability function phi_m(zeta)
! REPLACE THIS with your actual implementation (Businger, England-McNider, etc.)
!-------------------------------------------------------------------------------
    REAL, INTENT(IN)  :: ZETA
    REAL, INTENT(OUT) :: PHI_M
    
    ! Example: Businger stable form
    IF (ZETA .GT. 0.0) THEN
        PHI_M = 1.0 + 5.3 * ZETA
    ELSE
        PHI_M = (1.0 - 19.0 * ZETA)**0.25
    END IF
END SUBROUTINE GET_PHI_M

!-------------------------------------------------------------------------------
SUBROUTINE GET_PHI_H(ZETA, PHI_H)
!-------------------------------------------------------------------------------
! Placeholder for heat stability function phi_h(zeta)
! REPLACE THIS with your actual implementation
!-------------------------------------------------------------------------------
    REAL, INTENT(IN)  :: ZETA
    REAL, INTENT(OUT) :: PHI_H
    
    ! Example: Businger stable form
    IF (ZETA .GT. 0.0) THEN
        PHI_H = 0.74 + 8.0 * ZETA
    ELSE
        PHI_H = SQRT(1.0 - 11.6 * ZETA)
    END IF
END SUBROUTINE GET_PHI_H

!-------------------------------------------------------------------------------
SUBROUTINE COMPUTE_FC(B, ZETA, DELTA_Z, FC, ALPHA, P, Q)
!-------------------------------------------------------------------------------
! Computes correction factor fc using exponential form
! fc = exp[-alpha * (B-1) * (dz/dz_ref)^p * (zeta/zeta_ref)^q]
!
! INPUTS:
!   B          : Bias ratio
!   ZETA       : Nondimensional height z/L
!   DELTA_Z    : Layer thickness (m)
!   ALPHA, P, Q: Optional tuning parameters (use defaults if not provided)
! OUTPUT:
!   FC         : Correction factor (0.2 <= fc <= 1.0)
!-------------------------------------------------------------------------------
    REAL, INTENT(IN)  :: B, ZETA, DELTA_Z
    REAL, INTENT(OUT) :: FC
    REAL, INTENT(IN), OPTIONAL :: ALPHA, P, Q
    
    REAL :: ALPHA_USE, P_USE, Q_USE, EXPONENT
    
    ! Use provided parameters or defaults
    ALPHA_USE = ALPHA_DEFAULT
    P_USE     = P_DEFAULT
    Q_USE     = Q_DEFAULT
    IF (PRESENT(ALPHA)) ALPHA_USE = ALPHA
    IF (PRESENT(P))     P_USE     = P
    IF (PRESENT(Q))     Q_USE     = Q
    
    ! Apply threshold: no correction if bias is small
    IF (B .LE. B_THRESH) THEN
        FC = 1.0
        RETURN
    END IF
    
    ! Compute exponent
    EXPONENT = -ALPHA_USE * (B - 1.0) * &
               (DELTA_Z / DZ_REF)**P_USE * &
               (ABS(ZETA) / ZETA_REF)**Q_USE
    
    FC = EXP(EXPONENT)
    
    ! Apply floor
    FC = MAX(FC, FC_MIN)

END SUBROUTINE COMPUTE_FC

!-------------------------------------------------------------------------------
SUBROUTINE APPLY_FC_TO_PROFILE(Z, KM, KH, L, T, U, V, KZ, G, FIRST_CALL)
!-------------------------------------------------------------------------------
! Main driver: applies fc correction to entire K-profile
! Call this at the end of PROFILEK subroutine in your model
!
! INPUTS/OUTPUTS:
!   Z(KZ)      : Heights (m)
!   KM(KZ)     : Momentum diffusivity (m^2/s) - MODIFIED
!   KH(KZ)     : Heat diffusivity (m^2/s) - MODIFIED
!   L          : Monin-Obukhov length (m)
!   T(KZ)      : Potential temperature (K)
!   U(KZ)      : U-wind (m/s)
!   V(KZ)      : V-wind (m/s)
!   KZ         : Number of vertical levels
!   G          : Gravitational acceleration (m/s^2)
!   FIRST_CALL : Logical flag for initialization/diagnostics
!-------------------------------------------------------------------------------
    INTEGER, INTENT(IN)    :: KZ
    REAL, INTENT(INOUT)    :: KM(KZ), KH(KZ)
    REAL, INTENT(IN)       :: Z(KZ), L, T(KZ), U(KZ), V(KZ), G
    LOGICAL, INTENT(INOUT) :: FIRST_CALL
    
    INTEGER :: I
    REAL :: Z0, Z1, T0, T1, U0, U1, V0, V1, DELTA_Z
    REAL :: B, Z_G, FC, ZETA, KM_OLD, KH_OLD
    
    ! Open diagnostic file on first call
    IF (FIRST_CALL) THEN
        OPEN(UNIT=99, FILE='FC_DIAGNOSTICS.dat', STATUS='REPLACE')
        WRITE(99,'(A)') 'LEVEL  Z_G(m)  B      FC     KM_OLD  KM_NEW  KH_OLD  KH_NEW'
        FIRST_CALL = .FALSE.
    END IF
    
    ! Loop over layers (apply to interior levels only)
    DO I = 2, KZ-1
        Z0 = Z(I-1)
        Z1 = Z(I)
        T0 = T(I-1)
        T1 = T(I)
        U0 = U(I-1)
        U1 = U(I)
        V0 = V(I-1)
        V1 = V(I)
        DELTA_Z = Z1 - Z0
        
        ! Compute bias ratio B for this layer
        CALL COMPUTE_BIAS_RATIO(Z0, Z1, T0, T1, U0, U1, V0, V1, L, G, B, Z_G)
        
        ! Compute nondimensional height for this layer
        IF (ABS(L) .GT. 1.0E-10) THEN
            ZETA = Z_G / L
        ELSE
            ZETA = 0.0
        END IF
        
        ! Compute correction factor fc
        CALL COMPUTE_FC(B, ZETA, DELTA_Z, FC)
        
        ! Store old values for diagnostics
        KM_OLD = KM(I)
        KH_OLD = KH(I)
        
        ! Apply correction
        KM(I) = KM(I) * FC
        KH(I) = KH(I) * FC
        
        ! Write diagnostics every 10 levels (adjust as needed)
        IF (MOD(I, 10) .EQ. 0) THEN
            WRITE(99,'(I5,7F8.3)') I, Z_G, B, FC, KM_OLD, KM(I), KH_OLD, KH(I)
        END IF
    END DO

END SUBROUTINE APPLY_FC_TO_PROFILE

END MODULE FC_CORRECTION