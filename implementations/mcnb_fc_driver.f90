PROGRAM MCNB_FC_DRIVER
USE FC_CORRECTION
IMPLICIT NONE

! Simple demo: build a small synthetic profile, compute fc and print CSV
INTEGER, PARAMETER :: KZ = 6
REAL :: Z(KZ), U(KZ), V(KZ), T(KZ), KM(KZ), KH(KZ)
REAL :: L, G
INTEGER :: i, KZp

! Physical constants
G = 9.81
L = 50.0    ! demo Obukhov length (m)

! Synthetic profile (simple demo values)
Z = (/ 1.0, 2.0, 5.0, 10.0, 20.0, 40.0 /)
DO i = 1, KZ
    U(i) = 2.0 + 0.05*Z(i)     ! mild shear
    V(i) = 0.5
    T(i) = 270.0 + 0.1*Z(i)    ! small lapse
    KM(i) = 0.1                ! background K
    KH(i) = 0.05
END DO

! Disable diagnostic file output by default (notebook will parse stdout)
CALL SET_DIAGNOSTICS(.FALSE.)

! Apply correction (this will modify KM, KH in-place)
CALL APPLY_FC_TO_PROFILE(Z, KM, KH, L, T, U, V, KZ, G, FIRST_CALL=.TRUE.)

! Print CSV header to stdout for notebook to capture
PRINT *, 'level,z,z_g,km,k_h'
DO i = 2, KZ-1
    PRINT '(I0,1X,F6.2,1X,F6.2,1X,F8.4,1X,F8.4)', i, Z(i), SQRT(Z(i-1)*Z(i)), KM(i), KH(i)
END DO

END PROGRAM MCNB_FC_DRIVER
