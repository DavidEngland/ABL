module most_similarity
  !--------------------------------------------------------------------
  !  Monin–Obukhov similarity functions φ_m, φ_h
  !  Businger–Dyer type, tuned for use in WRF-style surface layer code.
  !
  !  Default set: Dyer (1974) / Brutsaert (1982)
  !    κ      = 0.40
  !    Pr0inv = 1.00   (φ_h neutral)
  !    b_m    = 16.0   (unstable)
  !    b_h    = 16.0   (unstable)
  !    β_m    = 5.0    (stable)
  !    β_h    = 5.0    (stable)
  !
  !  All reals are double precision (kind=dp).
  !--------------------------------------------------------------------
  implicit none
  private

  integer, parameter :: dp = selected_real_kind(15, 300)

  ! --- Public API
  public :: kappa_most
  public :: phi_m, phi_h
  public :: set_param_dyer1974
  public :: set_param_businger1971
  public :: set_param_hogstrom1996

  ! --- Core parameters (module variables, can be changed at runtime)
  real(dp) :: kappa    = 0.40_dp
  real(dp) :: pr0inv   = 1.00_dp   ! neutral φ_h
  real(dp) :: b_m      = 16.0_dp   ! unstable
  real(dp) :: b_h      = 16.0_dp   ! unstable
  real(dp) :: beta_m   = 5.0_dp    ! stable
  real(dp) :: beta_h   = 5.0_dp    ! stable

  ! --- ζ clipping (for numerical robustness)
  real(dp), parameter :: zeta_min = -5.0_dp   ! lower bound (strongly unstable)
  real(dp), parameter :: zeta_max =  5.0_dp   ! upper bound (very stable)

contains

  !====================================================================
  function kappa_most() result(k)
    ! Return von Kármán constant currently in use
    real(dp) :: k
    k = kappa
  end function kappa_most

  !====================================================================
  subroutine set_param_dyer1974()
    ! Dyer (1974) / Brutsaert (1982) canonical set
    kappa   = 0.40_dp
    pr0inv  = 1.00_dp
    b_m     = 16.0_dp
    b_h     = 16.0_dp
    beta_m  = 5.0_dp
    beta_h  = 5.0_dp
  end subroutine set_param_dyer1974

  !--------------------------------------------------------------------
  subroutine set_param_businger1971()
    ! Businger et al. (1971) Kansas set
    ! κ ≈ 0.35, Pr0inv ≈ 1.35, b_m = 15, b_h = 9, β_m = 4.7, β_h = 6.35
    kappa   = 0.35_dp
    pr0inv  = 1.35_dp
    b_m     = 15.0_dp
    b_h     =  9.0_dp
    beta_m  =  4.7_dp
    beta_h  =  6.35_dp
  end subroutine set_param_businger1971

  !--------------------------------------------------------------------
  subroutine set_param_hogstrom1996()
    ! Högström (1996) representative set
    ! κ = 0.40, Pr0inv ≈ 0.95–1.05 (take 1.00 here), b_m = 19, b_h = 11.6,
    ! β_m ≈ 5.3, β_h ≈ 5.3
    kappa   = 0.40_dp
    pr0inv  = 1.00_dp
    b_m     = 19.0_dp
    b_h     = 11.6_dp
    beta_m  =  5.3_dp
    beta_h  =  5.3_dp
  end subroutine set_param_hogstrom1996

  !====================================================================
  pure function phi_m(zeta) result(phi)
    ! Non-dimensional wind-speed gradient function φ_m(ζ)
    !
    ! Unstable (ζ < 0):
    !   φ_m = (1 - b_m ζ)^(-1/4)
    !
    ! Stable (ζ >= 0):
    !   φ_m = 1 + β_m ζ
    !
    real(dp), intent(in) :: zeta
    real(dp)             :: phi
    real(dp)             :: zz, arg

    ! Clip ζ for robustness
    zz = max(zeta_min, min(zeta_max, zeta))

    if (zz < 0.0_dp) then
       ! Unstable branch
       arg = 1.0_dp - b_m * zz
       ! Guard against negative/zero argument due to extreme ζ
       arg = max(arg, 1.0e-6_dp)
       phi = arg**(-0.25_dp)
    else
       ! Stable branch
       phi = 1.0_dp + beta_m * zz
    end if
  end function phi_m

  !====================================================================
  pure function phi_h(zeta) result(phi)
    ! Non-dimensional temperature gradient function φ_h(ζ)
    !
    ! Unstable (ζ < 0):
    !   φ_h = Pr0inv * (1 - b_h ζ)^(-1/2)
    !
    ! Stable (ζ >= 0):
    !   φ_h = Pr0inv + β_h ζ
    !
    real(dp), intent(in) :: zeta
    real(dp)             :: phi
    real(dp)             :: zz, arg

    ! Clip ζ for robustness
    zz = max(zeta_min, min(zeta_max, zeta))

    if (zz < 0.0_dp) then
       ! Unstable branch
       arg = 1.0_dp - b_h * zz
       arg = max(arg, 1.0e-6_dp)
       phi = pr0inv * arg**(-0.5_dp)
    else
       ! Stable branch
       phi = pr0inv + beta_h * zz
    end if
  end function phi_h

end module most_similarity