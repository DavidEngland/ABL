      MODULE RadParams
!## PARAMETERS
      integer, parameter :: nvx = 100
      integer, parameter :: nv1x = nvx + 1 
      integer, parameter :: ndfsx = nvx, mdfsx = nvx + 1 
      integer, parameter :: ndfs4x = 4*ndfsx, mbx = 18, mbsx = 6 
      integer, parameter :: mbirx = 12 , ncx=8, mby = 10 
      integer, parameter :: naer=21, nrh=8 , mxat=7 ,mxac=3

!#INPUTS
      logical, save :: fourssl,foursir
      integer, save :: nv,nv1,mb,mbs,mbir,nc,ndfs,mdfs,ndfs4,ndfs2

      integer       :: iaform,n_atau,ivd,ifg,nac
      integer,  dimension(mxac) :: itps
      real,dimension(mxat,mxac) :: a_wlis,a_taus
      real,dimension(nvx,mxac)  :: aprofs

      logical ,save :: edding, quadra, hemisp, mquadr

! naer ::MAX # of aerosol types( original 18, now add NAAPS smoke and dust
! mxat ::MAX # of wavelength dependent aerosol optical depths
! iaform ::Aerosol  1) tau visible only  2) mfrsr 3) wavelength dependent
! n_atau :: # of aerosol optical depths 'a_tau' at 'a_wli' wavelengths
! mxat   :: MAX # aerosol Tau / wavelengths 
! ifg    :: Aerosol Humidity Dependence 0)Taken from profile 1-8) Fixed 1: 0%; 2: 50%;  3: 70%; 4: 80%; 5:90%; 6: 95%; 7: 98%; and 8: 99%

! itps    :: 18 aerosol types ( maritime continental urban , 0.5um dust 1.0 2.0 4.0 8.0 insoluble soluble soot sea_salt_acc sea_salt_coarse dust_nuc dust_acc dust_coarse dust_transported sulfate

!  aprofs:: Aerosol Profile(s) 
!  nac  :: # Aerosol Constiuents    
!  mxac :: MAX # Aerosol Constiuents 

!  ************************************************************************
!  Aerosol modification.  atau or ataux is column optical depth, naer is
!  the number of aerosol types in the "block data aerosol" tables (see 
!  aerosol subroutines and block data).  nrh is the
!  number of relative humidities in the block data tables, iaform 
!  (aerosol formulation to use) is 1 for CERES, 2 for CAGEX.  You need
!  to set the aerosol taus in your calling program.  You also need to set
!  itp (aerosol type 1, 2, or 3, for maritime, continental, urban) and ifg,
!  the RH flag.  
!
!  ifg=0 will distribute according to RH in layer, ifg=1 - 
!  ifg=8 will distribute according to a set RH:  1: 0%; 2: 50%;
!             3: 70%; 4: 80%; 5:90%; 6: 95%; 7: 98%; and 8: 99%.
!
!  ivd is the vertical distribution flag, used to determine the vertical
!  distribution of aerosol optical depth. "0" corresponds to the Spinhirne
!  distribution.  Others may be added.
!
!!  ipr is the aerosol properties flag.  If set to zero, aerosol scattering 
!!  properties will be apportioned into the 18 Fu-Liou bands normally.  If 
!!  set to one, the first SW band will be divided into 10 sub-intervals, and 
!!  scattering properties for those intervals will be used.  The other 17 
!!  bands will be dealt with normally. 
!  ************************************************************************

      END MODULE RadParams
