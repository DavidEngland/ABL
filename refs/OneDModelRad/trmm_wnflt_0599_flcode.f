	MODULE TRMM_WINDOW
	integer ,dimension(11:13):: iwncld
        real,dimension(11:13):: fu_sr,fu_sf
	real,dimension(11:13):: adm, sat_flt_sr,sat_unf_sr
	real,dimension(11:13)::      sat_flt_sf,sat_unf_sf
        real :: 	 sat_flt_r,sat_unf_r, sat_f

	contains
!============================================================
	subroutine trmm_wnflt(COSVZA)
!fu_sr      = FU-liou Radiance in each window band
!fu_sf      = FU-liou Flux in each window band

! sat_flt_sr = emmulation of TRMM FILTERED RADIANCE in each window band
! sat_unf_sr = emmulation of TRMM UNFILTERED RADIANCE in each window band

! sat_flt_r = emmulation of TRMM FILTERED RADIANCE 
! sat_unf_r = emmulation of TRMM UNFILTERED RADIANCE 
! sat_f     = emmulation of TRMM Window Flux

	 sat_flt_r =0.0
	 sat_unf_r =0.0
	 sat_f     =0.0 
	BAND : do ib = 11,13

	adm(ib) = fu_sf(ib) /(3.14159*fu_sr(ib) )
	ff  = TRMM_FULIOU_FILTER(ib,0,iwncld(ib),fu_sr(ib),COSVZA)
	uff = TRMM_FULIOU_FILTER(ib,1,iwncld(ib),fu_sr(ib),COSVZA)

	sat_flt_sr(ib) =  fu_sr(ib) * ff
	sat_unf_sr(ib) = sat_flt_sr(ib) * uff
	sat_flt_r = sat_flt_r + sat_flt_sr(ib)
	sat_unf_r = sat_unf_r + sat_unf_sr(ib)
! Flux method 1 ( SPECTRAL ADM  )
	sat_f     = sat_f     + sat_unf_sr(ib) * (3.14159*adm(ib))


!! FLUX METHOD 2 ( DIFFUSE ANGLE Aprox)
!!	ff  = TRMM_FULIOU_FILTER(ib,0,iwncld(ib),fu_sr(ib),0.6)
!!	uff = TRMM_FULIOU_FILTER(ib,1,iwncld(ib),fu_sr(ib),0.6)
!!	sat_flt_sf(ib) =  fu_sf(ib) * ff
!!	sat_unf_sf(ib) = sat_flt_sf(ib) * uff
!!	sat_f     = sat_f     + sat_unf_sf(ib)
!!!
 
	enddo BAND	

!	print* , 'SCCWN',sum( sat_unf_sr ) / sum( sat_flt_sr )

	return
	end subroutine trmm_wnflt
!============================================================
	real function TRMM_FULIOU_FILTER   &
                   (ib,idir,icld,UNFILTERED_RADIANCE,COSVZA)
	parameter (ibw = 1) ! ibw=1 (844:1227cm-1) 
			    ! ibw=2 (847:1219cm-1) 
!!ib - FULIOU BAND !!
! 11 = 1100:1250 cm-1 Fu-Liou band
! 12 = 980:1100 cm-1  Fu-Liou band
! 13 = 800:980 cm-1   Fu-Liou band

!! idir - Filter or UNfilter
! 0 = FU-liou band to FILTERED BAND
! 1 = FILTERED BAND to TRMM WINDOW 

!! icld - Clear or Cloudy Sky
! 0 = Clear
! 1 = Cloudy 

! UNFILTERED_RADIANCE - Radiance[wm-2sr-1] corresponding 
!                       to the FU-LIOU band "ib" [11-13]

! COSVZA - Cosine of view zenith angle of toa radiance observation [0-1]

	real coefs(4,3,0:1,0:1,2)
	integer ibt(11:13) 
	data ibt/3,2,1/
	data coefs/ &
!! (844:1227cm-1)
      5.137e-01, 5.002e-04,-4.887e-03, 3.724e-05,&
      7.156e-01,-2.305e-04, 9.184e-04, 1.272e-05,&
      6.067e-01, 6.400e-04,-9.066e-03, 3.633e-05,&
      1.382e+00,-8.756e-05, 7.637e-04,-3.077e-06,&
      1.398e+00, 4.567e-04,-1.826e-03,-2.524e-05,&
      1.463e+00,-3.787e-04, 1.142e-03,-4.245e-05,&
      4.946e-01, 3.123e-03,-2.532e-04,-6.721e-05,&
      7.187e-01,-1.150e-03, 2.028e-04, 9.825e-05,&
      6.101e-01,-2.080e-03,-2.730e-03, 2.223e-04,&
      1.384e+00,-3.153e-04, 1.790e-04, 6.385e-06,&
      1.391e+00, 2.233e-03,-3.908e-04,-1.896e-04,&
      1.460e+00, 8.689e-04,-2.329e-04,-1.286e-04,&
!!! (847:1219cm-1)
      5.137e-01, 5.002e-04,-4.887e-03, 3.724e-05,&
      7.156e-01,-2.305e-04, 9.184e-04, 1.272e-05,&
      6.067e-01, 6.400e-04,-9.066e-03, 3.633e-05,&
      1.346e+00, 3.661e-05, 5.024e-04,-1.379e-06,&
      1.398e+00, 4.567e-04,-1.826e-03,-2.524e-05,&
      1.390e+00, 4.193e-04,-8.215e-03, 5.476e-05,&
      4.946e-01, 3.123e-03,-2.532e-04,-6.721e-05,&
      7.187e-01,-1.150e-03, 2.028e-04, 9.825e-05,&
      6.101e-01,-2.080e-03,-2.730e-03, 2.223e-04,&
      1.344e+00, 2.939e-04, 1.753e-04,-1.015e-05,&
      1.391e+00, 2.233e-03,-3.908e-04,-1.896e-04,&
      1.395e+00,-2.624e-03,-2.174e-03, 2.697e-04/

	if ( ib < 11 .or. ib > 13) then
	TRMM_FULIOU_FILTER=1.0000
	return
	endif

	if ( (idir == 0 .or. idir == 1) .and.   &
             (icld == 0 .or. icld == 1) )then
	else
	stop ' ERROR idir or icld NE to 0 or 1'
	endif

	icld_use = icld
	if ( icld_use== 1 )then
	if ( ib == 11 .and. UNFILTERED_RADIANCE >  8.0) icld_use =0
	if ( ib == 12 .and. UNFILTERED_RADIANCE >  8.0) icld_use =0
	if ( ib == 13 .and. UNFILTERED_RADIANCE > 20.0) icld_use =0
	endif

	TRMM_FULIOU_FILTER =    &
      				coefs(1,ibt(ib),idir,icld_use,ibw) +   &
       UNFILTERED_RADIANCE    * coefs(2,ibt(ib),idir,icld_use,ibw) +    &
       COSVZA                 * coefs(3,ibt(ib),idir,icld_use,ibw) +   &
       UNFILTERED_RADIANCE**2 * coefs(4,ibt(ib),idir,icld_use,ibw)  

	return
	end function TRMM_FULIOU_FILTER

	end MODULE TRMM_WINDOW
