! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -c 6-24-98 (1)
      subroutine rad ( as,as1, u0, ss, pts, ee, ur, ITriErr )
! +++ PAUSE FIX      subroutine rad ( as, u0, ss, pts, ee , ur )
! 6-24-98 (1)
! **********************************************************************
! MODIFIED 0598 Tho compute Window RADIANCES ( ONLY WITH 2-Stream OPTION)
! MODIFIED 0598 A) Outputs Window (800-1250cm-1) Flux Profiles., (F. Rose) 
!		B) Includes new Window K's from Dave Kratz,  (F. Rose) 
! 		C) H2O continuum absorbtion extended to (5-2200cm-1) 
!		   according to CKD_2.1, (F. Rose) 
! MODIFIED 4/2/97  to include CKD continuum for LW and associated 
!                  parameterizations, (F. Rose) 
! MODIFIED 4/1/97 to apportion 1st SW band into 10 sub-intervals 
!                  for aerosols. Removed commented out sections. 
!                  from older versions (T. Alberta)
! MODIFIED 10/26/96 to include aerosols and direct/diffuse. (T. Alberta)
!
! MODIFIED 10/25/96 - changed logicals for 2 and 4 stream
!                     options (T. Alberta)
! MODIFIED 10/23/96 - many cosmetic changes. (T. Alberta)
! MODIFIED 10/23/96 - many cosmetic changes. (T. Alberta)
!
! MODIFIED 9/96 to include variable levels after compilation. (F. Rose)
!
! In this radiation scheme,  six  and  12 bands are selected for solar 
! and thermal IR regions, respectively. The spectral division is below: 
! 0.2 - 0.7 um, 0.7 - 1.3 um, 1.3 - 1.9 um, 1.9 - 2.5 um, 2.5 -3.5 um,
! 3.5 - 4.0 um, and 2200 - 1900 cm**-1, 1900 - 1700 cm**-1, 1700 -1400
! cm**-1,  1400 - 1250 cm**-1,  1250 - 1100 cm**-1, 1100 - 980 cm**-1,
! 980 - 800 cm**-1,  800 - 670 cm**-1,  670 - 540 cm**-1, 540 - 400 cm
! **-1,  400 - 280 cm**-1,  280 - 0 cm**-1,  where  the index  for the
! spectral band ( ib = 1, 2, ..., 18 ) is defined.
!
!                                                 ***********
!                                                 * Common  *
!             **********************              * block   *  ******
!             *  INPUT PARAMETERS  *              * or call *  *TYPE*
!             **********************              ***********  ******
!   as(mbs)   Surface albedo, mbs = 6                CALL       REAL
!   as1(10)   Surface albedo in 10 vis_subbands      CALL       REAL
!             (Spectral reflectances for SW bands)
!   u0        cosine of solar zenith angle           CALL       REAL
!   ss        solar constant (W/m**2)                CALL       REAL
!   pts       surface temperature (K)                CALL       REAL
!   ee(mbir)  IR surface emissivity, mbir=12         CALL       REAL
!   pp(nv1x)  atmospheric pressure (millibars)      atmos       REAL
!   pt(nv1x)  atmospheric temperature (K)           atmos       REAL
!   ph(nv1x)  water vapor mixing ratio (kg/kg)      atmos       REAL
!   po(nv1x)  ozone mixing ratio (kg/kg)            atmos       REAL
!   pre(nvx)  effective radius of water cloud (um)  clouds      REAL
!   plwc(nvx) liquid water content (g/m**3)         clouds      REAL
!   pde(nvx)  effective diameter of ice cloud (um)  clouds      REAL
!   piwc(nvx) ice water content (g/m**3)            clouds      REAL
!   prwc(nvx) rain water content (g/m**3)           rains       REAL
!   pgwc(nvx) graupel water content (g/m**3)        graups      REAL
!   umco2     concentration of CO2 (ppmv)           umcon       REAL
!   umch4     concentration of CH4 (ppmv)           umcon       REAL
!   umn2o     concentration of N2O (ppmv)           umcon       REAL
!   fourssl   if true, run four stream for solar    tsfslog    LOGICAL
!             if false, two stream is run
!   foursir   if true, run four stream for IR       tsfslog    LOGICAL
!             if false, run two/four stream for IR  
!   nv        number of LAYERS                      levels     INTEGER
!   itp       aerosol type (1, 2, 3 - see header)   aer_tau    INTEGER  
!   ivd       Aerosol optical depth vertical        aer_tau    INTEGER
!             distribution (see header)
!   ipr       Aerosol properties flag (see header)  aer_tau    INTEGER  
!   irobckd   LW continuum option                   irobckd    INTEGER 
!             1: Roberts continuum
!             2: Exact CKD_2.1 continuum
!             3: No continuum
!             4: Parameterized CKD_2.1 continuum
!   a_tau(5)  MFRSR-derived aerosol optical depths  aer_tau     REAL
!             (For CAGEX only)
!   a_taux    Column aerosol optical depth          aer_tau     REAL
!             (For CERES only)
!
! NEW FOR 0598  NEW K's For WINDOW band 11,12,13 
!
!   idkfr     Window K's Option			   /dkfrwn/     INTEGER  
!	      0: orig Fu_liou 
!	      1:Kratz Window K's :: SLOW but ACCURATE
!             2: Hybred          :: FAST   SUGGEST idkfr=2
!
! H20 Continuum absorbtion now for entire LW (5-2200cm-1)
!    iwtas    1:Linear				    /cont_tas/   INTEGER
!             2:log ,
!             3:linear w/Plank wgt   :: SUGGEST iwtas=3
!             4: Log w/plank wgt 
! Fu suggests use of linear averaging of continuum spectral Taus weighted by
! plank function 
!
! cfc_conc(3)  Window K's now include CFC's         /cfcs/      REAL
!	      ( F11,F12,F22 ) in ppv
!
! Suggested data cfc_conc/0.268e-09 , 0.503e-09 ,0.105e-09/
!  
! NEW FOR 0698 RADIANCE
! Input:
!	ur - Cosine of View Zenith Angle 
! Output:: 
!	fiurt(nv1x)  - TOTAL LW (0-2200cm-1) radiance [Wm-2sr-1]
!	fiurw(nv1x)  - WINDOW (800-1250cm-1) radiance [Wm-2sr-1]
! 
!   Check these common blocks in this program and in the header file
!   for proper format.
!   
!
! Note:  (1)  as(mbs) and ee(mbir) consider the substantial wavelength
!             dependence of surface albedos and emissivities.
!        (2)  For CO2, CH4 and N2O, uniform mixing is assumed  through
!             the atmosphere. The  concentrations  can be changed
!             through 'common /umcon/ umco2, umch4, umn2o '.
!        (3)  nvx, nv1x, ndfsx, mdfsx, ndfs4x, mbx, mbsx, mbirx,  
!             and  ncx  are given through 'rad_0598.h'.  These
!             variables are used for dimensioning only!
!        (4)  nv1 and 1 are the surface and top levels, respectively.
!
!                       **********************
!                       *  OUTPUT PARAMETERS  *
!                       **********************
!              fds(nv1x)   downward solar flux ( W / m ** 2 )
!              fus(nv1x)   upward solar flux ( W / m **2 )
!              dts(nvx)    solar heating rate ( K / day )
!              fdir(nv1x)  downward IR flux ( W / m ** 2 )
!              fuir(nv1x)  upward IR flux ( W / m **2 )
!              dtir(nvx)   IR heating rate ( K / day )
!              fd(nv1x)    downward net flux ( W / m ** 2 )
!              fu(nv1x)    upward net flux ( W / m **2 )
!              fdsdr(nv1x) downward direct SW flux ( W / m **2 )
!              fdsdf(nv1x) downward diffuse SW flux ( W / m **2 )
!              dt(nvx)     net heating rate ( K / day )
!              fdwn(nv1x)  downward WINDOW flux ( W / m ** 2 )
!              fuwn(nv1x)  upward WINDOW flux ( W / m **2 )
!
! Note:  Solar, IR, and net represent 0.2 - 0.4 um, 2200 - 0 cm**-1,
!        and  entire spectral regions, respectively.
! 	 Window represents (800-1250cm-1)
! 6-24-98 (2)	
!                  ********************************* 
!                  * OUTPUT FOR IR WINDOW RADIANCE *
!                  *********************************
!              fiurw(nv1) upward radiance at ur in IR window (W/m**2/Sr)
!              fiurt(nv1) upward radiance at ur in TOTAL LW (W/m**2/Sr)
!
! Note: The IR window is defined as the spectral interval between
!       800 to 1250 cm**-1.
! 6-24-98 (2)
!
! Fu 07-08-98
! The improved parameterization of cirrus radiative properties in
! both solar and IR spectra (Fu 1996; Fu et al. 1998) has been 
! incorporated into the radiation model.  Note that the definition
! of the generalized effective size used in the new parameterization
! (Eq. 3.10 or Eq. 2.3 in Fu 1996) is different from the mean
! effective size defined in Eq. 2.1 of Fu and Liou (1993).  Now 
! you can make choice between the two versions of cirrus para-
! meterization through the logical variable "fl93i".  Use appropriate
! effective sizes of ice clouds through input "pde" for the two
! different versions of parameterization.
! Fu 07-08-98
!
! *********************************************************************
	use TRMM_WINDOW
! **********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
! Fu 07-08-98
        logical fl93i
! Fu 07-08-98
      integer i,ib,ig,mbn,kg(mbx),irobckd,iac
      real pp,pt,ph,po,pre,plwc,pde,piwc,prwc,pgwc,umco2,umch4,umn2o
      real fds,fus,dts,fdir,fuir,dtir,fd,fu,dt,fu1,fd1,bf,bs
      real as(mbsx),as1(10), ee(mbirx),u0,ss,pts,f0,fuq1,fuq2,xx,asx
      real hk

! Chou arrays :5/99
	integer npc,ntop,nbot
	real dn1(nv1x),up1(nv1x),u0d
	data u0d / 0.50 /
! Chou arrays :5/99

!---------- 10/29/96 (1)
!---------- 4/1/97 (1)
      real ctau, dtau, tae, wae, wwae 
	real aprop
      common /aer_proc/ ctau(18), dtau(10)
!---------- 4/1/97 (1)
      common /aer/ tae(nvx,mxac), wae(nvx,mxac), wwae(nvx,4,mxac)
      common /aerpout/aprop(nvx,mxac,3)
!---------- 10/29/96 (1)

!---------- 10/28/96 (1)
      real fdsdr,fdsdf,ffdr,ffdf
      common /radiat/ fds(nv1x), fus(nv1x), dts(nvx),&
                     fdir(nv1x), fuir(nv1x), dtir(nvx),&
                     fd(nv1x), fu(nv1x), dt(nvx),  &
                     fdsdr(nv1x), fdsdf(nv1x)
      common /dirdiff/ffdr(nv1x),ffdf(nv1x)

!---------- 10/28/96 (1)

      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /clouds/ pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
      common /rains/ prwc(nvx)
      common /graups/ pgwc(nvx)
      common /umcon/ umco2, umch4, umn2o 
      common /dfsout/ fu1(nv1x), fd1(nv1x)
      common /planci/ bf(nv1x), bs
!---------- 4/2/97 (1)
        common /robckd/ irobckd
!---------- 4/2/97 (1)
	real otau,tw,ti
      	common /outtau/ otau(mbx) !! 7-21-96
	common /wat/ tw(nvx)
	common /ic/  ti(nvx)


!!! NEW FOR 0598 ----------------------------------------------------------


  	real cfc_conc(3)
	real fuwn,fdwn
	integer idkfr,iwtas
! +++ PAUSE FIX 
        INTEGER ITriErr
!
        common /wndow/ fuwn(nv1x),fdwn(nv1x)
	common /dkfrwn/ idkfr  
	common /cfcs/ cfc_conc
        common /cont_tas/ iwtas
! 6-24-98 (4a)	
	real fiurt,fiurw,fiur,ur
  	common /radiance/ fiurt(nv1x),fiurw(nv1x),fiur(nv1x) 
	real trwn_flt_r, trwn_unf_r, trwn_f
	common /TRMMWNOUT/ trwn_flt_r, trwn_unf_r, trwn_f 
! 6-24-98 (4a)
	integer isolar_spectrum
	real hk1, fk1o3,sol_spect,fk1h2o

        common /band1/ hk1(10), fk1o3(10),sol_spect(0:7),fk1h2o(10)
 	common /select_solar_spectra/ isolar_spectrum

	logical lchou,lband6a,lpar,lray
	common/CHOU/lchou,lband6a,lpar,lray

!   uvfu(nv1,10) upward flux within 10 UV-VIS subintervals (W/m**2)
!   uvfd(nv1,10) downward flux within 10 UV-VIS subintervals (W/m**2)
      real uvfu,uvfd,uvdir,uvdif
      common /uvflux/ uvfu(nv1x,10),uvfd(nv1x,10),   &
                     uvdir(nv1x,10),uvdif(nv1x,10)
	real swfu,swfd,swdir,swdif
	common /swflux/ swfu(nv1x,6),swfd(nv1x,6),   &
                     swdir(nv1x,6),swdif(nv1x,6)
	real rlwfu,rlwfd,sbf,sbs
	common /lwflux/ rlwfu(nv1x,7:18),rlwfd(nv1x,7:18),   &
       sbf(nv1x,7:18),sbs(7:18)

!  ******************************************************************
!  kg(mb) is the number of intervals to perform the g-quadrature in
!  each band to consider the nongray gaseous absorption.  In total,
!  we need to perform 121 spectral calculations in  the  scattering
!  problem for each atmospheric profile.
!  ******************************************************************
      data kg / 10, 8, 12, 7, 12, 5, &
                 2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8 /


!! Depending on K's option # of Solver loops needed changes...
	if ( idkfr == 0) kg(11:13) = (/3,5,2/)
	if ( idkfr == 1) kg(11:13) = (/80,40,10/)
	if ( idkfr == 2) kg(11:13) = (/5,5,2/)

!  *********************************************************************
!  The following variables are declared here, instead of in the calling
!  program, since any change in these would require modifications
!  the code anyway.  A check is inserted here, to make sure the 
!  number of layers given by the user, nv, is less than or equal
!  to the number of levels nvx, given in the header file rad_0598.h.
!  nvx is used for dimensioning arrays.  Uncomment to use.
!  *********************************************************************
      nv1=nv+1        ! Number of levels
      mb=18           ! Number of bands
      mbs=6           ! number of shortwave bands
      mbir=12         ! number of longwave bands
      nc=8            ! number of cloud types
      ndfs=nv         ! number of layers
      mdfs=nv1        ! number of levels
      ndfs4=4*ndfs    ! number of layers * 4
      ndfs2=2*ndfs    ! number of layers * 2



      if (nv.gt.nvx) then
        print *,'ARRAY ERROR: number of levels specified (nv) '
        print *,'is greater than number allowed (nvx)!!!'
        stop
      endif
!c      ITriErr = 0
!  ******************************************************************
!  ******************** 10/23/96 (end) ******************************
!  ******************************************************************

! Fu 07-08-98
 	fl93i = .false. !! False = New 98Ice cld !! True Old93 ice
! Fu 07-08-98

      f0 = 1.0 / 3.14159
      do 10 i = 1, nv1
       fds(i) = 0.0
       fus(i) = 0.0
       fdir(i) = 0.0
       fuir(i) = 0.0
	fuwn(i)=0.0
	fdwn(i)=0.0
!---------- 10/28/96 (2)
       fdsdr(i) = 0.0        
       fdsdf(i) = 0.0  
       if (i.lt.nv1) then
         dts(i)=0.      
         dtir(i)=0.      
         dt(i)=0.  
       endif    
!---------- 10/28/96 (2)
! 6-24-98 (5)	
           fiurw(i) = 0.0
	   fiurt(i) = 0.0
! 6-24-98 (5)	
10     continue

	swfu=0.0
	swfd=0.0
	swdir=0.0
	swdif=0.0
	rlwfu=0.0
	rlwfd=0.0
	sbf=0.0 ;sbs=0
      call thicks
      call rayle2
!
! note 0.03 is the limit in surface albedo paramterization
! less than 0.03, the albedo values could be incorrect.
! see leaft 2 for details


!      if ( u0 .le. 1.0e-4 ) then
   
      if ( u0 .le. 1.0e-4) then
!        print*,'NO LONG WAVE'
        mbn = mbs + 1    ! Longwave only if sun is too low
      else
        mbn = 1          ! Shortwave and longwave
      endif


!---------- 10/29/96 (2)
      call aerosol_init
!---------- 10/29/96 (2)

!TRMM (begin)
	fu_sr=0.0
	fu_sf=0.0
!TRMM (end)

	call ckd1_init( isolar_spectrum,lray )

      do 20 ib = mbn,mb !!!!!!!!!!!!!!!!!!!!!!!!!! mb

           if ( fl93i ) then
	      call ice ( ib )
	   else
              call icenew ( ib )
           endif


! Use Yong Hu's Water cloud optical properties for SW bands
	if ( ib <=6 ) then
	  call water_hu ( ib ) ! CALL Yong Hu's WATER CLOUD OPTICS For SW bands ib=1:6
	else
 	  call water ( ib )
	endif

       call rain ( ib )
       call graup ( ib )

!---------- 4/1/97 (3)
! No more ipr option
       if (ib.ne.1) then
       call aerosolxy (ib,'x')
       ctau(ib)=0.
       do i=1,nv
	do iac=1,nac
        ctau(ib)=ctau(ib)+tae(i,iac)
	end do
        end do
       endif
!---------- 4/1/97 (3)

! --- 7-21-97 -- Vertically Integrated Cloud Water/Ice optical Depth
	otau(ib)=0.0
	do i=1,nv
        otau(ib)=otau(ib)+(tw(i)+ti(i))
	enddo
! --- 7-21-97
	if (ib ==1 .and. lchou ) then
! Chou scheme3 
	ntop=0
	nbot=0
	npc=0
	 do i=1,nv
	  if  ( (tw(i)+ti(i)) .gt. 0.0 )then
	   ntop =i
	   npc= i
	   exit
	  endif
	 enddo
	 do i=nv,1,-1
	  if  ( (tw(i)+ti(i)) .gt. 0.0 )then
	   nbot =i+1
	   npc= i
	   exit
	  endif
	 enddo
! Chou scheme3
	endif

! MOVED i       call rayle ( ib, u0  )
!---------- 4/2/97 (2)
!      call gascon ( ib )  ! REPLACED WITH CONDITIONAL BELOW.
!
      if (irobckd .eq.1) then
        call gascon ( ib )  !! OLD ROBERTS CONTINUUM
      elseif (irobckd .eq.2) then
        call gascon_ckd ( ib ) !! EXACT CKD_2.1 (SLOW)
      elseif (irobckd .eq.3) then
        call gascon_off      !! Option to Turn off Continuum Absorbtion
      elseif (irobckd .eq.4) then
        call gascon_ckd_parm(ib) !! Parameterized CKD_2.1 Cont. Abs.
      else
        stop &
       'Set Continuum Model 1=Roberts 2=CKD_2.1 3=NONE, 4=PARM CKD_2.1'
      endif
			
!---------- 4/2/97 (2)

       if ( ib .gt. mbs ) then
         call planck ( ib, pts )
       endif
				
       do 30 ig = 1, kg(ib)
!---------- 4/1/97 (4)
! No more ipr option
        if (ib.eq.1) then
          call aerosolxy (ig,'y')

	if(ig ==9) then
	 do i=1,nv
	 do iac=1,nac
	 aprop(i,iac,1)= tae(i,iac)
	 aprop(i,iac,2)= wae(i,iac)
	 aprop(i,iac,3)= wwae(i,1,iac) *0.3333
	 enddo;enddo
	endif

          dtau(ig)=0.
          do i=1,nv
	   do iac=1,nac
           dtau(ig)=dtau(ig)+tae(i,iac)
	   end do
           end do
        endif
!---------- 4/1/97 (4)
	call rayle ( ib, u0, ig , lray)
        call gases ( ib, ig, hk )
        call comscp 

!  ---------------------------
!  10/25/96 -- 11/4/95 (begin)
!  ---------------------------
        if ( ib .le. mbs ) then

	asx = as(ib)

!! If as1 array is filled with good albedo for band 1 use it !
	if ( ib == 1        .and.    &
           as1(ig) .ge. 0.0 .and.   &
           as1(ig) .le. 1.0 ) asx = as1(ig)

          if ( fourssl ) then
		
            call qfts ( ib, asx, u0, f0 )  !! FOUR STREAM SOLIR
		
          else
            quadra = .false.
            hemisp = .false.
            edding = .true.
! +++ PAUSE FIX           call qftsts ( ib, asx, u0, f0 ) ! TWO STREAM SOLIR
            call qftsts ( ib, asx, u0, f0, ITriErr ) ! TWO STREAM SOLIR
	            IF (ITriErr .NE. 0) RETURN
          endif
	
          do 40 i = 1, nv1
           fds(i) = fds(i) + fd1(i) * hk
           fus(i) = fus(i) + fu1(i) * hk
!---------- 10/28/96 (3)
           fdsdr(i) = fdsdr(i) + ffdr(i) * hk
           fdsdf(i) = fdsdf(i) + ffdf(i) * hk
!---------- 10/28/96 (3)
	if ( ib == 1 ) then
	 uvfu(i,ig)  = fu1(i)  * hk
	 uvfd(i,ig)  = fd1(i)  * hk
	 uvdir(i,ig) = ffdr(i) * hk
	 uvdif(i,ig) = ffdf(i) * hk
	endif
	if ( ib <= 6 ) then
	 swfu(i,ib)  = swfu(i,ib)  +fu1(i)  * hk
	 swfd(i,ib)  = swfd(i,ib)  +fd1(i)  * hk
	 swdir(i,ib) = swdir(i,ib) +ffdr(i) * hk
	 swdif(i,ib) = swdif(i,ib) +ffdf(i) * hk
	endif
 40        continue

        else

          if ( foursir ) then
            call qfti ( ib, ee(ib-mbs) ) !FOUR STREAM IR
          else
            quadra = .false.
            edding = .false.
            hemisp = .true.
            mquadr = .false.
! 6-24-98 (6)	
!                  call qftisf ( ib, ee(ib-mbs) )!TWO-FOUR STREAM COMB IR
! +++ PAUSE FIX                   call qftisf ( ib, ee(ib-mbs), ur )
                   call qftisf ( ib, ee(ib-mbs), ur, ITriErr )
! 6-24-98 (6)

! +++ PAUSE FIXc               call qftits ( ib, ee(ib-mbs) ) !TWO-STREAM IR
!               call qftits ( ib, ee(ib-mbs), ITriErr ) !TWO-STREAM IR
          endif
!  -------------------------
!  10/25/96 -- 11/4/95 (end)
!  -------------------------

          do i = 1, nv1
           fdir(i) = fdir(i) + fd1(i) * hk
           fuir(i) = fuir(i) + fu1(i) * hk
	   fiurt(i) = fiurt(i) + fiur(i) * hk * 0.1591549
           end do

	if ( ib >= 7 .and. ib <=18) then
   	 do i = 1, nv1
           rlwfu(i,ib) =  rlwfu(i,ib) + fu1(i) * hk
           rlwfd(i,ib) =  rlwfd(i,ib) + fd1(i) * hk
	   sbf(i,ib) = bf(i)
	   sbs(ib) = bs
  	enddo
	endif

	if ( ib.eq.11 .or. ib.eq.12 .or. ib.eq.13) then
		do i=1,nv1
	        fdwn(i) = fdwn(i) + fd1(i) * hk
                fuwn(i) = fuwn(i) + fu1(i) * hk
! 6-24-98 (7)	
           fiurw(i) = fiurw(i) + fiur(i) * hk * 0.1591549
! 6-24-98 (7)
		enddo
!TRMM (begin)
! GET Spectral Flux and Radiance for Filtering to TRMM WINDOW
          fu_sf(ib) = fu_sf(ib) + fu1(1)  * hk
          fu_sr(ib) = fu_sr(ib) + fiur(1) * hk * 0.1591549
!TRMM (end)      
	endif 


        endif

30      continue  
20     continue

!  ------------------------------------------------------------------
!  In this model, we used the solar spectral irradiance determined by
!  Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained 
!  in the spectral region 0.2 - 4.0 um.
!
!  fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
!  hk of 0.03.
!  ------------------------------------------------------------------
      if ( lband6a ) then
!	fuq1 = ss / ( 1340.0 +11.0 )
	fuq1 = ss /  sol_spect(0)  
	else
	fuq1 = ss / ( sol_spect(0) - sol_spect(7) )
       endif


!       fuq2 = bs * 0.03 * 3.14159 * ee(12)
! fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
! hk of 0.03.
        fuq2 = 0.0

      do 60 i = 1, nv1
       fds(i) = fds(i) * fuq1
       fus(i) = fus(i) * fuq1
       fuir(i) = fuir(i) + fuq2
       fd(i) = fds(i) + fdir(i)
       fu(i) = fus(i) + fuir(i)
!---------- 10/28/96 (4)
       fdsdr(i) = fdsdr(i) * fuq1
       fdsdf(i) = fdsdf(i) * fuq1
!---------- 10/28/96 (4)
60     continue

!! Correct Spectral Output
	 uvfu  = uvfu  * fuq1
	 uvfd  = uvfd  * fuq1
	 uvdir = uvdir * fuq1
	 uvdif = uvdif * fuq1
	
	 swfu  = swfu  * fuq1
	 swfd  = swfd  * fuq1
	 swdir = swdir * fuq1
	 swdif = swdif * fuq1

!---------------------------------------------------------------------
! Account for SW absorption by CO2 & O2 after CHOU : Journal Climate Feb90 & CLIRAD-SW NASA TM

!       Change by Nair 16 Oct 2002
!       Originally the logic was 'if (lchou ) then'. However if solar
!       radiation is zero, then the following code causes a lot of
!       problems. So the following section of code is not activated
!       if solar zenith angle > 90.0

	if (lchou .and. u0 .gt. 0.0) then

	do i =1,nv1
	 if ( fdsdr(i)/fds(i)  < 0.50 ) then
	  npc=i 
	if ( ntop==0 .or. nbot ==0 ) then
	  ntop = npc 
	  nbot = npc+1
	endif
	  exit
	 endif
	enddo

	if ( npc .ne.0 .and. npc < ntop) npc = ntop
	if ( npc .ne.0 .and. npc > nbot) npc = nbot

	call add_chou3(nv,umco2,u0,pp,pt,ph,fds,fus,dn1,up1,npc,ntop,nbot)

! Here we assume ALL of the CO2 & O2 absorption comes from the Direct Beam
!	fdsdr(1:nv1)=fdsdr(1:nv1)- (fds(1:nv1)-dn1(1:nv1) )

! Here we assume the CO2 & O2 absorption is split evenly between Direct & diffuse.
	fdsdr(1:nv1)=fdsdr(1:nv1)- (fds(1:nv1)-dn1(1:nv1) )   &
         *fdsdr(1:nv1)/fds(1:nv1)
	fdsdf(1:nv1)=fdsdf(1:nv1)- (fds(1:nv1)-dn1(1:nv1) )   &
         *fdsdf(1:nv1)/fds(1:nv1)

! DEBUG PRINT
!	do i=1,nv+1	
!	print'(A4,2i3,f8.1,3x,2(3f8.2,3x),4(2f8.2,3x))','CPRF',i,npc,pp(i),
!     &      fds(i),dn1(i),dn1(i)-fds(i),
!     &      fus(i),up1(i),up1(i)-fus(i)
!     ,      ,fdsdr(i),fdsdf(i)
!	enddo

!	print'(a30,4f8.4,i8)',
!     & 'CPRF: Toa ALBEDO (w/o, with)',fus(1)/fds(1),up1(1)/dn1(1)
!     &,u0,fus(nv+1)/fds(nv+1),npc

	fds(1:nv1)=dn1(1:nv1) 
	fus(1:nv1)=up1(1:nv1)

!	do i=1,nv1
!	if ( abs (fds(i)-fdsdr(i)-fdsdf(i) ) > 1.0) stop ' CHOU: DIR/DIFFUSE'
!	enddo

!End Chou absorption
	endif !lchou
!---------------------------------------------------------------------

     do 70 i = 1, nv

!      do 70 i = 2, nv

!      xx = fds(i) -fus(i) - fds(i+1) + fus(i+1)
!      dts(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
!      xx = fdir(i) -fuir(i) - fdir(i+1) + fuir(i+1)
!      dtir(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
!      dt(i) = dts(i) + dtir(i)

!      Nair: Changed units of heating rate from K/day to K/s

!       xx = fds(i-1) -fus(i-1) - fds(i) + fus(i)
!       dts(i) = 9.760956175E-5 * xx / ( pp(i) - pp(i-1) )
!       xx = fdir(i-1) -fuir(i-1) - fdir(i) + fuir(i)
!       dtir(i) =  9.760956175E-5 * xx / ( pp(i) - pp(i-1) )
!       dt(i) = dts(i) + dtir(i)
!       print*,'RADPROF',i,fuir(i),fdir(i)
      xx = fds(i) -fus(i) - fds(i+1) + fus(i+1)
      dts(i) = 9.760956175E-5* xx / ( pp(i+1) - pp(i) )
      xx = fdir(i) -fuir(i) - fdir(i+1) + fuir(i+1)
      dtir(i) = 9.760956175E-5* xx / ( pp(i+1) - pp(i) )
      dt(i) = dts(i) + dtir(i)
70     continue
        
!      do i = 1, nv1
!         print*,'KKKK',pp(i),fuir(i),fdir(i)
!      enddo

!TRMM (begin)
	iwncld(11:13) = 0
	where ( otau(11:13) > 1.0 ) iwncld(11:13) = 1

!       call trmm_wnflt(ur)

	 trwn_flt_r = sat_flt_r
	 trwn_unf_r = sat_unf_r
	 trwn_f     = sat_f

!TRMM (end)

      return
      end

