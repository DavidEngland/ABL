Subroutine oned_Raddrvr (p,theta,rhop,q,ozone,cosz, &
alb,em,tr,fthrd,rlong,raduptop,raddntop,delt)
! Code for computing radiative heating tendencies
! Stripped from code Fu-Liou scheme added to RAMS
! Driver interface for RAMS, Soubroutine FuLiou_Raddriv
! Inputs:
! p - Pressure Pa
! t - Temperature K
! rho - Density (Kg/m3)
! rv -  Water vapor mixing ratio kg/kg
! cw -  Cloud liquid water content, set to 0
! ci -  Cloud ice content, set to 0
! rr -  Rain water content, set to 0
! gr -  Graupel content, set to 0
! o3r - Ozone mixing ratio, gm/m3
! nz  - # of pressure levels
! cosz - Cosine of solar zenith angle 
! solconst - Solar Constant
! alb - Broadband albedo
! em  - Broadband emissivity
! sft - Surface Temperature (K) 
! Ouputs:
! rshort - Downwelling shortwave at the surface
! rlong - Downwelling longwave at the surface
! fthrd - Heating rate K/s 

!-------------------SET KZ HERE--------------------------------------------------------
	INTEGER, PARAMETER  :: KZ=100	! NUMBER OF VERTICAL LEVELS 
	
!---------------Must also be set in BL.incl	
	 real,dimension(Kz):: cw,ci,rr,gr,htr,p,theta,rhop,q,ozone,fthrd,T,oldT
        
        real cosz, solconst, alb, em, rshort, rlong,raduptop,raddntop
        integer nz
	nz=kz
	
        do i= 1, nz
           cw(i)=0.0
           ci(i) =0.00938567
           rr(i) = 0.0
           gr(i) =0.0
	   T(i)=Theta(i)*(p(i)/100000.)**.286
	  
        enddo
	
        em = 1.0
        call FuLiou_Raddriv (p,t,rhop,q,cw,ci,rr,gr,&
                             ozone, nz, cosz, solconst,&
                             alb,em,tr,&
                             rshort, rlong, fthrd,raduptop,raddntop)
			     
			     

end

!---------------------------------------------------------------------
! RAMS interface to Fu Liou radiative transfer code. Modified
! by Udaysankar S. Nair. 
!
! Version 0.0 10 OCT 2002
! 
! This is a brute force approach version to get Fu-Liou code to 
! interface with RAMS and get it working. I have not paid 
! much attention to efficiency of computation and memory use.
!
! subroutine FouLiou_Raddriv()
!
! Inputs to the subroutine
! 
! 
!    ===============================================================
!    The following variables are 1D arrays containing the specified 
!    variable value at vertical levels in an atmosperic column 
!    at a particular horizontal grid point location. 
!    ===============================================================
!        real p_rams(nzlev_rams_rad)     : Pressure (Pa)
!        real t_rams(nzlev_rams_rad)     : Temperature (K)
!        real rho_rams(nzlev_rams_rad)   : Density (kg/m3) 
!        real r_rams(nzlev_rams_rad)     : Water vapor (kg/kg)
!        real cldw_rams(nzlev_rams_rad)  : Cloud water (kg/kg)
!        real cldi_rams(nzlev_rams_rad)  : Cloud ice (kg/kg)
!        real rain_rams(nzlev_rams_rad)  : Rain (kg/kg)
!        real graup_rams(nzlev_rams_rad) : Graupel (kg/kg)
!        real o3_rams(nzlev_rams_rad)    : Ozone mixing ratio (g/m3)  
!    ================================================================
!     
!        real sza        : Solar zenith angle at the current RAMS 
!                          grid point
!
!        real nzlev_rams_rad : Number of RAMS vertical levels
!
!
! Outputs from the subroutine
!
!        real fup_rams(nzlev_rams_rad)  : Upwelling flux at each RAMS 
!                                     level
!
!        real fdn_rams(nzlev_rams_rad)  : Downwelling flux at each RAMS
!                                     level
!  
!--------------------------------------------------------------------
	subroutine FuLiou_Raddriv(p_rams,t_rams,rho_rams,&
          r_rams, cldw_rams, cldi_rams, rain_rams, graup_rams,&
          o3_rams, nzlev_rams_rad, cosz, solconst,albedt,em,surft,&
          rshort, rlong, fthrd,raduptop,raddntop)
 
        

	USE RadParams
	use TRMM_WINDOW
	real raduptop,raddntop

	real p_rams(nzlev_rams_rad),t_rams(nzlev_rams_rad), &
             rho_rams(nzlev_rams_rad), r_rams(nzlev_rams_rad), &
             cldw_rams(nzlev_rams_rad),cldi_rams(nzlev_rams_rad), &
             rain_rams(nzlev_rams_rad),graup_rams(nzlev_rams_rad), &
             o3_rams(nzlev_rams_rad),fthrd(nzlev_rams_rad), cosz,  &
             solconst,albedt,em,surft, rshort, rlong
!       integer itp0, igp, jgp
!************************************************************
! common data block for the type of atmosphere, the atmosphere
! could be tropical(ktrop.lay), mid-latitude summer (kmls.lay), 
! mid-lat winter (kmlw.lay), SubarcticSummer(ksas.lay), and 
! SubarcticWinter (ksaw.lay) pp(nv1x) 
! pp : atmospheric pressure (millibars)   
! pt : atmospheric temperature (K)     
! ph : water vapor mixing ratio (kg/kg) 
! po(nv1x)  ozone mixing ratio (kg/kg)   
!*************************************************************
	common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

!*************************************************************
! common data block for the tupe of cloud, the cloub could be
! ICE or WATER. For both of them, the effective radius of water 
! cloud (in um, pre),  liquid water content (g/m**3, plwc)
! and pde(nvx)  effective diameter of ice cloud (um) and  
! piwc(nvx) ice water content (g/m**3) is got from the data in 
! misc_subs.f
!*************************************************************
	common /clouds/ pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)

!*************************************************************
!**    rain water content (g/m**3), from common data block
!*************************************************************
	common /rains/ prwc(nvx)

!*************************************************************
!*pgwc(nvx) graupel water content (g/m**3) 
!*************************************************************
	common /graups/ pgwc(nvx)

!*************************************************************
! umco2     concentration of CO2 (ppmv)
! umch4     concentration of Ch4 (ppmv)
! umn2o    concentration of N2o (ppmv)
!*************************************************************
        common /umcon/ umco2, umch4, umn2o 

!*************************************************************
! fds(nv1x)   downward solar flux ( W / m ** 2 )
! fdsdr(nv1x) downward direct SW flux ( W / m **2 )
! fdsdf(nv1x) downward diffuse SW flux ( W / m **2 )
!*************************************************************
	common /radiat/ fds(nv1x), fus(nv1x), dts(nvx), &
                        fdir(nv1x), fuir(nv1x), dtir(nvx), &
     	                fd(nv1x), fu(nv1x), dt(nvx), &
                        fdsdr(nv1x), fdsdf(nv1x)  

	common /dfsout/ fu1(nv1x), fd1(nv1x)

!  ******************************************************************
!  bf and bs are the blackbody intensity function integrated over the
!  band ib at the nv1 levels and at the surface, respectively.    The
!  units of bf and bs are W/m**2/Sr. nd*10 is the band width from ve.
!  ******************************************************************
	common /planci/ bf(nv1x), bs

!*************************************************************
!  irobckd   LW continuum option
!*************************************************************
	common /robckd/ irobckd,no_rayle

!  ******************************************************************
!! WINDOW FLUX
!  ******************************************************************
        common /wndow/ fuwn(nv1x),fdwn(nv1x)

!********************************************************************
!  idkfr     Window K's Option
! idkfr=0 ! 1= DAVE KRATZ 0= Original Fu 2=HYBRED
!********************************************************************
	common /dkfrwn/ idkfr

!********************************************************************
!fc_conc(3)  Window K's now include CFC's
!*******************************************************************
        common /cfcs/ cfc_conc(3)

!***********************************************************************
!   1:Linea      2:log ,
!             3:linear w/Plank wgt   :: SUGGEST iwtas=3
!             4: Log w/plank wgt
!***********************************************************************
	common /cont_tas/ iwtas

!***********************************************************************
! fiurt(nv1x)  - TOTAL LW (0-2200cm-1) radiance [Wm-2sr-1]
!             fiurt(nv1) upward radiance at ur in TOTAL LW (W/m**2/Sr)
! fiurw(nv1x)  - WINDOW (800-1250cm-1) radiance [Wm-2sr-1]
!             fiurw(nv1) upward radiance at ur in IR window (W/m**2/Sr)
! fiurw(nv1) upward radiance at ur in IR window (W/m**2/Sr)
!              fiurt(nv1) upward radiance at ur in TOTAL LW (W/m**2/Sr)
! 6-24-98 (4a)	
!***********************************************************************
  	common /radiance/ fiurt(nv1x),fiurw(nv1x),fiur(nv1x) 

	common /TRMMWNOUT/ trwn_flt_r, trwn_unf_r, trwn_f 
! 6-24-98 (4a)
        common /select_solar_spectra/ isolar_spectrum
	character*10 cchou,cband6a,cpar,cray
	logical lchou,lband6a,lpar,lray
	common/CHOU/lchou,lband6a,lpar,lray

!***************************************************************
! uvfu(nv1,10) upward flux within 10 UV-VIS subintervals (W/m**2)
!***************************************************************
        common /uvflux/ uvfu(nv1x,10),uvfd(nv1x,10),   &
                     uvdir(nv1x,10),uvdif(nv1x,10)
	
	common /swflux/ swfu(nv1x,6),swfd(nv1x,6),   &
                     swdir(nv1x,6),swdif(nv1x,6)


	common /lwflux/ rlwfu(nv1x,7:18),rlwfd(nv1x,7:18),   &
       sbf(nv1x,7:18),sbs(7:18)

	common /aerpout/aprop(nvx,mxac,3)
        real szen
        common /solang/szen

	real px(nv1x),dtx(nvx),rh(nv1x)
	real sh_aer(3)

	character*20 state(17)
	character*20 aerosol_tau,ray_type
	character*30 cont_type

	real ,dimension(0:3) :: rhx,pwx,skinx
	real ,dimension(12) :: spece
	real outp(7,0:4)
	real cool_rate(nvx)
	data rhx  /1.0,1.0,1.0,1.0/
	data pwx  /1.0,2.0,1.0,2.0/
	data skinx/0.0, 0.0,0.0, 0.0/
        data spece/0.835, 0.916, 0.934, 0.923, 0.835, 0.877,&
               0.921, 0.926, 0.934, 0.934, 0.934, 0.934/  
!======================================================
! for surface, as: surface albedo at the solar 6 bands
!               as1: surface albedo at the10 intervals 
!		   at the first solar band
! ee : infrared
!=======================================================
	real*4 as(mbsx),as1(10), ee(mbirx), salb(mbsx)
	character*40 atm
	character*1 stream,pha,cwtas,cspec,caform,crh
	character*7 water_content,particle_radius
	character*7 cos_zen,sfc_alb,optical_depth
	character*7 co2_con,cdkfr,lwtas(0:4),cur
	character*30 aerosol_type(1:21) 
	character*2 ldkfr(0:2)
	character*2 cmlev
	character*4 ftag
	gpt(vtbar1,p0,p1) = 29.3* vtbar1* log(p0/p1)

	data lwtas/'Orig','LinTau ','LogTau','LinPLNK ','LogPlnk'/
	data ldkfr/'FU','DK','HY'/
! 18 aerosol types
	data aerosol_type/'Marine', 'continental', 'urban', '0.5_dust', &
        '1.0_dust', '2.0_dust', '4.0_dust', '8.0_dust', 'insoluble', &
        'soluble', 'soot', '(nucl)sea salt', '(coar)sea salt', &
        'nucleation dust', &
        'accumulation dust', 'coarse dust', 'transported dust',  &
        'sulfate droplets' , 'NAAPS DUST', 'NAAPS smoke','PRIDE DUST'/ 

101	format(1x, 6(1x, f7.5))
102	format(I2)

        ur = 0.5

        FOURSSL  = 0
        FOURSIR  = 0
        LCHOU    = 0
        LBAND6A  = 1
        LPAR     = 1
        LRAY     = 1
        ISOLSPEC = 1
        IROBCKD  = 1
        IDKFR    = 2
        IWTAS    = 3
        NAC      = 0
        IATYPE   = 1
        NATAU    = 1
        AWVL     = 0.55
        ATAU     = 0.0
        IVD      = 0
        IAFORM   = 1
        IARH     = 0
        

! print*,' IDKFR =',idkfr
! print*,' IWTAS =',iwtas
! print*,' LCHOU=',lchou
! print*,' LBAND6A',lband6a
! print*,' IFG', IFG
! print*,' IAFORM', IAFORM
! print*,' IROBCKD',IROBCKD
! print*,'IDKFR', IDKFR
! print*,' IVD',IVD
! print*,'ISOLAR_SPECTRUM',ISOLAR_SPECTRUM
! print*,' LPAR',lpar
! print*,' LRAY',lray 
! print*,'nzlev_rams_rad',nzlev_rams_rad

!========================================================================
!
! Transfer the atmospheric profile values from RAMS column to Fu-Liou
! grid.
!========================================================================
	
        levfl = nzlev_rams_rad -1
        do ilev = 2, nzlev_rams_rad 

!        levfl = nzlev_rams_rad
!        do ilev = 1, nzlev_rams_rad 
         
           pp(levfl) = p_rams(ilev)/100.0
           pt(levfl) = t_rams(ilev)
!           ph(levfl) = r_rams(ilev)
!          changed by nair 11 Jan 2003
           ph(levfl)  = r_rams(ilev)/rho_rams(ilev)
           po(levfl) = o3_rams(ilev)/(1000.0*rho_rams(ilev))
           pre(levfl) = 0.0
           plwc(levfl) = 0.0
           pde(levfl) = 0.0
           piwc(levfl) = 0.0
           prwc(levfl) = 0.0
           pgwc(levfl) = 0.0
!          print*,'p=',pp(levfl),'t=',pt(levfl),'ph=',ph(levfl),'po=',o3_rams(levfl)
           levfl = levfl -1
        enddo


        nv = nzlev_rams_rad-2
!        nv = nzlev_rams_rad-1
         
	nv1=nv+1

!! set up rad24a
	ndfs = nv
	mdfs =nv1
	ndfs4 = 4*ndfs 
	ndfs2=  2*ndfs
	mb   = 18
	mbs  = 6
	mbir = 12 
	nc   = 8 

	ss= solconst 

        umco2 = 370.0
        umch4 = 1.75
        umn2o = 0.31
        cfc_conc= (/0.268e-09 , 0.503e-09 ,0.105e-09/)

	do iiii = 1, mbsx
	as(iiii)=albedt
	enddo
	
	as1(1:10) =albedt


!       as1(1:iii) = albedt


!       ee(1:mbirx)=1.0
            
	do iiii = 1, mbirx
	ee(iiii)=em
	enddo
       

!        if(fourssl)then
!           print*,'FOUR STREAM CHOSEN FOR SW'
!        else
!           print*,'TWO  STREAM CHOSEN FOR SW'
!        endif
!           
!        if(foursir)then
!           print*,'FOUR STREAM CHOSEN FOR LW'
!        else
!           print*,'TWO  STREAM CHOSEN FOR LW'
!        endif

!        print*,'Number of aerosol types = ',nac
                  
!        print*,'UMCO2 = ',umco2
!        print*,'UMCH4 = ',umch4
!        print*,'UMN2O = ',umn2o
             
!        print*,'CFC(1)',cfc_conc(1)
!        print*,'CFC(2)',cfc_conc(2)
!        print*,'CFC(3)',cfc_conc(3)
                 
!        do i = 1, nac
!           print*,'Constituent = ', i,'Aerosol type =',aerosol_type(itps(i))
!           print*,'Number of wavelengths = ',N_ATAU
!           do j = 1, n_atau
!              print*,'WL =',a_wlis(j,i),'Tau =',a_taus(j,i)
!           enddo
!        enddo

        u0 = cosz
        szen = cosz
        pts = surft
           
	call rad ( as,as1, u0, ss, pts, ee , ur ,ITriErr)
             
        rshort = fds(nv1)
           
        rlong  = fdir(nv1)
	
	Raduptop = fuir(2)
	Raddntop = fdir(2)
       
        levfl = nv
!        print*,'DEBUG VALS',igp,jgp,u0
        do ilev = 2, nv1
           fthrd(ilev) = dt(levfl)
           levfl = levfl -1
!           print*,'XXXX',ilev,&
!           p_rams(ilev),t_rams(ilev),rho_rams(ilev),&
!           r_rams(ilev),o3_rams(ilev), fthrd(ilev),&
!           fds(ilev),fdir(ilev),fus(ilev),fuir(ilev)
        enddo
!       Make consistent with Wang's code. 10/18/07
        fthrd(1) = fthrd(2)
        
	end
