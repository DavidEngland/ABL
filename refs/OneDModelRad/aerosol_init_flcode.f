      subroutine aerosol_init
!
!                        8/14/95, 4/1/97 , 2/10/2000
!
!  **********************************************************************
!  Subroutine to create aerosol optical properties.  There are several
!  inputs and 6 outputs.  
!
!    INPUTS FROM COMMON BLOCKS OR HEADER FILE:
!
!    a_tau(nwi) :  The input column aerosol optical depth
!    (real)           (common block "aer_tau" - see header file).
!
!    a_wli(nwi) :  Wavelength in microns corresponding to aerosol tau in "a_tau"
!
!    aprof(# layers): The input aerosol optical depth profile - LAYERS
!    (real)           (common block "aer_prof").
!
!    itp:       Aerosol type, given in header file rad_0598.h.
!
!    ifg:       The table will compute vertical distributions based on
!    (integer)  relative humidity (see explanation below).  If ifg is
!               set to 0, each layer will have properties calculated
!               based on the relative humidity of that layer.  If ifg
!               is set equal to another integer (1 through the number of
!               relative humidities given in the block data "aerosol")
!               the routine will calculate a vertical profile of optical
!               properties based on the relative humidity corresponding
!               to the index given.  The indices are: 1: 0%; 2: 50%;
!               3: 70%; 4: 80%; 5:90%; 6: 95%; 7: 98%; and 8: 99%.
!               If the number of relative humidities changes, these
!               numbers will have to be modified.
!
!    ivd:       Vertical tau distribution flag.  If set to zero, the 
!               distribution is based on Jim Spinhirne's marine 
!               distribution formulation, and no user input is required.  
!               If set to one, the user's own vertical distribution is 
!               used, and must be present in the array aprof(nlayers).
!               NOTE: This vertical distribution is used as a weighting 
!               factor ONLY, to distribute input column optical depths!
!
!----------------------------------------------------------------------------
!    a_ssa, a_ext, a_asy:  Input single-scattering albedos, extinction
!           coefficients, and asymmetry parameters.  These variables 
!           are dimensioned (# of bands, # of relative humidities,
!           # of aerosol types). An x or y is appended on these 
!           variable names: if x, the numbers correspond to the 18 
!           original bands.  If y, the numbers are for the 10 
!           sub-intervals in the first shortwave band (.2-.7 microns).  
!           All of these variables come from the block data statements 
!           aerosol# (# corresponds to an integer, eg. aerosol1) and 
!           are in common blocks aer_optx and aer_opty.
!
!    nv,mb,pp,pt,ph,dz: number of layers, number of bands, and the
!           pressure, temperature, humidity and thickness profiles.
!           These are shared by several subroutines.
!
!    OUTPUTS:
!
!    a_tau1,a_ssa1,a_asy1:  The optical depth, single-scattering albedo,
!       and asymmetry parameter vertical profiles for 18 bands.  These
!       are dimensioned (nvx, 18)  These are in the common block
!       aer_initx, which is shared by the subroutine "aerosolx".  
!
!    a_tau2,a_ssa2,a_asy2:  Properties for SW band 1's 10 subintervals.  
!       These are dimensioned (nvx, 10)  These are in the common block
!       aer_inity, which is shared by the subroutine "aerosoly".  
!
!  **********************************************************************
      USE RadParams
!##      include 'rad_0698.h'
      implicit none
      integer iq,mtop,n,m,ict,ix,iy,irh,krh,iac,itp
      integer ib1
      real, dimension(mbx,nrh,naer) :: a_ssax,a_extx,a_asyx
      real, dimension(mby,nrh,naer) :: a_ssay,a_exty,a_asyy

      real, dimension(nvx) :: tauxxx
      real, dimension(nvx,mbx,mxac) :: a_tau1,a_ext1,a_ssa1,a_asy1
      real, dimension(nvx,mby,mxac) :: a_tau2,a_ext2,a_ssa2,a_asy2

      real ,dimension(nvx)  :: taux1,taux2,rh,ht1,rhp
      real sumxxx

      real,dimension(mxat) :: a_wli,a_tau
      real,dimension(nvx)  :: aprof
	real,dimension(nvx,mbx) :: wvd_x
	real,dimension(nvx,mby) :: wvd_y

      real pp,pt,ph,po,dz,p1,h1,z,sig,tp
      real rhx(nrh)
      real sclht(12), layht(12)
      real wts(4),tau3(2),tau3y(4)
      real aotf,wlf,sump,rirh
      real spinhirne_sig, spinhirne_tau
      real aero_dist,vht,colval,tsum,zt,zb,xfact,usrsclht,szen
      real coltau1(mbx),coltau2(mby),coltau3(mbx),coltau4(mby)
	
      common /aer_optx/ a_ssax,a_extx,a_asyx
      common /aer_opty/ a_ssay,a_exty,a_asyy
      common /aer_initx/ a_tau1,a_ssa1,a_asy1
      common /aer_inity/ a_tau2,a_ssa2,a_asy2
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /thick/ dz(nvx)

      common /tau_spline_aot/ aotf(15),wlf(15)
      common /nair_debug/aero_dist(nvx),vht(nvx),colval,tsum
      common /aero_vert/usrsclht
      common /solang/szen

      data rhx /0.,50.,70.,80.,90.,95.,98.,99./
      data wts /.23015,.28274,.25172,.23539/
!      data sclht /4.0,4.0,4.0,4.0,2.0,1.0,1.0,1.0,1.0,0.7,8.0,8.0/
      data sclht /8.0,8.0,8.0,8.0,2.0,1.0,1.0,1.0,1.0,0.7,8.0,8.0/
!      data sclht /8.0,8.0,8.0,8.0,2.0,8.0,8.0,8.0/
      data layht /2.0,2.0,2.0,2.0,6.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0/
      
      integer i 
      integer frst /0/
      integer initprof /0/

!      PRINT*,'AEROSOL INIT NEW'

!  Initialize.

	 if(frst .eq. 0)rh     = -9999.
	 a_ssa1 = 0. ; a_ext1 = 0. ; a_asy1 = 0. ; a_tau1 = 0.
	 a_ssa2 = 0. ; a_ext2 = 0. ; a_asy2 = 0. ; a_tau2 = 0.

	if ( nac < 0 .or. nac > mxac ) stop 'nac:# Aerosol Constituents'
	if (n_atau<0 .or.n_atau>mxat) stop 'n_atau:# Aerosol Tau / Wavelengths'
	if (ifg < 0 .or. ifg > 8) stop 'ifg: Aerosol RH% Flag'
	AEROSOL_CONSTITUENTS : do iac = 1,nac

	a_wli(1:n_atau) = a_wlis(1:n_atau,iac)
	a_tau(1:n_atau) = a_taus(1:n_atau,iac)
	aprof(1:nvx)  = aprofs(1:nvx,iac)
	itp	      = itps(iac)
!        write(*,*) 'itp = ', itp, 'naer = ', naer 
	if ( itp .lt. 1 .or. itp .gt. naer ) stop 'itp:Bad Aerosol Type'
!       print*,'CONSTITUENTS',iac,itp

! FOR Aerosol Optical Properties types that are constant with RH	 
	if (itp==1  .or. itp==2 .or. itp==3 .or. itp==4 .or. &
            itp==5 .or. itp==6 .or.itp==7 .or. itp==8 .or. itp & 
            == 10) then
!!       Has already been filled in Block data 
	else
	do krh=2,8
	 a_extx(1:mbx,krh,itp)= a_extx(1:mbx,1,itp)
	 a_ssax(1:mbx,krh,itp)= a_ssax(1:mbx,1,itp)
	 a_asyx(1:mbx,krh,itp)= a_asyx(1:mbx,1,itp)

	 a_exty(1:mby,krh,itp)= a_exty(1:mby,1,itp)
	 a_ssay(1:mby,krh,itp)= a_ssay(1:mby,1,itp)
	 a_asyy(1:mby,krh,itp)= a_asyy(1:mby,1,itp)

	enddo

	endif
!	if ( ifg .ne.0) print*,'CHECK',ifg,itp,a_ssax(1:mbx,ifg,itp)

!  ******************************************************************
!  Calculate heights at center of layer - find highest layer to place
!  aerosols (15 km) - calculate relative humidities of each layer as
!  needed.  Values of RH > 99% will be set equal to 99% to make table
!  lookup easier. "mtop" is the highest aerosol layer.
!  ******************************************************************
      z=0.
      m=nv
      iq=0
      do while (iq.eq.0 .and. m .ne. 0)
       ht1(m)=(z*2.+dz(m))/2.
       z=z+dz(m)
       if (z.gt.15.) then
         iq=1
         mtop=m
       endif
       p1=(pp(m)+pp(m+1))/2.
       tp=(pt(m)+pt(m+1))/2.
       h1=(ph(m)+ph(m+1))/2.
       call ql_rh(rh(m),tp,p1,h1)
       !print*,'RH= ',m,rh(m),ht1(m),dz(m),z,ph(m),nvx
       if (rh(m).gt.98.9) rh(m)=98.9
       if ((rh(m).lt..01).and.(rh(m).gt.-999.)) rh(m)=0.
       m=m-1
       end do

!  *************************************************************
!  Calculate vertical distribution of asymmetry, ss albedo and
!  extinction, based on aerosol type and relative humidity.  
!  If ifg is not equal to 0, parameters  will corresponds to a 
!  single RH, as described in header file. Loop 31 deals with 
!  the 18 original bands, loop 32 with the 10 band 1 subintervals.
!  *************************************************************
      do 30 m=mtop,nv
       do 31 n=1,mbx
        if (rh(m).eq.-9999.) then
          a_ext1(m,n,iac)=-9999.
          a_ssa1(m,n,iac)=-9999.
          a_asy1(m,n,iac)=-9999.
        else
          if (ifg.eq.0) then          ! Dependence on layer RH.
            ict=2
            do while (rh(m).ge.rhx(ict))
             ict=ict+1
             end do
            a_ext1(m,n,iac)=a_extx(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_extx(n,ict,itp)-a_extx(n,ict-1,itp))
            a_ssa1(m,n,iac)=a_ssax(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_ssax(n,ict,itp)-a_ssax(n,ict-1,itp))
            a_asy1(m,n,iac)=a_asyx(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_asyx(n,ict,itp)-a_asyx(n,ict-1,itp))
	  rhp(m) = rh(m)
!            if(m .eq. nv)then
!               print*,'EXT =',a_ext1(m,n,iac), a_ssa1(m,n,iac)
!               print*,'AS EXT =',a_ext1(m,n,iac)*(1.0-a_ssa1(m,n,iac))
!            endif
          else                        ! Dependence on prescribed RH.
            a_ext1(m,n,iac)=a_extx(n,ifg,itp)
            a_ssa1(m,n,iac)=a_ssax(n,ifg,itp)
            a_asy1(m,n,iac)=a_asyx(n,ifg,itp)
          endif
        endif
 31     continue
!-------------------------------------------
       do 32 n=1,mby
        if (rh(m).eq.-9999.) then
          a_ext2(m,n,iac)=-9999.
          a_ssa2(m,n,iac)=-9999.
          a_asy2(m,n,iac)=-9999.
        else
          if (ifg.eq.0) then          ! Dependence on layer RH.
            ict=2
            do while (rh(m).ge.rhx(ict))
             ict=ict+1
             end do
            a_ext2(m,n,iac)=a_exty(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_exty(n,ict,itp)-a_exty(n,ict-1,itp))
            a_ssa2(m,n,iac)=a_ssay(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_ssay(n,ict,itp)-a_ssay(n,ict-1,itp))
            a_asy2(m,n,iac)=a_asyy(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_asyy(n,ict,itp)-a_asyy(n,ict-1,itp))
          else                        ! Dependence on prescribed RH.
            a_ext2(m,n,iac)=a_exty(n,ifg,itp)
            a_ssa2(m,n,iac)=a_ssay(n,ifg,itp)
            a_asy2(m,n,iac)=a_asyy(n,ifg,itp)
          endif
        endif
 32     continue

 30    continue


        coltau1(1:mbx) =0.0
        coltau2(1:mby) =0.0
        coltau3(1:mbx) =0.0
        coltau4(1:mby) =0.0

        VERTICAL : do  m=mtop,nv
            vht(m) = ht1(m)
            zb = ht1(m)-(dz(m)*0.5)
            zt = ht1(m)+(dz(m)*0.5) 
                 
            if(zb .lt. layht(itp) .and. zt .gt. layht(itp))then
               zt = layht(itp)
            endif
              
            if(zt .le. layht(itp))then
               xfact = sclht(itp)*(exp(-(zb/sclht(itp)))-&
                                   exp(-(zt/sclht(itp))))
            else
               xfact = 0.0
            endif

            if(frst .eq. 0) then
               print*,'z1 z2 x bext tau',zb,zt,xfact,a_ext1(m,1,iac),&
               xfact*a_ext1(m,1,iac)
               if(usrsclht .le. 0.0)frst=1
            endif

            a_tau1(m,1:18,iac) = xfact*a_ext1(m,1:18,iac)
            a_tau2(m,1:10,iac) = xfact*a_ext2(m,1:10,iac)
             
            coltau1(1:mbx) = coltau1(1:mbx)+ a_tau1(m,1:mbx,iac)
            coltau2(1:mby) = coltau2(1:mby)+ a_tau2(m,1:mby,iac)

            aero_dist(m) = a_tau1(m,1,iac)

        enddo VERTICAL

        if(usrsclht .ge. -0.001) then
           if(frst .eq. 0) then
              print*,'Using scl ht of',usrsclht
              frst = 1
           endif
            VERTICAL2 : do  m=mtop,nv

                zb = ht1(m)-(dz(m)*0.5)
                zt = ht1(m)+(dz(m)*0.5) 
                 
                if(zb .lt. layht(itp) .and. zt .gt. layht(itp))then
                   zt = layht(itp)
                endif
              
                if(zt .le. layht(itp))then
                   xfact = usrsclht*(exp(-(zb/usrsclht))-&
                                     exp(-(zt/usrsclht)))
                else
                   xfact = 0.0
                endif

                a_tau1(m,1:18,iac) = xfact*a_ext1(m,1:18,iac)
                a_tau2(m,1:10,iac) = xfact*a_ext2(m,1:10,iac)

                coltau3(1:mbx) = coltau3(1:mbx)+ a_tau1(m,1:mbx,iac)
                coltau4(1:mbx) = coltau4(1:mbx)+ a_tau2(m,1:mbx,iac)
            
            enddo VERTICAL2

            VERTICAL3 : do  m=mtop,nv

                zb = ht1(m)-(dz(m)*0.5)
                zt = ht1(m)+(dz(m)*0.5) 
                 
                if(zb .lt. layht(itp) .and. zt .gt. layht(itp))then
                   zt = layht(itp)
                endif
              
                if(zt .le. layht(itp))then
                   xfact = usrsclht*(exp(-(zb/usrsclht))-&
                                     exp(-(zt/usrsclht)))
                else
                   xfact = 0.0
                endif
                if(szen .lt. 0.0 .and. usrsclht .gt. 99.0)xfact=0.0
                if(szen .gt. 0.0)xfact = xfact*0.5
                a_tau1(m,1:mbx,iac) = (coltau1(1:mbx)/coltau3(1:mbx))&
                                     *xfact*a_ext1(m,1:mbx,iac)
                a_tau2(m,1:mby,iac) = (coltau2(1:mby)/coltau4(1:mby))&
                                     *xfact*a_ext2(m,1:mby,iac)
                aero_dist(m) = a_tau1(m,1,iac)
                  
            enddo VERTICAL3
        endif 

	enddo AEROSOL_CONSTITUENTS
      
      return
      end

!===========================================================================
      subroutine aerosolxy ( ib,cmode )
! *********************************************************************
!                      Modified 2/14/00
!
! tae, wae, and wwae are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
! due to the Mie scattering of aerosols for a given layer. 
!
!  This subroutine is called for bands 2 - 18 (ib) 
!  or vis subbands 1-10 (ig)
! *********************************************************************
      USE RadParams
      implicit none
	character*1 cmode
      integer i,ib,iac
      real x1,x2,x3,x4,y1,y2,y3,y4,tae,wae,wwae
      real ,dimension(nvx,18,mxac) :: a_tau1,a_ssa1,a_asy1
      real ,dimension(nvx,10,mxac) :: a_tau2,a_ssa2,a_asy2
      common /aer_initx/ a_tau1,a_ssa1,a_asy1
      common /aer_inity/ a_tau2,a_ssa2,a_asy2

      common /aer/ tae(nvx,mxac), wae(nvx,mxac), wwae(nvx,4,mxac)
     
      AEROSOL_CONSTITUENTS  : do iac=1,nac

      LEVELS : do  i = 1, nv
       select case (cmode)
	case ('x')
         tae(i,iac) = a_tau1(i,ib,iac)
         wae(i,iac) = a_ssa1(i,ib,iac)
         x1     = a_asy1(i,ib,iac)
	case ('y')
         tae(i,iac) = a_tau2(i,ib,iac)
         wae(i,iac) = a_ssa2(i,ib,iac)
         x1     = a_asy2(i,ib,iac)
       end select

       x2 = x1 * x1
       x3 = x2 * x1
       x4 = x3 * x1
       y1 = 3.0 * x1
       y2 = 5.0 * x2
       y3 = 7.0 * x3
       y4 = 9.0 * x4
  
       wwae(i,1,iac) = y1
       wwae(i,2,iac) = y2
       wwae(i,3,iac) = y3
       wwae(i,4,iac) = y4

	enddo LEVELS
	enddo AEROSOL_CONSTITUENTS

      return
      end
!----------------------------------------------------------------
	real function spinhirne_sig(ht1)

       data sig0,a,ap,b,bp,f/0.025,0.4,2981.0,1.6,2.5,1.5e-7/

         t1=  sig0*(1+a)**2
         t4 = f*(1+ap)**2

         t2 = exp(ht1/b)
         t3 = (a+exp(ht1/b))**2
         t5 = exp(ht1/bp)
         t6 = (a+exp(ht1/bp))**2
         spinhirne_sig=(t1*t2/t3)+(t4*t5/t6)   ! scattering coefficient

	return 
	end
!---------------------------------------------
	real function spinhirne_tau(sig,ssa,dz)
	ext = sig / ssa
	spinhirne_tau = ext / dz
	return
	end
