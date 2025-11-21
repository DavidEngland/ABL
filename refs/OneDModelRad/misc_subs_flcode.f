      subroutine thicks 
! **********************************************
! dz is the thickness of a layer in units of km.
! **********************************************
      USE RadParams
      implicit none
      integer i
      real pp,pt,ph,po,dz
!##      include 'rad_0698.h'
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /thick/ dz(nvx)

      do 100 i = 1, nv
       dz(i) = 0.0146337 * ( pt(i) + pt(i+1) )    &
               * alog( pp(i+1) / pp(i) )
100    continue

      return
      end

      subroutine ice ( ib )
! *********************************************************************
! ti, wi, and wwi are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and
! 4) due to the scattering of ice clouds for a given layer.
! *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib,ibr
      real pre,plwc,pde,piwc,ap,bp,cps,dps,cpir,dz,ti
      real wi,wwi,fw1,fw2,fw3,fd,wf1,wf2,wf3,wf4,gg,x1,x2,x3,x4

      common /clouds/ pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
      common /ic1/ ap(3,mbx), bp(4,mbx), cps(4,4,mbsx), dps(4,mbsx),   &
                     cpir(4,mbirx)
      common /thick/ dz(nvx)
      common /ic/ ti(nvx), wi(nvx), wwi(nvx,4)

! ******************************************************************
! The constant 1000.0 below is to consider the units of dz(i) is km.
! ******************************************************************
      do 10 i = 1, nv
       if ( piwc(i) .lt. 1.0e-5 ) then
         ti(i) = 0.0
         wi(i) = 0.0
         wwi(i,1) = 0.0
         wwi(i,2) = 0.0
         wwi(i,3) = 0.0
         wwi(i,4) = 0.0
       else
         fw1 = pde(i)
         fw2 = fw1 * pde(i)
         fw3 = fw2 * pde(i)
         ti(i) = dz(i) * 1000.0 * piwc(i) * ( ap(1,ib) +   &
                 ap(2,ib) / fw1 + ap(3,ib) / fw2 )
         wi(i) = 1.0 - ( bp(1,ib) + bp(2,ib) * fw1 +   &
                 bp(3,ib) * fw2 + bp(4,ib) * fw3 )

         if ( ib .le. mbs ) then

           fd =  dps(1,ib) + dps(2,ib) * fw1 +   &
                 dps(3,ib) * fw2 + dps(4,ib) * fw3
           wf1 = cps(1,1,ib) + cps(2,1,ib) * fw1 +   &
                 cps(3,1,ib) * fw2 + cps(4,1,ib) * fw3
           wwi(i,1) = ( 1.0 - fd ) * wf1 + 3.0 * fd
           wf2 = cps(1,2,ib) + cps(2,2,ib) * fw1 +   &
                 cps(3,2,ib) * fw2 + cps(4,2,ib) * fw3
           wwi(i,2) = ( 1.0 - fd ) * wf2 + 5.0 * fd
           wf3 = cps(1,3,ib) + cps(2,3,ib) * fw1 +   &
                 cps(3,3,ib) * fw2 + cps(4,3,ib) * fw3
           wwi(i,3) = ( 1.0 - fd ) * wf3 + 7.0 * fd
           wf4 = cps(1,4,ib) + cps(2,4,ib) * fw1 +   &
                 cps(3,4,ib) * fw2 + cps(4,4,ib) * fw3
           wwi(i,4) = ( 1.0 - fd ) * wf4 + 9.0 * fd

         else

           ibr = ib - mbs
           gg = cpir(1,ibr) + cpir(2,ibr) * fw1 +   &
                cpir(3,ibr) * fw2 + cpir(4,ibr) * fw3
           x1 = gg
           x2 = x1 * gg
           x3 = x2 * gg
           x4 = x3 * gg
           wwi(i,1) = 3.0 * x1
           wwi(i,2) = 5.0 * x2
           wwi(i,3) = 7.0 * x3
           wwi(i,4) = 9.0 * x4

         endif
       endif
10     continue
      
      return
      end

      subroutine water ( ib )
! ******************************************************************
! tw, ww, and www are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and
! 4) due to the Mie scattering of water clouds for a given layer. 
! By using the mean single scattering properties of the eight drop
! size distributions in each spectral band, the single scattering
! properties of a water cloud with the given liquid water content
! and effective radius are obtained by interpolating (Eqs. 4.25 -
! 4.27 of Fu, 1991). 
! ******************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib,j
      real pre,plwc,pde,piwc,re,fl,bz,wz,gz,dz,tw,ww,www
      real gg,x1,x2,x3,x4
      common /clouds/ pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
      common /wat1/ re(ncx),fl(ncx),bz(ncx,mbx),wz(ncx,mbx),gz(ncx,mbx)
      common /thick/ dz(nvx)
      common /wat/ tw(nvx), ww(nvx), www(nvx,4)

      do 10 i = 1, nv
         if ( plwc(i) .lt. 1.0e-5 ) then
             tw(i) = 0.0
             ww(i) = 0.0
             www(i,1) = 0.0
             www(i,2) = 0.0
             www(i,3) = 0.0
             www(i,4) = 0.0
         else

           if ( pre(i) .lt. re(1) ) then
!            ---------------------------------------------------- 
!            A cloud with the effective radius smaller than 4.18 
!            um is assumed to have an effective radius of 4.18 um 
!            with respect to the single scattering properties. 
!            ---------------------------------------------------- 
             tw(i) = dz(i) * plwc(i) * bz(1,ib) / fl(1)
              ww(i) = wz(1,ib)
               x1 = gz(1,ib)
               x2 = x1 * gz(1,ib)
               x3 = x2 * gz(1,ib)
             x4 = x3 * gz(1,ib)
             www(i,1) = 3.0 * x1
             www(i,2) = 5.0 * x2
             www(i,3) = 7.0 * x3
             www(i,4) = 9.0 * x4

           elseif ( pre(i) .gt. re(nc) ) then
!            ---------------------------------------------------- 
!            A cloud with the effective radius larger than 31.23 
!            um is assumed to have an effective radius of 31.18 um 
!            with respect to the single scattering properties.  
!            ---------------------------------------------------- 
             tw(i) = dz(i) * plwc(i) * bz(nc,ib) / fl(nc)
             ww(i) = wz(nc,ib)
             x1 = gz(nc,ib)
               x2 = x1 * gz(nc,ib)
             x3 = x2 * gz(nc,ib)
               x4 = x3 * gz(nc,ib)
             www(i,1) = 3.0 * x1
             www(i,2) = 5.0 * x2
             www(i,3) = 7.0 * x3
             www(i,4) = 9.0 * x4

           else

             j = 1
             do while ((pre(i) .lt. re(j)).or.(pre(i) .gt. re(j+1)))
              j = j + 1
              end do

             tw(i) = dz(i) * plwc(i) * ( bz(j,ib) / fl(j) +    &
                   ( bz(j+1,ib) / fl(j+1) - bz(j,ib) / fl(j) ) /    &
                   ( 1.0 / re(j+1) - 1.0 / re(j) ) * ( 1.0 / pre(i)   &
                   - 1.0 / re(j) ) )
             ww(i) = wz(j,ib) + ( wz(j+1,ib) - wz(j,ib) ) /   &
                   ( re(j+1) - re(j) ) * ( pre(i) - re(j) )
             gg = gz(j,ib) + ( gz(j+1,ib) - gz(j,ib) ) /   &
               ( re(j+1) - re(j) ) * ( pre(i) - re(j) )
             x1 = gg
             x2 = x1 * gg
             x3 = x2 * gg
             x4 = x3 * gg
             www(i,1) = 3.0 * x1
             www(i,2) = 5.0 * x2
             www(i,3) = 7.0 * x3
             www(i,4) = 9.0 * x4
           endif
         endif

10      continue
      return
      end

      subroutine rayle2 
!  *******************************************************************
!  trp is P(mb)/T(K)*DZ(m) and the constant 14.6337=R(287)/g(9.806)/2.
!  *******************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,trp
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /ray2/ trp(nvx)

      do 100 i = 1, nv
       trp(i) = 14.6337 * ( pp(i) + pp(i+1) )   &
                * alog( pp(i+1) / pp(i) ) 
100    continue
      return
      end

      subroutine rayle ( ib, u0, ig ,lray)
!  ****************************************************************
!  tr, wr, and wwr are the optical depth, single scattering albedo,
!  and expansion coefficients of the phase function ( 1, 2, 3, and
!  4 ) due to the Rayleigh scattering for a given layer.
!  ****************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
	logical lray
      integer i,ib,ig
      real u0,ri,trp,tr,wr,wwr,x,rix,riy

      common /ray1/ ri(mbsx),rix(10),riy(10)
      common /ray2/ trp(nvx)
      common /ray/ tr(nvx), wr(nvx), wwr(nvx,4)

      if ( ib .le. mbs ) then
        if ( ib .eq. 1 ) then
           if (lray) then
	     x = riy(ig) !LOG & Lin
	    else
	     x = -3.902860e-6*u0*u0+6.120070e-6*u0+4.177440e-6
!	     x= rix(ig) !LIN 
	   endif
        else
          x = ri(ib)
        endif
!!!!	x=1.0E-09 !!! ZERO RAYLEIGH
        do i = 1, nv
         tr(i) = trp(i) * x
         wr(i) = 1.0
         wwr(i,1) = 0.0
         wwr(i,2) = 0.5
         wwr(i,3) = 0.0
         wwr(i,4) = 0.0
         end do
      else
        do i = 1, nv
         tr(i) = 0.0
         wr(i) = 0.0
         wwr(i,1) = 0.0
         wwr(i,2) = 0.0
         wwr(i,3) = 0.0
         wwr(i,4) = 0.0
         end do
      endif

      return
      end

      subroutine rain ( ib )
!  *******************************************************************
!  trn, wrn, and wwrn are the optical depth, single scattering albedo,
!  and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
!  due to the Mie scattering of rain for a given layer. 
!                        Jan. 19, 1993
!  *******************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib
      real prwc,rwc,brn,wrnf,grn,dz,trn,wrn,wwrn
      real x1,x2,x3,x4,y1,y2,y3,y4
      common /rains/ prwc(nvx)
      common /rai1/ rwc, brn(mbx), wrnf(mbx), grn(mbx)
      common /thick/ dz(nvx)
      common /rai/ trn(nvx), wrn(nvx), wwrn(nvx,4)

      x1 = grn(ib)
      x2 = x1 * grn(ib)
      x3 = x2 * grn(ib)
      x4 = x3 * grn(ib)
      y1 = 3.0 * x1
      y2 = 5.0 * x2
      y3 = 7.0 * x3
      y4 = 9.0 * x4

      do 10 i = 1, nv
        if ( prwc(i) .lt. 1.0e-5 ) then
          trn(i) = 0.0
          wrn(i) = 0.0
          wwrn(i,1) = 0.0
          wwrn(i,2) = 0.0
          wwrn(i,3) = 0.0
          wwrn(i,4) = 0.0
        else
          trn(i) = dz(i) * prwc(i) * brn(ib) / rwc
          wrn(i) = wrnf(ib)
          wwrn(i,1) = y1
          wwrn(i,2) = y2
          wwrn(i,3) = y3
          wwrn(i,4) = y4
        endif
10       continue
      return
      end

      subroutine graup ( ib )
!  *******************************************************************
!  tgr, wgr, and wwgr are the optical depth, single scattering albedo,
!  and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
!  due to the Mie scattering of graupel for a given layer. 
!                        Jan. 19, 1993
!  *******************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib
      real pgwc,gwc,bg,wgf,gg,dz,tgr,wgr,wwgr,x1,x2,x3,x4,y1,y2,y3,y4
      common /graups/ pgwc(nvx)
      common /gra1/ gwc, bg(mbx), wgf(mbx), gg(mbx)
      common /thick/ dz(nvx)
      common /gra/ tgr(nvx), wgr(nvx), wwgr(nvx,4)

      x1 = gg(ib)
      x2 = x1 * gg(ib)
      x3 = x2 * gg(ib)
      x4 = x3 * gg(ib)
      y1 = 3.0 * x1
      y2 = 5.0 * x2
      y3 = 7.0 * x3
      y4 = 9.0 * x4

      do 10 i = 1, nv
       if ( pgwc(i) .lt. 1.0e-5 ) then
         tgr(i) = 0.0
         wgr(i) = 0.0
         wwgr(i,1) = 0.0
         wwgr(i,2) = 0.0
         wwgr(i,3) = 0.0
         wwgr(i,4) = 0.0
       else
         tgr(i) = dz(i) * pgwc(i) * bg(ib) / gwc
         wgr(i) = wgf(ib)
         wwgr(i,1) = y1
         wwgr(i,2) = y2
         wwgr(i,3) = y3
         wwgr(i,4) = y4
       endif
10     continue

      return
      end



!---------- 4/1/97 (5) -- PREVIOUS 465 LINES

      subroutine ql_rh(rh,tl,pl,ql)

!      rh (0-100)
!      tl (K)
!      pl (mb)
!      ql (g/g) mixing ratio, defined as the water vapor
!         content per unit mass of dry air, w = 0.622*e/(p-e)
!      
      es=satvap(tl)      
      ws=0.622*es/(pl-es)
      rh= ql/ws *100.
      return
      end
      
      real function satvap(temp2)
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
!        toot = 273.16 / temp2
!        toto = 1 / toot
!        eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) /    &
!          log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
!        satvap = 10 ** eilog
        satvap = esif(temp)/100.0
      else
!        tsot = 373.16 / temp2
!        ewlog = -7.90298 * (tsot - 1) + 5.02808 *    &
!                   (log(tsot) / log(10.))
!        ewlog2 = ewlog - 1.3816e-07 *    &
!                   (10 ** (11.344 * (1 - (1 / tsot))) - 1)
!        ewlog3 = ewlog2 + .0081328 *    &
!                   (10 ** (-3.49149 * (tsot - 1)) - 1)
!        ewlog4 = ewlog3 + (log(1013.246) / log(10.))
!        satvap = 10 ** ewlog4
         satvap = eslf(temp)/100.0
      end if

      return
      end

      subroutine gascon ( ib )
!  ********************************************************************
!  tgm(nv) are the optical depthes due to water vapor continuum absorp-
!  tion in nv layers for a given band ib. We include continuum absorp-
!  tion in the 280 to 1250 cm**-1 region. vv(11)-vv(17) are the central
!  wavenumbers of each band in this region. 
! *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib
      real tgm,vv(18)

      common /con/ tgm(nvx)

      data vv / 10*0.0, 1175.0, 1040.0, 890.0, 735.0,    &
                      605.0, 470.0, 340.0, 0.0 /

      if ( ib .gt. 10 .and. ib .lt. 18 ) then
         call qopcon ( vv(ib), tgm )
      else
         do i = 1, nv
          tgm(i) = 0.0
          end do
      endif

      return
      end

      subroutine gases ( ib, ig, hk )
!  *****************************************************************
!  tg(nv) are the optical depthes due to nongray gaseous absorption, 
!  in nv layers for a given band ib and cumulative probability ig. 
!  *****************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib,ig
      real hk,pp,pt,ph,po,umco2,umch4,umn2o,hk1, fk1o3
      real hk2,c2h2o,hk3,c3h2o,hk4,c4h2o,hk5,c5h2o,hk6,c6h2o
      real hk7,c7h2o,hk8,c8h2o,hk9,c9h2o,hk10,c10h2o,c10ch4,c10n2o
      real hk11,c11h2o,c11ch4,c11n2o,hk12,c12o3,c12h2o,hk13,c13h2o
      real hk14,c14hca,c14hcb,hk15,c15hca,c15hcb,hk16,c16h2o
      real hk17,c17h2o,hk18,c18h2o,tg,fk
      real fkg(nv1x), fkga(nv1x), fkgb(nv1x), pq(nv1x)
      real tg1(nvx), tg2(nvx), tg3(nvx)
      real sol_spect,fk1h2o

      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /umcon/ umco2, umch4, umn2o 

      common /band1/ hk1(10), fk1o3(10),sol_spect(0:7),fk1h2o(10)
      common /band2/ hk2(8), c2h2o(3,11,8)
      common /band3/ hk3(12), c3h2o(3,11,12)
      common /band4/ hk4(7), c4h2o(3,11,7)
      common /band5/ hk5(12), c5h2o(3,11,12)
      common /band6/ hk6(5), c6h2o(3,11,5)
      common /band7/ hk7(2), c7h2o(3,19,2)
      common /band8/ hk8(3), c8h2o(3,19,3)
      common /band9/ hk9(4), c9h2o(3,19,4)
      common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
      common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
      common /band12/ hk12(5), c12o3(3,19,5), c12h2o(3,19)
      common /band13/ hk13(2), c13h2o(3,19,2)
      common /band14/ hk14(10), c14hca(3,19,10), c14hcb(3,19,10)
      common /band15/ hk15(12), c15hca(3,19,12), c15hcb(3,19,12)
      common /band16/ hk16(7), c16h2o(3,19,7)
      common /band17/ hk17(7), c17h2o(3,19,7)
      common /band18/ hk18(8), c18h2o(3,19,8)
      common /gas/ tg(nvx)
	integer idkfr
	common /dkfrwn/ idkfr
	logical lchou,lband6a,lpar,lray
	common/CHOU/lchou,lband6a,lpar,lray

      goto ( 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 ) ib 
      stop
!  ---------------------------------------------------------------------
!  In this band ( 50000 - 14500 cm**-1 ), we have considered the nongray
!  gaseous absorption of O3.    619.618 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
 1    fk = fk1o3(ig)
!!!	fk=1.0E-09 !ZERO OZONE
        call qopo3s ( fk, tg )
	
	if (lpar) call qoph2o_chou ( fk1h2o(ig)  , tg )
!      hk = 619.618 * hk1(ig)
       hk = sol_spect(1) * hk1(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 14500 - 7700 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.  484.295 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
 2    call qks ( c2h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
!      hk = 484.295 * hk2(ig)
       hk = sol_spect(2) * hk2(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 7700 - 5250 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 149.845 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
 3    call qks ( c3h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
!      hk = 149.845 * hk3(ig)
       hk = sol_spect(3) * hk3(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 5250 - 4000 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 48.7302 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
 4    call qks ( c4h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
!      hk = 48.7302 * hk4(ig)
       hk = sol_spect(4) * hk4(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 4000 - 2850 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 31.6576 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
 5    call qks ( c5h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
	if ( lband6a ) then
	 hk = ( sol_spect(5) + sol_spect(7) ) * hk5(ig)
	else
!	 hk = 31.6576 * hk5(ig)
	 hk = sol_spect(5) * hk5(ig)
	endif
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 2850 - 2500 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 5.79927 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
 6    call qks ( c6h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
!      hk = 5.79927 * hk6(ig)
       hk = sol_spect(6) * hk6(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 2200 - 1900 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
 7    call qki ( c7h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
      hk = hk7(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 1900 - 1700 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
 8    call qki ( c8h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
      hk = hk8(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 1700 - 1400 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
 9    call qki ( c9h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
      hk = hk9(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 1400 - 1250 cm**-1 ), we have considered the 
!  overlapping absorption of H2O, CH4, and N2O by approach one of 
!  Fu(1991).
!  ---------------------------------------------------------------------
 10   call qki ( c10h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg1 )
      call qki ( c10ch4, fkg )
      call qopch4 ( fkg, tg2 )
      call qki ( c10n2o, fkg )
      call qopn2o ( fkg, tg3 )
      do i = 1, nv
       tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
       end do
      hk = hk10(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 1250 - 1100 cm**-1 ), we have considered the 
!  overlapping absorption of H2O, CH4, and N2O by approach one of 
!  Fu(1991).
!  ---------------------------------------------------------------------
 11   continue
	if ( idkfr == 0) then

      call qki ( c11h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg1 )

      call qki ( c11ch4, fkg )
      call qopch4 ( fkg, tg2 )
      call qki ( c11n2o, fkg )
      call qopn2o ( fkg, tg3 )
      do i = 1, nv
       tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
       end do

       hk = hk11(ig)
	
	elseif( idkfr == 1) then
	call ck_dkfr_win_all(ib,ig, hk,tg)
	elseif( idkfr == 2) then
	call ck_dkfr_win_hyb(ib,ig, hk,tg)
	endif

      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 1100 - 980 cm**-1 ), we have considered the overlapping
!  absorption of H2O and O3 by approach one of Fu(1991).
!  ---------------------------------------------------------------------
 12   continue
	if ( idkfr == 0) then
      call qkio3 ( c12o3(1,1,ig), fkg )
      call qopo3i ( fkg, tg1 )
      call qki ( c12h2o, fkg )
      call qoph2o ( fkg, tg2 )	

       do i = 1, nv
       tg(i) = tg1(i) + tg2(i)
       end do

      hk = hk12(ig)

	elseif( idkfr == 1) then
	call ck_dkfr_win_all(ib,ig, hk,tg)
	elseif( idkfr == 2) then
	call ck_dkfr_win_hyb(ib,ig, hk,tg)
	endif

      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 980 - 800 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
 13   continue
	if ( idkfr == 0) then
      call qki ( c13h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
      hk = hk13(ig)
	elseif( idkfr == 1) then
        call ck_dkfr_win_all(ib,ig, hk,tg)
	elseif( idkfr == 2) then
	call ck_dkfr_win_hyb(ib,ig, hk,tg)
	endif
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 800 - 670 cm**-1), we have considered the overlapping
!  absorption of H2O and CO2 by approach two of Fu(1991).
!  ---------------------------------------------------------------------
 14   do i = 1, nv1
       if ( pp(i) .ge. 63.1 ) then
         pq(i) = ph(i)
       else
         pq(i) = 0.0
       endif
       end do
      call qki ( c14hca(1,1,ig), fkga )
      call qki ( c14hcb(1,1,ig), fkgb )
      do i = 1, nv1
       fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
       end do
      call qophc ( fkg, tg)
      hk = hk14(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 670 - 540 cm**-1), we have considered the overlapping
!  absorption of H2O and CO2 by approach two of Fu(1991).
!  ---------------------------------------------------------------------
 15   do i = 1, nv1
       if ( pp(i) .ge. 63.1 ) then
         pq(i) = ph(i)
       else
         pq(i) = 0.0
       endif
       end do
      call qki ( c15hca(1,1,ig), fkga )
      call qki ( c15hcb(1,1,ig), fkgb )
      do i = 1, nv1
       fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
       end do
      call qophc ( fkg, tg)
      hk = hk15(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 540 - 400 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
 16   call qki ( c16h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
      hk = hk16(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 400 - 280 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
 17   call qki ( c17h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
      hk = hk17(ig)
      goto 20

!  ---------------------------------------------------------------------
!  In this band ( 280 - 000 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
 18   call qki ( c18h2o(1,1,ig), fkg )
      call qoph2o ( fkg, tg )
      hk = hk18(ig)

 20   continue
      return
      end
! 338kfix
	subroutine qks ( coefks, fkg )
! *********************************************************************
! fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
! **-1 for a given cumulative probability in nv1 layers. coefks(3,11)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 11 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at conditions other than those eleven
! pressures is interpolated linearly with pressure (Fu, 1991).
! *********************************************************************
      USE RadParams
	
!##        include 'rad_0698.h'
	common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
	dimension coefks(3,11)
	dimension fkg(nv1x)
	dimension stanp(11)
	dimension tk(nv1x)
      	data stanp / 10.0, 15.8, 25.1, 39.8, 63.1, 100.0,   &
      	             158.0, 251.0, 398.0, 631.0, 1000.0 /
	do i = 1, nv1
           if ( pt(i) .le. 180.0 ) then
              tk(i) = 180.0
           elseif ( pt(i) .ge. 320.0 ) then
              tk(i) = 320.0
           else
              tk(i) = pt(i)
           endif
	enddo
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefks(1,1) + coefks(2,1) * ( tk(i) - 245.0 )   &
             + coefks(3,1) * ( tk(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(11) ) then
	     y1 = ( tk(i) - 245.0 ) * ( tk(i) - 245.0 )
	     x1 = exp ( coefks(1,10) + coefks(2,10) * ( tk(i) - 245.0 )   &
      	     + coefks(3,10) * y1 )
    	     x2 = exp ( coefks(1,11) + coefks(2,11) * ( tk(i) - 245.0 )   &
      	     + coefks(3,11) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(11) - stanp(10) )   &
      	     * ( pp(i) - stanp(10) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( tk(i) - 245.0 ) * ( tk(i) - 245.0 )
	     x1 = exp ( coefks(1,i1-1) + coefks(2,i1-1) * (tk(i)-245.0)   &
      	     + coefks(3,i1-1) * y1 )
	     x2 = exp ( coefks(1,i1) + coefks(2,i1) * ( tk(i) - 245.0 )   &
      	     + coefks(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )   &
      	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

	subroutine qki ( coefki, fkg )
! *********************************************************************
! fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
! **-1 for a given cumulative probability in nv1 layers. coefki(3,19)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 19 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).
! *********************************************************************
      USE RadParams
!##	include 'rad_0698.h'
	common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
	dimension coefki(3,19)
	dimension fkg(nv1x)
	dimension stanp(19)
	dimension tk(nv1x)
	data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51,    &
      	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,   &
      	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	do i = 1, nv1
           if ( pt(i) .le. 180.0 ) then
              tk(i) = 180.0
           elseif ( pt(i) .ge. 320.0 ) then
              tk(i) = 320.0
           else
              tk(i) = pt(i)
           endif
	enddo
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( tk(i) - 245.0 )   &
             + coefki(3,1) * ( tk(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( tk(i) - 245.0 ) * ( tk(i) - 245.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( tk(i) - 245.0 )   &
      	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( tk(i) - 245.0 )   &
      	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )   &
      	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( tk(i) - 245.0 ) * ( tk(i) - 245.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (tk(i)-245.0)   &
      	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( tk(i) - 245.0 )   &
      	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )   &
      	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

	subroutine qkio3 ( coefki, fkg )
! *********************************************************************
! fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
! **-1 for a given cumulative probability in nv1 layers. coefki(3,19)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 19 pressures by
!         ln k = a + b * ( t - 250 ) + c * ( t - 250 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).
! *********************************************************************
      USE RadParams
!##	include 'rad_0698.h'
	common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
	dimension coefki(3,19)
	dimension fkg(nv1x)
	dimension stanp(19)
	dimension tk(nv1x)
	data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51,    &
      	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,   &
      	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	do i = 1, nv1
           if ( pt(i) .le. 180.0 ) then
              tk(i) = 180.0
           elseif ( pt(i) .ge. 320.0 ) then
              tk(i) = 320.0
           else
              tk(i) = pt(i)
           endif
	enddo
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( tk(i) - 250.0 )   &
             + coefki(3,1) * ( tk(i) - 250.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( tk(i) - 250.0 ) * ( tk(i) - 250.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( tk(i) - 250.0 )   &
      	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( tk(i) - 250.0 )   &
      	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )   &
      	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( tk(i) - 250.0 ) * ( tk(i) - 250.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (tk(i)-250.0)   &
      	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( tk(i) - 250.0 )   &
      	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )   &
      	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

! 338kfix

      subroutine qopo3s ( fk, tg )
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,fk,tg(nvx),fq

      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      fq = 238.08 * fk
      do 10 i = 1, nv
       tg(i) = ( po(i) + po(i+1) ) * ( pp(i+1) - pp(i) ) * fq
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 476.16 * fk
!20      continue
! 476.16 = 2.24e4 / M * 10.0 / 9.8, where M = 48 for O3.

      return
      end

      subroutine qoph2o ( fkg, tg )
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,fkg(nv1x),tg(nvx)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      do 10 i = 1, nv
       tg(i) = ( fkg(i) * ph(i) + fkg(i+1) * ph(i+1) )   &
               * ( pp(i+1) - pp(i) ) * 634.9205
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 1269.841
!20      continue
! 1269.841 = 2.24e4 / M * 10.0 / 9.8, where M = 18 for H2O.

      return
      end


      subroutine qoph2o_chou ( fk, tg )
      USE RadParams
      implicit none
      integer i
      real pp,pt,ph,po,fk,tg(nvx)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

!        do  i = 1, nv
!       tg(i) = tg(i) +  fk* 1.02
!     & * (ph(i)+ph(i+1)) *0.5 
!     & * ( pp(i+1)-pp(i) )
!     & * (.5*(pp(i)+pp(i+1))/300.)**.8 
!     & * (1.+0.00135*(pt(i)-240.)) +1.e-11 
!	enddo

       do  i = 1, nv
       tg(i) = tg(i) + fk*0.51 *(ph(i)+ ph(i+1))* (pp(i+1)-pp(i)) 
       enddo 

      return
      end

      subroutine qopch4 ( fkg, tg )
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,fkg(nv1x),tg(nvx)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      do 10 i = 1, nv
       tg(i) = ( fkg(i)+fkg(i+1) ) * ( pp(i+1)-pp(i) ) * 6.3119e-4
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 1.26238e-3
!20      continue
! 1.26238e-3 = 2.24e4 / M * 10.0 / 9.8 * 1.6e-6 * M / 28.97, where 
! M = 16 for CH4.

      return
      end

      subroutine qopn2o ( fkg, tg )
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,fkg(nv1x),tg(nvx)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      do 10 i = 1, nv
       tg(i) = ( fkg(i)+fkg(i+1) ) * ( pp(i+1)-pp(i) ) * 1.10459e-4
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 2.20918e-4
!20      continue
! 2.20918e-4 = 2.24e4 / M * 10.0 / 9.8 * 0.28e-6 * M / 28.97, where
! M = 44 for N2O.

      return
      end

      subroutine qopo3i ( fkg, tg )
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,fkg(nv1x),tg(nvx)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      do 10 i = 1, nv
       tg(i) = ( fkg(i) * po(i) + fkg(i+1) * po(i+1) )   &
               * ( pp(i+1) - pp(i) ) * 238.08
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 476.16
!20      continue

      return
      end

      subroutine qophc ( fkg, tg )
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,fkg(nv1x),tg(nvx)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      do 10 i = 1, nv
       tg(i) = ( fkg(i) + fkg(i+1) ) * ( pp(i+1) - pp(i) ) * 0.5
10     continue
!  ------------------------
!  See page 86 of Fu (1991).
!  ------------------------

      return
      end

      subroutine qopcon ( vv, tg )
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,ff(nv1x),pe(nv1x),tg(nvx),x,y,z,r,s,vv,w

      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      x = 4.18
      y = 5577.8
      z = 0.00787
      r = 0.002
      s = ( x + y * exp ( - z * vv ) ) / 1013.25

      do 3 i = 1, nv1
       pe(i) = pp(i) * ph(i) / ( 0.622 + 0.378 * ph(i) )
       w = exp ( 1800.0 / pt(i) - 6.08108 )
       ff(i) = s * ( pe(i) + r * pp(i) ) * w
3      continue

      do 5 i = 1, nv
       tg(i) = ( ff(i) * ph(i) + ff(i+1) * ph(i+1) )*   &
               ( pp(i+1) - pp(i) ) * 0.5098835
5      continue

!      do 7 i = 1, nv
!       tg(i) = tg(i) * 10.0 / 9.80616
!7      continue

      return
      end

      subroutine comscp
!  *********************************************************************
!  This subroutine is used to  COMbine Single-Scattering Properties  due
!  to  ice crystals,  water droplets, and  Rayleigh molecules along with
!  H2O continuum absorption and nongray gaseous absorption.  See Section
!  3.4 of Fu (1991). wc, wc1, wc2, wc3, and wc4, are total (or combined)
!  single - scattering  albedo,  and   expansion   coefficients  of  the
!  phase function ( 1, 2, 3, and 4 ) in nv layers. tt(nv) are the normal
!  optical depth ( from the top of the atmosphere to a given level ) for
!  level 2 - level nv1( surface ). The single-scattering  properties  of
!  rain and graupel are also incorporated in ( Jan. 19, 1993 ).
!  The single-scattering properties of aerosols are incorporated in 
!  (10/29/96) based on earlier version (5/17/95).
!  *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
 
      integer i,j,iac
      real ti,wi,wwi,tw,ww,www,trn,wrn,wwrn,tgr,wgr,wwgr
      real tr,wr,wwr,tgm,tg,wc1,wc2,wc3,wc4,wc,tt,tc(nvx)
      real tis,tws,trns,tgrs,fw

!---------- 10/29/96 (4)
      real tae,wae,wwae,taes(0:4)
      common /aer/ tae(nvx,mxac), wae(nvx,mxac), wwae(nvx,4,mxac)
!---------- 10/29/96 (4)
      
      common /ic/ ti(nvx), wi(nvx), wwi(nvx,4)
      common /wat/ tw(nvx), ww(nvx), www(nvx,4)
      common /rai/ trn(nvx), wrn(nvx), wwrn(nvx,4)
      common /gra/ tgr(nvx), wgr(nvx), wwgr(nvx,4)
      common /ray/ tr(nvx), wr(nvx), wwr(nvx,4)
      common /con/ tgm(nvx)
      common /gas/ tg(nvx)
      common /dfsin/ wc1(nvx), wc2(nvx), wc3(nvx), wc4(nvx),    &
                           wc(nvx), tt(nvx)

      LEVELS : do  i = 1, nv

       tc(i) = ti(i) + tw(i) + tr(i) + tgm(i) + tg(i) +   &
              trn(i) +tgr(i)

	do iac = 1,nac
	tc(i) = tc(i) + tae(i,iac)
	enddo

!---------- 10/29/96 (5)

       tis = ti(i) * wi(i)
       tws = tw(i) * ww(i) 
       trns = trn(i) * wrn(i)
       tgrs = tgr(i) * wgr(i)

!---------- 10/29/96 (6)
	taes(0:4) = 0.0
	do iac = 1,nac
         taes(0)=taes(0)+tae(i,iac)*wae(i,iac)
	do j=1,4
	 taes(j)=taes(j)+tae(i,iac)*wae(i,iac)*wwae(i,j,iac)
	enddo
	enddo

        fw = tis + tws + tr(i) + trns + tgrs + taes(0)
!---------- 10/29/96 (6)
!WJ, tc(i) ==0 will lead wc(i) out of flow.
                                                                                
       if ( tc(i) .eq. 0 ) then
           tc(i) = 1.0e-15
       endif

       wc(i) =  fw / tc(i)

       if ( fw .lt. 1.0e-20 ) then
        wc1(i) = 0.0 ; wc2(i) = 0.0 ; wc3(i) = 0.0 ; wc4(i) = 0.0
       else
!---------- 10/29/96 (7)

!         wc1(i) = ( tis * wwi(i,1) + tws * www(i,1) + 
!     1      tr(i) * wwr(i,1) + trns * wwrn(i,1) + tgrs * wwgr(i,1) +
!     1      taes * wwae(i,1) )/fw

!        wc2(i) = ( tis * wwi(i,2) + tws * www(i,2) +
!     1      tr(i) * wwr(i,2) + trns * wwrn(i,2) + tgrs * wwgr(i,2) +
!     1       taes * wwae(i,2) )/fw

!         wc3(i) = ( tis * wwi(i,3) + tws * www(i,3) +
!     1      tr(i) * wwr(i,3) + trns * wwrn(i,3) + tgrs * wwgr(i,3) +
!     1       taes * wwae(i,3) )/fw

!         wc4(i) = ( tis * wwi(i,4) + tws * www(i,4) +
!     1      tr(i) * wwr(i,4) + trns * wwrn(i,4) + tgrs * wwgr(i,4) +
!     1       taes * wwae(i,4) )/fw
!----------------------------------------------------------------------
        wc1(i) = ( tis * wwi(i,1) + tws * www(i,1) +    &
            tr(i) * wwr(i,1) + trns * wwrn(i,1) + tgrs * wwgr(i,1) +   &
            taes(1) )/fw

         wc2(i) = ( tis * wwi(i,2) + tws * www(i,2) +   &
            tr(i) * wwr(i,2) + trns * wwrn(i,2) + tgrs * wwgr(i,2) +   &
             taes(2) )/fw

         wc3(i) = ( tis * wwi(i,3) + tws * www(i,3) +   &
            tr(i) * wwr(i,3) + trns * wwrn(i,3) + tgrs * wwgr(i,3) +   &
             taes(3) )/fw

         wc4(i) = ( tis * wwi(i,4) + tws * www(i,4) +   &
            tr(i) * wwr(i,4) + trns * wwrn(i,4) + tgrs * wwgr(i,4) +   &
             taes(4) )/fw

!---------- 10/29/96 (7)
       endif

	enddo LEVELS



      tt(1) = tc(1)

      do 20 i = 2, nv
       tt(i) = tt(i-1) + tc(i)
20     continue

      return
      end

      subroutine planck ( ib, pts )
!  ******************************************************************
!  bf and bs are the blackbody intensity function integrated over the
!  band ib at the nv1 levels and at the surface, respectively.    The
!  units of bf and bs are W/m**2/Sr. nd*10 is the band width from ve.
!  ******************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib,nv11,ibr,j,nd(mbirx)
      real pts,pp,pt,ph,po,bf,bs,ve(mbirx),bt(nv1x)
      real bts,v1,v2,w,fq1,fq2,x

      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /planci/ bf(nv1x), bs

      data ve / 2200.0, 1900.0, 1700.0, 1400.0, 1250.0, 1100.0,   &
                  980.0, 800.0, 670.0, 540.0, 400.0, 280.001 /
      data nd / 30, 20, 30, 15, 15, 12,   &
                  18, 13, 13, 14, 12, 28 /

      nv11 = nv1 + 1
      ibr = ib - mbs
      bts = 0.0

      do 10 i = 1, nv1
       bt(i) = 0.0
10     continue
      v1 = ve(ibr)
      do 20 j = 1, nd(ibr)
       v2 = v1 - 10.0
       w = ( v1 + v2 ) * 0.5
       fq1 = 1.19107e-8 * w * w * w
       fq2 = 1.43884 * w
       do 30 i = 1, nv11
        if ( i .eq. nv11 ) then
          x = fq1 / ( exp ( fq2 / pts ) - 1.0 )
          bts = bts + x
        else
          x = fq1 / ( exp ( fq2 / pt(i) ) - 1.0 )
          bt(i) = bt(i) + x
        endif
30      continue
       v1 = v2
20     continue
      do 40 i = 1, nv1
       bf(i) = bt(i) * 10.0
40     continue
      bs = bts * 10.0

      return
      end

      subroutine coeff1 
!  *********************************************************************
!  coefficient calculations for four first-order differential equations.
!  *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,j,ib
      real a,u,p0d,p1d,p2d,p3d,p11d,p22d,p33d,x,w0w,w1w,w2w,w3w,fw
      real w,w1,w2,w3,t0,t1,u0,f0,b,c(4,5),q1,q2,q3,fq
      common /dis/ a(4)
      common /pointlg/ u(4)
      common /legen/ p0d(4), p1d(4), p2d(4), p3d(4)
      common /legen1/ p11d(4,4), p22d(4,4), p33d(4,4)
      common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf1/ b(4,3)

      x = 0.5 * w
      w0w = x
      w1w = x * w1
      w2w = x * w2
      w3w = x * w3
      if ( ib .le. mbs ) then
        fw = u0 * u0
        q1 = - w1w * u0
        q2 = w2w * ( 1.5 * fw - 0.5 )
        q3 = - w3w * ( 2.5 * fw - 1.5 ) * u0
      endif
      fq = 0.5 * w0w

      do 10 i = 3, 4
       do 20 j = 1, 4
        c(i,j) = fq + w1w * p11d(i,j) +   &
                 w2w * p22d(i,j) + w3w * p33d(i,j) 
        if ( i .eq. j ) then 
          c(i,j) = ( c(i,j) - 1.0 ) / u(i)
        else
          c(i,j) = c(i,j) / u(i)
        endif
20     continue
10    continue

      do 30 i = 1, 4
       if ( ib .le. mbs ) then
         c(i,5) = w0w + q1 * p1d(i) +   &
                  q2 * p2d(i) + q3 * p3d(i) 
       else
         c(i,5) = 1.0
       endif
       c(i,5) = c(i,5) / u(i)
30     continue

      b(1,1) = c(4,4) - c(4,1)
      b(1,2) = c(4,4) + c(4,1)
      b(2,1) = c(4,3) - c(4,2)
      b(2,2) = c(4,3) + c(4,2)
      b(3,1) = c(3,4) - c(3,1)
      b(3,2) = c(3,4) + c(3,1)
      b(4,1) = c(3,3) - c(3,2)
      b(4,2) = c(3,3) + c(3,2)
      b(1,3) = c(4,5) - c(1,5)
      b(2,3) = c(3,5) - c(2,5)
      b(3,3) = c(3,5) + c(2,5)
      b(4,3) = c(4,5) + c(1,5)
      return
      end

      subroutine coeff2 
!  *****************************************************************
!  coefficient calculations for second order differential equations.
!  *****************************************************************
      implicit none
      integer ib
      real w,w1,w2,w3,t0,t1,u0,f0,b,a,d,fw1,fw2,fw3,fw4
      common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf1/ b(4,3)
      common /coedf2/ a(2,2,2), d(4)

      fw1 = b(1,1) * b(1,2)
      fw2 = b(2,1) * b(3,2)
      fw3 = b(3,1) * b(2,2)
      fw4 = b(4,1) * b(4,2)
      a(2,2,1) = fw1 + fw2
      a(2,1,1) = b(1,1) * b(2,2) + b(2,1) * b(4,2)
      a(1,2,1) = b(3,1) * b(1,2) + b(4,1) * b(3,2)
      a(1,1,1) = fw3 + fw4
      a(2,2,2) = fw1 + fw3
      a(2,1,2) = b(1,2) * b(2,1) + b(2,2) * b(4,1)
      a(1,2,2) = b(3,2) * b(1,1) + b(4,2) * b(3,1)
      a(1,1,2) = fw2 + fw4
      d(1) = b(3,2) * b(4,3) + b(4,2) * b(3,3) + b(2,3) / u0
      d(2) = b(1,2) * b(4,3) + b(2,2) * b(3,3) + b(1,3) / u0
      d(3) = b(3,1) * b(1,3) + b(4,1) * b(2,3) + b(3,3) / u0
      d(4) = b(1,1) * b(1,3) + b(2,1) * b(2,3) + b(4,3) / u0

      return
      end

      subroutine coeff4 
!  *****************************************************************
!  coefficient calculations for fourth-order differential equations.
!  *****************************************************************
      implicit none
      integer ib
      real w,w1,w2,w3,t0,t1,u0,f0,a,d,b1,c1,z,x
      common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf2/ a(2,2,2), d(4)
      common /coedf4/ b1, c1, z(4)
      x = u0 * u0
      b1 = a(2,2,1) + a(1,1,1)
      c1 = a(2,1,1) * a(1,2,1) - a(1,1,1) * a(2,2,1)
      z(1) = a(2,1,1) * d(3) + d(4) / x - a(1,1,1) * d(4)
      z(2) = a(1,2,1) * d(4) - a(2,2,1) *d(3) + d(3) / x
      z(3) = a(2,1,2) * d(1) + d(2) / x - a(1,1,2) * d(2)
      z(4) = a(1,2,2) * d(2) - a(2,2,2) * d(1) + d(1) / x
      return
      end

      subroutine coeffl  
!  ********************************
!  fk1 and fk2 are the eigenvalues.
!  ********************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer ib,i
      real w,w1,w2,w3,t0,t1,u0,f0,a,d,b,b1,c1,z,x,aa,zz,a1,z1,fk1,fk2
      real dt,fw,a2,b2,fw1,fw2,y,zx,fq0,fq1
      
      common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf1/ b(4,3)
      common /coedf2/ a(2,2,2), d(4)
      common /coedf4/ b1, c1, z(4)
      common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2

      dt = t1 - t0
      x = sqrt ( b1 * b1 + 4.0 * c1 )
      fk1 = sqrt ( ( b1 + x ) * 0.5 )
      fk2 = sqrt ( ( b1 - x ) * 0.5 )
      fw = u0 * u0
      x = 1.0 / ( fw * fw ) - b1 / fw - c1
!---------- 4/2/97 (4)
      if (abs (x) .lt. 1.0E-16) THEN
        if ( x .lt. 0.0) THEN
          x = -1.0E-6
        else
          x = 1.0E-6
        end if
      end if
!---------- 4/2/97 (4)

      fw = 0.5 * f0 / x
      z(1) = fw * z(1) 
      z(2) = fw * z(2) 
      z(3) = fw * z(3) 
      z(4) = fw * z(4) 
      z1(1) = 0.5 * ( z(1) + z(3) )
      z1(2) = 0.5 * ( z(2) + z(4) )
      z1(3) = 0.5 * ( z(2) - z(4) )
      z1(4) = 0.5 * ( z(1) - z(3) )
      a2 = ( fk1 * fk1 - a(2,2,1) ) / a(2,1,1)
      b2 = ( fk2 * fk2 - a(2,2,1) ) / a(2,1,1)
      x = b(1,1) * b(4,1) - b(3,1) * b(2,1)
      fw1 = fk1 / x
      fw2 = fk2 / x
      y = fw2 * ( b2 * b(2,1) - b(4,1) ) 
      zx = fw1 * ( a2 * b(2,1) - b(4,1) )
      a1(1,1) = 0.5 * ( 1 - y )
      a1(1,2) = 0.5 * ( 1 - zx )
      a1(1,3) = 0.5 * ( 1 + zx )
      a1(1,4) = 0.5 * ( 1 + y )
      y = fw2 * ( b(3,1) - b2 * b(1,1) ) 
      zx = fw1 * ( b(3,1) - a2 * b(1,1) ) 
      a1(2,1) = 0.5 * ( b2 - y )
      a1(2,2) = 0.5 * ( a2 - zx )
      a1(2,3) = 0.5 * ( a2 + zx )
      a1(2,4) = 0.5 * ( b2 + y )
      a1(3,1) = a1(2,4)
      a1(3,2) = a1(2,3)
      a1(3,3) = a1(2,2)
      a1(3,4) = a1(2,1)
      a1(4,1) = a1(1,4)
      a1(4,2) = a1(1,3)
      a1(4,3) = a1(1,2)
      a1(4,4) = a1(1,1)

      if ( ib .le. mbs ) then
        fq0 = exp ( - t0 / u0 )
        fq1 = exp ( - t1 / u0 )
      else
        fq0 = 1.0
        fq1 = exp ( - dt / u0 )
      endif
      x = exp ( - fk1 * dt )
      y = exp ( - fk2 * dt )

      do 40 i = 1, 4
       zz(i,1) = z1(i) * fq0
       zz(i,2) = z1(i) * fq1
       aa(i,1,1) = a1(i,1)
       aa(i,2,1) = a1(i,2)
       aa(i,3,1) = a1(i,3) * x
       aa(i,4,1) = a1(i,4) * y
       aa(i,3,2) = a1(i,3)
       aa(i,4,2) = a1(i,4)
       aa(i,1,2) = a1(i,1) * y
       aa(i,2,2) = a1(i,2) * x
40     continue
      return
      end

      subroutine coefft 
!  ********************************************************************
!  See the paper by Liou, Fu and Ackerman (1988) for the formulation of
!  the delta-four-stream approximation in a homogeneous layer.
!  ********************************************************************
      call coeff1 
      call coeff2 
      call coeff4 
      call coeffl 
      return
      end

      subroutine coefft0 
!  *****************************************************************
!  In the limits of no scattering ( Fu, 1991 ), fk1 = 1.0 / u(3) and
!  fk2 = 1.0 / u(4).
!  *****************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib,jj,j,k
      real u,w,w1,w2,w3,t0,t1,u0,f0,aa,zz,a1,z1,fk1,fk2,y,fw,dt,x

      common /pointlg/ u(4)
      common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
      fk1 = 4.7320545
      fk2 = 1.2679491
      y = exp ( - ( t1 - t0 ) / u0 )
      fw = 0.5 * f0

      do 10 i = 1, 4
       if ( ib .le. mbs ) then
         z1(i) = 0.0
         zz(i,1) = 0.0
         zz(i,2) = 0.0
       else
         jj = 5 - i

!!       Fix by Nair on 6/Dec/07.  Sudden breakdown 
!!       of radiation code, messes up everything
!!       without any warning.  Caused by division
!!       by zero error in the line 
!!       z1(i) = fw / ( 1.0 + u(jj) / u0 )
!!       Added the following quick fix, meed
!!       to revisit and figure out the exact problem
!!       Begin patch

         if(abs(abs(u(jj))-abs(u0)) .LT. 1E-9)then
            if(u(jj) .LT. 0.0)then
               u(jj) = u(jj) + 1E-8
            else
               u(jj) = u(jj) - 1E-8
            endif
            print*,'Detected possible DBZ, Attempting correction'
         endif
   
!!       End patch
         z1(i) = fw / ( 1.0 + u(jj) / u0 )
         zz(i,1) = z1(i) 
         zz(i,2) = z1(i) * y
       endif

!  ***************************************************************
!  ******************  10/24/96 **********************************
!  ***************************************************************
!  multiple references to 10 continue replaced by loops 11 and 12.
!  ***************************************************************
!  ***************************************************************
!  ***************************************************************
       do 11 j = 1, 4
        a1(i,j) = 0.0
        do 12 k = 1, 2
         aa(i,j,k) = 0.0
12       continue
11      continue
10     continue

      do 20 i = 1, 4
       j = 5 - i
       a1(i,j) = 1.0
20     continue

      dt = t1 - t0
      x = exp ( - fk1 * dt )
      y = exp ( - fk2 * dt )
      aa(1,4,1) = y
      aa(2,3,1) = x
      aa(3,2,1) = 1.0
      aa(4,1,1) = 1.0
      aa(1,4,2) = 1.0
      aa(2,3,2) = 1.0
      aa(3,2,2) = x
      aa(4,1,2) = y

      return
      end

      subroutine qccfe ( ib, asbs, ee )  
!  ********************************************************************
!  In the solar band  asbs is the surface albedo, while in the infrared
!  band asbs is  blackbody intensity emitted at the surface temperature
!  times surface emissivity.  In this subroutine, the delta-four-stream
!  is applied to nonhomogeneous atmospheres. See comments in subroutine
!  'qcfel' for array AB(13,4*n).
!  ********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib,j,k,n,n4,ibn,i8,kf,i1,i2,i3,j1,j2,j3,m1,m2,m18,m28
      real a,u,w1,w2,w3,w,t,u0,f0,wn,w1n,w2n,w3n,t0n,t1n
      real u0n,f0n,aa,zz,a1,z1,fk1t,fk2t,fk1,fk2,a4,z4,g4,ab,bx,xx
      real fu(4,4), wu(4),v1,v2,v3,asbs,ee,fw1,fw2

      common /dis/ a(4)
      common /pointlg/ u(4)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0(ndfsx), f0(ndfsx)
      common /coedfi/ ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
      common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4),   &
                        fk1t, fk2t
      common /qccfeo/ fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx),    &
                            z4(4,ndfsx), g4(4,ndfsx)
      common /qcfelc/ ab(13,ndfs4x), bx(ndfs4x), xx(ndfs4x)

      n = ndfs
      n4 = ndfs4

      do i = 1, n4
       do j = 1, 13
        ab(j,i) = 0.0
        end do
       end do

      ibn = ib
      wn = w(1)
      w1n = w1(1)
      w2n = w2(1)
      w3n = w3(1)
      t0n = 0.0
      t1n = t(1)
      u0n = u0(1)
      f0n = f0(1)
      if ( wn .ge. 0.999999 ) then
        wn = 0.999999
      endif
      if ( wn .le. 1.0e-4 ) then
        call coefft0 
        fk1(1) = fk1t
        fk2(1) = fk2t
      else
        call coefft 
        fk1(1) = fk1t
        fk2(1) = fk2t
      endif

      do 10 i = 1, 4
       z4(i,1) = z1(i)
       do 11 j = 1, 4
        a4(i,j,1) = a1(i,j)
11      continue
10     continue

      do 20 i = 1, 2
       bx(i) = - zz(i+2,1)
       i8 = i + 8
       do 21 j = 1, 4
        ab(i8-j,j) = aa(i+2,j,1)
21      continue
20     continue

      do 30 i = 1, 4
       wu(i) = zz(i,2)
       do 31 j = 1, 4
        fu(i,j) = aa(i,j,2)
31      continue
30     continue

      do 40 k = 2, n
       wn = w(k)
       w1n = w1(k)
       w2n = w2(k)
       w3n = w3(k)
       t0n = t(k-1)
       t1n = t(k)
       u0n = u0(k)
       f0n = f0(k)
       if ( wn .ge. 0.999999 ) then
         wn = 0.999999
       endif
       if ( wn .le. 1.0e-4 ) then
         call coefft0 
         fk1(k) = fk1t
         fk2(k) = fk2t
       else
         call coefft 
         fk1(k) = fk1t
         fk2(k) = fk2t
       endif

       do 50 i = 1, 4
        z4(i,k) = z1(i)
        do 51 j = 1, 4
         a4(i,j,k) = a1(i,j)
51       continue
50      continue

       kf = k + k + k + k
       i1 = kf - 5
       i2 = i1 + 3
       j1 = kf - 7
       j2 = j1 + 3
       i3 = 0
       do 60 i = i1, i2
        i3 = i3 + 1
        bx(i) = - wu(i3) + zz(i3,1)
        j3 = 0
        i8 = i + 8
        do 61 j = j1, j2
         j3 = j3 + 1
         ab(i8-j,j) = fu(i3,j3)
61      continue
        j3 = 0
        do 62 j = j2 + 1, j2 + 4
         j3 = j3 + 1
         ab(i8-j,j) = - aa(i3,j3,1)
62       continue
60      continue

       do 70 i = 1, 4
        wu(i) = zz(i,2)
        do 71 j = 1, 4
         fu(i,j) = aa(i,j,2)
71       continue
70      continue

40     continue

      if ( ib .le. mbs ) then
        v1 = 0.2113247 * asbs
        v2 = 0.7886753 * asbs
        v3 = asbs * u0(1) * f0(1) * exp ( - t(n) / u0(1) )
        m1 = n4 - 1
        m2 = n4
        m18 = m1 + 8
        m28 = m2 + 8
        fw1 = v1 * wu(3)
        fw2 = v2 * wu(4)
        bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
        bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
        do 80 j = 1, 4
         j1 = n4 - 4 + j
         fw1 = v1 * fu(3,j)
         fw2 = v2 * fu(4,j)
         ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
         ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
80       continue
      else
        v1 = 0.2113247 * ( 1.0 - ee )
        v2 = 0.7886753 * ( 1.0 - ee )
        v3 = asbs
        m1 = n4 - 1
        m2 = n4
        m18 = m1 + 8
        m28 = m2 + 8
        fw1 = v1 * wu(3)
        fw2 = v2 * wu(4)
        bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
        bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
        do 85 j = 1, 4
         j1 = n4 - 4 + j
         fw1 = v1 * fu(3,j)
         fw2 = v2 * fu(4,j)
         ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
         ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
85       continue
      endif

      call qcfel 
      do 90 k = 1, n
       j = k + k + k + k - 4
       do 91 i = 1, 4
        j = j + 1
        g4(i,k) = xx(j)
91      continue
90     continue
      return
      end

      subroutine qcfel 
!  **********************************************************************
!  1. `qcfel' is the abbreviation of ` qiu constants for each layer'.
!  2. The inhomogeneous atmosphere is divided into n adjacent homogeneous
!     layers where the  single scattering properties are constant in each
!     layer and allowed to vary from one to another. Delta-four-stream is
!     employed for each homogeneous layer. The boundary conditions at the
!     top and bottom of the atmosphere,  together with  continuity condi-
!     tions  at  layer interfaces lead to a system of algebraic equations
!     from which 4*n unknown constants in the problom can be solved.
!  3. This subroutine is used for solving the 4*n unknowns of A *X = B by
!     considering the fact that the coefficient matrix is a sparse matrix
!     with the precise pattern in this special problom.
!  4. The method is not different in principle from the general scheme of
!     Gaussian elimination with backsubstitution, but carefully optimized
!     so as to minimize arithmetic operations.  Partial  pivoting is used
!     to quarantee  method's numerical stability,  which will  not change
!     the basic pattern of sparsity of the matrix.
!  5. Scaling special problems so as to make  its nonzero matrix elements
!     have comparable magnitudes, which will ameliorate the stability.
!  6. a, b and x present A, B and X in A*X=B, respectively. and n4=4*n.
!  7. AB(13,4*n) is the matrix A in band storage, in rows 3 to 13; rows 1
!     and 2 and other unset elements should be set to zero on entry.
!  8. The jth column of A is stored in the jth column of the array AB  as
!     follows:
!             AB(8+i-j,j) = A(i,j) for max(1,j-5) <= i <= min(4*n,j+5).
!     Reversedly, we have
!            A(ii+jj-8,jj) = AB(ii,jj).
! **********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer n,n4,k,l,i,j,k44,m1,i0m1,i0,m18,i0f,n44,m,im1,n3,n2,n1
      integer m2,m3,m4,m1f,m28,m38,m48,ifq
      real ab,b,x,p,t,yy,xx

      common /qcfelc/ ab(13,ndfs4x), b(ndfs4x), x(ndfs4x)
      n = ndfs
      n4 = ndfs4

      do 200 k = 1, n - 1
       k44 = 4 * k - 4
       do 201 l= 1, 4
        m1 = k44 + l
        p = 0.0

        do 10 i = 8, 14 - l
         if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
           p = ab(i,m1)
           i0 = i
         endif
10       continue
        i0m1 = i0 + m1
        m18 = m1 + 8

!       **********************************************************
!       ****************** 10/24/96 ******************************
!       **********************************************************
!       Replaced "if i0=8 goto 20" statement with if "i0<>8 then".
!       line "20 replaced by the endif.
!       **********************************************************
!       **********************************************************
!       **********************************************************
        if ( i0 .ne. 8 ) then
          do 15 j = m1, m1 + 8 - l
           i0f = i0m1 - j
           m1f = m18 - j
           t = ab(i0f,j)
           ab(i0f,j) = ab(m1f,j)
           ab(m1f,j) = t
15         continue

          i0f = i0m1 - 8
          t = b(i0f)
          b(i0f) = b(m1)
          b(m1) = t
        endif

        yy = ab(8,m1)
        ab(8,m1) = 1.0
        do 25 j = m1 + 1, m1 + 8 - l
         m1f = m18 - j
         ab(m1f,j) = ab(m1f,j) / yy
25       continue
        b(m1) = b(m1) / yy

        do 30 i = 9, 14 - l
         xx = ab(i,m1)
         ab(i,m1) = 0.0
         im1 = i + m1
         do 31 j = m1 + 1, m1 + 8 - l
          ifq = im1 - j
          m1f = m18 - j
          ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
31        continue
         ifq = im1 - 8
         b(ifq) = b(ifq) - b(m1) * xx
30       continue

201     continue
200    continue

      n44 = n4 - 4
      do 300 l = 1, 3
       m1 = n44 + l
       p = 0.0
       do 40 i = 8, 12 - l
        if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
          p = ab(i,m1)
          i0 = i
        endif
40      continue
       i0m1 = i0 + m1
       m18 = m1 + 8

!      **********************************************************
!      ****************** 10/24/96 ******************************
!      **********************************************************
!      Replaced "if i0=8 goto 55" statement with if "i0<>8 then".
!      line "55 replaced by the endif.
!      **********************************************************
!      **********************************************************
!      **********************************************************
       if ( i0 .ne. 8 ) then
         do 50 j = m1, m1 + 4 - l
          i0f = i0m1 - j
          m1f = m18 - j
          t = ab(i0f,j)
          ab(i0f,j) = ab(m1f,j)
          ab(m1f,j) = t
50        continue
         i0f = i0m1 - 8
         t = b(i0f)
         b(i0f) = b(m1)
         b(m1) = t
       endif

       yy = ab(8,m1)
       ab(8,m1) = 1.0
       do 60 j = m1 + 1, m1 + 4 - l
        m1f = m18 - j
        ab(m1f,j) = ab(m1f,j) / yy
60      continue
       b(m1) = b(m1) / yy
       do 70 i = 9, 12 - l
        xx = ab(i,m1)
        ab(i,m1) = 0.0
        im1 = i + m1
        do 71 j = m1 + 1, m1 + 4 - l
         ifq = im1 - j
         m1f = m18 - j
         ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
71       continue
        ifq = im1 - 8
        b(ifq) = b(ifq) - b(m1) * xx
70      continue
300    continue

      yy = ab(8,n4)
      ab(8,n4) = 1.0
      b(n4) = b(n4) / yy
      n3 = n4 - 1
      n2 = n3 - 1
      n1 = n2 - 1
      x(n4) = b(n4)
      x(n3) = b(n3) - ab(7,n4) * x(n4)
      x(n2) = b(n2) - ab(7,n3) * x(n3) - ab(6,n4) * x(n4)
      x(n1) = b(n1) - ab(7,n2) * x(n2) - ab(6,n3) * x(n3) -   &
            ab(5,n4) * x(n4)
      do 80 k = 1, n - 1
       m4 = 4 * ( n - k )
       m3 = m4 - 1
       m2 = m3 - 1
       m1 = m2 - 1
       m48 = m4 + 8
       m38 = m3 + 8
       m28 = m2 + 8
       m18 = m1 + 8
       x(m4) = b(m4)
       do 81 m = m4 + 1, m4 + 4
        x(m4) = x(m4) - ab(m48-m,m) * x(m)
81      continue
       x(m3) = b(m3)
       do 82 m = m3 + 1, m3 + 5
        x(m3) = x(m3) - ab(m38-m,m) * x(m)
82      continue
       x(m2) = b(m2)
       do 83 m = m2 + 1, m2 + 6
        x(m2) = x(m2) - ab(m28-m,m) * x(m)
83      continue
       x(m1) = b(m1)
       do 84 m = m1 + 1, m1 + 7
        x(m1) = x(m1) - ab(m18-m,m) * x(m)
84      continue
80     continue
      return
      end

      subroutine adjust4
!  *****************************************************************
!  In this subroutine, we incorporate a delta-function adjustment to
!  account for the  forward  diffraction  peak in the context of the 
!  four-stream or two stream approximations ( Liou, Fu and Ackerman,
!  1988 ).  The w1(n), w2(n), w3(n), w(n), and t(n) are the adjusted
!  parameters.
!  *****************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,n
      real ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t,u0a, f0a
      real dtt(ndfsx), dt(ndfsx),tt0,f,fw

      common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),   &
                           ww(ndfsx), tt(ndfsx)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0a(ndfsx), f0a(ndfsx)

      n = ndfs
      tt0 = 0.0
      do i = 1, n
       f = ww4(i) / 9.0
       fw = 1.0 - f * ww(i) 
       w1(i) = ( ww1(i) - 3.0 * f ) / ( 1.0 - f )
       w2(i) = ( ww2(i) - 5.0 * f ) / ( 1.0 - f )
       w3(i) = ( ww3(i) - 7.0 * f ) / ( 1.0 - f )
       w(i) = ( 1.0 - f ) * ww(i) / fw
       dtt(i) = tt(i) - tt0
       tt0 = tt(i)
       dt(i) = dtt(i) * fw
       enddo

      t(1) = dt(1)
      do i = 2, n
       t(i) = dt(i) + t(i-1)
       enddo

      return
      end

      subroutine qfts ( ib, as, u0, f0 )
!  *******************************************************************
!  The delta-four-stream approximation for nonhomogeneous atmospheres
!  in the solar wavelengths (Fu, 1991). The input parameters are ndfs,
!  mdfs, and ndfs4 through 'para.file',  ib, as, u0, f0 for solar and
!  ib, bf, bs, ee for IR through arguments of  'qfts' and 'qfti', and
!  ww1(ndfs), ww2(ndfs), ww3(ndfs), ww4(ndfs), ww(ndfs), and tt(ndfs)
!  through common statement 'dfsin'.
!  *******************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer ib,n,m,i,k,jj,ii
      real a,u,ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t,u0a,f0a
      real fk1,fk2,a4,z4,g4,ffu,ffd,as,u0,f0,asbs
      real x(4), fi(4),ee,fw1,fw2,fw3,y,y1,fw4
!---------- 10/28/96 (5)
      real ffdr, ffdf, ydr
      common /dirdiff/ ffdr(nv1x), ffdf(nv1x)
!---------- 10/28/96 (5)

      common /dis/ a(4)
      common /pointlg/ u(4)
      common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),   &
                           ww(ndfsx), tt(ndfsx)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0a(ndfsx), f0a(ndfsx)
      common /qccfeo/ fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx),    &
                            z4(4,ndfsx), g4(4,ndfsx)
      common /dfsout/ ffu(mdfsx), ffd(mdfsx)

      n = ndfs
      m = mdfs
      ee = 0.0
      asbs = as
      call adjust4
      do 5 i = 1, n
       u0a(i) = u0
       f0a(i) = f0
5      continue

      call qccfe ( ib, asbs, ee ) 
      fw1 = 0.6638961
      fw2 = 2.4776962
      fw3 = u0 * 3.14159 * f0 
      do 10 i = 1, m
       if ( i .eq. 1 ) then
         x(1) = 1.0
         x(2) = 1.0
         x(3) = exp ( - fk1(1) * t(1) )
         x(4) = exp ( - fk2(1) * t(1) )
         k = 1
         y = 1.0
!---------- 10/28/96 (6a)
         ydr=1.0
!---------- 10/28/96 (6a)
       elseif ( i .eq. 2 ) then
         x(1) = exp ( - fk2(1) * t(1) )
         x(2) = exp ( - fk1(1) * t(1) )
         x(3) = 1.0
         x(4) = 1.0
         k = 1
         y = exp ( - t(1) / u0 )
!---------- 10/28/96 (6b)
         ydr = exp ( - tt(1) / u0 )
!---------- 10/28/96 (6b)
       else
         k = i - 1
         y1 = t(k) - t(k-1)
         x(1) = exp ( - fk2(k) * y1 )
         x(2) = exp ( - fk1(k) * y1 )
         x(3) = 1.0
         x(4) = 1.0
         y = exp ( - t(k) / u0 )
!---------- 10/28/96 (6c)
         ydr = exp ( - tt(k) / u0 )
!---------- 10/28/96 (6c)
       endif
       do 37 jj = 1, 4
        fi(jj) = z4(jj,k) * y
37      continue

       do 40 ii = 1, 4
        fw4 = g4(ii,k) * x(ii)
        do 41 jj = 1, 4
         fi(jj) = fi(jj) + a4(jj,ii,k) * fw4
41       continue
40      continue

       ffu(i)= fw1 * fi(2) + fw2 * fi(1) 
       ffd(i)= fw1 * fi(3) + fw2 * fi(4) + fw3 * y
!---------- 10/28/96 (7)
       ffdr(i) = fw3 * ydr
       ffdf(i) = ffd(i) - ffdr(i)
!---------- 10/28/96 (7)
10     continue
      return
      end

      subroutine qfti ( ib, ee )
! ******************************************************************
!  The exponential approximation for the Planck function in optical 
!  depth is used for the infrared (Fu, 1991). Since the direct solar 
!  radiation source has an exponential function form in terms of 
!  optical depth, the formulation of the delta-four-stream 
!  approximation for infrared wavelengths is the same as that for 
!  solar wavelengths. 
! ******************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,ib,n,k,ii,jj,m
      real a,u,ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t,u0,f0,fk1
      real fk2,a4,z4,g4,ffu,ffd,bf,bs,asbs,q1,q2,t0,fw1,fw2,xy
      real x(4), fi(4),fw3,ee,y1
	real timt0
      common /dis/ a(4)
      common /pointlg/ u(4)
      common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),   &
                           ww(ndfsx), tt(ndfsx)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0(ndfsx), f0(ndfsx)
      common /qccfeo/ fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx),    &
                            z4(4,ndfsx), g4(4,ndfsx)
      common /dfsout/ ffu(mdfsx), ffd(mdfsx)
      common /planci/ bf(nv1x), bs

      n = ndfs
      m = mdfs
      asbs = bs * ee
      call adjust4 
      t0 = 0.0

      do 3 i = 1, n
       q1 = alog ( bf(i+1) / bf(i) )

! Paul Stackhouse: July 5, 98
!          q2 = 1.0 / ( t(i) - t0 )
           timt0 = t(i) - t0 
           if (timt0 .lt. 1.0e-12) timt0 = 1.0e-12
           q2 = 1.0 / timt0
! Paul Stackhouse: July 5, 98

       f0(i) = 2.0 * ( 1.0 - w(i) ) * bf(i)
       if ( abs(q1) .le. 1.0e-10 ) then
         u0(i) = - 1.0e+10 / q2
       else
         u0(i) = - 1.0 / ( q1 * q2 )
       endif

!---------- 4/2/97 (5)
      if (abs(u0(i)) .gt. 4.25E+09) then
        if (u0(i) .lt. 0.0) then
          u0(i) = -4.25E+09
        else
          u0(i) = 4.25E+09
        end if
      end if
!---------- 4/2/97 (5)

       t0 = t(i)
3       continue

      call qccfe ( ib, asbs, ee ) 
      fw1 = 0.6638958
      fw2 = 2.4776962
      do 10 i = 1, m
       if ( i .eq. 1 ) then
         x(1) = 1.0
         x(2) = 1.0
         x(3) = exp ( - fk1(1) * t(1) )
         x(4) = exp ( - fk2(1) * t(1) )
         k = 1
         xy = 1.0
       elseif ( i .eq. 2 ) then
         x(1) = exp ( - fk2(1) * t(1) )
         x(2) = exp ( - fk1(1) * t(1) )
         x(3) = 1.0
         x(4) = 1.0
         k = 1
         xy =  exp ( - t(1) / u0(1) )
       else
         k = i - 1
         y1 = t(k) - t(k-1)
         x(1) = exp ( - fk2(k) * y1 )
         x(2) = exp ( - fk1(k) * y1 )
         x(3) = 1.0
         x(4) = 1.0
         xy =  exp ( - y1 / u0(k) )
       endif
       do 37 jj = 1, 4
        fi(jj) = z4(jj,k) * xy
37      continue
       do 40 ii = 1, 4
        fw3 = g4(ii,k) * x(ii)
        do 45 jj = 1, 4
         fi(jj) = fi(jj) + a4(jj,ii,k) * fw3
45       continue
40      continue
       ffu(i)= fw1 * fi(2) + fw2 * fi(1)
       ffd(i)= fw1 * fi(3) + fw2 * fi(4)
!      NAIR DEBUG
!       print*,'FFIR',fi(2),fi(1),fi(3),fi(4),ffu(i),ffd(i)
10     continue

      return
      end

      subroutine cfgts0 ( gamma1, gamma2, gamma3, gamma4, ugts1 )
!  *********************************************************************
!  This subroutine is used to calculate the Coefficients For Generalized
!  Two-Stream scheme. We can make choices between Eddington, quadrature
!  and  hemispheric  mean  schemes through  logical variables 'edding',
!  'quadra', and 'hemisp'.  The Eddington and quadrature schemes are 
!  discussed in detail by Liou (1992).  The hemispheric mean scheme is 
!  derived by assuming that the phase function is equal to 1 + g in the 
!  forward scattering hemisphere and 1 - g  in the backward scattering 
!  hemisphere where g is the asymmetry factor. The hemispheric mean is
!  only used for infrared wavelengths (Toon et al. 1989).
!
!   11/4/95
!
!  *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer ib
      real w,w1,w2,w3,t0,t1,u0,f0,x,y,gamma1,gamma2
      real gamma3,gamma4,z,ugts1

      common /coedfi / ib, w, w1, w2, w3, t0, t1, u0, f0
      if ( edding ) then
        x = 0.25 * w1
        y = w * x
        gamma1 = 1.75 - w - y
        gamma2 = - 0.25 + w - y
        gamma3 = 0.0
        gamma4 = 0.0
        if ( ib .le. mbs ) then
          gamma3 = 0.5 - x * u0
          gamma4 = 1.0 - gamma3
        endif
        ugts1 = 0.5
      endif

      if ( quadra ) then
        x = 0.866 * w
        y = 0.2887 * w1
        z = y * w
        gamma1 = 1.732 - x - z
        gamma2 = x - z
        gamma3 = 0.0
        gamma4 = 0.0
        if ( ib .le. mbs ) then
          gamma3 = 0.5 - y * u0
          gamma4 = 1.0 - gamma3
        endif
        ugts1 = 0.57735
      endif

      if ( hemisp ) then
        x = w * w1 / 3.0
        gamma1 = 2.0 - w - x
        gamma2 = w - x
        gamma3 = 0.0
        gamma4 = 0.0
        ugts1 = 0.5
      endif

      if ( mquadr ) then
        y = 0.2767 * w
        x = y + y + y
        z = y * w1
        gamma1 = 1.66 - x - z
        gamma2 = x - z
        gamma3 = 0.0
        gamma4 = 0.0
        ugts1 = 0.6024
      endif

      return
      end

      subroutine cfgts ( lamda, gamma, cadd0, cadd1, cmin0, cmin1,   &
                           g1g2, fkb )
!  *********************************************************
!  This subroutine is used to calculate the Coefficients For
!  Generalized Two-Stream scheme. 
!
!  11/4/95
!
!  ***********************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer ib
      real lamda,w,w1,w2,w3,t0,t1,u0,f0,gamma1,gamma2
      real gamma3,gamma4,ugts1,gamma,g1g2,alfa,beta,cadd0,cmin0
      real fw,cadd1,cmin1,fkb,x,z,fq

      common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0

      call cfgts0 ( gamma1, gamma2, gamma3, gamma4, ugts1 )
      lamda = sqrt ( ( gamma1 + gamma2 ) * ( gamma1 - gamma2 ) )
      gamma = gamma2 / ( gamma1 + lamda )
      g1g2 = gamma1 + gamma2
      fq = 1.0 / u0

      x = exp ( - fq * ( t1 - t0 ) )
      z = lamda * lamda - fq * fq

!---------- 4/2/97 (6)
      if (z .eq. 0.0) z = 0.001
!---------- 4/2/97 (6)

      fkb = z

      if ( ib .le. mbs ) then
        alfa = gamma3
        beta = gamma4
        fw = 3.1415927 * f0 * w * exp ( - fq * t0 )
        cadd0 = fw * ( ( gamma1 - fq ) * alfa +   &
                   beta * gamma2 ) / z
        cmin0 = fw * ( ( gamma1 + fq ) * beta +   &
                   alfa * gamma2 ) / z
      else
        fw = 3.1415927 * f0
        cadd0 = fw * ( g1g2 - fq ) / z
        cmin0 = fw * ( g1g2 + fq ) / z
      endif

      cadd1 = cadd0 * x
      cmin1 = cmin0 * x
      return
      end

! +++ PAUSE FIX      subroutine qccgts ( ib, asbs, ee )
      subroutine qccgts ( ib, asbs, ee, ITriErr )
!  ********************************************************************
!  In the solar band  asbs is the surface albedo, while in the infrared
!  band asbs is  blackbody intensity emitted at the surface temperature
!  times surface emissivity.  In this subroutine,  the generalized two-
!  stream is applied to nonhomogeneous atmospheres. ee is the IR
!  surface emissivity. 
!  ********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer ndfs2x,ib,k,ibn,k1,k2
      INTEGER ITriErr
      parameter ( ndfs2x = ndfsx * 2 )
      real lamdan,w1,w2,w3,w,t,u0,f0,wn,w1n,w2n
      real w3n,t0n,t1n,u0n, f0n,gamman,caddn,cminn,caddn0,cminn0
      real aa,bb,expn,g1g2n,fkbn,ee,asbs,wm1,wm2,rsfc,ssfc
      real a(ndfs2x), b(ndfs2x), c(ndfs2x), r(ndfs2x), u(ndfs2x)
      real xn(ndfsx), yn(ndfsx), zn(ndfsx),gam(ndfs2x)

      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0(ndfsx), f0(ndfsx)
      common /coedfi/ ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
      common / gtscoe / lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx),   &
                          cminn(ndfsx), caddn0(ndfsx), cminn0(ndfsx),   &
                      aa(ndfsx), bb(ndfsx), expn(ndfsx), g1g2n(ndfsx),   &
                          fkbn(ndfsx)

      ibn = ib
      do 40 k = 1, ndfs
       wn = w(k)
       w1n = w1(k)
       if ( k .eq. 1 ) then
         t0n = 0.0
       else
         t0n = t(k-1)
       endif
       t1n = t(k)
       u0n = u0(k)
       f0n = f0(k)
       if ( wn .ge. 0.999999 ) then
         wn = 0.999999
       endif
       call cfgts ( lamdan(k), gamman(k), caddn0(k), caddn(k),   &
                        cminn0(k), cminn(k), g1g2n(k) , fkbn(k))
       expn(k) = exp ( - lamdan(k) * ( t1n - t0n ) )
       xn(k) = gamman(k) * expn(k)
       yn(k) = ( expn(k) - gamman(k) ) / ( xn(k) - 1.0 )
       zn(k) = ( expn(k) + gamman(k) ) / ( xn(k) + 1.0 )
40     continue

      a(1) = 0.0
      b(1) = xn(1) + 1.0
      c(1) = xn(1) - 1.0
      r(1) = - cminn0(1)
      do 50 k = 1, ndfs - 1
       k1 = k + k
       k2 = k + k + 1
       a(k1) = 1.0 + xn(k) - yn(k+1) * ( gamman(k) + expn(k) )
       b(k1) = 1.0 - xn(k) - yn(k+1) * ( gamman(k) - expn(k) )
       c(k1) = yn(k+1) * ( 1.0 + xn(k+1) ) - expn(k+1) - gamman(k+1)
       r(k1) = caddn0(k+1) - caddn(k) - yn(k+1) *   &
                ( cminn0(k+1) - cminn(k) )
       a(k2) = gamman(k) - expn(k) - zn(k) * ( 1.0 - xn(k) )
       b(k2) = -1.0 - xn(k+1) + zn(k) * ( expn(k+1) + gamman(k+1) )
       c(k2) = zn(k) * ( expn(k+1) - gamman(k+1) ) - xn(k+1) + 1.0
       r(k2) = cminn0(k+1) - cminn(k) - zn(k) *   &
                ( caddn0(k+1) - caddn(k) )
50     continue
      if ( ib .le. mbs ) then
        rsfc = asbs
        ssfc = 3.1415927 * u0(1) * exp(-t(ndfs)/u0(1)) * rsfc * f0(1)
      else
        rsfc = 1.0 - ee
        ssfc = 3.1415927 * asbs
      endif

      wm1 = 1.0 - rsfc * gamman(ndfs)
      wm2 = xn(ndfs) - rsfc * expn(ndfs)
      a(ndfs2) = wm1 + wm2
      b(ndfs2) = wm1 - wm2
      c(ndfs2) = 0.0
      r(ndfs2) = rsfc * cminn(ndfs) - caddn(ndfs) + ssfc
! +++ PAUSE FIX      call tridag ( a, b, c, r, u, gam, ndfs2 )
      call tridag ( a, b, c, r, u, gam, ndfs2, ITriErr )
      IF (ITriErr .NE. 0) RETURN
      do 60 k = 1, ndfs
       k1 = k + k - 1
       k2 = k + k
       aa(k) = u(k1) + u(k2)
       bb(k) = u(k1) - u(k2)
60     continue
      return
      end

! +++ PAUSE FIX      subroutine tridag ( a, b, c, r, u, gam, n )
      subroutine tridag ( a, b, c, r, u, gam, n, ITriErr )
! *******************************************************************
!
!   | b1 c1 0  ...                |   | u1   |   | r1   |           
!   | a2 b2 c2 ...                |   | u2   |   | r2   |          
!   |          ...                | . | .    | = | .    |
!   |          ... an-1 bn-1 cn-1 |   | un-1 |   | rn-1 |          
!   |              0    an   bn   |   | un   |   | rn   |        
!
! This  subroutine solves for  a vector U of length N the tridiagonal
! linear set given by above equation. A, B, C and R are input vectors
! and are not modified (Numerical Recipes by Press et al. 1989).
! *******************************************************************
      implicit none
      integer n,j
      INTEGER ITriErr
      real gam(n), a(n), b(n), c(n), r(n), u(n),bet
      ITriErr = 0
! +++ PAUSE FIX      if ( b(1) .eq. 0. ) pause
      IF ( b(1) .EQ. 0. ) THEN
           ITriErr = 1
           RETURN
      END IF
!  *********************************************************
!  If this happens then you should rewrite your equations as
!  a set of order n-1, with u2 trivially eliminated.
!  *********************************************************
      bet = b(1)
      u(1) = r(1) / bet

!  ***************************************
!  Decomposition and forward substitution.
!  ***************************************
      do 11 j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j) * gam(j)
! +++ PAUSE FIX       if ( bet .eq. 0. ) pause
       IF ( bet .EQ. 0. ) THEN
           ITriErr = 2
           RETURN
       END IF
!      ---------------------------------------
!      Algorithm fails; see Numerical Recipes.
!      ---------------------------------------
       u(j) = ( r(j) - a(j) * u(j-1) ) / bet
11     continue

!  ****************
!  Backsubstitution
!  ****************
      do 12 j = n - 1, 1, -1
       u(j) = u(j) - gam(j+1) * u(j+1)
12     continue

      return
      end

! +++ PAUSE FIX      subroutine qftsts ( ib, as, u0, f0 )
      subroutine qftsts ( ib, as, u0, f0, ITriErr )
! **********************************************************************
! The generalized two stream approximation for nonhomgeneous atmospheres
! in  the  solar  wavelengths.  The  input  parameters are those through
! 'para.file', through argument of 'qftsts' and through common statement
! 'dfsin' and 'gtslog'.
! **********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      INTEGER ITriErr
      integer ib,n,m,i,k
      real ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t,u0a,f0a
      real lamdan,gamman,caddn,cminn,caddn0,cminn0,aa,bb,expn
      real g1g2n,fkbn,ffu,ffd,as,u0,f0,ee,asbs,fw3,xx
      real yy(ndfsx)

!---------- 10/28/96 (8)
      real ffdr,ffdf,yydr(ndfsx)
      common /dirdiff/ ffdr(nv1x), ffdf(nv1x)
!---------- 10/28/96 (8)
      common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),   &
                           ww(ndfsx), tt(ndfsx)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0a(ndfsx), f0a(ndfsx)
      common / gtscoe / lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx),   &
                          cminn(ndfsx), caddn0(ndfsx), cminn0(ndfsx),   &
                       aa(ndfsx), bb(ndfsx), expn(ndfsx), g1g2n(ndfsx),   &
                          fkbn(ndfsx)
      common /dfsout/ ffu(mdfsx), ffd(mdfsx)

      n = ndfs
      m = mdfs
      ee = 0.0
      asbs = as
      call adjust2

      do 5 i = 1, n
       u0a(i) = u0
       f0a(i) = f0
5      continue

! +++ PAUSE FIX      call qccgts ( ib, asbs, ee )
      call qccgts ( ib, asbs, ee, ITriErr )
      IF (ITriErr .NE. 0) RETURN
      fw3 = u0 * 3.1415927 * f0
      do k = 1, ndfs
       yy(k) = exp(-t(k)/u0)
!---------- 10/28/96 (9)
       yydr(k) = exp(-tt(k)/u0)
!---------- 10/28/96 (9)
       enddo

      xx = aa(1) * expn(1)
      ffu(1) = xx + gamman(1) * bb(1) + caddn0(1)
      ffd(1) = gamman(1) * xx + bb(1) + cminn0(1) + fw3
!---------- 10/28/96 (10)
      ffdr(1) =  fw3
      ffdf(1) = ffd(1) - ffdr(1)
!---------- 10/28/96 (10)
      do 10 i = 2, m
       k = i - 1
       xx = bb(k) * expn(k)
       ffu(i) = aa(k) + gamman(k) * xx + caddn(k)
       ffd(i) = gamman(k) * aa(k) + xx + cminn(k) + fw3 * yy(k)
!---------- 10/28/96 (11)
       ffdr(i) = fw3 * yydr(k)
       ffdf(i) = ffd(i) - ffdr(i)
!---------- 10/28/96 (11)
10     continue
      
      return
      end

! +++ PAUSE FIX      subroutine qftits ( ib, ee )
      subroutine qftits ( ib, ee, ITriErr )
! **********************************************************************
! The exponential approximation for the Planck function in optical depth
! is used for the infrared ( Fu, 1991). Since the direct solar radiation
! source has an exponential function form in terms of optical depth, the
! formulation of generalized two stream approximation for infrared  wave
! lengths is the same as that for solar wavelengths. 
! The generalized two stream approximation for nonhomgeneous atmospheres
! in the infrared wavelengths.  The  input  parameters are those through
! 'para.file', through argument of 'qftits' and through common statement
! 'dfsin', 'gtslog', and 'planci'.
! **********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      INTEGER ITriErr
      integer ib,n,m,i,k
      real lamdan,ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t
      real gamman,caddn,cminn,caddn0,cminn0,aa,bb,expn
      real g1g2n,fkbn,ffu,ffd,u0,f0,ee,asbs,xx,q1,q2,bf,bs,t0
	real timt0
      common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),   &
                           ww(ndfsx), tt(ndfsx)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0(ndfsx), f0(ndfsx)
      common / gtscoe / lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx),   &
                          cminn(ndfsx), caddn0(ndfsx), cminn0(ndfsx),   &
                       aa(ndfsx), bb(ndfsx), expn(ndfsx), g1g2n(ndfsx),   &
                          fkbn(ndfsx)
      common /dfsout/ ffu(mdfsx), ffd(mdfsx)
      common /planci/ bf(nv1x), bs
      n = ndfs
      m = mdfs
      asbs = bs * ee
      call adjust2
      t0 = 0.0

      do 3 i = 1, n
       q1 = alog ( bf(i+1) / bf(i) )

! Paul Stackhouse: July 5, 98
!          q2 = 1.0 / ( t(i) - t0 )
           timt0 = t(i) - t0 
           if (timt0 .lt. 1.0e-12) timt0 = 1.0e-12
           q2 = 1.0 / timt0
! Paul Stackhouse: July 5, 98

       if ( mquadr ) then
         f0(i) = 1.66 * ( 1.0 - w(i) ) * bf(i)
       else
         f0(i) = 2.0 * ( 1.0 - w(i) ) * bf(i)
       endif

       if ( abs(q1) .le. 1.0e-10 ) then
         u0(i) = - 1.0e+10 / q2
       else
         u0(i) = - 1.0 / ( q1 * q2 )
       endif
       t0 = t(i)
3      continue

! +++ PAUSE FIX      call qccgts ( ib, asbs, ee ) 
      call qccgts ( ib, asbs, ee, ITriErr ) 
      IF (ITriErr .NE. 0) RETURN
      xx = aa(1) * expn(1)
      ffu(1) = xx + gamman(1) * bb(1) + caddn0(1)
      ffd(1) = gamman(1) * xx + bb(1) + cminn0(1) 

      do 10 i = 2, m
       k = i - 1
       xx = bb(k) * expn(k)
       ffu(i) = aa(k) + gamman(k) * xx + caddn(k)
       ffd(i) = gamman(k) * aa(k) + xx + cminn(k)
10     continue
      return
      end

! 6-24-98 (8)	
!	subroutine qftisf ( ib, ee )
! +++ PAUSE FIX	subroutine qftisf ( ib, ee, ur )
	subroutine qftisf ( ib, ee, ur, ITriErr  )
! 6-24-98 (8)	
! *********************************************************************
! In this subroutine, the two- and four- stream combination  scheme  or
! the source function technique (Toon et al. 1989) is used to calculate
! the IR radiative fluxes. The exponential approximation for the Planck
! function in optical depth is used ( Fu, 1991).
! At IR wavelengths, the two-stream results are not exact in the limit 
! of no scattering. It also introduces large error in the case of sca-
! ttering. Since the no-scattering limit is of considerable significance
! at IR wavelengths, we have used  the source function technique  that
! would be exact in the limit of the pure absorption and would also en-
! hance the accuracy of the two-stream approach when scattering occurs
! in the IR wavelengths.
! Here, we use nq Gauss points to obtain the fluxes: when nq=2, we use
! double Gaussian quadrature as in Fu and Liou (1993) for  four-stream
! approximation; when nq = 3, we use the regular Gauss quadrature  but
! u1*w1+u2*w2+u3*w3=1.0.
! The upward radiance at the cosine of zenith angle, ur, in the IR
! atmospheric window 800-980 cm**-1, 980-1100 cm**-1, 1100-1250 cm**-1
! is calculated based on equation (2.25) in Fu et al. (1997). 6-24-98.
!
! *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      INTEGER ITriErr
      integer ib,nq,n,m,i,j,i1
      parameter ( nq = 2 )
      real ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t
      real lamdan,gamman,caddn,cminn,caddn0,cminn0,aa,bb,expn
      real g1g2n,fkbn,ffu,ffd,u0,f0,ee,asbs,xx,bf,bs,ugts1
      real t0,q1,q2,x,y,z,y1,yy
      real fg(ndfsx), fh(ndfsx), fj(ndfsx), fk(ndfsx)
      real alfa(ndfsx+1), beta(ndfsx)
      real fiu(mdfsx,nq), fid(mdfsx,nq),ub(ndfsx,nq)
      real fx(ndfsx,nq), fy(ndfsx), fz1(ndfsx,nq), fz2(ndfsx,nq)

      real fuq1(ndfsx), fuq2(ndfsx)
      real ug(nq), wg(nq), ugwg(nq)

      common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),   &
                           ww(ndfsx), tt(ndfsx)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0(ndfsx), f0(ndfsx)
      common / gtscoe / lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx),   &
                          cminn(ndfsx), caddn0(ndfsx), cminn0(ndfsx),   &
                       aa(ndfsx), bb(ndfsx), expn(ndfsx), g1g2n(ndfsx),   &
                          fkbn(ndfsx)
      common /dfsout/ ffu(mdfsx), ffd(mdfsx)
! 6-24-98 (9/10a)	
      real fiurt,fiurw,fiur,ur
	real timt0
      real  fxr(ndfsx), fz1r(ndfsx), fz2r(ndfsx), ubr(ndfsx)
  	common /radiance/ fiurt(nv1x),fiurw(nv1x),fiur(nv1x)  
! 6-24-98 (9/10a)
      common /planci/ bf(nv1x), bs

      data ug / 0.2113248, 0.7886752 /
      data wg / 0.5, 0.5 /
      data ugwg / 0.105662, 0.394338 /
      if ( mquadr ) then
        ugts1 = 0.60241
      else
        ugts1 = 0.5
      endif
      n = ndfs
      m = mdfs
      asbs = bs * ee
      call adjust2 
      t0 = 0.0

      do i = 1, n
       q1 = alog ( bf(i+1) / bf(i) )
! Paul Stackhouse: July 5, 98
!          q2 = 1.0 / ( t(i) - t0 )
           timt0 = t(i) - t0 
           if (timt0 .lt. 1.0e-12) timt0 = 1.0e-12
           q2 = 1.0 / timt0
! Paul Stackhouse: July 5, 98

       if ( mquadr ) then
         f0(i) = 1.66 * ( 1.0 - w(i) ) * bf(i)
       else
         f0(i) = 2.0 * ( 1.0 - w(i) ) * bf(i)
       endif
       if ( abs(q1) .le. 1.0e-10 ) then
         u0(i) = - 1.0e+10 / q2
       else
         u0(i) = - 1.0 / ( q1 * q2 )
       endif
       t0 = t(i)
       beta(i) = - 1.0 / u0(i)
       enddo

! +++ PAUSE FIX      call qccgts ( ib, asbs, ee ) 
      call qccgts ( ib, asbs, ee, ITriErr )
      IF (ITriErr .NE. 0) RETURN 
      do i = 1, n
       x = ( 1.0 - w(i) ) * w(i) / fkbn(i) / ugts1
       y1 = w1(i) / 3.0
       y = g1g2n(i)
       z = -y1 * beta(i)
       fuq1(i) = x * ( y - z ) + 1.0 - w(i)
       fuq2(i) = x * ( y + z ) + 1.0 - w(i)
       enddo

      do i = 1, n + 1
       alfa(i) = 6.2832 * bf(i)
       enddo

      do i = 1, n
       x = lamdan(i) * ugts1
       y = gamman(i) * ( 1.0 + x )
       z = 1.0 - x
       fg(i) = aa(i) * ( z + z )
       fh(i) = bb(i) * ( y + y )
       fj(i) = aa(i) * ( y + y )
       fk(i) = bb(i) * ( z + z )
       enddo

      do j = 1, nq
       fid(1,j) = 0.0
       enddo

      do j = 1, nq
       t0 = 0.0
       do i = 2, mdfs
        i1 = i - 1
        fx(i1,j) = exp ( - ( t(i1) - t0 ) / ug(j) )
        fy(i1) = expn(i1)
        xx = lamdan(i1) * ug(j)
        fz1(i1,j) = ( 1.0 - fx(i1,j) * fy(i1) ) / ( xx + 1.0 )
        fz2(i1,j) = ( fx(i1,j) - fy(i1) ) / ( xx - 1.0 )
        ub(i1,j) = ug(j) * beta(i1)
!!      NAIR Fix on 6Dec2007
!!      Origianl statement as follows
!!---------- 4/2/97 (7)
!!        if (ub(i1,j) .eq. 1.0) ub(i1,j) = 1.001
!!---------- 4/2/97 (7)
!!      Replaced by the following patch
        if (abs(ub(i1,j)) .eq. 1.0) ub(i1,j) = 1.001* &
         (abs(ub(i1,j))/ub(i1,j))
!!      End patch 


        fid(i,j) = fid(i1,j) * fx(i1,j) + fj(i1) * fz1(i1,j) +   &
                      fk(i1) * fz2(i1,j) +    &
                      fuq2(i1) / ( ub(i1,j) + 1.0 ) *   &
                      ( alfa(i) - alfa(i1) * fx(i1,j) )
        t0 = t(i1)
        enddo
       enddo

      yy = 0.0
      do j = 1, nq
       yy = yy + ugwg(j) * fid(mdfs,j) 
       enddo

      xx = yy * ( 1.0 - ee ) * 2.0 + 6.2831854 * ee * bs 
      do j = 1, nq
       fiu(mdfs,j) = xx
       enddo
! 6-24-98 (11)	
	fiur(mdfs) = xx
! 6-24-98 (11)
      do j = 1, nq
       do i = mdfs - 1, 1, -1
               if ( ub(i,j) .eq. 1.0 ) then
                  ub(i,j) = 1.001
               endif

        fiu(i,j) = fiu(i+1,j) * fx(i,j) + fg(i) * fz2(i,j) +   &
                      fh(i) * fz1(i,j) +    &
                      fuq1(i) / ( ub(i,j) - 1.0 ) *   &
                      ( alfa(i+1) * fx(i,j) - alfa(i) )
        enddo
       enddo

      do i = 1, mdfs
       ffu(i) = 0.0
       ffd(i) = 0.0
       enddo

      do i = 1, mdfs
       do j = 1, nq
        ffu(i) = ffu(i) + ugwg(j) * fiu(i,j) 
        ffd(i) = ffd(i) + ugwg(j) * fid(i,j)
        enddo
       enddo
! 6-24-98 (12) 
!!	if ( ib .eq. 11 .or. ib .eq. 12 .or. ib .eq. 13 ) then
	if ( ur .eq. 0.5 ) then
           ur = 0.500001
	endif
	t0 = 0.0
	do i = 2, mdfs
	   i1 = i - 1
	   fxr(i1) = exp ( - ( t(i1) - t0 ) / ur )
	   xx = lamdan(i1) * ur
           fz1r(i1) = ( 1.0 - fxr(i1) * fy(i1) ) / ( xx + 1.0 )
           fz2r(i1) = ( fxr(i1) - fy(i1) ) / ( xx - 1.0 )
	   ubr(i1) = ur * beta(i1)
	   t0 = t(i1)
	enddo
	do i = mdfs - 1, 1, -1
           fiur(i) = fiur(i+1) * fxr(i) + fg(i) * fz2r(i) +   &
                     fh(i) * fz1r(i) +    &
                     fuq1(i) / ( ubr(i) - 1.0 ) *   &
                     ( alfa(i+1) * fxr(i) - alfa(i) ) 
	enddo
!	else
!	do i = 1, mdfs
!           fiur(i) = 0.0
!	enddo
!	endif
! 6-24-98 (12) 
      return
      end

      subroutine adjust2 
!  *****************************************************************
!  In this subroutine, we incorporate a delta-function adjustment to
!  account for the  forward  diffraction  peak in the context of the 
!  two-stream approximation ( Liou, Fu and Ackerman, 1988 ).  w1(n),
!  w(n), and t(n) are the adjusted parameters.
!  *****************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,n
      real ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t,u0a
      real f0a,dtt(ndfsx), dt(ndfsx),tt0,f,fw
      common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),   &
                           ww(ndfsx), tt(ndfsx)
      common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),    &
                            t(ndfsx), u0a(ndfsx), f0a(ndfsx)

      n = ndfs
      tt0 = 0.0
      do 10 i = 1, n
       f = ww2(i) / 5.0
       fw = 1.0 - f * ww(i) 
       w1(i) = ( ww1(i) - 3.0 * f ) / ( 1.0 - f )
       w(i) = ( 1.0 - f ) * ww(i) / fw
       dtt(i) = tt(i) - tt0
       tt0 = tt(i)
       dt(i) = dtt(i) * fw
10     continue
      t(1) = dt(1)

      do 20 i = 2, n
       t(i) = dt(i) + t(i-1)
20     continue

      return
      end

!---------- 4/2/97 (3) -- NEXT 801 LINES
      subroutine gascon_ckd_parm( ib )
!  *******************************************
!          4/2/97 - F.Rose
!  Parameterized CKD_2.1 continuum absorption.
!  *******************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer m,ib
      real pp,pt,ph,po,tgm,dz,dz_quick
      real dp,amnt,patm,temp,parm_ckd
      integer iflb(18)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /con/ tgm(nvx)
      
      data iflb /6*0,12,11,10,9,8,7,6,5,4,3,2,1/

      do m=1,nv
      tgm(m)=0.0
      end do

      if( iflb(ib) .eq. 0) return

      do m=1,nv
       dp     =   pp(m+1) - pp(m)
       amnt   = 1.02*dp*(ph(m)+ph(m+1) ) *0.5
       if(amnt .lt. 0.0)print*,'PROBLEMMMMMM'
       patm   = (( pp(m)   + pp(m+1) ) *0.5 ) /1013.25
       temp   = ( pt(m)   + pt(m+1) ) *0.5
	dz = dz_quick(pp(m),pp(m+1),pt(m),pt(m+1),ph(m),ph(m+1) )
       tgm(m) = parm_ckd(iflb(ib),amnt,patm,temp,dz)
       end do

      return
      end
! --------------------------------------------
	function dz_quick(p1,p2,t1,t2,h1,h2)
	temp   =  ( t1   + t2 ) *0.5
	hum    =  ( h1   + h2 ) *0.5
	tv     = temp*(1+.61*hum)
	dz_quick =29.3*tv*log(p2/p1)*0.001
	return
	end
! --------------------------------------------
!*********************************************************************
      function parm_ckd(iband,amnt,patm,temp,dz)
! Parameterization of CKD_2.1 continuum over Fu-Liou Bands
! Input:
! iband  =  integer (1-8) where
!	  Band 1 ='  5:280cm-1'
!         Band 2 ='280:400cm-1'
!         Band 3 ='400:540cm-1'
!         Band 4 ='540:670cm-1'
!         Band 5 ='670:800cm-1'
!         Band 6 ='800:980cm-1'
!         Band 7 ='980:1100cm-1'
!         Band 8 ='1100:1250cm-1'
!         Band 9 ='1250:1400cm-1'
!         Band10 ='1400:1700cm-1'
!         Band11 ='1700:1900cm-1'
!         Band12 ='1900:2200cm-1'
! amnt = h2O ammount (g/cm**2)
! patm= pressure (atm)
! temp = temperature (k)
! dz = pathlength (Km)
! Output:
! clough_parm = parameterized CKD_2.1optical depth for band
!234567890123456789012345678901234567890123456789012345678901234567890

        parameter (ncoef=7,nband=12)
	real aa(ncoef,nband) ,aa1,aa2,aa3,aa4
	common /pcont_aa0/ aa0(ncoef,nband) ! ORIGINAL 280:1400
	common /pcont_aa1/ aa1(ncoef,nband)! lin no plank 
	common /pcont_aa2/ aa2(ncoef,nband)! log no plank
	common /pcont_aa3/ aa3(ncoef,nband)! lin Plank wgt
	common /pcont_aa4/ aa4(ncoef,nband)! log Plank wgt
	common /cont_tas/ iwtas
        if(iwtas <0 .or. iwtas >4) stop 'iwtas'
	if ( iwtas == 0) aa= aa0
	if ( iwtas == 1) aa= aa1
	if ( iwtas == 2) aa= aa2
	if ( iwtas == 3) aa= aa3
	if ( iwtas == 4) aa= aa4

       ph2o = amnt *(8.314d+07 *temp )/   &
                    (dz*1.0d+05*18.01534 *1.01325d+06)

	if (iwtas == 0) patmx = patm
	if (iwtas >  0) patmx = log(patm)


 	tau_log	    = aa(1,iband)	      +   &
      		      aa(2,iband)* log(amnt)  +   &
      		      aa(3,iband)* temp       +   &
                      aa(4,iband)* patmx       +   &
      		      aa(5,iband)* (ph2o)     +   &
      		      aa(6,iband)* amnt       +   &
      		      aa(7,iband)* log(ph2o) 

    
	parm_ckd = exp ( tau_log )

	return
	end

      subroutine gascon_off
!  ******************************
!          4/2/97 - F.Rose
! Continuum absorption turned off.
!  ******************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i
      real tgm
      common /con/ tgm(nvx)

      do i = 1, nv
       tgm(i) = 0.0
       end do

      return
      end

      subroutine gascon_ckd ( ib )
!  ********************************************************************
!  Exact CKD_2.1 continuum absorption.
!
!  tgm(nv) are the optical depthes due to water vapor continuum absorp-
!  tion in nv layers for a given band ib. We include continuum absorp-
!  tion in the 280 to 1400 cm**-1 region.
! *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      integer i,m,ib
      real pp,pt,ph,po
      real tgm
      integer iwtas
      common /con/ tgm(nvx)
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
	common /cont_tas/ iwtas

      integer iflb(18)
      real dp,p,t,pres,sum,dvabs
      real war,wn2,wo2,wco2,wh2o
      integer ii,iv1,iv2,idv,nptabs
	integer ib1,ie1
	real fe,fr
      real*8  v1abs,v2abs,vbnds(12,2)
      real absrb(1000),wnbr(1000)
	real dz_quick ,dz
       data vbnds/5,  280,400,540,670,800,980,1100,1250,1400,1700,1900,   &
                  280,400,540,670,800,980,1100,1250,1400,1700,1900,2200/

      data dvabs /1.0/
         if(iwtas <0 .or. iwtas >4) stop 'iwtas'

      if (iwtas == 0)iflb =(/0,0,0,0,0,0, 0, 0, 0,9,8,7,6,5,4,3,2,0/)
      if (iwtas > 0) iflb=(/0,0,0,0,0,0, 12,11,10,9,8,7,6,5,4,3,2,1/)
    
      do m=1,nv
      tgm(m)=0.0
      enddo

      if( iflb(ib) .eq. 0) return

      v1abs= vbnds(iflb(ib),1)
      v2abs= vbnds(iflb(ib),2)
      iv1=v1abs
      iv2=v2abs
      idv=dvabs

      ii=0
      do i=iv1,iv2,idv
       ii=ii+1
       wnbr(ii)=i
      enddo
      nptabs=ii

      do m=1,nv
       dp  =   pp(m+1) - pp(m)
       pres= ( pp(m) + pp(m+1) ) *0.5
       t   = ( pt(m) + pt(m+1) ) *0.5
	dz = dz_quick(pp(m),pp(m+1),pt(m),pt(m+1),ph(m),ph(m+1) )
       if(dp .lt. 0.0)print*,'PROBLEMMMMMM'
       wh2o =1.02*dp*(ph(m)+ph(m+1))*0.5   &
           *(6.023e+23/18.01534) ! g/g -> g/cm2 -->molecules/cm2

       p    = pres/1013.25
       war  = p*1.01325e+06/(1.3805e-16*t)*0.00934*1.0e+05*dz
       wn2  = p*1.01325e+06/(1.3805e-16*t)*0.78084*1.0e+05*dz
       wo2  = p*1.01325e+06/(1.3805e-16*t)*0.20948*1.0e+05*dz
       wco2 = p*1.01325e+06/(1.3805e-16*t)*3.50e-04*1.0e+05*dz

      do ii=1,1000
       absrb(ii)=0.0
       enddo

      call clough(pres,t,war,wn2,wo2,wco2,wh2o,   &
                    v1abs,v2abs,dvabs,nptabs,absrb)
 
      ii=0
      sum=0
      do i=iv1,iv2,idv
      ib1=i
      ie1=i+idv
      if ( iwtas == 3 .or.iwtas == 4)    &
       call subplank(t,iv1,iv2,idv,ib1,ie1,fe,fr)

       ii=ii+1
      if ( iwtas == 0)sum=sum+log(absrb(ii))       !! log noplank

      if ( iwtas == 1)sum=sum+    absrb(ii)        !! lin noplank
      if ( iwtas == 2)sum=sum+log(absrb(ii))       !! log noplank
      if ( iwtas == 3)sum = sum +absrb(ii)*fr      !! lin Plank wgt
      if ( iwtas == 4)sum = sum +log(absrb(ii))*fr !! log Plank wgt
       enddo
      if ( iwtas == 0)tgm(m) = exp(sum/float(ii))  !! log noplank

      if ( iwtas == 1)tgm(m) = sum/float(ii)       !! lin noplank
      if ( iwtas == 2)tgm(m) = exp(sum/float(ii))  !! log noplank
      if ( iwtas == 3)tgm(m) = sum                 !! lin Plank wgt
      if ( iwtas == 4)tgm(m) = exp(sum)            !! log Plank wgt
      enddo

      return
      end

!  ************************************************************
!  The Following routines are for the exact CKD_2.1 solution
!  From S.A.Clough with modifications for LWONLY by Dave Kratz.
!  ************************************************************
      subroutine clough(pres,temp,wa,wn2,wo2,wco2,wh2o,v1abs,v2abs,   &
                        dvabs,nptabs,absrb)
      implicit real*8 (v)
      dimension absrb(*)
      dimension wk(60)
      data xlosmt/2.68675e+19/
      data wk/60*0/

!  *********************************************************************
!   this program calculates the continuum optical depth
!         for an homogeneous layer
!
!   the following quantities must be specified:
!
!          pressure                   pave (mb)
!          temperature                tave ( k)
!          column amount
!            nitrogen                 wn2    (molec/cm**2)
!            oxygen                   wk(7)  (molec/cm**2)
!            carbon dioxide           wk(2)  (molec/cm**2)
!            water vapor              wk(1)  (molec/cm**2)
!          number of molecules        nmol
!          beginning wavenumber       v1abs (cm-1)
!          ending wavenumber          v2abs (cm-1)
!          sampling interval          dvabs (cm-1)
!          number of values           nptabs
!
!   the results are in array absorb
!
!   note that for an atmospheric layer:
!
!     wtot   = xlosmt * (pave/1013) * (273/tave) * (path length)
!     wbroad = the column amount for all species not explicitly provided
!     wk(m)  = (volume mixing ratio) * (column of dry air)
!
!      note that continua are included for
!
!             m=1    water vapor    air & self         0 - 20,000 cm-1
!             wn2    nitrogen       air             2020 -   2800 cm-1
!
!  The following is an example for a 1 km path (see cntnout for results)
!  *********************************************************************
      pave = pres
      tave = temp
      wk(7) = wo2
      wk(2) = wco2
      wk(1) =wh2o
      wbroad=wn2+wa
      nmol = 7
      jrad=1

      call contnm(jrad,pave,tave,wk,wbroad,nmol,v1abs,v2abs,dvabs,   &
                  nptabs,absrb)

      return
      end

      subroutine xint (v1a,v2a,dva,a,afact,vft,dvr3,r3,n1r3,n2r3)
!  ****************************************************************
!     this subroutine interpolates the a array stored
!     from v1a to v2a in increments of dva using a multiplicative
!     factor afact, into the r3 array from location n1r3 to n2r3 in
!     increments of dvr3
!  ****************************************************************
      implicit double precision (v)
      dimension a(*),r3(*)
      data onepl/1.001/, onemi/0.999/, argmin/34./

      recdva = 1./dva
      ilo = (v1a+dva-vft)/dvr3+1.+onemi
      ilo = max(ilo,n1r3)
      ihi = (v2a-dva-vft)/dvr3+onemi
      ihi = min(ihi,n2r3)

      do 10 i = ilo, ihi
       vi = vft+dvr3*float(i-1)
       j = (vi-v1a)*recdva+onepl
       vj = v1a+dva*float(j-1)
       p = recdva*(vi-vj)
       c = (3.-2.*p)*p*p
       b = 0.5*p*(1.-p)
       b1 = b*(1.-p)
       b2 = b*p
       conti = -a(j-1)*b1+a(j)*(1.-c+b2)+a(j+1)*(c+b1)-a(j+2)*b2
       r3(i) = r3(i)+conti*afact
 10    continue

      return
      end

      function radfn (vi,xkt)

!
!     function radfn calculates the radiation term for the line shape
!  *********************************************************************
!
!               last modification:    12 august 1991
!
!                  implementation:    r.d. worsham
!
!             algorithm revisions:    s.a. clough
!                                     r.d. worsham
!                                     j.l. moncet
!
!
!                     atmospheric and environmental research inc.
!                     840 memorial drive,  cambridge, ma   02139
!
!----------------------------------------------------------------------
!
!               work supported by:    the arm program
!                                     office of energy research
!                                     department of energy
!
!
!      source of original routine:    afgl line-by-line model
!
!                                             fascod3
!
!
!  *********************************************************************
      implicit double precision (v)
      data onepl/1.001/, onemi/0.999/, argmin/34./

!  *******************************************
!  In the small xviokt region 0.5 is required.
!  *******************************************
      xvi = vi
      if (xkt.gt.0.0) then
        xviokt = xvi/xkt
        if (xviokt.le.0.01) then
          radfn = 0.5*xviokt*xvi
        elseif (xviokt.le.10.0) then
          expvkt = exp(-xviokt)
          radfn = xvi*(1.-expvkt)/(1.+expvkt)
        else
          radfn = xvi
        endif
      else
        radfn = xvi
      endif

      return
      end

      subroutine contnm(jrad,pave,tave,wk,wbroad,nmol,v1abs,v2abs,   &
                        dvabs,nptabs,absrb)
!  *********************************************
!  Subroutine contnm contains the continuum data
!  which is interpolated into the array absrb
!  *********************************************
      implicit double precision (v)
      dimension absrb(*)
      dimension wk(60)
      dimension csh2o(2050),cfh2o(2050)
      dimension c(2050),c0(2050),c1(2050),c2(2050)
      dimension xfac(0:50)

      data pi /3.1415927/
      data planck / 6.626176e-27 /,boltz / 1.380662e-16 /,   &
           clight / 2.99792458e10 /,avog / 6.022045e23 /
      data p0 / 1013. /,t0 / 296. /
      data xlosmt / 2.68675e+19 /

!  ****************************************************************
!  These are self-continuum modification factors from 700-1200 cm-1
!  ****************************************************************

      data (xfac(i),i=0,50)/   &
          1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,   &
          1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,   &
          1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,   &
          1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,   &
          1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,   &
          1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,   &
          1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,   &
          1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,   &
          1.00000,1.00000,1.00000/

      radcn1=2.*planck*clight*clight*1.e-07
      radcn2=planck*clight/boltz

!  ************************************
!  Assign sccs version number to module
!      hvrcnt = '3.3'
!  ************************************
      rhoave = (pave/p0)*(t0/tave)
      xkt = tave/radcn2
      wtot = wbroad
      do 10 m = 1, nmol
       wtot = wtot+wk(m)
 10    continue
      wtot = wtot*1.e-20
      patm = pave/p0

!  ******************************
!  ********    water    ********
!  ******************************
      call sl296 (v1c,v2c,dvc,nptc,c0,v1abs,v2abs,dvabs,nptabs)
      call sl260 (v1c,v2c,dvc,nptc,c1,v1abs,v2abs,dvabs,nptabs)
      call frn296 (v1c,v2c,dvc,nptc,c2,v1abs,v2abs,dvabs,nptabs)

      w1 = wk(1)*1.0e-20
!  *************************************************************
!  The factor of 1.e-20 is handled this way to avoid underflows.
!  *************************************************************
      ph2o = patm*(w1/wtot)
      rh2o = ph2o*(t0/tave)
      pfrgn = patm-ph2o
      rfrgn = pfrgn*(t0/tave)
      xkt = tave/radcn2
      tfac = (tave-t0)/(260.-t0)

!  ******************************
!  ********     self     ********
!  ******************************
      alpha2 = 200.**2
      alphs2= 120.**2
      betas = 5.e-06
      v0s=1310.
      factrs= 0.15

!  ******************************
!  ********    foreign   ********
!  ******************************
      hwsqf= 330.**2
      betaf = 8.0e-11
      v0f =1130.
      factrf = 0.97

      v0f2 =1900.
      hwsqf2 = 150.**2
      beta2 = 3.e-06


      v1h=v1c
      dvh=dvc
      npth=nptc

      do 20 j = 1, nptc
       vj = v1c+dvc*float(j-1)
       vs2 = (vj-v0s)**2
       sh2o = 0.
       if (c0(j).gt.0.) then
         sh2o = c0(j)*(c1(j)/c0(j))**tfac
         sfac = 1.
         if (vj.ge.700. .and.  vj.le.1200.) then
           jfac = (vj-700.)/10. + 0.00001
           sfac = xfac(jfac)
         endif

!     ----------------------------------------------------------------
!     correction to self continuum (1 sept 85); factor of 0.78 at 1000
!                             and  .......
!     ----------------------------------------------------------------
       sh2o = sfac * sh2o*(1.-0.2333*(alpha2/((vj-1050.)**2+alpha2))) *   &
                        (1.-factrs*(alphs2/(vs2+(betas*vs2**2)+alphs2)))
       endif

!     -------------------------------
!     Correction to foreign continuum
!     -------------------------------
       vf2 = (vj-v0f)**2
       vf6 = vf2 * vf2 * vf2
       fscal  = (1.-factrf*(hwsqf/(vf2+(betaf*vf6)+hwsqf)))
       vf2 = (vj-v0f2)**2
       vf4 = vf2*vf2
       fscal = fscal* (1.- 0.6*(hwsqf2/(vf2 + beta2*vf4 + hwsqf2)))

       c2(j)=c2(j)*fscal
       c(j) = w1*(sh2o*rh2o+c2(j)*rfrgn)

       csh2o(j)=1.e-20 * sh2o
       cfh2o(j)=1.e-20 * c2(j)

!     ---------------
!     Radiation field
!     ---------------
       if (jrad.eq.1) c(j) = c(j)*radfn(vj,xkt)
 20    continue

      call xint (v1c,v2c,dvc,c,1.0,v1abs,dvabs,absrb,1,nptabs)

      return
      end


      subroutine sl296 (v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
!  ******************************************
!               06/28/82
!      units of (cm**3/mol) * 1.e-20
!  ******************************************
      implicit double precision (v)
      dimension s(305)
      dimension c(*)
      data v1s,v2s,dvs,npts / -20.0, 3020.0, 10.0, 305/

      data (s(i),i=1,2)/   &
           1.1109e-01 ,1.0573e-01/
      data (s(i),i=3,52)/   &
           1.0162e-01, 1.0573e-01, 1.1109e-01, 1.2574e-01, 1.3499e-01,   &
           1.4327e-01, 1.5065e-01, 1.5164e-01, 1.5022e-01, 1.3677e-01,   &
           1.3115e-01, 1.2253e-01, 1.1271e-01, 1.0070e-01, 8.7495e-02,   &
           8.0118e-02, 6.9940e-02, 6.2034e-02, 5.6051e-02, 4.7663e-02,   &
           4.2450e-02, 3.6690e-02, 3.3441e-02, 3.0711e-02, 2.5205e-02,   &
           2.2113e-02, 1.8880e-02, 1.6653e-02, 1.4626e-02, 1.2065e-02,   &
           1.0709e-02, 9.1783e-03, 7.7274e-03, 6.7302e-03, 5.6164e-03,   &
           4.9089e-03, 4.1497e-03, 3.5823e-03, 3.1124e-03, 2.6414e-03,   &
           2.3167e-03, 2.0156e-03, 1.7829e-03, 1.5666e-03, 1.3928e-03,   &
           1.2338e-03, 1.0932e-03, 9.7939e-04, 8.8241e-04, 7.9173e-04/
      data (s(i),i=53,102)/   &
           7.1296e-04, 6.4179e-04, 5.8031e-04, 5.2647e-04, 4.7762e-04,   &
           4.3349e-04, 3.9355e-04, 3.5887e-04, 3.2723e-04, 2.9919e-04,   &
           2.7363e-04, 2.5013e-04, 2.2876e-04, 2.0924e-04, 1.9193e-04,   &
           1.7618e-04, 1.6188e-04, 1.4891e-04, 1.3717e-04, 1.2647e-04,   &
           1.1671e-04, 1.0786e-04, 9.9785e-05, 9.2350e-05, 8.5539e-05,   &
           7.9377e-05, 7.3781e-05, 6.8677e-05, 6.3993e-05, 5.9705e-05,   &
           5.5788e-05, 5.2196e-05, 4.8899e-05, 4.5865e-05, 4.3079e-05,   &
           4.0526e-05, 3.8182e-05, 3.6025e-05, 3.4038e-05, 3.2203e-05,   &
           3.0511e-05, 2.8949e-05, 2.7505e-05, 2.6170e-05, 2.4933e-05,   &
           2.3786e-05, 2.2722e-05, 2.1736e-05, 2.0819e-05, 1.9968e-05/
      data (s(i),i=103,152)/   &
           1.9178e-05, 1.8442e-05, 1.7760e-05, 1.7127e-05, 1.6541e-05,   &
           1.5997e-05, 1.5495e-05, 1.5034e-05, 1.4614e-05, 1.4230e-05,   &
           1.3883e-05, 1.3578e-05, 1.3304e-05, 1.3069e-05, 1.2876e-05,   &
           1.2732e-05, 1.2626e-05, 1.2556e-05, 1.2544e-05, 1.2604e-05,   &
           1.2719e-05, 1.2883e-05, 1.3164e-05, 1.3581e-05, 1.4187e-05,   &
           1.4866e-05, 1.5669e-05, 1.6717e-05, 1.8148e-05, 2.0268e-05,   &
           2.2456e-05, 2.5582e-05, 2.9183e-05, 3.3612e-05, 3.9996e-05,   &
           4.6829e-05, 5.5055e-05, 6.5897e-05, 7.5360e-05, 8.7213e-05,   &
           1.0046e-04, 1.1496e-04, 1.2943e-04, 1.5049e-04, 1.6973e-04,   &
           1.8711e-04, 2.0286e-04, 2.2823e-04, 2.6780e-04, 2.8766e-04/
      data (s(i),i=153,202)/   &
           3.1164e-04, 3.3640e-04, 3.6884e-04, 3.9159e-04, 3.8712e-04,   &
           3.7433e-04, 3.4503e-04, 3.1003e-04, 2.8027e-04, 2.5253e-04,   &
           2.3408e-04, 2.2836e-04, 2.4442e-04, 2.7521e-04, 2.9048e-04,   &
           3.0489e-04, 3.2646e-04, 3.3880e-04, 3.3492e-04, 3.0987e-04,   &
           2.9482e-04, 2.8711e-04, 2.6068e-04, 2.2683e-04, 1.9996e-04,   &
           1.7788e-04, 1.6101e-04, 1.3911e-04, 1.2013e-04, 1.0544e-04,   &
           9.4224e-05, 8.1256e-05, 7.3667e-05, 6.2233e-05, 5.5906e-05,   &
           5.1619e-05, 4.5140e-05, 4.0273e-05, 3.3268e-05, 3.0258e-05,   &
           2.6440e-05, 2.3103e-05, 2.0749e-05, 1.8258e-05, 1.6459e-05,   &
           1.4097e-05, 1.2052e-05, 1.0759e-05, 9.1400e-06, 8.1432e-06/
      data (s(i),i=203,252)/   &
           7.1460e-06, 6.4006e-06, 5.6995e-06, 4.9372e-06, 4.4455e-06,   &
           3.9033e-06, 3.4740e-06, 3.1269e-06, 2.8059e-06, 2.5558e-06,   &
           2.2919e-06, 2.0846e-06, 1.8983e-06, 1.7329e-06, 1.5929e-06,   &
           1.4631e-06, 1.3513e-06, 1.2461e-06, 1.1519e-06, 1.0682e-06,   &
           9.9256e-07, 9.2505e-07, 8.6367e-07, 8.0857e-07, 7.5674e-07,   &
           7.0934e-07, 6.6580e-07, 6.2580e-07, 5.8853e-07, 5.5333e-07,   &
           5.2143e-07, 4.9169e-07, 4.6431e-07, 4.3898e-07, 4.1564e-07,   &
           3.9405e-07, 3.7403e-07, 3.5544e-07, 3.3819e-07, 3.2212e-07,   &
           3.0714e-07, 2.9313e-07, 2.8003e-07, 2.6777e-07, 2.5628e-07,   &
           2.4551e-07, 2.3540e-07, 2.2591e-07, 2.1701e-07, 2.0866e-07/
      data (s(i),i=253,302)/   &
           2.0082e-07, 1.9349e-07, 1.8665e-07, 1.8027e-07, 1.7439e-07,   &
           1.6894e-07, 1.6400e-07, 1.5953e-07, 1.5557e-07, 1.5195e-07,   &
           1.4888e-07, 1.4603e-07, 1.4337e-07, 1.4093e-07, 1.3828e-07,   &
           1.3569e-07, 1.3270e-07, 1.2984e-07, 1.2714e-07, 1.2541e-07,   &
           1.2399e-07, 1.2102e-07, 1.1878e-07, 1.1728e-07, 1.1644e-07,   &
           1.1491e-07, 1.1305e-07, 1.1235e-07, 1.1228e-07, 1.1224e-07,   &
           1.1191e-07, 1.1151e-07, 1.1098e-07, 1.1068e-07, 1.1109e-07,   &
           1.1213e-07, 1.1431e-07, 1.1826e-07, 1.2322e-07, 1.3025e-07,   &
           1.4066e-07, 1.5657e-07, 1.7214e-07, 1.9449e-07, 2.2662e-07,   &
           2.6953e-07, 3.1723e-07, 3.7028e-07, 4.4482e-07, 5.3852e-07/
      data (s(i),i=303,305)/   &
           6.2639e-07, 7.2175e-07, 7.7626e-07/

      dvc = dvs
      v1c = v1abs-dvc
      v2c = v2abs+dvc

      i1 = (v1c-v1s)/dvs
      if (v1c.lt.v1s) i1 = i1-1

      v1c = v1s+dvs*float(i1)
      i2 = (v2c-v1s)/dvs
      nptc = i2-i1+3
      v2c = v1c+dvs*float(nptc-1)
      do 10 j = 1, nptc
       i = i1+j
       c(j) = 0.
       if ((i.lt.1).or.(i.gt.npts)) go to 10
       c(j) = s(i)
 10    continue

      return
      end

      subroutine sl260 (v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
!  **************************************
!               06/28/82
!      units of (cm**3/mol) * 1.e-20
!  **************************************
      implicit double precision (v)
      dimension s(305)
      dimension c(*)

      data v1s,v2s,dvs,npts / -20.0, 3020.0, 10.0, 305/

      data (s(i),i=1,2)/   &
           1.7750e-01, 1.7045e-01/
      data (s(i),i=3,52)/   &
           1.6457e-01, 1.7045e-01, 1.7750e-01, 2.0036e-01, 2.1347e-01,   &
           2.2454e-01, 2.3428e-01, 2.3399e-01, 2.3022e-01, 2.0724e-01,   &
           1.9712e-01, 1.8317e-01, 1.6724e-01, 1.4780e-01, 1.2757e-01,   &
           1.1626e-01, 1.0098e-01, 8.9033e-02, 7.9770e-02, 6.7416e-02,   &
           5.9588e-02, 5.1117e-02, 4.6218e-02, 4.2179e-02, 3.4372e-02,   &
           2.9863e-02, 2.5252e-02, 2.2075e-02, 1.9209e-02, 1.5816e-02,   &
           1.3932e-02, 1.1943e-02, 1.0079e-02, 8.7667e-03, 7.4094e-03,   &
           6.4967e-03, 5.5711e-03, 4.8444e-03, 4.2552e-03, 3.6953e-03,   &
           3.2824e-03, 2.9124e-03, 2.6102e-03, 2.3370e-03, 2.1100e-03,   &
           1.9008e-03, 1.7145e-03, 1.5573e-03, 1.4206e-03, 1.2931e-03/
      data (s(i),i=53,102)/   &
           1.1803e-03, 1.0774e-03, 9.8616e-04, 9.0496e-04, 8.3071e-04,   &
           7.6319e-04, 7.0149e-04, 6.4637e-04, 5.9566e-04, 5.4987e-04,   &
           5.0768e-04, 4.6880e-04, 4.3317e-04, 4.0037e-04, 3.7064e-04,   &
           3.4325e-04, 3.1809e-04, 2.9501e-04, 2.7382e-04, 2.5430e-04,   &
           2.3630e-04, 2.1977e-04, 2.0452e-04, 1.9042e-04, 1.7740e-04,   &
           1.6544e-04, 1.5442e-04, 1.4425e-04, 1.3486e-04, 1.2618e-04,   &
           1.1817e-04, 1.1076e-04, 1.0391e-04, 9.7563e-05, 9.1696e-05,   &
           8.6272e-05, 8.1253e-05, 7.6607e-05, 7.2302e-05, 6.8311e-05,   &
           6.4613e-05, 6.1183e-05, 5.8001e-05, 5.5048e-05, 5.2307e-05,   &
           4.9761e-05, 4.7395e-05, 4.5197e-05, 4.3155e-05, 4.1256e-05/
      data (s(i),i=103,152)/   &
           3.9491e-05, 3.7849e-05, 3.6324e-05, 3.4908e-05, 3.3594e-05,   &
           3.2374e-05, 3.1244e-05, 3.0201e-05, 2.9240e-05, 2.8356e-05,   &
           2.7547e-05, 2.6814e-05, 2.6147e-05, 2.5551e-05, 2.5029e-05,   &
           2.4582e-05, 2.4203e-05, 2.3891e-05, 2.3663e-05, 2.3531e-05,   &
           2.3483e-05, 2.3516e-05, 2.3694e-05, 2.4032e-05, 2.4579e-05,   &
           2.5234e-05, 2.6032e-05, 2.7119e-05, 2.8631e-05, 3.0848e-05,   &
           3.3262e-05, 3.6635e-05, 4.0732e-05, 4.5923e-05, 5.3373e-05,   &
           6.1875e-05, 7.2031e-05, 8.5980e-05, 9.8642e-05, 1.1469e-04,   &
           1.3327e-04, 1.5390e-04, 1.7513e-04, 2.0665e-04, 2.3609e-04,   &
           2.6220e-04, 2.8677e-04, 3.2590e-04, 3.8624e-04, 4.1570e-04/
      data (s(i),i=153,202)/   &
           4.5207e-04, 4.9336e-04, 5.4500e-04, 5.8258e-04, 5.8086e-04,   &
           5.6977e-04, 5.3085e-04, 4.8020e-04, 4.3915e-04, 4.0343e-04,   &
           3.7853e-04, 3.7025e-04, 3.9637e-04, 4.4675e-04, 4.7072e-04,   &
           4.9022e-04, 5.2076e-04, 5.3676e-04, 5.2755e-04, 4.8244e-04,   &
           4.5473e-04, 4.3952e-04, 3.9614e-04, 3.4086e-04, 2.9733e-04,   &
           2.6367e-04, 2.3767e-04, 2.0427e-04, 1.7595e-04, 1.5493e-04,   &
           1.3851e-04, 1.1874e-04, 1.0735e-04, 9.0490e-05, 8.1149e-05,   &
           7.4788e-05, 6.5438e-05, 5.8248e-05, 4.8076e-05, 4.3488e-05,   &
           3.7856e-05, 3.3034e-05, 2.9592e-05, 2.6088e-05, 2.3497e-05,   &
           2.0279e-05, 1.7526e-05, 1.5714e-05, 1.3553e-05, 1.2145e-05/
      data (s(i),i=203,252)/   &
           1.0802e-05, 9.7681e-06, 8.8196e-06, 7.8291e-06, 7.1335e-06,   &
           6.4234e-06, 5.8391e-06, 5.3532e-06, 4.9079e-06, 4.5378e-06,   &
           4.1716e-06, 3.8649e-06, 3.5893e-06, 3.3406e-06, 3.1199e-06,   &
           2.9172e-06, 2.7348e-06, 2.5644e-06, 2.4086e-06, 2.2664e-06,   &
           2.1359e-06, 2.0159e-06, 1.9051e-06, 1.8031e-06, 1.7074e-06,   &
           1.6185e-06, 1.5356e-06, 1.4584e-06, 1.3861e-06, 1.3179e-06,   &
           1.2545e-06, 1.1951e-06, 1.1395e-06, 1.0873e-06, 1.0384e-06,   &
           9.9250e-07, 9.4935e-07, 9.0873e-07, 8.7050e-07, 8.3446e-07,   &
           8.0046e-07, 7.6834e-07, 7.3800e-07, 7.0931e-07, 6.8217e-07,   &
           6.5648e-07, 6.3214e-07, 6.0909e-07, 5.8725e-07, 5.6655e-07/
      data (s(i),i=253,302)/   &
           5.4693e-07, 5.2835e-07, 5.1077e-07, 4.9416e-07, 4.7853e-07,   &
           4.6381e-07, 4.5007e-07, 4.3728e-07, 4.2550e-07, 4.1450e-07,   &
           4.0459e-07, 3.9532e-07, 3.8662e-07, 3.7855e-07, 3.7041e-07,   &
           3.6254e-07, 3.5420e-07, 3.4617e-07, 3.3838e-07, 3.3212e-07,   &
           3.2655e-07, 3.1865e-07, 3.1203e-07, 3.0670e-07, 3.0252e-07,   &
           2.9749e-07, 2.9184e-07, 2.8795e-07, 2.8501e-07, 2.8202e-07,   &
           2.7856e-07, 2.7509e-07, 2.7152e-07, 2.6844e-07, 2.6642e-07,   &
           2.6548e-07, 2.6617e-07, 2.6916e-07, 2.7372e-07, 2.8094e-07,   &
           2.9236e-07, 3.1035e-07, 3.2854e-07, 3.5481e-07, 3.9377e-07,   &
           4.4692e-07, 5.0761e-07, 5.7715e-07, 6.7725e-07, 8.0668e-07/
      data (s(i),i=303,305)/   &
           9.3716e-07, 1.0797e-06, 1.1689e-06/

      dvc = dvs
      v1c = v1abs-dvc
      v2c = v2abs+dvc

      i1 = (v1c-v1s)/dvs
      if (v1c.lt.v1s) i1 = i1-1

      v1c = v1s+dvs*float(i1)
      i2 = (v2c-v1s)/dvs
      nptc = i2-i1+3
      v2c = v1c+dvs*float(nptc-1)
      do 10 j = 1, nptc
       i = i1+j
       c(j) = 0.
       if ((i.lt.1).or.(i.gt.npts)) go to 10
       c(j) = s(i)
 10    continue

      return
      end

      subroutine frn296 (v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
!  **********************************
!               06/28/82
!      units of (cm**3/mol)*1.e-20
!  **********************************
      implicit double precision (v)
      dimension s(305)
      dimension c(*)

      data v1s,v2s,dvs,npts / -20.0, 3020.0, 10.0, 305/

      data (s(i),i=1,2)/   &
           1.2859e-02, 1.1715e-02/
      data (s(i),i=3,52)/   &
           1.1038e-02, 1.1715e-02, 1.2859e-02, 1.5326e-02, 1.6999e-02,   &
           1.8321e-02, 1.9402e-02, 1.9570e-02, 1.9432e-02, 1.7572e-02,   &
           1.6760e-02, 1.5480e-02, 1.3984e-02, 1.2266e-02, 1.0467e-02,   &
           9.4526e-03, 8.0485e-03, 6.9484e-03, 6.1416e-03, 5.0941e-03,   &
           4.4836e-03, 3.8133e-03, 3.4608e-03, 3.1487e-03, 2.4555e-03,   &
           2.0977e-03, 1.7266e-03, 1.4920e-03, 1.2709e-03, 9.8081e-04,   &
           8.5063e-04, 6.8822e-04, 5.3809e-04, 4.4679e-04, 3.3774e-04,   &
           2.7979e-04, 2.1047e-04, 1.6511e-04, 1.2993e-04, 9.3033e-05,   &
           7.4360e-05, 5.6428e-05, 4.5442e-05, 3.4575e-05, 2.7903e-05,   &
           2.1374e-05, 1.6075e-05, 1.3022e-05, 1.0962e-05, 8.5959e-06/
      data (s(i),i=53,102)/   &
           6.9125e-06, 5.3808e-06, 4.3586e-06, 3.6394e-06, 2.9552e-06,   &
           2.3547e-06, 1.8463e-06, 1.6036e-06, 1.3483e-06, 1.1968e-06,   &
           1.0333e-06, 8.4484e-07, 6.7195e-07, 5.0947e-07, 4.2343e-07,   &
           3.4453e-07, 2.7830e-07, 2.3063e-07, 1.9951e-07, 1.7087e-07,   &
           1.4393e-07, 1.2575e-07, 1.0750e-07, 8.2325e-08, 5.7524e-08,   &
           4.4482e-08, 3.8106e-08, 3.4315e-08, 2.9422e-08, 2.5069e-08,   &
           2.2402e-08, 1.9349e-08, 1.6152e-08, 1.2208e-08, 8.9660e-09,   &
           7.1322e-09, 6.1028e-09, 5.2938e-09, 4.5350e-09, 3.4977e-09,   &
           2.9511e-09, 2.4734e-09, 2.0508e-09, 1.8507e-09, 1.6373e-09,   &
           1.5171e-09, 1.3071e-09, 1.2462e-09, 1.2148e-09, 1.2590e-09/
      data (s(i),i=103,152)/   &
           1.3153e-09, 1.3301e-09, 1.4483e-09, 1.6944e-09, 2.0559e-09,   &
           2.2954e-09, 2.6221e-09, 3.2606e-09, 4.2392e-09, 5.2171e-09,   &
           6.2553e-09, 8.2548e-09, 9.5842e-09, 1.1280e-08, 1.3628e-08,   &
           1.7635e-08, 2.1576e-08, 2.4835e-08, 3.0014e-08, 3.8485e-08,   &
           4.7440e-08, 5.5202e-08, 7.0897e-08, 9.6578e-08, 1.3976e-07,   &
           1.8391e-07, 2.3207e-07, 2.9960e-07, 4.0408e-07, 5.9260e-07,   &
           7.8487e-07, 1.0947e-06, 1.4676e-06, 1.9325e-06, 2.6587e-06,   &
           3.4534e-06, 4.4376e-06, 5.8061e-06, 7.0141e-06, 8.4937e-06,   &
           1.0186e-05, 1.2034e-05, 1.3837e-05, 1.6595e-05, 1.9259e-05,   &
           2.1620e-05, 2.3681e-05, 2.7064e-05, 3.2510e-05, 3.5460e-05/
      data (s(i),i=153,202)/   &
           3.9109e-05, 4.2891e-05, 4.7757e-05, 5.0981e-05, 5.0527e-05,   &
           4.8618e-05, 4.4001e-05, 3.7982e-05, 3.2667e-05, 2.7794e-05,   &
           2.4910e-05, 2.4375e-05, 2.7316e-05, 3.2579e-05, 3.5499e-05,   &
           3.8010e-05, 4.1353e-05, 4.3323e-05, 4.3004e-05, 3.9790e-05,   &
           3.7718e-05, 3.6360e-05, 3.2386e-05, 2.7409e-05, 2.3626e-05,   &
           2.0631e-05, 1.8371e-05, 1.5445e-05, 1.2989e-05, 1.1098e-05,   &
           9.6552e-06, 8.0649e-06, 7.2365e-06, 5.9137e-06, 5.2759e-06,   &
           4.8860e-06, 4.1321e-06, 3.5918e-06, 2.7640e-06, 2.4892e-06,   &
           2.1018e-06, 1.7848e-06, 1.5855e-06, 1.3569e-06, 1.1986e-06,   &
           9.4693e-07, 7.4097e-07, 6.3443e-07, 4.8131e-07, 4.0942e-07/
      data (s(i),i=203,252)/   &
           3.3316e-07, 2.8488e-07, 2.3461e-07, 1.7397e-07, 1.4684e-07,   &
           1.0953e-07, 8.5396e-08, 6.9261e-08, 5.4001e-08, 4.5430e-08,   &
           3.2791e-08, 2.5995e-08, 2.0225e-08, 1.5710e-08, 1.3027e-08,   &
           1.0229e-08, 8.5277e-09, 6.5249e-09, 5.0117e-09, 3.9906e-09,   &
           3.2332e-09, 2.7847e-09, 2.4570e-09, 2.3359e-09, 2.0599e-09,   &
           1.8436e-09, 1.6559e-09, 1.4910e-09, 1.2794e-09, 9.8229e-10,   &
           8.0054e-10, 6.0769e-10, 4.5646e-10, 3.3111e-10, 2.4428e-10,   &
           1.8007e-10, 1.3291e-10, 9.7974e-11, 7.8271e-11, 6.3833e-11,   &
           5.4425e-11, 4.6471e-11, 4.0209e-11, 3.5227e-11, 3.1212e-11,   &
           2.8840e-11, 2.7762e-11, 2.7935e-11, 3.2012e-11, 3.9525e-11/
      data (s(i),i=253,302)/   &
           5.0303e-11, 6.8027e-11, 9.3954e-11, 1.2986e-10, 1.8478e-10,   &
           2.5331e-10, 3.4827e-10, 4.6968e-10, 6.2380e-10, 7.9106e-10,   &
           1.0026e-09, 1.2102e-09, 1.4146e-09, 1.6154e-09, 1.7510e-09,   &
           1.8575e-09, 1.8742e-09, 1.8700e-09, 1.8582e-09, 1.9657e-09,   &
           2.1204e-09, 2.0381e-09, 2.0122e-09, 2.0436e-09, 2.1213e-09,   &
           2.0742e-09, 1.9870e-09, 2.0465e-09, 2.1556e-09, 2.2222e-09,   &
           2.1977e-09, 2.1047e-09, 1.9334e-09, 1.7357e-09, 1.5754e-09,   &
           1.4398e-09, 1.4018e-09, 1.5459e-09, 1.7576e-09, 2.1645e-09,   &
           2.9480e-09, 4.4439e-09, 5.8341e-09, 8.0757e-09, 1.1658e-08,   &
           1.6793e-08, 2.2694e-08, 2.9468e-08, 3.9278e-08, 5.2145e-08/
      data (s(i),i=303,305)/   &
           6.4378e-08, 7.7947e-08, 8.5321e-08/

      dvc = dvs
      v1c = v1abs-dvc
      v2c = v2abs+dvc

      i1 = (v1c-v1s)/dvs
      if (v1c.lt.v1s) i1 = i1-1

      v1c = v1s+dvs*float(i1)
      i2 = (v2c-v1s)/dvs
      nptc = i2-i1+3
      v2c = v1c+dvs*float(nptc-1)

      do 10 j = 1, nptc
       i = i1+j
       c(j) = 0.
       if ((i.ge.1).and.(i.le.npts)) then
         c(j) = s(i)
       endif
 10    continue

      return
      end
!---------- 4/2/97 (3) -- PREVIOUS 801 LINES
!==================================================
	subroutine subplank(temp,ib,ie,id,ib1,ie1,fe,fr)
!	data id_fine /1/
!        integer ib0,ie0,id_fine,id0,nnn
!        real fluxb,temp0
!	save ib0,ie0,id0,id_fine,nnn,fluxb,temp0
        integer,save:: ib0,ie0,id_fine,id0,nnn
        real, save:: fluxb,temp0 
        id_fine = 1
  
	if (ib0.ne.ib .or. ie0.ne.ie .or.    &
            id0.ne.id .or.temp.ne.temp0)then
	Print*,'initalize subplank'
	 ib0=ib
	 ie0=ie
	 id0=id
	 temp0=temp
	 call plank(temp,ib,ie,id_fine,fluxb)
	 nnn= (ie-ib)/float(id)  
	endif

	call plank(temp,ib1,ie1,id_fine,flux)
	fr=(flux/fluxb)
	fe=fr*nnn
!	print*,nnn,flux,fluxb,fe

	return
	end
!==================================================
	subroutine plank(t,ib,ie,id,fluxo)
	real*8 c1,c2,u,Flux,el,h,c,rk
	real*4 fluxo
	c=2.99793E+10
	h=6.62620E-27
	rk=1.38062E-16
	C1=3.74E-16
	C2=1.44E-02

	c1=3.741951e-16
	C2=1.439E-02

	rw=1.0e-06

	Flux=0
	do 1 ic=ib,ie-id,id
	um1=1.0e+04/float(ic)
	um2=1.0e+04/float(ic+id)
	dum = um1-um2

	rcm=ic+id/2.0 !center

	um=1.0e+04/rcm !!float(ic) ! WAVENUMBER to Micron

	upw=um/rcm     !!float(ic)

	u=um*1.0E-06
	el=  c1/( u**5 *(exp( c2/(u*t) ) -1.0 ) )

!	rw=rwm*1.0e-06
	rw= dum*1.0e-6

	el2=el*rw

	Flux=Flux+el2
	umo=um
1	continue

	fluxo=flux
	return
	end
! Fu 07-08-98	
	block data ice1new
! *********************************************************************
! Following Fu (1996; J. Climate) and Fu et al. (1998; J. Climate),
! ap is the empirical coefficients of Eq. (3.9a) of Fu (1996) and
! Eq. (3.1) of Fu et al. (1998) to calculate the extiction coefficient
! (1/m).  bps is for the single scattering albedo in the solar bands
! (3.9b in Fu) and bpir is for the absorption coefficient (1/m) in the
! IR bands (3.2 in Fu et al.).  cp is the empirical coefficients of
! Eq. (3.9c) in Fu or Eq. (3.3) in Fu et al. to compute the asymmetry
! factor of the phase function.  dps is the empirical coefficients of
! Eq. (3.9d) of Fu to calculate the forward delta-fraction in the 
! solar bands.  The units of generalized effective size and ice water
! content are um and g/m**3, respectively, in these equations.
! *********************************************************************
!	include 'para.file'
      USE RadParams
!##      include 'rad_0698.h'
	common /ic1new/ ap(3,mbx), bps(4,mbsx), bpir(4, mbirx),   &
                        cp(4,mbx), dps(4,mbsx)
        data ap /    &
                 -2.9172062e-05,  2.5192544e+00,  0.0,   &
                 -2.2948980e-05,  2.5212550e+00,  0.0,   &
                 -2.9772840e-04,  2.5400320e+00,  0.0,   &
                  4.2668223e-04,  2.4933372e+00,  0.0,   &
                  4.3226531e-04,  2.4642946e+00,  0.0,   &
                  9.5918990e-05,  2.5232218e+00,  0.0,   &
      		 -2.308881e-03, 2.814002e+00, 1.072211e+00,   &
      		 -2.465236e-03, 2.833187e+00,-4.227573e-01,   &
      		 -3.034573e-03, 2.900043e+00,-1.849911e+00,   &
      		 -4.936610e-03, 3.087764e+00,-3.884262e+00,   &
      		 -8.178608e-03, 3.401245e+00,-8.812820e+00,   &
      		 -8.372696e-03, 3.455018e+00,-1.516692e+01,   &
      		 -1.691632e-03, 2.765756e+00,-8.331033e+00,   &
      		 -4.159424e-03, 3.047325e+00,-5.061568e+00,   &
      		 -9.524174e-03, 3.587742e+00,-1.068895e+01,   &
      		 -1.334860e-02, 4.043808e+00,-2.171029e+01,   &
      		  3.325756e-03, 2.601360e+00,-1.909602e+01,   &
      		  4.919685e-03, 2.327741e+00,-1.390858e+01 /
	data bps /    &
        1.3540265e-07,  9.9282217e-08, -7.3843168e-11,  3.3111862e-13,   &
       -2.1458450e-06,  2.1984010e-05, -4.4225520e-09,  1.0711940e-11,   &
        1.4027890e-04,  1.3919010e-03, -5.1005610e-06,  1.4032930e-08,   &
        5.7801650e-03,  2.4420420e-03, -1.1985030e-05,  3.3878720e-08,   &
        2.7122737e-01,  1.9809794e-03, -1.5071269e-05,  5.0103900e-08,   &
        1.6215025e-01,  6.3734393e-03, -5.7740959e-05,  1.9109300e-07 /
	data bpir /   &
      	4.346482e-01, 1.721457e-02,-1.623227e-04, 5.561523e-07,   &
        7.428957e-01, 1.279601e-02,-1.391803e-04, 5.180104e-07,   &
        8.862434e-01, 1.226538e-02,-1.523076e-04, 6.000892e-07,   &
        7.152274e-01, 1.621734e-02,-1.868544e-04, 7.078738e-07,   &
        5.874323e-01, 1.876628e-02,-2.045834e-04, 7.510080e-07,   &
        5.409536e-01, 1.949649e-02,-2.050908e-04, 7.364680e-07,   &
        1.195515e+00, 3.350616e-03,-5.266996e-05, 2.233377e-07,   &
        1.466481e+00,-2.129226e-03,-1.361630e-05, 1.193649e-07,   &
        9.551440e-01, 1.309792e-02,-1.793694e-04, 7.313392e-07,   &
        3.003701e-01, 2.051529e-02,-1.931684e-04, 6.583031e-07,   &
        2.005578e-01, 2.132614e-02,-1.751052e-04, 5.355885e-07,   &
        8.869787e-01, 2.118409e-02,-2.781429e-04, 1.094562e-06 /
	data cp /    &
        7.4812728e-01,  9.5684492e-04, -1.1151708e-06, -8.1557303e-09,    &
        7.5212480e-01,  1.1045100e-03, -2.9157100e-06, -1.3429900e-09,   &
        7.5320460e-01,  1.8845180e-03, -9.7571460e-06,  2.2428270e-08,   &
        7.7381780e-01,  2.2260760e-03, -1.4052790e-05,  3.7896870e-08,   &
        8.7020490e-01,  1.6645530e-03, -1.4886030e-05,  4.9867270e-08,   &
        7.4212060e-01,  5.2621900e-03, -5.0877550e-05,  1.7307870e-07,   &
        7.962716e-01, 3.003488e-03,-2.082376e-05, 5.366545e-08,   &
        8.472918e-01, 2.559953e-03,-2.182660e-05, 6.879977e-08,   &
        8.741665e-01, 2.455409e-03,-2.456935e-05, 8.641223e-08,   &
        8.522816e-01, 2.523627e-03,-2.149196e-05, 6.685067e-08,   &
        8.609604e-01, 2.200445e-03,-1.748105e-05, 5.176616e-08,   &
        8.906280e-01, 1.903269e-03,-1.733552e-05, 5.855071e-08,   &
        8.663385e-01, 2.797934e-03,-3.187011e-05, 1.217209e-07,   &
        7.984021e-01, 3.977117e-03,-4.471984e-05, 1.694919e-07,   &
        7.363466e-01, 4.798266e-03,-4.513292e-05, 1.525774e-07,   &
        7.260484e-01, 2.664334e-03,-1.251136e-05, 2.243377e-08,   &
        6.891414e-01, 6.192281e-03,-6.459514e-05, 2.436963e-07,   &
        4.949276e-01, 1.186174e-02,-1.267629e-04, 4.603574e-07 /
	data dps /    &
        1.1572963e-01,  2.5648064e-04,  1.9131293e-06, -1.2460341e-08,   &
        1.1360752e-01,  2.4156171e-04,  2.0185942e-06, -1.2876106e-08,   &
        1.1241170e-01, -1.7635186e-07,  2.1499248e-06, -1.2949304e-08,   &
        1.0855775e-01, -3.2496217e-04,  3.4207304e-06, -1.6247759e-08,   &
        5.7783360e-02, -4.1158260e-04,  4.2361240e-06, -1.7204950e-08,   &
        1.1367129e-01, -1.9711061e-03,  1.6078010e-05, -5.1736898e-08 /
! **********************************************************************
	end

	subroutine icenew ( ib )
! *********************************************************************
! ti, wi, and wwi are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and
! 4) due to the scattering of ice clouds for a given layer.
! *********************************************************************
!	include 'para.file'
      USE RadParams
!##        include 'rad_0698.h'
!c	common /clouds/ pre(nv) , plwc(nv) , pdge(nv), piwc(nv)
        common /clouds/ pre(nvx), plwc(nvx), pdge(nvx), piwc(nvx)
	common /ic1new/ ap(3,mbx), bps(4,mbsx), bpir(4,mbirx),    &
                        cp(4,mbx), dps(4,mbsx)
	common /thick/ dz(nvx)
	common /ic/ ti(nvx), wi(nvx), wwi(nvx,4)
	do 10 i = 1, nv
	   if ( piwc(i) .lt. 1.0e-5 ) then
	     ti(i) = 0.0
	     wi(i) = 0.0
	     wwi(i,1) = 0.0
	     wwi(i,2) = 0.0
	     wwi(i,3) = 0.0
	     wwi(i,4) = 0.0
           else
	     fw1 = pdge(i)
	     fw2 = fw1 * pdge(i)
	     fw3 = fw2 * pdge(i)
	     if ( ib .le. mbs ) then
	        tau = dz(i) * 1000.0 * piwc(i) * ( ap(1,ib) +   &
      	             ap(2,ib) / fw1 )
	        omega = 1.0 - ( bps(1,ib) + bps(2,ib) * fw1 +   &
      	                bps(3,ib) * fw2 + bps(4,ib) * fw3 )
	        asy = cp(1,ib) + cp(2,ib) * fw1 +   &
                      cp(3,ib) * fw2 + cp(4,ib) * fw3
	        fd = dps(1,ib) + dps(2,ib) * fw1 +   &
                     dps(3,ib) * fw2 + dps(4,ib) * fw3
	        f = 0.5 / omega + fd
	        fw = f * omega
                ti(i) = ( 1.0 - fw ) * tau
                wi(i) = ( 1.0 - f ) * omega / ( 1.0 - fw )
                gg = ( asy - f ) / ( 1.0 - f )
	        x1 = gg
                x2 = x1 * gg
                x3 = x2 * gg
                x4 = x3 * gg
                wwi(i,1) = 3.0 * x1
	        wwi(i,2) = 5.0 * x2
                wwi(i,3) = 7.0 * x3
                wwi(i,4) = 9.0 * x4
             else
	        ibr = ib - mbs
	        betae = piwc(i) * ( ap(1,ib) +   &
      	                ap(2,ib) / fw1 + ap(3,ib) / fw2 )
	        betaa = piwc(i) / fw1 * ( bpir(1,ibr) + bpir(2,ibr) *    &
                        fw1 + bpir(3,ibr) * fw2 + bpir(4,ibr) * fw3 )
	        asy = cp(1,ib) + cp(2,ib) * fw1 +   &
                      cp(3,ib) * fw2 + cp(4,ib) * fw3
                ti(i) = dz(i) * 1000.0 * betae
                wi(i) = 1.0 - betaa / betae
                gg = asy
	        x1 = gg
                x2 = x1 * gg
                x3 = x2 * gg
                x4 = x3 * gg
                wwi(i,1) = 3.0 * x1
	        wwi(i,2) = 5.0 * x2
                wwi(i,3) = 7.0 * x3
                wwi(i,4) = 9.0 * x4
	     endif
           endif
10	continue
	return
	end
! Fu 07-08-98	
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
      block data
! **********************************************************************
! Double-Gauss quadratures and weights (Sykes, 1951).
! **********************************************************************
      implicit none
      real a,u
      common /dis/ a(4)
      common /pointlg/ u(4)
      data a / 0.5, 0.5, 0.5, 0.5 /
      data u / -0.7886752, -0.2113247, 0.2113247, 0.7886752 /
      end
      
      block data ice1
! ****************************************************************
! ap and bp are empirical coefficients of Eqs. (2.9) and (2.10) to
! calculate the extiction coefficient (1/m) and single scattering 
! albedo, cps and dps are empirical coefficients of Eq. (2.13) to
! compute the expansion coefficients of the phase function (1, 2, 
! 3, 4) in the solar bands, cpir is the empirical coefficients of 
! Eq. (2.15) to calculate the asymmetry factor in the IR bands (Fu
! and Liou, 1992). The units of mean effective size and ice water
! content are um and g/m*m*m, respectively, in these equations.
! ****************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      real ap,bp,cps,dps,cpir
      common /ic1/ ap(3,mbx), bp(4,mbx), cps(4,4,mbsx), dps(4,mbsx),   &
                     cpir(4,mbirx)
      data ap / -6.656e-3,          3.686,           0.00,   &
                  -6.656e-3,          3.686,           0.00,   &
                  -6.656e-3,          3.686,           0.00,   &
                  -6.656e-3,          3.686,           0.00,   &
                  -6.656e-3,          3.686,           0.00,   &
                  -6.656e-3,          3.686,           0.00,   &
                  -7.770e-3,          3.734,          11.85,   &
                  -8.088e-3,          3.717,          17.17,   &
                  -8.441e-3,          3.715,          19.48,   &
                  -9.061e-3,          3.741,          26.48,   &
                  -9.609e-3,          3.768,          34.11,   &
                  -1.153e-2,          4.109,          17.32,   &
                  -8.294e-3,          3.925,          1.315,   &
                  -1.026e-2,          4.105,          16.36,   &
                  -1.151e-2,          4.182,          31.13,   &
                  -1.704e-2,          4.830,          16.27,   &
                  -1.741e-2,          5.541,         -58.42,   &
                  -7.752e-3,          4.624,         -42.01 /
      data bp / .10998E-05, -.26101E-07,  .10896E-08, -.47387E-11,   &
                  .20208E-04,  .96483E-05,  .83009E-07, -.32217E-09,   &
                  .13590E-03,  .73453E-03,  .28281E-05, -.18272E-07,   &
                 -.16598E-02,  .20933E-02, -.13977E-05, -.18703E-07,   &
                  .46180E+00,  .24471E-03, -.27839E-05,  .10379E-07,   &
                  .42362E-01,  .86425E-02, -.75519E-04,  .24056E-06,   &
                  .19960E+00,  .37800E-02, -.14910E-04,  .00000E+00,   &
                  .30140E+00,  .26390E-02, -.11160E-04,  .00000E+00,   &
                  .39080E+00,  .12720E-02, -.55640E-05,  .00000E+00,   &
                  .31050E+00,  .26030E-02, -.11390E-04,  .00000E+00,   &
                  .20370E+00,  .42470E-02, -.18100E-04,  .00000E+00,   &
                  .23070E+00,  .38300E-02, -.16160E-04,  .00000E+00,   &
                  .56310E+00, -.14340E-02,  .62980E-05,  .00000E+00,   &
                  .52070E+00, -.97780E-03,  .37250E-05,  .00000E+00,   &
                  .32540E+00,  .34340E-02, -.30810E-04,  .91430E-07,   &
                  .10280E+00,  .50190E-02, -.20240E-04,  .00000E+00,   &
                  .39640E+00, -.31550E-02,  .64170E-04, -.29790E-06,   &
                  .80790E+00, -.70040E-02,  .52090E-04, -.14250E-06 /
      data cps / .22110E+01, -.10398E-02,  .65199E-04, -.34498E-06,   &
                   .32201E+01,  .94227E-03,  .80947E-04, -.47428E-06,   &
                   .41610E+01,  .74396E-03,  .82690E-04, -.45251E-06,   &
                   .51379E+01,  .51545E-02,  .11881E-04, -.15556E-06,   &
                   .22151E+01, -.77982E-03,  .63750E-04, -.34466E-06,   &
                   .31727E+01,  .15597E-02,  .82021E-04, -.49665E-06,   &
                   .40672E+01,  .25800E-02,  .71550E-04, -.43051E-06,   &
                   .49882E+01,  .86489E-02, -.18318E-04, -.59275E-07,   &
                   .22376E+01,  .10293E-02,  .50842E-04, -.30135E-06,   &
                   .31549E+01,  .47115E-02,  .70684E-04, -.47622E-06,   &
                   .39917E+01,  .82830E-02,  .53927E-04, -.41778E-06,   &
                   .48496E+01,  .15998E-01, -.39320E-04, -.43862E-07,   &
                   .23012E+01,  .33854E-02,  .23528E-04, -.20068E-06,   &
                   .31730E+01,  .93439E-02,  .36367E-04, -.38390E-06,   &
                   .39298E+01,  .16424E-01,  .10502E-04, -.35086E-06,   &
                   .47226E+01,  .25872E-01, -.77542E-04, -.21999E-07,   &
                   .27975E+01,  .29741E-02, -.32344E-04,  .11636E-06,   &
                   .43532E+01,  .11234E-01, -.12081E-03,  .43435E-06,   &
                   .56835E+01,  .24681E-01, -.26480E-03,  .95314E-06,   &
                   .68271E+01,  .42788E-01, -.45615E-03,  .16368E-05,   &
                   .19655E+01,  .20094E-01, -.17067E-03,  .50806E-06,   &
                   .28803E+01,  .36091E-01, -.28365E-03,  .79656E-06,   &
                   .34613E+01,  .58525E-01, -.46455E-03,  .13444E-05,   &
                   .39568E+01,  .81480E-01, -.64777E-03,  .19022E-05 /
      data dps / .12495E+00, -.43582E-03,  .14092E-04, -.69565E-07,   &
                   .12363E+00, -.44419E-03,  .14038E-04, -.68851E-07,   &
                   .12117E+00, -.48474E-03,  .12495E-04, -.62411E-07,   &
                   .11581E+00, -.55031E-03,  .98776E-05, -.50193E-07,   &
                  -.15968E-03,  .10115E-04, -.12472E-06,  .48667E-09,    &
                   .13830E+00, -.18921E-02,  .12030E-04, -.31698E-07 /
      data cpir / .79550,     2.524e-3,    -1.022e-5,     0.000e+0,   &
                    .86010,     1.599e-3,    -6.465e-6,     0.000e+0,   &
                    .89150,     1.060e-3,    -4.171e-6,     0.000e+0,   &
                    .87650,     1.198e-3,    -4.485e-6,     0.000e+0,   &
                    .88150,     9.858e-4,    -3.116e-6,     0.000e+0,   &
                    .91670,     5.499e-4,    -1.507e-6,     0.000e+0,   &
                    .90920,     9.295e-4,    -3.877e-6,     0.000e+0,   &
                    .84540,     1.429e-3,    -5.859e-6,     0.000e+0,   &
                    .76780,     2.571e-3,    -1.041e-5,     0.000e+0,   &
                    .72900,     2.132e-3,    -5.584e-6,     0.000e+0,   &
                    .70240,     4.581e-3,    -3.054e-5,     6.684e-8,   &
                    .22920,     1.724e-2,    -1.573e-4,     4.995e-7 /
      end

      block data water1
! *********************************************************************
! bz, wz and gz are the extinction coefficient(1/km), single scattering
! albedo and asymmetry factor for the water clouds (St II, Sc I, St I,
! As, Ns, Sc II, Cu, and Cb) in different bands.   re is the effective 
! radius and fl is the liquid water content (LWC).  See Tables 4.2-4.4 
! of Fu (1991).
! *********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      real re,fl,bz,wz,gz

      common /wat1/ re(ncx),fl(ncx),bz(ncx,mbx),wz(ncx,mbx),gz(ncx,mbx)
      data re /  4.18,  5.36,  5.89,  6.16,    &
                   9.27,  9.84, 12.10, 31.23 /
      data fl / 0.05, 0.14, 0.22, 0.28,   &
                  0.50, 0.47, 1.00, 2.50 /
      data bz /  15.11,  40.25,  59.81,  72.43,   &
                   83.69,  73.99, 128.17, 120.91,   &
                   15.74,  41.70,  61.52,  74.47,   &
                   85.78,  75.59, 130.46, 121.84,   &
                   16.38,  43.52,  64.84,  77.97,   &
                   87.31,  77.36, 134.30, 124.06,   &
                   17.57,  45.78,  66.44,  80.15,   &
                   90.49,  79.90, 137.56, 125.92,   &
                   18.19,  46.63,  69.39,  82.20,   &
                   91.46,  79.99, 138.21, 126.08,   &
                   21.30,  51.88,  77.77,  87.02,   &
                   94.91,  83.55, 143.46, 128.45,   &
                   22.44,  57.35,  84.41, 103.50,   &
                  103.49,  84.17, 152.77, 132.07,   &
                   18.32,  52.69,  76.67, 100.31,   &
                  105.46,  92.86, 157.82, 133.03,   &
                   17.27,  50.44,  74.18,  96.76,   &
                  105.32,  95.25, 158.07, 134.48,   &
                   13.73,  44.90,  67.70,  90.85,   &
                  109.16, 105.48, 163.11, 136.21,   &
                   10.30,  36.28,  57.23,  76.43,   &
                  106.45, 104.90, 161.73, 136.62,   &
                    7.16,  26.40,  43.51,  57.24,   &
                   92.55,  90.55, 149.10, 135.13,   &
                    6.39,  21.00,  33.81,  43.36,   &
                   66.90,  63.58, 113.83, 125.65,   &
                   10.33,  30.87,  47.63,  60.33,   &
                   79.54,  73.92, 127.46, 128.21,   &
                   11.86,  35.64,  54.81,  69.85,   &
                   90.39,  84.16, 142.49, 135.25,   &
                   10.27,  33.08,  51.81,  67.26,   &
                   93.24,  88.60, 148.71, 140.42,   &
                    6.72,  24.09,  39.42,  51.68,   &
                   83.34,  80.72, 140.14, 143.57,   &
                    3.92,  14.76,  25.32,  32.63,   &
                   60.85,  58.81, 112.30, 145.62 /
      data wz / .999999, .999999, .999999, .999999,   &
                  .999998, .999999, .999998, .999997,   &
                  .999753, .999700, .999667, .999646,   &
                  .999492, .999470, .999344, .998667,   &
                  .995914, .994967, .994379, .993842,   &
                  .991385, .990753, .988908, .974831,   &
                  .983761, .978981, .976568, .974700,   &
                  .963466, .959934, .953865, .897690,   &
                  .702949, .683241, .679723, .669045,   &
                  .642616, .632996, .629776, .588820,   &
                  .947343, .929619, .924806, .914557,   &
                  .877169, .867047, .853661, .737426,   &
                  .919356, .896274, .885924, .881097,   &
                  .812772, .781637, .775418, .637341,   &
                  .874717, .861122, .847850, .851677,   &
                  .787171, .772952, .753143, .618656,   &
                  .764750, .752410, .736529, .743435,   &
                  .671272, .659392, .639492, .549941,   &
                  .807536, .808700, .795994, .805489,   &
                  .750577, .755524, .709472, .571989,   &
                  .753346, .772026, .767273, .777079,   &
                  .751264, .760973, .712536, .568286,   &
                  .632722, .676332, .684631, .693552,   &
                  .707986, .717724, .682430, .552867,   &
                  .288885, .348489, .371653, .380367,   &
                  .454540, .465769, .475409, .493881,   &
                  .261827, .306283, .321340, .333051,   &
                  .392917, .406876, .417450, .484593,   &
                  .295804, .339929, .352494, .365502,   &
                  .416229, .430369, .435267, .491356,   &
                  .301214, .354746, .369346, .381906,   &
                  .433602, .447397, .447406, .486968,   &
                  .243714, .318761, .344642, .352770,   &
                  .427906, .438979, .445972, .477264,   &
                  .109012, .187230, .226849, .224976,   &
                  .331382, .335917, .374882, .457067 /
      data gz / .838, .839, .844, .847,   &
                  .849, .860, .853, .859,   &
                  .809, .810, .819, .823,   &
                  .823, .849, .833, .843,   &
                  .774, .787, .781, .792,   &
                  .812, .836, .815, .833,   &
                  .801, .802, .793, .793,   &
                  .814, .829, .818, .832,   &
                  .877, .873, .879, .880,   &
                  .885, .899, .891, .908,   &
                  .783, .769, .777, .756,   &
                  .764, .776, .770, .797,   &
                  .818, .805, .824, .830,   &
                  .815, .801, .820, .845,   &
                  .810, .802, .826, .840,   &
                  .829, .853, .840, .868,   &
                  .774, .766, .799, .818,   &
                  .815, .869, .834, .869,   &
                  .734, .728, .767, .797,   &
                  .796, .871, .818, .854,   &
                  .693, .688, .736, .772,   &
                  .780, .880, .808, .846,   &
                  .643, .646, .698, .741,   &
                  .759, .882, .793, .839,   &
                  .564, .582, .637, .690,   &
                  .719, .871, .764, .819,   &
                  .466, .494, .546, .609,   &
                  .651, .823, .701, .766,   &
                  .375, .410, .455, .525,   &
                  .583, .773, .637, .710,   &
                  .262, .301, .334, .406,   &
                  .485, .695, .545, .631,   &
                  .144, .181, .200, .256,   &
                  .352, .562, .413, .517,   &
                  .060, .077, .088, .112,   &
                  .181, .310, .222, .327 /
      end

      block data rayle1
!  *************************************************************
!  ri is the coefficient in Eq.(4.8) of Fu (1991) to compute the 
!  optical depth due to Rayleigh scattering in the solar bands.
!  *************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      real ri,rix,riy
      
      common /ray1/ ri(mbsx),rix(10),riy(10)

      data ri / 0.9022e-5, 0.5282e-6, 0.5722e-7,   &
                      0.1433e-7, 0.4526e-8, 0.1529e-8 /
!!	data rix /1.92E-04,1.23E-04,7.37E-05,4.62E-05,3.58E-05,
!!     &            2.43E-05,1.30E-05,6.52E-06,3.46E-06,1.76E-06/

!!!  Improved Solar Spectral Energy weighting ( for vis subband )
!!!	data rix /1.77E-04,1.19E-04,6.20E-05,4.48E-05,3.40E-05,
!!!     &              2.28E-05,1.18E-05,6.25E-06,3.40E-06,1.79E-06/

!!! RECOMPUTED FROM 111 SUBBAND FU CODE.
!!!  ENERGY WEIGHTED - LINEAR Averaging
	data rix / 1.775E-04,1.222E-04,6.606E-05,4.588E-05,3.524E-05,   &
                   2.420E-05,1.244E-05,6.510E-06,3.475E-06,1.792E-06/

!!!  ENERGY WEIGHTED - LOGARITHMETIC Averaging
!!!	data riy / 1.761E-04,1.214E-04,6.512E-05,4.583E-05,3.510E-05,
!!!     &             2.401E-05,1.209E-05,6.437E-06,3.396E-06,1.767E-06/

!!!  ENERGY WEIGHTED - LOG(1-6) Lin(7-10) Averaging
!!!	data riy / 1.761E-04,1.214E-04,6.512E-05,4.583E-05,3.510E-05,
!!!     &             2.401E-05, 1.244E-05,6.510E-06,3.475E-06,1.792E-06/

!!!  ENERGY WEIGHTED - LOG(1-7) Lin(8-10) Averaging
	data riy / 1.761E-04,1.214E-04,6.512E-05,4.583E-05,3.510E-05,   &
                   2.401E-05,1.209E-05, 6.510E-06,3.475E-06,1.792E-06/

      end

      block data rain1
!  ********************************************************************
!  brn,  wrnf and  grn  are  the extinction coefficient (1/km),  single
!  scattering  albedo  and  asymmetry  factor  for  the rain.  The size
!  distribution of  rain  is  in the form of a truncated constant-slope 
!  gamma function (Manton and Cotton, 1977)  where rmin = 60 um, rmax =
!  1800 um,  rc = 162 um,  density of water = 1 g/cm**3, and rain water
!  content (rwc) = 0.5 g/m**3.
!                        Jan. 19, 1993
!  ********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      real rwc,brn,wrnf,grn
      common /rai1/ rwc, brn(mbx), wrnf(mbx), grn(mbx)
      data rwc / 0.5 /
      data brn /  1.5377, 1.5377, 1.5379, 1.5385, 1.5396, 1.5417,   &
                    1.5454, 1.5478, 1.5512, 1.5559, 1.5600, 1.5642,   &
                    1.5647, 1.5741, 1.5862, 1.5993, 1.6149, 1.6765 /
      data wrnf /.999932, .97096, .74627, .56719, .53023, .53815,   &
                    .53233, .52884, .53192, .52969, .52716, .52321,   &
                    .51904, .53859, .55169, .55488, .55334, .55218 /
      data grn / .88323, .89067, .92835, .96626, .97553, .96626,   &
                   .97226, .97663, .97216, .97467, .97745, .98156,   &
                   .98584, .96374, .94218, .93266, .92990, .90729 /
      end

      block data graup1
!  ********************************************************************
!  bg,  wgf  and  gg are  the  extinction  coefficient (1/km),   single
!  scattering  albedo  and  asymmetry  factor for the graupel. The size
!  distribution of graupel is in the form of a truncated constant-slope 
!  gamma function (Manton and Cotton, 1977)  where rmin = 60 um, rmax =
!  5000 um, rc = 500 um, density of graupel = 0.6 g/cm**3, and  graupel
!  water content (gwc) = 0.5 g/m**3.
!                        Jan. 19, 1993
!  ********************************************************************
      USE RadParams
      implicit none
!##      include 'rad_0698.h'
      real gwc,bg,wgf,gg
      common /gra1/ gwc, bg(mbx), wgf(mbx), gg(mbx)
      data gwc / 0.5 /
      data bg /  0.83939,0.83940,0.83940,0.83941,0.83946,0.83951,   &
                   0.83967,0.83979,0.83995,0.84029,0.84058,0.84097,   &
                   0.84143,0.84286,0.84418,0.84825,0.85421,0.87477 /
      data wgf / 0.999911,0.97115,0.56192,0.53156,0.52579,0.53846,   &
                    0.53296,0.53017,0.53182,0.53180,0.52959,0.52446,   &
                    0.52342,0.54914,0.55258,0.54307,0.53160,0.55474 /
      data gg /  0.89218,0.89940,0.96820,0.97816,0.98141,0.96373,   &
                   0.97173,0.97559,0.97330,0.97327,0.97626,0.98274,   &
                   0.98396,0.94673,0.94213,0.95539,0.97097,0.93183 /
      end

      block data legend
! **********************************************************************
! p0d(4), p1d(4), p2d(4), and p3d(4) are Legendre polynomials p0(x), 
! p1(x), p2(x), and p3(x) when x = u(1), u(2), u(3), and u(4).
! **********************************************************************
      implicit none
      real p0d, p1d, p2d, p3d
      common /legen/ p0d(4), p1d(4), p2d(4), p3d(4)
      data p0d /  .100000E+01,  .100000E+01,  .100000E+01, .100000E+01 /
      data p1d / -.788675E+00, -.211325E+00,  .211325E+00, .788675E+00 /
      data p2d /  .433013E+00, -.433013E+00, -.433013E+00, .433013E+00 /
      data p3d / -.433940E-01,  .293394E+00, -.293394E+00, .433940E-01 /
      end

      block data legenf
! *********************************************************************
! p11d(4,4), p22d(4,4), and p33d(4,4) are defined as 0.5*p1d(i)*p1d(j),
! 0.5*p2d(i)*p2d(j), and 0.5*p3d(i)*p3d(j), respectively.
! *********************************************************************
      implicit none
      real p11d, p22d, p33d
      common /legen1/ p11d(4,4), p22d(4,4), p33d(4,4)
      data p11d / .311004E+00, .833334E-01,-.833334E-01,-.311004E+00,   &
                  .833334E-01, .223291E-01,-.223291E-01,-.833334E-01,   &
                 -.833334E-01,-.223291E-01, .223291E-01, .833334E-01,   &
                 -.311004E+00,-.833334E-01, .833334E-01, .311004E+00 /
      data p22d / .937501E-01,-.937501E-01,-.937501E-01, .937501E-01,   &
                 -.937501E-01, .937501E-01, .937501E-01,-.937501E-01,   &
                 -.937501E-01, .937501E-01, .937501E-01,-.937501E-01,   &
                  .937501E-01,-.937501E-01,-.937501E-01, .937501E-01 /
      data p33d / .941520E-03,-.636577E-02, .636577E-02,-.941520E-03,   &
                 -.636577E-02, .430400E-01,-.430400E-01, .636577E-02,   &
                  .636577E-02,-.430400E-01, .430400E-01,-.636577E-02,   &
                 -.941520E-03, .636577E-02,-.636577E-02, .941520E-03 /
      end
      
!---------- 4/1/97 (6) -- Replaces old ckd1 block data.
      subroutine ckd1_init( isolar_spectrum,lray )
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 to 
! one. fko3 is the corresponding ozone absorption coefficient in units
! of (cm-atm)**-1 (Fu, 1991). The spectral region is from 50000 cm**-1 
! to 14500 cm**-1.
! *********************************************************************
!      common /band1/ hk(10), fko3(10)
!      data hk / .24, .16, .24, .28, .03, 
!     1            .016, .01, .008, .008, .008 /
!      data fko3 / .2204e-08,.1207e-01,.4537e-01,.1032e+00,.1740e+00,
!     1           .1210e+01,.7367e+01,.2050e+02,.8100e+02,.2410e+03 /
! *********************************************************************
! The treatment of O3 absorption here is modified to consider  the
! dependence of aerosol single-scattering properties on wavelength
! in the first band.  We divide this band into 10 subintervals  as
! 175.439 - 224.719 - 243.902 - 285.714 - 298.507 - 322.500 -
! 357.500 - 437.500 - 497.500 - 595.000 - 692.500 following Grant
! (personal communication, 1996).
! *********************************************************************
!---------- 4/1/97 (6)
      common /band1/ hk(10), fko3(10),sol_spect(0:7),fkh2o(10)

	real ,dimension(10):: fko3_1,fko3_2,fko3_3,fko3_4
	logical lray
!!!   OLD FU Ozone Coefficients
      fko3_1 = (/ .2751e+02,.1564e+03,.1746e+03,.2684e+02,.1962e+01,   &
                  .7402e-01,.7295e-03,.1239e-01,.8011e-01,.8011e-01 /)

!!! RECOMPUTED FROM 111 SUBBAND FU CODE.
!!! ENERGY WEIGHTED - Linear Average
!!!    fko3_2 = ( /3.87E+01,1.66E+02,1.94E+02,2.93E+01,3.22E+00,
!!!     1            8.38E-02,7.38E-04,1.25E-02,8.45E-02,7.88E-02/)   
!!! ENERGY WEIGHTED - LOG Average
!!!	fko3_3 =( /3.28E+01,1.56E+02,1.75E+02,2.70E+01,2.19E+00,
!!!   &		   2.98E-02,1.94E-04,1.02E-02,7.76E-02,7.19E-02/)
!!! ENERGY WEIGHTED - LOG(1-6) Lin(7-10) Average
!!!	data fko3_4 /3.28E+01,1.56E+02,1.75E+02,2.70E+01,2.19E+00,
!!!     &	     2.98E-02,7.38E-04,1.25E-02,8.45E-02,7.88E-02/

!!! ENERGY WEIGHTED - LOG(1-5) Lin(6-10) Average
	 fko3_4 = (/3.28E+01,1.56E+02,1.75E+02,2.70E+01,2.19E+00,   &
      		    8.38E-02,7.38E-04,1.25E-02,8.45E-02,7.88E-02/)

! CHOU
!       data fkh2o / 6*0.0,4*0.00075 /
! Q.Fu
        fkh2o=(/0.,0.,0.,0.,0.,0.,0.,7.63E-05,1.90E-03,2.57E-03/)
	
! SOLAR SPECTRA Choice::

	if (lray) then
	 fko3=fko3_4  !!! LOG & Lin wgt
	else
	 fko3=fko3_1 !! OLD
	endif



	select case (isolar_spectrum)
	case default
        sol_spect =(/1351.0, 619.618, 484.295, 149.845, 48.7302,& 
                    31.6576, 5.79927,11.0 /)
        hk =(/ .0012221, .0014711, .0100210, .0097280, .0240560,   &
                .0528740, .1766475, .1897695, .2870026, .2472082 /)
!--------------------------------------------------------------------
      	case (1)      
        sol_spect = (/1373.013,626.929,493.654,155.256, 51.183, &
                     28.516,  5.550, 11.924/)
         hk  = (/0.00117,0.00145,0.01060,0.01060,0.02618,   &
         0.05375,0.17448,0.18855,0.28872,0.24451/)
 !-------------------------------------------------------- 
      	 case (2)!!!  Mod3.7_cebchkur     
         sol_spect=(/1362.119,   &
         626.734,484.953,153.697, 51.119, 28.420,  5.542,   &
          11.910/)
          hk = (/0.00125,0.00154,0.01040,0.01009,0.02481,   &
                0.05374,0.17567,0.18909,0.28927,0.24413/)
 !-------------------------------------------------------- 
    	 case (3)!!!  Mod3.7_chkur        
        sol_spect =(/1359.751,   &
         624.364,484.953,153.697, 51.119, 28.420,  5.542,   &
          11.910/)
         hk = (/0.00128,0.00158,0.01108,0.01120,0.02336,   &
         0.05156,0.17472,0.18980,0.29037,0.24506/)
 !-------------------------------------------------------- 
      	case (4)!!!  Mod3.7_Newkur       
       sol_spect =(/1367.997,   &
         630.440,487.123,153.697, 51.119, 28.420,  5.542,   &
          11.910/)
        hk =(/0.00127,0.00142,0.01142,0.01087,0.02440,   &
         0.05532,0.17325,0.19032,0.28889,0.24283/)
 !-------------------------------------------------------- 
 	case (5)!!!  Mod3.7_thkur        
        sol_spect =(/1376.231,   &
         634.524,491.276,153.697, 51.119, 28.420,  5.542,   &
          11.910/)
        hk =(/   &
         0.00125,0.00152,0.01005,0.00970,0.02423,   &
         0.05515,0.17862,0.19247,0.28592,0.24111/)
 !-------------------------------------------------------- 
	 case (6)!!!
!SPECIAL FU_111 subband integration of vis-subband
      sol_spect =(/1351.0 , &
             619.618, 484.295, 149.845, 48.7302, 31.6576, 5.79927,   &
             11.0 /)
       hk =(/&
         0.001113,0.001452,0.010599,0.010586,0.026237,   &
         0.053743,0.174724,0.188548,0.295757,0.237242/)
	end select
!	print*,'SOLSPEC',isolar_spectrum,sol_spect(0:7)

      end

      block data ckd2
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
! pressures,  and eight cumulative probabilities  ( Fu,  1991 ). The
! spectral region is from 14500 to 7700 cm**-1.
! *********************************************************************
      common /band2/ hk(8), coeh2o(3,11,8)
      data hk / .71, .11, .06, .06, .04, .016, .0034, .0006 /
!   .343849E+03    .532724E+02    .290577E+02    .290577E+02    .193718E+02
!   .774872E+01    .164660E+01    .290577E+00
      data ( ( ( coeh2o(k,j,i), i = 1, 8 ), j = 1, 11 ), k = 1, 3 ) /   &
      -.1735E+02,-.1407E+02,-.1268E+02,-.1131E+02,-.9261E+01,-.6666E+01,   &
      -.3937E+01,-.5448E+00,-.1690E+02,-.1365E+02,-.1232E+02,-.1101E+02,   &
      -.9058E+01,-.6574E+01,-.3914E+01,-.5529E+00,-.1643E+02,-.1323E+02,   &
      -.1195E+02,-.1068E+02,-.8840E+01,-.6475E+01,-.3889E+01,-.6143E+00,   &
      -.1598E+02,-.1282E+02,-.1157E+02,-.1035E+02,-.8598E+01,-.6339E+01,   &
      -.3848E+01,-.6636E+00,-.1551E+02,-.1241E+02,-.1119E+02,-.1001E+02,   &
      -.8342E+01,-.6178E+01,-.3788E+01,-.8181E+00,-.1506E+02,-.1201E+02,   &
      -.1082E+02,-.9692E+01,-.8073E+01,-.6017E+01,-.3703E+01,-.9003E+00,   &
      -.1446E+02,-.1154E+02,-.1042E+02,-.9332E+01,-.7810E+01,-.5846E+01,   &
      -.3576E+01,-.1083E+01,-.1394E+02,-.1112E+02,-.1005E+02,-.8992E+01,   &
      -.7548E+01,-.5674E+01,-.3477E+01,-.1266E+01,-.1351E+02,-.1076E+02,   &
      -.9722E+01,-.8702E+01,-.7334E+01,-.5531E+01,-.3401E+01,-.1524E+01,   &
      -.1311E+02,-.1044E+02,-.9422E+01,-.8423E+01,-.7117E+01,-.5383E+01,   &
      -.3410E+01,-.1785E+01,-.1274E+02,-.1015E+02,-.9162E+01,-.8190E+01,   &
      -.6949E+01,-.5236E+01,-.3477E+01,-.2082E+01, .2407E-02, .2847E-02,   &
       .3768E-02, .4626E-02, .5631E-02, .4542E-02, .3475E-02,-.3085E-02,   &
       .2428E-02, .2805E-02, .3412E-02, .3893E-02, .4773E-02, .3998E-02,   &
       .2742E-02,-.2556E-02, .2428E-02, .2721E-02, .3077E-02, .3161E-02,   &
       .4019E-02, .3224E-02, .2512E-02,-.1884E-02, .2449E-02, .2617E-02,   &
       .2763E-02, .2658E-02, .3286E-02, .2617E-02, .1989E-02,-.1740E-02,   &
       .2512E-02, .2470E-02, .2470E-02, .2282E-02, .2512E-02, .1926E-02,   &
       .1465E-02,-.2612E-02, .2554E-02, .2303E-02, .2303E-02, .1842E-02,   &
       .2030E-02, .1340E-02, .1068E-02,-.1413E-02, .2449E-02, .2198E-02,   &
       .2030E-02, .1465E-02, .1528E-02, .9838E-03, .1005E-02,-.1099E-02,   &
       .2868E-02, .2198E-02, .1968E-02, .1382E-02, .1172E-02, .5652E-03,   &
       .6070E-03,-.1662E-02, .3077E-02, .2219E-02, .1800E-02, .1277E-02,   &
       .1005E-02, .3349E-03, .2512E-03,-.1195E-02, .3182E-02, .2219E-02,   &
       .1758E-02, .1172E-02, .7326E-03, .4815E-03, .6280E-04,-.1880E-02,   &
       .3265E-02, .2114E-02, .1696E-02, .1298E-02, .4187E-03, .4187E-03,   &
      -.3768E-03,-.1467E-02,-.1180E-04,-.1294E-04,-.1142E-04,-.7232E-05,   &
      -.8754E-05,-.1484E-04,-.8373E-05, .1028E-04,-.1218E-04,-.1142E-04,   &
      -.9515E-05,-.1522E-05,-.9134E-05,-.1484E-04,-.3425E-05, .1142E-06,   &
      -.1294E-04,-.9895E-05,-.7231E-05,-.4187E-05,-.7612E-05,-.3806E-05,   &
       .1522E-05,-.3882E-05,-.1256E-04,-.8754E-05,-.7612E-05,-.6470E-05,   &
      -.4948E-05,-.3425E-05, .4948E-05,-.1054E-04,-.1370E-04,-.6089E-05,   &
      -.8373E-05,-.5709E-05,-.3045E-05,-.3806E-05, .5328E-05, .8678E-05,   &
      -.1370E-04,-.6851E-05,-.8373E-05,-.1522E-05,-.3425E-05, .0000E+00,   &
       .1256E-04,-.1572E-04,-.1484E-04,-.7231E-05,-.7992E-05,-.4567E-05,   &
      -.2664E-05,-.3807E-06,-.1522E-05, .2169E-05,-.1713E-04,-.9515E-05,   &
      -.6089E-05,-.6851E-05,-.3045E-05,-.1142E-05, .1903E-05, .9363E-05,   &
      -.1560E-04,-.9134E-05,-.5328E-05,-.4948E-05, .0000E+00, .7611E-06,   &
      -.6851E-05, .1252E-04,-.1522E-04,-.8373E-05,-.6089E-05,-.6089E-05,   &
      -.3805E-06,-.1142E-05,-.3807E-06, .2512E-05,-.1599E-04,-.7231E-05,   &
      -.5709E-05,-.4567E-05, .1522E-05,-.2284E-05,-.3941E-10, .5290E-05/
      end
      
      block data ckd3
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
! pressures,  and twelve cumulative probabilities ( Fu,  1991 ). The
! spectral region is from 7700 to 5250 cm**-1.
! *********************************************************************
      common /band3/ hk(12), coeh2o(3,11,12)
      data hk / .34, .11, .1, .09, .12, .1,   &
                  .06, .04, .026, .01, .0035, .0005 /
!   .509474E+02    .164830E+02    .149845E+02    .134861E+02    .179814E+02
!   .149845E+02    .899071E+01    .599381E+01    .389597E+01    .149845E+01
!   .524458E+00    .749226E-01
      data ( ( ( coeh2o(k,j,i), i = 1, 12 ), j = 1, 11 ), k = 1, 3 ) /   &
      -.1900E+02,-.1515E+02,-.1344E+02,-.1224E+02,-.1081E+02,-.9337E+01,   &
      -.7965E+01,-.6585E+01,-.4578E+01,-.2247E+01, .1747E+00, .3083E+01,   &
      -.1854E+02,-.1471E+02,-.1300E+02,-.1181E+02,-.1039E+02,-.8927E+01,   &
      -.7576E+01,-.6238E+01,-.4317E+01,-.2119E+01, .1888E+00, .3033E+01,   &
      -.1808E+02,-.1426E+02,-.1257E+02,-.1137E+02,-.9966E+01,-.8513E+01,   &
      -.7177E+01,-.5885E+01,-.4053E+01,-.1977E+01, .2245E+00, .3005E+01,   &
      -.1763E+02,-.1381E+02,-.1213E+02,-.1094E+02,-.9542E+01,-.8094E+01,   &
      -.6779E+01,-.5524E+01,-.3788E+01,-.1796E+01, .2961E+00, .2828E+01,   &
      -.1716E+02,-.1337E+02,-.1170E+02,-.1051E+02,-.9116E+01,-.7677E+01,   &
      -.6381E+01,-.5153E+01,-.3493E+01,-.1607E+01, .3850E+00, .2660E+01,   &
      -.1670E+02,-.1295E+02,-.1127E+02,-.1008E+02,-.8690E+01,-.7265E+01,   &
      -.5991E+01,-.4799E+01,-.3212E+01,-.1438E+01, .4582E+00, .2588E+01,   &
      -.1596E+02,-.1231E+02,-.1067E+02,-.9501E+01,-.8151E+01,-.6793E+01,   &
      -.5588E+01,-.4458E+01,-.2940E+01,-.1257E+01, .4888E+00, .2260E+01,   &
      -.1530E+02,-.1184E+02,-.1017E+02,-.8992E+01,-.7661E+01,-.6369E+01,   &
      -.5213E+01,-.4145E+01,-.2701E+01,-.1108E+01, .4239E+00, .1974E+01,   &
      -.1481E+02,-.1144E+02,-.9756E+01,-.8573E+01,-.7255E+01,-.5994E+01,   &
      -.4868E+01,-.3829E+01,-.2485E+01,-.9738E+00, .3343E+00, .1667E+01,   &
      -.1439E+02,-.1108E+02,-.9360E+01,-.8183E+01,-.6885E+01,-.5646E+01,   &
      -.4559E+01,-.3555E+01,-.2314E+01,-.8904E+00, .2169E+00, .1289E+01,   &
      -.1402E+02,-.1073E+02,-.8987E+01,-.7817E+01,-.6551E+01,-.5335E+01,   &
      -.4278E+01,-.3316E+01,-.2147E+01,-.8695E+00, .1587E-01, .8658E+00,   &
       .1132E-01, .8855E-02, .6698E-02, .5296E-02, .4396E-02, .3370E-02,   &
       .3245E-02, .4145E-02, .4731E-02, .4756E-02, .3116E-02,-.2763E-02,   &
       .1135E-01, .8917E-02, .6657E-02, .5170E-02, .4207E-02, .3056E-02,   &
       .2868E-02, .3433E-02, .3726E-02, .4109E-02, .2836E-02,-.3119E-02,   &
       .1135E-01, .8980E-02, .6615E-02, .5045E-02, .4061E-02, .2847E-02,   &
       .2491E-02, .2847E-02, .2910E-02, .2671E-02, .2396E-02,-.3245E-02,   &
       .1135E-01, .9043E-02, .6594E-02, .4940E-02, .3914E-02, .2638E-02,   &
       .2156E-02, .2261E-02, .2051E-02, .1978E-02, .1566E-02,-.3203E-02,   &
       .1139E-01, .9085E-02, .6531E-02, .4835E-02, .3768E-02, .2428E-02,   &
       .1842E-02, .1612E-02, .1591E-02, .1279E-02, .7201E-03,-.2763E-02,   &
       .1143E-01, .9085E-02, .6447E-02, .4752E-02, .3684E-02, .2261E-02,   &
       .1570E-02, .1235E-02, .1151E-02, .7243E-03, .6489E-04,-.2240E-02,   &
       .1135E-01, .9001E-02, .5694E-02, .4438E-02, .3412E-02, .1968E-02,   &
       .1235E-02, .9420E-03, .8792E-03, .5045E-03,-.1821E-03,-.1936E-02,   &
       .1174E-01, .9273E-02, .5882E-02, .4689E-02, .3454E-02, .1947E-02,   &
       .1151E-02, .6070E-03, .6698E-03, .9420E-04,-.6740E-03,-.2707E-02,   &
       .1218E-01, .9336E-02, .6050E-02, .4731E-02, .3475E-02, .1863E-02,   &
       .1151E-02, .4605E-03, .3768E-03,-.1214E-03,-.4396E-03,-.1903E-02,   &
       .1235E-01, .9294E-02, .6029E-02, .4584E-02, .3370E-02, .1800E-02,   &
       .1068E-02, .2303E-03, .1675E-03,-.4501E-03,-.7571E-03,-.1149E-02,   &
       .1233E-01, .9315E-02, .6029E-02, .4438E-02, .3203E-02, .1842E-02,   &
       .9629E-03, .0000E+00,-.2198E-03,-.5338E-03,-.9721E-03,-.7661E-03,   &
      -.3692E-04,-.3844E-04,-.2588E-04,-.1180E-04,-.1066E-04,-.3426E-05,   &
      -.2664E-05, .7611E-06, .6089E-05,-.4568E-06,-.2077E-04,-.1142E-04,   &
      -.3730E-04,-.3806E-04,-.2360E-04,-.1256E-04,-.1180E-04,-.4567E-05,   &
      -.3425E-05,-.2284E-05,-.1522E-05,-.4225E-05,-.9940E-05,-.4187E-05,   &
      -.3501E-04,-.3844E-04,-.2131E-04,-.1256E-04,-.9896E-05,-.3806E-05,   &
      -.4186E-05, .7612E-06,-.1903E-05, .4110E-05, .1789E-05,-.2169E-04,   &
      -.3425E-04,-.3882E-04,-.1941E-04,-.1294E-04,-.9515E-05,-.4567E-05,   &
      -.4186E-05, .1522E-05,-.4187E-10, .4605E-05,-.2588E-05, .6470E-05,   &
      -.3501E-04,-.3730E-04,-.1751E-04,-.1332E-04,-.1066E-04,-.3806E-05,   &
      -.4567E-05,-.1142E-05,-.3045E-05, .1104E-05,-.1058E-04, .2816E-04,   &
      -.3578E-04,-.3501E-04,-.1751E-04,-.1332E-04,-.1218E-04,-.3806E-05,   &
      -.3425E-05,-.3806E-06,-.4187E-05,-.6090E-06,-.6965E-05,-.3463E-04,   &
      -.3578E-04,-.3349E-04,-.1675E-04,-.9895E-05,-.9515E-05,-.6090E-05,   &
      -.6470E-05,-.3807E-06,-.5328E-05,-.4186E-06,-.3996E-05, .2074E-04,   &
      -.3540E-04,-.3083E-04,-.1789E-04,-.9896E-05,-.1104E-04,-.6470E-05,   &
      -.5709E-05, .3425E-05,-.4567E-05, .3463E-05, .5633E-05,-.3159E-05,   &
      -.3730E-04,-.2740E-04,-.1484E-04,-.1066E-04,-.1142E-04,-.6470E-05,   &
      -.6470E-05, .1522E-05,-.1522E-05,-.3045E-05, .3197E-05,-.1039E-04,   &
      -.3425E-04,-.2284E-04,-.1370E-04,-.1028E-04,-.1104E-04,-.8373E-05,   &
      -.4948E-05, .1903E-05,-.7612E-06,-.1104E-05, .2455E-05,-.3805E-07,   &
      -.3235E-04,-.2093E-04,-.1294E-04,-.1142E-04,-.1180E-04,-.6851E-05,   &
      -.3045E-05,-.7611E-06, .1256E-05,-.7231E-06, .9924E-05, .3578E-05/
      end
      
      block data ckd4
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
! pressures,  and seven cumulative probabilities ( Fu,  1991 ). The
! spectral region is from 5250 to 4000 cm**-1.
! *********************************************************************
      common /band4/ hk(7), coeh2o(3,11,7)
      data hk / .52, .21, .11, .1, .04, .015, .005 /
!   .253397E+02    .102333E+02    .536032E+01    .487302E+01    .194921E+01
!   .730953E+00    .243651E+00
      data ( ( ( coeh2o(k,j,i), i = 1, 7 ), j = 1, 11 ), k = 1, 3 ) /   &
      -.1722E+02,-.1402E+02,-.1202E+02,-.1001E+02,-.7702E+01,-.5273E+01,   &
      -.6530E+00,-.1677E+02,-.1359E+02,-.1164E+02,-.9662E+01,-.7419E+01,   &
      -.5001E+01,-.6040E+00,-.1630E+02,-.1316E+02,-.1125E+02,-.9303E+01,   &
      -.7092E+01,-.4750E+01,-.5715E+00,-.1584E+02,-.1274E+02,-.1086E+02,   &
      -.8939E+01,-.6751E+01,-.4458E+01,-.4928E+00,-.1538E+02,-.1232E+02,   &
      -.1048E+02,-.8579E+01,-.6399E+01,-.4191E+01,-.4683E+00,-.1493E+02,   &
      -.1192E+02,-.1011E+02,-.8241E+01,-.6065E+01,-.3910E+01,-.4310E+00,   &
      -.1440E+02,-.1145E+02,-.9643E+01,-.7873E+01,-.5710E+01,-.3668E+01,   &
      -.3304E+00,-.1391E+02,-.1104E+02,-.9238E+01,-.7479E+01,-.5367E+01,   &
      -.3387E+01,-.3604E+00,-.1348E+02,-.1069E+02,-.8918E+01,-.7122E+01,   &
      -.5086E+01,-.3152E+01,-.3030E+00,-.1310E+02,-.1037E+02,-.8626E+01,   &
      -.6790E+01,-.4815E+01,-.2945E+01,-.4789E+00,-.1275E+02,-.1011E+02,   &
      -.8347E+01,-.6484E+01,-.4584E+01,-.2788E+01,-.5807E+00, .7934E-02,   &
       .9231E-02, .1005E-01, .9043E-02, .8164E-02, .8980E-02, .6403E-02,   &
       .7954E-02, .9169E-02, .9797E-02, .8687E-02, .7724E-02, .7954E-02,   &
       .6652E-02, .7954E-02, .9043E-02, .9608E-02, .8499E-02, .7347E-02,   &
       .7473E-02, .6382E-02, .7996E-02, .8980E-02, .9378E-02, .8289E-02,   &
       .7264E-02, .6594E-02, .6674E-02, .8059E-02, .8938E-02, .9294E-02,   &
       .8227E-02, .7201E-02, .6678E-02, .7032E-02, .8122E-02, .8896E-02,   &
       .9189E-02, .8038E-02, .7033E-02, .5987E-02, .5475E-02, .8268E-02,   &
       .9064E-02, .8792E-02, .7975E-02, .6573E-02, .5087E-02, .4657E-02,   &
       .8541E-02, .8980E-02, .9085E-02, .7996E-02, .6133E-02, .4501E-02,   &
       .3860E-02, .8813E-02, .9043E-02, .9294E-02, .8122E-02, .5861E-02,   &
       .4354E-02, .3964E-02, .8875E-02, .8834E-02, .9797E-02, .8164E-02,   &
       .5463E-02, .4417E-02, .3270E-02, .8938E-02, .8771E-02, .1005E-01,   &
       .8247E-02, .5589E-02, .4835E-02, .3033E-02,-.1484E-04,-.2169E-04,   &
      -.2436E-04,-.2588E-04,-.1142E-04,-.1142E-05,-.1519E-04,-.1522E-04,   &
      -.2055E-04,-.2131E-04,-.2398E-04,-.4948E-05,-.1675E-04,-.3593E-04,   &
      -.1522E-04,-.2055E-04,-.1865E-04,-.2207E-04,-.4948E-05,-.1180E-04,   &
      -.1237E-04,-.1598E-04,-.2017E-04,-.1903E-04,-.2284E-04,-.1028E-04,   &
      -.1865E-04,-.2381E-04,-.1713E-04,-.2017E-04,-.1827E-04,-.2169E-04,   &
      -.1218E-04,-.9515E-05,-.2415E-04,-.1827E-04,-.2093E-04,-.1637E-04,   &
      -.1827E-04,-.9134E-05,-.8373E-05,-.1243E-04,-.1560E-04,-.1865E-04,   &
      -.1599E-04,-.1256E-04,-.1066E-04,-.1142E-05,-.2181E-04,-.1675E-04,   &
      -.1560E-04,-.1522E-04,-.1675E-04,-.1865E-04,-.1865E-04,-.9522E-05,   &
      -.1332E-04,-.1370E-04,-.1446E-04,-.2055E-04,-.1142E-04,-.2512E-04,   &
      -.3343E-04,-.1294E-04,-.1294E-04,-.1751E-04,-.2512E-04,-.1560E-04,   &
      -.2854E-04,-.7003E-05,-.8753E-05,-.1028E-04,-.1751E-04,-.2512E-04,   &
      -.1713E-04,-.1713E-04,-.1245E-04 /
      end
      
      block data ckd5
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
! pressures,  and twelve cumulative probabilities ( Fu,  1991 ). The
! spectral region is from 4000 to 2850 cm**-1.
! *********************************************************************
      common /band5/ hk(12), coeh2o(3,11,12)
      data hk / .13, .14, .13, .16, .18, .14,    &
                      .07, .02, .016, .008, .004, .002 /
!   .411549E+01    .443207E+01    .411549E+01    .506522E+01    .569837E+01
!   .443207E+01    .221603E+01    .633153E+00    .506522E+00    .253261E+00
!   .126631E+00    .633153E-01
      data ( ( ( coeh2o(k,j,i), i = 1, 12 ), j = 1, 11 ), k = 1, 3 ) /   &
      -.1499E+02,-.1267E+02,-.1118E+02,-.9696E+01,-.7992E+01,-.6323E+01,   &
      -.4414E+01,-.2961E+01,-.1715E+01,-.1406E+00, .1612E+01, .3689E+01,   &
      -.1454E+02,-.1223E+02,-.1075E+02,-.9277E+01,-.7576E+01,-.5915E+01,   &
      -.4043E+01,-.2630E+01,-.1449E+01, .2314E-01, .1708E+01, .3744E+01,   &
      -.1408E+02,-.1178E+02,-.1031E+02,-.8851E+01,-.7154E+01,-.5503E+01,   &
      -.3666E+01,-.2288E+01,-.1141E+01, .2772E+00, .1819E+01, .3788E+01,   &
      -.1363E+02,-.1134E+02,-.9876E+01,-.8423E+01,-.6733E+01,-.5091E+01,   &
      -.3286E+01,-.1938E+01,-.8649E+00, .5349E+00, .1969E+01, .3795E+01,   &
      -.1318E+02,-.1091E+02,-.9452E+01,-.8004E+01,-.6309E+01,-.4677E+01,   &
      -.2904E+01,-.1595E+01,-.5641E+00, .7592E+00, .2109E+01, .3783E+01,   &
      -.1275E+02,-.1048E+02,-.9028E+01,-.7585E+01,-.5892E+01,-.4267E+01,   &
      -.2524E+01,-.1274E+01,-.2782E+00, .9376E+00, .2257E+01, .3714E+01,   &
      -.1180E+02,-.9887E+01,-.8492E+01,-.7014E+01,-.5390E+01,-.3834E+01,   &
      -.2156E+01,-.9775E+00,-.3129E-01, .1151E+01, .2330E+01, .3592E+01,   &
      -.1114E+02,-.9367E+01,-.8002E+01,-.6514E+01,-.4928E+01,-.3435E+01,   &
      -.1835E+01,-.7064E+00, .2153E+00, .1309E+01, .2422E+01, .3488E+01,   &
      -.1074E+02,-.8941E+01,-.7582E+01,-.6116E+01,-.4536E+01,-.3072E+01,   &
      -.1521E+01,-.4651E+00, .4053E+00, .1465E+01, .2374E+01, .3260E+01,   &
      -.1041E+02,-.8545E+01,-.7180E+01,-.5745E+01,-.4177E+01,-.2735E+01,   &
      -.1245E+01,-.2356E+00, .5786E+00, .1516E+01, .2263E+01, .3074E+01,   &
      -.1008E+02,-.8149E+01,-.6804E+01,-.5409E+01,-.3855E+01,-.2427E+01,   &
      -.9857E+00,-.4939E-01, .7060E+00, .1483E+01, .2159E+01, .2745E+01,   &
       .9985E-02, .8373E-02, .7431E-02, .6866E-02, .4584E-02, .2952E-02,   &
       .3098E-02, .3768E-02, .4013E-02, .3960E-02, .3228E-02, .3203E-02,   &
       .1007E-01, .8436E-02, .7368E-02, .6657E-02, .4375E-02, .2617E-02,   &
       .2742E-02, .3286E-02, .3192E-02, .2992E-02, .2612E-02, .1968E-02,   &
       .1019E-01, .8457E-02, .7264E-02, .6426E-02, .4187E-02, .2365E-02,   &
       .2324E-02, .2614E-02, .2736E-02, .2068E-02, .2085E-02, .1005E-02,   &
       .1028E-01, .8478E-02, .7138E-02, .6259E-02, .3998E-02, .2156E-02,   &
       .1926E-02, .1953E-02, .2250E-02, .1844E-02, .1869E-02,-.6489E-03,   &
       .1030E-01, .8478E-02, .7033E-02, .6112E-02, .3852E-02, .1989E-02,   &
       .1716E-02, .1763E-02, .1432E-02, .1193E-02, .1306E-02,-.5861E-03,   &
       .1042E-01, .8499E-02, .6887E-02, .5987E-02, .3768E-02, .1800E-02,   &
       .1549E-02, .1712E-02, .1287E-02, .7389E-03, .7222E-03,-.1130E-02,   &
       .8227E-02, .7201E-02, .6866E-02, .5903E-02, .3412E-02, .1591E-02,   &
       .1402E-02, .1346E-02, .1041E-02, .8185E-03, .3349E-03,-.4815E-03,   &
       .8268E-02, .6992E-02, .7159E-02, .6384E-02, .3286E-02, .1591E-02,   &
       .1271E-02, .1202E-02, .9187E-03, .6531E-03,-.4187E-03,-.7954E-03,   &
       .8478E-02, .7159E-02, .7117E-02, .6447E-02, .3349E-02, .1528E-02,   &
       .9964E-03, .9210E-03, .6112E-03, .6259E-03,-.3768E-03,-.1298E-02,   &
       .8520E-02, .7075E-02, .7096E-02, .6405E-02, .3245E-02, .1528E-02,   &
       .1011E-02, .7877E-03, .7536E-03, .9001E-04,-.6719E-03,-.1026E-02,   &
       .8561E-02, .6950E-02, .7033E-02, .6280E-02, .2993E-02, .1528E-02,   &
       .6698E-03, .5847E-03, .2847E-03,-.6280E-04,-.9420E-03,-.1444E-02,   &
      -.1408E-04,-.2664E-04,-.1180E-04,-.1903E-04,-.9515E-05, .3806E-06,   &
      -.6851E-05,-.3806E-05,-.4834E-05,-.3239E-05,-.2284E-05,-.1028E-04,   &
      -.1484E-04,-.2550E-04,-.1142E-04,-.1827E-04,-.9515E-05, .3805E-06,   &
      -.4948E-05, .3806E-06,-.2664E-06, .1058E-04,-.1012E-04,-.1142E-04,   &
      -.1560E-04,-.2512E-04,-.1256E-04,-.1865E-04,-.9134E-05, .1142E-05,   &
      -.3425E-05, .2474E-05,-.9781E-05,-.1519E-05,-.7916E-05,-.1294E-04,   &
      -.1560E-04,-.2474E-04,-.1180E-04,-.2017E-04,-.7992E-05, .3805E-06,   &
      -.2283E-05,-.4453E-05,-.1180E-05,-.5138E-05,-.4453E-05,-.3425E-05,   &
      -.1522E-04,-.2550E-04,-.9896E-05,-.1903E-04,-.9134E-05,-.1142E-05,   &
      -.7611E-06,-.5252E-05,-.4567E-06,-.4643E-05,-.4567E-06,-.4567E-05,   &
      -.1294E-04,-.2512E-04,-.1028E-04,-.2055E-04,-.9896E-05,-.4567E-05,   &
      -.2284E-05,-.5100E-05,-.4339E-06,-.9515E-06,-.1252E-04,-.7612E-06,   &
      -.2246E-04,-.1370E-04,-.1066E-04,-.1598E-04,-.8754E-05,-.5328E-05,   &
      -.6622E-05,-.5138E-05,-.8754E-07,-.9515E-06, .6090E-05, .4187E-05,   &
      -.3463E-04,-.1599E-04,-.1218E-04,-.2093E-04,-.9515E-05,-.4567E-05,   &
      -.1104E-05,-.1903E-05,-.1488E-05,-.3730E-05,-.4567E-05, .3045E-05,   &
      -.3463E-04,-.1675E-04,-.1294E-04,-.1979E-04,-.1066E-04,-.4187E-05,   &
      -.4034E-05,-.2893E-05,-.2588E-05,-.9401E-05, .2284E-05, .3045E-05,   &
      -.2778E-04,-.1522E-04,-.1560E-04,-.1751E-04,-.1256E-04,-.5709E-05,   &
      -.2474E-05,-.2577E-05,-.2284E-05,-.4187E-06, .7650E-05,-.3425E-05,   &
      -.3083E-04,-.1827E-04,-.1370E-04,-.1751E-04,-.1104E-04,-.9515E-05,   &
      -.6318E-05,-.4358E-05,-.7613E-07, .4643E-05, .4415E-05, .1028E-04/
      end
      
      block data ckd6
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
! pressures,  and  five  cumulative probabilities ( Fu,  1991 ). The
! spectral region is from 2850 to 2500 cm**-1.
! *********************************************************************
      common /band6/ hk(5), coeh2o(3,11,5)
      data hk / .3, .2, .2, .2, .1 /
!   .173978E+01    .115985E+01    .115985E+01    .115985E+01    .579927E+00
      data ( ( ( coeh2o(k,j,i), i = 1, 5 ), j = 1, 11 ), k = 1, 3 ) /   &
      -.1905E+02,-.1602E+02,-.1472E+02,-.1307E+02,-.1024E+02,-.1823E+02,   &
      -.1555E+02,-.1427E+02,-.1266E+02,-.9938E+01,-.1749E+02,-.1508E+02,   &
      -.1381E+02,-.1225E+02,-.9641E+01,-.1684E+02,-.1462E+02,-.1337E+02,   &
      -.1185E+02,-.9367E+01,-.1630E+02,-.1417E+02,-.1294E+02,-.1145E+02,   &
      -.9123E+01,-.1578E+02,-.1373E+02,-.1251E+02,-.1108E+02,-.8881E+01,   &
      -.1517E+02,-.1327E+02,-.1209E+02,-.1072E+02,-.8653E+01,-.1463E+02,   &
      -.1284E+02,-.1169E+02,-.1040E+02,-.8453E+01,-.1421E+02,-.1244E+02,   &
      -.1133E+02,-.1014E+02,-.8312E+01,-.1382E+02,-.1207E+02,-.1100E+02,   &
      -.9887E+01,-.8220E+01,-.1348E+02,-.1173E+02,-.1071E+02,-.9685E+01,   &
      -.8220E+01, .1024E-01, .1842E-02, .6908E-03, .1737E-02, .3517E-02,   &
       .8394E-02, .2072E-02, .8164E-03, .1716E-02, .2805E-02, .8143E-02,   &
       .2240E-02, .9001E-03, .1570E-02, .1800E-02, .8227E-02, .2386E-02,   &
       .9420E-03, .1486E-02, .1068E-02, .8373E-02, .2533E-02, .9210E-03,   &
       .1319E-02, .9420E-03, .8394E-02, .2700E-02, .9629E-03, .1026E-02,   &
       .5233E-03, .8917E-02, .2575E-02, .8792E-03, .7536E-03, .4187E-03,   &
       .9378E-02, .2617E-02, .7955E-03, .6070E-03, .4815E-03, .9797E-02,   &
       .2638E-02, .6908E-03, .5233E-03, .6280E-03, .1009E-01, .2638E-02,   &
       .4815E-03, .2931E-03, .4815E-03, .1036E-01, .2428E-02, .3140E-03,   &
       .3977E-03, .2093E-03,-.5366E-04,-.1522E-04,-.5709E-05,-.2664E-05,   &
       .3806E-05,-.4301E-04,-.1484E-04,-.4948E-05,-.7610E-06, .7610E-06,   &
      -.3920E-04,-.1484E-04,-.4948E-05, .3804E-06,-.3806E-05,-.3920E-04,   &
      -.1522E-04,-.4948E-05, .3425E-05, .1903E-05,-.3806E-04,-.1484E-04,   &
      -.3045E-05, .2664E-05, .7993E-05,-.4148E-04,-.1408E-04,-.3806E-05,   &
       .4187E-05, .7993E-05,-.5481E-04,-.1180E-04,-.3045E-05, .3045E-05,   &
       .2284E-05,-.5709E-04,-.1104E-04,-.2283E-05,-.2664E-05,-.1142E-05,   &
      -.6090E-04,-.1218E-04,-.2664E-05, .3804E-06, .3045E-05,-.6698E-04,   &
      -.1218E-04,-.2664E-05, .1523E-05,-.1142E-05,-.6508E-04,-.1218E-04,   &
      -.3425E-05, .1903E-05, .7612E-06 /
      end

      block data ckd7
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  two  cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 2200 to 1900 cm**-1.
! *********************************************************************
      common /band7/ hk(2), coeh2o(3,19,2)
      data hk / 0.7, 0.3 /
      data ( ( ( coeh2o(k,j,i), i = 1, 2 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2008E+02,-.1467E+02,-.2004E+02,-.1426E+02,-.2001E+02,-.1386E+02,   &
      -.1998E+02,-.1345E+02,-.1995E+02,-.1304E+02,-.1992E+02,-.1263E+02,   &
      -.1989E+02,-.1223E+02,-.1986E+02,-.1183E+02,-.1984E+02,-.1143E+02,   &
      -.1758E+02,-.1038E+02,-.1602E+02,-.9480E+01,-.1469E+02,-.8752E+01,   &
      -.1349E+02,-.8218E+01,-.1255E+02,-.7677E+01,-.1174E+02,-.7184E+01,   &
      -.1110E+02,-.6735E+01,-.1056E+02,-.6332E+01,-.1019E+02,-.5975E+01,   &
      -.9874E+01,-.5644E+01, .2533E-02, .2269E-01, .2575E-02, .2263E-01,   &
       .2554E-02, .2267E-01, .2491E-02, .2250E-01, .2449E-02, .2244E-01,   &
       .2344E-02, .2234E-01, .2219E-02, .2208E-01, .5694E-02, .2190E-01,   &
       .9650E-02, .2162E-01, .3286E-01, .1848E-01, .2987E-01, .1578E-01,   &
       .2527E-01, .1465E-01, .2175E-01, .1386E-01, .2056E-01, .1235E-01,   &
       .1963E-01, .1116E-01, .1926E-01, .1040E-01, .2014E-01, .1040E-01,   &
       .2024E-01, .1042E-01, .1972E-01, .1080E-01,-.8754E-05,-.6698E-04,   &
      -.1104E-04,-.6432E-04,-.1142E-04,-.6051E-04,-.1180E-04,-.6128E-04,   &
      -.1180E-04,-.6242E-04,-.1218E-04,-.6280E-04,-.1218E-04,-.6204E-04,   &
       .5328E-04,-.5709E-04, .1275E-03,-.5214E-04,-.1370E-03,-.4148E-04,   &
      -.1100E-03,-.3045E-04,-.9248E-04,-.3197E-04,-.7346E-04,-.2436E-04,   &
      -.5100E-04,-.2131E-04,-.5861E-04,-.2550E-04,-.5328E-04,-.3311E-04,   &
      -.6090E-04,-.4225E-04,-.5443E-04,-.4415E-04,-.4034E-04,-.4339E-04/
      end

      block data ckd8
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  three cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 1900 to 1700 cm**-1.
! *********************************************************************
      common /band8/ hk(3), coeh2o(3,19,3)
       data hk / 0.2, 0.7, 0.1 /
      data ( ( ( coeh2o(k,j,i), i = 1, 3 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2283E+02,-.1639E+02,-.6155E+01,-.2237E+02,-.1595E+02,-.5775E+01,   &
      -.2191E+02,-.1551E+02,-.5381E+01,-.2145E+02,-.1507E+02,-.5004E+01,   &
      -.2099E+02,-.1463E+02,-.4617E+01,-.2053E+02,-.1419E+02,-.4218E+01,   &
      -.2025E+02,-.1375E+02,-.3806E+01,-.2021E+02,-.1330E+02,-.3403E+01,   &
      -.2018E+02,-.1287E+02,-.2993E+01,-.1998E+02,-.1091E+02,-.2586E+01,   &
      -.1744E+02,-.9171E+01,-.2162E+01,-.1490E+02,-.7642E+01,-.1763E+01,   &
      -.1303E+02,-.6526E+01,-.1373E+01,-.1113E+02,-.5846E+01,-.9699E+00,   &
      -.9814E+01,-.5280E+01,-.5955E+00,-.8582E+01,-.4787E+01,-.2510E+00,   &
      -.8020E+01,-.4350E+01, .2770E-01,-.7571E+01,-.3942E+01, .2406E+00,   &
      -.7140E+01,-.3537E+01, .3567E+00, .3722E-01, .1505E-01, .6615E-02,   &
       .3722E-01, .1518E-01, .5840E-02, .3720E-01, .1526E-01, .5170E-02,   &
       .3399E-01, .1530E-01, .4773E-02, .3012E-01, .1551E-01, .4333E-02,   &
       .2625E-01, .1553E-01, .3956E-02, .2240E-01, .1562E-01, .3454E-02,   &
       .1846E-01, .1574E-01, .3161E-02, .1446E-01, .1572E-01, .3098E-02,   &
       .5924E-02, .8875E-02, .2658E-02, .2204E-01, .7096E-02, .2504E-02,   &
       .1591E-01, .5233E-02, .2292E-02, .8855E-02, .4249E-02, .2190E-02,   &
       .5422E-02, .3496E-02, .2041E-02, .4919E-02, .3621E-02, .2200E-02,   &
       .6657E-02, .3663E-02, .2248E-02, .8645E-02, .3852E-02, .2118E-02,   &
       .8771E-02, .3873E-02, .2176E-02, .9043E-02, .3747E-02, .2079E-02,   &
      -.1568E-03,-.4681E-04, .4567E-05,-.1568E-03,-.4605E-04,-.3425E-05,   &
      -.1572E-03,-.4605E-04,-.1104E-04,-.2154E-03,-.4453E-04,-.6851E-05,   &
      -.2843E-03,-.4225E-04,-.7231E-05,-.3562E-03,-.4110E-04,-.7231E-05,   &
      -.3692E-03,-.4110E-04,-.1028E-04,-.3007E-03,-.4263E-04,-.6470E-05,   &
      -.2325E-03,-.3996E-04,-.8373E-05,-.5290E-04,-.7612E-05,-.4948E-05,   &
      -.7422E-04,-.1256E-04,-.8449E-05,-.3501E-04,-.1446E-04,-.4834E-05,   &
       .4529E-04,-.2246E-04,-.2893E-05, .6470E-05,-.1789E-04,-.7498E-05,   &
      -.4948E-05,-.1713E-04,-.8183E-05,-.5481E-04,-.1713E-04,-.1447E-04,   &
      -.4986E-04,-.1903E-04,-.1353E-04,-.5138E-04,-.1484E-04,-.1147E-04,   &
      -.5328E-04,-.1560E-04,-.6588E-05/
      end

      block data ckd9
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  four cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 1700 to 1400 cm**-1.
! *********************************************************************
      common /band9/ hk(4), coeh2o(3,19,4)
       data hk / 0.22, 0.51, 0.22, 0.05 /
      data ( ( ( coeh2o(k,j,i), i = 1, 4 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2066E+02,-.1464E+02,-.8301E+01,-.3548E+01,-.2025E+02,-.1419E+02,   &
      -.7905E+01,-.3260E+01,-.2019E+02,-.1374E+02,-.7495E+01,-.2927E+01,   &
      -.2013E+02,-.1329E+02,-.7078E+01,-.2584E+01,-.2007E+02,-.1284E+02,   &
      -.6675E+01,-.2247E+01,-.2001E+02,-.1239E+02,-.6268E+01,-.1890E+01,   &
      -.1996E+02,-.1194E+02,-.5853E+01,-.1530E+01,-.1991E+02,-.1150E+02,   &
      -.5441E+01,-.1133E+01,-.1987E+02,-.1105E+02,-.5022E+01,-.7447E+00,   &
      -.1575E+02,-.9657E+01,-.4191E+01,-.3728E+00,-.1329E+02,-.8133E+01,   &
      -.3638E+01, .1616E-01,-.1181E+02,-.6675E+01,-.3178E+01, .4083E+00,   &
      -.1036E+02,-.5655E+01,-.2731E+01, .7953E+00,-.8628E+01,-.4990E+01,   &
      -.2303E+01, .1153E+01,-.7223E+01,-.4453E+01,-.1877E+01, .1454E+01,   &
      -.6567E+01,-.3974E+01,-.1461E+01, .1663E+01,-.6077E+01,-.3551E+01,   &
      -.1071E+01, .1800E+01,-.5651E+01,-.3136E+01,-.7005E+00, .1809E+01,   &
      -.5241E+01,-.2726E+01,-.3859E+00, .1781E+01, .1315E-01, .4542E-02,   &
       .3496E-02, .4877E-02, .9650E-02, .4542E-02, .3098E-02, .3956E-02,   &
       .6154E-02, .4626E-02, .2763E-02, .3077E-02, .2658E-02, .4626E-02,   &
       .2512E-02, .2261E-02, .2658E-02, .4689E-02, .2219E-02, .1405E-02,   &
       .2700E-02, .4752E-02, .1926E-02, .7473E-03, .2658E-02, .4773E-02,   &
       .1737E-02, .5066E-03, .4668E-02, .4815E-02, .1507E-02, .1842E-03,   &
       .8541E-02, .4794E-02, .1382E-02,-.2156E-03, .1022E-01, .2198E-02,   &
       .3977E-03,-.2910E-03, .5484E-02, .6698E-03, .0000E+00,-.2339E-03,   &
       .3349E-02, .1068E-02,-.2512E-03,-.4228E-03, .1884E-02, .2093E-03,   &
      -.3977E-03,-.6405E-03,-.8373E-04,-.5233E-03,-.4124E-03,-.5945E-03,   &
       .7536E-03,-.6698E-03,-.4919E-03,-.4794E-03, .3600E-02,-.4605E-03,   &
      -.4375E-03,-.3517E-03, .3873E-02,-.5861E-03,-.3203E-03,-.4689E-03,   &
       .3935E-02,-.7326E-03,-.2072E-03,-.4228E-03, .4124E-02,-.8582E-03,   &
      -.4187E-04,-.5945E-03,-.8525E-04, .1865E-04,-.1142E-05, .2664E-05,   &
      -.1313E-03, .1865E-04, .0000E+00, .1256E-04,-.6470E-04, .1865E-04,   &
      -.3045E-05, .8754E-05, .3805E-06, .1789E-04,-.6851E-05, .5328E-05,   &
       .1142E-05, .1827E-04,-.6090E-05, .4148E-05, .1142E-05, .1865E-04,   &
      -.3806E-05,-.3768E-05,-.1903E-05, .1751E-04,-.4948E-05, .3121E-05,   &
       .3159E-04, .1979E-04,-.3045E-05,-.9896E-06, .1005E-03, .1789E-04,   &
      -.6089E-05,-.1865E-05,-.2207E-04, .1941E-04, .1903E-05, .2322E-05,   &
      -.1675E-04, .6090E-05,-.7611E-06, .4397E-05, .3425E-04, .3806E-06,   &
       .1522E-05, .3806E-05, .4796E-04, .1522E-05,-.3806E-06, .3654E-05,   &
      -.6851E-05, .2664E-05,-.3920E-05,-.6850E-06,-.1370E-04, .5328E-05,   &
      -.6584E-05,-.8716E-05,-.8374E-10, .1522E-05,-.6356E-05, .1294E-05,   &
      -.9515E-05, .7612E-06,-.3235E-05,-.1066E-05,-.7612E-05, .1142E-05,   &
      -.4529E-05, .3730E-05,-.2664E-05,-.3806E-06,-.3501E-05,-.5328E-06/
      end

      block data ckd10
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  four cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 1400 to 1250 cm**-1. coech4 and coen2o
! are the coefficients to calculate the CH4 and N2O absorption coe-
! fficients in units of (cm-atm)**-1 at three temperature, nineteen
! pressures, and one cumulative probability (Fu, 1991), respectively.
! *********************************************************************
      common /band10/hk(4), coeh2o(3,19,4), coech4(3,19), coen2o(3,19)
       data hk / 0.28, 0.42, 0.25, 0.05 /
      data ( ( ( coeh2o(k,j,i), i = 1, 4 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2023E+02,-.1641E+02,-.1171E+02,-.6090E+01,-.2016E+02,-.1595E+02,   &
      -.1133E+02,-.5867E+01,-.2011E+02,-.1550E+02,-.1095E+02,-.5660E+01,   &
      -.2005E+02,-.1504E+02,-.1055E+02,-.5407E+01,-.2001E+02,-.1459E+02,   &
      -.1015E+02,-.5137E+01,-.1997E+02,-.1413E+02,-.9749E+01,-.4852E+01,   &
      -.1993E+02,-.1367E+02,-.9337E+01,-.4534E+01,-.1990E+02,-.1321E+02,   &
      -.8920E+01,-.4211E+01,-.1987E+02,-.1276E+02,-.8506E+01,-.3889E+01,   &
      -.1645E+02,-.1179E+02,-.7711E+01,-.3613E+01,-.1442E+02,-.1081E+02,   &
      -.6942E+01,-.3316E+01,-.1308E+02,-.9950E+01,-.6344E+01,-.2950E+01,   &
      -.1212E+02,-.9217E+01,-.5904E+01,-.2577E+01,-.1131E+02,-.8559E+01,   &
      -.5519E+01,-.2256E+01,-.1064E+02,-.7962E+01,-.5183E+01,-.1929E+01,   &
      -.1013E+02,-.7447E+01,-.4833E+01,-.1643E+01,-.9712E+01,-.7071E+01,   &
      -.4485E+01,-.1410E+01,-.9305E+01,-.6760E+01,-.4145E+01,-.1249E+01,   &
      -.8966E+01,-.6477E+01,-.3820E+01,-.1114E+01, .7913E-02, .8206E-02,   &
       .1509E-01, .1869E-01, .4228E-02, .8247E-02, .1467E-01, .1783E-01,   &
       .2010E-02, .8227E-02, .1442E-01, .1687E-01, .1947E-02, .8289E-02,   &
       .1394E-01, .1568E-01, .1863E-02, .8289E-02, .1346E-01, .1484E-01,   &
       .1842E-02, .8415E-02, .1310E-01, .1400E-01, .1800E-02, .8457E-02,   &
       .1275E-01, .1377E-01, .1696E-02, .8478E-02, .1220E-01, .1321E-01,   &
       .1842E-02, .8478E-02, .1189E-01, .1250E-01, .1409E-01, .8624E-02,   &
       .1254E-01, .1214E-01, .9043E-02, .1045E-01, .1225E-01, .1260E-01,   &
       .8561E-02, .1202E-01, .1181E-01, .1296E-01, .1114E-01, .1235E-01,   &
       .1191E-01, .1330E-01, .1199E-01, .1271E-01, .1195E-01, .1371E-01,   &
       .1415E-01, .1315E-01, .1218E-01, .1361E-01, .1478E-01, .1338E-01,   &
       .1296E-01, .1306E-01, .1518E-01, .1375E-01, .1365E-01, .1334E-01,   &
       .1530E-01, .1411E-01, .1392E-01, .1327E-01, .1547E-01, .1507E-01,   &
       .1390E-01, .1264E-01,-.1089E-03,-.2740E-04,-.2017E-04,-.5519E-04,   &
      -.4491E-04,-.2740E-04,-.1408E-04,-.5937E-04,-.6090E-05,-.2702E-04,   &
      -.6470E-05,-.4719E-04,-.7232E-05,-.2740E-04,-.6089E-05,-.4910E-04,   &
      -.7231E-05,-.2969E-04,-.4186E-05,-.5366E-04,-.6090E-05,-.3045E-04,   &
      -.2284E-05,-.4986E-04,-.4568E-05,-.3121E-04,-.4948E-05,-.5100E-04,   &
      -.3426E-05,-.3007E-04,-.7993E-05,-.4910E-04, .1522E-05,-.2931E-04,   &
      -.9896E-05,-.5366E-04,-.5823E-04,-.1599E-04,-.1713E-04,-.4110E-04,   &
      -.3121E-04,-.1713E-04,-.3159E-04,-.3578E-04,-.3996E-04,-.1598E-04,   &
      -.3958E-04,-.4605E-04,-.3349E-04,-.1751E-04,-.3844E-04,-.5576E-04,   &
      -.2626E-04,-.2474E-04,-.3920E-04,-.4464E-04,-.1979E-04,-.3045E-04,   &
      -.3958E-04,-.5336E-04,-.2893E-04,-.3616E-04,-.3996E-04,-.4754E-04,   &
      -.2398E-04,-.3083E-04,-.4415E-04,-.5119E-04,-.2702E-04,-.2664E-04,   &
      -.4605E-04,-.4038E-04,-.2398E-04,-.2360E-04,-.4948E-04,-.5149E-04/
      data ( ( coech4(k,j), j = 1, 19 ), k = 1, 3 ) /   &
      -.8909E+01,-.8464E+01,-.8018E+01,-.7573E+01,-.7133E+01,-.6687E+01,   &
      -.6240E+01,-.5803E+01,-.5377E+01,-.4534E+01,-.3983E+01,-.3502E+01,   &
      -.3062E+01,-.2648E+01,-.2265E+01,-.1896E+01,-.1568E+01,-.1234E+01,   &
      -.9298E+00, .9629E-03, .9838E-03, .1088E-02, .1172E-02, .1256E-02,   &
       .1402E-02, .1528E-02, .1633E-02, .1716E-02, .4815E-03,-.3977E-03,   &
      -.5652E-03,-.5024E-03,-.4605E-03,-.4563E-03,-.4438E-03,-.4521E-03,   &
      -.4312E-03,-.3789E-03,-.1294E-04,-.1408E-04,-.1522E-04,-.1675E-04,   &
      -.1751E-04,-.1941E-04,-.2246E-04,-.2207E-04,-.1827E-04,-.1256E-04,   &
      -.9515E-05,-.6470E-05,-.3045E-05,-.3806E-05,-.2055E-05,-.3730E-05,   &
      -.7612E-06,-.3806E-05, .1256E-05/
      data ( ( coen2o(k,j), j = 1, 19 ), k = 1, 3 ) /   &
      -.7863E+01,-.7412E+01,-.6963E+01,-.6514E+01,-.6065E+01,-.5611E+01,   &
      -.5167E+01,-.4720E+01,-.4283E+01,-.3454E+01,-.2858E+01,-.2404E+01,   &
      -.1922E+01,-.1491E+01,-.1097E+01,-.7177E+00,-.3548E+00, .1218E-01,   &
       .3088E+00, .4459E-02, .4542E-02, .4668E-02, .4752E-02, .4815E-02,   &
       .4919E-02, .5087E-02, .5254E-02, .5296E-02, .2324E-02, .2093E-02,   &
       .2294E-02, .2125E-02, .2058E-02, .1920E-02, .1786E-02, .1689E-02,   &
       .1788E-02, .2144E-02,-.7231E-05,-.7231E-05,-.7231E-05,-.6470E-05,   &
      -.6851E-05,-.7231E-05,-.5709E-05,-.6470E-05,-.4186E-05, .8754E-05,   &
      -.7612E-05,-.9134E-06,-.8640E-05,-.8487E-05,-.8259E-05,-.9553E-05,   &
      -.8107E-05,-.1654E-04,-.1858E-04/
      end

      block data ckd11
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and three cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 1250 to 1100 cm**-1. coech4 and coen2o
! are the coefficients to calculate the CH4 and N2O absorption coe-
! fficients in units of (cm-atm)**-1 at three temperature, nineteen
! pressures, and one cumulative probability (Fu, 1991), respectively.
! *********************************************************************
      common /band11/hk(3), coeh2o(3,19,3), coech4(3,19), coen2o(3,19)
       data hk / 0.80, 0.15, 0.05 /
      data ( ( ( coeh2o(k,j,i), i = 1, 3 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2005E+02,-.1548E+02,-.1021E+02,-.2001E+02,-.1504E+02,-.1001E+02,   &
      -.1997E+02,-.1459E+02,-.9814E+01,-.1993E+02,-.1416E+02,-.9595E+01,   &
      -.1989E+02,-.1373E+02,-.9349E+01,-.1985E+02,-.1328E+02,-.9072E+01,   &
      -.1982E+02,-.1286E+02,-.8833E+01,-.1957E+02,-.1243E+02,-.8566E+01,   &
      -.1911E+02,-.1200E+02,-.8276E+01,-.1743E+02,-.1134E+02,-.7958E+01,   &
      -.1625E+02,-.1078E+02,-.7629E+01,-.1524E+02,-.1036E+02,-.7334E+01,   &
      -.1429E+02,-.9970E+01,-.7051E+01,-.1348E+02,-.9620E+01,-.6749E+01,   &
      -.1282E+02,-.9270E+01,-.6505E+01,-.1229E+02,-.8932E+01,-.6277E+01,   &
      -.1186E+02,-.8628E+01,-.6120E+01,-.1148E+02,-.8345E+01,-.6049E+01,   &
      -.1112E+02,-.8066E+01,-.5906E+01, .1842E-02, .2131E-01, .3033E-01,   &
       .1905E-02, .2137E-01, .2841E-01, .1926E-02, .2135E-01, .2696E-01,   &
       .1926E-02, .2133E-01, .2514E-01, .1884E-02, .2154E-01, .2401E-01,   &
       .5589E-02, .2156E-01, .2321E-01, .9483E-02, .2156E-01, .2210E-01,   &
       .1333E-01, .2150E-01, .2133E-01, .1725E-01, .2154E-01, .2074E-01,   &
       .2254E-01, .1999E-01, .2005E-01, .2118E-01, .1926E-01, .1978E-01,   &
       .1936E-01, .1920E-01, .1963E-01, .1905E-01, .1911E-01, .1934E-01,   &
       .1909E-01, .1903E-01, .1920E-01, .1922E-01, .1901E-01, .1899E-01,   &
       .1934E-01, .1930E-01, .1974E-01, .1966E-01, .1909E-01, .2014E-01,   &
       .1976E-01, .1905E-01, .1984E-01, .1963E-01, .1940E-01, .1897E-01,   &
      -.1522E-05,-.6013E-04,-.5062E-04,-.2665E-05,-.6204E-04,-.5519E-04,   &
      -.3806E-05,-.6394E-04,-.5633E-04,-.4567E-05,-.6280E-04,-.5214E-04,   &
      -.6090E-05,-.6128E-04,-.5290E-04, .6051E-04,-.6242E-04,-.5823E-04,   &
       .1313E-03,-.6013E-04,-.5176E-04, .1336E-03,-.5747E-04,-.4072E-04,   &
       .6318E-04,-.5671E-04,-.3996E-04,-.5595E-04,-.3996E-04,-.4263E-04,   &
      -.3958E-04,-.4719E-04,-.4453E-04,-.3387E-04,-.5138E-04,-.5100E-04,   &
      -.5252E-04,-.4986E-04,-.4491E-04,-.5100E-04,-.4453E-04,-.4529E-04,   &
      -.5176E-04,-.4795E-04,-.4453E-04,-.5557E-04,-.5176E-04,-.5062E-04,   &
      -.5747E-04,-.4795E-04,-.5633E-04,-.5709E-04,-.4643E-04,-.3806E-04,   &
      -.5481E-04,-.5671E-04,-.4948E-04/
      data ( ( coech4(k,j), j = 1, 19 ), k = 1, 3 ) /   &
      -.1207E+02,-.1162E+02,-.1116E+02,-.1070E+02,-.1024E+02,-.9777E+01,   &
      -.9319E+01,-.8858E+01,-.8398E+01,-.7384E+01,-.6643E+01,-.6081E+01,   &
      -.5602E+01,-.5188E+01,-.4822E+01,-.4479E+01,-.4184E+01,-.3884E+01,   &
      -.3627E+01, .1036E-01, .1036E-01, .1040E-01, .1040E-01, .1045E-01,   &
       .1047E-01, .1049E-01, .1055E-01, .1059E-01, .1059E-01, .1026E-01,   &
       .1011E-01, .1024E-01, .1049E-01, .1072E-01, .1089E-01, .1109E-01,   &
       .1153E-01, .1191E-01,-.4910E-04,-.4834E-04,-.4910E-04,-.4910E-04,   &
      -.4910E-04,-.4872E-04,-.4834E-04,-.4948E-04,-.5100E-04,-.5633E-04,   &
      -.6166E-04,-.5595E-04,-.5366E-04,-.5366E-04,-.5328E-04,-.5328E-04,   &
      -.4948E-04,-.5519E-04,-.5595E-04/
      data ( ( coen2o(k,j), j = 1, 19 ), k = 1, 3 ) /   &
      -.9461E+01,-.9003E+01,-.8543E+01,-.8084E+01,-.7629E+01,-.7166E+01,   &
      -.6707E+01,-.6249E+01,-.5793E+01,-.5312E+01,-.4847E+01,-.4393E+01,   &
      -.3974E+01,-.3587E+01,-.3231E+01,-.2885E+01,-.2602E+01,-.2358E+01,   &
      -.2108E+01, .4710E-02, .4752E-02, .4773E-02, .4773E-02, .4815E-02,   &
       .4877E-02, .4898E-02, .4982E-02, .5066E-02, .5296E-02, .5149E-02,   &
       .5129E-02, .5024E-02, .4752E-02, .4501E-02, .4270E-02, .4019E-02,   &
       .3646E-02, .2759E-02,-.1484E-04,-.1408E-04,-.1446E-04,-.1446E-04,   &
      -.1522E-04,-.1560E-04,-.1522E-04,-.1522E-04,-.1598E-04,-.1484E-04,   &
      -.9895E-05,-.1028E-04,-.7612E-05,-.1903E-05, .1903E-05, .0000E+00,   &
       .2283E-05, .6166E-05,-.2740E-05/
      end

      block data ckd12
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeo3 is the coefficient to calculate the ozone absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  five cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 1100 to  980 cm**-1.    coeh2o is the
! coefficient to calculate the H2O absorption coefficient in units
! of (cm-atm)**-1 at three temperature, nineteen pressures, and one
! cumulative probability ( Fu, 1991 ).
! *********************************************************************
      common /band12/ hk(5), coeo3(3,19,5), coeh2o(3,19)
      data hk / 0.45, 0.30, 0.2, 0.04, 0.01 /
      data ( ( ( coeo3(k,j,i), i = 1, 5 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.6590E+01,-.3912E+01,-.8513E+00, .2731E+01, .5515E+01,-.6157E+01,   &
      -.3583E+01,-.7292E+00, .2740E+01, .5508E+01,-.5731E+01,-.3242E+01,   &
      -.5800E+00, .2782E+01, .5485E+01,-.5301E+01,-.2901E+01,-.4131E+00,   &
       .2805E+01, .5455E+01,-.4879E+01,-.2551E+01,-.2288E+00, .2878E+01,   &
       .5416E+01,-.4449E+01,-.2201E+01,-.2228E-01, .3000E+01, .5374E+01,   &
      -.4018E+01,-.1843E+01, .2055E+00, .3143E+01, .5342E+01,-.3615E+01,   &
      -.1502E+01, .4561E+00, .3288E+01, .5204E+01,-.3228E+01,-.1172E+01,   &
       .7099E+00, .3396E+01, .5077E+01,-.2828E+01,-.8499E+00, .9664E+00,   &
       .3463E+01, .4893E+01,-.2480E+01,-.5393E+00, .1229E+01, .3493E+01,   &
       .4656E+01,-.2181E+01,-.2653E+00, .1504E+01, .3456E+01, .4398E+01,   &
      -.1950E+01,-.1469E-01, .1735E+01, .3387E+01, .4115E+01,-.1788E+01,   &
       .2517E+00, .1919E+01, .3251E+01, .3832E+01,-.1677E+01, .5027E+00,   &
       .2032E+01, .3088E+01, .3581E+01,-.1637E+01, .7373E+00, .2100E+01,   &
       .2910E+01, .3364E+01,-.1650E+01, .9383E+00, .2123E+01, .2793E+01,   &
       .3150E+01,-.1658E+01, .1091E+01, .2112E+01, .2683E+01, .3021E+01,   &
      -.1654E+01, .1163E+01, .2099E+01, .2602E+01, .2871E+01, .9498E-02,   &
       .8894E-02, .1161E-01, .8828E-02,-.1669E-02, .9613E-02, .8347E-02,   &
       .1053E-01, .8462E-02,-.1612E-02, .9700E-02, .7829E-02, .9101E-02,   &
       .7915E-02,-.1439E-02, .9815E-02, .7167E-02, .7981E-02, .7282E-02,   &
      -.1094E-02, .9671E-02, .6764E-02, .6930E-02, .5613E-02,-.8347E-03,   &
       .9613E-02, .6312E-02, .6225E-02, .4145E-02,-.1295E-02, .9728E-02,   &
       .6099E-02, .5293E-02, .2965E-02,-.1756E-02, .9844E-02, .5915E-02,   &
       .4496E-02, .1871E-02,-.2044E-02, .9930E-02, .5817E-02, .3509E-02,   &
       .1324E-02,-.2044E-02, .9988E-02, .5535E-02, .2711E-02, .6620E-03,   &
      -.1813E-02, .1034E-01, .5247E-02, .1926E-02,-.2303E-03,-.1842E-02,   &
       .1058E-01, .4795E-02, .1197E-02,-.9498E-03,-.2216E-02, .1084E-01,   &
       .4414E-02, .6188E-03,-.1123E-02,-.2303E-02, .1079E-01, .3926E-02,   &
       .1756E-03,-.1497E-02,-.2274E-02, .1039E-01, .3425E-02,-.1900E-03,   &
      -.1353E-02,-.2389E-02, .9815E-02, .2769E-02,-.6620E-03,-.1756E-02,   &
      -.1785E-02, .9818E-02, .2444E-02,-.1016E-02,-.1410E-02,-.1698E-02,   &
       .1074E-01, .3218E-02,-.1235E-02,-.1900E-02,-.2533E-02, .1145E-01,   &
       .3684E-02,-.1364E-02,-.1353E-02,-.1957E-02,-.4030E-04,-.2375E-04,   &
      -.3814E-05,-.4943E-04,-.3166E-04,-.3742E-04,-.1871E-04,-.1137E-04,   &
      -.4317E-04,-.2878E-04,-.3526E-04,-.2015E-04,-.1295E-04,-.4821E-04,   &
      -.2303E-04,-.3382E-04,-.2087E-04,-.1519E-04,-.2231E-04,-.1871E-04,   &
      -.3454E-04,-.2087E-04,-.8109E-05,-.6476E-05,-.1511E-04,-.3454E-04,   &
      -.1820E-04,-.1269E-05,-.1439E-04,-.5037E-05,-.4173E-04,-.2598E-04,   &
       .6645E-05,-.1943E-04,-.2087E-04,-.3454E-04,-.2267E-04, .2159E-05,   &
      -.2231E-04,-.2159E-05,-.2950E-04,-.2080E-04, .2159E-06,-.4317E-05,   &
       .1799E-04,-.3670E-04,-.1590E-04,-.4461E-05,-.9354E-05,-.3598E-05,   &
      -.3216E-04,-.1475E-04,-.2231E-05,-.1295E-04,-.2878E-05,-.3576E-04,   &
      -.7347E-05,-.1022E-04,-.2159E-05,-.7915E-05,-.3015E-04,-.5230E-05,   &
      -.5109E-05,-.6476E-05,-.7196E-05,-.2331E-04,-.1079E-04,-.4102E-05,   &
       .1439E-05,-.1223E-04,-.2216E-04,-.1094E-04,-.5325E-05,-.7196E-06,   &
      -.1655E-04,-.1036E-04,-.7627E-05,-.2878E-05, .5037E-05,-.1295E-04,   &
       .1029E-04,-.1346E-04,-.4821E-05,-.7915E-05, .7915E-05, .2835E-04,   &
      -.2893E-04,-.1367E-05,-.7196E-05,-.1871E-04, .3965E-04,-.3310E-04,   &
      -.3310E-05,-.7195E-06, .2303E-04/
      data ( ( coeh2o(k,j), j = 1, 19 ), k = 1, 3 ) /   &
      -.1984E+02,-.1983E+02,-.1982E+02,-.1981E+02,-.1963E+02,-.1917E+02,   &
      -.1871E+02,-.1825E+02,-.1779E+02,-.1639E+02,-.1545E+02,-.1484E+02,   &
      -.1433E+02,-.1387E+02,-.1345E+02,-.1305E+02,-.1268E+02,-.1231E+02,   &
      -.1196E+02, .6071E-03, .2072E-02, .6196E-02, .1030E-01, .1436E-01,   &
       .1846E-01, .2259E-01, .2667E-01, .2993E-01, .2878E-01, .2803E-01,   &
       .2851E-01, .2864E-01, .2874E-01, .2862E-01, .2859E-01, .2853E-01,   &
       .2868E-01, .2887E-01,-.3808E-06, .2474E-04, .9895E-04, .1728E-03,   &
       .1911E-03, .1165E-03, .4225E-04,-.3121E-04,-.8982E-04,-.9553E-04,   &
      -.9705E-04,-.9591E-04,-.9287E-04,-.9172E-04,-.9096E-04,-.9134E-04,   &
      -.9248E-04,-.1050E-03,-.1031E-03/
      end

      block data ckd13
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  two  cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 980 to 800 cm**-1.
! *********************************************************************
      common /band13/ hk(2), coeh2o(3,19,2)
       data hk / 0.95, 0.05 /
      data ( ( ( coeh2o(k,j,i), i = 1, 2 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.1992E+02,-.1446E+02,-.1992E+02,-.1405E+02,-.1991E+02,-.1363E+02,   &
      -.1990E+02,-.1322E+02,-.1989E+02,-.1282E+02,-.1989E+02,-.1242E+02,   &
      -.1988E+02,-.1201E+02,-.1987E+02,-.1159E+02,-.1986E+02,-.1119E+02,   &
      -.1982E+02,-.1079E+02,-.1817E+02,-.1039E+02,-.1659E+02,-.1000E+02,   &
      -.1537E+02,-.9623E+01,-.1460E+02,-.9266E+01,-.1406E+02,-.8959E+01,   &
      -.1354E+02,-.8676E+01,-.1309E+02,-.8411E+01,-.1267E+02,-.8232E+01,   &
      -.1229E+02,-.8094E+01, .5024E-03, .3199E-01, .5652E-03, .3199E-01,   &
       .6071E-03, .3211E-01, .6489E-03, .3199E-01, .6699E-03, .3178E-01,   &
       .6908E-03, .3157E-01, .6908E-03, .3109E-01, .6698E-03, .3075E-01,   &
       .6698E-03, .3054E-01, .1474E-01, .3000E-01, .3085E-01, .2960E-01,   &
       .3659E-01, .2935E-01, .3016E-01, .2920E-01, .2834E-01, .2895E-01,   &
       .2780E-01, .2870E-01, .2753E-01, .2843E-01, .2755E-01, .2820E-01,   &
       .2765E-01, .2732E-01, .2769E-01, .2705E-01, .6299E-09,-.7993E-04,   &
      -.3802E-06,-.7992E-04,-.3802E-06,-.8525E-04,-.3808E-06,-.8449E-04,   &
      -.7610E-06,-.7764E-04,-.1142E-05,-.7231E-04,-.1142E-05,-.7345E-04,   &
      -.2284E-05,-.8259E-04,-.2284E-05,-.8031E-04, .2436E-03,-.7878E-04,   &
       .7612E-05,-.8525E-04,-.1248E-03,-.9439E-04,-.9477E-04,-.9172E-04,   &
      -.8982E-04,-.8640E-04,-.7916E-04,-.6813E-04,-.7574E-04,-.6090E-04,   &
      -.7612E-04,-.7117E-04,-.7498E-04,-.7041E-04,-.7269E-04,-.7992E-04/
      end

      block data ckd14
! **********************************************************************
! hk is the interval in the g (cumulative probability) space from 0
! to one. coehca and coehcb are the coefficients to calculate the
! H2O and CO2 overlapping absorption coefficients in units of (cm-
! atm)**-1 at three temperature, nineteen pressures, and ten cumu-
! lative probabilities (Fu, 1991). The spectral region is from 800
! to 670 cm**-1.
! **********************************************************************
      common /band14/ hk(10), coehca(3,19,10), coehcb(3,19,10)
       data hk / .3,.3,.2,.12,.06,.012,.004,.0025,.0011,.0004 /
      data ( ( ( coehca(k,j,i), i = 1, 10 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.1847E+02,-.1399E+02,-.1106E+02,-.8539E+01,-.5852E+01,-.3295E+01,   &
      -.1208E+01,-.6272E-01, .2055E+01, .6071E+01,-.1801E+02,-.1357E+02,   &
      -.1067E+02,-.8171E+01,-.5562E+01,-.3071E+01,-.1073E+01, .1033E+00,   &
       .2055E+01, .6071E+01,-.1755E+02,-.1314E+02,-.1027E+02,-.7798E+01,   &
      -.5224E+01,-.2823E+01,-.9280E+00, .2723E+00, .2165E+01, .5969E+01,   &
      -.1709E+02,-.1272E+02,-.9868E+01,-.7404E+01,-.4880E+01,-.2569E+01,   &
      -.6908E+00, .4453E+00, .2241E+01, .5969E+01,-.1663E+02,-.1230E+02,   &
      -.9467E+01,-.7013E+01,-.4535E+01,-.2297E+01,-.4408E+00, .6353E+00,   &
       .2359E+01, .5969E+01,-.1617E+02,-.1188E+02,-.9050E+01,-.6619E+01,   &
      -.4160E+01,-.1967E+01,-.1687E+00, .8213E+00, .2421E+01, .5969E+01,   &
      -.1571E+02,-.1147E+02,-.8629E+01,-.6230E+01,-.3771E+01,-.1648E+01,   &
       .1573E+00, .1019E+01, .2511E+01, .5884E+01,-.1525E+02,-.1106E+02,   &
      -.8215E+01,-.5841E+01,-.3393E+01,-.1331E+01, .4013E+00, .1198E+01,   &
       .2654E+01, .5794E+01,-.1480E+02,-.1066E+02,-.7800E+01,-.5454E+01,   &
      -.3032E+01,-.9870E+00, .6323E+00, .1373E+01, .2905E+01, .5647E+01,   &
      -.1402E+02,-.9693E+01,-.7206E+01,-.4846E+01,-.2656E+01,-.6540E+00,   &
       .8323E+00, .1530E+01, .3211E+01, .5355E+01,-.1343E+02,-.9060E+01,   &
      -.6596E+01,-.4399E+01,-.2294E+01,-.3519E+00, .9823E+00, .1673E+01,   &
       .3420E+01, .5083E+01,-.1279E+02,-.8611E+01,-.5785E+01,-.4010E+01,   &
      -.1936E+01,-.1177E+00, .1134E+01, .1974E+01, .3591E+01, .4770E+01,   &
      -.1230E+02,-.8174E+01,-.5298E+01,-.3611E+01,-.1607E+01, .3636E-01,   &
       .1433E+01, .2260E+01, .3539E+01, .4439E+01,-.1192E+02,-.7763E+01,   &
      -.4946E+01,-.3228E+01,-.1321E+01, .1991E+00, .1720E+01, .2420E+01,   &
       .3383E+01, .4041E+01,-.1154E+02,-.7377E+01,-.4576E+01,-.2851E+01,   &
      -.1093E+01, .4430E+00, .1896E+01, .2462E+01, .3122E+01, .3620E+01,   &
      -.1118E+02,-.7003E+01,-.4210E+01,-.2524E+01,-.8973E+00, .7490E+00,   &
       .1966E+01, .2363E+01, .2818E+01, .3182E+01,-.1080E+02,-.6677E+01,   &
      -.3872E+01,-.2264E+01,-.6846E+00, .9392E+00, .1867E+01, .2138E+01,   &
       .2505E+01, .2738E+01,-.1031E+02,-.6353E+01,-.3596E+01,-.1938E+01,   &
      -.4537E+00, .1015E+01, .1659E+01, .1830E+01, .2142E+01, .2287E+01,   &
      -.9695E+01,-.5977E+01,-.3427E+01,-.1596E+01,-.1979E+00, .9458E+00,   &
       .1363E+01, .1545E+01, .1743E+01, .1832E+01, .3628E-01, .2728E-01,   &
       .2213E-01, .1656E-01, .1507E-01, .1564E-01, .1623E-01, .1419E-01,   &
       .1455E-01, .1089E-02, .3632E-01, .2740E-01, .2164E-01, .1606E-01,   &
       .1369E-01, .1418E-01, .1444E-01, .1275E-01, .1331E-01, .9210E-03,   &
       .3636E-01, .2746E-01, .2114E-01, .1557E-01, .1239E-01, .1285E-01,   &
       .1237E-01, .1141E-01, .1141E-01, .9210E-03, .3640E-01, .2748E-01,   &
       .2064E-01, .1516E-01, .1141E-01, .1125E-01, .1092E-01, .1026E-01,   &
       .1011E-01,-.5652E-03, .3646E-01, .2746E-01, .2024E-01, .1478E-01,   &
       .1036E-01, .9688E-02, .9610E-02, .9305E-02, .9399E-02,-.6489E-03,   &
       .3651E-01, .2734E-01, .1984E-01, .1438E-01, .9436E-02, .8486E-02,   &
       .8214E-02, .8995E-02, .7892E-02,-.8582E-03, .3655E-01, .2723E-01,   &
       .1951E-01, .1402E-01, .8716E-02, .7433E-02, .7169E-02, .8072E-02,   &
       .5443E-02,-.1172E-02, .3659E-01, .2709E-01, .1911E-01, .1379E-01,   &
       .8107E-02, .6818E-02, .6818E-02, .7033E-02, .3056E-02,-.1047E-02,   &
       .3670E-01, .2698E-01, .1890E-01, .1363E-01, .7502E-02, .6371E-02,   &
       .6558E-02, .6489E-02,-.5652E-03,-.1340E-02, .3592E-01, .2238E-01,   &
       .1804E-01, .1007E-01, .6730E-02, .5512E-02, .6194E-02, .4375E-02,   &
      -.1109E-02,-.3559E-03, .3609E-01, .2242E-01, .1526E-01, .8582E-02,   &
       .6284E-02, .5809E-02, .4501E-02, .9420E-03,-.9001E-03,-.1005E-02,   &
       .3703E-01, .2196E-01, .1281E-01, .7860E-02, .5861E-02, .5842E-02,   &
       .1800E-02,-.1591E-02,-.1235E-02,-.9420E-03, .3728E-01, .2114E-01,   &
       .1347E-01, .6678E-02, .5449E-02, .4837E-02,-.1084E-02,-.1361E-02,   &
      -.6699E-03,-.1256E-03, .3683E-01, .2061E-01, .1350E-01, .6133E-02,   &
       .5449E-02, .2111E-02,-.1386E-02,-.1235E-02,-.5652E-03,-.8373E-04,   &
       .3656E-01, .1988E-01, .1348E-01, .5441E-02, .5149E-02,-.8813E-03,   &
      -.1116E-02,-.8373E-03,-.3140E-03,-.6280E-04, .3669E-01, .1934E-01,   &
       .1363E-01, .5035E-02, .3585E-02,-.1250E-02,-.9357E-03,-.8227E-03,   &
      -.3140E-03,-.4187E-04, .3618E-01, .1856E-01, .1390E-01, .3836E-02,   &
       .1470E-02,-.1096E-02,-.8080E-03,-.4480E-03,-.2093E-03,-.2093E-04,   &
       .3416E-01, .1741E-01, .1431E-01, .1951E-02,-.2923E-04,-.9422E-03,   &
      -.4576E-03,-.2395E-03,-.1565E-03,-.2799E-04, .3219E-01, .1674E-01,   &
       .1516E-01, .6652E-03,-.5051E-03,-.7052E-03,-.2002E-03,-.2135E-03,   &
      -.7633E-04,-.7300E-04,-.1290E-03,-.9934E-04,-.5595E-04,-.3996E-04,   &
       .1294E-04,-.9134E-05, .1294E-05,-.3121E-05,-.4757E-04,-.1979E-04,   &
      -.1305E-03,-.9629E-04,-.5481E-04,-.4301E-04, .1827E-04,-.9363E-05,   &
       .1777E-04,-.2185E-04,-.1903E-04,-.1675E-04,-.1313E-03,-.9439E-04,   &
      -.5404E-04,-.4263E-04, .9134E-05,-.1020E-04, .3524E-04,-.2599E-04,   &
      -.2093E-04, .1675E-04,-.1313E-03,-.9172E-04,-.5252E-04,-.4567E-04,   &
       .4186E-05,-.3920E-05, .2552E-04,-.2059E-04,-.2246E-04,-.1028E-04,   &
      -.1324E-03,-.9210E-04,-.5138E-04,-.4491E-04, .6470E-05,-.2131E-05,   &
       .1496E-04,-.1572E-04,-.3311E-04,-.8754E-05,-.1324E-03,-.9058E-04,   &
      -.5328E-04,-.4225E-04, .1827E-05,-.8411E-06, .4719E-05,-.6813E-05,   &
      -.2474E-04,-.1256E-04,-.1340E-03,-.8868E-04,-.5633E-04,-.4187E-04,   &
      -.4415E-05, .6055E-05,-.1648E-04,-.1507E-04, .1979E-04,-.2131E-04,   &
      -.1340E-03,-.8373E-04,-.5899E-04,-.3920E-04,-.4072E-05, .1491E-04,   &
      -.9781E-05,-.5328E-05, .3578E-04,-.1979E-04,-.1321E-03,-.7954E-04,   &
      -.5899E-04,-.4072E-04, .1066E-05, .5728E-05,-.5138E-05,-.8373E-05,   &
       .2626E-04,-.2436E-04,-.1363E-03,-.6432E-04,-.5176E-04,-.3083E-04,   &
       .2169E-05,-.8944E-05, .3159E-05, .6470E-05,-.4187E-05, .4948E-05,   &
      -.1302E-03,-.7802E-04,-.3311E-04,-.1903E-04, .5328E-05,-.1884E-04,   &
       .1408E-04, .3311E-04, .1142E-05,-.7613E-06,-.1473E-03,-.6737E-04,   &
      -.7536E-04,-.1085E-04,-.1903E-05,-.1458E-04, .4034E-04,-.3941E-10,   &
      -.7992E-05, .2664E-05,-.1361E-03,-.5709E-04,-.8550E-04,-.5709E-05,   &
      -.8640E-05, .6523E-05, .1903E-05,-.8221E-05,-.3045E-05,-.9134E-05,   &
      -.1329E-03,-.5529E-04,-.7107E-04, .2664E-05,-.9020E-05, .3320E-04,   &
      -.2131E-05,-.4187E-05,-.7231E-05,-.3806E-05,-.1278E-03,-.5247E-04,   &
      -.6465E-04, .3806E-05,-.6091E-05, .1245E-04,-.3844E-05,-.6090E-05,   &
      -.8754E-05,-.2664E-05,-.1321E-03,-.5632E-04,-.5897E-04, .1012E-04,   &
       .1168E-04,-.4196E-06,-.8411E-05,-.8868E-05,-.1484E-04,-.1522E-05,   &
      -.1252E-03,-.4907E-04,-.5932E-04, .3245E-04, .1996E-04,-.3325E-05,   &
      -.5785E-05,-.6394E-05,-.6851E-05,-.1142E-05,-.1093E-03,-.4731E-04,   &
      -.6761E-04, .1808E-04, .1754E-04,-.5079E-05,-.5809E-05,-.5649E-05,   &
      -.3988E-05,-.5849E-06,-.1151E-03,-.4965E-04,-.7163E-04, .7839E-05,   &
       .5505E-05,-.6084E-05,-.3344E-05,-.3894E-05,-.1391E-05,-.1327E-05/
      data ( ( ( coehcb(k,j,i), i = 1, 10 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.9398E+01,-.5678E+01,-.3606E+01,-.2192E+01, .2104E+01, .3044E+01,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.9094E+01,-.5422E+01,   &
      -.3448E+01,-.1650E+01, .2046E+01, .2749E+01,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.8760E+01,-.5270E+01,-.3329E+01,-.1147E+01,   &
       .2112E+01, .2709E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.8537E+01,-.5152E+01,-.3129E+01,-.9544E+00, .2254E+01, .2771E+01,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.8176E+01,-.4936E+01,   &
      -.2680E+01,-.9259E+00, .2247E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.7836E+01,-.4676E+01,-.2378E+01,-.3550E+00,   &
       .1396E+01, .1976E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.7419E+01,-.4122E+01,-.2407E+01,-.1204E-01, .1744E+01,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.7124E+01,-.3727E+01,   &
      -.2160E+01, .6158E+00, .1953E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.6823E+01,-.3324E+01,-.1748E+01,-.9806E-01,   &
       .2319E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.5957E+01,-.3017E+01,-.1647E+01, .1398E+01,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.5115E+01,-.2290E+01,   &
      -.5273E+00, .5662E+00, .1459E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.4162E+01,-.1453E+01, .1116E+00,-.4587E+02,   &
       .9569E+00,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.3611E+01,-.9744E+00,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.3075E+01,-.4176E+00,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.3469E+01,-.9395E+00, .5092E+00, .6200E+00,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.3808E+01,-.1505E+01, .3901E+00, .6264E+00,-.1155E+01,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4058E+01,-.1818E+01,   &
       .2693E+00, .7087E+00, .3820E+00,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.4587E+02,-.4587E+02,-.4262E+01,-.2097E+01,-.5711E-01, .5681E+00,   &
       .1310E+01, .7371E+00,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.3997E+01,-.1784E+01, .4388E-01, .5167E+00, .6930E+00,-.6906E+00,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02, .2944E-01, .2723E-01,   &
       .1854E-01, .2023E-01, .2254E-01, .3059E-02, .4788E+00, .3059E-02,   &
       .3059E-02, .3059E-02, .3080E-01, .2549E-01, .1547E-01, .2225E-01,   &
       .2107E-01, .3059E-02, .4737E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .3269E-01, .2656E-01, .2125E-01, .2179E-01, .2162E-01, .4589E+00,   &
       .4643E+00, .3059E-02, .3059E-02, .3059E-02, .3322E-01, .2476E-01,   &
       .2075E-01, .2139E-01, .1907E-01, .4501E+00, .4441E+00, .3059E-02,   &
       .3059E-02, .3059E-02, .3387E-01, .2182E-01, .2665E-01, .1841E-01,   &
       .2506E-01, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3532E-01, .2091E-01, .1995E-01, .2067E-01, .1949E-01, .4491E+00,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3468E-01, .2075E-01,   &
       .2587E-01, .1401E-01, .8646E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3059E-02, .3059E-02, .3666E-01, .2430E-01, .1919E-01, .2007E-01,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3613E-01, .2147E-01, .1892E-01, .1361E-01, .3059E-02, .4506E+00,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3129E-01, .1954E-01,   &
       .2442E-01, .1011E-01, .4420E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .3059E-02, .3059E-02, .3177E-01, .2101E-01, .1526E-01, .4376E+00,   &
       .4379E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .2887E-01, .2044E-01, .1285E-01, .3059E-02,-.4862E-03, .3059E-02,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .2759E-01, .2114E-01,   &
       .4303E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3059E-02, .3059E-02, .2880E-01, .1690E-01,-.4187E+00, .3059E-02,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .2852E-01, .2255E-01, .2184E-01, .4334E+00, .4217E+00, .3059E-02,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .2840E-01, .2136E-01,   &
       .1644E-01, .2812E-01, .4358E+00, .4288E+00, .3059E-02, .3059E-02,   &
       .3059E-02, .3059E-02, .2809E-01, .2173E-01, .1708E-01, .3346E-01,   &
       .4225E-01, .4419E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .2702E-01, .2260E-01, .1607E-01, .2720E-01, .3982E-01, .4452E+00,   &
       .4365E+00, .4345E+00, .4432E+00, .4623E+00, .2684E-01, .2328E-01,   &
       .2099E-01, .3040E-01, .3867E-01, .4389E+00, .3132E-01, .3158E-01,   &
       .4083E-01, .4580E+00,-.1581E-03,-.9707E-04,-.1250E-03, .2580E-03,   &
       .7378E-04,-.1617E-01, .8646E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.1319E-03,-.9528E-04,-.1710E-03, .7118E-04, .2076E-04,-.1608E-01,   &
       .8552E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.1721E-03,-.4680E-04,   &
      -.5522E-04,-.6242E-04, .4517E-04,-.7777E-02, .8382E-02,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.1482E-03,-.4208E-04,-.5216E-04,-.6514E-04,   &
      -.8378E-04,-.7956E-02, .8013E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.1501E-03,-.4002E-04,-.1664E-03, .2272E-04,-.1888E-03,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1201E-03,-.4709E-04,   &
      -.5371E-04,-.1574E-03, .1854E-03,-.7712E-02,-.4656E-05,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.1333E-03,-.1062E-03, .5785E-04,-.4150E-04,   &
      -.5717E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.1212E-03,-.8524E-04,-.5895E-04,-.2884E-03,-.1581E-01,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.8148E-04,-.9361E-04,   &
      -.2873E-03, .1883E-03,-.1594E-01, .8133E-02,-.4656E-05,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.1221E-03,-.1430E-04, .6335E-04,-.2581E-03,   &
       .7977E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.9257E-04,-.5008E-04, .6389E-04,-.7455E-02,-.7745E-02,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1186E-03,-.9037E-04,   &
      -.7461E-04,-.4656E-05, .1168E-03,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.8513E-04,-.5708E-04, .7763E-02,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.1124E-03,-.1228E-03, .7663E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1015E-03,-.8369E-04,   &
      -.2167E-03,-.7548E-02, .7608E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.4656E-05,-.4656E-05,-.1049E-03,-.6414E-04,-.1384E-03,-.1644E-03,   &
      -.6919E-02, .7736E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.1008E-03,-.7047E-04,-.1276E-03,-.2445E-03,-.1860E-03, .7975E-02,   &
      -.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.9629E-04,-.1007E-03,   &
      -.1127E-03,-.1527E-03,-.3238E-03,-.7373E-02, .7877E-02, .7840E-02,   &
       .7997E-02, .8345E-02,-.8800E-04,-.1072E-03,-.1046E-03,-.1777E-03,   &
      -.2146E-03,-.7016E-02, .1516E-01, .1532E-01, .1509E-01, .8268E-02/
      end

      block data ckd15
! **********************************************************************
! hk is the interval in the g (cumulative probability) space from 0
! to one. coehca and coehcb are the coefficients to calculate the
! H2O and CO2 overlapping absorption coefficients in units of (cm-
! atm)**-1 at three temperatures, nineteen pressures, and 12 cumu-
! lative probabilities (Fu, 1991). The spectral region is from 670
! to 540 cm**-1.
! **********************************************************************
      common /band15/ hk(12), coehca(3,19,12), coehcb(3,19,12)
       data hk /.24,.36,.18,.1,.05,.02,.016,.012,.01,.006,.0039,.0021/
      data ( ( ( coehca(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 1, 2 ) /   &
      -.1921E+02,-.1363E+02,-.1080E+02,-.8392E+01,-.6776E+01,-.5696E+01,   &
      -.4572E+01,-.3752E+01,-.2382E+01,-.1110E+01, .6803E+00, .3259E+01,   &
      -.1875E+02,-.1321E+02,-.1040E+02,-.8026E+01,-.6449E+01,-.5401E+01,   &
      -.4316E+01,-.3498E+01,-.2141E+01,-.9439E+00, .8103E+00, .3314E+01,   &
      -.1829E+02,-.1278E+02,-.1000E+02,-.7646E+01,-.6089E+01,-.5085E+01,   &
      -.4047E+01,-.3217E+01,-.1872E+01,-.7106E+00, .9573E+00, .3390E+01,   &
      -.1783E+02,-.1236E+02,-.9596E+01,-.7264E+01,-.5735E+01,-.4740E+01,   &
      -.3743E+01,-.2882E+01,-.1587E+01,-.4714E+00, .1120E+01, .3425E+01,   &
      -.1737E+02,-.1195E+02,-.9193E+01,-.6877E+01,-.5371E+01,-.4404E+01,   &
      -.3405E+01,-.2574E+01,-.1298E+01,-.1747E+00, .1327E+01, .3547E+01,   &
      -.1691E+02,-.1153E+02,-.8776E+01,-.6490E+01,-.4993E+01,-.4049E+01,   &
      -.3039E+01,-.2256E+01,-.1012E+01, .1103E+00, .1530E+01, .3651E+01,   &
      -.1644E+02,-.1112E+02,-.8360E+01,-.6105E+01,-.4623E+01,-.3688E+01,   &
      -.2694E+01,-.1915E+01,-.6855E+00, .3993E+00, .1714E+01, .3950E+01,   &
      -.1598E+02,-.1073E+02,-.7943E+01,-.5723E+01,-.4236E+01,-.3314E+01,   &
      -.2338E+01,-.1596E+01,-.3583E+00, .6963E+00, .1868E+01, .4127E+01,   &
      -.1553E+02,-.1034E+02,-.7542E+01,-.5357E+01,-.3856E+01,-.2942E+01,   &
      -.1986E+01,-.1299E+01,-.5472E-01, .9443E+00, .2149E+01, .4261E+01,   &
      -.1485E+02,-.9661E+01,-.7008E+01,-.4830E+01,-.3458E+01,-.2566E+01,   &
      -.1658E+01,-.9639E+00, .2083E+00, .1182E+01, .2458E+01, .4452E+01,   &
      -.1427E+02,-.9166E+01,-.6373E+01,-.4404E+01,-.3073E+01,-.2209E+01,   &
      -.1349E+01,-.6648E+00, .4023E+00, .1452E+01, .2739E+01, .4466E+01,   &
      -.1380E+02,-.8726E+01,-.5772E+01,-.3982E+01,-.2732E+01,-.1874E+01,   &
      -.1052E+01,-.4403E+00, .5763E+00, .1792E+01, .2999E+01, .4335E+01,   &
      -.1305E+02,-.8270E+01,-.5304E+01,-.3586E+01,-.2392E+01,-.1568E+01,   &
      -.8299E+00,-.2650E+00, .8584E+00, .2062E+01, .3141E+01, .4168E+01,   &
      -.1269E+02,-.7900E+01,-.4956E+01,-.3205E+01,-.2065E+01,-.1332E+01,   &
      -.6415E+00,-.7921E-01, .1170E+01, .2269E+01, .3198E+01, .4066E+01,   &
      -.1227E+02,-.7536E+01,-.4576E+01,-.2859E+01,-.1815E+01,-.1139E+01,   &
      -.4520E+00, .2272E+00, .1371E+01, .2351E+01, .3150E+01, .3935E+01,   &
      -.1186E+02,-.7159E+01,-.4223E+01,-.2538E+01,-.1619E+01,-.9324E+00,   &
      -.1566E+00, .5151E+00, .1520E+01, .2339E+01, .3132E+01, .3880E+01,   &
      -.1120E+02,-.6777E+01,-.3919E+01,-.2330E+01,-.1387E+01,-.6737E+00,   &
       .1108E+00, .6991E+00, .1531E+01, .2163E+01, .3150E+01, .3767E+01,   &
      -.9973E+01,-.6279E+01,-.3638E+01,-.2048E+01,-.1098E+01,-.4407E+00,   &
       .3043E+00, .7797E+00, .1424E+01, .2002E+01, .3122E+01, .3611E+01,   &
      -.8483E+01,-.5607E+01,-.3357E+01,-.1744E+01,-.8884E+00,-.2264E+00,   &
       .3800E+00, .7504E+00, .1245E+01, .2032E+01, .3097E+01, .3546E+01,   &
       .3762E-01, .2372E-01, .1643E-01, .1208E-01, .1170E-01, .1164E-01,   &
       .1214E-01, .1161E-01, .1028E-01, .9185E-02, .7712E-02, .1001E-01,   &
       .3762E-01, .2382E-01, .1593E-01, .1145E-01, .1059E-01, .1049E-01,   &
       .1080E-01, .1057E-01, .8894E-02, .7807E-02, .7132E-02, .1032E-01,   &
       .3764E-01, .2386E-01, .1555E-01, .1080E-01, .9692E-02, .9231E-02,   &
       .9585E-02, .9644E-02, .7711E-02, .6443E-02, .6223E-02, .9922E-02,   &
       .3764E-01, .2395E-01, .1516E-01, .1028E-01, .8917E-02, .8415E-02,   &
       .8457E-02, .8777E-02, .6436E-02, .5428E-02, .5499E-02, .8017E-02,   &
       .3768E-01, .2399E-01, .1482E-01, .9692E-02, .8247E-02, .7640E-02,   &
       .7582E-02, .7783E-02, .5432E-02, .4482E-02, .4919E-02, .5903E-02,   &
       .3770E-01, .2401E-01, .1449E-01, .9252E-02, .7620E-02, .6678E-02,   &
       .6845E-02, .6925E-02, .4939E-02, .3471E-02, .4124E-02, .3873E-02,   &
       .3776E-01, .2395E-01, .1419E-01, .8959E-02, .7096E-02, .6184E-02,   &
       .6110E-02, .6075E-02, .4419E-02, .2891E-02, .3056E-02, .1214E-02,   &
       .3780E-01, .2391E-01, .1392E-01, .8687E-02, .6573E-02, .5733E-02,   &
       .5359E-02, .5009E-02, .4034E-02, .2755E-02, .1968E-02,-.4187E-04,   &
       .3791E-01, .2382E-01, .1373E-01, .8561E-02, .6060E-02, .5120E-02,   &
       .4618E-02, .4713E-02, .3965E-02, .2481E-02, .8164E-03,-.1088E-02,   &
       .3843E-01, .2148E-01, .1302E-01, .6384E-02, .5256E-02, .4260E-02,   &
       .4077E-02, .4181E-02, .4132E-02, .2135E-02,-.2931E-03,-.1151E-02,   &
       .3896E-01, .2081E-01, .1097E-01, .5568E-02, .4475E-02, .3795E-02,   &
       .3828E-02, .3996E-02, .3766E-02, .1193E-02,-.1089E-02,-.9420E-03,   &
       .3973E-01, .2024E-01, .9943E-02, .4815E-02, .3820E-02, .3663E-02,   &
       .3568E-02, .3881E-02, .2859E-02, .6698E-03,-.1549E-02,-.6280E-03,   &
       .3635E-01, .1963E-01, .1061E-01, .3812E-02, .3509E-02, .3429E-02,   &
       .3693E-02, .3316E-02, .1120E-02, .6552E-03,-.1193E-02,-.1109E-02,   &
       .3631E-01, .1893E-01, .1056E-01, .3172E-02, .3378E-02, .3164E-02,   &
       .2751E-02, .1722E-02, .1112E-02, .4354E-03,-.7327E-03,-.1319E-02,   &
       .3500E-01, .1828E-01, .1050E-01, .2831E-02, .2784E-02, .2564E-02,   &
       .1469E-02, .7739E-03, .1209E-02, .7913E-03,-.2512E-03,-.1758E-02,   &
       .3352E-01, .1763E-01, .1045E-01, .2401E-02, .1928E-02, .1340E-02,   &
       .3753E-03, .5794E-03, .9060E-03, .1042E-02, .1465E-03,-.2533E-02,   &
       .2880E-01, .1729E-01, .1077E-01, .1347E-02, .1194E-02,-.1191E-03,   &
       .2828E-03, .6606E-03, .9743E-03, .1002E-02, .0000E+00,-.3140E-02,   &
       .2040E-01, .1585E-01, .1165E-01, .3871E-05, .1509E-04,-.1046E-02,   &
       .2444E-03, .4359E-03, .1041E-02, .2429E-02,-.1721E-03,-.2786E-02,   &
       .1737E-01, .1560E-01, .1240E-01,-.2139E-03,-.1025E-02,-.1248E-02,   &
      -.6934E-04, .1649E-03, .4062E-03, .1554E-02,-.4179E-03,-.7795E-03/
      data ( ( ( coehca(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 3, 3 ) /   &
      -.1488E-03,-.9248E-04,-.2322E-04,-.4187E-05, .1104E-04, .9895E-05,   &
      -.2283E-05, .2512E-05,-.9058E-05, .8449E-05, .8297E-05,-.3882E-04,   &
      -.1488E-03,-.9058E-04,-.2398E-04,-.5709E-05, .1218E-04, .1180E-04,   &
       .1522E-05, .6927E-05,-.1161E-04, .1714E-04,-.4948E-06,-.3540E-04,   &
      -.1500E-03,-.8830E-04,-.2474E-04,-.8373E-05, .6470E-05, .7992E-05,   &
       .9096E-05, .6737E-05,-.1485E-04, .1873E-04,-.4948E-06,-.4491E-04,   &
      -.1500E-03,-.8601E-04,-.2664E-04,-.1028E-04, .6851E-05, .6851E-05,   &
       .1294E-04,-.2550E-05,-.1520E-04, .2310E-04, .4948E-06,-.2017E-04,   &
      -.1507E-03,-.8373E-04,-.2664E-04,-.1256E-04, .4567E-05, .1028E-04,   &
       .9210E-05,-.2131E-05,-.6995E-05, .7498E-05,-.1104E-04,-.2284E-05,   &
      -.1519E-03,-.8183E-04,-.2816E-04,-.1142E-04, .7611E-06, .7231E-05,   &
       .1751E-05,-.7612E-06, .8312E-05, .2436E-05,-.7231E-05, .2398E-04,   &
      -.1530E-03,-.7992E-04,-.2893E-04,-.9896E-05, .3806E-06, .8906E-05,   &
       .3159E-05,-.5328E-05, .3692E-05,-.2093E-05,-.6851E-05,-.3045E-05,   &
      -.1538E-03,-.7536E-04,-.3007E-04,-.8754E-05,-.3045E-05, .5138E-05,   &
       .9134E-06,-.1979E-06, .1560E-05,-.1507E-04, .2284E-04, .9895E-05,   &
      -.1541E-03,-.7688E-04,-.2969E-04,-.5709E-05,-.3996E-05, .1142E-05,   &
      -.8373E-06, .1235E-04,-.7079E-05,-.6737E-05, .1028E-04, .3578E-04,   &
      -.1560E-03,-.6851E-04,-.1903E-04,-.4187E-05,-.4605E-05,-.1142E-06,   &
       .3878E-05, .3597E-05,-.9591E-05, .5328E-05, .7612E-05,-.4948E-05,   &
      -.1587E-03,-.6546E-04,-.2740E-04,-.7612E-06,-.3578E-05, .1713E-05,   &
       .6064E-05,-.9781E-05, .1408E-05, .5709E-05, .8373E-05,-.1256E-04,   &
      -.1484E-03,-.5823E-04,-.4301E-04,-.1522E-05, .7498E-05,-.5328E-06,   &
      -.7855E-05,-.1599E-05, .1964E-04,-.2284E-05, .7882E-10, .5328E-05,   &
      -.1238E-03,-.5700E-04,-.5266E-04, .3286E-05, .4910E-05,-.8602E-05,   &
       .6090E-06, .8454E-05, .1256E-05,-.4072E-05,-.1903E-05, .6470E-05,   &
      -.1155E-03,-.5231E-04,-.4396E-04, .3626E-05,-.7051E-05,-.1743E-05,   &
       .9667E-05, .2064E-04,-.2778E-05,-.6546E-05,-.4948E-05, .1903E-05,   &
      -.1024E-03,-.5129E-04,-.4506E-04, .7943E-06, .3074E-06, .3243E-05,   &
       .2754E-04,-.1479E-05, .1661E-05,-.2969E-05,-.1066E-04, .7612E-06,   &
      -.8473E-04,-.5418E-04,-.4674E-04,-.3418E-05, .9460E-05, .1151E-04,   &
       .5714E-05,-.1069E-04,-.2022E-05,-.9061E-05,-.1104E-04,-.3083E-04,   &
      -.4283E-04,-.5037E-04,-.4476E-04, .1951E-04, .8922E-05, .1296E-04,   &
      -.4053E-05,-.4355E-05,-.2355E-05,-.5004E-05,-.1218E-04,-.1522E-04,   &
       .6411E-05,-.5937E-04,-.5331E-04, .1934E-04, .5284E-05, .1129E-04,   &
      -.2166E-05,-.1484E-06,-.5407E-05,-.1364E-04,-.3115E-05, .3004E-04,   &
      -.5074E-04,-.6256E-04,-.5097E-04, .2218E-04, .1228E-04,-.1160E-05,   &
      -.1105E-05, .1618E-06,-.6089E-05,-.4216E-06,-.5314E-05, .7903E-05/
      data ( ( ( coehcb(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 1, 2 ) /   &
      -.9593E+01,-.4078E+01,-.2812E+01,-.6506E+00,-.4123E+00, .2055E+01,   &
       .4097E+01, .4671E+01, .4639E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.9276E+01,-.3757E+01,-.2467E+01,-.5784E+00, .8833E-01, .2232E+01,   &
       .3826E+01, .4723E+01, .4942E+01, .5135E+01,-.4587E+02,-.4587E+02,   &
      -.8968E+01,-.3508E+01,-.2116E+01,-.1363E+00, .1662E+00, .2424E+01,   &
       .4220E+01, .4513E+01, .1375E+01, .4601E+01,-.4587E+02,-.4587E+02,   &
      -.8662E+01,-.3164E+01,-.1722E+01, .5178E-01, .7288E+00, .2411E+01,   &
       .3805E+01, .4766E+01, .4342E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.8292E+01,-.2799E+01,-.1359E+01, .3271E+00, .1650E+01, .2395E+01,   &
       .4192E+01, .4758E+01, .2470E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.7812E+01,-.2404E+01,-.1085E+01, .7167E+00, .2202E+01, .2922E+01,   &
       .4322E+01, .4591E+01, .4186E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.7441E+01,-.2066E+01,-.7142E+00, .1057E+01, .2524E+01, .2946E+01,   &
       .4220E+01, .3607E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.7191E+01,-.1745E+01,-.3487E+00, .1453E+01, .2739E+01, .3660E+01,   &
       .4114E+01, .3245E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.6895E+01,-.1326E+01,-.3500E+00, .1647E+01, .2899E+01, .4023E+01,   &
       .3361E+01, .3360E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.5876E+01,-.9573E+00, .2014E+00, .2130E+01, .3493E+01, .4088E+01,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.4429E+01,-.3417E+00, .1204E+01, .2780E+01, .3843E+01, .3099E+01,   &
      -.4587E+02, .3605E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.3122E+01, .2697E+00, .1866E+01, .3526E+01, .3569E+01, .1025E+01,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.2284E+01, .8186E+00, .2754E+01, .3206E+01, .3704E+01,-.4587E+02,   &
      -.4587E+02, .4625E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.1711E+01, .1220E+01, .3248E+01,-.4587E+02, .2565E+01, .3297E+01,   &
      -.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.1758E+01, .7970E+00, .2758E+01, .2926E+01, .2613E+01, .1974E+01,   &
      -.4587E+02, .2310E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.1737E+01, .3499E+00, .2246E+01, .2673E+01, .3308E+01, .3463E+01,   &
       .3103E+01, .2611E+01, .2178E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.1559E+01, .2215E+00, .1875E+01, .2500E+01, .3346E+01, .3585E+01,   &
       .3946E+01, .3533E+01, .3205E+01,-.4587E+02,-.4587E+02,-.4587E+02,   &
      -.1601E+01, .5060E-01, .1275E+01, .2176E+01, .3081E+01, .3649E+01,   &
       .3940E+01, .4106E+01, .4112E+01, .4349E+01, .2292E+01,-.4587E+02,   &
      -.1222E+01, .3199E+00, .1642E+01, .2380E+01, .3254E+01, .3534E+01,   &
       .3687E+01, .3717E+01, .3402E+01, .3868E+01,-.4587E+02,-.4587E+02,   &
       .2967E-01, .1697E-01, .1795E-01, .1387E-01, .2032E-01, .1187E-01,   &
       .2560E-01, .1044E-01,-.4560E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .2998E-01, .1586E-01, .1786E-01, .1521E-01, .1710E-01, .1061E-01,   &
       .2030E-01, .1158E-01, .4452E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .2993E-01, .1551E-01, .1481E-01, .9846E-02, .2443E-01, .1150E-01,   &
       .1865E-01, .1376E-01, .4617E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .3035E-01, .1417E-01, .1438E-01, .1511E-01, .1901E-01, .8582E-02,   &
       .1746E-01, .1450E-01, .4523E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .2970E-01, .1347E-01, .1322E-01, .1252E-01, .1665E-01, .1037E-01,   &
       .1320E-01, .1199E-01, .4436E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .2949E-01, .1291E-01, .1671E-01, .1111E-01, .1400E-01, .1318E-01,   &
       .1060E-01, .1046E-01, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3004E-01, .1300E-01, .1413E-01, .9085E-02, .9764E-02, .2260E-01,   &
       .9778E-02, .4671E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3086E-01, .1436E-01, .1205E-01, .1081E-01, .4681E-02, .1479E-01,   &
       .1888E-01, .3494E-01, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3094E-01, .1500E-01, .1457E-01, .1060E-01, .8319E-02, .8983E-02,   &
       .3791E-01, .2232E-01, .4631E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .3158E-01, .1585E-01, .1292E-01, .6531E-02, .1383E-01, .4605E+00,   &
       .4662E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .3182E-01, .1586E-01, .8724E-02, .5798E-02, .2454E-01, .4607E+00,   &
       .4560E+00, .4511E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .2369E-01, .1606E-01, .5477E-02, .1228E-01, .4579E+00, .4561E+00,   &
       .4497E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .2190E-01, .1779E-01, .6267E-02, .4535E+00, .4533E+00, .3059E-02,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .2100E-01, .1653E-01, .7449E-02, .4543E+00, .4472E+00, .4439E+00,   &
       .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,   &
       .1864E-01, .1771E-01, .7040E-02, .2877E-01, .3381E-01, .2691E-01,   &
       .4466E+00, .3059E-02, .4613E+00, .3059E-02, .3059E-02, .3059E-02,   &
       .1637E-01, .1641E-01, .8424E-02, .1318E-01, .2060E-01, .3426E-01,   &
       .4122E-01, .4621E+00, .4555E+00, .4525E+00, .3059E-02, .3059E-02,   &
       .1607E-01, .1452E-01, .8013E-02, .1213E-01, .1482E-01, .2125E-01,   &
       .3379E-01, .3562E-01, .4619E+00, .4569E+00, .3059E-02, .3059E-02,   &
       .1698E-01, .1538E-01, .6616E-02, .1147E-01, .1217E-01, .1696E-01,   &
       .1871E-01, .2273E-01, .4513E-01, .4702E+00, .4617E+00, .4553E+00,   &
       .1700E-01, .1547E-01, .6456E-02, .1324E-01, .1502E-01, .2095E-01,   &
       .2547E-01, .2823E-01, .4107E-01, .4676E+00, .4583E+00, .4498E+00/
      data ( ( ( coehcb(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 3, 3 ) /   &
      -.6747E-05,-.2483E-04, .6575E-04, .1026E-03, .3888E-03,-.8519E-04,   &
      -.1629E-03,-.1808E-04,-.8355E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.2270E-04,-.3427E-04, .5118E-04, .1218E-03, .1245E-03,-.1245E-03,   &
       .3841E-05,-.4151E-04,-.8763E-02,-.1687E-01,-.4656E-05,-.4656E-05,   &
      -.4557E-04,-.3023E-04, .2286E-04, .5656E-04, .4113E-04,-.1407E-03,   &
      -.1301E-03, .8503E-04,-.7284E-02,-.1669E-01,-.4656E-05,-.4656E-05,   &
      -.5325E-04,-.5309E-04,-.1246E-04, .2244E-04, .5136E-04,-.1272E-03,   &
       .4217E-04,-.1749E-04,-.8435E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.6857E-04,-.7217E-04, .1740E-05, .3653E-04,-.1490E-03,-.4090E-04,   &
      -.2376E-04, .2047E-04,-.7974E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.1232E-03,-.9826E-04,-.2849E-04, .1703E-04,-.1895E-03,-.3363E-03,   &
       .7102E-04,-.1838E-05,-.1655E-01,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.9896E-04,-.5127E-04,-.2704E-04,-.1218E-04,-.1207E-03,-.5883E-04,   &
       .6893E-04,-.7924E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.7837E-04,-.4980E-04, .6902E-05,-.1072E-03,-.4051E-04,-.1991E-05,   &
      -.1173E-03,-.5195E-04,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.8136E-04,-.8102E-04, .1254E-03,-.4658E-04, .3173E-04,-.4461E-05,   &
      -.1558E-03,-.2036E-03, .8360E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.2232E-04,-.6411E-04, .9486E-04,-.2322E-03,-.8282E-04,-.8202E-02,   &
       .8416E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.1398E-03,-.7165E-04,-.4258E-04,-.3970E-04,-.2839E-03,-.7873E-02,   &
       .8231E-02,-.8213E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.6754E-04,-.7469E-04,-.6898E-04,-.1702E-03,-.8079E-02,-.7270E-02,   &
       .8116E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.2396E-04,-.2361E-04,-.8664E-04,-.8038E-02,-.8207E-02,-.4656E-05,   &
      -.4656E-05,-.1670E-01,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.5479E-04,-.7593E-04,-.1005E-03, .8199E-02,-.7942E-02,-.8244E-02,   &
      -.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.3806E-04,-.5825E-04,-.1003E-03,-.2925E-03,-.1506E-03, .3148E-04,   &
       .8060E-02,-.1593E-01, .8327E-02,-.4656E-05,-.4656E-05,-.4656E-05,   &
      -.4706E-04,-.3630E-04,-.7811E-04,-.6881E-04,-.1822E-03,-.3091E-03,   &
      -.3033E-03,-.7684E-02,-.7663E-02, .8167E-02,-.4656E-05,-.4656E-05,   &
      -.7669E-04,-.4610E-04,-.8063E-04,-.7250E-04,-.1094E-03,-.1241E-03,   &
      -.2944E-03,-.1736E-03,-.7886E-02, .8248E-02,-.4656E-05,-.4656E-05,   &
      -.7138E-04,-.4545E-04,-.3653E-04,-.6075E-04,-.4528E-04,-.1077E-03,   &
      -.1119E-03,-.1657E-03,-.4695E-03,-.8112E-02,-.7587E-02, .8217E-02,   &
      -.6812E-04,-.4558E-04,-.6739E-04,-.8861E-04,-.9386E-04,-.1334E-03,   &
      -.2007E-03,-.2179E-03,-.1650E-03,-.8001E-02, .8273E-02, .8118E-02/
      end

      block data ckd16
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  seven cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 540 to 400 cm**-1.
! *********************************************************************
      common /band16/ hk(7), coeh2o(3,19,7)
       data hk / .12, .24, .24, .20, .12, .06, .02 /
      data ( ( ( coeh2o(k,j,i), i = 1, 7 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2344E+02,-.2016E+02,-.1986E+02,-.1655E+02,-.1243E+02,-.8437E+01,   &
      -.4858E+01,-.2298E+02,-.2014E+02,-.1984E+02,-.1609E+02,-.1198E+02,   &
      -.8020E+01,-.4548E+01,-.2252E+02,-.2012E+02,-.1981E+02,-.1564E+02,   &
      -.1153E+02,-.7596E+01,-.4239E+01,-.2206E+02,-.2009E+02,-.1957E+02,   &
      -.1517E+02,-.1111E+02,-.7161E+01,-.3871E+01,-.2160E+02,-.2007E+02,   &
      -.1911E+02,-.1472E+02,-.1065E+02,-.6721E+01,-.3479E+01,-.2113E+02,   &
      -.2005E+02,-.1865E+02,-.1426E+02,-.1021E+02,-.6302E+01,-.3081E+01,   &
      -.2067E+02,-.2003E+02,-.1819E+02,-.1379E+02,-.9765E+01,-.5883E+01,   &
      -.2678E+01,-.2026E+02,-.2001E+02,-.1773E+02,-.1333E+02,-.9332E+01,   &
      -.5443E+01,-.2253E+01,-.2024E+02,-.1999E+02,-.1727E+02,-.1288E+02,   &
      -.8897E+01,-.5029E+01,-.1858E+01,-.2026E+02,-.1959E+02,-.1481E+02,   &
      -.1147E+02,-.7477E+01,-.4555E+01,-.1464E+01,-.2022E+02,-.1632E+02,   &
      -.1305E+02,-.9885E+01,-.6689E+01,-.4108E+01,-.1068E+01,-.1936E+02,   &
      -.1438E+02,-.1163E+02,-.8499E+01,-.6146E+01,-.3673E+01,-.6816E+00,   &
      -.1675E+02,-.1281E+02,-.1020E+02,-.7716E+01,-.5678E+01,-.3256E+01,   &
      -.3125E+00,-.1510E+02,-.1124E+02,-.8821E+01,-.7140E+01,-.5243E+01,   &
      -.2851E+01,-.2560E-01,-.1334E+02,-.9708E+01,-.8061E+01,-.6611E+01,   &
      -.4842E+01,-.2459E+01, .1711E+00,-.1155E+02,-.8798E+01,-.7440E+01,   &
      -.6123E+01,-.4439E+01,-.2089E+01, .2480E+00,-.1020E+02,-.8154E+01,   &
      -.6945E+01,-.5681E+01,-.4055E+01,-.1737E+01, .2390E+00,-.9464E+01,   &
      -.7677E+01,-.6512E+01,-.5284E+01,-.3707E+01,-.1453E+01, .2015E+00,   &
      -.9033E+01,-.7246E+01,-.6093E+01,-.4882E+01,-.3346E+01,-.1264E+01,   &
       .1033E+00, .4658E-01, .5840E-02, .4626E-02, .2688E-01, .2395E-01,   &
       .1804E-01, .2074E-01, .4660E-01, .1884E-02, .8561E-02, .2690E-01,   &
       .2403E-01, .1788E-01, .1934E-01, .4660E-01, .1800E-02, .1252E-01,   &
       .2694E-01, .2393E-01, .1786E-01, .1825E-01, .4660E-01, .1779E-02,   &
       .1649E-01, .2696E-01, .2397E-01, .1779E-01, .1765E-01, .4348E-01,   &
       .1758E-02, .2043E-01, .2696E-01, .2393E-01, .1748E-01, .1675E-01,   &
       .3944E-01, .1737E-02, .2445E-01, .2698E-01, .2384E-01, .1752E-01,   &
       .1549E-01, .3538E-01, .1654E-02, .2847E-01, .2702E-01, .2384E-01,   &
       .1714E-01, .1565E-01, .3127E-01, .1570E-02, .3245E-01, .2705E-01,   &
       .2374E-01, .1712E-01, .1514E-01, .2715E-01, .1444E-02, .3540E-01,   &
       .2711E-01, .2363E-01, .1702E-01, .1446E-01, .2960E-01, .1760E-01,   &
       .2977E-01, .2397E-01, .2087E-01, .1618E-01, .1445E-01, .2466E-01,   &
       .3039E-01, .2428E-01, .2217E-01, .1821E-01, .1593E-01, .1463E-01,   &
       .2640E-01, .2545E-01, .2231E-01, .2060E-01, .1773E-01, .1555E-01,   &
       .1473E-01, .3456E-01, .2135E-01, .2030E-01, .1844E-01, .1740E-01,   &
       .1559E-01, .1428E-01, .3203E-01, .2047E-01, .1809E-01, .1760E-01,   &
       .1725E-01, .1545E-01, .1541E-01, .2137E-01, .1857E-01, .1616E-01,   &
       .1698E-01, .1700E-01, .1537E-01, .1636E-01, .1338E-01, .1518E-01,   &
       .1580E-01, .1658E-01, .1710E-01, .1518E-01, .1513E-01, .1570E-01,   &
       .1614E-01, .1603E-01, .1673E-01, .1706E-01, .1497E-01, .1439E-01,   &
       .1987E-01, .1731E-01, .1601E-01, .1675E-01, .1681E-01, .1535E-01,   &
       .1425E-01, .2018E-01, .1723E-01, .1597E-01, .1691E-01, .1666E-01,   &
       .1509E-01, .1446E-01,-.2873E-03,-.8031E-04, .4225E-04,-.9287E-04,   &
      -.6013E-04,-.4339E-04,-.2474E-04,-.2862E-03,-.8372E-05, .1146E-03,   &
      -.9248E-04,-.6166E-04,-.3882E-04,-.1827E-04,-.2870E-03,-.6851E-05,   &
       .1865E-03,-.9172E-04,-.6128E-04,-.3616E-04,-.7612E-05,-.2877E-03,   &
      -.7231E-05, .1880E-03,-.9287E-04,-.5671E-04,-.4110E-04,-.1104E-04,   &
      -.3429E-03,-.7612E-05, .1149E-03,-.9287E-04,-.6356E-04,-.4529E-04,   &
      -.2436E-04,-.4187E-03,-.7992E-05, .4339E-04,-.9325E-04,-.6280E-04,   &
      -.4225E-04,-.3197E-04,-.4925E-03,-.8754E-05,-.2740E-04,-.9477E-04,   &
      -.6432E-04,-.3768E-04,-.3361E-04,-.5511E-03,-.8753E-05,-.9972E-04,   &
      -.9515E-04,-.6394E-04,-.3806E-04,-.3787E-04,-.4792E-03,-.1028E-04,   &
      -.1534E-03,-.9477E-04,-.6356E-04,-.3616E-04,-.2923E-04,-.5070E-03,   &
       .1922E-03,-.1028E-03,-.5823E-04,-.7954E-04,-.2550E-04,-.3893E-04,   &
      -.3776E-03,-.1043E-03,-.7993E-04,-.7422E-04,-.4948E-04,-.3007E-04,   &
      -.3863E-04, .8335E-04,-.5709E-04,-.6090E-04,-.7840E-04,-.3692E-04,   &
      -.3007E-04,-.4251E-04,-.6204E-04,-.4872E-04,-.3806E-04,-.4681E-04,   &
      -.3463E-04,-.3007E-04,-.4312E-04,-.1142E-04,-.5176E-04,-.5024E-04,   &
      -.3007E-04,-.3730E-04,-.3037E-04,-.3888E-04, .2550E-04,-.6508E-04,   &
      -.2512E-04,-.3083E-04,-.3197E-04,-.3041E-04,-.3750E-04, .1484E-04,   &
      -.1941E-04,-.2626E-04,-.3349E-04,-.3463E-04,-.2896E-04,-.1716E-04,   &
      -.7231E-04,-.3920E-04,-.2893E-04,-.3540E-04,-.3311E-04,-.3734E-04,   &
      -.2550E-05,-.7650E-04,-.3159E-04,-.2778E-04,-.3121E-04,-.2169E-04,   &
      -.4365E-04,-.1546E-04,-.7916E-04,-.2931E-04,-.2854E-04,-.3654E-04,   &
      -.1979E-04,-.4811E-04,-.1435E-04/
      end

      block data ckd17
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and  seven cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 400 to 280 cm**-1.
! *********************************************************************
      common /band17/ hk(7), coeh2o(3,19,7)
       data hk / .12, .26, .22, .20, .10, .085, .015 /
      data ( ( ( coeh2o(k,j,i), i = 1, 7 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2255E+02,-.2000E+02,-.1703E+02,-.1282E+02,-.9215E+01,-.5938E+01,   &
      -.2009E+01,-.2209E+02,-.1997E+02,-.1657E+02,-.1236E+02,-.8764E+01,   &
      -.5499E+01,-.1582E+01,-.2163E+02,-.1993E+02,-.1611E+02,-.1191E+02,   &
      -.8324E+01,-.5061E+01,-.1170E+01,-.2117E+02,-.1990E+02,-.1565E+02,   &
      -.1146E+02,-.7889E+01,-.4631E+01,-.7737E+00,-.2071E+02,-.1987E+02,   &
      -.1519E+02,-.1100E+02,-.7440E+01,-.4179E+01,-.3719E+00,-.2026E+02,   &
      -.1985E+02,-.1473E+02,-.1054E+02,-.6995E+01,-.3721E+01, .0000E+00,   &
      -.2024E+02,-.1982E+02,-.1426E+02,-.1009E+02,-.6549E+01,-.3284E+01,   &
       .4053E+00,-.2022E+02,-.1980E+02,-.1381E+02,-.9639E+01,-.6097E+01,   &
      -.2821E+01, .8375E+00,-.2021E+02,-.1933E+02,-.1335E+02,-.9187E+01,   &
      -.5653E+01,-.2379E+01, .1272E+01,-.2010E+02,-.1503E+02,-.1125E+02,   &
      -.7665E+01,-.4492E+01,-.1893E+01, .1642E+01,-.1747E+02,-.1278E+02,   &
      -.9547E+01,-.6120E+01,-.3756E+01,-.1443E+01, .1995E+01,-.1529E+02,   &
      -.1095E+02,-.8107E+01,-.5036E+01,-.3182E+01,-.1032E+01, .2429E+01,   &
      -.1370E+02,-.9303E+01,-.6691E+01,-.4357E+01,-.2683E+01,-.6173E+00,   &
       .2805E+01,-.1150E+02,-.7859E+01,-.5618E+01,-.3843E+01,-.2234E+01,   &
      -.2171E+00, .2973E+01,-.9590E+01,-.6537E+01,-.4886E+01,-.3355E+01,   &
      -.1805E+01, .1615E+00, .3157E+01,-.7530E+01,-.5699E+01,-.4306E+01,   &
      -.2892E+01,-.1388E+01, .5448E+00, .3155E+01,-.6758E+01,-.5112E+01,   &
      -.3809E+01,-.2464E+01,-.9947E+00, .8713E+00, .3203E+01,-.6245E+01,   &
      -.4610E+01,-.3376E+01,-.2058E+01,-.6166E+00, .1073E+01, .3109E+01,   &
      -.5777E+01,-.4175E+01,-.2963E+01,-.1671E+01,-.2556E+00, .1241E+01,   &
       .3014E+01, .4264E-01, .1968E-02, .1863E-01, .1436E-01, .1101E-01,   &
       .1055E-01, .1281E-01, .4264E-01, .1989E-02, .1861E-01, .1438E-01,   &
       .1095E-01, .1030E-01, .1211E-01, .3996E-01, .1968E-02, .1861E-01,   &
       .1434E-01, .1103E-01, .1019E-01, .1160E-01, .3600E-01, .1947E-02,   &
       .1861E-01, .1442E-01, .1086E-01, .1003E-01, .1157E-01, .3203E-01,   &
       .5756E-02, .1861E-01, .1444E-01, .1080E-01, .9922E-02, .1151E-01,   &
       .2801E-01, .9713E-02, .1859E-01, .1446E-01, .1070E-01, .9880E-02,   &
       .1066E-01, .2393E-01, .1369E-01, .1859E-01, .1451E-01, .1057E-01,   &
       .9880E-02, .1072E-01, .1987E-01, .1767E-01, .1863E-01, .1451E-01,   &
       .1040E-01, .9880E-02, .1057E-01, .1572E-01, .2169E-01, .1863E-01,   &
       .1442E-01, .1022E-01, .9742E-02, .1036E-01, .3391E-02, .1884E-01,   &
       .1566E-01, .1105E-01, .1011E-01, .1001E-01, .1017E-01, .1982E-01,   &
       .1444E-01, .1189E-01, .1030E-01, .9859E-02, .9861E-02, .1038E-01,   &
       .1748E-01, .1321E-01, .9922E-02, .1068E-01, .1013E-01, .9937E-02,   &
       .9958E-02, .1346E-01, .9943E-02, .9566E-02, .1097E-01, .9815E-02,   &
       .9964E-02, .1059E-01, .9817E-02, .7159E-02, .8687E-02, .1114E-01,   &
       .1007E-01, .1014E-01, .1058E-01, .3370E-02, .7264E-02, .9378E-02,   &
       .1112E-01, .9767E-02, .1016E-01, .1101E-01, .2993E-02, .8017E-02,   &
       .9566E-02, .1116E-01, .9738E-02, .1025E-01, .1086E-01, .8331E-02,   &
       .8771E-02, .1001E-01, .1117E-01, .9847E-02, .1076E-01, .1084E-01,   &
       .7850E-02, .9378E-02, .1001E-01, .1105E-01, .9964E-02, .1113E-01,   &
       .1168E-01, .8038E-02, .9336E-02, .9817E-02, .1096E-01, .1024E-01,   &
       .1175E-01, .1107E-01,-.2188E-03,-.2283E-05,-.8069E-04,-.4415E-04,   &
      -.2284E-04,-.4491E-04,-.4518E-04,-.2196E-03,-.2665E-05,-.8107E-04,   &
      -.4301E-04,-.2398E-04,-.4795E-04,-.4693E-04,-.2683E-03,-.3045E-05,   &
      -.8107E-04,-.4301E-04,-.2246E-04,-.4757E-04,-.4152E-04,-.3403E-03,   &
      -.4187E-05,-.8031E-04,-.3996E-04,-.1865E-04,-.4301E-04,-.4350E-04,   &
      -.4118E-03, .6584E-04,-.8107E-04,-.4034E-04,-.1903E-04,-.4643E-04,   &
      -.4834E-04,-.4803E-03, .1378E-03,-.8069E-04,-.4072E-04,-.1713E-04,   &
      -.5176E-04,-.3460E-04,-.4099E-03, .2101E-03,-.8069E-04,-.3920E-04,   &
      -.1713E-04,-.5024E-04,-.3524E-04,-.3391E-03, .2809E-03,-.7992E-04,   &
      -.3616E-04,-.2017E-04,-.5633E-04,-.4886E-04,-.2668E-03, .2078E-03,   &
      -.8069E-04,-.3768E-04,-.2131E-04,-.5580E-04,-.5454E-04,-.2207E-04,   &
      -.8601E-04,-.4643E-04,-.2436E-04,-.4148E-04,-.5458E-04,-.4579E-04,   &
      -.5138E-04,-.2893E-04,-.3273E-04,-.3882E-04,-.3920E-04,-.5035E-04,   &
      -.3170E-04,-.2169E-04,-.3007E-04,-.2740E-04,-.5328E-04,-.4491E-04,   &
      -.4403E-04,-.6383E-04, .4834E-04,-.2702E-04,-.4453E-04,-.4339E-04,   &
      -.4457E-04,-.4551E-04,-.8133E-04, .3768E-04,-.7611E-06,-.2626E-04,   &
      -.4643E-04,-.4305E-04,-.4840E-04,-.5149E-04, .7193E-04,-.2169E-04,   &
      -.4491E-04,-.3996E-04,-.4483E-04,-.4487E-04,-.6698E-04,-.4834E-04,   &
      -.3463E-04,-.4986E-04,-.4377E-04,-.4514E-04,-.5377E-04,-.2626E-04,   &
      -.4187E-04,-.3692E-04,-.5100E-04,-.4651E-04,-.4392E-04,-.5386E-04,   &
      -.4643E-04,-.4301E-04,-.3578E-04,-.5176E-04,-.4594E-04,-.4551E-04,   &
      -.3920E-04,-.3425E-04,-.4491E-04,-.3654E-04,-.5138E-04,-.4377E-04,   &
      -.5614E-04,-.5758E-04,-.3600E-04/
      end

      block data ckd18
! *********************************************************************
! hk is the interval in the g (cumulative probability) space from 0 
! to one. coeh2o is the coefficient to calculate the H2O absorption
! coefficient in units of (cm-atm)**-1 at there temperatures, nine-
! teen pressures, and eight cumulative probabilities ( Fu,  1991 ).
! The spectral region is from 280 to 0 cm**-1.
! *********************************************************************
      common /band18/ hk(8), coeh2o(3,19,8)

!       data hk / .07, .1, .2, .25, .2, .1, .03, .02 /
        data hk / .10, .1, .2, .25, .2, .1, .03, .02 /
      data ( ( ( coeh2o(k,j,i), i = 1, 8 ), j = 1, 19 ), k = 1, 3 ) /   &
      -.2121E+02,-.2002E+02,-.1676E+02,-.1274E+02,-.8780E+01,-.5167E+01,   &
      -.2692E+01,-.6275E+00,-.2075E+02,-.1996E+02,-.1630E+02,-.1228E+02,   &
      -.8324E+01,-.4718E+01,-.2260E+01,-.2303E+00,-.2029E+02,-.1990E+02,   &
      -.1584E+02,-.1182E+02,-.7868E+01,-.4269E+01,-.1806E+01, .1645E+00,   &
      -.2022E+02,-.1985E+02,-.1538E+02,-.1136E+02,-.7417E+01,-.3820E+01,   &
      -.1373E+01, .5657E+00,-.2018E+02,-.1981E+02,-.1492E+02,-.1090E+02,   &
      -.6965E+01,-.3369E+01,-.9319E+00, .9577E+00,-.2013E+02,-.1937E+02,   &
      -.1446E+02,-.1044E+02,-.6512E+01,-.2917E+01,-.4928E+00, .1376E+01,   &
      -.2009E+02,-.1891E+02,-.1400E+02,-.9984E+01,-.6063E+01,-.2466E+01,   &
      -.6887E-01, .1768E+01,-.2006E+02,-.1845E+02,-.1354E+02,-.9530E+01,   &
      -.5618E+01,-.2024E+01, .3615E+00, .2196E+01,-.2003E+02,-.1800E+02,   &
      -.1308E+02,-.9075E+01,-.5174E+01,-.1593E+01, .7820E+00, .2600E+01,   &
      -.1827E+02,-.1464E+02,-.1097E+02,-.7525E+01,-.3733E+01,-.1077E+01,   &
       .1204E+01, .3014E+01,-.1525E+02,-.1210E+02,-.9275E+01,-.5876E+01,   &
      -.2768E+01,-.6286E+00, .1622E+01, .3394E+01,-.1298E+02,-.1060E+02,   &
      -.7764E+01,-.4462E+01,-.2154E+01,-.2001E+00, .2034E+01, .3756E+01,   &
      -.1157E+02,-.8941E+01,-.5984E+01,-.3509E+01,-.1651E+01, .2279E+00,   &
       .2422E+01, .4066E+01,-.9986E+01,-.7062E+01,-.4794E+01,-.2818E+01,   &
      -.1196E+01, .6394E+00, .2791E+01, .4283E+01,-.8064E+01,-.5512E+01,   &
      -.3933E+01,-.2274E+01,-.7559E+00, .1036E+01, .3085E+01, .4444E+01,   &
      -.6440E+01,-.4863E+01,-.3219E+01,-.1791E+01,-.3279E+00, .1427E+01,   &
       .3304E+01, .4527E+01,-.5902E+01,-.4207E+01,-.2756E+01,-.1350E+01,   &
       .7686E-01, .1776E+01, .3475E+01, .4550E+01,-.5439E+01,-.3739E+01,   &
      -.2330E+01,-.9233E+00, .4612E+00, .2066E+01, .3564E+01, .4502E+01,   &
      -.5006E+01,-.3316E+01,-.1906E+01,-.5066E+00, .8352E+00, .2272E+01,   &
       .3587E+01, .4419E+01, .2338E-01, .1968E-02, .9503E-02, .3412E-02,   &
       .6280E-03,-.1109E-02,-.1089E-02,-.1026E-02, .1972E-01, .2093E-02,   &
       .9503E-02, .3391E-02, .6489E-03,-.1172E-02,-.1164E-02,-.1158E-02,   &
       .1603E-01, .3328E-02, .9524E-02, .3391E-02, .6489E-03,-.1277E-02,   &
      -.1229E-02,-.1296E-02, .1229E-01, .7138E-02, .9524E-02, .3370E-02,   &
       .6070E-03,-.1319E-02,-.1264E-02,-.1610E-02, .8478E-02, .1095E-01,   &
       .9566E-02, .3412E-02, .5652E-03,-.1382E-02,-.1266E-02,-.1566E-02,   &
       .4563E-02, .1480E-01, .9566E-02, .3412E-02, .5443E-03,-.1423E-02,   &
      -.1199E-02,-.1679E-02, .2261E-02, .1865E-01, .9608E-02, .3454E-02,   &
       .4815E-03,-.1423E-02,-.1296E-02,-.1555E-02, .2198E-02, .2250E-01,   &
       .9671E-02, .3412E-02, .4187E-03,-.1426E-02,-.1472E-02,-.1800E-02,   &
       .2072E-02, .2600E-01, .9734E-02, .3433E-02, .3977E-03,-.1428E-02,   &
      -.1541E-02,-.1591E-02, .1987E-01, .8645E-02, .6280E-02, .1298E-02,   &
      -.1151E-02,-.1509E-02,-.1662E-02,-.1570E-02, .4668E-02, .8373E-02,   &
       .3956E-02,-.4187E-04,-.1968E-02,-.1624E-02,-.1700E-02,-.1947E-02,   &
       .9231E-02, .5694E-02, .1444E-02,-.2512E-03,-.1827E-02,-.1662E-02,   &
      -.1576E-02,-.1633E-02, .8666E-02, .3077E-02,-.1737E-02,-.1277E-02,   &
      -.1507E-02,-.1757E-02,-.1612E-02,-.1612E-02, .8164E-03,-.4375E-02,   &
      -.1884E-02,-.1277E-02,-.1564E-02,-.1853E-02,-.1591E-02,-.1486E-02,   &
      -.1486E-02,-.2596E-02,-.1633E-02,-.1539E-02,-.1662E-02,-.1846E-02,   &
      -.1423E-02,-.1277E-02,-.1423E-02,-.2617E-02,-.1005E-02,-.1379E-02,   &
      -.1687E-02,-.1905E-02,-.1528E-02,-.1298E-02,-.1675E-03,-.1947E-02,   &
      -.5024E-03,-.1325E-02,-.1696E-02,-.1698E-02,-.1486E-02,-.1277E-02,   &
       .1047E-03,-.1109E-02,-.5861E-03,-.1363E-02,-.1620E-02,-.1666E-02,   &
      -.1507E-02,-.9210E-03, .1047E-03,-.1047E-02,-.8394E-03,-.1342E-02,   &
      -.1591E-02,-.1323E-02,-.1340E-02,-.9420E-03,-.1085E-03, .2283E-05,   &
      -.4719E-04,-.3807E-06,-.1522E-05,-.3425E-05,-.7612E-06, .1751E-05,   &
      -.1766E-03, .1523E-05,-.4719E-04,-.7609E-06,-.3807E-06,-.3045E-05,   &
       .1599E-05, .8723E-05,-.2443E-03, .1941E-04,-.4757E-04,-.1522E-05,   &
      -.3806E-06,-.1903E-05,-.2778E-05, .1294E-04,-.1838E-03, .8563E-04,   &
      -.4757E-04,-.1903E-05, .1142E-05,-.2664E-05,-.6090E-06, .1321E-04,   &
      -.1161E-03, .1526E-03,-.4757E-04,-.2664E-05,-.3805E-06,-.3806E-05,   &
      -.2093E-05, .2253E-04,-.4795E-04, .9248E-04,-.4757E-04,-.1903E-05,   &
       .0000E+00,-.3045E-05,-.7992E-06, .1393E-04,-.9134E-05, .2246E-04,   &
      -.4834E-04,-.2664E-05, .3804E-06,-.5328E-05,-.1510E-05, .1465E-04,   &
      -.1028E-04,-.4757E-04,-.4948E-04,-.1142E-05, .7614E-06,-.4910E-05,   &
      -.5709E-06, .1477E-04,-.1256E-04,-.1066E-03,-.4910E-04,-.1523E-05,   &
      -.3805E-06,-.3121E-05,-.2512E-05, .1142E-04,-.7878E-04,-.2664E-05,   &
      -.8373E-05,-.7612E-06, .1104E-04,-.3311E-05,-.1979E-05, .5709E-05,   &
      -.2626E-04,-.4872E-04,-.3808E-06,-.2283E-05, .2284E-05,-.3349E-05,   &
      -.4034E-05, .7231E-05,-.4910E-04, .1599E-04, .1256E-04,-.7612E-05,   &
       .1180E-05,-.1815E-05,-.7193E-05, .3045E-05, .1576E-09, .6470E-05,   &
      -.1408E-04,-.1903E-05, .1522E-05,-.4746E-05,-.4948E-05, .3806E-06,   &
       .9020E-04, .5214E-04, .6090E-05,-.1104E-04, .1180E-05,-.2778E-05,   &
      -.6090E-05,-.2664E-05,-.6737E-04,-.1218E-04,-.3806E-05,-.5214E-05,   &
      -.1066E-05,-.1294E-05,-.3045E-05,-.2664E-05,-.4643E-04, .1713E-04,   &
      -.1218E-04,-.6204E-05,-.2360E-05,-.1979E-05,-.1903E-05,-.3806E-05,   &
      -.3045E-04,-.1256E-04,-.9134E-05,-.6508E-05,-.1027E-05,-.7993E-06,   &
      -.1142E-05,-.7992E-05,-.3616E-04,-.1028E-04,-.1066E-04,-.6051E-05,   &
       .1066E-05,-.1751E-05,-.2284E-05,-.2284E-05,-.3920E-04,-.9895E-05,   &
      -.1321E-04,-.3844E-05,-.2055E-05,-.2512E-05,-.3806E-05,-.3425E-05/
      end

	block data cont_coef     
	parameter (ncoef=7,nband=12)
	real aa0,aa1,aa2,aa3,aa4

	common /pcont_aa0/ aa0(ncoef,nband)
	common /pcont_aa1/ aa1(ncoef,nband)
	common /pcont_aa2/ aa2(ncoef,nband)
	common /pcont_aa3/ aa3(ncoef,nband)
	common /pcont_aa4/ aa4(ncoef,nband)

	data aa0/& !! ORIGINAL !! pres atm
       -1.00e+01, 0.000e+00, 0.000e+00, 0.000e+00,   &
       0.000e+00, 0.000e+00, 0.000e+00,   &
       3.457e+00, 9.977e-01,-6.202e-03, 1.822e+00,   &
       8.200e+00, 8.343e-03, 1.413e-02,   &
       2.327e+00, 1.011e+00,-7.050e-03, 1.368e+00,   &
       1.508e+01,-1.747e-03, 1.271e-01,   &
       6.139e+00, 1.015e+00,-1.505e-02, 5.499e-01,   &
       6.284e+00,-7.911e-03, 5.933e-01,   &
       9.234e+00, 1.003e+00,-2.254e-02, 1.204e-01,   &
       3.237e-01,-1.697e-03, 9.105e-01,   &
       9.426e+00, 9.999e-01,-2.438e-02, 2.410e-02,   &
      -7.441e-01, 4.155e-05, 9.880e-01,   &
       8.408e+00, 1.000e+00,-2.319e-02, 2.111e-02,   &
      -8.137e-01, 6.080e-06, 9.933e-01,   &
       7.285e+00, 1.001e+00,-2.049e-02, 5.285e-02,   &
      -5.057e-01,-3.719e-04, 9.647e-01,   &
       1.096e+00, 1.008e+00,-6.852e-03, 1.021e+00,   &
       8.523e+00,-5.127e-03, 2.513e-01,   &
       -1.00e+01, 0.000e+00, 0.000e+00, 0.000e+00,   &
       0.000e+00, 0.000e+00, 0.000e+00,   &
       -1.00e+01, 0.000e+00, 0.000e+00, 0.000e+00,   &
       0.000e+00, 0.000e+00, 0.000e+00,   &
       -1.00e+01, 0.000e+00, 0.000e+00, 0.000e+00,   &
       0.000e+00, 0.000e+00, 0.000e+00   &
      /

	data aa1/& !ckdfu.fuliou.lin.out
       6.762e+00, 9.996e-01,-7.259e-03, 9.977e-01,   &
       1.603e+01, 1.709e-01,-1.394e-03,   &
       5.261e+00, 9.971e-01,-6.219e-03, 9.868e-01,   &
       1.561e+01, 1.275e-01, 2.789e-03,   &
       3.235e+00, 9.872e-01,-6.404e-03, 8.690e-01,   &
       2.107e+01, 1.222e-01, 6.908e-02,   &
       5.942e+00, 9.982e-01,-1.368e-02, 3.916e-01,   &
       9.076e+00, 2.730e-02, 5.256e-01,   &
       9.217e+00, 1.003e+00,-2.221e-02, 8.940e-02,   &
       6.127e-01,-6.341e-04, 8.978e-01,   &
       9.478e+00, 1.000e+00,-2.441e-02, 1.837e-02,   &
      -7.174e-01,-8.151e-06, 9.865e-01,   &
       8.451e+00, 1.000e+00,-2.323e-02, 1.559e-02,   &
      -8.067e-01, 5.546e-06, 9.932e-01,   &
       7.281e+00, 1.001e+00,-2.035e-02, 4.216e-02,   &
      -4.622e-01,-4.440e-04, 9.615e-01,   &
       1.081e+00, 9.983e-01,-4.259e-03, 8.882e-01,   &
       9.883e+00, 1.679e-02, 6.666e-02,   &
       3.738e+00, 9.973e-01,-3.824e-03, 9.842e-01,   &
       8.171e+00, 5.865e-02, 4.520e-03,   &
       2.879e+00, 9.953e-01,-3.841e-03, 9.721e-01,   &
       7.980e+00, 5.344e-02, 1.106e-02,   &
       7.045e-01, 1.008e+00,-5.968e-03, 7.708e-01,   &
       7.788e+00,-5.452e-03, 1.683e-01   &
      /
	data aa2/& !ckdfu.fuliou.log.out
       6.405e+00, 9.986e-01,-7.277e-03, 9.963e-01,   &
       1.304e+01, 1.588e-01,-3.569e-04,   &
       5.039e+00, 9.967e-01,-6.107e-03, 9.820e-01,   &
       1.819e+01, 1.188e-01, 4.073e-03,   &
       3.333e+00, 9.832e-01,-6.811e-03, 8.141e-01,   &
       2.107e+01, 1.416e-01, 1.102e-01,   &
       6.498e+00, 1.000e+00,-1.486e-02, 3.513e-01,   &
       7.720e+00, 2.083e-02, 5.794e-01,   &
       9.328e+00, 1.003e+00,-2.251e-02, 8.079e-02,   &
       4.519e-01,-5.421e-04, 9.084e-01,   &
       9.448e+00, 9.999e-01,-2.438e-02, 1.711e-02,   &
      -7.349e-01, 5.084e-05, 9.879e-01,   &
       8.424e+00, 1.000e+00,-2.318e-02, 1.598e-02,   &
      -8.060e-01,-2.452e-05, 9.931e-01,   &
       7.343e+00, 1.001e+00,-2.051e-02, 3.936e-02,   &
      -5.032e-01,-3.829e-04, 9.655e-01,   &
       1.962e+00, 1.004e+00,-6.678e-03, 7.064e-01,   &
       9.220e+00, 5.062e-03, 2.399e-01,   &
       3.616e+00, 9.972e-01,-3.827e-03, 9.825e-01,   &
       8.823e+00, 5.525e-02, 5.027e-03,   &
       2.481e+00, 9.947e-01,-3.944e-03, 9.568e-01,   &
       9.687e+00, 5.039e-02, 1.942e-02,   &
       2.141e+00, 1.009e+00,-1.016e-02, 5.728e-01,   &
       5.100e+00,-5.114e-03, 3.688e-01   &
      /

	data aa3/& !ckdfu.fuliou.lin.plnk.out
       6.899e+00, 9.989e-01,-7.230e-03, 9.979e-01,   &
       1.680e+01, 2.287e-01,-9.470e-04,   &
       5.300e+00, 9.971e-01,-6.528e-03, 9.864e-01,   &
       1.615e+01, 1.283e-01, 2.595e-03,   &
       3.351e+00, 9.867e-01,-6.864e-03, 8.681e-01,   &
       2.138e+01, 1.249e-01, 6.981e-02,   &
       5.981e+00, 9.983e-01,-1.384e-02, 3.944e-01,   &
       9.316e+00, 2.791e-02, 5.220e-01,   &
       9.253e+00, 1.003e+00,-2.230e-02, 9.128e-02,   &
       6.481e-01,-6.383e-04, 8.963e-01,   &
       9.555e+00, 1.000e+00,-2.460e-02, 1.921e-02,   &
      -6.885e-01,-3.569e-05, 9.855e-01,   &
       8.495e+00, 1.000e+00,-2.333e-02, 1.484e-02,   &
      -7.828e-01,-1.622e-06, 9.927e-01,   &
       7.328e+00, 1.001e+00,-2.046e-02, 3.898e-02,   &
      -5.016e-01,-3.945e-04, 9.649e-01,   &
       8.158e-01, 9.974e-01,-3.582e-03, 8.775e-01,   &
       9.817e+00, 1.744e-02, 7.603e-02,   &
       3.482e+00, 9.971e-01,-3.205e-03, 9.827e-01,   &
       8.142e+00, 5.638e-02, 5.275e-03,   &
       3.262e+00, 9.969e-01,-4.671e-03, 9.760e-01,   &
       8.702e+00, 4.119e-02, 8.187e-03,   &
       1.006e+00, 1.008e+00,-6.560e-03, 8.004e-01,   &
       8.480e+00,-5.247e-03, 1.390e-01   &
      /
	data aa4/& !ckdfu.fuliou.log.plnk.out
       6.883e+00, 1.007e+00,-7.231e-03, 1.006e+00,   &
       1.711e+01, 2.401e-01,-8.842e-04,   &
       5.128e+00, 1.006e+00,-6.531e-03, 9.905e-01,   &
       1.907e+01, 1.219e-01, 4.032e-03,   &
       3.529e+00, 9.903e-01,-7.494e-03, 8.184e-01,   &
       2.150e+01, 1.434e-01, 1.137e-01,   &
       6.623e+00, 1.008e+00,-1.525e-02, 3.557e-01,   &
       7.867e+00, 2.121e-02, 5.826e-01,   &
       9.438e+00, 1.010e+00,-2.281e-02, 8.279e-02,   &
       4.923e-01,-4.263e-04, 9.138e-01,   &
       9.574e+00, 1.005e+00,-2.472e-02, 1.761e-02,   &
      -7.078e-01, 1.032e-04, 9.917e-01,   &
       8.530e+00, 1.007e+00,-2.347e-02, 1.612e-02,   &
      -7.979e-01, 1.613e-05, 1.000e+00,   &
       7.429e+00, 1.006e+00,-2.076e-02, 3.668e-02,   &
      -5.353e-01,-2.511e-04, 9.738e-01,   &
       1.838e+00, 1.013e+00,-6.037e-03, 6.632e-01,   &
       8.190e+00, 1.126e-04, 2.812e-01,   &
       3.284e+00, 9.994e-01,-3.042e-03, 9.830e-01,   &
       8.944e+00, 4.588e-02, 5.890e-03,   &
       2.977e+00, 9.987e-01,-5.050e-03, 9.663e-01,   &
       9.910e+00, 4.875e-02, 1.525e-02,   &
       2.177e+00, 1.010e+00,-1.035e-02, 6.566e-01,   &
       6.613e+00,-5.513e-03, 2.837e-01   &
      /
      end


!===================================================================

!     ******************************************************************

function esif(t)

!     This function calculates the ice saturation vapor pressure as a
!     function of Celsius temperature

implicit none
real esif,x,t,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6114327e+03,c1= .5027041e+02,c2= .1875982e+01)
parameter (c3= .4158303e-01,c4= .5992408e-03,c5= .5743775e-05)
parameter (c6= .3566847e-07,c7= .1306802e-09,c8= .2152144e-12)

x=max(-80.,t)
esif=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end

!     ******************************************************************
function eslf(t)

!     This function calculates the liquid saturation vapor pressure as a
!     function of Celcius temperature

implicit none
real eslf,x,t,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01)
parameter (c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05)
parameter (c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13)

x=max(-80.,t)
eslf=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end

!     ******************************************************************

