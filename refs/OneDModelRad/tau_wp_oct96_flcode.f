!=========================================================================
	function tau_wp(idir,rin,radius,ipha)
! NEW (10-30-96)MINNIS VERSION OF VISIBLE OPTICAL DEPTH TO 
! WATER/ICE PATH RELATION

! IDIR= 1        TAU ----> WATER PATH
! IDIR=-1 WATER PATH ----> TAU
! RIN ( TAU if IDIR=1) ( LWP/IWP if IDIR=-1)
! RADIUS ( radius in um WATER droplets)
!         ( effective diameter in um for ICE crystals)
! IPHA  ( 1= WATER) ( 2=ICE)

	real*4 iwpcoefs(3),qextcoefs(3)
!c OLD	data iwpcoefs/0.1473,0.2833,-2.594E-4/
	data iwpcoefs/2.453e-1,1.2196e-3,-3.4745e-6/ !NEW
	data qextcoefs/2.416, -0.1854,.0209/   
! ICE	
	if(ipha.eq.2)then
! OLD	veae= iwpcoefs(1)+
!     &        iwpcoefs(2)*radius +
!     &        iwpcoefs(3)*radius*radius
	veae= iwpcoefs(1)*radius +   &
              iwpcoefs(2)*radius *radius+   &
              iwpcoefs(3)*radius*radius*radius

	val=veae
	elseif(ipha.eq.1) then
! WATER
	 alr=alog(radius)
	 qext = qextcoefs(1) + qextcoefs(2)*alr +qextcoefs(3)*alr*alr
	 val =4.*radius/3./qext
	else
! OTHER
!        print*, ipha ne 1or2
	 tau_wp=-9999.
	 return
	endif

	if(idir.eq. 1)tau_wp=rin*val !TAU ----> WATER PATH
	if(idir.eq.-1)tau_wp=rin/val !WATER PATH ----> TAU

	return
	end
	
