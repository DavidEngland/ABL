	program model

IMPLICIT NONE
!-------------------------------------------------------------------------------
!
!		ONE-DIMENSIONAL BOUNDARY LAYER MODEL (1DBLM)
!	
!-------------------------------------------------------------------------------

	INCLUDE 'BL.incl'
	INCLUDE 'ASSIM.incl'
        Real  cb,cosz,alb,em,rlongdnsfc	
	Integer IRadswitch
	REAL FLUXIN,cp,rhoa
!-------------------------
!Set atmospheric radiation on IRadswitch =1
!			   off IRadswitch=0
	
	IRadswitch =1	
	FirstCE=1
	Numheatingavg=0
	Nstepsavg=0
	Do j=1,Kz
	Radheating(j)=0
	Diffheating(j)=0
	End do
	
!---------------------------------------------------------
!!			MODEL INITIALIZATION
!-------------------------------------------------------------------------------
        write (*,*) "stepping into INITIALIZE...."
	CALL INITIALIZE(Z, T, U, V, Q, P, rhop,ozone,TIME, MTIME, TTIME, &
	                TINC,Tincflux,Tnextflux, Z0, 	&
                        ZM, DELZ, DELZM, TEMPOUT, DELT, T0, TSTEPS, FIRST, &
			Firstflux,	&
			FLUXTYPE, TIME0, STABFUNC, UG, VG, VONKAR, LAT, IMO, IDY)
        
	
	
	
	ALLOCATE(OBS_INSOL(TSTEPS))
        ALLOCATE(OBS_TRADT(TSTEPS))
        ALLOCATE(TR_MODEL(TSTEPS))
        ALLOCATE(TG_MODEL(TSTEPS))
!       CALL EKMAN(UGARR, VGARR, Z, KZ, U, V, KM, KH, LAT)
	CALL ASSIM_INIT(IRAD, ITRAD, OBS_INSOL, OBS_TRADT, TSTEPS, ADJ_FLAG)
        write (*,*) "stepping into SURFACEINIT...."
	Write(*,*) 'Before call surfaceinit  T(1)', T(1)
	CALL SURFACEINIT (USTAR, TSTAR, QSTAR, L, HO, LH, TIME, T0, TR, TAERO,	&
                         T(1), Q(1), IRAD, LAT, IMO, IDY, TSTEP, DELT, TSTEPS,  &
                         OBS_INSOL, U(1), V(1), Z(1), Z0, Q0, FLUXTYPE, 	&
			 STABFUNC, RIHALF,        &
                         XIL,cb,rlongdnsfc,Iradswitch,ColumnEnergy)

	write(*,*) "AFTER SURFACE INIT: ustar, tstar, tr, taero, t0, t1"
	write(*,*) ustar, tstar, tr, taero, t0, t(1)
	
	
	If (Iradswitch.EQ.1) then
		
		Cosz =-3.0
		
		
		Call InitRADIATION (Delt,P,T,Q,RHOP,ozone,cosz, alb,em,Tr,rlongdnsfc)
		Write (*,*) 'Initial Down Longwave', rlongdnsfc
		
	End If
	Write(*,*) 'After Init Radiation'
	Write(*,*) 'T(1)', T(1)
!Get initial pressure consistent with hydrosatic equation	
	   
	Call Pressure (P,T,Zm)
	
	Write(*,*) 'After Pressure'
	Write(*,*) 'T(1)', T(1)

	
	Call ColumnenergyCalc(FirstCE,T,p,Columnenergy,InitialColumnE,Radiativeloss,Fthrd,delt,zm)
	Write(*,*) 'Init Column Energy', Columnenergy
	
	Write(*,*) 'T(1)', T(1)
!
!-------------------------------------------------------------------------------
! 		       BEGIN MODEL INTEGRATION
!-------------------------------------------------------------------------------

	RESETT 	= 24.0*3600.
	TCOUNT  = TEMPOUT		! TEMPERATURE OUTPUT COUNTER 
	TSTEP	= 1			! ACTUAL TIME STEP NUMBER
	MTIME  	= 0			! MODEL TIME (STARTS AT 0)
	IT     	= 0
        write (*,*) "stepping into the loop...."
1	IT  	= IT + 1
        
        
!-----------------------------3 HOURLY OUTPUT-------------------------------------
!---------------- Write profile output every 3 hours
!-----------------Set TINC = 3 hours in Bl.dat
!	Write(*,*) 'Time,TTIME, Mtime,Tnext,Tnextflux,TINCFLUx', &
!	Time,Mtime,TTIME, Tnext,Tnextflux,Tincflux
        IF (MTIME .GE. TNEXT) THEN
        TNEXT = MTIME + TINC
	write(*,*) 'Numheatingavg',Numheatingavg
	Do j=1,Kz
	Diffheating(j)=Diffheating(j)/numheatingavg
	radheating(j)=radheating(j)/numheatingavg
	end do
	Numheatingavg=0
	
	Write(*,*)  'Calling Store Mtime', Mtime
 	CALL STORE(FIRST, TIME, TR, TAERO, T0, L, T, Q, U, V, Z, KM, KH, &
	Diffheating,Radheating)
        Do j=1,kz
	Diffheating(j)=0
	Radheating(j)=0
	End Do
	Write (*,*) 'after store'
!	write(*,*) "ustar, tstar"
!	write(*,*) ustar, tstar
	write(*,*) "t0, tr, taero, t1"
        write(*,*) t0, tr, taero, t(1)
!	write(*,*)
	ENDIF
	
! ---------------Write high frequency flux data nominally at 5 miniute
! ----------------set TINCFlux in Bl.dat
       If (Mtime.GE.TNEXTFLUX) then
!----------------Compute average over the Time increment of output
!		Write(*,*) 'Time, Mtime,Tnextflux,TINCFLUx', &
!	         Time,Mtime,Tnextflux,Tincflux
	XiSavg =Xisavg/Nstepsavg 
	Xilavg =Xilavg/Nstepsavg
	HOAVG =HOavg/Nstepsavg
	LHAVG =LHAVG/Nstepsavg
	GHFAVG =GHFAVG/Nstepsavg
	LWUAVG =LWUAVG/Nstepsavg
	Taeroavg = taeroavg/Nstepsavg
	Travg = Travg/Nstepsavg
!   Don't average To or columnenergy	
!	T0avg=T0avg/Nstepsavg
!	Columnenergyavg=Columnenergyavg/Nstepsavg
	Difftime=Nstepsavg*delt
	Radiativelossavg=Radiativelossavg/Nstepsavg
	
       
       Call STOREflux(FIRSTflux, TIME, TRAVG, TAEROAVG, T0AVG, 	&
		    Z0, Q0, XISAVG, XILAVG, HOAVG, LHAVG, GHFAVG, GHFRAVG,  &
		   
       LWUAVG,cb,Difftime,columnEnergyavg,Raduptop,Raddntop,radiativelossavg,p)
       
       TnextFlux= Mtime+Tincflux
       
	XiSavg =0
	Xilavg =0
	HOAVG =0
	LHAVG = 0
	GHFAVG =0
	GHFRAVG =0
	LWUAVG =0
	Taeroavg = 0
	Travg = 0
	T0avg=0
	Nstepsavg=0
	columnenergyavg=0
	radiativelossavg=0
	End If
	
       

!-------------------------ASSIMILATION SECTION----------------------------------

	CALL SAVE(TR_MODEL, TG_MODEL, TSTEPS, TR, T0, TSTEP)
	IF(MTIME .GE. ASSIM_TIME)THEN
        	CURRENT_ST = ASSIM_TIME / DELT
                ASSIM_TIME = MTIME + ASSIM_INC

                IF(ITRAD .EQ. 1)THEN
		WRITE(*,*)'NOT YET IMPLEMENTED'
                IF(TIME/(3600.) .GE. ASSIM_ST .AND. TIME/(3600.) .LE. ASSIM_EN)THEN

!               CALL MA_ADJ(OBS_TRADT, TG_MODEL, TR_MODEL, TIME, TTIME, ADJ_MA, &
!               		T_STEPS, LH, TSTEP, ADJ_FLAG, CURRENT_ST, 	&
!				ASSIM_TIME, DELT, ASSIM_INC)
                ADJ_FLAG = 1
                ENDIF
                ENDIF
        ENDIF

!--------------------------END ASSIMILATION SECTION-----------------------------

	MTIME= IT*DELT
	TSTEP= TSTEP + 1
	TIME = TIME0 + MTIME
	IF ( TIME .GE. RESETT) THEN
	TIME  = TIME  - RESETT
	TIME0 = TIME0 - RESETT
        IDY   = IDY+1
      	END IF
	IF (MTIME .LT. TTIME) THEN

!-------------------------TIME STEP MODEL INTEGRATION---------------------------

	CALL GET_FLUX (TR, TAERO, T0, T(1), Q(1), LAT, IMO, IDY, IRAD, DELT, TSTEPS,    &
                        XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U(1), V(1), Z(1), 	&
			Z0, ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 		&
			STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFR, QSTAR, &
			rlongdnsfc,Iradswitch)

        Write(*,*) 'After GetFlux'
	write(*,*) "AFTER FLUX: ustar, tstar, tr, taero, t0, t1"
	write(*,*) ustar, tstar, tr, taero, t0, t(1)

	CALL PROFILEK(U, V, T, Q, Z, DELZ, KZ, STABFUNC, ZI, LAT, T0, TSTAR, 		&
			USTAR, L, DELT, VONKAR, KM, KH)
        Write(*,*) 'After Call Profile'
	Write(*,*) 'T(1)', T(1)
	CALL CENTERDIFF(U, KM, DELZ, DELZM, DELT, KZ, USTAR, UFLUX, UG)	
		U(KZ)=UG						! UG is top boundary
	Write(*,*) 'After CenterDiff U'
	CALL CENTERDIFF(V, KM, DELZ, DELZM, DELT, KZ, USTAR, VFLUX, VG)	
		V(KZ)=VG						! VG is top boundary
	
	
	Do j=1,kZ
	Dfrate(j)=T(j)
	EndDo
	CALL DoubleCENTERDIFF(T, KH, DELZ, DELZM, DELT, KZ, USTAR, TSTAR, T(KZ))
!	CALL CENTERDIFF(T, KH, DELZ, DELZM, DELT, KZ, USTAR, 0, T(KZ))
	cp=1004
	Rhoa=1.1
	FluxIn=Rhoa*cp*ustar*tstar
	Write(*,*) ho,FluxIn
	Do j=1,kz
	dfrate(j)=-(dfrate(j)-T(j))/delt
	End Do
	Write(*,*) 'After DoubleCenterDiff'
	Write(*,*) 'T(1)', T(1)
	
	CALL CENTERDIFF(Q, KH, DELZ, DELZM, DELT, KZ, USTAR, QSTAR, Q(KZ))

	CALL CORIOLIS(U, V, UG, VG, Z, DELZ, ZI, LAT, DELT, KZ)
! ---------------------
!       Shortwave not turned on
!----------------------------
	Call Pressure (P,T,Zm)
	Write(*,*) 'After Pressure'
	Write(*,*) 'T(1)', T(1)
	Cosz=-3.
	If (Iradswitch.EQ.1) then

	Call RADIATION (Delt,P,T,Q,RHOP,ozone,cosz, &
	alb,em,Tr,fthrd,rlongdnsfc,raduptop,raddntop)
	Write(*,*) 'After Radiation'
	Write(*,*) 'T(1)', T(1)
	EndIF
	
	Call ColumnenergyCalc(FirstCE,T,p,Columnenergy,InitialColumnE,Radiativeloss,Fthrd,delt,zm)
	Write(*,*) 'After Columnenergy'
	Write(*,*) 'T(1)', T(1)
	
	Call SumHeatingProf(fthrd, dfRate,Numheatingavg,Radheating,Diffheating)
	
	!      Sum components for average later
			
	Call SumFlux(Nstepsavg,Tr,TrAVG,TAERO, TAEROAVG, &
	T0,T0AVG,XIS,XISAVG,XIL,XILAVG,HO,HOAVG,LH,LHavg,LWU,LWUAVG,   &
	GHF,GHFAVG,ColumnEnergy,Columnenergyavg,radiativeloss,radiativelossavg)

	
!------------------------END TIME STEP MODEL INTEGRATION------------------------
!	If (Tcount.GE.20) then
!	STOP
!	End If 
	IF(TCOUNT.EQ.TEMPOUT)THEN
	        WRITE(16,'(9F12.4)') TIME, T0, TR, TAERO, T(1),HO,L,Ustar,Tstar
		TCOUNT = 1
	ELSE
		TCOUNT = TCOUNT + 1	
	ENDIF
	
	GO TO 1

	ENDIF	

!-------------------------TIME STEP MODEL INTEGRATION---------------------------




	end program model




!-------------------------------------------------------------------------------
!
!				SUBROUTINES
!
!
!	INITIALIZE	! Initializes vertical profiles
!	ASSIM_INIT	! Initializes assimilation of insolation and/or Skin T Tend.
!	SURFACEINIT	! Initializes land surface model and three-T system
!	STORE		! Stores variables for model output
!
!	COMPUTEZ	! Computes vertical grid variables
!	EKMAN		! Provides U- and V- profiles if needed
!	SWNET		! Calculates SWnet term or uses observed insolation
!	STABLEFLUX	! Calculates surface similarity flux terms for stable BL	
!	GETSFCSPECHUM	! Calculates surface specific humidity	
!	GETQSTAR	! Calculates surface similarity humidity flux variable 
!	SAVE		! Saves the Skin and Ground T every assimilation interval
!	BISECT		! Bisect method Root Finder for diagnostic sfc. energy eqn.
!	SFCFLUX		! Calls Businger Form similarity relationships:
!	GETUSTAR	! Calculates Friction Velocity
!	GETTSTAR	! Calculates Friction Potential Temperature Scale
!	GETL		! Calculates Monin-Obukhov length
!	BLHEIGHT	! Calculates the BL height during unstable regimes
!
!	FLUX		! Solves TR, TAERO, TG, and surface similarity fluxes each dt
!	PROFILEK	! Calculates a diffusion (K-Profile) Profile for vertical diffusion
!	CENTERDIFF	! Centered in space Finite Differencing Solver
!	CORIOLIS	! Accounts for Coriolis force through PGF term

	SUBROUTINE INITIALIZE(Z, T, U, V, Q, P, rhop,ozone, TIME, MTIME, TTIME, TINC, &
				Tincflux,TNextflux,Z0,	 &
				      ZM, DELZ, DELZM, TEMPOUT, DELT, T0, 	 &
				      TSTEPS, FIRST, Firstflux, FLUXTYPE, TIME0, STABFUNC,	 &
				      UG, VG, VONKAR, LAT, IMO, IDY)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine initialized the model using profiles of
!	potential temperature, specific humidity, u- and v-component 
!	wind, height of each level, and the pressure of each level.
!	Profiles are read in from INPUT.DAT, e.g., 
!	Z(m)    TH(K)    U(m/s)  V(m/s)  Q(g/kg)  P(kPa)
!	2       300.00  -1.00   1.00    20.00   970.0
!	4       300.50  -1.00   1.00    20.50   969.5
!-------------------------------------------------------------------------------

	INCLUDE 'BL.incl'
	CHARACTER fname*124

	INCLUDE 'BL.DAT'		! Get geostrophic wind from this file
!------	READ IN Z, T, U, V, Q, AND P--------------------------------------------

!	WRITE(*,*)"Please Enter the Name of the Input File" // "(e.g., INPUT.DAT) : "
!	READ(*,*)fname
	fname='CASESAD991023003.txt'
	OPEN(unit=11, file=fname, status='old')
	READ(11,*)
	Write(*,*) 'Z(i), p(i),T(i), U(i), V(i)'
	DO I = 1, KZ
		READ(11,*)Z(I), p(i),T(I), U(I), V(I), rhop (i), Q(I), ozone(i)
		Write(*,*) Z(i), p(i),T(i), U(i), V(i) 
!                Q(I) = Q(I)*1.E-3	! Convert to kg/kg
!		P(I) = P(I)*100.	! Convert to Pa

!********* INSERTED FOR NSF RUNS *************
		U(I) = UG
		V(I) = VG 
!********* END INSERT FOR NSF RUNS ***********
	ENDDO
	WRITE(*,*)'Set U = Ug, V = Vg, for NSF Runs (Remove?)'
	WRITE(*,*)"Vertical Profiles Read, Setting up Vertical Grid...."
	
	

!------ SET UP TIME VARIABLES---------------------------------------------------

	WRITE(*,*)"Model start time is  :", TIME/3600., "UTC"
	WRITE(*,*)"Local (Central) time :", (TIME/3600.)-5.0
	MTIME = 0.0
	WRITE(*,*)"Model time is        :", MTIME
        TSTEPS = INT(TTIME/DELT)
	
!------ SET UP VERTICAL GRID----------------------------------------------------

	CALL COMPUTEZ(Z0, Z, ZM, DELZ, DELZM, KZ)
!-------------Open File for Model Input Parameters	
	OPEN(unit=18, file='ModelPrintInput.txt', status='replace')
	Write(18,'A12, A40') 'Snd File,   ', 'CASESAD991023003.txt'

	Write(18,'A6,F10.2,A6,F10.2') 'Ug=,',UG, ',Vg =,', VG
	
	Write (18,'A6,F10.3') 'Z0 =,',Z0
	
	Write(18,'A6,F10.3') 'Lat =,', Lat
	
	Write (18,'A14,I4,A20,I4') 'Fluxtype =,', Fluxtype, ', Stability Type =,', Stabfunc


! ---------------Open Surface high frequency Temperature and flux file	
	OPEN(unit=16, file= ' SFC_TEMPFlux.TXT', status='replace')
	
	end subroutine initialize

	SUBROUTINE ASSIM_INIT(IRAD, ITRAD, OBS_INSOL, OBS_TRADT, TSTEPS, ADJ_FLAG)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine reads in the observation data and places it
!	into the appropriate arrays.  If assimilation of Skin Temp or 
!	Insolation is turned OFF this subroutine does nothing.
!-------------------------------------------------------------------------------

	INTEGER				::	IRAD, ITRAD, I, TSTEPS, ADJ_FLAG
	REAL, DIMENSION(TSTEPS)		::	OBS_INSOL, OBS_TRADT

	IF(IRAD .EQ. 0) THEN
	   ADJ_FLAG = 0
           RETURN
        END IF

	IF(IRAD .EQ. 1)THEN

	OPEN(unit=1001, file='INSOL_IN')
	READ(1001, '(F15.8)')OBS_INSOL
	OPEN(unit=1002, file='TRAD_TEND_IN')
	READ(1002, '(F15.8)')OBS_TRADT
	
	ENDIF

	end subroutine assim_init


   SUBROUTINE STORE(FIRST, TIME, TR, TAERO, T0, L, T, Q, U, V, Z, KM, KH, &
                    Diffheating,Radheating)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!       This subroutine stores model variables at the output
!	interval and writes them out to files
!-------------------------------------------------------------------------------

	INCLUDE 'BL.incl'

	REAL			::	NULL
        
        Integer Itimeprint,jj
!---------OPEN THE OUTPUT FILE ONLY THE FIRST TIME THROUGH----------------------

	IF (FIRST.EQ.1) THEN
		FIRST = 0
	        OPEN(unit=12, FILE='MODELOUT.dat', status='replace')
! 	        WRITE(12,*)'TIME = ,', TIME/3600
!		WRITE(12,*)'   HGHT(m) level(j)  TEMP(TH)  MIXR(kg/kg)  &
!		U(m/s)   V(m/s) &
!		        KH       Turbulent   Radiational'
!	        WRITE(12,'(6F8.4,3F8.6)')(z(j),P(j), T(J), Q(J), U(J), V(J),&
!		KH(J),Diffheating(j),Radheating(j), J=1,KZ)
 	        

		WRITE(*,*)'Begining Model Integration'
	
	ENDIF

	
	WRITE(*,'(5F12.4)')TIME/3600., TR, TAERO, T(1), RIHALF

        WRITE(12,*)'TIME =,', TIME/3600.
	WRITE(12,*)'   HGHT(m) Press(mb) TEMP(TH)  MIXR(kg/kg)   U(m/s)   V(m/s) &
		   KM        KH       Turbulent Heating   Radiational Heating'
	
	WRITE(12,'(8(F12.6,A1),F12.6')(z(j),',',P(j)/1000.,',', T(J),',', Q(J),',', &
	U(J),',', V(J), ',',KH(J),',',Diffheating(j),',',Radheating(j), J=1,KZ)
 	        
      
 
        


          end subroutine store


   SUBROUTINE STOREflux(FIRSTflux, TIME, TRAVG, TAEROAVG, T0AVG, 	&
		    Z0, Q0, XISAVG, XILAVG, HOAVG, LHAVG, GHFAVG, GHFR, &
		   LWUAVG,cb,Difftime,columnEnergyavg,Raduptop,Raddntop, &
		   radiativelossavg,p)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!       This subroutine stores model variables at the output
!	interval and writes them out to files
!-------------------------------------------------------------------------------

	INCLUDE 'BL.incl'

	REAL			::	NULL
        REal Cb
        Integer Itimeprint,jj
!---------OPEN THE OUTPUT FILE ONLY THE FIRST TIME THROUGH----------------------

	IF (FIRSTFlux.EQ.1) THEN
		FIRSTflux = 0
 	 	Write (*,*) 'Opening Flux File'
 	        OPEN(unit=13, FILE='FLUXES.dat', status='replace')
		
		WRITE(13,'A120')'Time(hr)  XILAVG*DT   LWUAVG*DT   HOAVG*DT  &
		  GHFAVG*Dt  RADUPTOP    RADDNTOP  cb*t0  Column Energy  &
		  Radiativeloss'

		
	ENDIF


        WRITE(13,'(7(F12.3,A1), 2(f14.1,A1),2(F12.3,A1))') Time/3600.,',',XILAVG*Difftime,',', &
	LWUAVG*Difftime,',', HOAVG*difftime,',', GHFAVG*difftime,',', &
        Raduptop*difftime,',',raddntop*difftime,',',cb*t0AVG,',',columnenergyavg, &
	',', radiativelossavg*Difftime
        

      end subroutine storeflux

!
!	1	DDD	BBBB	L	MM   MM
!      11	D  D	B  B	L	M  M  M
!	1	D  D	BBB	L	M  M  M
!	1	D  D	B  B	L	M     M
!     11111	DDD	BBBB	LLLLL	M     M
!

	SUBROUTINE COMPUTEZ(Z0, Z, ZM, DELZ, DELZM, KZ)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine uses the vertical height grid provided 
!	in the input.dat file to create the mid-level array and 
!	calculate and create the delta-z and delta-zm arrays.
!-------------------------------------------------------------------------------
	INTEGER			::	KZ, I
	REAL, DIMENSION(KZ)	::	Z, ZM, DELZ, DELZM 	 
	REAL			::	Z0

	DO I = 1, KZ-1
		ZM(I)    = 0.5*(Z(I) + Z(I+1))
	ENDDO
		ZM(KZ)   = Z(KZ)
		DELZ(1)  = Z(1) - Z0
		DELZM(1) = Z0 + 0.5*( Z(1) - Z0 )
	DO I = 2, KZ
		DELZ(I)  = Z(I)  - Z(I-1)
		DELZM(I) = ZM(I) - ZM(I-1)
	ENDDO

	WRITE(*,*)"Vertical Grid Complete, Continuing Initialization..."

	end subroutine computez		


   SUBROUTINE SAVE(TR_MODEL, TG_MODEL, TSTEPS, TR, TG, IT)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!       This subroutine saves model variables every time step
!-------------------------------------------------------------------------------
	INTEGER                 ::      ITI, TSTEPS
	REAL                    ::      TR, TG, IT
	REAL, DIMENSION(TSTEPS) ::      TR_MODEL, TG_MODEL
	
	ITI = INT(IT)

        TR_MODEL(ITI) = TR
        TG_MODEL(ITI) = TG

end subroutine save



   SUBROUTINE EKMAN(UGARR, VGARR, Z, KZ, U, V, KM, KH, PHI)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	SETS UP A U AND V PROFILE IF ONE IS NOT AVAILABLE		                                                                           !
!-------------------------------------------------------------------------------

	INTEGER                 ::      KZ, I
	REAL, DIMENSION(KZ)     ::      UGARR, VGARR, U, V, KM, KH, Z
	REAL                    ::      C, UG, VG, PI, OMEGA, F, GAM, PHI, ZI

!       Here C = Km = Kh

	C       = 1.0
	UG      = UGARR(1)
	VG      = VGARR(1)
	PI      = ATAN2(0.0,-1.0)
	OMEGA   = 2.0*PI/(24.0*3600.0)
	F       = 2.0*OMEGA*SIND(PHI)
	GAM     = SQRT(ABS(F)/(2.0*C))
		
	ZI      = PI * SQRT(2.0/ABS(F))
	DO I = 1, KZ

        	IF(Z(I) .LT. ZI)THEN

        	U(I) = UG-EXP(-GAM*Z(I))*( UG*COS(GAM*Z(I))+VG*SIN(GAM*Z(I)) )
        	V(I) = VG-EXP(-GAM*Z(I))*( VG*COS(GAM*Z(I))-UG*SIN(GAM*Z(I)) )

        	ELSE

        	U(I) = UG
        	V(I) = VG

        	ENDIF

        	KM(I) = C
        	KH(I) = C

        	WRITE(*,*)U(I), V(I)
	END DO

end subroutine ekman

	SUBROUTINE SWNET(LAT, IMO, IDY, IRAD, TSTEP, TIME, DELT, TSTEPS, XIS, OBS_INSOL)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine calculates the Net shortwave radiation
!   		PHI 	= LATITUDE(IN RADIANS)
!   		N	= DAY OF THE YEAR
!   		HRA	= HOUR ANGLE
!   		COSZ	= COS OF ZENITH ANGLE
!   		AR	= A**2/R**2
!   		RD	= GAS CONST FOR DRY AIR
!-------------------------------------------------------------------------------

	INCLUDE 'SFC.incl'

	INTEGER, DIMENSION(12)          ::      MO
	INTEGER                         ::      IMO, IDY, I, IRAD, TSTEPI, TSTEPS
	REAL                            ::      PI, DTR, PHI, DELT, TIME0
	REAL                            ::      AQ31, AQ32, AQ33, SCAT1
	REAL                            ::      TSTEP
	REAL, DIMENSION(TSTEPS)         ::      OBS_INSOL
	DATA MO /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /

	INCLUDE 'SFC.DAT'
	

	PI      = 3.1415926
	DTR     = PI/180.
	PHI     = LAT*DTR
	XKR     = 1.18*OMEGA
	TSTEPI  = INT(TSTEP)

	IF (IRAD .EQ. 1)THEN

!               XIS = OBS_INSOL(TSTEPI)*(1.-ALB)	
                XIS = OBS_INSOL(TSTEPI)
	
	ELSE

!--------------------------- Atmospheric absorption-----------------------------

	SCAT1 = 0.000949 * (PRESS/100.)
	SCAT  = SCAT1 + 0.051
	AQ31  = 2.9 * PW
	AQ32  = ( 1.0 + (141.5*PW) ) ** (0.645)
	AQ33  = 5.925 * PW
	AQ3   = AQ31 / ( AQ32 + AQ33 )

!----------------------SWD (Shortwave Downward Radiation)-----------------------				

	N=(MO(IMO)+IDY)*1.
	DEC=23.45*SIND(360.*((284.+N)/365.))
	HR= ACOSD(-TAND(DEC)*TAND(LAT))
	THR=TIME/3600.
	THR=THR-6.0               ! Converts to LST
	IF(THR.LT.12.)THEN
	HRA=-15.*(12.-THR)
	ELSE IF(THR.GE.12.)THEN
	HRA=15.*(12.-THR)
	END IF
	COSZ=COSD(DEC)*COSD(LAT)*COSD(HRA)+SIND(DEC)*SIND(LAT)
	DO=2.*PI*N/365.
	AR=1.000110 +.034221 *COS(DO)+.001280*SIN(DO)+.000719*COS(2.*DO)+.000077*SIN(2.*DO)

	IF(ACOSD(COSZ).LT.(90.))THEN
                        XIS = (1.- ALB) * SO * AR * COSZ * (SCAT-AQ3)
        ELSE
                        XIS = 0.0
        ENDIF

	ENDIF

!***********NSF TESTING************REMOVE EVENTUALY
	XIS = 0.0
!********************************

	end subroutine swnet

   SUBROUTINE STABLEFLUX(T1, TAERO, U, V, Z, Z0, USTAR, TSTAR, L, VONKAR, STABFUNC, RIHALF)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine uses the Ri Formulations of surface fluxes
!	proposed by England and McNider (1995; BL-Met)
!-------------------------------------------------------------------------------
	INTEGER	::	STABFUNC
	Double Precision T1
	REAL    ::       TAERO, U, V, Z, Z0, USTAR, TSTAR, L, VONKAR
	REAL    ::      TGRAD, ZGRAD, TREF, GOVTH, WND, RI, RIC, FUNC, U2, V2, SHEAR
	REAL	::	RIHALF, H1, USTAR2, USTARTHSTAR, A, B
	REAL	::	LOGZOZ0, FUNC2

	TGRAD  = T1-TAERO
	ZGRAD  = Z-Z0
	TREF   = (T1+TAERO)/2.0
	H1     = ( 0.5*(Z0+Z) ) * LOG(Z/Z0)
	GOVTH  = 9.81/TREF
	U2     = U*U
	V2     = V*V
	WND    = ( U2 + V2 ) ** (0.5)
	RIHALF = GOVTH * H1 * (TGRAD/(WND*WND))
	RIC    = 0.25
	LOGZOZ0= LOG(Z/Z0)

	IF(RIHALF .GE. 0.0)THEN
	  IF(RIHALF.GE.RIC)THEN
		FUNC2   = 0.0
		USTAR   = 0.0
		TSTAR   = 0.0
		L       = 1.0
	  ELSE
		IF (STABFUNC .EQ. 1)THEN
		  FUNC2   = ( ( (RIC-RIHALF) / RIC ) ** (2.0) )         
		ENDIF
		IF (STABFUNC .EQ. 2)THEN
		  A = 0.61
		  B = 2.61
		  FUNC2   = 1.0 /  ( ((RIHALF**A)  + 1.0) ** (B) )
		ENDIF		
		IF (STABFUNC .EQ. 3)THEN
		  A     = 0.70
		  B     = 6.27
		  FUNC2   = 1.0 /  ( ((RIHALF**A)  + 1.0) ** (B) )
		ENDIF

		FUNC = SQRT ( FUNC2)

!		USTAR2 = ( (VONKAR*VONKAR) * FUNC2 * (WND*WND) ) / (LOGZOZ0*LOGZOZ0)
!		USTAR  = SQRT(USTAR2)
		USTAR  = (VONKAR * FUNC * WND ) / (LOGZOZ0)

!		USTARTHSTAR = ( (VONKAR*VONKAR) * FUNC2 * WND * TGRAD ) / (LOGZOZ0*LOGZOZ0)
!		TSTAR  = USTARTHSTAR/USTAR
		TSTAR  = (VONKAR * FUNC * TGRAD ) / (LOGZOZ0)
		IF ( TSTAR .GT. 0.) THEN
		  L = 0.5 * (T1+Taero) * (USTAR*USTAR) / (VONKAR*TSTAR*GOVTH)
		ELSE
		  L = 1.E+05
		END IF
	  ENDIF

	ELSE	! special case for unstable condition, should not be here

		USTAR = ( (VONKAR * WND) ) / (LOGZOZ0)
         	TSTAR = VONKAR * (TGRAD)   / (LOGZOZ0)
		L     = 0.5 * (T1+Taero) * (USTAR*USTAR) / (VONKAR*TSTAR*GOVTH)

	END IF


	
end subroutine stableflux

   SUBROUTINE GETSFCSPECHUM(QA, TA, PRESS)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!       This subroutine calculates the specific humidity at
!	the roughness height as a function of the aerodynamic 
!	temperature and surface pressure
!	McCumber formulation (Pielke, 2nd ed. p. 420)
!-------------------------------------------------------------------------------
	REAL            ::      QA, TA, PRESS
	REAL            ::      QVS, ES

	ES  = 6.1078* EXP(17.269*( (TA-273.16)/(TA-35.86) ) )
	QVS = 0.622* (ES / ( (PRESS/100.)-(0.378*ES) ) )
	QA  = QVS

end subroutine getsfcspechum

   SUBROUTINE GETQSTAR(L, QSTAR, Z, Q, QA, Z0, VONKAR, stability)

!-------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------
!       This subroutine calculates the friction specific humidity, QSTAR
!-----------------------------------------------------------------------------------------------

	REAL            ::      L, QSTAR, T, Z, T0, Q, QA, Z0, VONKAR
	REAL            ::      PHIH, SIH, DENOM, ZOL, stability
	REAL		::	LOGZOZ0

	LOGZOZ0 = LOG(Z/Z0)
	IF (L .NE. 0.) THEN
	  ZOL = Z/L
	ELSE
	  ZOL = 1.E3 * stability
	END IF

	IF(stability.LT.0.0 .AND. ZOL .LT. 0.)THEN
!          ZOL = Z/L
           PHIH = (1.0 -11.6*ZOL)**0.25
           SIH  = 2.0*LOG(0.5*(1.+(PHIH*PHIH)))
           IF (SIH .LT. LOGZOZ0) THEN
	     DENOM = 0.95*(LOGZOZ0 - SIH)
             QSTAR = VONKAR * (Q-QA) / DENOM
	   ELSE
             QSTAR = 1.E+03 * VONKAR * (Q-QA) / LOGZOZ0
	   END IF

        ELSE IF (stability.GT.0.0 .AND. ZOL .GT. 0.)THEN
!           ZOL = Z / L
!           IF(ZOL.GT.2.0)THEN
!             ZOL = 2.0
!           ENDIF

           SIH   = 8.0 * ZOL
           DENOM = 0.95*(LOGZOZ0 + SIH)
           QSTAR = VONKAR * (Q-QA) / DENOM

      	ELSE

           QSTAR = VONKAR * (Q-QA) / LOGZOZ0

      	ENDIF

end subroutine getqstar


   SUBROUTINE SFCFLUX(L, USTAR, TSTAR, T1, Z, Z0, VONKAR, TAERO, stability, U, V, T0, G)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	Uses Businger flux formulations
!-------------------------------------------------------------------------------
	INTEGER	::	LIMIT, IT
	Double Precision T1
	REAL    ::      L, USTAR, TSTAR, Z, Z0, TOLL, OLDL
	REAL    ::      VONKAR, TAERO, stability, U, V, T0, G

	LIMIT	= 500
	IT	= 0
	TOLL    = (1.e-3)

	stability = T1 - TAERO
	IF(stability.LT.0.0.AND.L.GT.0.0)  L = -1.0*L	
	IF(stability.GT.0.0.AND.L.LT.0.0)  L = -1.0*L	

7       OLDL = L
		CALL GETTSTAR(L, TSTAR, T1, Z, T0, Z0, VONKAR, TAERO, stability)
                CALL GETUSTAR(L, USTAR, U, V, Z, Z0, VONKAR, stability)
                CALL GETL(L, Z, USTAR, TSTAR, T1, VONKAR, G, TAERO)
	IT = IT + 1
        IF ( IT .GT. LIMIT )STOP'ERROR IN SFCFLUX - LIMIT EXCEEDED'
!	IF ( L .EQ. 100000.0)STOP'ERROR IN SFC FLUX - L TOO HIGH'
        IF (ABS((L-OLDL)/L) .GT. TOLL) GO TO 7
	
end subroutine sfcflux

   SUBROUTINE GETUSTAR(L, USTAR, U, V, Z, Z0, VONKAR, stability)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine calculates the friction velocity, USTAR
!-------------------------------------------------------------------------------

	REAL	::	L, USTAR, U, V, Z, Z0, VONKAR, stability
	REAL	::	SIM, DENOM, PHIM, PI, SPEED, T, T0, ZOL
	REAL	::	LOGZOZ0

	PI      = ATAN2(0.0,-1.0)
	SPEED   = SQRT(U*U+V*V)
	LOGZOZ0 = LOG(Z/Z0)
	IF (L .NE. 0.) THEN
	  ZOL = Z/L
	ELSE
	  ZOL = 1.E03*stability
	END IF

	IF(stability.LT.0.0 .AND. ZOL .LT. 0.)THEN
!        	ZOL = Z/L
         	PHIM = (1.0 - 19.0*ZOL)**0.25
         	SIM  = 2.0*LOG(0.5*(1.+PHIM)) + LOG(0.5*(1.+PHIM*PHIM)) -              &
               		2.0*ATAN(PHIM) + 0.5*PI
         	IF (SIM .LT. LOGZOZ0 ) THEN
		  DENOM = LOGZOZ0 - SIM
        	  USTAR = VONKAR * SPEED / DENOM
		ELSE
        	  USTAR = 1.E+03 * VONKAR * SPEED / LOGZOZ0
		END IF
        ELSE IF (stability.GT.0.0 .AND. ZOL .GT. 0.)THEN
!        	ZOL = Z/L
         	SIM   = 5.3 * ZOL
		DENOM = LOGZOZ0 + SIM
        	USTAR = VONKAR * SPEED / DENOM
      	ELSE
        	USTAR = VONKAR * SPEED / LOGZOZ0
      	END IF

	IF (USTAR .LE. 0.) THEN
	  write(*,*)
	  write(*,*) " Z, L, ZOL, LOGZOZ0: ", Z, L, ZOL, LOGZOZ0
	  write(*,*) " stability, PHIM, SIM, ustar: ", stability, PHIM, SIM, ustar
	  write(*,*) " u, v, speed: ", u, v, speed
	  STOP'GETUSTAR: SOMETHING IS WRONG HERE, CHECK THE CODE'
	END IF

end subroutine getustar

   SUBROUTINE GETTSTAR(L, TSTAR, T, Z, T0, Z0, VONKAR, TAERO, stability)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!       This subroutine calculates the friction temperature, TSTAR
!-------------------------------------------------------------------------------

	REAL	::	L, TSTAR, Z, T0, Z0, VONKAR, TAERO
	REAL	::	PHIH, SIH, DENOM, ZOL, stability
	REAL	::	LOGZOZ0
        Double Precision T
	STABILITY = T - TAERO
	LOGZOZ0 = LOG(Z/Z0)
	IF (L .NE. 0.) THEN
	  ZOL = Z/L
	ELSE
	  ZOL = 1.E03*STABILITY
	END IF

	IF(stability.LT.0.0 .AND. ZOL .LT. 0.)THEN
!		ZOL = Z/L
		PHIH = (1.0 -11.6*ZOL)**0.25
         	SIH  = 2.0*LOG(0.5*(1.+(PHIH*PHIH)))
         	IF (SIH .LT. LOGZOZ0) THEN
		  DENOM = 0.95*(LOGZOZ0 - SIH)
         	  TSTAR = VONKAR * (T-TAERO) / DENOM
		ELSE
        	  TSTAR = 1.E+03 * VONKAR * (T-TAERO) / LOGZOZ0
		END IF
        ELSE IF (stability.GT.0.0 .AND. ZOL .GT. 0.)THEN
!        	ZOL = Z/L
		SIH   = 8.0 * ZOL
		DENOM = 0.95*(LOGZOZ0 + SIH )
         	TSTAR = VONKAR * (T-TAERO) / DENOM
  	ELSE
         	TSTAR = VONKAR * (T-TAERO) / LOGZOZ0
      	END IF

end subroutine gettstar

   SUBROUTINE GETL(L, Z, USTAR, TSTAR, T1, VONKAR, G, TAERO)

IMPLICIT NONE
!--------------------------------------------------------------------------------
!       This subroutine calculates the Monin Obukov Length scale L
!--------------------------------------------------------------------------------
        Double Precision T1
	REAL	::	USTAR, TSTAR, L, Z, TAERO, VONKAR, G
	REAL	::	DENOM, USTAR2, TBAR

	TBAR  = (T1 + TAERO) / 2.0

	IF (TSTAR.NE.0.0)THEN
        	USTAR2 =  USTAR*USTAR
        	DENOM  =  VONKAR * G * TSTAR
        	L      = TBAR * USTAR2 / DENOM
      	ELSE
        	L      = 1.E5
      	END IF

end subroutine getl

   SUBROUTINE BLHEIGHT(LAT, T, USTAR, TSTAR, DELT, ZI)

IMPLICIT NONE
!--------------------------------------------------------------------------------
!       This subroutine calculates the boundary layer height as given
!	in Pielke 2nd edition p. 192
!--------------------------------------------------------------------------------
        Double Precision T
	REAL            :: LAT, T0, USTAR, TSTAR, DELT, ZI
	REAL            :: G, GOVTH, WSTAR, DZI, F
	REAL            :: DUM1, DUM2, DUM3, DUM4, DUM5, DUM6
	REAL            :: UST3, UST2, WST3, WST2, ZI2OVTH, OMEGA

	OMEGA	= 7.292e-5
	G	= 9.81
	GOVTH	= G / T

	F	= 2. * OMEGA * SIND(LAT)		! Coriolis Term

	IF(TSTAR .GT. 0.0)THEN				! For the rare case were
		ZI = (0.33*USTAR)/ABS(F)		! near neutral conditions
	RETURN						! sneak in.
	ENDIF

!	IF(ZI .EQ. 0.0)THEN
!		ZI = (0.33*USTAR)/ABS(F)		! Initial BL height
!	ENDIF
		WSTAR   = 	(-1.0 * GOVTH * USTAR * TSTAR * ZI)**(.33)
		UST2    = 	USTAR*USTAR
		UST3    = 	UST2*USTAR
		WST2    = 	WSTAR*WSTAR
		WST3    = 	WST2*WSTAR
		ZI2OVTH = 	ZI*ZI / T
		DUM1    =       1.1*UST3
                DUM2    =       3.3*UST2*F*ZI
                DUM3    =       5. / 1000.              ! 5/1000 ~ lapse rate above ZI
                DUM4    =       G*ZI2OVTH*DUM3
                DUM5    =       9.0*WST2
                DUM6    =       7.2*UST2
                DZI     =       (1.8*(WST3+DUM1-DUM2)) / (DUM4+DUM5+DUM6)
                ZI      =       ZI + (DZI*DELT)

end subroutine blheight

!       1       DDD     BBBB    L       MM   MM
!      11       D  D    B  B    L       M  M  M
!       1       D  D    BBB     L       M  M  M
!       1       D  D    B  B    L       M     M
!     11111     DDD     BBBB    LLLLL   M     M
!


   SUBROUTINE PROFILEK(U, V, T, Q, Z, DELZ, KZ, STABFUNC, ZI, LAT, T0, 		&
			TSTAR, USTAR, L, DELT, VONKAR, KM, KH)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!   This subroutine produces the necessary K-profile for vertical diffusion
!-------------------------------------------------------------------------------

	INTEGER			::	KZ, STABFUNC, I, IFLAG
	Double Precision T(KZ)
	REAL, DIMENSION(KZ)	::	U, V, Q, Z, DELZ, RL, KM, KH
	REAL			::	HS, A, B, FUNC, RI
	REAL			::	VONKAR, ZDIFF, KZI, LAT, T0, USTAR, TSTAR
	REAL			::	DELT, ZI, L, TOL, PHIM, PHIH, KMHS, KHHS
	REAL			::	TLAYER, RIC, UDIFF, VDIFF, DENOM, DT
	REAL			::	GOVTH, ZM, SHEAR, G, WEIGHT, PHIMI, PHIHI
	REAL			::	KMDIFF, KHDIFF, HS2, HS3, HS4, HS5
	REAL			::	TERM1, TERM2, KMDER, KHDER
	
!------ Miscellaneous Settings--------------------------------------------------
	HS	= Z(1)		! Height of Surface Layer, taken as 1st model level
	KZI	= 0.000001	! Background Diffusion coefficient
	G	= 9.81		! Gravity
	RIC	= 0.25		! Critical Richardson number
	TOL	= 1.E-6		! Tolerance value for minimum Denominator value
	ZDIFF	= ZI - HS	! Depth from surface layer height to top of BL

!------ Unstable K-profile------------------------------------------------------
	IF (L .LT. 0.0) THEN
	
		CALL BLHEIGHT(LAT, T(1), USTAR, TSTAR, DELT, ZI)

		ZDIFF = ZI  - HS
		
		PHIMI = (1.0 - 19.0*HS/L)**0.25		! Assuming VonKarman = 0.4
		PHIHI = SQRT(1.0 - 11.6*HS/L)		! Inverse form

!------ Begin Forming O'Brian Profile components -------------------------------

		KMHS  = VONKAR*HS*USTAR*PHIMI		! Surface Layer diffusion
        	KHHS  = VONKAR*HS*USTAR*PHIHI		
	        KMDIFF= KMHS - KZI
         	KHDIFF= KHHS - KZI

         	HS2   = HS*HS
         	HS3   = HS2*HS
         	HS4   = HS3*HS
         	HS5   = HS4*HS

         	TERM1 = 4.0*HS3 - 75.0*HS4/L
         	TERM2 = 1.0/(HS4 - 15.0*HS5/L)**0.75
         	KMDER = 0.25 * VONKAR * USTAR * TERM1 * TERM2
         	TERM1 = 2.0*HS - 27.0*HS2/L
         	TERM2 = 1.0 / SQRT(HS2 - 9.0*HS3/L)
         	KHDER = 0.5 * VONKAR * USTAR * TERM1 * TERM2		

!------ Calculate diffusion coefficient for each layer -------------------------

		DO I = 1, KZ-1

			ZM = Z(I)+0.5*DELZ(I+1)		! Mid-level height
			IF(ZM .LT. ZI)THEN		! Only fill to BL height
			
!------ Fill Polynomial (O'Brian 1970) -----------------------------------------
			KM(I+1) = KZI + ( (ZI-ZM)*(ZI-ZM)/(ZDIFF*ZDIFF) ) * ( KMDIFF + (ZM - HS) *  &
		                       ( KMDER + 2.0*KMDIFF/ZDIFF) )
			KH(I+1) = KZI + ( (ZI-ZM)*(ZI-ZM)/(ZDIFF*ZDIFF) ) * ( KHDIFF + (ZM - HS) *  &
		                       ( KHDER + 2.0*KHDIFF/ZDIFF) )
			
			IF(KM(I+1).LT.KZI)THEN
                	KM(I+1) = KZI				! Don't let Km, Kh get too small
                	ENDIF
                	IF(KH(I+1).LT.KZI)THEN
                	KH(I+1) = KZI
                	ENDIF
			
!------ Set Km, Kh above BL height equal to background value--------------------
			ELSE
				KM(I+1) = KZI
				KH(I+1) = KZI
			ENDIF	
		ENDDO

	ENDIF

!------ Stable K-profile--------------------------------------------------------
	IF (L .GE. 0.0) THEN

	IFLAG 	= 1		! Logical flag for determination of nocturnal BL height
	
	RL(1) 	= VONKAR*( Z(1) / 2.0 )
	PHIM	= 1.0 + 5.3*HS/L
	PHIH	= 0.74 + 8.0*HS/L
	KMHS	= VONKAR*HS*USTAR / PHIM
	KHHS	= VONKAR*HS*USTAR / PHIH

		DO I = 1, KZ-1

!------ Calculate Richardson Number for each layer-------------------------------
			TLAYER		= ( (T(I+1)+T(I)) / 2.0 )	! Temperature of layer
			ZM      	= Z(I)+0.5*DELZ(I+1)		! Mid-level height
			UDIFF       	= ( U(I+1)-U(I) ) / DELZ(I+1)	
			VDIFF       	= ( V(I+1)-V(I) ) / DELZ(I+1)
            		DENOM       	= (UDIFF*UDIFF + VDIFF*VDIFF)
            		DT          	= (T(I+1) - T(I)) / DELZ(I+1)
            		GOVTH       	= G / TLAYER

			IF(DENOM.LT.TOL)THEN
                		DENOM 	= TOL				! For near shear-less situations
            		ENDIF

				RI    	= GOVTH  * DT / DENOM		! Richardson Number Calculation

			IF(RI .LT. 0.0)THEN
				RI	= 0.0
			ENDIF

!------ Calculate Mixing Length--------------------------------------------------
			IF (ZM .LT. 200. ) THEN
				RL(I+1) = VONKAR*ZM
			ELSE
				RL(I+1) = 70.
			ENDIF

!------ Calculate Layer Wind Shear-----------------------------------------------
			SHEAR		= SQRT(UDIFF*UDIFF + VDIFF*VDIFF)
!------ Set up Stability Function------------------------------------------------
	                IF (STABFUNC .EQ. 1)THEN
        	        	FUNC   	= ( (RIC-RI) / RIC ) ** (2.0)
        	        ENDIF
        	        IF (STABFUNC .EQ. 2)THEN
        	        	A 	= 0.61
        	        	B 	= 2.61
        	        	FUNC   	= 1.0 /  ( ((RI**A)  + 1.0) ** (B) )
        	        ENDIF
        	        IF (STABFUNC .EQ. 3)THEN
        	        	A     	= 0.70
        	        	B     	= 6.27
        	        	FUNC   	= 1.0 /  ( ((RI**A)  + 1.0) ** (B) )
        	        ENDIF

!------ Calculate Diffusion Profile (K-Profile)--------------------------------- 						
			IF (RI .LE. RIC) THEN
 	            		IF (DT .GT. 0.0) THEN
               				KM(I+1) = 1.1 * FUNC * RL(I+1)*RL(I+1) * SHEAR
               				KH(I+1) = 1.1 * FUNC * RL(I+1)*RL(I+1) * SHEAR
            			ELSE IF (DT .LE. 0.0)THEN
               				KM(I+1) = ( (1.0 - 18.0*RI)**(0.5) ) * RL(I+1)*RL(I+1) * SHEAR
               				KH(I+1) = ( (1.0 - 18.0*RI)**(0.5) ) * RL(I+1)*RL(I+1) * SHEAR
            			ENDIF
                		IF(KM(I+1).LT.KZI)THEN
                			KM(I+1) = KZI
                		ENDIF
                		IF(KH(I+1).LT.KZI)THEN
                			KH(I+1) = KZI
                		ENDIF
			ELSE
				IF ( IFLAG .EQ. 1 ) THEN
					ZI    = Z(I)
					IFLAG = 0
				END IF
				 	KM(I+1) = KZI
               				KH(I+1) = KZI	
			ENDIF
		ENDDO
	ENDIF
!------ Smooth Profile using weighting------------------------------------------
! Turn off smoothing by comment 
!	KM(1) = KZI
!        KH(1) = KZI

!	DO I = 2, KZ-1
!       	WEIGHT = RL(I) / DELZ(I)
!	                IF (WEIGHT .GT. 1.0) WEIGHT = 1.0
!                	KM(I) = WEIGHT*KM(I-1) + KM(I) + WEIGHT*KM(I+1)
!                        KM(I) = KM(I) / (1.0 + 2.0*WEIGHT)
!                        KH(I) = WEIGHT*KH(I-1) + KH(I) + WEIGHT*KH(I+1)
!                        KH(I) = KH(I) / (1.0 + 2.0*WEIGHT)
!        ENDDO

end subroutine profilek

   SUBROUTINE CENTERDIFF(F, K, DELZ, DELZM, DELT, KZ, USTAR, PHISTAR, FTOP)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!   This subroutine performs vertical diffusion using a centered-difference 
!   solver.
!-------------------------------------------------------------------------------

	INTEGER			::	I, KZ
	REAL, DIMENSION(KZ)	::	F, K, DELZ, DELZM
	REAL			::	HOVZ, FLOW, FUP, USTAR, PHISTAR, FTOP, DELT

	DO I = 1, KZ

		IF(I.GT.1.AND.I.LT.KZ)THEN

                	HOVZ = DELT / DELZM(I)
                	FUP  = K(I+1) * ( ( F(I+1)-F(I) ) / DELZ(I+1) )
                	FLOW = K(I) * ( ( F(I) - F(I-1) ) / DELZ(I) )
                	F(I) = F(I) + HOVZ*( FUP - FLOW )
        
		ENDIF

        	IF(I.EQ.1)THEN

                	HOVZ = DELT / DELZM(I)
	                FUP  = K(2) * ( ( F(2)-F(1) ) / DELZ(2) )
                	F(1)  = F(1) + HOVZ *  ( FUP  - USTAR*PHISTAR )

	        ENDIF

        	IF(I.EQ.KZ)THEN

                	HOVZ = DELT / DELZM(I)
                	FUP  = K(KZ) * ( ( FTOP - F(KZ) ) / DELZM(KZ) )
	
                     		
                	FLOW = K(KZ) * ( ( F(KZ) - F(KZ-1) ) / DELZM(KZ) )
                	F(KZ) = F(KZ) + HOVZ * ( FUP - FLOW )
        	ENDIF
	
	ENDDO

end subroutine centerdiff
  SUBROUTINE doubleCENTERDIFF(F, K, DELZ, DELZM, DELT, KZ, USTAR, PHISTAR, FTOP)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!   This subroutine performs vertical diffusion using a centered-difference 
!   solver.
!-------------------------------------------------------------------------------

	INTEGER			::	I, KZ
	Double Precision F(Kz)
	REAL, DIMENSION(KZ)	::	 K, DELZ, DELZM
	REAL			::	HOVZ, FLOW, FUP, USTAR, PHISTAR, FTOP, DELT

	DO I = 1, KZ

		IF(I.GT.1.AND.I.LT.KZ)THEN

                	HOVZ = DELT / DELZM(I)
                	FUP  = K(I+1) * ( ( F(I+1)-F(I) ) / DELZ(I+1) )
                	FLOW = K(I) * ( ( F(I) - F(I-1) ) / DELZ(I) )
                	F(I) = F(I) + HOVZ*( FUP - FLOW )
        
		ENDIF

        	IF(I.EQ.1)THEN

                	HOVZ = DELT / DELZM(I)
	                FUP  = K(2) * ( ( F(2)-F(1) ) / DELZ(2) )
                	F(1)  = F(1) + HOVZ *  ( FUP  - USTAR*PHISTAR )

	        ENDIF

        	IF(I.EQ.KZ)THEN

                	HOVZ = DELT / DELZM(I)
                	FUP  = K(KZ) * ( ( FTOP - F(KZ) ) / DELZM(KZ) )
	
                     		
                	FLOW = K(KZ) * ( ( F(KZ) - F(KZ-1) ) / DELZM(KZ) )
                	F(KZ) = F(KZ) + HOVZ * ( FUP - FLOW )
        	ENDIF
	
	ENDDO

end subroutine doublecenterdiff

   SUBROUTINE CORIOLIS(U, V, UG, VG, Z, DELZ, ZI, LAT, DELT, KZ)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!   This subroutine calculates the effect of the coriolis parameter (PGF term)
!   and is the second part of the operator splitting
!-------------------------------------------------------------------------------
	INTEGER			::	KZ, I
	REAL, DIMENSION(KZ)	::	U, V, Z, DELZ
	REAL			::	UG, VG, DELT, LAT, OMEGA, F, ZI

	OMEGA   = 7.292E-5		! Angular Velocity of Earth, 2*PI/(24*3600) (1/sec)

	F   = 2. * OMEGA * SIND(LAT)

	IF (ZI .LT. Z(2)) ZI=Z(2)

	DO I=1,KZ

            U(I) = U(I) + DELT * (  F * (V(I) - VG) )
            V(I) = V(I) + DELT * ( -F * (U(I) - UG) )

      	END DO

end subroutine coriolis


!=========================================================================
!
! The following is the driver for calculating surface fluxes
!
!=========================================================================
!

!       1       DDD     BBBB    L       MM   MM
!      11       D  D    B  B    L       M  M  M
!       1       D  D    BBB     L       M  M  M
!       1       D  D    B  B    L       M     M
!     11111     DDD     BBBB    LLLLL   M     M
!

   SUBROUTINE FLUX   (TR, TAERO, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS, &
			XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, Z0,	&
			ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 	&
			STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFR,        &
			QSTAR,rlongdnsfc,Iradswitch)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine solves for the following variables:
!	TAERO		= Aerodynamic Temperature using Zilitinkevich(1970)
!	USTAR, TSTAR	= Friction velocity and temperature (U*, T*) used
!			  in flux calculations
!	HO, LH		= Sensible, Latent Heat Flux
!	UFLUX, VFLUX	= w'u', w'v' (Pielke v1, p.140)
!-------------------------------------------------------------------------------
	Integer                 ::      Iradswitch
	REal                    ::      rlongdnsfc
	INTEGER			:: 	LIMIT, IMO, IDY, FLUXTYPE, IT, ITB
	INTEGER			::	ITA, ITR, ADJ_FLAG, STABFUNC, TSTEPS
	REAL			:: 	TOLBISECT, TOLTR, TOLTAERO
	REAL			::	TSTEP, T0, Q1, IRAD, DELT, VISC, CAPG, GHFR
	REAL			::	OLDTAERO, stability, U1, V1, Z1, G, RIHALF
	REAL			::	ADJ_MA, LHV_T, OLDTR, TIME0, ANGLE, UFLUX, VFLUX
	REAL			::	oldftr, newftr, deltr, ghfrold
	REAL, DIMENSION(TSTEPS)	::	OBS_INSOL
	
	LOGICAL			::	CONVERGED
	INTEGER			::	ICNT, LIMIT_TA
	REAL			::	OLDL, TOLL
	REAL			::	Zilitinkevich
        REAL Tactual,Poo
	INCLUDE 'SFC.incl'
	INCLUDE 'SFC.DAT'


!----------SET THE TOLERANCE FOR THE VARIOUS ITERATIONS-------------------------

	LIMIT	  = 500
	TOLL      = (3.e-3)
	LIMIT_TA  = 10
	TOLTAERO  = (3.e-1)

!---------CALCULATE SURFACE UTILITY VARIABLES (VISCOCITY, HEAT CAPACITY)--------

	VISC    = (1.32+0.009*(TAERO-273.15))*1.e-5
	Poo =100000
!CHECK RHO        
        TActual = T1*(press/poo)**(rd/cp)
	RHO	= press/(RD*Tactual)
!  use constant density
        RHO=1.1
        
!	CHI     = SQRT(RHOS*CS*XLAMBDA) / (4.18*10000.)
!	CAPGB   = 3.293e6 * CHI
!	DELHS   = (XLAMBDA) / ( CAPGB*(1.18*(7.292e-5)) )
	CAPG    = (XLAMBDA) / DELHS
	Cb      = Rhos*cs*delhs
	G	= 9.81
	C1	= 0.45

!---------Calculate SWnet radiation---------------------------------------------

        CALL SWNET(LAT, IMO, IDY, IRAD, TSTEP, TIME, DELT, TSTEPS, XIS, OBS_INSOL)

!---------Calculate Incoming LW radiation---------------------------------------


        If (Iradswitch.EQ.0) then
		XIL = (QC+(1.-QC)*.67*((1670*Q1)**.08))*SIGMA*T1**4.
        else
	        Xil = rlongdnsfc
		
	End If 
                Xil=Xil+GHGADD
	

	
!------ Update Ground heat flux and outgoing long wave radiation ---------------

	GHF = Capg*(T0-THM)
	LWU = EMISS*SIGMA*TR**4

!---------Calculate Skin Temperature--------------------------------------------
!
! The following calculates all the temperature and fluxes at the 
! surface interface until a balanced is achieved.
!
! We calculate 3 temperatures:
!		Taero:	Aerodynamic temperature
!		TR:	Skin radiative temperature
!		T0:	Ground temperature (TG)
!
! We iterate computing these temperatures and the resulting
! surface fluxes (u*u*, u*T*) until the fluxes converge.
! Since Moni-Obukhov length (L) indicates the variation of T* and u*
!
!		L = (Tbar.u*^2)/(k.g.T*)
!
! then, convergence of L is an indication of convergence of
! surface fluxes.
!
! Thus the following block of code comprises two loops:
!
!	1) the outer loop, updates TAERO until TAERO converges
!	2) the inner loop, updates L until L converges
!
!--------------------------------------------------------------------------------
 
!--------------------------------------------------------------------------------
!
! First calculate TR by allowing delTR/delt = 0, using bisection method
!
!--------------------------------------------------------------------------------
!
!
!------ START OF OUTER LOOP (for TAERO to converge)
!
!	CALL GETL(L, Z1, USTAR, TSTAR, T1, VONKAR, G, TAERO)
        ITA=0
	ITR=0
9	CONTINUE


!--------------------------------------------------------------------------------
!
! Now calculate surface fluxes until Monin_Obukhov length (L) converges
!
!--------------------------------------------------------------------------------
!


!------ Check for convergence in L
!
        ICNT = 0
	stability = T1 - TAERO
	CONVERGED = .FALSE.
	DO WHILE (.NOT. CONVERGED)

	   IF (FLUXTYPE .EQ. 1) THEN		! Businger Form for BL Fluxes
		CALL SFCFLUX(L, USTAR, TSTAR, T1, Z1, Z0, VONKAR, TAERO, stability, U1, V1, T0, G)
	   ENDIF

	   IF (FLUXTYPE .EQ. 2) THEN		! England Form for Stable BL Fluxes

		IF(stability.LT.0.0)THEN

			CALL SFCFLUX(L, USTAR, TSTAR, T1, Z1, Z0, VONKAR, TAERO, stability, U1, V1, T0, G)
		ELSE
			CALL STABLEFLUX(T1, TAERO, U1, V1, Z1, Z0, USTAR, TSTAR, L, VONKAR, STABFUNC, RIHALF)		

		ENDIF
	   ENDIF

!
!----- Check for convergence in L
!
	   IF (L == 0.0) THEN
		IF (OLDL == 0.0) CONVERGED = .TRUE.
	   ELSE
		IF (ABS((L-OLDL)/L) .LT. TOLL) CONVERGED = .TRUE.
	   END IF

	   OLDL=L
           ICNT=ICNT+1
	   IF (ICNT .GT. LIMIT) THEN
		STOP'ICNT EXCEEDED LIMIT, L DID NOT CONVERGE'
	   END IF

	END DO


!----- Calculate QSTAR                        ---------------------------------

        CALL GETSFCSPECHUM(QA, TAERO, PRESS)
        CALL GETQSTAR(L, QSTAR, Z1, Q1, QA, Z0, VONKAR, T1-TAERO)
        poo=100000.
	TActual = T1*(press/poo)**(rd/cp)
	RHO	= Press/(RD*Tactual)
!       use constant density
	Rho=1.1	
!----- Calculate Sensible and Latent Heat Flux ---------------------------------

	HO 	 = -CP  * RHO    * USTAR * TSTAR
	LHV_T = (2.501 - (0.00237*(TAERO-273.15)) ) * 1e6
        IF(ADJ_FLAG .EQ. 1)THEN
		LH = -RHO * ADJ_MA * USTAR * QSTAR * LHV_T
        ELSE
        	LH = -RHO * MAVAIL * USTAR * QSTAR * LHV_T
        ENDIF

!--------------------------------------------------------------------------------
!
! Now calculate aerodynamic temperature by performing Zilitinkevich adjustment
!
!     Taero = Tz0 = TR + .0962 (T*/k) (u*.Z0/v)^.45
!
!--------------------------------------------------------------------------------
!

	Zilitinkevich = 0.0962*(TSTAR/VONKAR)*( (USTAR*Z0/VISC)**C1 )
	IF (ABS(Zilitinkevich) .LT. 1.E-3) THEN
	  TAERO = 0.5 * (TR+T1)
	ELSE
	  TAERO = TR + Zilitinkevich
	END IF
!
! Adjust TAERO to a reasonable temperature
!   if Zilitinkevich adjustment was too much
!
!	IF (TR .GT. T1) TAERO = MIN (MAX(TAERO, T1), TR)
!	IF (TR .LT. T1) TAERO = MAX (MIN(TAERO, T1), TR)

	IF (TR .GT. T1) THEN
	  IF (TAERO .GT. TR) THEN
	    TAERO = TR
	  ELSE
	    IF (TAERO .LT. T1) TAERO = T1
	  END IF
	ELSE
	  IF (TAERO .LT. TR) THEN
	    TAERO = TR
	  ELSE
	    IF (TAERO .GT. T1) TAERO = T1
	  END IF
	END IF

!	write(*,*) "TAERO, ustar, tstar: ", taero, ustar, tstar

!
!----- Check for convergence in TA 
!
	ITA = ITA + 1
	IF (ITA .LT. LIMIT_TA)THEN
		IF (ABS((TAERO - OLDTAERO)/TAERO) .GT. TOLTAERO) THEN
	   	  OLDTAERO = TAERO
		  GO TO 9
		END IF
	END IF
!        TAERO = 0.5 * (TAERO + OLDTAERO)


!------ Update Ground heat flux and outgoing long wave radiation ---------------

	GHF = Capg*(T0-THM)
	LWU = EMISS*SIGMA*TR**4


end subroutine flux




!================================================================================
!
! The following routine uses bisection root finding method to calculate TR
! and in the process calculates all fluxes
!
!================================================================================

!       1       DDD     BBBB    L       MM   MM
!      11       D  D    B  B    L       M  M  M
!       1       D  D    BBB     L       M  M  M
!       1       D  D    B  B    L       M     M
!     11111     DDD     BBBB    LLLLL   M     M
!

   SUBROUTINE GET_TR (TR, TAERO, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS, &
			XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, Z0,	&
			ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 	&
			STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFR, QSTAR)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine solves for the following variables:
!	TR		= Skin Temperature using root finding
!	TAERO		= Aerodynamic Temperature using Zilitinkevich(1970)
!	T0 (TG)		= Ground Temperature using prognostic sfc. energy eqn.
!	HO, LH		= Sensible, Latent Heat Flux
!	UFLUX, VFLUX	= w'u', w'v' (Pielke v1, p.140)
!-------------------------------------------------------------------------------
	INTEGER                 ::      IRadswitch
	Real 			::      rLongdnsfc
	INTEGER			:: 	LIMIT, IMO, IDY, FLUXTYPE, IT, ITB
	INTEGER			::	ITA, ITR, ADJ_FLAG, STABFUNC, TSTEPS
	REAL			:: 	TOLBISECT, TOLTR, TOLTAERO
	REAL			::	TSTEP, T0, Q1, IRAD, DELT, VISC, CAPG, GHFR
	REAL			::	OLDTAERO, stability, U1, V1, Z1, G, RIHALF
	REAL			::	ADJ_MA, LHV_T, OLDTR, TIME0, ANGLE, UFLUX, VFLUX
	REAL			::	oldftr, newftr, deltr, ghfrold
	REAL, DIMENSION(TSTEPS)	::	OBS_INSOL
	
	INTEGER			::	NINC
	REAL			::	STEP
	REAL			::	GHFRA, GHFRB, FTRA, FTRB, TRA, TRB
	REAL			::	TAEROA, TAEROB
	REAL			::	TRX, GHFRX, FTRX, TAEROX

	LOGICAL, SAVE		::	LFIRST=.TRUE.
	LOGICAL			::	CONVERGED
	INTEGER			::	ICNT
	REAL			::	OLDL, LTOL, TOLFTR
        Real   Tactual,poo
	INCLUDE 'SFC.incl'
	INCLUDE 'SFC.DAT'

!----------SET THE TOLERANCE FOR THE VARIOUS ITERATIONS-------------------------

	TOLBISECT = (3.e-7)	
	TOLTR	  = (1.e-5)
	TOLTAERO  = (3.e-5)
	TOLFTR	  = (1.e-1)
	LTOL      = (3.e-3)

!---------CALCULATE SURFACE UTILITY VARIABLES (VISCOCITY, HEAT CAPACITY)--------
        Poo=100000.
	VISC    = (1.32+0.009*(TAERO-273.15))*1.e-5
        TActual = T1*(press/poo)**(rd/cp)
	RHO	= Press/(RD*Tactual)
!  use contant density
        Rho=1.1
		
!       CHI     = SQRT(RHOS*CS*XLAMBDA) / (4.18*10000.)
!	CAPGB   = 3.293e6 * CHI
!	DELHS   = (XLAMBDA) / ( CAPGB*(1.18*(7.292e-5)) )
	CAPG    = (XLAMBDA) / DELHS
	Cb      = Rhos*cs*delhs
	
	G	= 9.81
	C1	= 0.45


!--------- USING BISECTION METHOD TO FIND TR

	NINC    = 0
	IF (LFIRST) THEN
	  LFIRST =.FALSE.
	  N	 = 30.
	  LIMIT	 = 120
	ELSE
	  N	 = 20.
	  LIMIT	 = 100
	END IF
	
	
3	CONTINUE

	SIGMA	= 5.67e-8
	ITR	= 0
	ITA	= 0
	STEP	= 2.0 * N / float(LIMIT)
	TRA	= TR - N
	TAEROA	= TRA
	TRB	= TRA + STEP
	TAEROB	= TRB

!
! Step through until zero line is crossed
!

	GHFRA	= ( CAPG*(TRA-T0) )
	CALL FLUX  (TRA, TAEROA, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS, &
                      XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, Z0, 	&
		      ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 		&
		      STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFRA, QSTAR,  &
		      rlongdnsfc,Iradswitch)
	FTRA	= ( XIS + XIL - (EMISS*SIGMA*TRA**4.) - HO - LH - GHFRA )
        Write (*,*) 'A FTRA, XIL, H0' , FTRA, XIL, HO 

4	CONTINUE

	GHFRB	= ( CAPG*(TRB-T0) )
	CALL FLUX  (TRB, TAEROB, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS, &
                      XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, Z0, 	&
		      ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 		&
		      STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFRB, QSTAR,  &
		      rlongdnsfc,Iradswitch)
	FTRB	= ( XIS + XIL - (EMISS*SIGMA*TRB**4.) - HO - LH - GHFRB )
        Write (*,*) 'B  FTRB, XIL, H0' ,FTRB, XIL, HO                                                 
	 
	IF (FTRA*FTRB .GT. 0.0)THEN		! Continue stepping

		TRA     = TRB
		TAEROA  = TAEROB
		GHFRA	= GHFRB
		FTRA	= FTRB

		TRB     = TRA + STEP
		TAEROB  = TRB

		ITA = ITA + 1
		IF(ITA .GT. LIMIT)THEN
			N = 2.*N
			LIMIT = 2*LIMIT
                        NINC=NINC+1 
			IF (NINC .GT. 2 ) THEN
			  Write(*,*)  'T(1)', T1
			  WRITE(*,*) 'BRACKET FOR CONVERGENCE: ', TR-.5*N, TR+.5*N
			  WRITE(*,*) 'TR NOT CONVERGING in BRACKET'
			  WRITE(*,*) 'TRA, TRB, STEP: ', TRA, TRB, STEP
			  WRITE(*,*) 'TR, T0: ', TR, T0
			  write(*,*) 'XIS, XIL, H0, LH: ',XIS, XIL, -HO, -LH
			  write(*,*) 'LWA: ',- (EMISS*SIGMA*TRA**4.)
			  write(*,*) 'GHFA, GHFB, FTRA, FTRB: ',-GHFRA, -GHFRB, FTRA, FTRB
			  STOP
			ELSE
!			  WRITE(*,*)'INCREASING N TO ', N
			  GO TO 3
			END IF
		ENDIF
		GO TO 4				! Continue stepping
	ENDIF					! TRA and TRB surround Root

!
! Zero line was crossed, now zoom-in on the root
! 

	ITR=0
	LIMIT = 500
	CONVERGED = .FALSE.
	DO WHILE (.NOT. CONVERGED)
	  TRX    = 0.5*(TRA+TRB)
	  TAEROX = TRX
	  GHFRX	 = ( CAPG*(TRX-T0) )
	  CALL FLUX (TRX, TAEROX, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS, &
                          XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, Z0,  &
		        ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 	 &
		        STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFRX, QSTAR,&
			rlongdnsfc,Iradswitch )
	  FTRX	= ( XIS + XIL - (EMISS*SIGMA*TRX**4.) - HO - LH - GHFRX )

	  IF (ABS(FTRX) .LT. TOLFTR .OR. ITR .GT. LIMIT) THEN 		! CONVERGED
	    CONVERGED = .TRUE.
	  ELSE
	    ITR = ITR + 1
	    IF ( FTRX*FTRB .LT. 0 ) THEN
	    
		TRA     = TRX
		TAEROA  = TAEROX
		FTRA    = FTRX

	    ELSE					! Close enough, check with Tolerance
		
	        TRB     = TRX
		TAEROB  = TAEROX
		FTRB    = FTRX
	
	    ENDIF
	  END IF



	END DO

!------ CONVERGED, SET TR and TAERO

	TR    = TRX
	TAERO = TAEROX
!	TR    = 0.5*(TRA+TRB)
!	TAERO = 0.5*(TAEROA+TAEROB)

!------ Update Ground heat flux and outgoing long wave radiation ---------------

	GHF = Capg*(T0-THM)
	LWU = EMISS*SIGMA*TR**4

!	WRITE (*,*) "RETURNING FROM GET_TR"
!	WRITE (*,*) "TR, TAERO, FTR: ", TR, TAERO, FTRA, FTRB
!	WRITE (*,*) "XIS, XIL, H0, GHF: ", XIS, XIL, -HO, GHFRB 
	
	

end subroutine GET_TR






	SUBROUTINE SURFACEINIT (USTAR, TSTAR, QSTAR, L, HO, LH, TIME, T0, TR, TAERO,	&
                        	T1, Q1, IRAD, LAT, IMO, IDY, TSTEP, DELT, TSTEPS,	&
				OBS_INSOL, U1, V1, Z1, Z0, QA, FLUXTYPE, STABFUNC, 	&
				RIHALF, XIL,cb,rlongdnsfc,Iradswitch)

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine initializes the three temeprature system used
!	in 1DBLM, as well as the initial similarity fluxes, ustar, tstar, 
!	and the Monin Length, L, using the initial state estimate of the
!	surface energy budget.
!-------------------------------------------------------------------------------

	INCLUDE 'SFC.incl'
	Integer                 ::      Iradswitch
	Real                    ::      rlongdnsfc
	INTEGER			::	IRAD, IMO, IDY, IT, ITA, LIMIT, TSTEPS
	INTEGER			::	FLUXTYPE, STABFUNC
	REAL			::	TSTEP, DELT
	REAL			::	TOL, TOLTR, VISC, Q1, U1, V1, Z1
	REAL			::	OLDTAERO, OLDL, stability, RIHALF
	REAL, DIMENSION(TSTEPS)	::	OBS_INSOL
	REAL			::	T0	! ground temperature (TG)
	REAL			::	OLDTR, OLDSTABILITY
	REAL			::	CAPG, G


	INTEGER			::	ADJ_FLAG
	REAL			::	TIME0


	REAL			::	UFLUX, VFLUX, GHFR
        REAL poo, Tactual
	INCLUDE 'SFC.DAT'
	

	IF (IRAD .EQ. 1) THEN
		WRITE(18,*)'Using Observed Solar Insolation...'
	ENDIF
	IF (FLUXTYPE .EQ. 2)THEN
		WRITE(18,'A50')'Using Ri Formulation for Stable BL Fluxes...'
	ELSE
		WRITE(18,'A50')'Using Businger formulations for BL Fluxes...'
	ENDIF

!---------CALCULATE SURFACE UTILITY VARIABLES (VISCOCITY, HEAT CAPACITY)--------
        Poo     =  100000
	VISC    = (1.32+0.009*(T1-273.15))*1.e-5
	TActual = T1*(press/poo)**(rd/cp)
	RHO	= Press/(RD*Tactual)
! use constant density
        Rho=1.1
	      
	
!	CHI     = SQRT(RHOS*CS*XLAMBDA) / (4.18*10000.)
!	CAPGB   = 3.293e6 * CHI
!	DELHS   = (XLAMBDA) / ( CAPGB*(1.18*(7.292e-5)) )
	CAPG    = (XLAMBDA) / DELHS
	Cb      = Rhos*cs*delhs
	
	G	= 9.81
	C1	= 0.45

	IT      = 0
	ITA     = 0
	LIMIT   = 1000
	TOL     = 1.e-3
	TOLTR   = 3.e-5	
	Write(18,'A25,F12.3,A25,F12.3 ') 'Soil Heat Capacity =,', CS, ',Soil Conductivity,', XLAMbDA
	Write(18,'A25, F12.3, A25,F12.3')'Soil Density,', RHOS, ',Slab Depth,', delhs
	Write(18,'A15,F12.4')'GHGGADD=,' ,GHGADD 
	Write (18,'A40,I4') 'Radiation on (1) off (0),', Iradswitch

!----- Calculate SWnet radiation------------------------------------------------

	CALL SWNET(LAT, IMO, IDY, IRAD, TSTEP, TIME, DELT, TSTEPS, XIS, OBS_INSOL)

!----- Calculate Incoming LW radiation------------------------------------------

1	If (Iradswitch.EQ.0) then
		XIL = (QC+(1.-QC)*.67*((1670*Q1)**.08))*SIGMA*T1**4.
        else
	        Xil = rlongdnsfc
		
	End If 
                Xil=Xil+GHGADD
        
        

	
!----- Calculate Skin and Aerodynamic Temperature ------------------------------

				
	
	TR    = T0
	TAERO = TR
        Write (*,*) 'Before call to get _TR  T1' , T1
	CALL GET_TR (TR, TAERO, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS,    	&
                        XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, 	&
			Z0, ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 		&
			STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFR, QSTAR)

	
	end subroutine surfaceinit







!       1       DDD     BBBB    L       MM   MM
!      11       D  D    B  B    L       M  M  M
!       1       D  D    BBB     L       M  M  M
!       1       D  D    B  B    L       M     M
!     11111     DDD     BBBB    LLLLL   M     M
!

   SUBROUTINE GET_FLUX (TR, TAERO, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS, 	&
			XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, Z0,	&
			ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 	&
			STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFR,        &
	
			QSTAR,rlongdnsfc,Iradswitch )

IMPLICIT NONE
!-------------------------------------------------------------------------------
!	This subroutine solves for the following variables:
!	TR		= Skin Temperature using root finding
!	TAERO		= Aerodynamic Temperature using Zilitinkevich(1970)
!	T0 (TG)		= Ground Temperature using prognostic sfc. energy eqn.
!	HO, LH		= Sensible, Latent Heat Flux
!	UFLUX, VFLUX	= w'u', w'v' (Pielke v1, p.140)
!-------------------------------------------------------------------------------
        Integer                 ::      Iradswitch
	Real                    ::      rlongdnsfc
	INTEGER			:: 	LIMIT, IMO, IDY, FLUXTYPE, IT, ITB
	INTEGER			::	ITA, ITR, ADJ_FLAG, STABFUNC, TSTEPS
	REAL			:: 	TOLBISECT, TOLTR, TOLTAERO
	REAL			::	TSTEP, T0, Q1, IRAD, DELT, VISC, CAPG, GHFR
	REAL			::	OLDTAERO, stability, U1, V1, Z1, G, RIHALF
	REAL			::	ADJ_MA, LHV_T, OLDTR, TIME0, ANGLE, UFLUX, VFLUX
	REAL			::	oldftr, newftr, deltr, ghfrold
	REAL, DIMENSION(TSTEPS)	::	OBS_INSOL
	
	LOGICAL, SAVE		::	LFIRST=.TRUE.
	LOGICAL			::	CONVERGED
	INTEGER			::	ICNT
	REAL			::	OLDL, LTOL
	Real Tactual,Poo

	INCLUDE 'SFC.incl'
	INCLUDE 'SFC.DAT'

!----------SET THE TOLERANCE FOR THE VARIOUS ITERATIONS-------------------------

	LIMIT	  = 500
	TOLBISECT = (3.e-7)	
	TOLTR	  = (1.e-5)
	TOLTAERO  = (3.e-5)
	LTOL      = (3.e-3)

!---------CALCULATE SURFACE UTILITY VARIABLES (VISCOCITY, HEAT CAPACITY)--------
        Poo=100000.
	VISC    = (1.32+0.009*(TAERO-273.15))*1.e-5
        TActual = T1*(press/poo)**(rd/cp)
	RHO	= Press/(RD*Tactual)
!use constant density
        Rho=1.1
	
!	CHI     = SQRT(RHOS*CS*XLAMBDA) / (4.18*10000.)
!	CAPGB   = 3.293e6 * CHI
!	DELHS   = (XLAMBDA) / ( CAPGB*(1.18*(7.292e-5)) )
	CAPG    = (XLAMBDA) / DELHS
	Cb      = Rhos*cs*delhs

	
	G	= 9.81
	C1	= 0.45
!---------Calculate SWnet radiation---------------------------------------------

        CALL SWNET(LAT, IMO, IDY, IRAD, TSTEP, TIME, DELT, TSTEPS, XIS, OBS_INSOL)

!---------Calculate Incoming LW radiation---------------------------------------


        If (Iradswitch.EQ.0) then
		XIL = (QC+(1.-QC)*.67*((1670*Q1)**.08))*SIGMA*T1**4.
        else
	        Xil = rlongdnsfc
	End if
        	Xil=Xil+GHGADD
        
        
	
	
!---------Calculate Skin Temperature--------------------------------------------
!
! The following calculates all the temperature and fluxes at the 
! surface interface until a balanced is achieved.
!
! We calculate 3 temperatures:
!		Taero:	Aerodynamic temperature
!		TR:	Skin radiative temperature
!		T0:	Ground temperature (TG)
!
! We iterate computing these temperatures and the resulting
! surface fluxes (u*u*, u*T*) until the fluxes converge.
! Since Moni-Obukhov length (L) indicates the variation of T* and u*
!
!		L = (Tbar.u*^2)/(k.g.T*)
!
! then, convergence of L is an indication of convergence of
! surface fluxes.
!
! Thus the following subroutine, GET_TR, comprises two loops:
!
!	1) the outer loop, updates TR until TR converges
!	2) the inner loop, updates L until L converges
!
!--------------------------------------------------------------------------------
 
!--------------------------------------------------------------------------------
!
! First calculate TR by allowing delTR/delt = 0, using bisection method
!
!--------------------------------------------------------------------------------
!
	CALL GET_TR (TR, TAERO, T0, T1, Q1, LAT, IMO, IDY, IRAD, DELT, TSTEPS,  &
                        XIS, XIL, OBS_INSOL, TSTEP, FLUXTYPE, U1, V1, Z1, 	&
			Z0, ADJ_FLAG, TIME0, HO, LH, L, TSTAR, USTAR, TIME, 	&
			STABFUNC, RIHALF, UFLUX, VFLUX, LWU, GHF, GHFR, QSTAR)


!------ Update Ground Temperature, TG-------------------------------------------

	GHF = Capg*(T0-THM)
	LWU = EMISS*SIGMA*TR**4
	T0 = T0 + (DELT) * ( (XIS+XIL-LWU-HO-LH-GHF)/Cb )

!------ Calculate w'u' and w'v' using eqn 7.24 in Pielke p 140------------------

	ANGLE = ATAN2(V1,U1)
	UFLUX = USTAR * COS(ANGLE)
	VFLUX = USTAR * SIN(ANGLE)

end subroutine get_flux

Subroutine RADIATION (Delt,P,T,Q,RHOP,ozone,cosz,   &
            alb,em,Tr,fthrd,rlongdnsfc,raduptop,raddntop)
! Calculates Radiative Heating in the Atmosphere
!  Using call to FLIOU radiative driver

Implicit None
 Real  cb,cosz,alb,em,rlongdnsfc
	
include 'BL.incl'
Real oldtheta(kz)
   

   
   
   
   Call oned_Raddrvr(p,t,rhop,q,ozone,cosz,   &
          alb,em,tr,fthrd,rlongdnsfc,raduptop,raddntop,delt)
   Do i =1,kz-1
   T(i)=T(i)+fthrd(i)*delt
   end do
   
  
   

end Subroutine Radiation

Subroutine InitRADIATION (Delt,P,T,Q,RHOP,ozone,cosz, alb,em,Tr,rlongdnsfc)

! Calculates Radiative Heating in the Atmosphere
!  Using call to FLIOU radiative driver
Implicit None
 Real  cb,cosz,alb,em,rlongdnsfc	
include 'BL.incl'
 Double Precision oldtheta(kz)


   Do I =1,Kz
   oldtheta(i)=t(i)
   end do
   
    Call oned_Raddrvr(p,t,rhop,q,ozone,cosz,   &
    alb,em,tr,fthrd,rlongdnsfc,raduptop,raddntop)
!   discard heating rates on initial call   
       Do i =1,kz
   t(i)=oldtheta(i)
   end do
   

end Subroutine InitRadiation 

Subroutine ColumnEnergycalc (FirstCE,T,p,Columnenergy,InitialcolumnE,Radiativeloss,fthrd,delt,zm)
Implicit None
include 'BL.incl'
Real cp,rhoa,Rd,Tactual

ColumnEnergy=0.
	CP	=	1004.
	Press=1000*100
	RD	=	287.

	Columnenergy=0.0
	Radiativeloss=0.0
	Do j= 1,85
! sum energy  #################CORRECT ENERGY

        TActual=T(j)*(p(j)/press)**(rd/cp)
	RHOA	= P(j)/(RD*Tactual)
        Rhoa=1.1
	 ColumnEnergy= ColumnEnergy + rhoA*Cp*(.5*(T(j+1)+T(j)))*(zm(j+1)-zm(j))
	Radiativeloss=Radiativeloss+ rhoA*Cp*fthrd(j)*(zm(j+1)-zm(j))
	
	End Do
	
	If (FirstCE.EQ.1) then
	InitialColumnE=ColumnEnergy
	FirstCE=0
	End If
	 write(*,*) 'Columnenergy', columnenergy,InitialcolumnE
	Columnenergy=columnenergy-InitialColumnE
	
       
End Subroutine ColumnEnergyCalc

	Subroutine SumFlux (Nstepsavg,Tr,TrAVG,TAERO, TAEROAVG, &
	        T0,T0AVG,XIS,XISAVG,XIL,XILAVG,HO,HOAVG,LH,LHAVG,LWU,LWUAVG, &
	        GHF,GHFAVG,ColumnEnergy,Columnenergyavg,radiativeloss,radiativelossavg)
	Implicit None
	Include 'BL.incl'
	
!       Sum components
	Nstepsavg=Nstepsavg+1
	
	TRAVG=TRAVG+TR
	TAEROAVG=TAEROAVG+TAERO
!	REM Don't average
	T0AVG=T0
	XISAVG=XISAVG+XIS
	XILAVG=XILAVG+XIL
	HOAVG= HOAVG+HO
	LHAVG=LHAVG+LH
	LWUAVG=LWUAVG+LWU
	GHFAVG=GHFAVG+GHF
!	REM Don't average 
	Columnenergyavg=columnenergy
	radiativelossavg=radiativelossavg+radiativeloss
	
	End Subroutine SumFlux
	
	Subroutine SumHeatingProf(fthrd,  &
	dfRate,Numheatingavg,Radheating,Diffheating)
	Implicit None
	Include 'BL.incl'

	Do j=1,kz
	RadHeating(j)=Radheating(j)+fthrd(j)
	DifFheating(j)=Diffheating(j)+dfrate(j)
	End Do

	Numheatingavg=Numheatingavg+1
	End Subroutine SumHeatingProf
	
	
	
	
	Subroutine Pressure (P,T,zm)
	Implicit None
	Include 'BL.incl'
	Real cp,rhoa,Rd,Thetaavg
	Integer kzm
	Real Pip(kz)
	

	G=9.81
	CP	=	1004.
	Press=1000*100
	RD	=	287.
	
	Pip(1)=cp*(p(1)/Press)**(rd/cp)
	
	Kzm=kz-1
	Do J=1,Kzm

	Thetaavg=(.5*(t(j)+t(j+1)))
	PiP(j+1)= Pip(j)-(g/thetaavg)*(zm(j+1)-zm(j))
	P(j+1)= Press*(pip(j+1)/cp)**(Cp/rd)
! 	Write(*,*) 'Pressure',j+1, pip(j+1), p(j+1)
	End Do
	
	
	
	
	
End subroutine Pressure
	



