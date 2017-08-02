!----------------------------------------------------------------------!
! SUBROUTINE MODELS                                                    !
! Author(s): Daniel van Kraalingen  (original)                         !
! Date     : 5-Jul-1993, Version: 1.1                                  !
!          Version august, 503                                        !
!                                                                      !
! This version used in ORYZA2000 FSEWin model;                         !
! Date     : November 2002; Adapted by B.A.M. Bouman                   !
! Purpose  : This subroutine is the interface routine between the FSE  !
!            driver and the simulation models. This routine is called  !
!            by the FSE driver at each new task at each time step. It  !
!            can be used by the user to specify calls to the different !
!            models that have to be simulated.                         !
!                                                                      !
! FORMAL PARAMETERS:  I=input, O=output                                !
! name   type meaning (unit)                                     class !
! ----   ---- ---------------                                    ----- !
! ITASK   I4  Task that subroutine should perform (-)               I  !
! IUNITD  I4  Unit number that is used for input files (-)          I  !
! IUNITO  I4  Unit number that is used for output file (-)          I  !
! IUNITL  I4  Unit number that is used for log file (-)             I  !
! FILEIT  C*  Name of timer file (-)                                I  !
! FILEI1  C*  Name of input file no. 1 (-)                          I  !
! FILEI2  C*  Name of input file no. 2 (-)                          I  !
! FILEI3  C*  Name of input file no. 3 (-)                          I  !
! FILEI4  C*  Name of input file no. 4 (-)                          I  !
! FILEI5  C*  Name of input file no. 5 (-)                          I  !
! OUTPUT  L4  Flag to indicate if output should be done (-)         I  !
! TERMNL  L4  Flag to indicate if simulation is to stop (-)        I/O !
! DOY     R4  Day number since 1 January (day of year) (d)          I  !
! IDOY    I4  Day number within year of simulation (d)              I  !
! YEAR    R4  Year of simulation (y)                                I  !
! IYEAR   I4  Year of simulation (y)                                I  !
! STTIME  R4  Start day of simulation (d)                           I  !
! TIME    R4  Time of simulation (d)                                I  !
! DELT    R4  Time interval of integration (d)                      I  !
! LAT     R4  Latitude of site (dec.degr.)                          I  !
! WSTAT   C*  Status code from weather system (-)                   I  !
! WTRTER  L4  Flag whether weather can be used by model (-)         O  !
! RDD     R4  Daily shortwave radiation (kJ.m-2.d-1)                I  !
! TMMN    R4  Daily minimum temperature (degrees C)                 I  !
! TMMX    R4  Daily maximum temperature (degrees C)                 I  !
! VP      R4  Early morning vapour pressure (kPa)                   I  !
! WN      R4  Average wind speed (m.s-1)                            I  !
! RAIN    R4  Daily amount of rainfall (mm.d-1)                     I  !
!                                                                      !
! Subroutines called: ET, WSTRESS, WNOSTRESS, IRRIG, PADDY             !
!                     NSOIL, NNOSTRESS, NCROP, ORYZA1                  !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE MODELS (ITASK , IUNITD, IUNITO, IUNITL,                 &
                         FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5, &
                         OUTPUT, TERMNL,                                 &
                         DOY   , IDOY  , YEAR  , IYEAR , STTIME,         &
                         TIME  , DELT ,  LAT   , WSTAT , WTRTER,         &
                         RDD   , TMMN  , TMMX  , VP    , WN    , RAIN)

	  USE public_module		!VARIABLES
      USE Module_OutDat
      USE RADIATION_C
      
      IMPLICIT NONE

!-----Formal parameters
      INTEGER   ITASK, IUNITD, IUNITO, IUNITL, IDOY, IYEAR

      REAL      DOY  , YEAR  , TIME  , DELT  , LAT , STTIME
      REAL      RDD  , TMMN  , TMMX  , VP    , WN  , RAIN
      CHARACTER (*) FILEIT, FILEI1, FILEI2, FILEI3, &
                    FILEI4, FILEI5
      LOGICAL   OUTPUT, TERMNL, WTRTER
      CHARACTER (*) WSTAT
!-----Local variables
      INTEGER        I, I1, NL
      CHARACTER (6)  WUSED
      CHARACTER (80) ESTAB  , ETMOD   , PRODENV, RUNMODE, NITROENV, NUTRIENT, PRODENV1
      CHARACTER (80) WATBAL 
      INTEGER        CROPSTA, DTFSECMP, EMYR   , EMD    , IDATE , SBDUR, GERMDELY
      REAL           ANGA   , ANGB    , DAE    , DVS    , ETD   , EVSC
      REAL           FAOF   , IR      , LAI    , LDSTRS , LESTRS, LRSTRS
      REAL           PCEW   , RAINCU  , SWR    , TKLT   , TMDA  , TRC  
      REAL           TRW    , WL0     , ZRT    , ZRTMS	, NRT
      REAL           LAIROL , CPEW    ,SLA

      INTEGER      NLXM
      PARAMETER   (NLXM=10)
      REAL         MSKPA(NLXM), TKL(NLXM) , TRWL(NLXM)
      REAL         WCAD(NLXM) , WCWP(NLXM)  , WCFC(NLXM), WCST(NLXM)
      REAL         WCLQT(NLXM)
!     Parameters for climate change
      INTEGER IMXX 
      PARAMETER (IMXX = 800)
      INTEGER ILTMCT, ILTMAXCT, ILRADCT, ILRAINCT, ILVAPPCT, ILWINDCT, CCYEAR,GERMIN
      REAL TMCTB(IMXX), TMINCTB(IMXX), TMAXCTB(IMXX), TCOR,Y1,LINT2  
      REAL RADCTB(IMXX), RAINCTB(IMXX), VAPPCTB(IMXX), WINDCTB(IMXX)
      LOGICAL      GIVEN, GERMCAL,ISGERMDELY

      REAL LLV, DLDR, WLVG, WST, WSO, GSO, GGR, GST, GLV, GRT, PLTR
      REAL TNSOIL, NACR, NSLLV, NFLV, RNSTRS
      REAL ANSO, ANLV, ANST, ANLD, ANRT
!	  REAL IROOTD		!THE INITIAL ROOT DEPTH REGARDLESS OF PLANTING METHODS FROM EXPERIMENTAL FILE. TAOLI, APRIL 14, 509
      REAL TTEMPD,TTEMPN,TCHANG,SHOUR,TSYEAR, EHOUR, TEYEAR, SDAY, EDAY
	  REAL HU2, XTAV, XTMIN, XTMAX, XHU, XTAVD, TMAXO, TMINO
	  INTEGER ISTEMC ,CONTRM
	  LOGICAL RDINQR, TEMPC, ISSUBOLD
	  REAL		  SURRES		!THE THICKNESS OF SURFACE RESIDUE (mm), TAOLI, 16 JULY 2009
	  REAL SOILT(0:10)			!FOR SOIL TEMPERATURE, AND SURFACE TEMPERATURE SOILT(0), TAOLI, 21JULY 2009
      SAVE         !&#@TAOLI

      DATA WUSED /'------'/, GIVEN /.FALSE./

!     avoid compiler warnings for not using these variables
      FILEI4 = ' '
      FILEI5 = ' '

!==============================================================!
! Initialization section: ITASK = 1                            !
!==============================================================!
      IF (ITASK.EQ.1) THEN

!--------Read data from the experimental data file
         CALL RDINIT (IUNITD, IUNITL, FILEIT)
         CALL RDSCHA ('RUNMODE' , RUNMODE)
         CALL RDSCHA ('ESTAB'   , ESTAB  )
         CALL RDSCHA ('ETMOD'   , ETMOD  )
         IF(RDINQR('GERMIN')) THEN
            CALL RDSINT('GERMIN',GERMIN)
         ELSE
             GERMIN=0
         END IF
         CALL RDSCHA ('PRODENV' , PRODENV)
         PRODENV1 = PRODENV     !PRODENV1 is specifically the switch of water module under submergence and none-submergence condition
         CALL RDSCHA ('NITROENV' , NITROENV)
         CALL UPPERC (PRODENV)
         CALL UPPERC (NITROENV)
         CALL RDSCHA ('WATBAL' , WATBAL)          
         CALL UPPERC (WATBAL)        
         IF(INDEX(NITROENV,'NITROGEN BALANCE').GT.0) THEN
		    CALL RDSCHA ('NUTRIENT' , NUTRIENT)
		    CALL UPPERC (NUTRIENT)
		    pv%NBALANCE = .TRUE.
		 ELSE
		    pv%NBALANCE = .FALSE.
		 END IF 
         !--------TAOLI, June 10, 509
		 IF(RDINQR('ISTEMC')) THEN
			CALL RDSINT('ISTEMC',ISTEMC )
			IF(ISTEMC.GT.0) THEN		 
				TEMPC = .TRUE.
			ELSE
				TEMPC = .FALSE.
			ENDIF
			IF(TEMPC) THEN
			    CALL RDSREA('TTEMPN',TTEMPN )
                CALL RDSREA('TTEMPD',TTEMPD )
			    CALL RDSREA('TCHANG',TCHANG )
			    CALL RDSREA('SHOUR ',SHOUR  )
			    CALL RDSREA('EHOUR ',EHOUR  )
                CALL RDSREA('TSYEAR',TSYEAR )
			    CALL RDSREA('TEYEAR',TEYEAR )
			    CALL RDSREA('SDAY  ',SDAY   )
			    CALL RDSREA('EDAY  ',EDAY   )
			    CALL RDSINT('CONTRM',CONTRM )
			ENDIF
		 ENDIF
         CALL UPPERC (RUNMODE)
         CALL UPPERC (ESTAB  )
         CALL UPPERC (ETMOD  ) 	 
         CALL RDSREA ('ANGA' , ANGA)
         CALL RDSREA ('ANGB' , ANGB)
         IF (RUNMODE.EQ.'EXPERIMENT') THEN
            CALL RDSINT('EMYR' , EMYR)
            CALL RDSINT('EMD'  , EMD )
         ELSE IF (RUNMODE.EQ.'EXPLORATION') THEN
            IF(PV%KILLSOIL) THEN
                SWR = IYEAR
                EMYR = IYEAR
                EMD = NINT(STTIME)
			ELSE
			    CALL RDSINT('EMYR' , EMYR)
                CALL RDSINT('EMD'  , EMD )
			END IF
         ELSE
            CALL FATALERR  &
            ('MODELS','unknown name for RUNMODE')
         END IF
         CALL RDSINT('SBDUR  ',SBDUR)
         CALL RDSREA('FAOF'   ,FAOF )
!        Read climate change parameters in 
         ILTMCT = 0; ILTMAXCT = 0;ILRADCT=0; ILRAINCT=0; ILVAPPCT=0; ILWINDCT=0
		 IF(RDINQR('TMCTB')) THEN
			CALL RDAREA('TMCTB ',TMCTB,IMXX,ILTMCT)
		 ELSEIF(RDINQR('TMAXCTB')) THEN
			CALL RDAREA('TMINCTB ',TMINCTB,IMXX,ILTMAXCT)
			CALL RDAREA('TMAXCTB ',TMAXCTB,IMXX,ILTMAXCT)
         ENDIF
         IF(RDINQR('RADCTB')) THEN
			CALL RDAREA('RADCTB ',RADCTB,IMXX,ILRADCT)
         END IF
         IF(RDINQR('RAINCTB')) THEN
			CALL RDAREA('RAINCTB ',RAINCTB,IMXX,ILRAINCT)
         END IF
         IF(RDINQR('VAPPCTB')) THEN
			CALL RDAREA('VAPPCTB ',VAPPCTB,IMXX,ILVAPPCT)
         END IF
         IF(RDINQR('WINDCTB')) THEN
			CALL RDAREA('WINDCTB ',WINDCTB,IMXX,ILWINDCT)
         END IF
         IF(RDINQR('CCYEAR')) THEN
			CALL RDSINT('CCYEAR',CCYEAR)
         ELSE
             CCYEAR = IYEAR
         END IF
         !!--If the diffusion radiation change is also considered in the climate change
         IF(RDINQR('XFRDIF')) THEN
             CALL RDSINT('XFRDIF',XFRDIF)
             ILFRDIF = 0
             IF(XFRDIF.GT.0) THEN
                 IF(RDINQR('FRDIFCTB')) THEN
                    CALL RDAREA('FRDIFCTB ',FRDIFCTB,IMXX,ILFRDIF)
                 END IF
                 IF(ILFRDIF.LE.0) THEN 
                    CALL FATALERR ('MODELS','Did not provide the table for difussion radiation change!')
                 END IF
             END IF
         ELSE
             XFRDIF = 0
         END IF
!!-----------------------------------------------------------------------------------------
         CLOSE (IUNITD)

!--------Initialize variables
         CROPSTA  = 0; WL0 = 0.; RAINCU = 0.
         DO I=1,NLXM
            WCLQT(I) = 0.3
            WCST(I)  = 0.3
         END DO

!--------Write and check water production situation setting
         IF (PRODENV.EQ.'POTENTIAL') THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','Rice grown in potential water production situation'
         ELSE IF (PRODENV.EQ.'WATER BALANCE') THEN
!           Select water balance model
            IF (WATBAL.EQ.'PADDY') THEN
               WRITE (IUNITO,'(A,T7,A)')'*','Water balance PADDY used'
            ELSE IF (WATBAL.EQ.'SAWAH') THEN
               WRITE (IUNITO,'(A,T7,A)')'*','Water balance SAWAH used'
            ELSE IF (WATBAL.EQ.'SAHEL') THEN
               WRITE (IUNITO,'(A,T7,A)')'*','Water balance SAHEL used'
            ELSE IF (WATBAL.EQ.'LOWBAL') THEN
               WRITE (IUNITO,'(A,T7,A)')'*','Water balance LOWBAL used'
            ELSE IF (WATBAL.EQ.'SOILPF') THEN
               WRITE (IUNITO,'(A,T7,A)')'*','Water balance SOILPF used'
            ELSE
               CALL FATALERR  &
                  ('MODELS','unknown name for soil water balance')
            END IF
            WUSED(6:6) = 'U'
         ELSE
            CALL FATALERR  &
               ('MODELS','unknown name for production situation')
         END IF

!--------Choose and check info on nitrogen production situation setting
         IF (NITROENV.EQ.'POTENTIAL') THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','Rice grown in potential N production situation'
         ELSE IF (NITROENV.EQ.'NITROGEN BALANCE') THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','Crop and soil nitrogen balance used'
         ELSE
            CALL FATALERR ('MODELS','unknown name for NITROENV')
         END IF

!--------Send warning if water and nitrogen limitations are combined
         IF (NITROENV.EQ.'NITROGEN BALANCE' .AND. PRODENV.EQ.'WATER BALANCE') THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','****************************************************************'
            WRITE (IUNITO,'(A,T7,A)') &
              '*','WARNING: Combined water and nitrogen limitations not validated!!'
            WRITE (IUNITO,'(A,T7,A)') &
              '*','****************************************************************'
         END IF

!--------Write information about RUNMODE to output file
         IF (RUNMODE.EQ.'EXPERIMENT') THEN
             WRITE (IUNITO,'(A,T7,A)') &
             '*','ORYZA model runs to simulate an experiment'
         ELSE IF (RUNMODE.EQ.'EXPLORATION') THEN
             WRITE (IUNITO,'(A,T7,A)') &
              '*','ORYZA model runs for exploration'
         ELSE
            CALL FATALERR  &
               ('MODELS','unknown name for RUNMODE')
         END IF

!--------Choose and check establishment setting
         IF (ESTAB.EQ.'TRANSPLANT') THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','Rice crop is transplanted'
         ELSE IF (ESTAB.EQ.'DIRECT-SEED') THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','Rice crop is direct-seeded'
!BB: set SBDUR to 0. if direct seeded:
            SBDUR = 0.
         ELSE
            CALL FATALERR &
               ('MODELS','unknown name for establishment')
         END IF

!--------Choose and check evapotranspiration modules
         IF (INDEX(ETMOD,'PENMAN').GT.0) THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','SETPMD: Penman evapotranspiration'
            WUSED(1:5) = 'UUUUU'
         ELSE IF (INDEX(ETMOD,'MAKKINK').GT.0) THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','SETMKD: Makkink evapotranspiration'
            WUSED(1:5) = 'UUU--'
         ELSE IF (INDEX(ETMOD,'PRIESTLY TAYLOR').GT.0) THEN
            WRITE (IUNITO,'(A,T7,A)') &
              '*','SETPTD: Priestley Taylor evapotranspiration'
            WUSED(1:5) = 'UUU--'
         ELSE
            CALL FATALERR &
               ('MODELS','unknown module name for evapotranspiration')
         END IF

!--------Check weather data for ORYZA crop model
         WUSED(1:3) = 'UUU'


!--------Write log messages to output file
         WRITE (IUNITO,'(A,76A1)') '*',('=',I1=1,76)
         WRITE (IUNITO,'(A)') '*'
         WRITE (IUNITO,'(A)') '* ORYYZA200 driver info:'
         WRITE (IUNITO,'(A,T7,A,I5,A,I4,A)') &
           '*','Year:',IYEAR,', day:',IDOY,', System start'
         PV%GERMRATE=0.0;PV%GERMDAYS=0.0;GERMCAL=.FALSE.;PV%TOOWET=0;PV%TOODRY=0; GERMDELY=0; ISGERMDELY=.TRUE.
      END IF

!==============================================================!
! Here ended initialization section ITASK = 1                  !
!==============================================================!

!-----Check weather data availability
      IF (ITASK.EQ.1.OR.ITASK.EQ.2.OR.ITASK.EQ.4) THEN
         IF (WSTAT(6:6).EQ.'4') THEN
            RAIN       = 0.
            WSTAT(6:6) = '1'
            IF (.NOT.GIVEN) THEN
               WRITE (IUNITL,'(2A)') ' Rain not available,', &
                 ' value set to zero, (patch DvK, Jan 1995)'
               GIVEN = .TRUE.
            END IF
         END IF
!--------Check whether there is an error in the I1-th weather variable
         DO I1=1,6
            IF (WUSED(I1:I1).EQ.'U' .AND. &
                WSTAT(I1:I1).EQ.'4') THEN
               WTRTER = .TRUE.
               TERMNL = .TRUE.
               RETURN
            END IF
         END DO
      END IF
 !------Climate change PARAMETERS, THIS SECTION IS MOVED FROM SUBROUTINE ORYZA1, TAOLI 21SEP 2013
      IF(ITASK.EQ.2) THEN
         !------Climate change PARAMETERS, THIS SECTION IS MOVED FROM SUBROUTINE ORYZA1, TAOLI 21SEP 2013
        IF((ILTMCT.GT.0).and.(ILTMAXCT.LE.0)) THEN
			TCOR = LINT2('TMCTB',TMCTB,ILTMCT,DOY)
			TMMX = TMMX+TCOR
			TMMN = TMMN+TCOR
		ELSEIF((ILTMAXCT.GT.0).AND.(ILTMCT.LE.0)) THEN	!Added by TaoLi, 13 Oct. 2010
			TCOR = LINT2('TMAXCTB',TMAXCTB,ILTMAXCT,DOY)
			TMMX = TMMX+TCOR
			TCOR = LINT2('TMINCTB',TMINCTB,ILTMAXCT,DOY)
			TMMN = TMMN+TCOR		!Added by TaoLi, 13 Oct. 2010
        ENDIF
            !Protect for ensuring TMMX>TMMN
            IF(TMMN.GT.TMMX) THEN
                Y1 = TMMN;TMMN = TMMX;TMMX = Y1;Y1=0.0
            END IF
        Y1=MAX(0, -CCYEAR+IYEAR+1)
        IF(ILRAINCT.GT.0) THEN
            TCOR = LINT2('RAINCTB',RAINCTB,ILRAINCT,DOY)
            RAIN=RAIN*(1.0+TCOR/100.0)**Y1
        END IF
        IF(ILRADCT.GT.0) THEN
            TCOR = LINT2('RADCTB',RADCTB,ILRADCT,DOY)
            RDD=RDD*(1.0+TCOR/100.0)**Y1
        END IF
        IF((ILFRDIF.GT.0).AND.(XFRDIF.GT.0)) THEN
            FRDIFVC = LINT2('FRDIFCTB',FRDIFCTB,ILFRDIF,DOY) !Get the given value of change
            IF(XFRDIF.EQ.1) THEN                   !If change is in percentage year by year
                FRDIFVC = (1.0+FRDIFVC/100.0)**Y1  !Calculated change multiplication number
            END IF  
        ELSE
            XFRDIF = 0
        END IF
        IF(ILVAPPCT.GT.0) THEN
            TCOR = LINT2('VAPPCTB',VAPPCTB,ILVAPPCT,DOY)
            VP=VP*(1.0+TCOR/100.0)**Y1
        END IF
        IF(ILWINDCT.GT.0) THEN
            TCOR = LINT2('WINDCTB',WINDCTB,ILWINDCT,DOY)
            WN=WN*(1.0+TCOR/100.0)**Y1
        END IF
!----------To exchange the max and min temperature if TMIN>TMAX to ensure TMAX>TMIN when climate change, TAOLI, 21Sept 2013  
        IF(TMMN.GT.TMMX) THEN
            Y1 = TMMN;TMMN = TMMX;TMMX = Y1;Y1=0.0
        END IF
        IF(((TMMN+TMMX)*0.5.LT.12.0).AND.(ITASK.EQ.2)) THEN
!            pv%TCOMPEN = 12.0-(TMMN+TMMX)*0.5   !reset back to 12.
!            TMMN = TMMN+pv%TCOMPEN
!            TMMX=TMMX+pv%TCOMPEN
!        else
            pv%TCOMPEN = 0.0
        END IF
!-----Calculate average temperature
      TMDA = (TMMX+TMMN)/2.;TMAXO=TMMX;TMINO=TMMN
!-----added by TAOLI, 18 June 2009. If artificial temperature management is applied
		IF((TEMPC).AND.((Year*1000+DOY).GE.(TSYEAR*1000+SDAY)).AND.&
            ((YEAR*1000+DOY).LE.(TEYEAR*1000+EDAY))) THEN
			CALL TSHIFT(TMMX,TMMN,12.0, SHOUR, EHOUR, TTEMPD,TTEMPN,TCHANG, &
				ISTEMC,CONTRM, 8.0,30.0,42.0,HU2, XTAV, XTMIN, XTMAX, XHU, XTAVD)
		ELSE
			XTAV=0.0; XTMIN=0.0; XTMAX=0.0; XHU=0.0; XTAVD=0.0
		ENDIF
	   TMDA = TMDA + XTAV;TMMX=TMMX+XTMAX;TMMN=TMMN+XTMIN
    !------end section
      END IF
!-----Calculate potential soil evaporation and transpiration
       CALL ET2(ITASK,OUTPUT, ANGA,  ANGB,  RDD,   TMDA,    VP,  WN,   LAT, &
                    IDOY,  DELT, ETMOD, CROPSTA, ESTAB, NL,  FAOF, WL0, &
                    WCLQT, WCST,  LAI,   EVSC,    ETD, TRC)
	!-----Calculate drought stress factors
      IF (PRODENV.EQ.'WATER BALANCE') THEN
	     CALL WSTRESS2 (ITASK,  DELT,   OUTPUT, IUNITD, IUNITL, FILEI1, FILEIT, &
                       TRC,    ZRT,    TKL,    NL,    CROPSTA, &
                       WCLQT,  WCWP,   MSKPA,                 &
                       TRW,    TRWL,   LRSTRS, LDSTRS, LESTRS, PCEW, CPEW)
      ELSE IF (PRODENV.EQ.'POTENTIAL') THEN
	     IF((ITASK.EQ.1).AND.(NITROENV.EQ.'NITROGEN BALANCE')) THEN        !READ in soil layers, thickness information, surface water
	        IF(LEN_TRIM(FILEI2).GT.0) THEN
	            GOTO 1000
	        ELSE
	            CALL FATALERR ('MODELS',  &
                'Input file for soil is missing!')
	        END IF
	     ELSE
	        CALL POTENTIAL_SOIL(ITASK, IUNITD, IUNITL, FILEI2, OUTPUT, NITROENV, &
                      NL, TKL, TKLT, WCAD, WCWP, WCFC, WCST, WCLQT, WL0)  
         END IF
500         TRW = TRC; ZRTMS = TKLT
!!		 CALL WNOSTRESS(NL,TRW, TRWL, LRSTRS, LDSTRS, LESTRS, PCEW, CPEW)	   
		 CALL WNOSTRESS(NL,TRW,TRWL,ZRT,TKL,LRSTRS,LDSTRS,LESTRS,PCEW,CPEW)
      END IF
		 IF ((PRODENV.EQ.'WATER BALANCE').OR.(NITROENV.EQ.'NITROGEN BALANCE')) THEN
			 PV%PROOT_NUTRIENT =	.TRUE.
		 ELSE
			 PV%PROOT_NUTRIENT = .FALSE.
         ENDIF
   
!-----Call the crop growth module
        CALL ORYZA1(ITASK,  IUNITD, IUNITL, FILEI1, FILEI2,FILEIT, &
                        OUTPUT, TERMNL, IDOY  , DOY, IYEAR, YEAR,EMD, EMYR, &
                        TIME,   DELT,   LAT,    RDD,    TMINO,   TMAXO, &
                        NFLV,   NSLLV,  NRT, RNSTRS,                 &
                        ESTAB,  TKLT,   ZRTMS,  CROPSTA, &
                        LRSTRS, LDSTRS, LESTRS, PCEW,  CPEW, TRC, &
                        DAE,    SLA, LAI,    LAIROL, ZRT,    DVS, &
                        LLV,    DLDR, WLVG, WST, WSO, GSO, GGR, GST, GLV, GRT, &
                        PLTR,WCLQT, WL0)
        
!-----Call the nitrogen crop demand and soil supply modules
       IF (NITROENV.EQ.'NITROGEN BALANCE') THEN
           CALL NCROP3 (ITASK, IUNITD, IUNITL, FILEI1, FILEI2, FILEIT, DELT, TIME, OUTPUT, &
                       TERMNL, DVS, LLV, DLDR, WLVG, WST, WSO, GSO, GST, GLV, GRT, &
                       PLTR, LAI, SLA, CROPSTA, TNSOIL, NACR, ANSO, ANLV, ANST, ANLD, &
                       ANRT, NFLV, NSLLV,NRT, RNSTRS)                
      ELSE IF (NITROENV.EQ.'POTENTIAL') THEN
          CALL NNOSTRESS2(DELT, IUNITD, IUNITL, ITASK, FILEI1, FILEIT, &
                           CROPSTA, DVS, WLVG, LAI, SLA, NFLV, NSLLV, RNSTRS)
      END IF

!-----Call the water balance module
      IF (PRODENV.EQ.'WATER BALANCE') THEN
!--      First, the irrigation subroutine
         CALL IRRIG (ITASK, IUNITD, IUNITL,  FILEIT, OUTPUT, &
                        DOY,   DELT,   CROPSTA, WL0, DAE, &
                        DVS, NL,    WCLQT,  MSKPA,   IR)   
!        Then the soil water balans module
1000      IF (WATBAL.EQ.'PADDY') THEN
            CALL PADDY (ITASK, IUNITD, IUNITL, FILEI2, OUTPUT, &
                      DOY,    DELT,   TIME,   CROPSTA,  ESTAB, &
                      RAIN,   EVSC,   TRWL,   TRW,      IR, &
                      NL,     ZRTMS,  TKL,    TKLT, &
                      WCAD,   WCWP,   WCFC,   WCST,     WCLQT, &
                      WL0,    MSKPA)
         ELSE IF (WATBAL.EQ.'SAWAH') THEN
            CALL SAWAH  (ITASK, IUNITD, IUNITL, FILEI2, OUTPUT, &
                       ESTAB,  CROPSTA, &
                       DOY,   DELT,   TIME,   RAIN,   IR,  EVSC, &
                       TRC,    TRWL, NL,    ZRTMS,  TKL,    TKLT, &
                       WCAD,   WCWP,   WCFC,   WCST, WCLQT, &
                       WL0,    MSKPA)
         ELSE IF (WATBAL.EQ.'SAHEL') THEN
            CALL DRSAHE (ITASK , IUNITD, IUNITL, FILEI2 ,       &
                        DELT  , OUTPUT, TERMNL, ESTAB, CROPSTA,       &         
                        NL    , EVSC  , RAIN  , IR, ZRT,   &
                        TRWL  , TKL   , TKLT  , ZRTMS , WL0 ,  &
                        WCAD  , WCWP  , WCFC  , WCST ,         &
                        WCLQT , MSKPA)
         ELSE IF (WATBAL.EQ.'LOWBAL') THEN
            CALL LOWBAL (ITASK, IUNITD, IUNITO, &
                         FILEI2, OUTPUT, DELT, TIME, ESTAB, CROPSTA,IR,&
                         TRWL, EVSC, RAIN, NL, TKL, TKLT, ZRTMS,    &
                         WCWP, WCFC, WCST, WCLQT, WL0, MSKPA) 
         ELSE IF (WATBAL.EQ.'SOILPF') THEN
           CALL SOILPF(ITASK, IUNITD, IUNITO, FILEI2, OUTPUT, & 
                       DOY, DVS, NL, TKL, TKLT, ZRTMS, &
                       WCWP, WCFC, WCST, WCLQT, WL0, MSKPA)
         END IF
         IF((ITASK.EQ.1).AND.(PRODENV.EQ.'POTENTIAL')) GOTO 500
!---  No water balance in potential situation
      ELSE IF (PRODENV.EQ.'POTENTIAL') THEN
      
      END IF

	  IF ((PRODENV.EQ.'WATER BALANCE').OR.(NITROENV.EQ.'NITROGEN BALANCE')) THEN
	  !-----added by TAOLI, 16 July 2009, 
			 CALL SOILTEMP(ITASK,  IUNITD, IUNITL, FILEI1, FILEI2,FILEIT, &
                        OUTPUT, TERMNL, IDOY, DOY, DELT, TMMN, TMMX, max(0.0,LAI),EVSC+TRC, &
						ETD, DVS, WCLQT, WL0)
	  ENDIF
	  IF (NITROENV.EQ.'NITROGEN BALANCE') THEN
		IF(NUTRIENT .EQ. 'GENERAL SOM') THEN
		     CALL SoilCarbonNitrogen(ITASK, IUNITD, IUNITL, FILEI2, FILEIT,OUTPUT, &
                   DOY,    DELT,   DAE,   DVS, CROPSTA,  RAIN,  NL,     ZRTMS,  &
		           TKL,    WCAD,   WCWP,   WCFC,   WCST,  WL0,   WCLQT, pv%Psoiltx)
		ELSEIF(NUTRIENT.EQ.'APSIM SOILN') THEN
		     CALL SoilN2 (ITASK, IUNITD, IUNITL, FILEI1, FILEI2, FILEIT, OUTPUT, &
                   DOY,    DELT, DAE)		
			 CALL soiln2_pond (ITASK, IUNITD, IUNITL, FILEI1, FILEI2, FILEIT)
		ELSEIF(NUTRIENT.EQ.'FIXED SUPPLY') THEN     !added by TAOLI, 23June 2012
		     CALL NSOIL3(ITASK, IUNITD, IUNITL, FILEI2,FILEIT,OUTPUT, DELT, DAE, &
                    DVS, NACR, TKL, NL)		
		ENDIF
      ENDIF
 !!---Updait the accumulative of rain
      
!+-------------------------------------------------------------------------------------
!+	use the variables available here to update global variable values, TAOLI
	IF((ITASK.EQ.1).AND.(PV%CRUN.GT.PV%STARTRUN).AND.(.NOT.PV%KILLSOIL)) THEN
		 NL=pv%pnl;			pv%pevap = 0.0;		pv%ptrans = 0.0 
		 pv%petp = 0.0;		pv%peta = 0.0;		
		 pv%plai = 0.0;     pv%pdae = DAE
		 CALL SOIL_RESERVE("GET_VALUES")
		 WL0= pv%pwl0
		 do i = 1, nl
			 tkl(i)=pv%pdlayer(i) /1000.0;		wclqt(i)= pv%pswc(i)	
			 wcst(i)=pv%pwcst(i);				wcfc(i)= pv%pwcfc(i)
			 wcwp(i)=pv%pwcwp(i);				wcad(i) = pv%pwcad(i)			 	
			 IF(PV%PBD(I).LE.0.0) THEN
				 pv%pbd(i)= 2.65*(1-PV%PWCST(I))
			 ENDIF
		 enddo
		 DEALLOCATE(TSI)
	ELSE
	     pv%pnl = nl;			pv%pevap = evsc;	pv%ptrans = trw 
	     pv%petp = trc + evsc;	pv%peta = trc;     pv%pdae = DAE
	     pv%pwl0 = wl0;			pv%plai = lai
	     do i = 1, nl
		     pv%pdlayer(i) = tkl(i)*1000.0;		pv%pswc(i) = wclqt(i)	
		     pv%pwcst(i) = wcst(i);				pv%pwcfc(i) = wcfc(i)
		     pv%pwcwp(i) = wcwp(i);				pv%pwcad(i) = wcad(i);	
		     IF(PV%PBD(I).LE.0.0) THEN
			     pv%pbd(i)= 2.65*(1-PV%PWCST(I))
		     ENDIF
	     enddo
	END IF
!+	added by TAOLI, 13 Sept 2009
!+-------------------------------------------------------------------------------------

!==============================================================!
! Output writing only at ITASK = 2                             !
!==============================================================!
      IF (ITASK.EQ.2) THEN
         IF (OUTPUT) THEN
            CALL OUTDAT (2, 0, 'YEAR', YEAR)
            CALL OUTDAT (2, 0, 'DOY' , DOY)
            CALL OUTDAT (2, 0, 'CROPSTA ', REAL (CROPSTA))
            CALL OUTDAT (2, 0, 'ETD' , ETD)
            CALL OUTDAT (2, 0, 'TRC' , TRC)
            CALL OUTDAT (2, 0, 'EVSC', EVSC)
            CALL OUTDAT (2, 0, 'RAIN', RAIN)
            !CALL OUTDAT (2, 0, 'RAINCU',   RAINCU)
         END IF
         IF((INDEX(ESTAB,'DIRECT-SEED').GT.0).AND.(GERMIN.EQ.1)) THEN   !ADDITIONAL SECTION TO ESTIMATE THE GERMINATE RATE
             IF((PV%PSOWYEAR.EQ.YEAR).AND.(PV%PSOWDAY.EQ.DOY)) THEN
                 GERMCAL=.TRUE.
             END IF
             IF(GERMCAL.AND.(CROPSTA.LT.1)) THEN
                 IF(((RAIN.GE.(EVSC+1.0)).OR.(MSKPA(1).LE.5.0)).AND.(ISGERMDELY)) THEN
                     ISGERMDELY = .FALSE.
                 END IF
                 IF(.NOT. ISGERMDELY) THEN
                     IF(WL0.GE.50.0)THEN
                         PV%GERMRATE = PV%GERMRATE + 0.00; PV%GERMDAYS = PV%GERMDAYS + 1.0; PV%TOOWET = PV%TOOWET + 1; PV%TOODRY = 0
                     ELSEIF((WL0.GT.0.0).AND.(WL0.LT.50.0)) THEN
                         PV%TOOWET = 0; PV%GERMDAYS = PV%GERMDAYS + 1.0; PV%TOODRY = 0
                         PV%GERMRATE = PV%GERMRATE !+ 0.40*((30.0-WL0)/30.0)**2
                     ELSE
                         PV%TOOWET = 0; PV%GERMDAYS = PV%GERMDAYS + 1.0
                         IF(MSKPA(1).LE.10.0) THEN
                              PV%TOODRY = 0; PV%GERMRATE = PV%GERMRATE + 0.70 + 0.03*MSKPA(1)
                         ELSEIF(MSKPA(1).GT.70.0) THEN
                             PV%GERMRATE = PV%GERMRATE; PV%TOODRY = PV%TOODRY + 1
                         !ELSEIF((MSKPA(1).GT.10.0).AND.(MSKPA(1).LE.30.0)) THEN
                         !    PV%GERMRATE = PV%GERMRATE+1.0; PV%TOODRY = 0
                         ELSEIF((MSKPA(1).GT.10.0).AND.(MSKPA(1).LE.70.0)) THEN
                             PV%TOODRY = 0; PV%GERMRATE = PV%GERMRATE + 1.0-1.0/60.0*(MSKPA(1)-10.0) 
                         END IF                     
                     END IF
                 ELSE                     
                     IF(MOD(REAL(EMYR),4.0).EQ.0.0) THEN
                         IF((EMD+1.0).GT.366) THEN
                             EMD = (EMD + 1) - 366; EMYR = EMYR + 1
                         ELSE
                             EMD = EMD + 1
                         END IF
                     ELSE
                         IF((EMD+1.0).GT.365) THEN
                             EMD = (EMD + 1) - 365; EMYR = EMYR + 1
                         ELSE
                             EMD = EMD + 1
                         END IF
                     END IF
                     GERMDELY = GERMDELY + 1; PV%GERMRATE=PV%GERMRATE - 0.02
                 End If
             END IF
         ELSE
             PV%GERMRATE=1.0;PV%GERMDAYS=1.0;PV%TOOWET=0;PV%TOODRY=0
         END IF
         CALL OUTDAT (2, 0, 'GERMIN', PV%GERMRATE)
         IF(GERMDELY.GE.10) THEN
             TERMNL = .TRUE.
         END IF
      END IF

!-----Set CROPSTA: 0=before sowing; 1=day of sowing; 2=in seedbed;
!                  3=day of transplanting; 4=main growth period
      IF (ITASK.EQ.1 .OR. ITASK.EQ.3) THEN
         IF (CROPSTA .EQ. 3) CROPSTA = 4
         IF (CROPSTA .EQ. 2) THEN
            IF (DAE .EQ. REAL(SBDUR)) CROPSTA = 3
         END IF
         IF (CROPSTA .EQ. 1) THEN
            IF (ESTAB.EQ.'TRANSPLANT') THEN
               CROPSTA = 2
            ELSE IF (ESTAB.EQ.'DIRECT-SEED') THEN
               CROPSTA = 4
            END IF
        END IF
        IF (CROPSTA .EQ. 0) THEN
           !IDATE = DTFSECMP(EMYR, EMD, IYEAR, IDOY)
           !IF (IDATE .EQ. 0) THEN
            IF((EMYR.EQ.IYEAR).AND.(EMD.EQ.IDOY)) THEN
              CROPSTA = 1
           !ELSE IF (IDATE .EQ. 1) THEN
            ELSEIF((IYEAR.GT.EMYR).OR.((EMYR.EQ.IYEAR).AND.(EMD.LT.IDOY))) THEN
              CALL FATALERR ('MODELS',  &
                'Time past supplied sowing date or year')
           END IF
        END IF

      END IF
!============================================================*
!-----Integration section
!============================================================*

      IF (ITASK .EQ. 3) THEN
!-------Summation of some state variables
        RAINCU = RAINCU + RAIN

      END IF

!==============================================================!
! Terminal calculations at ITASK = 4                           !
!==============================================================!
      IF (ITASK.EQ.4) THEN
!         WRITE (IUNITO,'(A)') '*'
!         WRITE (IUNITO,'(A)') '* FSE driver info:'
         WRITE (IUNITO,'(A,T7,A,I5,A,I4,A)') &
           '*','Year:',IYEAR,', day:',IDOY,', System end'
!---     Terminal output
         IF(LEN_TRIM(PV%OPSTRING).GT.1) THEN 
            IF(INDEX(PV%OPSTRING,'EMD').GT.0) CALL OPSTOR ('EMD', 1.0*EMD)
            IF(INDEX(PV%OPSTRING,'DAE').GT.0) CALL OPSTOR ('DAE', DAE)
         ELSE
           CALL OPSTOR ('EMD', 1.0*EMD)
           CALL OPSTOR ('DAE', DAE)
         END IF
      END IF
!==============================================================!
! End of section ITASK = 4                                     !
!==============================================================!

      RETURN
      END





