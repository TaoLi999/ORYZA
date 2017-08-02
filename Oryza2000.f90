!----------------------------------------------------------------------!
! Program ORYZA2000: Simulation model for potential and water- and     !
!                    nitrogen-production of rice.                      !
! Version 4.0                                                          !
! Date      : December 2001                                            !
! Programmed: B.A.M. Bouman                                            !
! History   : Adapted from ORYZA1 (1995), and ORYZA_W (1995) models    !
!----------------------------------------------------------------------!
      PROGRAM ORYZA2000
!      USE DFLIB
!       USE IFNLS   !can read multiple language input
!      INCLUDE 'link_fnl_shared_hpc.h'
      INCLUDE 'for_iosdef.for'
	  character(128) FILEIC
	  CHARACTER(128) THEFILE
	  CHARACTER*20 INITINP
	  INTEGER I, IUNITC
      LOGICAL I_OPEN
            
      DO I=1,600
          INQUIRE(I,OPENED=I_OPEN)
          IF(I_OPEN) CLOSE(I)
      END DO

      OPEN(UNIT=1, FILE='control.dat', STATUS='OLD',IOSTAT=IERR, ERR=100)
      100  IF (IERR .EQ. FOR$IOS_FILNOTFOU) THEN
        OPEN(UNIT=1, FILE=' ', STATUS = 'OLD')
      END IF

      !---------------------------------------------
      !Version information
      WRITE(*,*) '||=============================================='
      WRITE(*,*) '|| ORYZA version 3.                             '
      WRITE(*,*) '|| Modified by Dr. Tao Li, 14 Aug 2013.         '
      WRITE(*,*) '|| Owner: International Rice Research Institute.'
      WRITE(*,*) '|| Any difficulty for running, contact Dr. Tao Li'
      WRITE(*,*) '|| at t.li@irri.org                             '
      WRITE(*,*) '||=============================================='
      FILEC=''
      DO  
	  READ(1,'(A)',END=20) THEFILE
        THEFILE =TRIM(ADJUSTL(THEFILE))
      IF(LEN_TRIM(THEFILE).GT.0) THEN  
       I = INDEX(THEFILE, "=")  
        IF(I.GT.11) THEN
            INITINP = trim(THEFILE(1:I))
            CALL UPPERC(INITINP)            
            IF(index(trim(INITINP),'*CONTROLFILE').gt.0) THEN
                FILEIC = TRIM(THEFILE((I+1):LEN_TRIM(THEFILE)))
                CLOSE(1)
                GOTO 20
            ELSEIF(INDEX(INITINP,'CONTROLFILE').GT.0) THEN
                 FILEIC= trim(THEFILE((i+1):LEN_TRIM(THEFILE)))
                 FILEIC = TRIM(FILEIC(2:LEN_TRIM(FILEIC)))
                 I=INDEX(FILEIC,"!")
                 IF(I.GT.0)  FILEIC = TRIM(FILEIC(1:(I-1)))
                 I=INDEX(FILEIC,"*")
                 IF(I.GT.0)  FILEIC = TRIM(FILEIC(1:(I-1)))
                 I=LEN_TRIM(FILEIC) !INDEX(FILEIC,"'")
                 FILEIC = FILEIC(2:(I-1))
                 CLOSE(1)
                 GOTO 20
            END IF
        END IF
      END IF
    END DO
20      IF(LEN_TRIM(FILEIC).GT.0) THEN
        CALL FSE(FILEIC)
      ELSE      !Added by TAOLI, 17DEC 2012
40       WRITE(*,'(A)') 'There is no control file or path provided!'
        stop
      END IF
    END

!----------------------------------------------------------------------*
!                                                                      *
!                                                                      *
!              FORTRAN Simulation Environment (FSE 2.1)                *
!                            March, 1995                               *
!                                                                      *
! THIS VERSION WITH MINR ADAPTATIONS FOR MODEL ORYZA 4.0 (BOUMAN)      *
!                                                                      *
!     FSE 2.1 is a simulation environment suited for simulation of     *
!     biological processes in time, such as crop and vegetation growth,*
!     insect population development etc.                               *
!                                                                      *
!     The MAIN program, subroutine FSE and subroutine MODELS are       *
!     programmed by D.W.G. van Kraalingen, DLO Institute for           *
!     Agrobiological and Soil Fertility Research (AB-DLO),             *
!     PO Box 14, 6700 AA, Wageningen, The Netherlands (e-mail:         *
!     d.w.g.van.kraalingen@ab.agro.nl).                                *
!                                                                      *
!     FSE version 2.1 is described in:                                 *
!     Kraalingen, D.W.G. van 1995. The FSE system for crop simulation, *
!     Simulation Report no. 23, Institute for Agrobiological and Soil  *
!     Fertility Research and Dept. of Theoretical Production Ecology,  *
!     Agricultural University Wageningen.                              *
!                                                                      *
!     Data files needed for FSE 2.1:                                   *
!          (excluding data files used by models called from MODELS):   *
!        - CONTROL.DAT (contains file names to be used),               *
!        - timer file whose name is specified in CONTROL.DAT,          *
!        - optionally, a rerun file whose name is specified in         *
!          CONTROL.DAT,                                                *
!        - weather data files as specified in timer file               *
!     Object libraries needed for FSE 2.1:                             *
!        - TTUTIL (at least version 3.2)                               *
!        - WEATHER (at least version from 17-Jan-1990)                 *
!     Differences with standard version:                               *
!        - common block FSECM1 activated                               *
!        - calls to OP* subroutines activated                          *
!        - RDFROM logical set to false (unused rerun variables not     *
!          fatal)                                                      *
!        - calls to OBS* routine activated                             *
!----------------------------------------------------------------------*

      SUBROUTINE FSE (FILEIC)

      !USE CHART
	  USE public_module	!VARIABLES
	  use Module_OutDat
             
      IMPLICIT NONE

!-----Standard declarations for simulation and output control

      CHARACTER*(*) FILEIC
      INTEGER       ITASK   , INSETS, ISET  , IPFORM, IL 
      LOGICAL       OUTPUT  , TERMNL, RDINQR, STRUNF, ENDRNF
      CHARACTER (1) COPINF, DELTMP
	  CHARACTER(3) theLEN
      INTEGER       INPRS   , STRUN , ENDRUN
      LOGICAL ISEXIS, SAMESOIL, KILLSOIL1
      INTEGER        IMNPRS
      PARAMETER     (IMNPRS=200)
      CHARACTER (11) PRSEL(IMNPRS)
	  
!-----Declarations for time control
      INTEGER   IDOY, IYEAR
      REAL      DELT, DOY, FINTIM, PRDEL, STTIME, TIME, YEAR

!-----Declarations for weather system
      INTEGER   IFLAG    , ISTAT1, ISTAT2 , ISTN, UNITWL 
      !THE UNITWL was added to close weather log file correctly, TAOLI, 29 JAN 2013
      REAL      ANGA     , ANGB  , ELEV   , LAT , LONG, RDD
      REAL      TMMN     , TMMX  , VP     , WN  , RAIN
      LOGICAL   WTRMES   , WTRTER
      CHARACTER (80) WTRDIR,ESTAB
      CHARACTER (7)  CNTR
      CHARACTER (6)  WSTAT
      CHARACTER (1)  DUMMY
	  CHARACTER (4)	 MYOUT1, MYOUT2, SOILKILL, TOPHENOL
      CHARACTER(256) MYOP
	  CHARACTER ISTNC*9,MULTIY*3,ISTNCSTR*14
!-----Declarations for file names and units
      INTEGER   IUNITR, IUNITD, IUNITO, IUNITL, IUNITC, TOTALGROWTHDAYS, TOTALDAE,EMD,EMYR,DAYSDIFF
      CHARACTER (128) FILEON, FILEOL
      CHARACTER (128)  FILEIR, FILEIT
      CHARACTER (128) FILEI1, FILEI2, FILEI3, FILEI4, FILEI5
      CHARACTER (128) FILEIX
!-----Declarations for observation data facility
      INTEGER   INOD , IOD, I
      LOGICAL I_OPEN,I_EXIST

      INTEGER   IMNOD
      PARAMETER (IMNOD=100)
      INTEGER   IOBSD(IMNOD)

!-----For communication with OBSSYS routine
      COMMON /FSECM1/ YEAR,DOY,IUNITD,IUNITL,TERMNL

!     File name for control file and empty strings for input files 1-5.
!     WTRMES flags any messages from the weather system

!      DATA FILEIC /'CONTROL.DAT'/
      DATA FILEI1 /' '/, FILEI2 /' '/, FILEI3 /' '/
      DATA FILEI4 /' '/, FILEI5 /' '/
      DATA WTRMES /.FALSE./
!     FILEI1 for crop parameter; FILEI2 for soil parameter; FILEI3 for soil salinity dynamics; FILEI4 for submergence

      DATA STRUNF /.FALSE./, ENDRNF /.FALSE./

!-----Unit numbers for control file (C), data files (D),
!     output file (O), log file (L) and rerun file (R).
      IUNITC = 10
      IUNITD = 20
      IUNITO = 30
      IUNITL = 40
      IUNITR = 50
!-----Multilple year weather control
      ISTNCSTR = ""
!-----get control file, added by TAOLI, 1 Oct, 2010
	  ALLOCATE(pv)
	  PV%KILLSOIL = .TRUE.;PV%OPSTRING = ''
      CALL UPPERC(FILEIC)
!-----Open control file and read names of normal output file, log file
!     and rerun file (these files cannot be used in reruns)
      CALL RDINIT (IUNITC,0, FILEIC)
      CALL RDSCHA ('FILEON', FILEON)
      CALL RDSCHA ('FILEOL', FILEOL)
!     Check whether the log and res files are still opened
      FILEIX = ADJUSTL(FILEON)
      I=INDEX(FILEIX,".")
      FILEIX = FILEIX(1:I)
      FILEIX = TRIM(FILEIX)//".BIN"
      INQUIRE (FILE=FILEON, OPENED=I_OPEN,EXIST = I_EXIST,NUMBER=I)
      IF(I_OPEN) THEN
          CLOSE(I,STATUS ='DELETE')
      ELSE IF(I_EXIST) THEN
          OPEN(IUNITD,FILE =FILEON)
          CLOSE(IUNITD,STATUS ='DELETE')
      END IF
      INQUIRE (FILE=FILEIX, OPENED=I_OPEN,EXIST =I_EXIST,NUMBER=I)
      IF(I_OPEN) THEN
          CLOSE(I,STATUS ='DELETE')
      ELSE IF (I_EXIST) THEN
          OPEN(IUNITD,FILE=FILEIX)
          CLOSE(IUNITD,STATUS ='DELETE')
      END IF
      INQUIRE (FILE=FILEOL, OPENED=I_OPEN,EXIST=I_EXIST, NUMBER=I)
      IF(I_OPEN) THEN
          CLOSE(I,STATUS ='DELETE')
      ELSE IF(I_EXIST) THEN
          OPEN(IUNITD,FILE=FILEOL)
          CLOSE(IUNITD,STATUS ='DELETE')
      END IF
      INQUIRE (FILE='WEATHER.LOG', OPENED=I_OPEN,EXIST=I_EXIST, NUMBER=I)
      IF(I_OPEN) THEN
          CLOSE(I,STATUS ='DELETE')
      ELSE IF(I_EXIST) THEN
          OPEN(IUNITD,FILE='WEATHER.LOG')
          CLOSE(IUNITD,STATUS ='DELETE')
      END IF      
      
      CALL RDSCHA ('FILEIR', FILEIR)
!-----DETERMINE WHETHER THE SOIL WILL BE CONTINUOUESLY SIMULATION FOR FOLLOWING CROP SEASON
	  IF(RDINQR('SOILKILL')) THEN
			CALL RDSCHA ('SOILKILL', SOILKILL)
			CALL UPPERC(SOILKILL)
			IF(INDEX(SOILKILL,"YES").GT.0) THEN
				KILLSOIL1 = .TRUE.
			ELSE
				KILLSOIL1 = .FALSE.
			END IF					
	  ELSE
			KILLSOIL1 = .TRUE.	
	  ENDIF

!     check if start run number was found, if there, read it
      IF (RDINQR('STRUN'))  THEN
         CALL RDSINT ('STRUN',STRUN)
         STRUNF = .TRUE.
      END IF
!     check if end run number was found, if there, read it
      IF (RDINQR('ENDRUN')) THEN
         CALL RDSINT ('ENDRUN',ENDRUN)
         ENDRNF = .TRUE.
      END IF
      CLOSE (IUNITC)

!-----Open output file and possibly a log file
      CALL FOPENS  (IUNITO, FILEON, 'NEW', 'DEL')
      IF (FILEOL.NE.FILEON) THEN
         CALL FOPENS  (IUNITL, FILEOL, 'NEW', 'DEL')
      ELSE
         IUNITL = IUNITO
      END IF

!     initialization of logfile for processing of end_of_run values
!      CALL OPINIT

!-----See if rerun file is present, and if so read the number of rerun
!     sets from rerun file

      CALL RDSETS (IUNITR, IUNITL, FILEIR, INSETS)

!-----Initialise logfile for end-of-year state values
!&      CALL OPINIT
           CALL OPINIT(INSETS+1)

!======================================================================*
!======================================================================*
!                                                                      *
!                   Main loop and reruns begin here                    *
!                                                                      *
!======================================================================*
!======================================================================*

      IF (.NOT.ENDRNF) THEN
!        no end run was found in control.dat file
         ENDRUN = INSETS
      ELSE
         ENDRUN = MAX (ENDRUN, 0)
         ENDRUN = MIN (ENDRUN, INSETS)
      END IF

      IF (.NOT.STRUNF) THEN
!        no start run was found in control.dat file
         STRUN = 0
      ELSE
         STRUN = MAX (STRUN, 0)
         STRUN = MIN (STRUN, ENDRUN)
      END IF
!@@@  prepareing a  file with unit 600 for daily layer output (temporally for test ORYZA3 and DSSAT-ORYZA)
!@@@      open(UNIT=600, FILE='DailyLayer.txt')
 !     CALL ChartInit (IUNITO+2)

      DO ISET=STRUN,ENDRUN
      WRITE (*,'(A)') '   ORYZA3: Initialize model'
      
      !@@@ Temporal section for testing DSSAT-ORYZA --- Start of section
      
      !@@@ WRITE(FILEI5,'(I5)') ISET
      !@@@WRITE(600, '(A)') 'RERUN '//TRIM(FILEI5) 
       !@@@ Temporal section for testing DSSAT-ORYZA --- end of section     
 !-----Select data set
      CALL RDFROM (ISET, .FALSE.)
      PV%TOTALGROWTHDAYS = 0

      !CALL ChartSetRunID (ISET)


!======================================================================*
!                                                                      *
!                        Initialization section                        *
!                                                                      *
!======================================================================*

      ITASK  = 1
      TERMNL = .FALSE.
      WTRTER = .FALSE.

!-----Read names of timer file and input files 1-5 from control
!     file (these files can be used in reruns)
      CALL RDINIT (IUNITC,IUNITL,FILEIC)
      CALL RDSCHA ('FILEIT', FILEIT)
!      CALL TABTOSPACE(FILEIT)        !Check if TAB exist, then replace it with space
      IF (RDINQR ('FILEI1')) THEN
          CALL RDSCHA ('FILEI1', FILEI1)
!          CALL TABTOSPACE(FILEI1)        !Check if TAB exist, then replace it with space
      END IF
      IF (RDINQR ('FILEI2')) THEN
          CALL RDSCHA ('FILEI2', FILEI2)
!          CALL TABTOSPACE(FILEI2)        !Check if TAB exist, then replace it with space
      END IF
      IF (RDINQR ('FILEI3')) THEN
          CALL RDSCHA ('FILEI3', FILEI3)
!          CALL TABTOSPACE(FILEI3)        !Check if TAB exist, then replace it with space
      END IF
      IF (RDINQR ('FILEI4')) THEN
        CALL RDSCHA ('FILEI4', FILEI4)
!        CALL TABTOSPACE(FILEI4)        !Check if TAB exist, then replace it with space
      END IF
      IF (RDINQR ('FILEI5')) THEN
          CALL RDSCHA ('FILEI5', FILEI5)
!          CALL TABTOSPACE(FILEI5)        !Check if TAB exist, then replace it with space
      END IF
      CALL RDSREA ('PRDEL' , PRDEL )
      CALL RDSINT ('IPFORM', IPFORM)
      CALL RDSCHA ('COPINF', COPINF)
      CALL RDSCHA ('DELTMP', DELTMP )
      CALL RDSINT ('IFLAG' , IFLAG)
!-----If rotation in the same soil, then force the IYEAR and STTIME to be the end time of last rerun
	  	IF((.NOT.KILLSOIL1).AND.(ISET.GT.STRUN)) THEN
		  	pv%CRUN=ISET		!ALL SOIL INFORMATION WILL NE KEPT	
		  	IYEAR = PV%PYEAR; STTIME = REAL(PV%PDOY)
		  	PV%KILLSOIL = .FALSE.; PV%OPSTRING=MYOP	  
	  	ELSEIF((KILLSOIL1).AND.(ISET.GT.STRUN)) THEN
	        ALLOCATE(pv)
		  	CALL SET_ZERO_PV
		  	pv%CRUN=ISET; PV%KILLSOIL = .TRUE.
		  	PV%OPSTRING='';PV%ISPHENAD=.FALSE.
	  	ELSE
		  	CALL SET_ZERO_PV
		  	pv%CRUN=ISET;PV%STARTRUN = STRUN; PV%KILLSOIL = KILLSOIL1; PV%OPSTRING=MYOP
	  	END IF
!-----determine whether the soil information is correct for rotation in the same soil and soil process continue
	    IF((ISET.EQ.STRUN).AND.(.NOT.KILLSOIL1)) THEN
			CALL CHECK_RERUN(FILEIR,(endrun-strun+1), FILEI2,SAMESOIL)
			IF(.NOT.SAMESOIL) THEN
				PV%KILLSOIL = .TRUE.;KILLSOIL1 =.TRUE.
				CALL FATALERR ('Rerun data file','Try rotation in same soil, but different soils are given!')
			ELSE
			    PV%KILLSOIL = .FALSE.;KILLSOIL1 =.FALSE.
			END IF
		END IF
!     GET THE OUTPUT VARIABLE STRING FOR OP.DAT IF IT IS EXISTED
!     THE OPSTRING SHOULD BE 'XXX,AAA,BBB,CCCC'
          IF(RDINQR('OPSTRING')) THEN
            CALL RDSCHA('OPSTRING',PV%OPSTRING) 
          ELSE
            PV%OPSTRING = ''
          END IF
          IF(RDINQR('ISPHENAD')) THEN
            CALL RDSCHA('ISPHENAD',TOPHENOL) 
            IF(INDEX(TOPHENOL,'YES').GT.0) THEN
                PV%ISPHENAD=.TRUE.
            ELSE
                PV%ISPHENAD=.FALSE.
            END IF       
          ELSE
            PV%ISPHENAD=.FALSE.
          END IF
		
!-----See if observation data variable exists, if so read it
      INOD = 0
      IF (RDINQR('IOBSD')) THEN
         CALL RDAINT ('IOBSD' , IOBSD, IMNOD, INOD)
         IF (IOBSD(1).EQ.0) INOD = 0
      END IF

!-----See if variable with print selection exists, if so read it
      INPRS = 0
      IF (RDINQR('PRSEL')) CALL RDACHA ('PRSEL',PRSEL,IMNPRS,INPRS)

      CLOSE (IUNITC)

!-----Read time, control and weather variables from timer file
      CALL RDINIT (IUNITD  , IUNITL, FILEIT)
      CALL RDSCHA ('ESTAB'  , ESTAB)
      CALL RDSINT ('EMD' , EMD )
      CALL RDSINT('EMYR',EMYR  )
      CALL RDSREA ('FINTIM', FINTIM)
      IF(RDINQR('DELT')) THEN
      	CALL RDSREA ('DELT'  , DELT  )
      ELSE
        DELT=1.0
      ENDIF
      CALL RDSREA ('STTIME', STTIME)
      CALL RDSINT ('IYEAR' , IYEAR )
      IF(INDEX(ESTAB,"DIRECT-SEED").GT.0) THEN
          IF(DAYSDIFF(EMYR,EMD,IYEAR,INT(STTIME)).LT.9) THEN
              IF((EMD-9).LE.0) THEN
                  PV%PSOWYEAR=EMYR-1;I=PV%PSOWYEAR
                  IF(MOD(REAL(I),4.0).EQ.0.0) THEN
                      PV%PSOWDAY = EMD+366-9
                  ELSE
                      PV%PSOWDAY = EMD+365-9
                  END IF
              ELSE
                  PV%PSOWYEAR=EMYR;I=PV%PSOWYEAR
                  IF(MOD(REAL(I),4.0).EQ.0.0) THEN
                      PV%PSOWDAY = EMD-9
                  ELSE
                      PV%PSOWDAY = EMD-9
                  END IF
              END IF
              IYEAR=PV%PSOWYEAR;STTIME=REAL(PV%PSOWDAY)
          ELSE
              IF((EMD-9).LE.0) THEN
                  PV%PSOWYEAR=EMYR-1;I=PV%PSOWYEAR
                  IF(MOD(REAL(I),4.0).EQ.0.0) THEN
                      PV%PSOWDAY = EMD+366-9
                  ELSE
                      PV%PSOWDAY = EMD+365-9
                  END IF
              ELSE
                  PV%PSOWYEAR=EMYR;I=PV%PSOWYEAR
                  IF(MOD(REAL(I),4.0).EQ.0) THEN
                      PV%PSOWDAY = EMD-9
                  ELSE
                      PV%PSOWDAY = EMD-9
                  END IF
              END IF
          END IF
      END IF
      CALL RDSINT ('ISTN'  , ISTN  )
      CALL RDSCHA ('WTRDIR', WTRDIR)
      CALL RDSCHA ('CNTR'  , CNTR)
      IF(RDINQR('MULTIY')) THEN
            CALL RDSCHA ('MULTIY', MULTIY)            
      END IF
      CLOSE (IUNITD)

      !Modified by dr. TAOLI 18Aug 2014, for large scale modeling management,
      !the weather station information can be provided in control file,
      !the information of weather station in experiment file will be overwriten
      !bu those in control file
       CALL RDINIT (IUNITC,IUNITL,FILEIC)
        IF (RDINQR('ISTN')) THEN
            CALL RDSINT ('ISTN'  , ISTN  )
        END IF
        IF (RDINQR('WTRDIR')) THEN
            CALL RDSCHA ('WTRDIR', WTRDIR)
        END IF
        IF (RDINQR('CNTR')) THEN
            CALL RDSCHA ('CNTR'  , CNTR)
        END IF
		IF(RDINQR('MULTIY')) THEN
            CALL RDSCHA ('MULTIY', MULTIY)            
		END IF
      CLOSE (IUNITC)  !End the update, TAOLI 18Aug 2014
      I = LEN_TRIM(WTRDIR)
      IF(I.GT.2) THEN
          IF(INDEX(WTRDIR,"/").GT.0) THEN
              IF(WTRDIR(I:I).NE."/") WTRDIR = ADJUSTL(TRIM(WTRDIR))//"/"
          ELSE
              IF(WTRDIR(I:I).NE."\") WTRDIR = ADJUSTL(TRIM(WTRDIR))//"\"
          END IF              
      END IF
	  IF(INDEX(MULTIY,'YES').GT.0) THEN
		 WRITE(ISTNC,'(I9)') ISTN
         IF(TRIM(ADJUSTL(CNTR))//TRIM(ADJUSTL(ISTNC)).NE.TRIM(ADJUSTL(ISTNCSTR))) THEN
             FILEI5 =TRIM(ADJUSTL(WTRDIR))//TRIM(ADJUSTL(CNTR))//TRIM(ADJUSTL(ISTNC))//".cli"
		    CALL WEATHER_SPLIT(FILEI5)
            ISTNCSTR = TRIM(ADJUSTL(CNTR))//TRIM(ADJUSTL(ISTNC))
         END IF
      END IF
!-----Read the information of maximum growth duration from crop file
      CALL RDINIT (IUNITC,IUNITL,FILEI1)
      IF(RDINQR('TOTALDAE')) THEN
          CALL RDSREA('TOTALDAE',TOTALDAE)
          TOTALDAE=TOTALDAE*1.15   !The critical total growth duration is 15% more than given maximum growth duration
      ELSE
          TOTALDAE=250.0
      END IF
      CLOSE(IUNITC)
!---- If STTIME is not a whole day value, transfrom:
      STTIME = 1.0 *INT(STTIME)

!-----Initialize TIMER and OUTDAT routines
      CALL TIMER2 (ITASK, STTIME, DELT, PRDEL, FINTIM, &
                   IYEAR, TIME  , DOY , IDOY , TERMNL, OUTPUT)
      YEAR = REAL (IYEAR)
      CALL OUTDAT (ITASK, IUNITO, 'TIME', TIME)
     
!      CALL ChartInitialGroup

!-----Open weather file and read station information and return
!     weather data for start day of simulation.
!     Check status of weather system, WTRMES flags if warnings or errors
!     have occurred during the whole simulation. WTRTER flags if the run
!     should be terminated because of missing weather

      CALL STINFO (IFLAG , WTRDIR, ' ', CNTR, ISTN, IYEAR, &
                   ISTAT1, LONG  , LAT, ELEV, ANGA, ANGB, UNITWL)
      CALL WEATHR (IDOY  , ISTAT2, RDD, TMMN, TMMX, VP, WN, RAIN)
      IF (ISTAT1.NE.0.OR.ISTAT2.NE.0) WTRMES = .TRUE.
      WSTAT  = '444444'
      IF (ABS (ISTAT2).GE.111111) THEN
         WRITE (WSTAT,'(I6)') ABS (ISTAT2)
      ELSE IF (ISTAT2.EQ.0) THEN
         WSTAT = '111111'
      END IF
!     Test the adjustment of phenology paramters BASE, OPTIMAL AND MAXIMUM due to 
!     the change of elevation and latitude. TAOLI 9FEB 2014
      IF(PV%ISPHENAD) THEN
        CALL PHENAD(LAT,ELEV, DOY) !TAOLI 9FEB 2014
      END IF
!-----initialize OBSSYS routine
      IF (ITASK.EQ.1) CALL OBSINI

!-----Conversion of total daily radiation from kJ/m2/d to J/m2/d
      RDD = max(3000.,RDD*1000.)
!-----Call routine that handles the different models
      CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL, &
                   FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5, &
                   OUTPUT, TERMNL,  &
                   DOY   , IDOY  , YEAR  , IYEAR , STTIME,  &
                   TIME  , DELT  , &
                   LAT   , WSTAT , WTRTER, &
                   RDD   , TMMN  , TMMX  , VP    , WN, RAIN)


!======================================================================*
!                                                                      *
!                      Dynamic simulation section                      *
!                                                                      *
!======================================================================*

      WRITE (*,'(A)') '   ORYZA3: DYNAMIC loop'

      DO WHILE (.NOT.TERMNL)

!----------------------------------------------------------------------*
!                     Integration of rates section                     *
!----------------------------------------------------------------------*

      IF (ITASK.EQ.2) THEN

!--------Carry out integration only when previous task was rate
!        calculation

         ITASK = 3

!--------Call routine that handles the different models
         CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL, &
                      FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5, &
                      OUTPUT, TERMNL,  &
                      DOY   , IDOY  , YEAR  , IYEAR , STTIME,  &
                      TIME  , DELT  , &
                      LAT   , WSTAT , WTRTER, &
                      RDD   , TMMN  , TMMX  , VP    , WN, RAIN)

!--------Turn on output when TERMNL logical is set to .TRUE.
         
         IF(PV%TOTALGROWTHDAYS.GT.TOTALDAE) THEN
             WRITE(*,*) "Total growth days more than the critical total growth duration. Program will be terminated!"
             TERMNL = .TRUE.
         END IF
         IF (TERMNL.AND.PRDEL.GT.0.) OUTPUT = .TRUE.
      END IF

!----------------------------------------------------------------------*
!               Calculation of driving variables section               *
!----------------------------------------------------------------------*

      ITASK = 2
!+---------------------------------------------------------------------*
!	update model global variable values, TAOLI, 13 Sept. 2009
!+		setting doy and year into public module for other models
	
		pv%pdoy = int(doy); pv%pyear = int(year)
		pv%PRad = rdd/1000.0		!PRad is in Mj/d
		pv%pmaxt = tmmx;		pv%pminT= tmmn;		pv%prain = rain
		pv%pwind = wn;		pv%pPressur = vp
!+		add by TAOLI, 13 Sept, 509	

!-----Write time of output to screen and file
        CALL OUTDAT (2, 0, 'TIME', TIME)
            !CALL ChartNewGroup

      IF (OUTPUT) THEN
         IF (ISET.EQ.0) THEN
 !           WRITE (*,'(13X,A,I5,A,F7.2)') &
 !             'Default set, Year:', IYEAR, ', Day:', DOY
         ELSE
 !           WRITE (*,'(13X,A,I3,A,I5,A,F7.2)') &
 !             'Rerun set:', ISET, ', Year:', IYEAR, ', Day:', DOY
         END IF
      END IF

!-----Get weather data for new day and flag messages
      CALL STINFO (IFLAG , WTRDIR, ' ', CNTR, ISTN, IYEAR, &
                   ISTAT1, LONG  , LAT, ELEV, ANGA, ANGB, UNITWL)
      CALL WEATHR (IDOY, ISTAT2, RDD, TMMN, TMMX, VP, WN, RAIN)
      IF (ISTAT1.NE.0.OR.ISTAT2.NE.0) WTRMES = .TRUE.
      WSTAT  = '444444'
      IF (ABS (ISTAT2).GE.111111) THEN
         WRITE (WSTAT,'(I6)') ABS (ISTAT2)
      ELSE IF (ISTAT2.EQ.0) THEN
         WSTAT = '111111'
      END IF

!-----Conversion of total daily radiation from kJ/m2/d to J/m2/d
      RDD = max(3000.,RDD*1000.)
!----------------------------------------------------------------------*
!               Calculation of rates and output section                *
!----------------------------------------------------------------------*

!-----Call routine that handles the different models
      CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL, &
                   FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5, &
                   OUTPUT, TERMNL,  &
                   DOY   , IDOY  , YEAR  , IYEAR , STTIME,  &
                   TIME  , DELT  , &
                   LAT   , WSTAT , WTRTER, &
                   RDD   , TMMN  , TMMX  , VP    , WN, RAIN)

      IF (TERMNL.AND..NOT.OUTPUT.AND.PRDEL.GT.0.) THEN
!--------Call model routine again if TERMNL is switched on while
!        OUTPUT was off (this call is necessary to get output to file
!        when a finish condition was reached and output generation
!        was off)
         IF (ISET.EQ.0) THEN
            WRITE (*,'(13X,A,I5,A,F7.2)') &
              'Default set, Year:', IYEAR, ', Day:', DOY
         ELSE
            WRITE (*,'(13X,A,I3,A,I7,A,F7.2)') &
              'Rerun set:', ISET, ', Year:', IYEAR, ', Day:', DOY 
         END IF
         OUTPUT = .TRUE.
         CALL OUTDAT (2, 0, 'TIME', TIME)
             !CALL ChartNewGroup
    !        patch dvk
         CALL OUTDAT (2, 0, 'ISET', REAL (ISET))
         CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL, &
                      FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5, &
                      OUTPUT, TERMNL,  &
                      DOY   , IDOY  , YEAR  , IYEAR , STTIME, &
                      TIME  , DELT  , &
                      LAT   , WSTAT , WTRTER, &
                      RDD   , TMMN  , TMMX  , VP    , WN, RAIN)
      END IF

!----------------------------------------------------------------------*
!                             Time update                              *
!----------------------------------------------------------------------*

!-----Check for FINTIM, OUTPUT and observation days
      CALL TIMER2 (ITASK, STTIME, DELT, PRDEL, FINTIM, &
                   IYEAR, TIME  , DOY , IDOY , TERMNL, OUTPUT)
      YEAR = REAL (IYEAR)
      DO IOD=1,INOD,2
         IF (IYEAR.EQ.IOBSD(IOD).AND.IDOY.EQ.IOBSD(IOD+1)) &
             OUTPUT = .TRUE.
      END DO

      END DO
 
!======================================================================*
!                                                                      *
!                           Terminal section                           *
!                                                                      *
!======================================================================*

      ITASK = 4

      WRITE (*,'(A)') '   ORYZA3: Terminate model'

      !CALL ChartTerminalGroup

!-----Call routine that handles the different models
      CALL MODELS (ITASK , IUNITD, IUNITO, IUNITL, &
                   FILEIT, FILEI1, FILEI2, FILEI3, FILEI4, FILEI5, &
                   OUTPUT, TERMNL,  &
                   DOY   , IDOY  , YEAR  , IYEAR , STTIME,  &
                   TIME  , DELT  , &
                   LAT   , WSTAT , WTRTER, &
                   RDD   , TMMN  , TMMX  , VP    , WN, RAIN)

!-----Generate output file dependent on option from timer file
      IF (IPFORM.GE.4) THEN
         IF (INPRS.EQ.0) THEN
            CALL OUTDAT (IPFORM, 0, 'Simulation results',0.)
         ELSE
!           Selection of output variables was in timer file
!           write tables according to output selection array PRSEL
            CALL OUTSEL (PRSEL,IMNPRS,INPRS,IPFORM,'Simulation results')
         END IF
      END IF

      IF (WTRTER) THEN
         WRITE (*,'(/,A,/,/,/)') &
           ' The run was terminated due to missing weather'
         WRITE (IUNITO,'(/,A,/,/,/)') &
           ' The run was terminated due to missing weather'
         IF (IUNITO.NE.IUNITL) WRITE (IUNITL,'(/,A,/,/,/)') &
           ' The run was terminated due to missing weather'
      END IF

!-----Delete temporary output file dependent on switch from timer file
      IF (DELTMP.EQ.'Y'.OR.DELTMP.EQ.'y') CALL OUTDAT (99, 0, ' ', 0.)
	    IF((.NOT.KILLSOIL1).AND.(ISET.GE.STRUN)) THEN
		    PV%PDOY = DOY+1; PV%PYEAR = YEAR; PV%KILLSOIL = .FALSE.
		    IF(MOD(PV%PYEAR,4).GT.0.0) THEN
			    IF(PV%PDOY.GT.365) THEN
				    PV%PDOY=PV%PDOY-365
				    PV%PYEAR = PV%PYEAR+1
			    END IF
		    ELSE
			    IF(PV%PDOY.GT.366) THEN
				    PV%PDOY=PV%PDOY-366
				    PV%PYEAR = PV%PYEAR+1
			    END IF
		    END IF
		    ALLOCATE(TSI)
		    CALL SOIL_RESERVE("SET_VALUES")				
	    ELSEIF(ISET.EQ.ENDRUN) THEN
		    DEALLOCATE(pv)
	    ELSE
    !				DEALLOCATE(pv)	
	    END IF
    END DO

      IF (INSETS.GT.0) CLOSE (IUNITR)

!-----If input files should be copied to the output file,
!     copy rerun file (if present) and timer file and if there, input
!     files 1-5

      IF (COPINF.EQ.'Y'.OR.COPINF.EQ.'y') THEN
         IF (INSETS.GT.0) CALL COPFL2 (IUNITR, FILEIR, IUNITO, .TRUE.)
         CALL COPFL2 (IUNITD, FILEIT, IUNITO, .TRUE.)
         IF (FILEI1.NE.' ') CALL COPFL2 (IUNITD, FILEI1, IUNITO, .TRUE.)
         IF (FILEI2.NE.' ') CALL COPFL2 (IUNITD, FILEI2, IUNITO, .TRUE.)
         IF (FILEI3.NE.' ') CALL COPFL2 (IUNITD, FILEI3, IUNITO, .TRUE.)
         IF (FILEI4.NE.' ') CALL COPFL2 (IUNITD, FILEI4, IUNITO, .TRUE.)
         IF (FILEI5.NE.' ') CALL COPFL2 (IUNITD, FILEI5, IUNITO, .TRUE.)
      END IF

!-----Delete all .TMP files that were created by the RD* routines
!     during simulation
      CALL RDDTMP (IUNITD)

!-----Write to screen which files contain what
      IL = LEN_TRIM (FILEON)
      WRITE (*,'(/,3A)') ' File: ',FILEON(1:IL), &
        ' contains simulation results'
      WRITE (*,'(2A)') ' File: WEATHER.LOG', &
        ' contains messages from the weather system'
      IL = LEN_TRIM (FILEOL)
      WRITE (*,'(3A,/)') ' File: ',FILEOL(1:IL), &
        ' contains messages from the rest of the model'

!-----Write message to screen and output file if warnings and/or errors
!     have occurred from the weather system, pause and wait for return
!     from user to make sure he has seen this message

      IF (WTRMES) THEN

!         WRITE (*,'(/,A,/,A,/,A)') ' WARNING from ORYZA3:', &
!           ' Please check your model and weather log files!'
!         WRITE (IUNITO,'(A,/,A,/,A)') ' WARNING from ORYZA3:', &
!           ' Please check your model and weather log files!'

!         WRITE (*,'(A)') ' Press <Enter>'
!         READ  (*,'(A)') DUMMY

      END IF

!-----Close output file and temporary file of OUTDAT
      CLOSE (IUNITO)
      CLOSE (IUNITO+1)

!-----Close log file (if used)
      IF (FILEOL.NE.FILEON) CLOSE (IUNITL)

!-----Close log file of weather system
      IF((UNITWL.NE.91).AND.(UNITWL.NE.100)) THEN
        CLOSE(UNITWL)  !ADDED BY TAOLI TO CLOSE WEATHER FILE CORRECTLY
      END IF
      CLOSE (91)
	  CLOSE(100) 
!-----Write end_of_run values to file
      CALL OPWRITE (IUNITO,(INSETS+1))
      close(IUNITO)	 
      !@@@ Temporal section for testing DSSAT-ORYZA --- Start of section      
      !@@@CLOSE(600) 
       !@@@ Temporal section for testing DSSAT-ORYZA --- end of section   
!      call copy_op_weather(fileon)  !TAOLI 20 Aug 2013, Move the OP.DAT and weather.log to the same folder of RES.DAT, it
                                  !It is commented out for adoptation of Windows interface 
      RETURN
      END



	  SUBROUTINE SET_ZERO_PV

	  USE PUBLIC_MODULE
	
		pv%CRUN=0;				pv%PROOT_NUTRIENT=.FALSE.
		pv%pond_active ='NO';	pv%PYear=0;pv%Pdoy=0
		pv%Pdae=-1;				pv%Pdat=0
		pv%Pno3=0.0;			pv%Pnh4=0.0
		pv%PmineralizedN =0.0
		pv%Purea=0.0;			pv%pond_no3=0.0
		pv%pond_nh4=0.0;		pv%pond_urea=0.0
		pv%Psoc=0.0;			pv%Pson=0.0
		pv%pph=7.0;				pv%PFOM_type=0.0
		pv%PFOM_C=0.0;			pv%PFOM_N=0.0
		pv%Pnl=0;				pv%Pdlayer=0.0
		pv%Pbd=0.0;				pv%Psand=0.0
		pv%Pclay=0.0;			pv%Pkst=0.0
		pv%Pwcst=0.0;			pv%Pwcfc=0.0
		pv%Pwcwp=0.0;			pv%Pwcad=0.0
		pv%PplowDepth=0.0;		pv%Psoiltx=0.0
		pv%Pwl0=0.0;			pv%Pswc=0.0
		pv%Pwflux=0.0;			pv%PTRWL=0.0
		pv%Prunoff=0.0;			pv%Pdrain=0.0
		pv%Pirrig=0.0;			pv%Plai=0.001
		pv%PResNum=0;			pv%PResName=''
		pv%PResType=0;			pv%PResC=0.0
		pv%PResN=0.0;			pv%dlt_res_c_biom=0.0
		pv%dlt_res_c_hum=0.0;	pv%ProotD=0.0; pv%PRootAD = 0.0
		pv%ProotDen=0.0;		pv%PSROOTL=0.0
		pv%PRMINT=0.0;			pv%PROPTT=0.0
		pv%PRTBS=0.0;			pv%PRCNL=0.0
		pv%PMAXD=0.0;			pv%PSODT=0.0
		pv%PmaxT=0.0;			pv%PminT=0.0
		pv%PRad=0.0;			pv%PRain=0.0
		pv%PPressur=0.0;		pv%Pwind=0.0
		pv%PETp=0.0;			pv%PETa=0.0
		pv%Pevap=0.0;			pv%Ptrans=0.0
		pv%PDt=0.0;				pv%PRdn=0.0		
		PV%OPSTRING = ''		
	  
	  END SUBROUTINE SET_ZERO_PV


!-------------------------------------------------------------------------------------
!This routine will check whether the soil file is same in the rerun file
!The soil is same, if: 1) no soil file presented in rerun file;
!                      2) soil file presented in rerun file, but the full name is same
!It will be used when the cotinue rerun but soil process is carried on under:
!     1) with different rice varieties as rotation
!     2) with different field management
!     3) with different cropping system
!Developed: Dr. Tao Li, April 15, 2011
!-------------------------------------------------------------------------------------

SUBROUTINE CHECK_RERUN(RERUNFILE,ISETS,INSOILFILE,SAMESOIL)

    CHARACTER*(*) RERUNFILE, INSOILFILE
    LOGICAL SAMESOIL

    !LOCAL VARIABLES
    INTEGER I, II, IK
    INTEGER ISETS
    CHARACTER FILES*128, ST*1
    CHARACTER*(128), ALLOCATABLE:: SOILFILE(:)

    ALLOCATE(SOILFILE(ISETS))
    ! OPEN RERUN FILE
    
    OPEN(UNIT =500, FILE =RERUNFILE)
    II = 0
    DO WHILE (.NOT.EOF(500))
        READ(500,'(A)') FILES
        IF(LEN_TRIM(FILES).GT.0) THEN
            CALL UPPERC(FILES)
            FILES=ADJUSTL(TRIM(FILES))
            ST = FILES(1:1)
            IF((ST.EQ."*").OR.(ST.EQ."!")) THEN
            ELSE
                I = INDEX(FILES,"!")
                IF(I.GT.0) THEN
                    FILES = TRIM(FILES(1:(I-1)))
                END IF
                IF(INDEX(FILES,"FILEI2").GT.0) THEN
                    I=INDEX(FILES,"'")
                    IF(I.GT.0) FILES = ADJUSTL(TRIM(FILES((i+1):LEN_TRIM(FILES))))
                    I=INDEX(FILES,"'")
                    IF(I.GT.0) FILES = TRIM(FILES(1:(I-1)))
                    II= II+1
                    SOILFILE(II) = FILES
                END IF            
            END IF        
        END IF
    END DO
    
    CLOSE(500)
    IF(II.GT.0) THEN
        SAMESOIL = .TRUE.
        DO I=1, II
            DO IK=1, II
                IF(I.NE.IK) THEN
                    IF(SOILFILE(I).NE.SOILFILE(IK)) THEN
                        SAMESOIL =.FALSE.
                        EXIT
                    END IF
                END IF
            END DO
            IF(.NOT.SAMESOIL) EXIT
        END DO
    ELSE
        SAMESOIL = .TRUE.   
    END IF
    IF((SAMESOIL).AND.(II.GT.0)) THEN
        INSOILFILE = ADJUSTL(INSOILFILE)
        IF(LEN_TRIM(INSOILFILE).GT.0) THEN
            CALL UPPERC(SOILFILE(1))
            CALL UPPERC(INSOILFILE)
            IF (SOILFILE(1).EQ.INSOILFILE) THEN
				SAMESOIL = .TRUE.
			ELSE
				SAMESOIL =.FALSE.
			END IF
        ELSE
            SAMESOIL = .TRUE.
        END IF
	ELSEIF((II.EQ.0).AND.(LEN_TRIM(INSOILFILE).EQ.0)) THEN
		SAMESOIL = .FALSE.
    END IF
    DEALLOCATE(SOILFILE)
END SUBROUTINE CHECK_RERUN

!--------------------------------------------------------------------
!This subroutine will temporary STORE or SEND the soil physical and chemical 
!information if the soil is used for rotation running
!Developed" Dr. Tao Li, 18 April 2011
!
!-------------------------------------------------------------------- 
SUBROUTINE SOIL_RESERVE(ACTIONS)
	
	USE PUBLIC_module
	CHARACTER*(*) ACTIONS

	IF(INDEX(ACTIONS, "SET_VALUES").GT.0) THEN
		TSI%Xsno3= PV%PNO3;			TSI%Xsnh4=PV%PNH4;		TSI%Xurea=PV%PUREA
		TSI%Xpond_no3=PV%POND_NH4;	TSI%Xpond_nh4=PV%POND_NH4    
		TSI%Xpond_urea=PV%POND_UREA;TSI%Xsoc=PV%PSOC;		TSI%Xson=PV%PSON
		TSI%Xph=PV%PPH;				TSI%XFOM_type=PV%PFOM_type	
		TSI%XFOM_C=PV%PFOM_C;		TSI%XFOM_N=PV%PFOM_N;   TSI%XNL = PV%PNL
		TSI%Xdlayer=PV%PDLAYER;		TSI%Xbd=PV%PBD;			TSI%Xsand=PV%PSAND
		TSI%Xclay=PV%PCLAY;			TSI%Xkst=PV%PKST;		tSI%Xwcst=PV%PWCST		
		TSI%Xwcfc=PV%PWCFC;			TSI%Xwcwp=PV%PWCWP;		TSI%Xwcad=PV%PWCAD
		TSI%XplowDepth=PV%PPLOWDEPTH;						TSI%Xsoilt=PV%Psoiltx
		TSI%Xwl0 = PV%PWL0;			TSI%Xswc=PV%PSWC;		TSI%XResNum=PV%PRESNUM
		TSI%XResName=PV%PRESNAME;	TSI%XResType=PV%PRESTYPE
		TSI%XResC=PV%PRESC;			TSI%XResN=PV%PRESN
		TSI%Xdlt_res_c_biom=PV%dlt_res_c_biom;	TSI%Xdlt_res_c_hum=PV%dlt_res_c_hum     
	ELSEIF(INDEX(ACTIONS, "GET_VALUES").GT.0) THEN
		PV%Pno3= TSI%XSNO3;			PV%Pnh4=TSI%XSNH4;		PV%Purea=TSI%XUREA
		PV%pond_no3=TSI%XPOND_NH4;	PV%pond_nh4=TSI%XPOND_NH4    
		PV%pond_urea=TSI%XPOND_UREA;PV%Psoc=TSI%XSOC;		PV%Pson=TSI%XSON
		PV%Pph=TSI%XPH;				PV%PFOM_type=TSI%XFOM_type	
		PV%PFOM_C=TSI%XFOM_C;		PV%PFOM_N=TSI%XFOM_N;   PV%PNL = TSI%XNL
		PV%Pdlayer=TSI%XDLAYER;		PV%Pbd=TSI%XBD;			PV%Psand=TSI%XSAND
		PV%Pclay=TSI%XCLAY;			PV%Pkst=TSI%XKST;		PV%Pwcst=TSI%XWCST		
		PV%Pwcfc=TSI%XWCFC;			PV%Pwcwp=TSI%XWCWP;		PV%Pwcad=TSI%XWCAD
		PV%PplowDepth=TSI%XPLOWDEPTH;						PV%Psoiltx=TSI%XSOILT
		PV%Pwl0 = TSI%XWL0;			PV%Pswc=TSI%XSWC;		PV%PResNum=TSI%XRESNUM
		PV%PResName=TSI%XRESNAME;	PV%PResType=TSI%XRESTYPE
		PV%PResC=TSI%XRESC;			PV%PResN=TSI%XRESN
		PV%dlt_res_c_biom=TSI%Xdlt_res_c_biom;	PV%dlt_res_c_hum=TSI%Xdlt_res_c_hum    
	END IF
END SUBROUTINE

subroutine copy_OP_weather(fileop)
    USE IFPORT
    USE IFCORE
    character*(*) fileop
    character*(256) file,file1, file2, file3
    CHARACTER(3)        drive
     CHARACTER(256)      dir
     CHARACTER(256)      name
     CHARACTER(256)      ext     
    integer k, k1, k2
    logical kk
    
    if(len_trim(fileop).gt.0) then
        k=SPLITPATHQQ(fileop, drive, dir, name, ext)
        if(k.gt.0) then
            file = dir(1:k)
            k1= getcwd(file2)
            if(k1.eq.0) then
                file2=trim(file2)//'\op.dat'
                file1 = trim(drive)//trim(file)//'op.dat'
                CALL UPPERC(FILE2)
                CALL UPPERC(FILE1)
                IF(TRIM(FILE1).NE.TRIM(FILE2)) THEN
                    INQUIRE (FILE = file1, EXIST = kk)
                    if(kk) then
                        k2 = DELFILESQQ(file1)
                        !call DELFIL(file1, kk)
                    end if 
                    kk =RENAMEFILEQQ(file2, file1) 
                END IF                
            end if
             !move weather.log
            k1= getcwd(file3)
            if(k1.eq.0) then
                file3=trim(file3)//'\weather.log'
                file1 = trim(drive)//trim(file)//'weather.log'
                CALL UPPERC(FILE3)
                CALL UPPERC(FILE1)
                IF(TRIM(FILE1).NE.TRIM(FILE3)) THEN
                    INQUIRE (FILE = file1, EXIST = kk)
                    if(kk) then
                        k2 = DELFILESQQ(file1)
                        !call DELFIL(file1, kk)
                    end if 
                    kk =RENAMEFILEQQ(file3, file1) 
                END IF
        end if
    end if
        
    end if
end subroutine copy_OP_weather

SUBROUTINE PHENAD(LAT,ELEV, DOY)
!This subroutine is used the adjust the BASE, OPTIMAL AND MAXIMUM temperature of phenology
!development of a large spatial scale modeling which cover wide geolocations 
    USE public_module
    REAL LAT, ELEVV,DOY
    !LOCAL VARIABLES
    REAL D, W, X, L, X0, L0, W0

    !calculate latitude adjustment
    L= ABS(LAT);DOY=100
    L = L*3.1415926/180.0 !CONVERT TO RADIANS
    L0 = 23.5*3.1415926/180.0
    D=0.409*SIN(2.0*3.1415926/366*DOY-1.309)
    W = ACOS(-TAN(L)*TAN(D))
    W0 = ACOS(-TAN(L0)*TAN(D))
    X = W*SIN(L)*SIN(D)+COS(L)*COS(D)*SIN(W)
    X0 = W0*SIN(L0)*SIN(D)+COS(L0)*COS(D)*SIN(W0)
    PV%PLAALA = 0.12*(X0-X)/X0
    !CALCULATE ALTITUDE ADJUSTMENT
    IF(ELEV.GT.100.0) THEN
        PV%PLAALA =PV%PLAALA + (ELEV-100.0)/100.0*0.649
    END IF 
    END SUBROUTINE PHENAD
    
    INTEGER FUNCTION DAYSDIFF(YEAR1,DAY1,YEAR2,DAY2)
    INTEGER YEAR1,DAY1,YEAR2,DAY2
    
    IF(YEAR2.GT.YEAR1) THEN
        IF(MOD(REAL(YEAR1),4.0).EQ.0.0) THEN
            DAY2=DAY2+366
            DAYSDIFF=DAY2-DAY1
        ELSE
            DAY2=DAY2+365
            DAYSDIFF=DAY2-DAY1
        END IF
    ELSEIF(YEAR2.EQ.YEAR1) THEN
        DAYSDIFF=DAY2-DAY1    
    ELSE
        IF(MOD(REAL(YEAR2),4.0).EQ.0.0) THEN
            DAY1=DAY1+366
            DAYSDIFF=DAY2-DAY1
        ELSE
            DAY1=DAY1+365
            DAYSDIFF=DAY2-DAY1
        END IF
    END IF
    RETURN
END FUNCTION DAYSDIFF
