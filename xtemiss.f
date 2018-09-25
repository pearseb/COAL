c Commented out carbonaceous tracers as not used in Mk3L for the moment.
c BP  2010/08/06
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: xtemiss.f,v $
c Revision 1.10  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.9  2001/06/04 02:26:54  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.8  2001/03/07 04:28:58  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.7  2001/02/22 05:56:37  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.6  2000/12/08 03:58:52  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.5  2000/11/16 01:19:01  rot032
c Code to read corrected volcano data.
c
c Revision 1.4  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.3  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.2  2000/06/20 02:27:09  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.1  1999/12/02 02:56:11  rot032
c Initial revision
c

      SUBROUTINE XTEMISS
     &   (TWODT, KROW, PXTM1, P1MXTM1, TSM1M, SEAICEM, G3X01, APHP1, !Inputs
     &    LOLAND, PFOREST, PSNOW, PWLMX, WLM1M, WSMXM, WSM1M, ustar, !Inputs
     &    XTE, PXTEMS, so2dd, so4dd, sem)                          !Outputs
C
C    THIS ROUTINE CALCULATES THE LOWER BOUNDARY CONDITIONS
C    FOR VDIFF DEPENDING ON THE SURFACE EMISSION AND THE
C    DRY DEPOSITION FLUX.
C
C    JOHANN FEICHTER          UNI/HAMBURG         08/91
C    MODIFIED  U. SCHLESE    DKRZ-HAMBURG        JAN-95
C    Adapted for CSIRO GCM by Leon Rotstayn, 12/99
C
C    PURPOSE
C   ---------
C    THE LOWER BOUNDARY CONDITION FOR CALCULATING THE
C    TURBULENT EXCHANGE IN THE BOUNDARY LAYER IS
C    DETERMINED BY THE EMISSION AND THE DRY DEPOSITION FLUX.
C
C    INTERFACE
C   ------------
C    *XTEMISS* IS CALLED FROM *VDIFF*
C    Called from RADIN in CSIRO GCM, prior to CALL HVERTMX.
C

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f'

C Argument list
      REAL TWODT
      INTEGER KROW
      REAL PXTM1(KLON,KLEV,NTRAC) !Tracer mixing ratio at t-1 [kg/kg]
      REAL P1MXTM1(KLON)          !Density of air in surface layer
      real TSM1M(KLON)   !Surface temp
      real SEAICEM(KLON) !Sea-ice fraction
      real G3X01(KLON)          !10m wind (corrected to neutral for Nightingale scheme)
      real APHP1(KLON,KLEV+1)   !P at half levels at current timestep
      LOGICAL LOLAND(KLON)      !Land flag
      REAL PFOREST(KLON)        !Fractional vegetation cover
      REAL PSNOW(KLON)          !Snow depth [m]
      real ustar(klon)          !Friction velocty (m/s)
c Land-surface details needed to specify dry deposition velocity
      REAL PWLMX(KLON)   !maximum skin reservoir content of plant [ditto]
      real WLM1M(KLON)   !skin reservoir content of plant [mm for CSIRO GCM, not m]
      real WSMXM(KLON)   !field capacity of soil [ditto]
      real WSM1M(KLON)   !surface wetness [vol fraction for CSIRO GCM, not m]
      real XTE(KLON,KLEV,NTRAC) !Tracer tendencies (kg/kg/s)
      REAL PXTEMS(KLON,NTRAC)     !Sfc. flux of tracer passed to vertical mixing [kg/m2/s]
c Some diagnostics
      real so2dd(klon)
      real so4dd(klon)
      real sem(klon)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f' !Input TMELT from seaice scheme
      include 'ECAMEAN.f' !Input annual-mean source fields
      include 'ECFIELD.f' !Input big array of monthly means FIELD(LN2,NUMFL2,LAT)
      include 'FEWFLAGS.f'!Input debug stuff
      include 'TIMEX.f'   !Input MONTH

C Local work arrays and variables
C***  real SNMLTM(KLON)  !Area of melting snow?
c Annual-mean fields
      real ZVOLCANO(KLON)
      integer IVOLCHEIGHT(KLON)
      real ZBCSOURCE(KLON)
      real ZOCSOURCE(KLON)
      REAL ZVDRD(KLON,2)
      REAL ZMAXVDRY(KLON)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c      pcvdifts=1. !A time-stepping ratio used in ECHAM, not needed in CSIRO
      g=grav
      pxtems(:,:)=0.
      xte(:,:,:)=0.

c Indices in FIELD array need to change for 18-level CSIRO model, as 
c original code assumes 19 levels!

c These indices are 4 less than originally assumed

      iso2b=74
      iso2a1=75
      iso2a2=76
      idmso=77  
      ibcbi=78
      iocbi=79
c      iocna=116
      iocna=80

      do i=1,klon
        zvolcano(i)=vso2(i,krow)
        ivolcheight(i)=nint(hvolc(i,krow))
c        ivolcheight(i)=klev+1-nint(hvolc(i,krow)) !Because K values read in are in CSIRO sense.
        zbcsource(i)=bcsrc(i,krow)
        zocsource(i)=ocsrc(i,krow)
      enddo

c******************************************************************************

      ZTMST=TWODT
C
C     M WATER EQUIVALENT  CRITICAL SNOW HEIGHT (FROM *SURF*)
      ZSNCRI=0.025
C
C     COEFFICIENTS FOR ZVDRD = FUNCTION OF SOIL MOISTURE
C
      ZVWC2=(0.8E-2 - 0.2E-2)/(1. - 0.9)
      ZVW02=ZVWC2-0.8E-2
      ZVWC4=(0.2E-2 - 0.025E-2)/(1. - 0.9)
      ZVW04=ZVWC4-0.2E-2     
C
C --------------------------------------------------------------
C
C*     1.   SURFACE EMISSION.
C           ------- --------
C
C   CALCULATE DMS EMISSIONS FOLLOWING LISS+MERLIVAT
C   DMS SEAWATER CONC. FROM KETTLE ET AL.
      DO 112 JL=1,KLON
      IF (LOLAND(JL)) THEN
        ZDMSCON=0.
      ELSE
        ZDMSCON=FIELD(JL,idmso,KROW)*(1.-SEAICEM(JL))
c        ZDMSCON=FIELD(JL,idmso,KROW) !Sens test
        ZSST=TSM1M(JL)-273.15

c Andreae's Sc formula is not happy above 35 degrees, which can occur in q-flux runs!

        zsst=min(zsst, 35.)
C
C  G3X01:  10-M WINDS
C
        ZZSPEED=G3X01(JL)

c Nightingale (2000) scheme (J. Biogeochem. Cycles, 14, 373-387)
c For this scheme, zzspeed is the 10m wind adjusted to neutral stability.
c The formula for ScDMS from Saltzman et al (1993) is given by Kettle & Andreae (ref below)

C***        VpCO2 = 0.222*zzspeed**2 + 0.333*zzspeed
C***        ScCO2=600.
C***        ScDMS = 2674 - 147.12*zsst + 3.726*zsst**2 - 0.038*zsst**3 !Sc for DMS (Saltzman)
C***        if(zzspeed.lt.3.6)then
C***          zVdms = VpCO2 * (ScCO2/ScDMS)**(2./3)
C***        else
C***          zVdms = VpCO2 * sqrt(ScCO2/ScDMS)
C***        endif        

c Erickson (1993) scheme as described by Kettle and Andreae (2000), JGR 105 (D22), 26793 - 26808

C***        Fw = min (1., 0.2*ustar(jl)**3 )  !Whitecap fraction
C***        Vpw = 475.  !Piston velocity for Radon over whitecap area (cm/hr)
C***        Vpn = 9.58  !Piston velocity for Radon over non-whitecap area (cm/hr)
C***        VpRn = Fw * Vpw + (1. - Fw) * Vpn  !Piston velocity for Radon
C***        ScRn = 3147.3 - 201.9*zsst + 5.5*zsst**2 - 0.055*zsst**3  !Sc for Radon
C***        ScDMS = 2674 - 147.12*zsst + 3.726*zsst**2 - 0.038*zsst**3  !Sc for DMS
C***        if(zzspeed.lt.3.6)then
C***          zVdms = VpRn * (ScRn/ScDMS)**(2./3)
C***        else
C***          zVdms = VpRn * sqrt(ScRn/ScDMS)
C***        endif        

c Original ECHAM scheme...
        ZSCHMIDT=3652.047271-246.99*ZSST+8.536397*ZSST*ZSST
     *          -0.124397*ZSST*ZSST*ZSST  !Andreae's formula
        IF(ZZSPEED.GT.3.6.AND.ZZSPEED.LE.13.) THEN
            ZKW=2.85*ZZSPEED-9.65
            ZVDMS=ZKW*(ZSCHMIDT/600.)**(-0.5)
        ELSEIF(ZZSPEED.LE.3.6) THEN
            ZKW=0.17*ZZSPEED
            ZVDMS=ZKW*(ZSCHMIDT/600.)**(-2./3.)
        ELSE
            ZKW=5.9*ZZSPEED-49.3
            ZVDMS=ZKW*(ZSCHMIDT/600.)**(-0.5)
        ENDIF
C
        ZDMSEMISS=ZDMSCON*ZVDMS*32.064E-11/3600.
C   NANOMOL/LTR*CM/HOUR --> KG/M**2/SEC
C***        PXTEMS(JL,ITRACSO2-1)=PXTEMS(JL,ITRACSO2-1)+ZDMSEMISS
c Apply these terms directly, rather than as a surface flux via hvertmx.
        jk=klev
        xte(jl,jk,itracso2-1)=xte(jl,jk,itracso2-1)+
     *                 zdmsemiss*
     *                 g/(aphp1(jl,jk+1)-aphp1(jl,jk))
      ENDIF
C
C   EMISSIONS ON" FIELD":
c Subtract 4 from all of these, except for OC natural (subtract 42)
C                   
C                   NO 77 DMSTERR; 78 = SO2 BIOMASS BURN;  
C                   79= SO2ANTHR.LEVEL19; 80 = LEVEL18 (GEIA)
C                   81 = DMS SEAWATER CONC. (KETTLE ET AL.)
C                   82 = BC BIOMASS BURNING
C                   83 = OC BIOMASS BURNING
C                  122 = OC NATURAL
C   ANTHROPOGENIC SO2 AND BIOMASS BURNING

c Half of biomass emissions are in first layer and half in next layer up.

      PXTEMS(JL,ITRACSO2)=FIELD(JL,iso2b,KROW)*0.5
     &     +FIELD(JL,iso2a1,KROW)*0.97 !Formerly 78,79
      PXTEMS(JL,ITRACSO2+1)=FIELD(JL,iso2a1,KROW)*0.03

c Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via hvertmx.
c Also assume that 30% of the "surface" emissions really occur in the 2nd layer,
c either due to unresolved subgrid topography or the buoyancy of the smokestack plume. The
c 30% is added to the 2nd layer below.

      jk=klev
      gdp1=g/(aphp1(jl,jk+1)-aphp1(jl,jk))
      xte(jl,jk,itracso2)=xte(jl,jk,itracso2)
     &                      +pxtems(jl,itracso2)*gdp1*0.7
      xte(jl,jk,itracso2+1)=xte(jl,jk,itracso2+1)
     &                      +pxtems(jl,itracso2+1)*gdp1*0.7
  112 CONTINUE

C BP commented out this portion for debugging Mk3L because ntrac=3 only
C      if(ntrac.gt.3)then !Do carbonaceous aerosols
C       do 113 jl=1,klon
CC
Cc BC: 80% hydrophobic and 20% hydrophilic
C
C       PXTEMS(JL,ITRACBC)=0.8*(ZBCSOURCE(JL)+FIELD(JL,ibcbi,KROW))  !Formerly 82
C       PXTEMS(JL,ITRACBC+1)=0.2*(ZBCSOURCE(JL)+FIELD(JL,ibcbi,KROW))!Formerly 82 
C
Cc OC, 50% hydrophobic and 50% hydrophilic.
C
C       PXTEMS(JL,ITRACOC)=0.5*(ZOCSOURCE(JL)+FIELD(JL,iocbi,KROW)+  !Formerly 83
C     >  FIELD(JL,iocna,KROW))                                       !Formerly 122
C       PXTEMS(JL,ITRACOC+1)=0.5*(ZOCSOURCE(JL)+FIELD(JL,iocbi,KROW)+!Formerly 83
C     >  FIELD(JL,iocna,KROW))                                       !Formerly 122
C
Cc Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via hvertmx.
C
C      jk=klev
C      gdp1=g/(aphp1(jl,jk+1)-aphp1(jl,jk))
C      xte(jl,jk,itracbc)=xte(jl,jk,itracbc)+pxtems(jl,itracbc)*gdp1
C      xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)
C     &                      +pxtems(jl,itracbc+1)*gdp1
C      xte(jl,jk,itracoc)=xte(jl,jk,itracoc)+pxtems(jl,itracoc)*gdp1
C      xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)
C     &                      +pxtems(jl,itracoc+1)*gdp1
C 113   continue
C      endif !ntrac.gt.3    

C
C
C   END OF SURFACE EMISSIONS
C
C   SEA SALT CONCENTRATION
C
CLDR      CALL XTSEASALT( KLON, KLEV, G3XM40, G3X01, P1MXTM1, LOLAND)
C
C  EMISSION OF ANTHROPOGENIC SO2 IN THE NEXT HIGHER LEVEL PLUS BIOMASS BURNING
      JT=ITRACSO2
      DO 116 JL=1,KLON
C
        JK=KLEV-1
        XTE(JL,JK,JT)=XTE(JL,JK,JT)+
     *                (FIELD(JL,iso2a2,KROW)*0.97  !100% of the "above 100m" SO2 emission
     *                +pxtems(jl,jt)*0.3           !30% of the "surface" SO2 emission
     *                +FIELD(JL,iso2b,KROW)*0.5)   !50% of the biomass SO2 emission
     *                *G/(APHP1(JL,JK+1)-APHP1(JL,JK))
        XTE(JL,JK,JT+1)=XTE(JL,JK,JT+1)+
     *                (FIELD(JL,iso2a2,KROW)*0.03  !100% of the "above 100m" SO4 emission
     *                +pxtems(jl,jt+1)*0.3 )       !30% of the "surface" SO4 emission
     *                *G/(APHP1(JL,JK+1)-APHP1(JL,JK))

c 3rd layer

C***        JK=KLEV-2
C***        XTE(JL,JK,JT)=XTE(JL,JK,JT)+
C***     *                 FIELD(JL,iso2a2,KROW)*0.5*0.97  !50% of the "above 100m" SO2 emission
C***     *                *G/(APHP1(JL,JK+1)-APHP1(JL,JK))
C***        XTE(JL,JK,JT+1)=XTE(JL,JK,JT+1)+
C***     *                 FIELD(JL,iso2a2,KROW)*0.5*0.03  !50% of the "above 100m" SO4 emission
C***     *                *G/(APHP1(JL,JK+1)-APHP1(JL,JK))
C
C    VOLCANIC BACKGROUND EMISSIONS 
C
C   3 EMISSION LEVELS: 
C    1. PRE-INTRA ERUPTION IN LEVEL IVOLC-HEIGHT (=TOP OF VOLCANO)
C    2. POST-EXTRA ERUPTION IN LEVEL 15 -16 (CA 550-1736M)
C    3. EXPLOSIVE ERUPTION IN LEVEL 10 - 11 (CA 5000-7900M)
C
C       ZVOLCEMI  TOTAL EMISSION FROM VOLCANOES IN TG/YR
        ZVOLCEMI=8.
        ZVOLCEMI1=ZVOLCEMI*0.36
        ZVOLCEMI2=ZVOLCEMI*0.36
        ZVOLCEMI3=ZVOLCEMI*0.28
        JK=IVOLCHEIGHT(JL)
        IF(JK.GT.0) THEN
        ZDP=G/(APHP1(JL,JK+1)-APHP1(JL,JK))
        ZDP10=G/(APHP1(JL,11)-APHP1(JL,10))
        ZDP11=G/(APHP1(JL,12)-APHP1(JL,11))
        ZDP15=G/(APHP1(JL,16)-APHP1(JL,15))
        ZDP16=G/(APHP1(JL,17)-APHP1(JL,16))
        XTE(JL,JK,JT)=XTE(JL,JK,JT)+
     *                ZVOLCEMI1*ZVOLCANO(JL)*ZDP
        XTE(JL,16,JT)=XTE(JL,16,JT)+
     *                0.5*ZVOLCEMI2*ZVOLCANO(JL)*ZDP16
        XTE(JL,15,JT)=XTE(JL,15,JT)+
     *                0.5*ZVOLCEMI2*ZVOLCANO(JL)*ZDP15
        XTE(JL,11,JT)=XTE(JL,11,JT)+
     *                0.5*ZVOLCEMI3*ZVOLCANO(JL)*ZDP11
        XTE(JL,10,JT)=XTE(JL,10,JT)+
     *                0.5*ZVOLCEMI3*ZVOLCANO(JL)*ZDP10
        ENDIF       
  116 CONTINUE
C
C   --------------------------------------------------------------
C
C*      2.    DRY DEPOSITION.
C             --- ----------
      DO 205 JL=1,KLON
       ZMAXVDRY(JL)=(APHP1(JL,NLEV+1)-APHP1(JL,NLEV))/
     *              (G*P1MXTM1(JL)*ZTMST)
  205 CONTINUE
C
C      DRY DEPOSITION OF SO2, SO4
C
C
        ZVDPHOBIC=0.025E-2
        DO 206 JL=1,KLON
C     -  SEA -
         IF(.NOT.LOLAND(JL)) THEN
C         - SEA ICE -
C           - MELTING/NOT MELTING SEAICE-
              IF(TSM1M(JL).GE.(TMELT-0.1)) THEN
                ZVD2ICE=0.8E-2
                ZVD4ICE=0.2E-2
              ELSE
                ZVD2ICE=0.1E-2
                ZVD4ICE=0.025E-2
              ENDIF
c           ZVDRD(JL,1)=(1.-SEAICEM(JL))*1.0E-2+SEAICEM(JL)*ZVD2ICE
           ZVDRD(JL,1)=(1.-SEAICEM(JL))*0.8E-2+SEAICEM(JL)*ZVD2ICE !So leads agree with ocean
           ZVDRD(JL,2)=(1.-SEAICEM(JL))*0.2E-2+SEAICEM(JL)*ZVD4ICE
         ELSE
C      - LAND -
C        - NON-FOREST AREAS -
C         -  SNOW/NO SNOW -
             IF(PSNOW(JL).GT.ZSNCRI) THEN
C            - MELTING/NOT MELTING SNOW -
c               IF(SNMLTM(JL).GT.0.) THEN
               if(tsm1m(jl).ge.tmelt) then !This is a simplification of above line
                 ZVD2NOF=0.8E-2
                 ZVD4NOF=0.2E-2
               ELSE
                 ZVD2NOF=0.1E-2
                 ZVD4NOF=0.025E-2
               ENDIF
             ELSE
C           -  FROZEN/NOT FROZEN SOIL -
               IF(TSM1M(JL).LE.TMELT) THEN
                 ZVD2NOF=0.2E-2
                 ZVD4NOF=0.025E-2
               ELSE
C            - WET/DRY -
C               - COMPLETELY WET -
                 IF((WLM1M(JL)/PWLMX(JL)).GE.0.01
     *                             .OR.WSM1M(JL).EQ.WSMXM(JL)) THEN
                   ZVD2NOF=0.8E-2
                   ZVD4NOF=0.2E-2
                 ELSE
C                  - DRY -
                   IF((WSM1M(JL)/WSMXM(JL)).LT.0.9) THEN
                     ZVD2NOF=0.2E-2
                     ZVD4NOF=0.025E-2
                   ELSE
C                  - PARTLY WET -
                     ZVD2NOF=ZVWC2*(WSM1M(JL)/WSMXM(JL))-ZVW02
                     ZVD4NOF=ZVWC4*(WSM1M(JL)/WSMXM(JL))-ZVW04
                   ENDIF
                 ENDIF
               ENDIF
             ENDIF
             ZVDRD(JL,1)=PFOREST(JL)*0.8E-2+(1.-PFOREST(JL))*ZVD2NOF
             ZVDRD(JL,2)=PFOREST(JL)*0.2E-2+(1.-PFOREST(JL))*ZVD4NOF
           ENDIF
C
C Apply lower and upper bounds.

        ZVDRD(JL,1)=AMIN1(ZVDRD(JL,1),ZMAXVDRY(JL)) !SO2
        ZVDRD(JL,2)=AMIN1(ZVDRD(JL,2),ZMAXVDRY(JL)) !aerosols
  206 CONTINUE
C

c Sulfur emission diagnostic (hard-coded for 3 sulfur variables)

      pxtems(:,:)=0. !Zero this because fluxes are now passed in thru xte(:,:,:)

      sem(:)=pxtems(:,1)+pxtems(:,2)+pxtems(:,3) ! At surface

      do jt=1,3
        do jk=1,klev
          sem(:)=sem(:)+xte(:,jk,jt)*(aphp1(:,jk+1)-aphp1(:,jk))/g  !Above surface
        enddo
      enddo

      if(debug)then
        tdt=2*dt
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, In xtemiss, before dry dep'
          write(25,9)'del(SO2) ',(tdt*xte(mg,k,itracso2),k=nl,1,-1)
          write(25,9)'del(SO4) ',(tdt*xte(mg,k,itracso2+1),k=nl,1,-1)
          write(25,1)'sem*tdt ',tdt*sem(mg),
     &                'DMS surf emission ',pxtems(mg,1)*tdt,
     &                'SO2 surf emission ',pxtems(mg,2)*tdt,
     &                'SO4 surf emission ',pxtems(mg,3)*tdt
        endif
      endif

C
C
C      ZVDRD   DRY DEPOSITION VELOCITY IN M/S
C      ZVDRD(JL,1)  FOR SO2 GAS
C      ZVDRD(JL,2)  FOR AEROSOLS
C
      JT=ITRACSO2
      jk=klev
      ptmst=2*dt
      DO 211 JL=1,KLON
      gdp=g/(aphp1(jl,jk+1)-aphp1(jl,jk))
c      xt=pxtm1(jl,klev,itracso2)+0.5*ptmst*xte(jl,jk,itracso2) !Time-centered
      xt=pxtm1(jl,klev,itracso2)+ptmst*xte(jl,jk,itracso2) !Forward in time
      zhilso2=p1mxtm1(jl)*xt*zvdrd(JL,1)
      xt=pxtm1(jl,klev,itracso2+1)+ptmst*xte(jl,jk,itracso2+1) !Ditto
      zhilso4t=p1mxtm1(jl)*xt*zvdrd(jl,2)
      so2dd(jl)=zhilso2
      so4dd(jl)=zhilso4t
C***      PXTEMS(JL,JT)=PXTEMS(JL,JT)-ZHILSO2
C***      PXTEMS(JL,JT+1)=PXTEMS(JL,JT+1)-ZHILSO4T
      xte(jl,jk,jt)=xte(jl,jk,jt)-zhilso2*gdp
      xte(jl,jk,jt+1)=xte(jl,jk,jt+1)-zhilso4t*gdp
  211 CONTINUE

C BP commented out this portion for debugging Mk3L because ntrac=3 only
C      if(ntrac.gt.3)then
C        do jl=1,klon
C          gdp=g/(aphp1(jl,jk+1)-aphp1(jl,jk))
C          ZHILBCO=P1MXTM1(JL)*PXTM1(JL,KLEV,ITRACBC)*ZVDPHOBIC
C          ZHILBCY=P1MXTM1(JL)*PXTM1(JL,KLEV,ITRACBC+1)*ZVDRD(JL,2)
C          ZHILOCO=P1MXTM1(JL)*PXTM1(JL,KLEV,ITRACOC)*ZVDPHOBIC
C          ZHILOCY=P1MXTM1(JL)*PXTM1(JL,KLEV,ITRACOC+1)*ZVDRD(JL,2)
CC***          PXTEMS(JL,ITRACBC)=PXTEMS(JL,ITRACBC)-PCVDIFTS*ZHILBCO
CC***          PXTEMS(JL,ITRACBC+1)=PXTEMS(JL,ITRACBC+1)-PCVDIFTS*ZHILBCY
CC***          PXTEMS(JL,ITRACOC)=PXTEMS(JL,ITRACOC)-PCVDIFTS*ZHILOCO
CC***          PXTEMS(JL,ITRACOC+1)=PXTEMS(JL,ITRACOC+1)-PCVDIFTS*ZHILOCY
C          xte(jl,jk,itracbc)=xte(jl,jk,itracbc)-zhilbco*gdp
C          xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)-zhilbcy*gdp
C          xte(jl,jk,itracoc)=xte(jl,jk,itracoc)-zhiloco*gdp
C          xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)-zhilocy*gdp
C        enddo
C      endif

      if(debug)then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, end xtemiss'
          write(25,9)'aphp1 ',(aphp1(mg,k),k=nl+1,1,-1)
          write(25,9)'del(SO2) ',(tdt*xte(mg,k,itracso2),k=nl,1,-1)
          write(25,9)'del(SO4) ',(tdt*xte(mg,k,itracso2+1),k=nl,1,-1)
        endif
      endif
      
 1    format(4(a,g12.5))
 9    format(a,30g10.3)

      RETURN
      END
