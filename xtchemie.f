c Commented out carbonaceous tracers as not used in Mk3L for the moment.
c BP  2010/08/06
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: xtchemie.f,v $
c Revision 1.11  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.10  2001/06/04 02:27:01  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.9  2001/03/07 04:29:02  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.8  2001/02/28 04:40:57  rot032
c Bring sulfur cycle to run N58.
c
c Revision 1.7  2001/02/22 05:56:38  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.6  2000/12/08 03:58:54  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.5  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.4  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.3  2000/06/20 07:35:44  rot032
c Remove redundant declarations.
c
c Revision 1.2  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.1  2000/02/01 04:59:51  rot032
c Initial revision
c
      SUBROUTINE XTCHEMIE
     * (KTOP, KROW, PTMST, PDPP1, PMRATEP, PFPREC,  !Inputs
     *  PCLCOVER, PMLWC, PRHOP1, PTP1, sg, xtm1, pfevap,  !Inputs
     &  pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstay,pqfsed,plambs,  !Inputs
     &  prscav,pclcon,pccw,pfconv,pcscav,xtu,  !Inputs
     &  conwd,so2wd,so4wd,            !In and Out
     &  xte,so2oh,so2h2,so2o3,dmsoh,dmsn3)  !Outputs

c Inputs

c KROW: latitude index (LG in main program)      
c PTMST: timestep (seconds; tdt in main program)
c PDPP1: delta P (SI units)
c PMRATEP: precip formation rate (kg/kg/s)
c PFPREC: precipitation flux (entering from above) (kg/m2/s)
c PFEVAP: precipitation flux evaporating in layer k (kg/m2/s)
c PCLCOVER: cloud fraction (input; don't pass in cfrac though)
c PMLWC: liquid-water mixing ratio (kg/kg)
c PRHOP1: density of air (kg/m3)
c PTP1: temperature (K)
c sg: net solar radiation at ground (W/m2; used to determine if daytime)
c xtm1: tracer mixing ratio (kg/kg)
c xtu: tracer mixing ratio in convective updraught (kg/kg)
c pcscav: Convective tracer scavenging fraction (b/w 0 and 1)
c
c Outputs

c xte: tracer tendency (kg/kg/s)
C
C**** *XTCHEMIE*  CALCULATES DRY AND WET CHEMISTRY
C
C      J. FEICHTER             UNI HAMBURG    30/06/92
C
C      PURPOSE
C      ---------
C      THIS ROUTINE COMPUTES THE OXIDATION AND THE WET SCAVENGING
C      OF CHEMICAL SPECIES.
C
C**    INTERFACE
C      -----------
C      *XTCHEMIE*   IS CALLED FROM  progcld in CSIRO GCM
C
C      EXTERNALS
C      ------------
C          *XTWETDEP*  CALCULATES THE WET DEPOSITION
C

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f'

C Argument list
      INTEGER KTOP
      INTEGER KROW
      REAL PTMST
      REAL PDPP1(KLON,KLEV)
      REAL PMRATEP(KLON,KLEV)
      REAL PFPREC(KLON,KLEV)
      REAL PFEVAP(KLON,KLEV)
      REAL PCLCOVER(KLON,KLEV)
      REAL PMLWC(KLON,KLEV)
      REAL PRHOP1(KLON,KLEV)
      REAL PTP1(KLON,KLEV)
      real sg(klon)
      REAL XTM1(KLON,KLEV,NTRAC)
      real pcscav(klon,klev)    !Convective tracer scavenging fraction (b/w 0 and 1)
      real xtu(klon,klev,ntrac)
      REAL XTE(KLON,KLEV,NTRAC)
      real pfsnow(klon,klev)
      real pfconv(klon,klev)
      real pfsubl(klon,klev)
      real pcfcover(klon,klev)
      real pmiwc(klon,klev)
      real pmaccr(klon,klev)
      real pfmelt(klon,klev)
      real pfstay(klon,klev)
      real pqfsed(klon,klev)
      real plambs(klon,klev)
      real prscav(klon,klev)
      real pclcon(klon,klev)
      real pccw(klon,klev)
      real conwd(klon,ntrac)

      real dmsoh(klon) !Diagnostic output
      real dmsn3(klon) !Diagnostic output
      real so2wd(klon) !Diagnostic output
      real so4wd(klon) !Diagnostic output
      real so2oh(klon) !Diagnostic output
      real so2h2(klon) !Diagnostic output
      real so2o3(klon) !Diagnostic output

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'ECFIELD.f' !Input big array of monthly means ZOXIDANT(LN2,NUMFL2,LAT)
      include 'FEWFLAGS.f'

C Local work arrays and variables
      real so2oh3d(klon,klev),dmsoh3d(klon,klev),dmsn33d(klon,klev)
      real stratwd(klon) !Diagnostic output
C
      REAL ZXTP10(KLON,KLEV),          ZXTP1C(KLON,KLEV),
     *     ZHENRY(KLON,KLEV),
     *     ZSO4(KLON,KLEV),            ZRKH2O2(KLON,KLEV),
     *     ZSO4i(KLON,KLEV),           ZSO4C(KLON,KLEV),
     *     ZHENRYC(KLON,KLEV),         ZXTP1CON(KLON,KLEV)
      REAL ZZOH(KLON,KLEV),            ZZH2O2(KLON,KLEV),
     *     ZZO3(KLON,KLEV)
     *    ,ZZNO2(KLON,KLEV),           ZDXTE(KLON,KLEV,NTRAC)
      REAL ZDEP3D(KLON,KLEV),
     *     zlwcic(klon,klev)
c      real ziwcic(klon,klev)
      real zrevap(klon,klev)
c      real zso2ev(klon,klev)
      real xto(klon,klev,ntrac)
      integer ZRDAYL(klon)
      real zdayfac(2) !Two hemispheric values 
      parameter (nfastox=0) !1 for "fast" in-cloud oxidation; 0 for "slow"

C Local data, functions etc
C
C     DEFINE FUNCTION FOR CHANGING THE UNITS
C     FROM MASS-MIXING RATIO TO MOLECULES PER CM**3 AND VICE VERSA
      XTOC(X,Y)=X*6.022E+20/Y
      CTOX(X,Y)=Y/(6.022E+20*X)
C   X = DENSITY OF AIR, Y = MOL WEIGHT IN GRAMM
C
      ZFARR(ZK,ZH,ZTPQ)=ZK*EXP(ZH*ZTPQ)

C Start code : ----------------------------------------------------------

      if(debug)then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, start xtchemie'
          write(25,9)'pclcov ',(pclcover(mg,k),k=nl,1,-1)
          write(25,9)'DMS ',(xtm1(mg,k,1),k=nl,1,-1)
          write(25,9)'SO2 ',(xtm1(mg,k,2),k=nl,1,-1)
          write(25,9)'SO4 ',(xtm1(mg,k,3),k=nl,1,-1)
        endif
      endif

      dmsoh(:)=0.
      dmsn3(:)=0.
      so2oh(:)=0.
      so2h2(:)=0.
      so2o3(:)=0.
      so2oh3d(:,:)=0.
      dmsoh3d(:,:)=0.
      dmsn33d(:,:)=0.
      xte(:,:,:)=0.
      pg=grav
      pcons2=1./(ptmst*pg)
      pdtime=0.5*ptmst
      do i=1,klon
        if(sg(i).gt.0.)then
          zrdayl(i)=1
        else
          zrdayl(i)=0
        endif
      enddo

c Calculate xto, tracer mixing ratio outside convective updraughts
c Assumes pclcon < 1, but this shouldn't be a problem.
      
      do jt=1,ntrac
        do jk=1,klev
          do jl=1,klon
            xto(jl,jk,jt)=(xtm1(jl,jk,jt)-pclcon(jl,jk)*xtu(jl,jk,jt))
     &              /(1.-pclcon(jl,jk))
            xto(jl,jk,jt)=max(0.,xto(jl,jk,jt))
          enddo
        enddo
      enddo

C***c Simplicity: assume xto and xtu both equal grid-mean value.
C***      
C***      xto(:,:,:)=xtm1(:,:,:)
C***      xtu(:,:,:)=xtm1(:,:,:)

C
      CONTINUE
C
C   CALCULATE THE ZRDAYL (=0 --> NIGHT; =1 --> DAY) AND
C                 ZAMUO  =  ZENITH ANGLE
C
C    CONSTANTS
      PQTMST=1./PTMST
      ZMIN=1.E-20
      if(nfastox.eq.0)then
        NITER=5  !Slow in-cloud oxidation
      else 
        NITER=1  !Fast
      endif
C
C    REACTION RATE SO2-OH
      ZK2I=2.0E-12
      ZK2=4.0E-31
      ZK2F=0.45
C   REACTION RATE DMS-NO3
      ZK3=1.9E-13
C   MOLECULAR WEIGHTS IN G
      ZMOLGS=32.064
      ZMOLGAIR=28.84
C
      ZHPBASE=2.5E-06
      ZE1K=1.1E-02
      ZE1H=2300.
      ZE2K=1.23
      ZE2H=3020.
      ZE3K=1.2E-02
      ZE3H=2010.
      ZQ298=1./298.
      ZRGAS=8.2E-02
C
      ZAVO=6.022E+23
      ZNAMAIR=1.E-03*ZAVO/ZMOLGAIR
C
      CONTINUE

c Calculate in-cloud ql

      do jk=1,klev
        do jl=1,klon
          if(pclcover(jl,jk).gt.zmin)then
            zlwcic(jl,jk)=pmlwc(jl,jk)/pclcover(jl,jk)
          else
            zlwcic(jl,jk)=0.
          endif
c          if(pcfcover(jl,jk).gt.zmin)then
c            ziwcic(jl,jk)=pmiwc(jl,jk)/pcfcover(jl,jk)
c          else
c            ziwcic(jl,jk)=0.
c          endif
        enddo
      enddo

C  OXIDANT CONCENTRATIONS IN MOLECULE/CM**3
C
        DO 112 JK=1,KLEV
          DO 111 JL=1,KLON
           ZX=PRHOP1(JL,JK)*1.E-03
           JS1=JK
           JS2=KLEV+JK
           JS3=2*KLEV+JK
           JS4=3*KLEV+JK
           ZZOH(JL,JK)=ZOXIDANT(JL,JS1,KROW)
           ZZH2O2(JL,JK)=ZOXIDANT(JL,JS2,KROW)*ZX
           ZZO3(JL,JK)=ZOXIDANT(JL,JS3,KROW)*ZX
           ZZNO2(JL,JK)=ZOXIDANT(JL,JS4,KROW)*ZX
  111     CONTINUE
  112   CONTINUE
C
      CALL RESETR(ZHENRY,KLON*KLEV,0.)
      CALL RESETR(ZHENRYC,KLON*KLEV,0.)
      CALL RESETR(ZDXTE,KLON*KLEV*NTRAC,0.)
C
      CONTINUE

C
C   PROCESSES WHICH ARE DIFERENT INSIDE AND OUTSIDE OF CLOUDS
      CONTINUE
C
      JT=ITRACSO2+1
      DO 304 JK=KTOP,KLEV
        DO 302 JL=1,KLON
          ZSO4(JL,JK)=XTO(JL,JK,JT)
          ZSO4(JL,JK)=AMAX1(0.,ZSO4(JL,JK))
  302   CONTINUE
  304 CONTINUE
C
C
C   CALCULATE THE REACTION-RATES FOR SO2-H2O2
      DO 308 JK=KTOP,KLEV
       DO 306 JL=1,KLON
        IF(ZLWCIC(JL,JK).GT.ZMIN) THEN
         ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
         ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
         ZHP=ZHPBASE+ZSO4(JL,JK)*1000./(ZLWCIC(JL,JK)*ZMOLGS)
         ZQTP1=1./PTP1(JL,JK)-ZQ298
         ZRK=ZFARR(8.E+04,-3650.,ZQTP1)/(0.1+ZHP)
         ZRKE=ZRK/(ZLWCL*ZAVO)
C
         ZH_SO2=ZFARR(ZE2K,ZE2H,ZQTP1)
         ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
         ZP_SO2=ZH_SO2*ZPFAC
         ZF_SO2=ZP_SO2/(1.+ZP_SO2)
C
         ZH_H2O2=ZFARR(9.7E+04,6600.,ZQTP1)
         ZP_H2O2=ZH_H2O2*ZPFAC
         ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)
C
         ZRKH2O2(JL,JK)=ZRKE*ZF_SO2*ZF_H2O2
        ELSE
         ZRKH2O2(JL,JK)=0.
        ENDIF
  306  CONTINUE
  308 CONTINUE
C
C
C   HETEROGENEOUS CHEMISTRY
      JT=ITRACSO2
      DO 322 JK=KTOP,KLEV
       DO 320 JL=1,KLON
        ZXTP1=XTO(JL,JK,JT)
        ZXTP10(JL,JK)=XTO(JL,JK,JT)
        ZXTP1C(JL,JK)=XTO(JL,JK,JT)
        IF(ZXTP1.GT.ZMIN.AND.ZLWCIC(JL,JK).GT.ZMIN) THEN
          X=PRHOP1(JL,JK)
C
            ZQTP1=1./PTP1(JL,JK)-ZQ298
            ZE1=ZFARR(ZE1K,ZE1H,ZQTP1)
            ZE2=ZFARR(ZE2K,ZE2H,ZQTP1)
            ZE3=ZFARR(ZE3K,ZE3H,ZQTP1)
C
            ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
C    ZLWCL = LWC IN L/CM**3
            ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
C   ZLWCV = LWC IN VOL/VOL
            ZFAC1=1./(ZLWCL*ZAVO)
C   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
            ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
C   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
            ZZA=ZE2*ZRKFAC
            ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
            ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
            ZPH_O3=ZE1*ZRKFAC
            ZF_O3=ZPH_O3/(1.+ZPH_O3)
            ZDT=PTMST/FLOAT(NITER)
C
            ZH2O2M=ZZH2O2(JL,JK)
            ZSO2M=ZXTP1*XTOC(X,ZMOLGS)
            ZSO4M=ZSO4(JL,JK)*XTOC(X,ZMOLGS)
C
            ZSUMH2O2=0.
            ZSUMO3=0.
C
             DO 309 JN=1,NITER
              ZQ=ZRKH2O2(JL,JK)*ZH2O2M
              ZSO2MH=(1-nfastox)*ZSO2M*EXP(-ZQ*ZDT) ! = zero if nfastox.eq.1
     &                +nfastox * max (0., zso2m - zh2o2m )
C
              ZDSO2H=ZSO2M-ZSO2MH
              ZH2O2M=ZH2O2M-ZDSO2H
              ZH2O2M=AMAX1(0.,ZH2O2M)
              ZSUMH2O2=ZSUMH2O2+ZDSO2H
C
              ZSO4M=ZSO4M+ZDSO2H
C   CALCULATE THE PH VALUE
              ZSO2L=ZSO2MH*ZFAC1
              ZSO4L=ZSO4M*ZFAC1
              ZZB=ZHPBASE+ZSO4L
              ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
              ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
              ZZP=0.5*ZZP
              ZZP2=ZZP*ZZP
              ZHP=-ZZP+SQRT(ZZP2-ZZQ)
              ZQHP=1./ZHP
C
C   CALCULATE THE REACTION RATE FOR SO2-O3
C
              ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
              ZHENEFF=1.+ZE3*ZQHP
              ZP_SO2=ZZA*ZHENEFF
              ZF_SO2=ZP_SO2/(1.+ZP_SO2)
              ZRKO3=ZA2*ZF_O3*ZF_SO2
C
              ZQ=ZZO3(JL,JK)*ZRKO3
              ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
              ZDSO2O=ZSO2MH-ZSO2MO
              ZSO4M=ZSO4M+ZDSO2O
              ZSO2M=ZSO2MO
              ZSUMO3=ZSUMO3+ZDSO2O
  309       CONTINUE  !End of iteration loop
C
            ZDSO2TOT=ZXTP1-ZSO2M*CTOX(X,ZMOLGS)
            ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
            ZXTP1C(JL,JK)=ZXTP1-ZDSO2TOT
            ZSO4(JL,JK)=ZSO4(JL,JK)+ZDSO2TOT
C
            ZHENRY(JL,JK)=ZF_SO2
c Diagnostic only...
            ZFAC=PQTMST*CTOX(X,ZMOLGS)*PCLCOVER(JL,JK)
            ZFAC1=ZFAC*PDPP1(JL,JK)/PG
            so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
            so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
         ENDIF
  320  CONTINUE
  322 CONTINUE


c Repeat the aqueous oxidation calculation for ice clouds.

      JT=ITRACSO2+1
      DO JK=KTOP,KLEV
        DO JL=1,KLON
          ZSO4i(JL,JK)=XTO(JL,JK,JT)
          ZSO4i(JL,JK)=AMAX1(0.,ZSO4i(JL,JK))
        ENDDO
      ENDDO

c******************************************************************************
C***c Comment out from here when not using oxidation in ice clouds...
C***C
C***C
C***C   CALCULATE THE REACTION-RATES FOR SO2-H2O2
C***      DO JK=KTOP,KLEV
C***       DO JL=1,KLON
C***        IF(ziwcic(JL,JK).GT.ZMIN) THEN
C***         ZLWCL=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-06
C***         ZLWCV=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-03
C***         ZHP=ZHPBASE+ZSO4i(JL,JK)*1000./(ziwcic(JL,JK)*ZMOLGS)
C***         ZQTP1=1./PTP1(JL,JK)-ZQ298
C***         ZRK=ZFARR(8.E+04,-3650.,ZQTP1)/(0.1+ZHP)
C***         ZRKE=ZRK/(ZLWCL*ZAVO)
C***C
C***         ZH_SO2=ZFARR(ZE2K,ZE2H,ZQTP1)
C***         ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
C***         ZP_SO2=ZH_SO2*ZPFAC
C***         ZF_SO2=ZP_SO2/(1.+ZP_SO2)
C***C
C***         ZH_H2O2=ZFARR(9.7E+04,6600.,ZQTP1)
C***         ZP_H2O2=ZH_H2O2*ZPFAC
C***         ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)
C***C
C***         ZRKH2O2(JL,JK)=ZRKE*ZF_SO2*ZF_H2O2
C***        ELSE
C***         ZRKH2O2(JL,JK)=0.
C***        ENDIF
C***       ENDDO
C***      ENDDO
C***C
C***C
C***C   HETEROGENEOUS CHEMISTRY
C***      JT=ITRACSO2
C***      DO JK=KTOP,KLEV
C***       DO JL=1,KLON
C***        ZXTP1=XTO(JL,JK,JT)
C***        IF(ZXTP1.GT.ZMIN.AND.ziwcic(JL,JK).GT.ZMIN) THEN
C***          X=PRHOP1(JL,JK)
C***C
C***            ZQTP1=1./PTP1(JL,JK)-ZQ298
C***            ZE1=ZFARR(ZE1K,ZE1H,ZQTP1)
C***            ZE2=ZFARR(ZE2K,ZE2H,ZQTP1)
C***            ZE3=ZFARR(ZE3K,ZE3H,ZQTP1)
C***C
C***            ZLWCL=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-06
C***C    ZLWCL = LWC IN L/CM**3
C***            ZLWCV=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-03
C***C   ZLWCV = LWC IN VOL/VOL
C***            ZFAC1=1./(ZLWCL*ZAVO)
C***C   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
C***            ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
C***C   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
C***            ZZA=ZE2*ZRKFAC
C***            ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
C***            ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
C***            ZPH_O3=ZE1*ZRKFAC
C***            ZF_O3=ZPH_O3/(1.+ZPH_O3)
C***            ZDT=PTMST/FLOAT(NITER)
C***C
C***            ZH2O2M=ZZH2O2(JL,JK)
C***            ZSO2M=ZXTP1*XTOC(X,ZMOLGS)
C***            ZSO4M=ZSO4i(JL,JK)*XTOC(X,ZMOLGS)
C***C
C***            ZSUMH2O2=0.
C***            ZSUMO3=0.
C***C
C***             DO JN=1,NITER
C***              ZQ=ZRKH2O2(JL,JK)*ZH2O2M
C***              ZSO2MH=(1-nfastox)*ZSO2M*EXP(-ZQ*ZDT) ! = zero if nfastox.eq.1
C***     &                +nfastox * max (0., zso2m - zh2o2m )
C***C
C***              ZDSO2H=ZSO2M-ZSO2MH
C***              ZH2O2M=ZH2O2M-ZDSO2H
C***              ZH2O2M=AMAX1(0.,ZH2O2M)
C***              ZSUMH2O2=ZSUMH2O2+ZDSO2H
C***C
C***              ZSO4M=ZSO4M+ZDSO2H
C***C   CALCULATE THE PH VALUE
C***              ZSO2L=ZSO2MH*ZFAC1
C***              ZSO4L=ZSO4M*ZFAC1
C***              ZZB=ZHPBASE+ZSO4L
C***              ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
C***              ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
C***              ZZP=0.5*ZZP
C***              ZZP2=ZZP*ZZP
C***              ZHP=-ZZP+SQRT(ZZP2-ZZQ)
C***              ZQHP=1./ZHP
C***C
C***C   CALCULATE THE REACTION RATE FOR SO2-O3
C***C
C***              ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
C***              ZHENEFF=1.+ZE3*ZQHP
C***              ZP_SO2=ZZA*ZHENEFF
C***              ZF_SO2=ZP_SO2/(1.+ZP_SO2)
C***              ZRKO3=ZA2*ZF_O3*ZF_SO2
C***C
C***              ZQ=ZZO3(JL,JK)*ZRKO3
C***              ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
C***              ZDSO2O=ZSO2MH-ZSO2MO
C***              ZSO4M=ZSO4M+ZDSO2O
C***              ZSO2M=ZSO2MO
C***              ZSUMO3=ZSUMO3+ZDSO2O
C***            ENDDO  !End of iteration loop
C***C
C***            ZDSO2TOT=ZXTP1-ZSO2M*CTOX(X,ZMOLGS)
C***            ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
C***
C***            ZXTP10(JL,JK)=ZXTP1-ZDSO2TOT
C***     &                   *pcfcover(jl,jk)/(1.-pclcover(jl,jk))
C***            ZSO4i(JL,JK)=ZSO4i(JL,JK)+ZDSO2TOT
C***     &                   *pcfcover(jl,jk)/(1.-pclcover(jl,jk))
C***c            ZHENRY(JL,JK)=ZF_SO2
C***c Diagnostic only...
C***            ZFAC=PQTMST*CTOX(X,ZMOLGS)*pcfcover(jl,jk)
C***            ZFAC1=ZFAC*PDPP1(JL,JK)/PG
C***            so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
C***            so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
C***         ENDIF
C***       ENDDO
C***      ENDDO
c******************************************************************************

c Repeat the aqueous oxidation calculation for convective clouds.

      JT=ITRACSO2+1
      DO JK=KTOP,KLEV
        DO JL=1,KLON
          ZXTP1CON(JL,JK)=XTU(JL,JK,JT-1)
          ZSO4C(JL,JK)=XTU(JL,JK,JT)
          ZSO4C(JL,JK)=AMAX1(0.,ZSO4C(JL,JK))
        ENDDO
      ENDDO

c Comment from here when not using convective oxidation...

C   CALCULATE THE REACTION-RATES FOR SO2-H2O2
      DO JK=KTOP,KLEV
       DO JL=1,KLON
        IF(PCCW(JL,JK).GT.ZMIN) THEN
         ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
         ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
         ZHP=ZHPBASE+ZSO4C(JL,JK)*1000./(PCCW(JL,JK)*ZMOLGS)
         ZQTP1=1./PTP1(JL,JK)-ZQ298
         ZRK=ZFARR(8.E+04,-3650.,ZQTP1)/(0.1+ZHP)
         ZRKE=ZRK/(ZLWCL*ZAVO)
C
         ZH_SO2=ZFARR(ZE2K,ZE2H,ZQTP1)
         ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
         ZP_SO2=ZH_SO2*ZPFAC
         ZF_SO2=ZP_SO2/(1.+ZP_SO2)
C
         ZH_H2O2=ZFARR(9.7E+04,6600.,ZQTP1)
         ZP_H2O2=ZH_H2O2*ZPFAC
         ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)
C
         ZRKH2O2(JL,JK)=ZRKE*ZF_SO2*ZF_H2O2
        ELSE
         ZRKH2O2(JL,JK)=0.
        ENDIF
       ENDDO
      ENDDO
C

      if(debug)then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, Before aqueous oxidation'
          write(25,9)'zso4c ',(zso4c(mg,k),k=nl,1,-1)
          write(25,9)'zxtp1con ',(zxtp1con(mg,k),k=nl,1,-1)
        endif
      endif

C
C   HETEROGENEOUS CHEMISTRY
      JT=ITRACSO2
      DO JK=KTOP,KLEV
       DO JL=1,KLON
        ZXTP1=XTU(JL,JK,JT)
        IF(ZXTP1.GT.ZMIN.AND.PCCW(JL,JK).GT.ZMIN) THEN
          X=PRHOP1(JL,JK)
C
            ZQTP1=1./PTP1(JL,JK)-ZQ298
            ZE1=ZFARR(ZE1K,ZE1H,ZQTP1)
            ZE2=ZFARR(ZE2K,ZE2H,ZQTP1)
            ZE3=ZFARR(ZE3K,ZE3H,ZQTP1)
C
            ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
C    ZLWCL = LWC IN L/CM**3
            ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
C   ZLWCV = LWC IN VOL/VOL
            ZFAC1=1./(ZLWCL*ZAVO)
C   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
            ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
C   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
            ZZA=ZE2*ZRKFAC
            ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
            ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
            ZPH_O3=ZE1*ZRKFAC
            ZF_O3=ZPH_O3/(1.+ZPH_O3)
            ZDT=PTMST/FLOAT(NITER)
C
            ZH2O2M=ZZH2O2(JL,JK)
            ZSO2M=ZXTP1*XTOC(X,ZMOLGS)
            ZSO4M=ZSO4C(JL,JK)*XTOC(X,ZMOLGS)
C
            ZSUMH2O2=0.
            ZSUMO3=0.
C
             DO JN=1,NITER
              ZQ=ZRKH2O2(JL,JK)*ZH2O2M
              ZSO2MH=(1-nfastox)*ZSO2M*EXP(-ZQ*ZDT) ! = zero if nfastox.eq.1
     &                +nfastox * max (0., zso2m - zh2o2m )
C
              ZDSO2H=ZSO2M-ZSO2MH
              ZH2O2M=ZH2O2M-ZDSO2H
              ZH2O2M=AMAX1(0.,ZH2O2M)
              ZSUMH2O2=ZSUMH2O2+ZDSO2H
C
              ZSO4M=ZSO4M+ZDSO2H
C   CALCULATE THE PH VALUE
              ZSO2L=ZSO2MH*ZFAC1
              ZSO4L=ZSO4M*ZFAC1
              ZZB=ZHPBASE+ZSO4L
              ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
              ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
              ZZP=0.5*ZZP
              ZZP2=ZZP*ZZP
              ZHP=-ZZP+SQRT(ZZP2-ZZQ)
              ZQHP=1./ZHP
C
C   CALCULATE THE REACTION RATE FOR SO2-O3
C
              ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
              ZHENEFF=1.+ZE3*ZQHP
              ZP_SO2=ZZA*ZHENEFF
              ZF_SO2=ZP_SO2/(1.+ZP_SO2)
              ZRKO3=ZA2*ZF_O3*ZF_SO2
C
              ZQ=ZZO3(JL,JK)*ZRKO3
              ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
              ZDSO2O=ZSO2MH-ZSO2MO
              ZSO4M=ZSO4M+ZDSO2O
              ZSO2M=ZSO2MO
              ZSUMO3=ZSUMO3+ZDSO2O
            ENDDO  !End of iteration loop
C
            ZDSO2TOT=ZXTP1-ZSO2M*CTOX(X,ZMOLGS)
            ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
            ZXTP1CON(JL,JK)=ZXTP1CON(JL,JK)-ZDSO2TOT
            ZSO4C(JL,JK)=ZSO4C(JL,JK)+ZDSO2TOT
            ZHENRYC(JL,JK)=ZF_SO2

c Diagnostic only...
            ZFAC=PQTMST*CTOX(X,ZMOLGS)*pclcon(jl,jk)
            ZFAC1=ZFAC*PDPP1(JL,JK)/PG
            so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
            so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
         ENDIF
       ENDDO
      ENDDO

      if(debug)then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, After aqueous oxidation'
          write(25,9)'zhenryc ',(zhenryc(mg,k),k=nl,1,-1)
          write(25,9)'zso4c ',(zso4c(mg,k),k=nl,1,-1)
          write(25,9)'zxtp1con ',(zxtp1con(mg,k),k=nl,1,-1)
        endif
      endif

c*******************************************************************************
C
C    CALCULATE THE WET DEPOSITION
C
C
      DO 335 JT=ITRACSO2,NTRAC
      CALL RESETR(ZDEP3D,KLON*KLEV,0.)
      IF (LWETDEP(JT)) THEN
      IF(JT.EQ.ITRACSO2+1) THEN
        DO 332 JK=KTOP,KLEV
          DO 330 JL=1,KLON
            ZXTP1C(JL,JK)=ZSO4(JL,JK)
            ZXTP10(JL,JK)=zso4i(JL,JK)
            ZXTP1CON(JL,JK)=zso4C(JL,JK)
  330     CONTINUE
  332   CONTINUE
        zhenry (:,:)=1.
c        zhenry (:,:)=0.6
        zhenryc(:,:)=1.
      ENDIF
      IF(JT.GE.ITRACBC+1) THEN
        DO 334 JK=KTOP,KLEV
          DO 333 JL=1,KLON
            ZXTP10(JL,JK)=XTO(JL,JK,JT)
            ZXTP1C(JL,JK)=XTO(JL,JK,JT)
            ZXTP1CON(JL,JK)=XTU(JL,JK,JT)
  333     CONTINUE
  334   CONTINUE
        zhenry (:,:)=1.
        zhenryc(:,:)=1.
      ENDIF

      if(debug.and.(jt.eq.itracso2.or.jt.eq.itracso2+1))then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, before xtwetdep, JT= ',jt
          write(25,9)'conwd ',conwd(mg,jt)*ptmst
          write(25,9)'zxtp10 ',(zxtp10(mg,k),k=nl,1,-1)
          write(25,9)'zxtp1c ',(zxtp1c(mg,k),k=nl,1,-1)
          write(25,9)'zxtp1con ',(zxtp1con(mg,k),k=nl,1,-1)
        endif
      endif

C
      JTRAC=JT
      CALL XTWETDEP( JTRAC,krow,
     *               PTMST, PCONS2, PDTIME,
     *               PDPP1,
     *               PMRATEP, PFPREC, PFEVAP,
     *               PCLCOVER, PRHOP1, ZHENRY, zhenryc, pmlwc,
     &               pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,
     &               pfstay,pqfsed,plambs,prscav,pfconv,pclcon,pcscav, !Inputs
     &               ZXTP10, ZXTP1C,ZDEP3D,conwd,zxtp1con,        !In and Out
     &               zrevap )

c Reevaporate the SO2 as SO4
c Convective bit only in this version...
      
C***      if(jt.eq.itracso2)then
C***        zso2ev(:,:)=zrevap(:,:) !Save this for SO4 step
C***        do jk=1,klev
C***          do jl=1,klon
C***            if(zso2ev(jl,jk).gt.zmin)then
C***              zxtp10(jl,jk)=zxtp10(jl,jk)-zso2ev(jl,jk)
C***     &                              /(1.-pclcover(jl,jk)-pclcon(jl,jk))
C***              zdep3d(jl,jk)=zdep3d(jl,jk)+zso2ev(jl,jk)
C***            endif
C***          enddo
C***        enddo
C***      elseif(jt.eq.itracso2+1)then
C***        do jk=1,klev
C***          do jl=1,klon
C***            if(zso2ev(jl,jk).gt.zmin)then
C***              zxtp10(jl,jk)=zxtp10(jl,jk)+zso2ev(jl,jk)
C***     &                              /(1.-pclcover(jl,jk)-pclcon(jl,jk))
C***              zdep3d(jl,jk)=zdep3d(jl,jk)-zso2ev(jl,jk)
C***            endif
C***          enddo
C***        enddo
C***      endif

C
C   CALCULATE NEW CHEMISTRY AND SCAVENGING TENDENCIES

      DO 342 JK=KTOP,KLEV
       DO 340 JL=1,KLON
          ZXTP1=(1.-pclcover(jl,jk)-pclcon(jl,jk))*ZXTP10(JL,JK)+
     *           PCLCOVER(JL,JK)*ZXTP1C(JL,JK)
     *           + pclcon(jl,jk)*zxtp1con(jl,jk)
          ZDXTE(JL,JK,JT)=(ZXTP1-XTM1(JL,JK,JT))*PQTMST  !Total tendency (Dep + chem)
  340  CONTINUE
  342 CONTINUE

c Note that stratwd as coded here includes the below-cloud convective scavenging/evaporation

      if(jt.eq.itracso2)then
c        so2ao(:)=so2h2(:)+so2o3(:)
        stratwd(:)=0.
        do jk=1,klev
          do jl=1,klon
            stratwd(jl)=stratwd(jl)+zdep3d(jl,jk)*pdpp1(jl,jk)
     &           /(grav*ptmst)
          enddo
        enddo
        so2wd(:)=so2wd(:)+stratwd(:)
      elseif(jt.eq.itracso2+1)then
        stratwd(:)=0.
        do jk=1,klev
          do jl=1,klon
            stratwd(jl)=stratwd(jl)+zdep3d(jl,jk)*pdpp1(jl,jk)
     &           /(grav*ptmst)
          enddo
        enddo
        so4wd(:)=so4wd(:)+stratwd(:)
      endif

      if(debug.and.(jt.eq.itracso2.or.jt.eq.itracso2+1))then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          tdt=ptmst
          write(25,'(a,3i3)')'IPASS=1, after xtwetdep, JT= ',jt
          write(25,9)'conwd, stratwd, so2wd, so4wd ', conwd(mg,jt)*tdt,
     &         stratwd(mg)*tdt,so2wd(mg)*tdt,so4wd(mg)*tdt
          write(25,9)'zdep3d ',(zdep3d(mg,k),k=nl,1,-1)
          write(25,9)'zdxte ',(zdxte(mg,k,jt)*ptmst,k=nl,1,-1)
          write(25,9)'pmratep',(pmratep(mg,k),k=nl,1,-1)
          write(25,9)'pfprec ',(pfprec(mg,k),k=nl,1,-1)
          write(25,9)'pfevap ',(pfevap(mg,k),k=nl,1,-1)
          write(25,9)'pfconv ',(pfconv(mg,k),k=nl,1,-1)
          write(25,9)'pclcov ',(pclcover(mg,k),k=nl,1,-1)
          write(25,9)'pmaccr',(pmaccr(mg,k),k=nl,1,-1)
          write(25,9)'pqfsed',(pqfsed(mg,k),k=nl,1,-1)
          write(25,9)'pfsnow ',(pfsnow(mg,k),k=nl,1,-1)
          write(25,9)'pfmelt ',(pfmelt(mg,k),k=nl,1,-1)
          write(25,9)'pfsubl ',(pfsubl(mg,k),k=nl,1,-1)
          write(25,9)'pfstay ',(pfstay(mg,k),k=nl,1,-1)
          write(25,9)'pcfcov ',(pcfcover(mg,k),k=nl,1,-1)
          write(25,9)'pclcon ',(pclcon(mg,k),k=nl,1,-1)
          write(25,9)'pccw   ',(pccw(mg,k),k=nl,1,-1)
        endif
 9    format(a,30g10.3)
      endif


C
      ENDIF
  335 CONTINUE
C   END OF TRACER LOOP
C
C    CHANGE THE TOTAL TENDENCIES
C
C  ZDXTE(ITRACSO2) = TENDENCY OF SO2
C  ZDXTE(ITRACSO2+1) = CHEMISTRY AND SCAVENGING TENDENCY OF TOTAL
C       SULFATE
C
        JT=ITRACSO2
        DO 352 JK=KTOP,KLEV
          DO 350 JL=1,KLON
           xte(jl,jk,jt)=xte(jl,jk,jt)+zdxte(jl,jk,jt)
           XTE(JL,JK,JT+1)=XTE(JL,JK,JT+1)+ZDXTE(JL,JK,JT+1)
  350 CONTINUE
  352 CONTINUE
      
C BP commented out this portion for debugging Mk3L because ntrac=3 only
C      if(ntrac.gt.3)then
C        DO 362 JK=KTOP,KLEV
C          DO 360 JL=1,KLON
C            XTE(JL,JK,ITRACBC+1)=XTE(JL,JK,ITRACBC+1)+
C     *                          ZDXTE(JL,JK,ITRACBC+1)
C            XTE(JL,JK,ITRACOC+1)=XTE(JL,JK,ITRACOC+1)+
C     *                          ZDXTE(JL,JK,ITRACOC+1)
C 360      CONTINUE
C 362   CONTINUE
C      endif

C   OXIDATION WITH OH
      CONTINUE
C   CALCULATE THE DAY-LENGTH
c Need to hack this because of irritating CSIRO coding! (NH+SH latitudes concatenated!)
      ZDAYL=0.
      DO 402 JL=1,lon
        IF(ZRDAYL(JL).EQ.1) THEN
          ZDAYL=ZDAYL+1.
        ENDIF
  402 CONTINUE
      ZDAYFAC(:)=0.
      ZNLON=FLOAT(lon)
      IF(ZDAYL.NE.0.) ZDAYFAC(1)=ZNLON/ZDAYL !NH
      IF(ZDAYL.NE.znlon) ZDAYFAC(2)=ZNLON/(znlon-ZDAYL) !SH
C

        JT=ITRACSO2
C   DAY-TIME CHEMISTRY
      DO 412 JK=1,KLEV
       DO 410 JL=1,KLON
        IF(ZRDAYL(JL).EQ.1) THEN
         ins=(jl-1)/lon + 1 !hemisphere index
         X=PRHOP1(JL,JK)
         ZXTP1SO2=XTM1(JL,JK,JT)+XTE(JL,JK,JT)*PTMST
           ZTK2=ZK2*(PTP1(JL,JK)/300.)**(-3.3)
           ZM=X*ZNAMAIR
           ZHIL=ZTK2*ZM/ZK2I
           ZEXP=ALOG10(ZHIL)
           ZEXP=1./(1.+ZEXP*ZEXP)
           ZTK23B=ZTK2*ZM/(1.+ZHIL)*ZK2F**ZEXP
           ZSO2=ZXTP1SO2*ZZOH(JL,JK)*ZTK23B*ZDAYFAC(ins)
           ZSO2=AMIN1(ZSO2,ZXTP1SO2*PQTMST)
           ZSO2=AMAX1(ZSO2,0.)
           XTE(JL,JK,JT)=XTE(JL,JK,JT)-ZSO2
           XTE(JL,JK,JT+1)=XTE(JL,JK,JT+1)+ZSO2
           so2oh3d(jl,jk)=zso2
C
         ZXTP1DMS=XTM1(JL,JK,JT-1)+XTE(JL,JK,JT-1)*PTMST
         IF(ZXTP1DMS.LE.ZMIN) THEN
           ZDMS=0.
         ELSE
          T=PTP1(JL,JK)
c          ZTK1=(T*EXP(-234./T)+8.46E-10*EXP(7230./T)+
c     *         2.68E-10*EXP(7810./T))/(1.04E+11*T+88.1*EXP(7460./T)) !Original
          ztk1=1.646e-10-1.850e-12*t+8.151e-15*t**2-1.253e-17*t**3 !Cubic fit good enough
          ztk1=max(ztk1,5.e-12) !Because cubic falls away for T > 300 K
          ztk1=2.0*ztk1         !This is the fudge factor to account for other oxidants
         ZDMS=ZXTP1DMS*ZZOH(JL,JK)*ZTK1*ZDAYFAC(ins)
         ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
         XTE(JL,JK,JT-1)=XTE(JL,JK,JT-1)-ZDMS
         XTE(JL,JK,JT)=XTE(JL,JK,JT)+ZDMS
         dmsoh3d(jl,jk)=zdms
        ENDIF
       ENDIF
  410  CONTINUE
  412 CONTINUE
C
C   NIGHT-TIME CHEMISTRY
      DO 416 JK=1,KLEV
       DO 414 JL=1,KLON
        IF(ZRDAYL(JL).NE.1) THEN
         X=PRHOP1(JL,JK)
         ZXTP1DMS=XTM1(JL,JK,JT-1)+  XTE(JL,JK,JT-1)*PTMST
         IF(ZXTP1DMS.LE.ZMIN) THEN
           ZDMS=0.
         ELSE
           ZTK3=ZK3*EXP(520./PTP1(JL,JK))
C    CALCULATE THE STEADY STATE CONCENTRATION OF NO3
         ZQT=1./PTP1(JL,JK)
         ZQT3=300.*ZQT
         ZRHOAIR=PRHOP1(JL,JK)*ZNAMAIR
         ZKNO2O3=1.2E-13*EXP(-2450.*ZQT)
         ZKN2O5AQ=0.1E-04
         ZRX1=2.2E-30*ZQT3**3.9*ZRHOAIR
         ZRX2=1.5E-12*ZQT3**0.7
         ZKNO2NO3=ZRX1/(1.+ZRX1/ZRX2)*0.6**
     *            (1./(1.+(ALOG10(ZRX1/ZRX2))**2))
         ZEQN2O5=4.E-27*EXP(10930.*ZQT)
         ZKN2O5=ZKNO2NO3/ZEQN2O5
C
         ZNO3=ZKNO2O3*(ZKN2O5+ZKN2O5AQ)*ZZNO2(JL,JK)*ZZO3(JL,JK)
         ZZQ=ZKNO2NO3*ZKN2O5AQ*ZZNO2(JL,JK)+(ZKN2O5+ZKN2O5AQ)*ZTK3
     *      *ZXTP1DMS*XTOC(X,ZMOLGS)
         IF(ZZQ.GT.0.) THEN
             ZNO3=ZNO3/ZZQ
         ELSE
             ZNO3=0.
         ENDIF
         ZDMS=ZXTP1DMS*ZNO3*ZTK3
         ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
         XTE(JL,JK,JT-1)=XTE(JL,JK,JT-1)-ZDMS
         XTE(JL,JK,JT)=XTE(JL,JK,JT)+ZDMS
         dmsn33d(jl,jk)=zdms
       ENDIF
      ENDIF
  414 CONTINUE
  416 CONTINUE

c Calculate tendency of SO2 due to oxidation by OH (diagnostic) and ox. tendencies of DMS

      do jk=1,klev
        do jl=1,klon
          so2oh(jl)=so2oh(jl)+so2oh3d(jl,jk)*pdpp1(jl,jk)/grav
          dmsoh(jl)=dmsoh(jl)+dmsoh3d(jl,jk)*pdpp1(jl,jk)/grav
          dmsn3(jl)=dmsn3(jl)+dmsn33d(jl,jk)*pdpp1(jl,jk)/grav
        enddo
      enddo

      if(debug)then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, end xtchemie'
          write(25,9)'del(DMS)',(ptmst*xte(mg,k,itracso2-1),k=nl,1,-1)
          write(25,9)'del(SO2)',(ptmst*xte(mg,k,itracso2),k=nl,1,-1)
          write(25,9)'del(SO4)',(ptmst*xte(mg,k,itracso2+1),k=nl,1,-1)
        endif
      endif

      RETURN
      END

c******************************************************************************

      subroutine resetr(arr,len,val)
      
      real arr(len), val
      
      do i=1,len
        arr(i)=val
      enddo

      return
      end
