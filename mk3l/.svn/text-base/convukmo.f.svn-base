c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c Changes to add machine type 'ALPH'.
c SJP 2001/11/22
c
c $Log: convukmo.f,v $
c Revision 1.16  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.15  2001/06/04 02:26:56  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.14  2001/03/07 04:28:59  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.13  2001/02/22 05:34:41  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.12  2000/11/14 03:11:38  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.11  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.10  2000/06/20 02:27:10  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.9  1999/11/09 03:24:01  rot032
c Allow different CCW in R21 and T63 models.
c
c Revision 1.8  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.7  1998/12/10  00:55:41  ldr
c HBG changes to V5-1-21
c
c Revision 1.6  1997/12/23  00:23:36  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.5  1997/12/17  23:22:52  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.4.1.1  1997/12/19  02:03:15  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.4  1997/10/03  06:28:46  ldr
c Corrections to HBG/MRD stuff for f77 compiler.
c
c Revision 1.3  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.2  1997/07/24  05:42:50  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.1  1997/06/13  06:02:04  ldr
c Initial revision
c
      subroutine convukmo(ipass,lg,tdt,pg,qsg,xso4,
     &                ttg,qtg,rhg,precc,preci,fscav,qlg,qfg,
     &                fmp,fmu,fmd,kta,kba,fluxc,rainx,qevc,CCW)
C
C     CSIRO interface to the UKMO set of convection routines
C
C
C  Tree structure for the UKMO convection
C  
C  The convection is driven through an interface routine convukmo
C   (convukmo.f) set up specifically to call the main UKMO 
C    routine CONVECU (convec2a.f)
C
C---------------------------------------------------------------------
C  Flowtrace Showing Caller/Callee Relationships
C---------------------------------------------------------------------
C  
C  convukmo (convukmo.f)
C     |
C     |
C*****************************************************************
CNOTE - To save on the number of subroutines, all of the routines below 
C     - are now kept in one file called "ukall.f"
C*****************************************************************
C     |
C     |
C  CONVECU -(convec2a.f)
C          |
C         QSATU    (qsat1a.f)
C         FLAG_WET (flagw1a.f)
C         LAYER_CN (laycn1a.f)
C         DQS_DTH  (dqsdth1a.f)
C         WHENIMD  (UKMO/CRAY based - very slow - do not use)
C         LIFT_PAR (lifpar1a.f)
C          |    |
C          |   QSATU    (qsat1a.f)
C          |   LATENT_H (latent1a.f)
C          |
C          |
C         CONVEC2  (conv21a.f)
C          |    |
C          |   PARCEL   (parcel1a.f)
C          |    |   |
C          |    |  WHENIMD
C          |    |  DETRAIN  (detrai1a.f)
C          |    |   |   |
C          |    |   |  THP_DET  (thpdet1a.f)
C          |    |   |  QSATU    (qsat1a.f)
C          |    |   |  THETAR   (thetar1a.f)
C          |    |   |   |   |
C          |    |   |   |  QSATU    (qsat1a.f)
C          |    |   |   |
C          |    |   |  DET_RATE (detrat1a.f)
C          |    |   |
C          |    |  TERM_CON (termco1a.f)
C          |    |  CLOUD_W  (cloudw1a.f)
C          |    |       |
C          |    |      CON_RAD  (conrad1a.f)
C          |    |
C          |   ENVIRON  (enviro1a.f)
C          |
C          |
C         DD_CALL  (ddcall2a.f)
C          |    |
C          |   WHENIMD
C          |   FLX_INIT (flxini2a.f)
C          |   LAYER_DD (layerd2a.f)
C          |   DD_INIT  (ddinit2a.f)
C          |    |   |
C          |    |   |
C          |    |  SATCAL   (satcal2a.f)
C          |    |       |
C          |    |      QSATU    (qsat1a.f)
C          |    |      DQS_DTH  (dqsdth1a.f)
C          |    |
C          |   DOWND    (downd2a.f)
C          |        |
C          |       WHENIMD
C          |       DDRAUGHT (ddraug2a.f) 
C          |        |   |
C          |        |  SATCAL   (satcal2a.f)
C          |        |   |   |
C          |        |   |  QSATU    (qsat1a.f)
C          |        |   |  DQS_DTH  (dqsdth1a.f)
C          |        |   |
C          |        |  CRS_FRZL (crsfrz2a.f)
C          |        |  QSATU    (qsat1a.f)
C          |        |  DEVAP    (devap2a.f)
C          |        |   |   |
C          |        |   |  EVP      (evp2a.f)
C          |        |   |
C          |        |  TERMDD   (termdd2a.f)
C          |        |  DD_ENV   (ddenv2a.f)
C          |        |
C          |       CHG_PHSE (chgphs2a.f)
C          |       PEVP_BCB (pevbcb2a.f)
C          |            |
C          |           QSATU    (qsat1a.f)
C          |           EVP      (evp2a.f)
C          |           SATCAL   (satcal2a.f)
C          |                |
C          |               QSATU    (qsat1a.f)
C          |               DQS_DTH  (dqsdth1a.f)
C          |
C          |
C         COR_ENGY (coreng1a.f)
C
C---------------------------------------------------------------------
C
      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'
      include 'FEWFLAGS.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'C_LHEAT.f'
CHBG  include 'C_G.f'
      include 'UKPARAMS.f'
      include 'ECPARM.f'

C Argument list
      integer ipass               ! leads indicator
      integer lg                  ! latitude row indicator
      real tdt                    ! timestep
      real pg(ln2)                ! IN surface pressure (mbs)
      real qsg(ln2,nl)            ! OUT Saturation moisture
      real xso4(ln2,nl)           ! IN so4 mixing ratio (=xtg(:,:,3)) [kg/kg]
      real ttg(ln2,nl)            ! INOUT  MODEL TEMPERATURE (K)
      real qtg(ln2,nl)            ! INOUT
      real rhg(ln2,nl)            ! OUT Relative Humidity
      real precc(ln2)
      real preci(ln2)
      real fscav(ln2,nl,ntrac)    !INOUT tracer scavenging fraction
      real qlg(ln2,nl)            ! Atmos liquid water
      real qfg(ln2,nl)            ! Atmos frozen water 
      real fmp(ln2,nl)
      real fmu(ln2,nl)            ! Out updraught flux (s**(-1)) +ve
      real fmd(ln2,nl)            ! Out downdraught flux (s**(-1)) +ve
      integer kta(ln2)
      integer kba(ln2)
      real fluxc(ln2,nl)
      real rainx(ln2)
      real qevc(ln2,nl) !For LDR cloud scheme
      real CCW(ln2,nl)            ! OUT CONVECTIVE CLOUD LIQUID WATER
                                  ! (KG/KG) ON MODEL LEVELS

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'LOGIMSL.f'         ! IN for land mask

C Global data blocks
      include 'CNSTA.f'
      include 'HYBARR.f'
      include 'PRINTT.f'
      include 'TIMEX.f'           ! IN for mins,mstep
      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel
      integer icln,ipcln
      real avcvrn,pavcvrn
      common/newcl/icln(ln2,2,lat),avcvrn(ln2,lat)
     &           ,ipcln(ln2,2,lat),pavcvrn(ln2,lat)
      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

C Local work arrays and variables
      REAL PSTAR(ln2)             ! IN SURFACE PRESSURE (PA)
      real deltm(ln2,nl)
C---------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT/LOCAL
C---------------------------------------------------------------------
C
      REAL EXNER(ln2,nl+1)        ! IN EXNER RATIO
C
      REAL AK(nl),                ! IN HYBRID CO-ORDINATE COEFFICIENTS
     &     BK(nl)                 !    DEFINE PRESSURE AT MID-POINT
                                  !    OF LAYER K
C
      REAL AKM12(nl+1),           ! IN HYBRID CO-ORDINATE COEFFICIENTS
     &     BKM12(nl+1)            !    TO DEFINE PRESSURE AT
                                  !    LEVEL K-1/2
C
      REAL DELAK(nl),             ! IN DIFFERENCE IN HYBRID CO-ORDINATE
     &     DELBK(nl)              !    COEFFICIENTS ACROSS LAYER K
C
C---------------------------------------------------------------------
C  VARIABLES WHICH ARE INPUT AND OUTPUT
C---------------------------------------------------------------------
C
      REAL TH(ln2,nl)             ! INOUT
                                  ! IN MODEL POTENTIAL TEMPERATURE (K)
                                  ! OUT MODEL POTENTIAL TEMPERATURE
                                  !     AFTER CONVECTION (K)
C
C
      REAL Q(ln2,nl)              ! INOUT
                                  ! IN MODEL MIXING RATIO (KG/KG)
                                  ! OUT MODEL MIXING RATIO AFTER
                                  !     AFTER CONVECTION (KG/KG)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL DTHBYDT(ln2,nl)        ! OUT INCREMENTS TO POTENTIAL
                                  !     TEMPERATURE DUE TO CONVECTION
                                  !     (K/S)
C
      REAL DQBYDT(ln2,nl)         ! OUT INCREMENTS TO MIXING RATIO
                                  !     DUE TO CONVECTION (KG/KG/S)
C
      REAL RAIN(ln2)              ! OUT SURFACE CONVECTIVE RAINFALL
                                  !     (KG/M**2/S)
C
      REAL SNOW(ln2)              ! OUT SURFACE CONVECTIVE SNOWFALL
                                  !     (KG/M**2/S)
C
      REAL PRECIP(ln2,nl)         ! AMOUNT OF PRECIPITATION
                                  ! FROM EACH LAYER (KG/M*:2/S)
C
      REAL CCA(ln2)               ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(ln2)           ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(ln2)           ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCLWP(ln2)             ! OUT CONDENSED WATER PATH (KG/M**2)
C
      REAL FLX(ln2,nl)            ! PARCEL MASSFLUX IN LAYER K (PA/S)
C
      REAL DDFLX(ln2,nl)          ! DOWNDRAUGHT MASSFLUX IN LAYER K (PA/S)
C
      real CCWI(ln2,nl)           ! OUT CONVECTIVE CLOUD ICE
                                  ! (KG/KG) ON MODEL LEVELS
C
      real CCWL(ln2,nl)           ! OUT CONVECTIVE CLOUD LIQUID
                                  ! (KG/KG) ON MODEL LEVELS
C
      real dryttg(ln2,nl)         ! Amount of dry adjustment
C
CD    real ttgold(ln2,nl)         ! For checking temp changes
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCAL
C----------------------------------------------------------------------
C
      integer mg                  ! EW point indicator
      integer k                   ! vertical level indicator
      integer ns                  ! Hemisphere indicator
      logical test
      real temxk,temxkp,temx,dta
      character*1 chcon
      real rnx,snx
CD    real sumt,sumq,conrev
CD    real sumccwl,sumccwi,sumdq,sumdrs
      logical firstcall
      data firstcall/.true./
      save firstcall
      save AK,BK,DELAK,DELBK,AKM12,BKM12

C Local data, functions etc
      real epsil
      parameter (epsil=0.622) ! Ratio molec wt of H2O vapour to dry air
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------


      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,i1)')'Before convukmo. IPASS = ',ipass
C***          write(25,99)'DMS ',(xtg(mg,k,1),k=1,nl)
C***          write(25,99)'SO2 ',(xtg(mg,k,2),k=1,nl)
C***          write(25,99)'SO4 ',(xtg(mg,k,3),k=1,nl)
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,99)'qtg ',(qtg(mg,k),k=1,nl)
          do k=1,nl
            qsg(mg,k)=qsat(100.0*prf(mg,k),ttg(mg,k))
            rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
          enddo
          write(25,91)'rhg ',(rhg(mg,k),k=1,nl)
          write(25,*)
        endif
      endif

C
C----------------------------------------------------------------------
C Before using th UKMO convection, ensure dry static stability
C----------------------------------------------------------------------
C
      do 10 k=2,nl
      do 10 mg=1,ln2
        deltm(mg,k)=-0.5*KAPPA*log(prf(mg,k)/prf(mg,k-1))
   10 continue

c---- Keep the initial temperatures so that the amount of dry adjustment
c---- (if any) can be computed. Dry instabilities are removed before
c---- any moist convection is computed to avoid numerical instablities.
c---- The original dry instabilities can then be replaced following moist
c---- convection adjustments if so desired.
      do k=1,nl
      do mg=1,ln2
        dryttg(mg,k)=ttg(mg,k)
      enddo
      enddo

c**** check for dry instability. (adjust with overkill).

CSJP  Former machine dependence at this point
         k=0
 6       k=k+1
         if (k.eq.nl) goto 9
 7       test=.true.
         do 8 mg=1,ln2
c**   dry adjustment of temps if (s(k)-s(k+1))/cp>0 .
c**   if adjustment occurs, then skip back down to
c**   level underneath to check if adjustment has
c**   created instability at that level.
            temxk=1.0-deltm(mg,k+1)
            temxkp=1.0+deltm(mg,k+1)
            temx=temxk*ttg(mg,k)-temxkp*ttg(mg,k+1)
            if (temx.le.0.0) go to 8
            test=.false.
            dta=0.1+temx*dprf(mg,k+1)/
     &      (dprf(mg,k+1)*temxk+dprf(mg,k)*temxkp)
            ttg(mg,k)=ttg(mg,k)-dta
            ttg(mg,k+1)=ttg(mg,k+1)+(dprf(mg,k)/dprf(mg,k+1))*dta
 8        continue
c--   check if any points adjusted.
c--   if so, and k.ne.1 , skip back 1 level.
         if (k.eq.1) go to 6
         if (test) go to 6
c--   adjustment has occured : skip back 1 level
         k=k-1
         go to 7

 9    continue 

      do k=1,nl
      do mg=1,ln2
c---- compute the temp adjustments due to dry instability removal
        dryttg(mg,k)=dryttg(mg,k)-ttg(mg,k)
      enddo
      enddo

C
C----------------------------------------------------------------------
C Set up Hybrid arrays to suit UKMO scheme (convert mbs to Pa)
C----------------------------------------------------------------------
C
      If(firstcall)Then
        do 20 k=1,nl
          AK(K)   = 100.0*anf(k)
          BK(K)   =       bnf(k)
          DELAK(K)=-100.0*dadnf(k)*dsk(k)
          DELBK(K)=      -dbdnf(k)*dsk(k)
   20   continue
        do 30 k=1,nl+1
          AKM12(K)= 100.0*anh(k)
          BKM12(K)=       bnh(k)
   30   continue
        firstcall=.false.
      Endif
C
C----------------------------------------------------------------------
C Set up surface pressure in Pa from mbs
C----------------------------------------------------------------------
C
      do 40 mg=1,ln2
        PSTAR(mg)=pg(mg)*100.0
   40 continue
C
C----------------------------------------------------------------------
C Set up potential temperature (K) and moisture
C Set mass fluxes to zero
C----------------------------------------------------------------------
C
      do 50 k=1,nl
      do 50 mg=1,ln2
CD      ttgold(mg,k)=ttg(mg,k)
        Q(mg,k)=qtg(mg,k)
        FLX(mg,k)=0.0
        DDFLX(mg,k)=0.0
        CCW(mg,k)=0.0
        CCWI(mg,k)=0.0
        CCWL(mg,k)=0.0
   50 continue
C
C----------------------------------------------------------------------
C Set up half level EXNER ratio : (p/P00)**KAPPA
C----------------------------------------------------------------------
C
      do 60 mg=1,ln2
        EXNER(mg,1)=(PSTAR(mg)/PREF)**KAPPA
   60 continue
      do 70 k=2,nl
      do 70 mg=1,ln2
        EXNER(mg,k)=(prh(mg,k)*100.0/PREF)**KAPPA
   70 continue
      do 80 mg=1,ln2
        if(prh(mg,nl+1).ne.0.0)then
          EXNER(mg,nl+1)=(prh(mg,nl+1)*100.0/PREF)**KAPPA
        else
          EXNER(mg,nl+1)=0.0
        endif
   80 continue
C
C----------------------------------------------------------------------
C  Reset convective cloud stats each day for pure restarts
C----------------------------------------------------------------------
C
      if((.not.SCM.and.(mod(mins,1440_8).eq.0_8).and.(ipass.eq.1))
     &     .or.(SCM.and.(nsteps.eq.0)))then
        do 90 mg=1,ln2
        icln(mg,1,lg)=0
        icln(mg,2,lg)=0
        avcvrn(mg,lg)=0.0
        ipcln(mg,1,lg)=0
        ipcln(mg,2,lg)=0
   90   pavcvrn(mg,lg)=0.0
      end if
C
C----------------------------------------------------------------------
C Call UKMO convection 
C----------------------------------------------------------------------
C
c     print *,'lg=',lg
      CALL CONVECU (qcloud,ipass,
     &              ttg,xso4,fscav,TH,Q,PSTAR,land,DTHBYDT,
     &              DQBYDT,RAIN,SNOW,CCA,ICCB,ICCT,CCLWP,EXNER,
     &              AK,BK,AKM12,BKM12,DELAK,DELBK,tdt,CCW,
     &              CCWI,CCWL,FLX,DDFLX,PRECIP)
c     if(lg.eq.lat)stop
C
C----------------------------------------------------------------------
C Diagnostics for checking energy conservation
C (Usually commented out) : CD
C----------------------------------------------------------------------
C
CD    print *,'Sum dT,dQ ',lg
CD    do 110 mg=1,ln2
CD    conrev=(100.0*pg(mg)/g/tdt)
CD    sumt=0.0
CD    sumdq=0.0
CD    sumccwi=0.0
CD    sumccwl=0.0
CD    do 100 k=1,nl
CD      sumt=sumt+(ttg(mg,k)-ttgold(mg,k))*dprf(mg,k)/pg(mg)
CD      sumccwi=sumccwi+CCWI(mg,k)*dprf(mg,k)/pg(mg)
CD      sumccwl=sumccwl+CCWL(mg,k)*dprf(mg,k)/pg(mg)
CD       sumdq=sumdq+(Q(mg,k)-qtg(mg,k))*dprf(mg,k)/pg(mg)
CD100 continue
CD       sumdrs=(RAIN(mg)+SNOW(mg))/conrev
CD   &          +sumccwl+sumccwi
CD    sumq=(LC*RAIN(mg)+(LC+LF)*SNOW(mg))/conrev/CP
CD    sumccwl=sumccwl*LC/CP
CD    sumccwi=sumccwi*(LC+LF)/CP
CD    sumq=sumq+sumccwl+sumccwi
CD    rnx=(0.5*tdt)*RAIN(mg)*(1440.0/mstep)
CD    snx=(0.5*tdt)*SNOW(mg)*(1440.0/mstep)
CD    if(sumq.gt.0.0)
CD   &write (6,105) mg,lg,sumq,sumt,(sumq/sumt),rnx,snx
CD   & ,sumccwl,sumccwi,sumdq,sumdrs
CD105 format (1x,2i3,' sumq,sumt ',2e13.6,2x,f6.3,2(1x,f6.2),
CD   & 4(1x,e13.6))
CD110 continue
CD    if(lg.eq.lat)stop
C
C----------------------------------------------------------------------
C Get back temperature (K) and moisture
C Set up mass flux fmp (sec**-1) for momentum mixing
C Get flux of rainfall from each layer fluxc (for qcloud) in KG/M**2
C----------------------------------------------------------------------
C
      do 120 k=1,nl
      do 120 mg=1,ln2
c---- Code to allow the vertical temp profile to retain the original 
c---- well mixed dry instabilities (i.e. replace dryttg)
        ttg(mg,k)=ttg(mg,k)+dryttg(mg,k)
        qtg(mg,k)=Q(mg,k)
        fmp(mg,k)=(FLX(mg,k)+DDFLX(mg,k))/(muf(mg,k)*100.0)
c        fmu(mg,k)=FLX(mg,k)/(muf(mg,k)*100.0)
        fmu(mg,k)=FLX(mg,k)/100. !In hPa/s
c        fmd(mg,k)=DDFLX(mg,k)/(muf(mg,k)*100.0)
        fmd(mg,k)=DDFLX(mg,k)/100. !In hPa/s
        fluxc(mg,k)=PRECIP(mg,k)*tdt
c       fluxc(mg,k)=0.0
        qevc(mg,k)=0.0
        qlg(mg,k)=qlg(mg,k)+CCWL(mg,k) ! Ejected from cloud top
        qfg(mg,k)=qfg(mg,k)+CCWI(mg,k) ! Ejected from cloud top
C***        if(.not.land(mg).and.lg.ge.20.and.mg.ge.25.and.mg.le.34)then
C***        if(ccwi(mg,k).gt.0)write(26,'(2g13.4)')
C***     $         fmu(mg,iccb(mg)),1.e3*ccwi(mg,k)
C***        endif
C***c        if(k.eq.icct(mg))write(27,'(2g13.3)')ttg(mg,1),1.e3*ccwi(mg,k)
  120 continue
C
C----------------------------------------------------------------------
C  RAIN & SNOW in (KG/M**2/S)=mms/sec
C  Convective rainfall : precc in mms/step,   rainx in mms/day
C  Convective snowfall : preci in mms/step
C----------------------------------------------------------------------
C
      do 130 mg=1,ln2
        rnx=(0.5*tdt)*RAIN(mg)
        snx=(0.5*tdt)*SNOW(mg)
        rainx(mg)=(rnx+snx)*(1440.0/mstep)
        preci(mg)=preci(mg)+snx
        precc(mg)=precc(mg)+(rnx+snx)
  130 continue

C
C----------------------------------------------------------------------
C           generation of convective cloud
C
C.... save the max kt and min kb over the radiation step (2 Hrs)
C.... for the Slingo 87 cloud formulation.
C.... icln will be reset to zero following the call to CLOUD
C----------------------------------------------------------------------
C
      do 140 mg=1,ln2
        kba(mg)=ICCB(mg)
        kta(mg)=min(nl-1,ICCT(mg))
  140 continue

c The scaling down of CCW has been moved to progcld (LDR, 3/01)
c Wish to retain updraft (larger) values for use in xtchemie.

        do k=1,nl
          do mg=1,ln2
           if(k.lt.kba(mg))CCW(mg,k)=0.0
          enddo
        enddo


      if(ipass.eq.1)then
        do 150 mg=1,ln2
         if(icln(mg,1,lg).eq.0)icln(mg,1,lg)=kba(mg)
         icln(mg,1,lg)=min(icln(mg,1,lg),kba(mg))
         icln(mg,2,lg)=max(icln(mg,2,lg),kta(mg))
  150   continue
      elseif(ipass.eq.2)then
        do 160 mg=1,ln2
         if(ipcln(mg,1,lg).eq.0)ipcln(mg,1,lg)=kba(mg)
         ipcln(mg,1,lg)=min(ipcln(mg,1,lg),kba(mg))
         ipcln(mg,2,lg)=max(ipcln(mg,2,lg),kta(mg))
  160   continue
      endif
C
C----------------------------------------------------------------------
C.... add convective cloud amount : use the avcvrn arrays for this
C.... => divide by number of steps when using
C----------------------------------------------------------------------
C
      do 170 mg=1,ln2
        if(qcloud)rainx(mg)=CCA(mg)
        avcvrn(mg,lg)=avcvrn(mg,lg) + CCA(mg)*(2.0-ipass)
        pavcvrn(mg,lg)=pavcvrn(mg,lg) + CCA(mg)*(ipass-1.0)
  170 continue
C
C----------------------------------------------------------------------
C.... set up moisture values following convection
C----------------------------------------------------------------------
C
      do 175 k=1,nl
        do 175 mg=1,ln2
          qsg(mg,k)=qsat(100.0*prf(mg,k),ttg(mg,k))
          rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
  175 continue

C
C----------------------------------------------------------------------
C.... Update convection map at end of day only
C.... (see initax.f for klcmc)
C----------------------------------------------------------------------
C
      if ((mod(mins+int(mstep, 8),1440_8).eq.0_8).and.cvrnm) then
      do 180 mg=1,ln2
         chcon=' '
         if (kta(mg).eq.0) go to 180
         chcon='L'
         if (kta(mg).le.klcmc) go to 180
         chcon='M'
         if (kba(mg).ge.klcmc) go to 180
         chcon='P'
 180     cnmlp(mg,lg)=chcon
      endif

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(3(a,i2))')'After convukmo. IPASS = ',ipass,' kba ',
     &         kba(mg),' kta ',kta(mg)
C***          write(25,99)'DMS ',(xtg(mg,k,1),k=1,nl)
C***          write(25,99)'SO2 ',(xtg(mg,k,2),k=1,nl)
C***          write(25,99)'SO4 ',(xtg(mg,k,3),k=1,nl)
          write(25,99)'fmp ',(fmp(mg,k),k=1,nl)
          write(25,99)'fmu ',(fmu(mg,k),k=1,nl)
          write(25,99)'fmd ',(fmd(mg,k),k=1,nl)
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,99)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,91)'rhg ',(rhg(mg,k),k=1,nl)
          write(25,99)'ccw ',(ccw(mg,k),k=1,nl)
          write(25,1)'precc ',precc(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,g10.3))
 91   format(a,30f10.3)
 99   format(a,30g10.3)

      return
      end
