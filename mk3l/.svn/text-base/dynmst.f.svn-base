c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Extracted the COMMON block /GIANT4/ to a header file. The array PHNF
c is transferred to a new COMMON block /GIANT5/.
c SJP 2009/03/12
c
c nsib is modified to lsm_type = "nsib " comparison
c AJA 2009/01/22
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
c $Log: dynmst.f,v $
c Revision 1.35  2001/02/28 04:36:37  rot032
c Further tidy ups from HBG
c
c Revision 1.34  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.33  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.32  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.31  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.30  1998/12/10  00:55:37  ldr
c HBG changes to V5-1-21
c
c Revision 1.29  1997/12/23  00:23:35  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.28  1997/12/17  23:22:46  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.27.1.1  1997/12/19  02:03:12  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.27  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.26  1996/10/24  01:02:41  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.25  1996/06/13  02:06:24  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.24  1996/03/21  03:18:37  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.23  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.22  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.19.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.21  1995/06/30  02:44:40  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.20  1994/08/08  17:21:08  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.19  93/12/17  15:32:16  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.18  93/10/15  14:17:03  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c 
c Revision 1.17  93/10/07  12:09:39  ldr
c Move machine parameter statement into new include file MACHINE.f.
c 
c Revision 1.16  93/09/06  12:28:12  ldr
c Use correct formula for qs and dqsdt (i.e. don't subtract es from p).
c 
c Revision 1.15  93/08/19  15:07:54  ldr
c Minor cosmetic changes.
c 
c Revision 1.14  93/08/10  16:13:45  ldr
c Merge of HBG changes to V4-3-28l with other changes since then.
c 
c Revision 1.13  93/08/03  11:25:43  ldr
c Replace intrinsic functions with generic forms where possible.
c 
c Revision 1.12.1.1  93/08/10  15:27:05  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c
c     INPUT/OUTPUT:
c     Input:   from common/gausl in GAUSL.f
c                  sia - cosine latitude
c
c              from common/giant4 in this subroutine
c                  psl - sea level pressure
c                  rpn - pure moisture    tpn - pure temp(real temp)
c
c              from common/dystb in this subroutine
c                  profile arrays:
c                  qpa, qpb, qpc - moisture
c                  tpa, tpb, tpc - temperature
c                  wpa, wpb, wpc - winds
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from common/timex in TIMEX.f
c                  mins  - current model time in mins
c                  mstep - timestep in minutes 
c
c     In/Out:  from arguments
c                  ns - hemisphere index   lg - latitude index
c
c              from common/bgrnd in BGRND.f
c                  ibeta  - array of no of additions of above surface values
c                  month arrays of:
c                  pmonth -  sea level pressure    qmonth - moisture
c                  tmonth -  temperature           umonth - zonal wind
c                  vmonth -  meridional wind
c
c              from common/giant4 in this subroutine
c                  upn - u*cos(phi)/erad    vpn - v*cos(phi)/erad
c
c              from common/hybrpr in HYBRPR.f
c                  prf -pressure at full levels 
c
c              from common/masiv4 in MASIV4.f
c                  savegrid - global surface variables
c
c              from common/work1x in this subroutine
c                  pn - latitude row of surface pressure 
c
      subroutine dynmst(lg)

      implicit none

!$OMP THREADPRIVATE ( /GIANT4/ )
!$OMP THREADPRIVATE ( /GIANT5/ )
!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /WORK1X/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'GIANT4.f'
      include 'GIANT5.f'
      include 'HYBRPR.f'

      real un,vn,pn,dvn,von,ten,cpn,pln,rmn
      common/work1x/un(ln2,nl),vn(ln2,nl),pn(ln2),dvn(ln2,nl)
     & ,von(ln2,nl),ten(ln2,nl),cpn(ln2),pln(ln2),rmn(ln2,nl)

C Global data blocks
      include 'BGRND.f'
      include 'FEWFLAGS.f'
      include 'GAUSL.f'
      include 'LSMI.f'
      include 'MASIV4.f'
      include 'TIMEX.f'

      real tpa,tpb,tpc,wpa,wpb,wpc,qpa,qpb,qpc,plev,fred
      common/dystb/tpa(nl),tpb(nl),tpc(nl),wpa(nl),wpb(nl),wpc(nl)
     &,qpa(nl),qpb(nl),qpc(nl),plev(nl),fred(nl)

C Local work arrays and variables
      real at(ln2,0:nlp),bt(ln2,0:nlp),au(ln2,0:nlp),bu(ln2,0:nlp)
     &    ,av(ln2,0:nlp),bv(ln2,0:nlp),aq(ln2,0:nlp),bq(ln2,0:nlp)
     &    ,ar(ln2,0:nlp),br(ln2,0:nlp),ag(ln2,0:nlp),bg(ln2,0:nlp)
     &    ,rhn(ln2,nl),tg(ln2),qg(ln2)
      real psig(ln2,0:nl),splev(ln2,nl),aspl(ln2,nl),asp3(ln2,nl)
      real tdayex(ln2,nl),udayex(ln2,nl),vdayex(ln2,nl)
      real prfl(ln2,0:nl),prf3(ln2,0:nl)
      integer jindx(ln2,nl)
      logical test2(ln2,nl)
      logical test(lon,nl)  ! (lon is correct)
      logical tesx

      integer j
      integer k
      integer lga
      integer ma
      integer mg
      integer nbx
      integer ns

      real conr
      real gplk
      real pk
      real pz4m
      real qpak
      real qpbk
      real qpck
      real qplk
      real qs
      real qsg
      real rplk
      real siai
      real sigz4m
      real snowdep
      real temp0
      real tpak
      real tpbk
      real tpck
      real tplk
      real uplk
      real vplk
      real wetfac
      real wgmax
      real wpak
      real wpbk
      real wpck
      real z4max

C Local data, functions etc
      include 'ESTABL.f'  !Contains statement function form of establ

C Start code : ----------------------------------------------------------


c**** The pressure levels (pre-set) to which the model dynamics variables
c**** of temperature, winds, water vapour (T,U,V,Q) and the derived 
c**** relative humidity (RH) are interpolated are set up when the data
c**** files are initialized. This is done in ncinit.f (IF statsflag=T)
c**** (If statsflag=F, dummy plev values set in initax.f).

c**** The data levels MAY be given by 1000mb*(Sigma-level-value) of
c**** the model being used (usual strategy) 
c****    OR
c**** they may be pre-set to any sensible values BUT NOTE THAT
c**** when specifying arbitrary data levels (must have Nl of them)
c**** then the required pressure levels given by plev()
c**** must be in order and decreasing in pressure.
c     CALL TIMER('DYNMST  ',3)

C**** ROUTINE TO DERIVE DATA ON PRESET PRESSURE LEVELS
C**** (FOR ACCUMULATING DYNAMICAL STATISTICS OVER EG 1 MONTH).
C**** PRESET LEVELS ARE given by PLEV(Nl).
C****
C**** NOTE : THE CURRENT CODING IN THIS ROUTINE ONLY ALLOWS
C**** (FOR THE PRESSURE LEVELS SELECTED) FOR THE
C****  BOTTOM (NBELOW) LEVELS TO INTERSECT MOUNTAINS.
C****
C**** ACCUMULATING T,U,V,Q,RH,GPH(Nl LEVS),and PMSL
C****
C**** T ETC INTERPOLATED USING A*LOG(SIGMA)+B PROFILE
C****  (SURFACE TEMP T* IS USED).
C**** U,V ETC INTERPOLATED USING A*SIGMA+B PROFILE
C****  (SURFACE VELOCITY SET AT 0.7*BL VEL).
C**** Q ETC INTERPOLATED USING A*SIGMA**3+B PROFILE
C
C   ***************************************************
C   * NOTE **** NOTE **** NOTE **** NOTE **** NOTE    *
C   *    ONLY USE THIS ROUTINE AFTER ROUTINE DYNMNL   *
C   ***************************************************
C
C**** REMOVE COS(LAT)/ERAD WEIGHTING FROM UPN,VPN :
      siai=erad/sia(lg)
      do 10 k=1,nl
      do 10 mg=1,ln2
      upn(mg,k)=upn(mg,k)*siai
   10 vpn(mg,k)=vpn(mg,k)*siai
C**** SET UP SIGMA LEVELS ACCORDING TO REQUIRED P-LEVS.
C**** FOR THE REQUIRED DYNAMICAL STATS PRESSURE LEVELS,
C**** IT IS ASSUMED THERE WILL BE MOUNTAINS STICKING UP
C**** THROUGH LEVELS K = 1 - NBELOW
      do 12 k=1,nl
      do 12 mg=1,ln2
      if(hybrid)then
        splev(mg,k)=plev(k)
        prfl(mg,k)=log(prf(mg,k))
        prf3(mg,k)=prf(mg,k)**3
      else
        splev(mg,k)=plev(k)/pn(mg)
      endif
      aspl(mg,k)=log(splev(mg,k))
   12 asp3(mg,k)=splev(mg,k)**3

C**** CHECK IF THE SIGMA LEVEL HERE IS > 1 , THEN NBELOW IS 
C**** TOO SMALL. INCREASE (ELSE GARBAGE VALUES AT POINTS).
      if(hybrid)then
        do k=nbelow+1,nl
        do mg=1,ln2
          test2(mg,k)=splev(mg,k).gt.pn(mg)
        enddo
        enddo
      else
        do k=nbelow+1,nl
        do mg=1,ln2
          test2(mg,k)=splev(mg,k).gt.1.0
        enddo
        enddo
      endif
      call checkl(test2(1,nbelow+1),(nl-nbelow)*ln2,tesx)
      if(tesx)go to 1313

C**** SET UP SURFACE VALUES (FOR P-LEV BETWEEN LOWEST MODEL
C****  LEVEL AND SURFACE)
      wgmax=0.36
      if(lsm_type .eq. "nsib ")wgmax=0.42
      do 14 mg=1,ln2
      tg(mg)=savegrid(mg,3,lg)
      qsg=qsat(100.0*pn(mg),tg(mg))
      wetfac=1.0
      snowdep=savegrid(mg,12,lg)
      if((snowdep.le.0.0).and.(imsl(mg,lg).eq.4))
     &wetfac=savegrid(mg,8,lg)/wgmax
      qg(mg)=wetfac*qsg
   14 continue

C---- Create RH arrays
      do k=1,nl
      do mg=1,ln2
       pk=100.0*prf(mg,k)
       qs=qsat(pk,tpn(mg,k))
       rhn(mg,k)=100.0*max(0.0,min(1.0,rpn(mg,k)/qs))
      enddo
      enddo

C**** SET UP INTERPOLATER ARRAYS
      IF(hybrid)THEN

      do mg=1,ln2
        prfl(mg,0)=log(pn(mg))
        prf3(mg,0)=pn(mg)**3
      enddo

CXXXX BETWEEN ATMOS LEVELS 1:NL ( LEVEL INDEX 2:NL )
      do 15 k=2,nl
      do 15 mg=1,ln2
       tpak=prfl(mg,k-1)-prfl(mg,k)
       tpbk=prfl(mg,k-1)/tpak
       tpck=prfl(mg,k)/tpak
       at(mg,k)=(tpn(mg,k-1)-tpn(mg,k))/tpak
       bt(mg,k)=tpn(mg,k)*tpbk-tpn(mg,k-1)*tpck
       ag(mg,k)=(phnf(mg,k-1)-phnf(mg,k))/tpak
       bg(mg,k)=phnf(mg,k)*tpbk-phnf(mg,k-1)*tpck
       wpak=prf(mg,k-1)-prf(mg,k)
       wpbk=prf(mg,k-1)/wpak
       wpck=prf(mg,k)/wpak
       au(mg,k)=(upn(mg,k-1)-upn(mg,k))/wpak
       bu(mg,k)=upn(mg,k)*wpbk-upn(mg,k-1)*wpck
       av(mg,k)=(vpn(mg,k-1)-vpn(mg,k))/wpak
       bv(mg,k)=vpn(mg,k)*wpbk-vpn(mg,k-1)*wpck
       ar(mg,k)=(rhn(mg,k-1)-rhn(mg,k))/wpak
       br(mg,k)=rhn(mg,k)*wpbk-rhn(mg,k-1)*wpck
       qpak=prf3(mg,k-1)-prf3(mg,k)
       qpbk=prf3(mg,k-1)/qpak
       qpck=prf3(mg,k)/qpak
       aq(mg,k)=(rpn(mg,k-1)-rpn(mg,k))/qpak
       bq(mg,k)=rpn(mg,k)*qpbk-rpn(mg,k-1)*qpck
   15 continue
CXXXX BETWEEN SURFACE AND LEVEL 1 ( LEVEL INDEX 1 )
      do 16 mg=1,ln2
       tpak=prfl(mg,0)-prfl(mg,1)
       tpbk=prfl(mg,0)/tpak
       tpck=prfl(mg,1)/tpak
       at(mg,1)=(tg(mg)-tpn(mg,1))/tpak
       bt(mg,1)=tpn(mg,1)*tpbk-tg(mg)*tpck
       ag(mg,1)=(phisf(mg)-phnf(mg,1))/tpak
       bg(mg,1)=phnf(mg,1)*tpbk-phisf(mg)*tpck
       wpak=pn(mg)-prf(mg,1)
       wpbk=pn(mg)/wpak
       wpck=prf(mg,1)/wpak
       au(mg,1)=-0.3*upn(mg,1)/wpak
       bu(mg,1)=upn(mg,1)*(wpbk-0.7*wpck)
       av(mg,1)=-0.3*vpn(mg,1)/wpak
       bv(mg,1)=vpn(mg,1)*(wpbk-0.7*wpck)
       ar(mg,1)=0.0
       br(mg,1)=rhn(mg,1)
       qpak=prf3(mg,0)-prf3(mg,1)
       qpbk=prf3(mg,0)/qpak
       qpck=prf3(mg,1)/qpak
       aq(mg,1)=(qg(mg)-rpn(mg,1))/qpak
       bq(mg,1)=rpn(mg,1)*qpbk-qg(mg)*qpck
   16 continue

      ELSE ! (.not.hybrid)

CXXXX BETWEEN ATMOS LEVELS 1:NL ( LEVEL INDEX 2:NL )
      do 151 k=2,nl
      do 151 mg=1,ln2
       at(mg,k)=(tpn(mg,k-1)-tpn(mg,k))/tpa(k)
       bt(mg,k)=tpn(mg,k)*tpb(k)-tpn(mg,k-1)*tpc(k)
       ag(mg,k)=(phnf(mg,k-1)-phnf(mg,k))/tpa(k)
       bg(mg,k)=phnf(mg,k)*tpb(k)-phnf(mg,k-1)*tpc(k)
       au(mg,k)=(upn(mg,k-1)-upn(mg,k))/wpa(k)
       bu(mg,k)=upn(mg,k)*wpb(k)-upn(mg,k-1)*wpc(k)
       av(mg,k)=(vpn(mg,k-1)-vpn(mg,k))/wpa(k)
       bv(mg,k)=vpn(mg,k)*wpb(k)-vpn(mg,k-1)*wpc(k)
       ar(mg,k)=(rhn(mg,k-1)-rhn(mg,k))/wpa(k)
       br(mg,k)=rhn(mg,k)*wpb(k)-rhn(mg,k-1)*wpc(k)
       aq(mg,k)=(rpn(mg,k-1)-rpn(mg,k))/qpa(k)
       bq(mg,k)=rpn(mg,k)*qpb(k)-rpn(mg,k-1)*qpc(k)
  151 continue
CXXXX BETWEEN SURFACE AND LEVEL 1 ( LEVEL INDEX 1 )
      do 161 mg=1,ln2
       at(mg,1)=(tg(mg)-tpn(mg,1))/tpa(1)
       bt(mg,1)=tpn(mg,1)*tpb(1)-tg(mg)*tpc(1)
       ag(mg,1)=(phisf(mg)-phnf(mg,1))/tpa(1)
       bg(mg,1)=phnf(mg,1)*tpb(1)-phisf(mg)*tpc(1)
       au(mg,1)=-0.3*upn(mg,1)/wpa(1)
       bu(mg,1)=upn(mg,1)*(wpb(1)-0.7*wpc(1))
       av(mg,1)=-0.3*vpn(mg,1)/wpa(1)
       bv(mg,1)=vpn(mg,1)*(wpb(1)-0.7*wpc(1))
       ar(mg,1)=0.0
       br(mg,1)=rhn(mg,1)
       aq(mg,1)=(qg(mg)-rpn(mg,1))/qpa(1)
       bq(mg,1)=rpn(mg,1)*qpb(1)-qg(mg)*qpc(1)
  161 continue

      ENDIF

      do 17 mg=1,ln2
CXXXX For bottom points
      at(mg,0)=0.0
      bt(mg,0)=0.0
      ag(mg,0)=0.0
      bg(mg,0)=0.0
      au(mg,0)=0.0
      bu(mg,0)=0.0
      av(mg,0)=0.0
      bv(mg,0)=0.0
      ar(mg,0)=0.0
      br(mg,0)=0.0
      aq(mg,0)=0.0
      bq(mg,0)=0.0
CXXXX For top points
      at(mg,nlp)=0.0
      bt(mg,nlp)=tpn(mg,nl)
      ag(mg,nlp)=0.0
      bg(mg,nlp)=phnf(mg,nl)
      au(mg,nlp)=0.0
      bu(mg,nlp)=upn(mg,nl)
      av(mg,nlp)=0.0
      bv(mg,nlp)=vpn(mg,nl)
      ar(mg,nlp)=0.0
      br(mg,nlp)=rhn(mg,nl)
      aq(mg,nlp)=0.0
      bq(mg,nlp)=rpn(mg,nl)
   17 continue
CXXXX
C**** SET UP PRESSURE AT EACH MODEL SIGMA LEVEL
      do 18 k=1,nl
      do 18 mg=1,ln2
   18 psig(mg,k)=prf(mg,k)

C****
C**** LOOP FOR DETERMINING BETWEEN WHICH
C****  SIGMA LEVELS THE REQUIRED PRESSURE LEVELS LIE.
C****
c---- set level index array first

CSJP  Former machine dependence at this point
      do 76 mg=1,ln2
   76 psig(mg,0)=pn(mg)
      do 78 k=1,nbelow
      do 78 mg=1,ln2
   78 jindx(mg,k)=0
      do 80 k=1,nl
      do 80 j=0,nl
      do 80 mg=1,ln2
   80 if(psig(mg,j).gt.plev(k))jindx(mg,k)=j+1

      do 610 k=1,nl
      do 610 mg=1,ln2
      j=jindx(mg,k)
      tplk=bt(mg,j)+at(mg,j)*aspl(mg,k)
      uplk=bu(mg,j)+au(mg,j)*splev(mg,k)
      vplk=bv(mg,j)+av(mg,j)*splev(mg,k)
      rplk=br(mg,j)+ar(mg,j)*splev(mg,k)
      qplk=bq(mg,j)+aq(mg,j)*asp3(mg,k)
      gplk=bg(mg,j)+ag(mg,j)*aspl(mg,k)
c**** add derived values at p-lev "k" to monthly registers
      tmonth(mg,lg,k)=tmonth(mg,lg,k)+tplk
      umonth(mg,lg,k)=umonth(mg,lg,k)+uplk
      vmonth(mg,lg,k)=vmonth(mg,lg,k)+vplk
      rmonth(mg,lg,k)=rmonth(mg,lg,k)+rplk
      qmonth(mg,lg,k)=qmonth(mg,lg,k)+qplk
      gmonth(mg,lg,k)=gmonth(mg,lg,k)+gplk/grav
      udayex(mg,k)=uplk  
      vdayex(mg,k)=vplk  
      tdayex(mg,k)=tplk  
  610 continue

c---- check for peculiar values - model going unstable?
      do ns=1,2
      ma=(ns-1)*lon

      do 6002 k=1,nl
      do 6002 mg=1,lon
 6002 test(mg,k)=abs(vdayex(mg+ma,k)).gt.100.
      call checkl(test,lon*nl,tesx)
      if(tesx)then
        do 6010 k=1,nl
        do 6010 mg=1,lon
        if(test(mg,k))
     &   write(6,*) 'strange V velocity : ',vdayex(mg+ma,k),mg,ns,lg,k
 6010 continue
      end if

      enddo

      do ns=1,2
      ma=(ns-1)*lon

      do 6020 k=1,nl
      do 6020 mg=1,lon
      test(mg,k)=(abs(tdayex(mg+ma,k)).lt.100.)
     &            .and.(tdayex(mg+ma,k).ne.0.) 
 6020 continue
      call checkl(test,lon*nl,tesx)
      if(tesx)then
        do 6030 k=1,nl
        do 6030 mg=1,lon
        if(test(mg,k)) 
     &   write(6,*) 'strange Temperature : ',tdayex(mg+ma,k),mg,ns,lg,k
 6030 continue
      end if

      enddo
c----

c**** add up points above ground level (for first nbelow p-levs
c****  in this version - could need more/less for different p-levs)
      do 122 k=1,nbelow
      do 122 mg=1,ln2
  122 ibeta(mg,lg,k)=ibeta(mg,lg,k)+min(1,jindx(mg,k))

C*    COLLECT MEAN SEA LEVEL PRESSURE PSL
      do 70 mg=1,ln2
   70 pmonth(mg,lg)=pmonth(mg,lg)+psl(mg)
c     CALL TIMER('DYNMST  ',4)

c Various diagnostics

      if(uvtflag.and.mod(mins+int(mstep, 8),1440_8).eq.0_8)then
        do ns=1,2
          ma=(ns-1)*lon
          write(61) real(lg),real(ns)
     &             ,((udayex(mg+ma,k),mg=1,lon),k=1,nl)
          write(62) real(lg),real(ns)
     &             ,((vdayex(mg+ma,k),mg=1,lon),k=1,nl)
          write(63) real(lg),real(ns)
     &             ,((tdayex(mg+ma,k),mg=1,lon),k=1,nl)
        enddo
      endif

      return

 1313 print *,'nbelow is not suitable for the selection of'
      print *,'pressure levels for dynamical stats in dynmst.f'
      print *,'Please increase nbelow in PARAMS.f'
c.... compute approx sigma level of highest mountain on the globe
c.... (Himalayas)
      z4max=-9999.0
      do 500 lga=1,lat
      do 500 mg=1,ln2
  500 z4max=max(z4max,savegrid(mg,1,lga))
      print *,'max mountain height=',z4max
      temp0=289.0
      conr=grav/(0.0065*rdry)
      sigz4m=(1.0-0.0065*z4max/temp0)**conr
      pz4m=1013.25*sigz4m + 20.0
      print *,'typical pressure +20 at top of mountain is ',pz4m
c.... check how many of the preset pressure levels are likely to be
c.... below this level
      nbx=0
      do 510 k=1,nl
  510 if(pz4m.lt.plev(k))nbx=nbx+1
      print *,'nbelow should probably be set to ',nbx

      stop
      end
