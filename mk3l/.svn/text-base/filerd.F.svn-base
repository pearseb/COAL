c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c $Log: filerd.f,v $
c Revision 1.98  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.97  2001/02/28 04:36:37  rot032
c Further tidy ups from HBG
c
c Revision 1.96  2001/02/22 05:34:38  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.95  2001/02/12 05:39:46  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.94  2000/12/08 03:58:52  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.93  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.92  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.91  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.90  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.89.1.1  2000/06/20 02:27:09  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.89  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.88  1998/12/10  01:08:39  ldr
c HBG changes to V5-1-21
c
c Revision 1.87  1998/05/27  02:07:37  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.86  1997/12/23  00:23:34  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.85  1997/12/17  23:22:45  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.84  1997/12/08  06:00:40  ldr
c Correct order of reading sdot, ochfa (needed for qflux runs).
c
c Revision 1.83.2.1  1997/12/19  02:03:11  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.83  1997/10/06  04:15:03  ldr
c Final corrections to V5-1.
c
c Revision 1.82  1997/07/24  05:42:49  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.81  1996/12/19  23:40:50  ldr
c EAK's snow albedo changes merged and debugged.
c
c Revision 1.80  1996/11/22  04:54:06  ldr
c Enable HBG's hybrid stuff in filerd9r also!
c
c Revision 1.79  1996/11/19  06:05:53  ldr
c Enable HBG's hybrid stuff.
c
c Revision 1.78  1996/10/24  01:02:43  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.77  1996/03/21  03:41:39  ldr
c Merge of TIE and LDR changes.
c
c Revision 1.76  1996/03/21  03:18:39  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.75.1.1  1996/03/21  03:29:07  ldr
c A couple of housekeeping tidy-ups from LDR.
c
c Revision 1.75  1996/02/19  04:09:47  ldr
c Generalize for 24 levels.
c
c Revision 1.74  1996/02/08  00:21:01  ldr
c Tidy up reading of qcloud stuff in filerd and remove obsolete slowice flag.
c
c Revision 1.73  1995/11/23  06:09:40  ldr
c Merge of LDR and EAK changes since V4-7-21l.
c
c Revision 1.72  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.71.1.1  1995/11/23  05:52:11  ldr
c Fix for T63 header and new nsstop=-2 option.
c
c Revision 1.71  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.70  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.63.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.69.1.1  1995/11/10  05:23:03  ldr
c EAK's changes to soil types (increase to 9 soil types).
c
c Revision 1.69  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.68  1995/02/27  01:40:34  ldr
c Bugfixes received from HBG in response to problems on Fujitsu.
c
c Revision 1.67  1994/12/22  23:51:06  ldr
c Mods to cloud scheme of V4-6-13l to produce faster version for HBG and
c to allow vertically subgrid cloud below inversions in BL. Set slowice=T
c for accurate (but slow) ice fall calculation. Also, Ri now includes hlfusionn
c
c Revision 1.66  94/08/08  17:21:12  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.65  94/08/08  13:11:17  ldr
c Implement memory-saving trick: NX for tracers and NC for qcloud.
c 
c Revision 1.64  94/08/04  10:35:19  ldr
c Add cloud fraction to RS file - only needed so that t-1 value can be 
c used in estimate of Ri.
c 
c Revision 1.63  94/05/13  14:55:26  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c 
c Revision 1.62  94/03/30  12:34:17  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.61  94/01/17  11:28:27  ldr
c Changes from Hal interpolate qfluxes and make them implicit. Also move
c read of qfluxes to datard so it is done every month.
c 
c Revision 1.60  94/01/04  17:18:05  ldr
c A little fix so that new qflux runs are initialized properly.
c 
c Revision 1.59  93/12/17  15:32:27  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.58  93/12/10  16:38:29  ldr
c New albedo and SiB data sets to go with new mask in RS file.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  hybrid - if true, use new hybrid vertical coordinate
c                  ifwds  - default F, only T when doing forward step from
c                           new restart file
c                  lcouple - flag to run model in coupled mode, default F
c                  nsstop  - no of time steps for model run
c                  qcloud  - T if using prognostic cloud scheme
c                  qflux   - set to T for qflux run, default F
c                  semice  - if T, use the Semtner sea-ice model
c                  sltrace - if T compute tracer transport by semi-
c                            Langrangian advection and by mixing and
c                            convection in the vertical
c
c              from common/files in FILES.f
c                  irfilename - name of input restart file
c
c              from common/timex in TIMEX.f
c                  mstep - timestep in minutes
c
c     Output:  from common/fldmri in FLDMRI.f
c                  fldm4x - array of zeroes to pad out common block for 
c                         - vectorization
c
c              from common/fldri in FLDRI.f
c                  fld4x - as for fldm4x, above
c
c              from common/freeze in FREEZE.f
c                  tfi - temperature at bottom of ice
c
c              from common/masiv1 in this subroutine
c                  ronmx- global array of measure of east-west surface stress
c                  sonmx- global array of measure of north-sth surface stress
c                  iphys ) - current phys time step and previous phys time
c                  iphysm)   step, these indices swap around
c
c              from common/masiv3 in MASIV3.f
c                  rgsav - incoming long-wave radiation
c                  sgamp - amplitude of net solar radiation at the ground
c                  sgsav - net solar radiation at ground
c
c              from common/masiv4 in MASIV4.f
c                  savegrid - global surface variables
c
c              from common/masiv5 in this subroutine
c                  statsice - global statistics relating to sea-ice scheme
c
c              from common/pmasiv3 in PMASIV3.f
c                  prgsav - rgsav over leads(see above)
c                  psgamp - sgamp over leads  "    "
c                  psgsav - sgsav over leads  "    "
c
c              from common/qcloud1 in QCLOUD1.f
c                  qfb  - cloud ice at t     qlb  - cloud liquid water at t
c                  qfbm - cloud ice at t-1   qlbm - cloud liquid water at t-1
c                  cfb  - cloud fraction
c
c              from common/qflxdat in QFLXDAT.f
c                  newqfluxrun - if T, use weight of zero at start of new 
c                                qflux run
c                  ochfa - ocean heat flux average
c
c              from common/positn in this subroutine
c                  latpos - latitude level of tracer particles
c                  lonpos - longitude level of tracer particles
c                  npts   - no. tracer particles
c                  sigpos - sigma level of tracer particles
c
c              from common/rmgrid in RMGRID.f
c                  rgt  - moisture tendency
c                  rmg  - pressure weighted moistures at current timestep
c                  rmmg - pressure weighted mositures at previous timestep
c
c              from common/surf in this subroutine
c                  aftfh  - reciprocal of air resistance for ground 
c                           covered with foliage
c                  aftfhg - reciprocal of air resistance for bare ground
c                  ocondx - precipitation at previous timestep
c                  opreci - frozen precipitation at previous timestep
c
c              from common/surf1 in this subroutine
c                  mc   - moisture depth on canopy
c                  tgf  - canopy temperature
c                  tgg  - bare ground temperature
c
c              from common/traceblk in TRACEBLK.f
c                  con  - concentration of atmospheric tracers at t
c                  conm - concentration of atmospheric tracers at t-1
c
c              from common/uvpgd
c                  sdot - d(sigma)/d(t)
c
c     In/Out:  from common/fldmri in FLDMRI.f
c                  pmi  - spectral surface pressure, imaginary part
c                  pmr  - spectral surface pressure, real part
c                  psimi - stream function, imaginary part
c                  psimr - stream function, real part 
c                  temi - spectral temp, imaginary part at t-1
c                  temr - spectral temp, real part at t-1
c                  xhimi - velocity potential at t-1, imaginary part
c                  xhimr - velocity potential at t-1, real part
c
c              from common/fldri in FLDRI.f
c                  pi - spectral pressure, imaginary part
c                  pr - spectral pressure, real part
c                  psii - spectral stream function, imaginary part
c                  psir - spectral stream function, real part
c                  tei - spectral temp, imaginary part
c                  ter - spectral temp, real part
c                  xhii - velocity potental, imaginary part
c                  xhir - velocity potental, real part
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from common/mountn in this subroutine
c                  phisi- spectral surface geopotential height,imaginary part
c                  phisr- spectral surface geopotential height, real part
c
c              from common/rshead in this subroutine
c                  header2 - second header on restart file
c
c              from common/timex in TIMEX.f
c                  ndays - day counter for model run
c
c              from argument
c                  exptyp - restart file header
c
c 
      subroutine filerd(exptyp)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f'
      include 'VERSION.f'
      integer nps
      parameter (nps=63*48)

C Argument list
      character*50 exptyp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'DIFM.f'
      include 'ECTRAC.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'FLDMRI.f'
      include 'FLDRI.f'
      include 'FREEZE.f'
      include 'LSMI.f'
      include 'MASIV3.f'
      include 'MASIV4.f'
      include 'PMASIV3.f'
      include 'QCLOUD1.f'
      include 'QFLXDAT.f'
      include 'RMGRID.f'
      include 'SURF1.f'
      include 'TIMEX.f'
      include 'TRACEBLK.f'
      include 'UVPGD.f'

      real ronmx,sonmx
       integer iphysm,iphys
      common/masiv1/ronmx(ln2,nl,lat,2),sonmx(ln2,nl,lat,2)
     &             ,iphysm,iphys

      real statsice
      common/masiv5/statsice(ln2,14,lat)

      real phisr,phisi
      common/mountn/phisr(lw,mw),phisi(lw,mw)

      real lonpos,latpos,sigpos
      integer npts
      common/positn/lonpos(nps),latpos(nps),sigpos(nps),npts

      character*50 header2
      common/rshead/header2

      real resvr1,resrem,srunoff,sresvr1
      common/rivers/resvr1(0:lon+1,0:lat2),resrem(lon,lat2)
     & ,srunoff(lon,lat2),sresvr1(lon,lat2)

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      real griddata(ln2,ngrid)
      real ron(ln2,nl),son(ln2,nl)
      real qin(ln2,lat)
      real stice(ln2,9)

      complex spec(nw)

      integer ierr
      integer ifiler
      integer k
      integer lg
      integer lwf
      integer mstepx
      integer mg
      integer ndaysf
      integer ntracef
      integer ntracers

      character char

C Local data, functions etc

C Start code : ----------------------------------------------------------
C
C FUNCTION : TO INITIALISE Mk3 MODEL FROM RESTART FILE.
C

c----  OPEN THE ATMOSPHERIC RESTART FILE
      ifiler=8  
      open(unit=ifiler,file=irfilename,form='unformatted',status='old'
     &     ,iostat=ierr)
      call errcheck(ierr,'restart file     ','filerd    ')
c
C**** READ HEADER FROM RESTART FILE
      read(ifiler)exptyp

      if(exptyp(1:3).ne.'Mk3')then
        write(6,*)
     &   'Error: Not a Mk3 restart file. Use Mk3file.f to convert'
        write(0,*)
     &   'Error: Not a Mk3 restart file. Use Mk3file.f to convert'
        stop
      endif

C---- Special fix for HBG T63 restart file
C____ Remove redundant bit in header 
      if(exptyp(12:15).eq.'MA97')exptyp(12:15)='    '
C---- 

      read(ifiler)ndays,mstepx

      if(.not.ncepagcm)then
      if(ifwds.or.nsstop.eq.-2)then
        ndays=0
        write(6,*)' Setting NDAYS to 0 in filerd (New model)'
      endif
      endif
      write(6,601)ndays,mstepx
  601 format(1h ,'ndays=',i9,' mstep=',i2)
      if(mstep.ne.mstepx)
     &  write(6,602)mstepx,mstep
  602   format(1x,'change of step from ',i2,' to ',i2,' mins')
      write(6,93)exptyp
  93  format(2x,'Restart file type is : ',a50)

c Second header with additional information

      header2='                                                  '
      read(ifiler)header2
      write(6,931)header2
 931  format(' Header2 on RS file is : ',a50)

#ifdef ALPH
      header2(1:4)="ALPH"
#else
      header2(1:4)="LINU"
#endif

      header2(15:24)=version

c**** Check for matching hybrid/sigma vertical coordinate & restart file
      if(hybrid.and.header2(32:36).ne.'HYBRD')then
        print*,'Warning: RS file not hybrid, but hybrid=T.'
        write(0,*)'Warning: RS file not hybrid, but hybrid=T.'
        header2(32:36)='HYBRD'
      elseif(.not.hybrid.and.header2(32:36).eq.'HYBRD')then
        print*,'Warning: RS file is hybrid, but hybrid=F.'
        write(0,*)'Warning: RS file is hybrid, but hybrid=F.'
        header2(32:36)='     '
      endif

c**** read spectral geopotential
      read (ifiler) spec
      call setr (spec,
     &           phisr,phisi)

C**** READ IN PSI(flux) (T then T-1)
      do 10 k=1,nl
         read (ifiler) spec
 10      call setr (spec,
     &              psir(1,1,k),psii(1,1,k))

      do 20 k=1,nl
         read (ifiler) spec
 20      call setr (spec,
     &              psimr(1,1,k),psimi(1,1,k))

C**** READ IN XHI(flux) (T then T-1)
      do 30 k=1,nl
         read (ifiler) spec
 30      call setr (spec,
     &                  xhir(1,1,k),xhii(1,1,k))
      do 40 k=1,nl
         read (ifiler) spec
 40      call setr (spec,
     &              xhimr(1,1,k),xhimi(1,1,k))

C**** READ IN TEMPS (P*T') (T then T-1)
      do 50 k=1,nl
         read (ifiler) spec
 50      call setr (spec,
     &              ter(1,1,k),tei(1,1,k))
      do 60 k=1,nl
         read (ifiler) spec
 60      call setr (spec,
     &              temr(1,1,k),temi(1,1,k))

C**** READ IN P* (T then T-1)
      read (ifiler) spec
      call setr (spec,
     &           prr,pri)
      read (ifiler) spec
      call setr (spec,
     &           pmr,pmi)

C**** READ IN GRID MOISTURE FIELDS PLUS LAST TENDENCY
      do 70 k=1,nl
      read (ifiler) qin
      do 71 lg=1,lat
      do 71 mg=1,ln2
   71 rmg(mg,k,lg)=qin(mg,lg)
   70 continue
      do 72 k=1,nl
      read (ifiler) qin
      do 73 lg=1,lat
      do 73 mg=1,ln2
   73 rmmg(mg,k,lg)=qin(mg,lg)
   72 continue
      do 74 k=1,nl
      read (ifiler) qin
      do 75 lg=1,lat
      do 75 mg=1,ln2
   75 rgt(mg,k,lg)=qin(mg,lg)
   74 continue

C**** GET SURFACE TYPE INDICATOR ARRAY
      read (ifiler) imsl

C**** SET UP SURFACE STRESS FILE IDENTIFIERS
      iphysm=1 ! Tau-1
      iphys=2  ! Tau
C****
C****     PUT GRID POINT DATA INTO STORAGE ARRAYS
C****
      do 90 lg= 1,lat
      read (ifiler) griddata
      do 91 k=1,ngrid
      do 91 mg=1,ln2
   91 savegrid(mg,k,lg)=griddata(mg,k)
      read(ifiler)ron,son
      do 94 k=1,nl
      do 94 mg=1,ln2
      ronmx(mg,k,lg,iphysm)=ron(mg,k)
   94 sonmx(mg,k,lg,iphysm)=son(mg,k)
 90   continue


c Read sea ice data
      
      tfi=271.3  !Temperature at bottom of ice
      if(exptyp(5:10).eq.'SEMICE')then !i.e. if restart file has ice data
        do 96 lg= 1,lat
        read (ifiler) stice
        do 95 k=1,9
        do 95 mg=1,ln2
 95     statsice(mg,k,lg)=stice(mg,k)
 96     continue
        if(.not.semice) then
          exptyp(5:10)='      '
          print*,'Warning: SEMICE false but RS file has SEMICE data'
        endif
      elseif(semice)then        !if restart file has no sea-ice model data
        print*,' Warning: Restart file has no sea-ice model data'
        call initice
        exptyp(5:10)='SEMICE'
      endif                     ! exptyp
        

c**** read dpsi,dxhi (spectral diffusion energy
c****    loss conversion to temp tendency terms)
c**** Split vertical levels now used for all resolutions
      
      do k=1,nl
        read (ifiler) spec
        call setr (spec,
     &             dpsir(1,1,k),dpsii(1,1,k))
      enddo
      do k=1,nl
        read (ifiler) spec
        call setr (spec,
     &             dxhir(1,1,k),dxhii(1,1,k))
      enddo

c Read in vertical velocity sdot, in case it is needed

      read (ifiler) sdot


c     QFLUX section
c     The running mean of the surface forcing
      newqfluxrun=.false.
      if ( header2(26:30).eq.'OCHFA' ) then
        read(ifiler) ochfa
        print*,' Read running mean of ocean forcing OCHFA for QFLUX run'
        if (.not. qflux ) then
          print*,'WARNING: RS file contains OCHFA data but qflux=F'
          header2(26:30)='     '
        endif
      elseif ( qflux ) then
        print*, ' WARNING - Running mean of ocean forcing ochfa is not'
        print*, ' present in restart file. It is initialised to zero.'
        header2(26:30)='OCHFA'
        newqfluxrun=.true.
          do lg=1,lat
            do mg=1,ln2
              ochfa(mg,lg) = 0.
          end do
        end do
      endif


c NSIB section
C***      if(nsib)then

c.... Mk3 restart data for sib
      SELECT CASE(lsm_type)
      CASE("nsib ")
        read (ifiler) wb
        read (ifiler) wbice
        read (ifiler) tggsl
        read (ifiler) tggsn
        read (ifiler) ssdnn
        read (ifiler) ssdn3
        read (ifiler) smass
        read (ifiler) gflux
        read (ifiler) sgflux
        read (ifiler) isflag
        read (ifiler) snage
        read (ifiler) osnowd
        read (ifiler) tgf
        read (ifiler) mc
        read (ifiler) aftfh
        read (ifiler) aftfhg 
        read (ifiler) ocondx
        read (ifiler) pmc
      CASE("cable")
        read (ifiler) wb
        read (ifiler) wbice
        read (ifiler) tggsl
        read (ifiler) tggsn
        read (ifiler) ssdnn
        read (ifiler) ssdn3
        read (ifiler) smass
        read (ifiler) gflux
        read (ifiler) sgflux
        read (ifiler) isflag
        read (ifiler) snage
        read (ifiler) osnowd
        read (ifiler) tgf
        read (ifiler) mc
        read (ifiler) aftfh
        read (ifiler) aftfhg 
        read (ifiler) ocondx
        read (ifiler) pmc
      CASE DEFAULT ! old land surface scheme removed
        PRINT *,'***  Invalid LAND SCHEME: Program Exiting'
        STOP
      END SELECT

      
c       do lg=1,lat
c       do mg=1,ln2
c        pmc(mg,lg)=0.0
c       enddo
c       enddo

c Set non-land values to missing value

      do 299 lg=1,lat
          do 297 mg=1,ln2
            if(imsl(mg,lg).ne.4)then
              tggsl(mg,1,lg)=-7777777
              tgf(mg,lg)=-7777777
              mc(mg,lg)=0.0
            endif
 297      continue
 299  continue

C***      endif
c End of NSIB section

c Read in radiation stuff

      read(ifiler)sgsav
      read(ifiler)rgsav
      read(ifiler)sgamp
      read(ifiler)psgsav
      read(ifiler)prgsav
      read(ifiler)psgamp

c******************************************************************************

c Cloud water and cloud ice

      if(qcloud.and.nc.eq.0)then
        print*,'ERROR: Need to set nc=1 in QCLOUD1.f'
        stop
      endif

c First do the header consistency check

      if(exptyp(47:49).eq.'QCL'.and..not.qcloud)then !RS file has qcloud data
        print*,'Warning: RS file has qcloud data but qcloud=F'
      elseif(exptyp(47:49).eq.'   '.and.qcloud)then !No data in RS file
        print*,'WARNING: Setting qcloud stuff to 0 in filerd'
        do lg=1,lat
          do k=1,nl
            do mg=1,ln2
              qlb(mg,k,lg)=0.
              qlbm(mg,k,lg)=0.
              qfb(mg,k,lg)=0.
              qfbm(mg,k,lg)=0.
              cfb(mg,k,lg)=0.
            enddo
          enddo
        enddo
        do lg=1,lat
            do mg=1,ln2
              opreci(mg,lg)=0.
            enddo
        enddo
      endif

c Now read the RS data, using the header as a guide

      if(exptyp(47:50).eq.'QCL4')then
        read(ifiler)qlb
        read(ifiler)qlbm
        read(ifiler)qfb
        read(ifiler)qfbm
        read(ifiler)cfb
        read(ifiler)opreci
      elseif(exptyp(47:49).eq.'QCL')then
        print*,'Sorry, these RS files are no longer supported'
        stop
      endif

      if(qcloud)then
        exptyp(47:50)='QCL4'
      else
        exptyp(47:50)='    '
      endif

c******************************************************************************

c Aerosol stuff

      if(coupled_aero.and.ne.eq.0)then
        print*,'ERROR: Need to set ne=1 in ECPARM.f'
        stop
      endif

c First do the header consistency check

      if(exptyp(12:14).eq.'AER'.and..not.coupled_aero)then !RS file has aero data
        print*,'Warning: RS file has aerosol data but coupled_aero = F'
        exptyp(12:14)='    '
      elseif(exptyp(12:14).eq.'    '.and.coupled_aero)then !No data in RS file
        print*,'WARNING: Setting aerosol stuff to 0 in filerd!'
        bigxt(:,:,:,:)=0.
        bigxtm(:,:,:,:)=0.
      endif

c Now read the RS data, using the header as a guide

      char=exptyp(15:15)
      if(char.ge.'0'.and.char.le.'9')then
        read(char,'(i1)')ntracers
        if(ntracers.ne.ntrac)then
          print*,'Error: No. of aerosol tracers incorrect in RS file.'
          stop
        endif
      elseif(char.eq.'O')then
        print*,'WARNING: No. of aerosol tracers not stated in RS file!'
        print*,'May be incompatible with model version!'
      endif

      if(exptyp(12:14).eq.'AER')then
        read(ifiler)bigxt
        read(ifiler)bigxtm
      endif

      if(coupled_aero)then
        exptyp(12:14)='AER'
        write(char,'(i1)')ntrac
        exptyp(15:15)=char
      endif

c******************************************************************************

C**** --- END OF DATA INPUT FROM RESTART FILE ---

      close(unit=ifiler)
      write(6,92)ndays,ifiler
  92  format(10x,'model initialised at day',i9,' from file ',i2)

c Input tracer fields con, conm

      if ( sltrace ) then
        open(7,file='trace.in',form='unformatted',status='OLD'
     &       ,iostat=ierr)
        call errcheck(ierr,'trace.in         ','filerd    ')
        read(7) ntracef,ndaysf,lwf
        if(ntracef.ne.ntrace)then
          print*,'Error: ntrace in tracer restart file incorrect'
          print*,' ntracef, ntrace',ntracef,ntrace
          stop
        elseif(ndaysf.ne.ndays)then
          print*,'Error: ndays in tracer restart file incorrect'
          print*,' ndaysf, ndays',ndaysf,ndays
ctemp     stop
        elseif(lwf.ne.lw)then
          print*,'Resolution of tracer input file incorrect'
          stop
        endif
        read(7)conm
        read(7)con
        close(7)
      endif

      do 438 k=1,4
      fldm4x(k)=0.0
  438 fld4x(k)=0.0

c******************************************************************************

      if(ltrace)then

c READ POSITION OF PARTICLES (PARTICLE TRACER EXPERIMENT)
        open(7,file='partpos',form='unformatted',status='old'
     &       ,iostat=ierr)
        call errcheck(ierr,'partpos          ','inital    ')
        read(7)npts
        print * ,' reading inital positions for ',npts,' particles'
        if(npts.gt.nps) stop 'ltrace - filerd'
        read(7)lonpos
        read(7)latpos
        read(7)sigpos
        close(7)

      endif

c******************************************************************************

      if(newriver)then

c Restart file for the new river routing scheme.
        open(7,file='river.in',form='unformatted',status='OLD'
     &       ,iostat=ierr)
        call errcheck(ierr,'river.in         ','filerd    ')
        read(7)ndaysf
        if(ndaysf.ne.ndays)then
          print*,'Error: ndays in river restart file incorrect'
          print*,' ndaysf, ndays',ndaysf,ndays
          stop
        endif
        read(7)((resvr1(mg,lg),mg=1,lon),lg=1,lat2) ! m of water
        close(7)
      endif

      return
      end
C---------------------------------------------------------------------
      subroutine setr (spec,
     &                 specr,speci)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      complex spec(nw) !1-D read field
      real specr(lw,mw)
      real speci(lw,mw)

C Local work arrays and variables
      integer ll
      integer llmax
      integer llmax1
      integer lx
      integer mm


C Start code : ----------------------------------------------------------

c To turn 1-D input field into triangular or rhomboidal format

      lx=0
      do mm=1,mw
        if(trunc.eq.'T')then
          llmax=lw-mm+1
        else
          llmax=lw
        endif
        do ll=1,llmax
          lx=lx+1
          specr(ll,mm)=real(spec(lx))
          speci(ll,mm)=aimag(spec(lx))
        enddo
        llmax1=llmax+1
        if(llmax1.le.lw)then
          do ll=llmax1,lw
            specr(ll,mm)=0.0
            speci(ll,mm)=0.0
          enddo
        endif
      enddo

      return
      end
