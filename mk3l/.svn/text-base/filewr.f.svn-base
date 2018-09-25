c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c nsib comparison is modified to lsm_type = "nsib "
c AJA 2009/01/22
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Removed the line "include 'OCEANPARS.f'", as none of the parameters defined
c in this header file are referenced.
c SJP 2007/05/31
c
c Modified for five-digit year numbers.
c SJP 2004/09/29
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c $Log: filewr.f,v $
c Revision 1.58  2001/02/28 04:36:37  rot032
c Further tidy ups from HBG
c
c Revision 1.57  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.56  2001/02/12 05:39:49  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.55  2000/08/30 04:18:59  rot032
c Add fname to call writecif.
c
c Revision 1.54  2000/08/28 05:51:45  rot032
c Add a routine to dump cif files.
c
c Revision 1.53  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.52  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.51.1.1  2000/06/20 02:27:09  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.51  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.50  1998/12/10  00:55:37  ldr
c HBG changes to V5-1-21
c
c Revision 1.49  1997/12/23  01:26:59  ldr
c Hal's restart fixes.
c
c Revision 1.48  1997/12/23  00:23:35  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.47  1997/12/17  23:22:47  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.46.3.1  1997/12/19  02:03:13  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.46  1997/01/07  03:47:51  ldr
c Make write of new NSIB variables conditional, so that NSSTOP < 0 works.
c
c Revision 1.45  1996/12/19  23:40:50  ldr
c EAK's snow albedo changes merged and debugged.
c
c Revision 1.44  1996/10/24  01:02:44  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.43  1996/03/21  03:18:41  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.42  1996/02/19  04:09:48  ldr
c Generalize for 24 levels.
c
c Revision 1.41  1995/11/23  06:09:40  ldr
c Merge of LDR and EAK changes since V4-7-21l.
c
c Revision 1.40  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.39.1.1  1995/11/23  05:52:11  ldr
c Fix for T63 header and new nsstop=-2 option.
c
c Revision 1.39  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.32.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.38.1.1  1995/11/10  05:23:03  ldr
c EAK's changes to soil types (increase to 9 soil types).
c
c Revision 1.38  1995/05/02  07:49:15  ldr
c Merge of HBG stuff in V4-5-43h3 into latest version.
c
c Revision 1.37  1995/03/14  05:02:17  ldr
c Merge of INS changes for new ncput into LDR's latest version.
c
c Revision 1.36  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.35.2.1  1995/03/14  04:58:52  ldr
c Changes for new version of ncput which handles reduced T63 arrays.
c
c Revision 1.35  94/08/08  17:21:17  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.34  94/08/04  17:07:00  ldr
c Merge of LDR and HBG changes.
c 
c Revision 1.33  94/08/04  10:36:17  ldr
c Add cloud fraction to RS file - only needed so that t-1 value can be 
c used in estimate of Ri.
c 
c Revision 1.32.1.1  94/08/04  16:55:22  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.32.1.2  1995/05/02  07:31:32  ldr
c Extra changes from HBG to files which already went in as V4-5-43h.
c
c Revision 1.32  94/05/13  14:55:41  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c 
c Revision 1.31  93/12/06  16:55:28  ldr
c Changes to V4-4-45l to get T63 running on SGI (Seca).
c 
c Revision 1.30  93/11/29  11:41:59  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.29  93/11/22  14:25:36  ldr
c Correction for 18 level version.
c 
c Revision 1.28.1.1  93/11/29  11:38:25  ldr
c Changes to V4-4-32l from HBG for coupled model
c
c     INPUT/OUTPUT
c     Input:  from common/fewflags in FEWFLAGS.f
c                  qcloud  - T if using prognostic cloud scheme
c                  qflux   - set to T for qflux run, default F
c                  semice  - if T, use the Semtner sea-ice model
c                  sltrace - if T compute tracer transport by semi-
c                            Langrangian advection and by mixing and
c                            convection in the vertical
c
c             from common/files in FILES.f
c                 orfilename - name of output restart file
c
c             from common/lsmi in LSMI.f
c                 imsl - land surface type indicator
c
c             from common/masiv1 in this subroutine
c                 ronmx - global array of measure of east-west surface stress
c                 sonmx - global array of measure of north-sth surface stress
c                 iphys ) - current phys time step and previous phys time
c                 iphysm)   step, these indices swap around
c
c             from common/masiv3 in MASIV3.f
c                 rgsav - incoming long-wave radiation
c                 sgamp - amplitude of net solar radiation at the ground
c                 sgsav - net solar radiation at ground
c
c             from common/masiv4 in MASIV4
c                 savegrid - global surface variables
c
c             from common/masiv5 in this subroutine
c                 statsice - global statistics relating to sea-ice scheme
c
c             from common/pmasiv3 in PMASIV3.f
c                 prgsav - rgsav over leads(see above)
c                 psgamp - sgamp over leads  "    "
c                 psgsav - sgsav over leads  "    "
c
c             from common/qcloud1 in QCLOUD1.f
c                 qfb  - cloud ice at t     qlb  - cloud liquid water at t
c                 qfbm - cloud ice at t-1   qlbm - cloud liquid water at t-1
c                 cfb  - cloud fraction
c
c             from common/qflxdat in QFLXDAT.f
c                 ochfa - ocean heat flux average 
c
c             from common/rmgrid in RMGRID.f
c                 rgt  - moisture tendency
c                 rmg  - pressure weighted moistures at current timestep
c                 rmmg - pressure weighted mositures at previous timestep
c
c             from common/rshead in this subroutine
c                 header2 - second header on restart file
c
c             from common/surf in this subroutine
c                 aftfh  - reciprocal of air resistance for ground 
c                          covered with foliage
c                 aftfhg - reciprocal of air resistance for bare ground
c                 ocondx - precipitation at previous timestep
c                 opreci - frozen precipitation at previous timestep
c
c             from common/surf1 in SURF1.f
c                 mc   - moisture depth on canopy
c                 tgf  - canopy temperature
c                 tgg  - bare ground temperature
c
c             from common/timex in TIMEX.f
c                 mstep  - timestep in minutes
c                 ndays  - day counter for model runay counter for model run
c                 nsteps - time step counter
c
c             from common/traceblk in TRACEBLK.f
c                 con  - concentration of atmospheric tracers at t
c                 conm - concentration of atmospheric tracers at t-1
c                 nslfirst - if 1, starting from data set with no backward
c                            step (mainly for use by IGW)
c
c             from common/uvpgd in UVPGD.f
c                 sdot - d(sigma)/d(t)
c
c             from argument
c                 exptyp - restart file header
c
c     In/Out: from common/files in FILES.f
c                 str - character string YYYYYMM.XXXXX where XXXXX is runtype
c
c             from common/fldmri in FLDMRI.f
c                 pmi  - spectral surface pressure, imaginary part
c                 pmr  - spectral surface pressure, real part
c                 psimi - stream function, imaginary part
c                 psimr - stream function, real part 
c                 temi - spectral temp, imaginary part at t-1
c                 temr - spectral temp, real part at t-1
c                 xhimi - velocity potential at t-1, imaginary part
c                 xhimr - velocity potential at t-1, real part
c
c             from common/fldri in FLDRI.f
c                 pi -  spectral pressure, imaginary part
c                 pr -  spectral pressure, real part
c                 psii - spectral stream function, imaginary part
c                 psir - spectral stream function, real part
c                 tei - spectral temp, imaginary part
c                 ter - spectral temp, real part
c                 xhii - velocity potental, imaginary part
c                 xhir - velocity potental, real part
c
c             from common/mountn in this subroutine
c                 phisi - spectral surface geopotential height,imaginary part
c                 phisr - spectral surface geopotential height,real part
c
c 
      subroutine filewr(exptyp)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f'

C Argument list
      character*50 exptyp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'DIFM.f'
      include 'ECTRAC.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'FLDRI.f'
      include 'FLDMRI.f'
      include 'LSMI.f'
      include 'MASIV3.f'
      include 'MASIV4.f'
      include 'PMASIV3.f'
      include 'QFLXDAT.f'
      include 'QCLOUD1.f'
      include 'RMGRID.f'
      include 'SURF1.f'
      include 'TIMEX.f'
      include 'TRACEBLK.f'      !Tracer arrays con, conm
      include 'UVPGD.f'

      real ronmx,sonmx
      integer iphysm,iphys
      common/masiv1/ronmx(ln2,nl,lat,2),sonmx(ln2,nl,lat,2)
     &             ,iphysm,iphys

      real statsice
      common/masiv5/statsice(ln2,14,lat)

      real phisr,phisi
      common/mountn/phisr(lw,mw),phisi(lw,mw)

      real resvr1,resrem,srunoff,sresvr1
      common/rivers/resvr1(0:lon+1,0:lat2),resrem(lon,lat2)
     & ,srunoff(lon,lat2),sresvr1(lon,lat2)

      character*50 header2
      common/rshead/header2

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      real griddata(ln2,ngrid)
      real ron(ln2,nl),son(ln2,nl)
      real qout(ln2,lat)
      real stice(ln2,9)

      complex spec(nw)

      integer icomp
      integer irestf
      integer k
      integer lg
      integer mg

      real consx

C Local data, functions etc

C Start code : ----------------------------------------------------------
C
C FUNCTION : TO ARCHIVE DAILY INFORMATION OF NL LEVEL MODEL.
C

C**** TO OUTPUT RESTART DATA TO RESTART FILE
C**** AND TO ACCUMULATE OCEAN DATA ON TAPE22.
C**** OVERWRITE IFF OUTPUT FILENAMES EQUAL INPUT FILENAMES (IN NAMELIST)

      irestf=18

c****  OPEN THE RESTART FILES

      open(unit=irestf,file=orfilename,form='unformatted')

C**** READ/WRITE HEADER

      write(irestf)exptyp
C****
      write (irestf) ndays,mstep
      write(irestf)header2

c**** output spectral geopotential
      call setw (phisr,phisi,spec)
      write (irestf) spec
C**** READ/WRITE PSI(flux)
      do 10 k=1,nl
         call setw (psir(1,1,k),psii(1,1,k),spec)
 10      write (irestf) spec
      do 12 k=1,nl
         call setw (psimr(1,1,k),psimi(1,1,k),spec)
 12      write (irestf) spec
C**** READ/WRITE XHI(flux)
      do 14 k=1,nl
         call setw (xhir(1,1,k),xhii(1,1,k),spec)
 14      write (irestf) spec
      do 16 k=1,nl
         call setw (xhimr(1,1,k),xhimi(1,1,k),spec)
 16      write (irestf) spec
C**** NOW SET UP TO WRITE (P*T') NOT (P*T)
      do 18 k=1,nl
         call setw (ter(1,1,k),tei(1,1,k),spec)
 18      write (irestf) spec
      do 20 k=1,nl
         call setw (temr(1,1,k),temi(1,1,k),spec)
 20      write (irestf) spec
C**** READ/WRITE P*
      call setw (prr,pri,spec)
      write (irestf) spec
      call setw (pmr,pmi,spec)
      write (irestf) spec
C**** OUTPUT GRID MOISTURE + LAST TENDENCY
      do 30 k=1,nl
      do 31 lg=1,lat
      do 31 mg=1,ln2
   31 qout(mg,lg)=rmg(mg,k,lg)
   30 write (irestf)qout
      do 32 k=1,nl
      do 33 lg=1,lat
      do 33 mg=1,ln2
   33 qout(mg,lg)=rmmg(mg,k,lg)
   32 write (irestf)qout
      do 34 k=1,nl
      do 35 lg=1,lat
      do 35 mg=1,ln2
   35 qout(mg,lg)=rgt(mg,k,lg)
   34 write (irestf)qout
C**** OUTPUT SURFACE TYPE INDICATOR ARRAY
      write (irestf) imsl
C****
C****     OUTPUT PHYSICAL STATISTICS
C****
      do 80 lg=1,lat
      do 85 k=1,ngrid
      do 85 mg=1,ln2
   85 griddata(mg,k)=savegrid(mg,k,lg)
C**** OUTPUT PHYSICAL RESTART DATA TO RESTART FILE.
      write (irestf) griddata
      do 86 k=1,nl
      do 86 mg=1,ln2
      ron(mg,k)=ronmx(mg,k,lg,iphysm)
   86 son(mg,k)=sonmx(mg,k,lg,iphysm)
      write(irestf)ron,son
 80   continue

c Extra stuff for SPO's sea-ice model

      if(semice)then
        do 91 lg=1,lat
        do 90 k=1,9
        do 90 mg=1,ln2
 90     stice(mg,k)=statsice(mg,k,lg)
 91     write(irestf)stice
      endif

c**** output dpsi,dxhi (spectral diffusion energy
c****    loss conversion to temp tendancy terms)
c**** Split vertical levels now standard for all resolutions

      do k=1,nl
        call setw (dpsir(1,1,k),dpsii(1,1,k),spec)
        write (irestf) spec
      enddo
      do k=1,nl
        call setw (dxhir(1,1,k),dxhii(1,1,k),spec)
        write (irestf) spec
      enddo

c Write vertical velocity (sdot)

      write (irestf) sdot
      if (qflux) write(irestf) ochfa

c Write nsib fields

      SELECT CASE(lsm_type)
      CASE("nsib ") 
        write(irestf) wb
        write(irestf) wbice
        write(irestf) tggsl
        write(irestf) tggsn
        write(irestf) ssdnn
        write(irestf) ssdn3
        write(irestf) smass
        write(irestf) gflux
        write(irestf) sgflux
        write(irestf) isflag
        write(irestf) snage
        write(irestf) osnowd
        write(irestf) tgf
        write(irestf) mc
        write(irestf) aftfh 
        write(irestf) aftfhg 
        write(irestf) ocondx
        write(irestf) pmc
      CASE("cable")
        write(irestf) wb
        write(irestf) wbice
        write(irestf) tggsl
        write(irestf) tggsn
        write(irestf) ssdnn
        write(irestf) ssdn3
        write(irestf) smass
        write(irestf) gflux
        write(irestf) sgflux
        write(irestf) isflag
        write(irestf) snage
        write(irestf) osnowd
        write(irestf) tgf
        write(irestf) mc
        write(irestf) aftfh 
        write(irestf) aftfhg 
        write(irestf) ocondx
        write(irestf) pmc
      CASE DEFAULT ! old land surface scheme removed
        PRINT *,'***  Invalid LAND SCHEME: Program Exiting'
        STOP
      END SELECT

c Add radiation arrays to restart file.

      write(irestf)sgsav
      write(irestf)rgsav
      write(irestf)sgamp
      write(irestf)psgsav
      write(irestf)prgsav
      write(irestf)psgamp

c Write cloud water and cloud ice

      if(qcloud)then
        write(irestf)qlb
        write(irestf)qlbm
        write(irestf)qfb
        write(irestf)qfbm
        write(irestf)cfb
        write(irestf)opreci
      endif

      if(coupled_aero)then
        write(irestf)bigxt
        write(irestf)bigxtm
      endif

C****
C**** --- END OF REWRITE OF RESTART FILE ---
C****
      endfile irestf
      close(unit=irestf)
      write(6,140)irestf,ndays
  140 format(t10,'restart file(',i2,') written at day',i9)

c Output tracer fields con, conm

      if ( sltrace ) then
        open(51,file='trace.out',form='unformatted')
        write(51) ntrace,ndays,lw
        write(51)conm
        write(51)con
        write(51)nslfirst
        close(51)
        consx=1.0/nsteps
c---- files tT??.., t.  With 9:1 compession at T63
        if(mw.eq.64)then
          icomp=9
        else
          icomp=1
        endif
        call dmpcif(0,str,consx,icomp)
c       call dmpcif(1,str,consx,icomp)
      endif

c******************************************************************************

      if(newriver)then
c Restart file for the new river routing scheme.
        open(51,file='river.out',form='unformatted')
         write(51)ndays
         write(51)((resvr1(mg,lg),mg=1,lon),lg=1,lat2)
        close(51)
      endif

      return
      end

c******************************************************************************

      subroutine setw (specr,speci,spec)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real specr(lw,mw)
      real speci(lw,mw)
      complex spec(nw)

C Local work arrays and variables
      integer ll
      integer llmax
      integer lx
      integer mm

C Start code : ----------------------------------------------------------

c To make 1-D output field from triangular or rhomboidal format

      lx=0
      do mm=1,mw
        if(trunc.eq.'T')then
          llmax=lw-mm+1
        else
          llmax=lw
        endif
        do ll=1,llmax
          lx=lx+1
          spec(lx)=cmplx(specr(ll,mm),speci(ll,mm))
        enddo
      enddo

      return
      end

c******************************************************************************

      subroutine dmpcif(ntime,str,consx,ind)

c Routine to dump tracer con field to files

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ntime
      character*13 str
      real consx
      integer ind ! compression factor (1 or 9)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'GAUSL.f'
      include 'TIMEX.f'
      include 'TRACEBLK.f'      !Tracer arrays con, conm, conmth

C Local work arrays and variables
      real tx(lon,lat2)
      character filename*13
      character*3 varname
c---- These are the 9:1 compression arrays
      real txmid(64,lat2)
      real txcomp(64,32)    ! this is the minimal field for plotting.
      equivalence (tx(1,1),txcomp(1,1))
c----

      integer k
      integer lat9
      integer lg
      integer lgns
      integer lid
      integer lon9
      integer lx
      integer mg
      integer mga
      integer mid
      integer nt
      integer ntrace1

      real sumw
      real wlg
      real wlgm
      real wlgp

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** ROUTINE TO COLLECT TRACER STATS FROM THE 
c****      conm(-1:lonx+2,-1:lat2x+2,nl,ntrace)
c**** METHOD OF STORAGE IN TRACEBLK INTO ARRAYS SUITABLE FOR
c**** PLOTTING,ETC : EG. (LON,LAT*2) WHERE 1,1 IS (SOUTH POLE,GM).
c----
c---- ind controls whether the full field is dumped (ind.eq.1), or
c----     if a compressed (area averaged) form is created (ind.eq.9).
c---- (ind=9 for T63 only). If T63, then
c---- this routine takes data from T63 in arrays (192,96) and may
c---- compress (average) to arrays (64,32). The use of "64" etc is
c---- NOT a hangover from the R21 model!!
c----
c        lat k+1  :  x - x - x  :              : X :
c        lat k    :  x - x - x  :  reduced to  : X :
c        lat k-1  :  x - x - x  :              : X :
c
c----      by 3 point averaging (simple) in EW direction followed
c----      by a 3 point weighted average in NS direction onto lat k
c----      i.e. 9:1 reduction in data storage.

      if(ind.eq.9)then
c---- check that T63 is in use
      if (mw.ne.64) then
         print *,' wrong resolution for dmpcif 9:1 compression'
         stop
      end if
      lon9=lon/3
      lat9=lat2/3
      endif

c---- Set up start of file name
      filename='tnnn          '

c     ntrace1=1
      ntrace1=2
      do 100 nt=ntrace1,ntrace
      do 100 k=1,nl

        if(ntime.eq.0)then ! original 'a' files now with Lower Case
          if(nt.eq.1) write(filename(2:2),'(a1)') 't'
          if(nt.eq.2) write(filename(2:2),'(a1)') 'f'
          if(nt.eq.3) write(filename(2:2),'(a1)') 'v'
        elseif(ntime.eq.1)then ! original 'e' files now with Upper Case
          if(nt.eq.1) write(filename(2:2),'(a1)') 'T'
          if(nt.eq.2) write(filename(2:2),'(a1)') 'F'
          if(nt.eq.3) write(filename(2:2),'(a1)') 'V'
        endif
        write(filename(3:4),'(i2.2)') k 
c       print *,' writing ',filename
c
c Work out filename:
c Character string STR*13 contains YYYYYMM.XXXXX where XXXXX is runtype
        filename(5:13)='_'//str(9:13)//'.nc'
        varname=filename(2:4)

        do lgns=1,lat2
          do mg=1,lon
            if(ntime.eq.0) tx(mg,lgns)=conmnth(mg,lgns,k,nt)*consx
            if(ntime.eq.1) tx(mg,lgns)=con(mg,lgns,k,nt)
          enddo
        enddo

c Transfer data to plotting array : (1,1) is SP,GM
c and Write data
      if (ind.eq.1) then
c---- write full field
        call nc_put(lon, lat2, filename, varname, iyear, month, tx)
      else
c---- form 9:1 compressed field and then write
         do 30 lg=1,lat2
            do 30 mg=1,lon-2,3
            mid=(mg/3)+1
            mga=mg-1
            if (mga.eq.0) mga=lon
            txmid(mid,lg)=(tx(mga,lg)+tx(mg,lg)+tx(mg+1,lg))/3.0
 30      continue

         do 40 lg=2,lat2-1,3
            lid=(lg/3)+1
            if(lg.lt.lat)then
              wlgp=w(lg+1)
              wlg =w(lg)
              wlgm=w(lg-1)
            else
              lx=lat2+1-lg
              wlgp=w(lx-1)
              wlg =w(lx)
              wlgm=w(lx+1)
            endif
            sumw=wlgm+wlg+wlgp
            do 40 mg=1,64
            txcomp(mg,lid)=(wlgm*txmid(mg,lg-1)+
     &       wlg*txmid(mg,lg)+wlgp*txmid(mg,lg+1))/sumw
 40      continue

        call nc_put(lon9, lat9, filename, varname, iyear, month,
     $              txcomp)
      end if

  100 continue

      return
      end
