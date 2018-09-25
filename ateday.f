c Removed unnecessary "include 'MACHINE.f'"
c SJP 2001/11/22
c
c $Log: ateday.f,v $
c Revision 1.18  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.17  2001/02/12 05:39:50  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.16  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.15  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.14  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.13  1996/10/24  01:03:24  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.12  1996/06/13  02:08:41  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.11  1995/08/29  01:30:09  ldr
c HBG's corrections to his V4-7-13h, to make results with hybrid=F agree
c with previous version.
c
c Revision 1.10  1995/08/14  05:29:54  ldr
c HBG's new improved printing routines.
c
c Revision 1.9  1995/05/04  00:07:43  ldr
c Declare dcflag as type logical to avoid warnings on Fujitsu.
c
c Revision 1.8  1994/12/15  06:23:52  ldr
c Remove the endfile(6) statement which is causing problems on the Cray now.
c
c Revision 1.7  1994/08/18  14:25:50  ldr
c Remove the spurious zeroing of cloud statistics which caused zeroes to
c be printed to stats files. Sneak this into V4-6 while nobody is looking.
c
c Revision 1.6  94/08/09  12:35:27  ldr
c Rationalize zonal mean diagnostics: Use conzp flag to control HBG's
c cloud map and introduce savezonu, savezonv for dump of zmean u,v fields.
c Also remove savezcls which is fairly irrelevant.
c 
c Revision 1.5  94/08/04  17:07:09  ldr
c Merge of LDR and HBG changes.
c 
c Revision 1.4  94/07/11  12:30:35  ldr
c Rearrange variables in common block /atdat to avoid misalignment warnings
c at 64 bit on SGI.
c 
c Revision 1.3.1.1  94/08/04  16:56:52  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.3  94/01/17  11:31:59  ldr
c Make annoying endfile statement work on Cray only to keep VP happy.
c 
c Revision 1.2  93/12/23  15:31:28  ldr
c Changes to V4-4-54l from HBG for coupled model of V4-5
c 
c Revision 1.1  93/08/10  15:03:35  ldr
c Initial revision
c 
c Revision 1.38  93/06/18  14:43:30  ldr
c 
c     INPUT/OUTPUT:
c     Input:   from common/atdat in this subroutine
c                  endofmonth - counter
c                  idaypsav   - interval in days b/t calls to printing 
c                               routines
c
c              from common/fewflags in FEWFLAGS.f
c                  mlomap   - if T, calls prsmap to print map of surface temp
c                             of sea-ice & mixed layer ocean points
c                  qflux    - set to T for qflux run, default F
c                  savegbrf - if T saves various global means
c
c              from common/glmean in GLMEAN.f
c                  sidbar    - mean sea ice depth
c                  sndbar    - mean snow depth
c                  tempbar   - mean temp (all levels)
c                  totkei    - total kinetic energy
c                  totkzi    - total kinetic energy, zonal average
c                  tot_water - global integral of moisture
c
c              from common/printt in PRINTT.f
c                  conzp - if T, print a zonal height map of zonal mean
c                          cloud fraction
c                  zavgp - if T, przav looks at next 8 flags to determine 
c                          which zonal mean, time mean fields should be printed
c
c     In/Out:  from common/atdat in this subroutine
c                  ndadd   - counts number of (day) additions to print files
c                  nrunday - counts days since start of this run
c
c              from common/fewflags in FEWFLAGS.f
c                  idayp   - interval in days b/t calls to printing routines
c                  lcouple - flag to run model in coupled mode, default F
c
c     Output:  from common/atm2oc in ATM2OC
c                  ittx - ocean model counter
c
      subroutine ateday

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Global data blocks
      include 'ATM2OC.f'  !Parameters passed back to ogcm
      include 'FEWFLAGS.f'
      include 'GLMEAN.f'     ! Common block glmean contains global means
      include 'PRINTT.f'

      real tdt
      integer ndadd,nrunday,idaypsav,mstepsav
      logical iday,endofmonth
      common/atdat/tdt,ndadd,nrunday,idaypsav,mstepsav,iday,endofmonth

C Local work arrays and variables
      logical lprint

C Local data, functions etc

C Start code : ----------------------------------------------------------

c     CALL TIMER('ATEDAY  ',1)
C****
C**** PERFORMS JOBS AT END OF DAY FOR DAR CLIMATE MODEL
C****

c IF END OF DAY,UPDATE SEA TEMPS
c IF END OF DAY PRINT DAILY RAIN AND EXTREME TEMPS
c AT PRESET INTERVALS CALL PRINTING ROUTINES
c NOTE THAT NRUNDAY IS DAYS SINCE START OF THIS RUN (NOT SINCE START
c OF YEAR)
c TO BE PRECISE:
c For non-coupled runs: Call printing routines if NRUNDAY is multiple of IDAYP
c or at end of run. Exception: Don't do print on day 30 of 31 day run,
c unless IDAYP = 1 - avoids problems with 31 day months if printing
c stats after 10,20 and 31 days.
c For coupled runs: Call printing routines at end of month only.

        ndadd=ndadd+1
        nrunday=nrunday+1
        call timet
        call prdaily
        call energy(.false.,1) !To calculate totkzi and totkei
        write(6,438)totkzi,totkei
 438    format(' Total kzon = ',e17.11,' Total keddy = ',e17.11)
        if(idaypsav.eq.0)idayp=incd
        lprint=mod(nrunday,idayp).eq.0.and.(incd.eq.1.or.nrunday.ne.30)
        if(.not.lcouple.and.(lprint.or.nrunday.eq.incd)
     &     .or. lcouple.and.endofmonth) then
          if(savegbrf) write(69,1969) tot_water,sndbar,sidbar,tempbar
1969      format(1x,2f7.2,f7.3,f7.2)
          if(zavgp.and.conzp) call prclza(ndadd)
          call przav(ndadd,ittx)
          call prtcl
          call prtt(ndadd)
          if(mlomap)call prsmap
	  if(qflux.or.lcouple)call prmlomap
          ndadd=0
          call zerost(1)
        endif
C****
C**** END OF 1 DAY OF ATMOSPHERIC MODEL
C****
c     CALL TIMER('ATEDAY  ',2)

      return
      end
C---------------------------------------------------------------------
      subroutine prclza(ndadd)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ndadd

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'GAUSL.f'
      include 'TIMEX.f'

      real hlat_a,clat_a
      common/latheat/hlat_a(ln2,lat,nl),clat_a(ln2,lat,nl)

C Local work arrays and variables
      real clatz(nl),clg(nl)
      character*5 chcl(nl)

      integer inl
      integer k
      integer ksk
      integer lg
      integer lx
      integer ma
      integer mg
      integer ns

      real avg
      real clatzm
      real clgm
      real consdr
      real consrad

C Local data, functions etc

C Start code : ----------------------------------------------------------

c------------ print the zonally averaged cloud/height map ------------

      consdr=mstep/1440.0/ndadd
c     Scaling for quantities accumulated once per radiation step
c     plus factor for zonal average
      consrad=consdr*nrad
      avg=100.0*consrad/lon

c--------------------- Check number of levels ------------------------
c---- skip levels printed if nl>24
      ksk=1+(nl-1)/24

      print *,'zonal/height cloud display'
      do 10 k=1,nl
      clg(k)=0.0
      inl=nint(sig(k)*1000.0)
   10 write (chcl(k),"(' C',i3.3)") inl
      write (6,105) (chcl(k),k=ksk,nl,ksk)
 105  format(4x,'cloud &',24a5)
      clgm=0.0

      do 20 ns=1,2
      do 20 lg=1,lat
      lx=(ns-1)*(lat+1)+lg*(3-2*ns)

      do 17 k=1,nl
      clatz(k)=0.0
      do 15 mg=1,lon
      ma=mg+(ns-1)*lon
   15 clatz(k)=clatz(k)+clat_a(ma,lx,k)
   17 clatz(k)=clatz(k)*avg
      clatzm=clatz(nl)
      clatz(nl)=0.0

      do 18 k=1,nl
   18 clg(k)=clg(k)+clatz(k)*w(lx)
      clgm=clgm+clatzm*w(lx)

   20 write(6,100)lx,clatzm,(clatz(k),k=ksk,nl,ksk)
  100 format(1x,i2,1x,f5.1,2x,24f5.1)

      do 22 k=1,nl
   22 clg(k)=clg(k)*0.5
      clgm=clgm*0.5
      write (6,105) (chcl(k),k=ksk,nl,ksk)
      write(6,100)lx,clgm,(clg(k),k=ksk,nl,ksk)

      return
      end
