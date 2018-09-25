c CABLE initialization is called when lsm_type = "cable"
c AJA 2009/01/29
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Modified to remove the reading of the file "c9details", and the subsequent
c writing of the data to standard output, as this is completely unnecessary.
c SJP 2004/01/05
c
c Modified so that HMREAD is only called if QFLUX=T, as a constant mixed-layer
c ocean depth of 100m is now used otherwise.
c SJP 2003/05/29
c
c $Log: atstart.f,v $
c Revision 1.15  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.14  2001/02/22 06:50:59  rot032
c Fix unit no.
c
c Revision 1.13  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.12  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.11  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.10  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.9.1.1  2000/06/20 02:27:09  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.9  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.8  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.7  1998/12/10  00:55:36  ldr
c HBG changes to V5-1-21
c
c Revision 1.6  1997/12/19  02:03:12  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.5  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.4  1996/10/24  01:03:32  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.3  1996/06/13  02:08:56  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.2  1996/04/17  01:32:52  ldr
c Changes from TIE to allow initialization of monthly mean netcdf files
c from within model.
c
c Revision 1.1  1996/02/19  03:53:59  ldr
c Initial revision
c
c     INPUT/OUTPUT:
c     Input:   from common/fewflags in FEWFLAGS.f
c                  idayp - interval in days b/t calls to printing routines
c                  idyn  - if T and semice=T, leads=T, use dynamical 
c                          sea-ice model
c
c              from common/glmean in GLMEAN.f
c                  totkei - total kinetic energy
c                  totkzi - total kinetic energy, zonal average
c
c              from common/timex in TIMEX.f
c                  dt - time step in seconds
c
c              from arguments
c                  exptyp - restart file header
c
c     Output:  from common/atm2oc in ATM2OC.f
c                  nmth - month counter for ocean routines
c
c              from common/atdat in this subroutine
c                  iday     - logical variable, T if end of day
c                  idaypsav - nterval in days b/t calls to printing routines
c                  ndadd    - counts number of (day) additions to print files
c                  nrunday  - counts days since start of this run
c                  tdt      - 2 X timestep (dt)
c
c              from common/glmean in GLMEAN.f
c                  grunof   - global mean runoff
c
c     In/Out:  from common/fewflags in FEWFLAGS.f
c                  lcouple - flag to run model in coupled mode, default F
c
c              from common/timex in TIMEX.f
c                  month  - month counter
c                  nsteps - time step counter
c

      subroutine atstart(exptyp)
      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'RDPARM.f' ! included to use RADISW.f

C Argument list
      character exptyp*50

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'ATM2OC.f'  !Parameters passed back to ogcm
      include 'DATICE1.f'
      include 'FEWFLAGS.f'
      include 'GLMEAN.f'   ! Common block glmean contains global means
      include 'TIMEX.f'
      
      include "LSMI.f"      ! the array imsl is made available  
      include "GAUSL.f"     ! latitude and longitude computation variables are made available       
      include "RADISW.f"    ! For CO2 concentration to pass to CABLE
 
      real tdt
      integer ndadd,nrunday,idaypsav,mstepsav
      logical iday,endofmonth
      common/atdat/tdt,ndadd,nrunday,idaypsav,mstepsav,iday,endofmonth

      integer idfrz
      common/icedvz/idfrz(lat,2)

C Local work arrays and variables
      integer lg
      integer m
      integer mms
      integer ns

C Local data, functions etc

C Start code : ----------------------------------------------------------


c     CALL TIMER('ATSTART ',3)
C****                                                         *
C**** INITIALIZE THE NL LEVEL DAR CLIMATE MODEL               *
C****                                                         *

      idaypsav=idayp  !idaypsav is the  original idayp

      mms=1
      go to 1400
          
c SET UP PARAMETERS (OPENFL now called from INITAL)

c----
c---- Enter here if starting new month after month m=1
      entry atsmon(m)
      mms=m
c----
c---- Note only some of the following calls are needed after the 
c---- initial call (i.e. when mms>1)
c----
 1400 continue
c SET UP GAUSSIAN VALUES & LEGENDRE POLYS
      If (mms.eq.1) call gauleg

c      Prelim call to icesetup to create blonu,blatu etc,
c       which are needed by ncinit called from inital
      call icesetup(1)

      call inital(mms)
      call initax

c       Initialize dynamical ice model
      if(idyn)then
         if(.not.lcouple)call flatset(month)
         call flatset1
      endif
c For DATICE1.f
      rhocw1=4.1e6
      if(.not.semice)rhocw1=0.0
      if(lcouple)then
        if(lw.eq.22)then !Old Mk2 treatment
          antk=0.15e-4
          arck=0.15e-4
        else
          antk=0.075e-4
          arck=0.075e-4
        endif
      else
        antk=0.15e-4
        arck=0.15e-4
      endif
      diz=25.
      do ns=1,2
      do lg=1,lat
        idfrz(lg,ns)=0
      enddo
      enddo

c.... Only call HMREAD if QFLUX=T
      if (qflux) call hmread(month)

C**
C*    GET RESTART DATA FROM DISK FILE AND READ DATA FILES
      if(mms.eq.1)then
        call filerd(exptyp)
        call datard
        if(coupled_aero) call ecread
      end if

      nsteps=0
      iday=.true.

      dt=mstep*60.0
      tdt=2.0*dt
      if(ncepstrt.eq.2)then
        dt=0.5*mstep*60.0
        tdt=2.0*dt
      endif
      grunof=0.0
c.... nrunday counts days so far this month
      nrunday=0
c.... ndadd counts number of (day) additions to print files
c.... (can be reset to zero each idayp days)
      ndadd=0
c--  link steps to ocean model
      nmth=month

      If (mms.eq.1) Then
c Initialize for Fels-Schwarzkopf Radiation Scheme

        call initfs(0) !Argument 0 means this is call at start of run only

c Initialization for dynamical ice model

        if(idyn) call icesetup(2)
        call icecon
      End If

c* update surface albedo and sib fields each month (see also filerd)
      if(mms.gt.1)then
        call datard
        if(coupled_aero) call ecread
      endif

c* Initialise CABLE
      if(mms.eq.1)then
        if (lsm_type .eq. "cable")
     &    call CABLE_init(imsl,lon,ln2,lat,sia,mstep,rrvco2,iyear,month
     &                                                           ,mins)
      end if

c* In coupled mode, read in the Surface Temp and Salinity correction data
      if(lcouple.and.fluxadj)call tmread(month)
c* Read in the runoff descent data if coupled of savefcor
      if(lcouple.or.savefcor) call landrun
C*    ZERO OUT STATISTIC ARRAYS
C*    AT BEGINNING OF EACH MONTH
      call zerost(0)

c WRITE INITIAL VALUE OF GLOBAL K.E. AND HEADER TO OUTPUT FILE

      call energy(.false.,0)   !To calculate totkzi,totkei
      write(6,438)totkzi,totkei
 438    format(' Total kzon = ',e17.11,' Total keddy = ',e17.11)
      write(6,*)
      write(6,21)
 21   format(1x,'    DAY',2x,'HR',1x,'  PMSL',1x,'SBAL',1x,'T(9)',
     &     1x,' T*  ',1x,' Tmn',' Tmx',' Kzn',' Ked',' LWG',
     &     ' LWt',' SWg',' SIt',' SOt',' CLFR',
     &     1x,'AL',1x,'HC',
     &     1x,'MC',1x,'LC',1x,'TC',1x,'EVP',1x,'RAI',1x,'HFLX',
     &     1x,'WATR',1x,'SNOD',1x,'SICD',1x,'TEMP',' K/yr',1x,'Umx')
c     CALL TIMER('ATSTART ',4)

      return
      end
