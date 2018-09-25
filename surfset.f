c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c nsib comparison is modified to lsm_type = "nsib "
c the if-else statements for lsm_type are modified into case statements
c AJA 2009/01/22
c
c Transferred COMMON blocks to separate header files, as follows:
c /AMM/  ->  AMM.f
c SJP 2007/05/29
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: surfset.f,v $
c Revision 1.56  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.55  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.54.1.1  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.54  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.53  2000/11/21 01:05:37  rot032
c Change ompl to sicef for coupled_aero.
c
c Revision 1.52  2000/11/17 04:35:05  rot032
c Correct def'n of ompl for aerosol model.
c
c Revision 1.51  2000/11/14 06:43:01  rot032
c Fixes to HBG code for kaos.
c
c Revision 1.49  2000/06/20 02:08:34  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.48  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.47  1998/12/10  00:55:51  ldr
c HBG changes to V5-1-21
c
c Revision 1.46  1998/05/26  05:29:40  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.45  1997/12/23  00:23:38  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.44  1997/12/17  23:22:54  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.43  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.42.2.2  1997/12/23  00:10:59  ldr
c Hal's restart fixes.
c
c Revision 1.42.2.1  1997/12/19  02:03:17  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.42  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.41  1996/12/19  23:40:50  ldr
c EAK's snow albedo changes merged and debugged.
c
c Revision 1.40  1996/11/27  04:05:14  ldr
c Reduce melting snow and melting ice albedos if nl.ge.18.
c
c Revision 1.39  1996/10/24  01:03:17  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.38.1.1  1996/12/05  03:42:34  ldr
c Fix to snow albedo from EAK.
c
c Revision 1.38  1996/03/21  03:19:07  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.37  1995/08/31  04:30:48  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.36  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.31.3.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.35  1995/05/02  07:49:15  ldr
c Merge of HBG stuff in V4-5-43h3 into latest version.
c
c Revision 1.34  1995/02/08  05:16:10  mrd
c Removed resetting of wg under snow. This affected the soil moisture budget.
c
c Revision 1.33  1994/08/08  13:16:28  ldr
c Push the debug prints down into individual routines.
c
c Revision 1.32  94/08/04  16:56:37  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.32.1.1  1995/05/02  07:31:32  ldr
c Extra changes from HBG to files which already went in as V4-5-43h.
c
c Revision 1.31  94/03/30  12:58:48  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.30  94/03/23  09:35:47  mrd
c Removed the extra zenith angle calculation by avsol by moving call to
c zenith in radin and adding zenith angle as an argument to surfset.
c 
c Revision 1.29.1.1  94/03/30  12:35:16  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.29  93/12/17  15:33:57  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.28  93/12/06  17:32:20  ldr
c Get rid of some irritating warnings on the VP.
c 
c Revision 1.27  93/11/03  11:44:34  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.26  93/08/19  15:10:58  ldr
c Minor cosmetic changes.
c 
c Revision 1.25  93/08/03  11:32:52  ldr
c Replace intrinsic functions with generic forms where possible.
c 
c Revision 1.24  93/07/27  17:13:48  ldr
c Changes to ensure coherent approach to pure restarts, monthly
c interpolation factors..
c 
c Revision 1.23  93/07/12  14:27:49  ldr
c Merge of SPO's changes to V4-3-10l with other changes since then.
c 
c Revision 1.22  93/07/06  16:08:44  ldr
c      New map for ice depth averaged over leads area. This requires extra array
c      in common block cloudm1. The mapping in surfset.f reqires pl(mg).
c 
c Revision 1.21.1.1  93/07/12  14:15:08  ldr
c Minor changes from SPO for coupled model V4-4.
c 
c Revision 1.21  93/06/22  18:11:16  ldr
c LDR's new scheme for handling IMSL under freezing melting conditions - 
c No need for special treatment on SGI with this scheme.
c 
c Revision 1.20  93/06/18  14:45:14  ldr
c Small seaice change from SPO.
c 
c Revision 1.19  93/06/16  14:28:09  ldr
c Increased snow/ice albedos (HBG).
c 
c Revision 1.18  93/03/12  10:07:27  ldr
c Minor SPO changes.
c 
c Revision 1.17  93/02/03  12:45:10  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
c Revision 1.16  92/12/09  15:24:06  ldr
c Put /datice in include file and change name to datice1.
c 
c Revision 1.15  92/12/07  15:34:10  ldr
c Increase NH snow albedo over ice from 0.7 to 0.8.
c 
c Revision 1.14  92/11/20  15:43:15  ldr
c Merge of HBG's 1.12.1.1 and 1.13.
c 
c Revision 1.13  92/10/22  11:42:54  ldr
c SPO's changes to allow moisture from leads to affect model (includes moving
c around the radiation scheme.) These changes include V4-0-1l changes relating
c to coupled model where applicable.
c 
c Revision 1.12.1.1  92/11/20  15:41:42  ldr
c Modified snowmelt albedos. (HBG)
c 
c Revision 1.12  92/09/01  12:01:28  ldr
c Replace 271.2 with tfi (temperature at bottom of ice) and put it in common
c block /freeze/. tfi differs b/w coupled and uncoupled runs.
c 
c 
c Revision 1.11  92/08/31  17:10:36  ldr
c Make these latest albedo changes conditional on semice.eq..true.
c 
c Revision 1.9  92/06/30  14:57:41  ldr
c SPO's changes to LDR's initial version of sea-ice stuff.
c 
c Revision 1.8  92/06/16  12:03:07  ldr
c Rewritten to pass physical variables as arguments rather than in common,
c in order to implement sea-ice model.
c 
c Revision 1.7  92/05/12  15:38:22  ldr
c Changes to allow write of surface fluxes needed to calculate flux
c correction for coupled model. (If statsflag is true.)
c 
c Revision 1.6  92/05/11  16:41:46  ldr
c Use monthly albedo data set for R21 too.
c 
c Revision 1.5  92/04/22  11:57:43  ldr
c Generalized for R21/R42 and Cray/SGI.
c 
c Revision 1.4  91/05/17  14:48:16  mrd
c Modfified for qflux model.
c 
c Revision 1.3  91/05/15  11:11:26  ldr
c Modified for coupled model.
c 
c Revision 1.2  91/03/13  13:00:47  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:38:10  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT
c     Input:   from common/compsst in this subroutine
c                  tsst1 - sea surface temp at beginning of current month
c                  tsst2 - sea surface temp at end of current month
c
c              from common/fewflags in FEWFLAGS.f
c                  debug    - flag to control column debugging
c                  insdebug - flag to control hemisphere debugging
c                  lgdebug  - flag to control latitude debugging
c                  mgdebug  - flag to control longitude debugging
c                  lcouple  - flag to run model in coupled mode, default F
c                  leads    - if T and semice=T allow leads in sea-ice model
c                  lsm_type - if "nsib ", use "New SIB" land surface scheme
c                  qflux    - set to T for qflux run, default F
c                  semice   - if T, use the Semtner sea-ice model
c
c              from common/freeze in FREEZE.f
c                  tfi - temperature at bottom of ice
c
c              from common/printt in PRINTT.f
c                  gwicm - if true, print maps of convection & rainfall
c
c              from common/timex in TIMEX.f
c                  mins  - current model time in mins
c                  mstep - timestep in minutes
c                  ratlm - fraction of month elapsed, runs 0-1
c
c              from common/floes
c                  pl - proportion of leads    siced - seaice thickness
c
c              from arguments
c                  ns    - hemisphere index      lg - latitude index
c                  alsf1  - surface albedo at start of month
c                  alsf2  - surface albedo at end of month
c                  tstar1 - sea surface temperature at start of month
c                  tstar2 - sea surface temperature at end of month
c                  cosz   - mean cosine of zenith angle
c                  z4     - grid elevations [non-spectral]
c                  wb - soil moisture at the previous time step
c
c     Output:  from arguments
c                  cie - internal conductivity of ice/depth
c                  cls - factor for latent heat of fusion
c                  gamms - density*specific heat*depth 
c
c     In/Out:  from common/datice in DATICE1.f
c                  diz - depth of surface layer of ocean
c                  hoice -heat flux from ocean to ice
c                  rhocw1 - const = 4.1e6
c
c              from common/logimsl in LOGIMSL.f
c                  cice - seaice
c                  land,sea - logical variables for surface type
c                  mlo  - mixed layer ocean
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from arguments
c                  bs  - internal conductivity of snow/depth
c                  ci - set to cie (see above)
c                  tintp - temperature at sea points
c                  tstar - surface temperature
c                  als - surface albedo        snowd - snow depth
c                  tg  - surface temperature   wetfac - soil westness
c 
c     In:      from common/surf1 in SURF1.f
c                  osnowd - snow depth at the previous dt
c                  sigmf - vegetation fraction 
c                  tgg - bare ground temperature 
c
c     Out:     from common/surf1 in SURF1.f
c                  ssdnn   - snow density 
c                  snage   - snow age 
c
c     In:      from common/vegdat
c                  vegt - vegetation type
c
c
      subroutine surfset(lg,z4,tstar1,tstar2,alsf1,alsf2,
     &                   wg2,he,snowd,siced,cosz,tstar,wg,tintp,
     &                   tg,als,cls,wetfac,cie,ci,bs,gamms,sicef)
 
      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /VEGDAT/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      real wgmax
      parameter (wgmax=0.36)

C Argument list
      integer lg
      real z4(ln2)
      real tstar1(ln2)
      real tstar2(ln2)
      real alsf1(ln2)
      real alsf2(ln2)
      real wg2(ln2)
      real he(ln2)
      real snowd(ln2)
      real siced(ln2)
      real cosz(ln2)
      real tstar(ln2)
      real wg(ln2)
      real tintp(ln2)
      real tg(ln2)
      real als(ln2)
      real cls(ln2)
      real wetfac(ln2)
      real cie(ln2)
      real ci(ln2)
      real bs(ln2)
      real gamms(ln2)
      real sicef(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)

      include 'LOGIMSL.f'

      real rsmin,z0m,sigmf,vegt
      common/vegdat/ rsmin(ln2),z0m(ln2),sigmf(ln2),vegt(ln2)

C Global data blocks
      include 'CHMAP.f'

      include 'DATICE1.f'
        real antarck(2)
        equivalence (antk,antarck)

      include 'FEWFLAGS.f'
      include 'FREEZE.f'
      include 'LSMI.f'
      include 'PRINTT.f'
      include 'QFLXDAT.f'
      include 'SURF1.f'
      include 'TIMEX.f'
      include 'AMM.f'

      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel

      real tsst1,tsst2,tdtm1,tdtm2
      common/compsst/tsst1(ln2,lat),tsst2(ln2,lat)
     &              ,tdtm1(ln2,lat),tdtm2(ln2,lat)

      real zisum
      integer jice
      common/iceprt/zisum(lat,2),jice(lat,2)

      real prhos,pcss,pcnsd,pbch,pssat,pswilt,phsbh,psfc,psucs
     & ,phyds,pcdr3,prlaim,prlais,pslveg,pscveg
      integer ip2bp3,ipbp2
      common/psoilvegpr/prhos(ln2,lat),pcss(ln2,lat),
     &pcnsd(ln2,lat),pbch(ln2,lat),pssat(ln2,lat),
     &pswilt(ln2,lat),phsbh(ln2,lat),ip2bp3(ln2,lat),
     &ipbp2(ln2,lat),psfc(ln2,lat),psucs(ln2,lat),
     &phyds(ln2,lat),pcdr3(ln2,lat),prlaim(ln2,lat),
     &prlais(ln2,lat),pslveg(ln2,lat),pscveg(ln2,lat)

      character*50 header2
      common/rshead/header2

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      real sifeet(ln2)
      logical notsnalb

      integer jicex
      integer kx
      integer ma
      integer mg
      integer ns

      real ar1
      real ar2
      real ar3
      real alvo
      real aliro
      real alvd
      real alv
      real alird
      real alir
      real alss
      real al2
      real cc
      real cczen
      real celciusfrz
      real dtau
      real dnsnow
      real fcls
      real fage
      real fzen
      real fzenm
      real ftsoil
      real hfrz
      real hoicex
      real reset
      real snalb
      real saveals
      real sicft
      real snr
      real snrat
      real sstocean
      real tratlm
      real tintpx
      real ttbg
      real talb
      real tsigmfx
      real tsoil
      real wgsat
      real zisumx

C Local data, functions etc
      character*1 cint(12)
      data cint/'0','1','2','3','4','5','6','7','8','9','t','t'/

C Start code : ----------------------------------------------------------

c     surfset1 abandons leapfrog tg, wg; include correction to wg for snow
c     note: cls still not being used (latent heat correction for snow)
c
 
c.... hl = el (defined in inital.f) : define fcls = additional latent
c.... heat factor for sublimation (1.13344)
      fcls=1.0+(hlf/hl)
C
C     SET SURFACE PARAMETER VALUES
C
c      wgsat=0.36
      do 359 mg=1,ln2
      cie(mg)=0.0
  359 ci(mg)=0.0
      do 360 mg=1,ln2
c**** set logical surface type arrays
      land(mg)=imsl(mg,lg).ge.4
      sea(mg)=imsl(mg,lg).eq.3
      mlo(mg)=iabs(imsl(mg,lg)).eq.2
  360 cice(mg)=imsl(mg,lg).eq.1
c     For the qflux ocean set the sea points to be mlo
c Now for coupled runs as well set mlo points to be sea points
      if ( qflux) then
        do mg=1,ln2
          mlo(mg)=mlo(mg).or.sea(mg)
	  sea(mg) = .false.
        end do
      end if
      if ( lcouple) then
        do mg=1,ln2
          sea(mg)=mlo(mg).or.sea(mg)
	  mlo(mg) = .false.
        end do
      end if
c**** set up surface values for different surfaces (s,l,i,m)
c Use monthly albedo data set for R42 version (4/92)
c Now using monthly albedos for both R42 and R21 (5/92)

      do 361 mg=1,ln2
        als(mg)=alsf1(mg)*(1.0-ratlm)+alsf2(mg)*ratlm
        wetfac(mg)=1.0
        cls(mg)=1.0
 361  continue 

c set up temp at sea points
      if(lcouple)then
        do 3600 mg=1,ln2
 3600   tintp(mg)=tstar1(mg)*(ammratlm-1.0)-tstar2(mg)*ammratlm
      else
        do 362 mg=1,ln2
        if(.not.land(mg))
     &  tintp(mg)=tstar1(mg)*(ratlm-1.0)-tstar2(mg)*ratlm
  362   continue
      endif

c Calc. heat flux from ocean to ice assuming heat diffusion.

c These values now set in atstart.f
c     rhocw1=4.1e6
c     if(.not.semice)rhocw1=0.0
c     antarck(1)=0.15e-4
c     antarck(2)=0.15e-4
c     diz=25.

c Use IGW's qflux code with variable MLO depth          

      if (qflux) then
        do 3609 mg=1,ln2
        hoice(mg,lg)=0.0
        tratlm=1.-ratlm
        if(cice(mg)) hoice(mg,lg)=tratlm**2 * gm1cur(mg,lg)
     &       +tratlm*gm2cur(mg,lg) + gm3cur(mg,lg)
 3609   continue
      else
        do 3610 ns=1,2
        do 3610 mg=1,lon
        ma=mg+(ns-1)*lon
        hoice(ma,lg)=0.0
        if(cice(ma)) hoice(ma,lg)=
     &    rhocw1*diz*antarck(3-ns)*(tintp(ma)-tfi)/((0.5*diz)**2)
 3610   continue
      endif

c this is old code for Mark2 model initailization problems
      if(lcouple)then
        do 3612 mg=1,ln2
        if(cice(mg))then
          tintpx=tsst1(mg,lg)*(ratlm-1.0)-tsst2(mg,lg)*ratlm
          reset=max(tintp(mg)-tintpx-3.0,0.0)
          hoicex=rhocw1*diz*reset/(20.0*86400.0)
          hoice(mg,lg)=hoice(mg,lg)+hoicex
        endif
 3612   continue
C ..................................................... begin .. ach . 13/11/98

c   If the actual ocean model top level temperature (sstocean) is less
c   than freezing, then add an additional ice to ocean heat flux
c   proportional to the deficit of this temperature below freezing.
c   Use a heat transfer coefficient hfrz, typically 20-200 times
c   larger than antarck above. This simulates formation of sea ice from
c   the ocean as the water temperature drops below freezing, and thereby
c   prevents the ocean temperature from falling much below freezing
c   while conserving heat flux.
c   It also seems a numerically gentle way of supporting ice freezing
c   from very cold water (compared to the total instant freeze and reset
c   method).
c
c   Notes: 1. "sstocean" uses temperature from mean of ocean arrays T and TA
c             as is done above for tstar1 and tstar2
c now linked to  tintp(mg) than full sstocean from ocean model
c          2. A "hfrz" of 6.0e-4 is equivalent to a heat exchange rate of
c             400 W/m^2 per deg C. (This is found adequate to keep the
c             surface water temperature above about -2.1 in most
c             applications.)
c          3. An "hfrz" of 210.e-4 is equivalent to the heat exchange
c             associated with an instant reset of ocean model temp to
c             freezing.  An "hfrz" larger than this may lead to numerical
c             instability.
c still assumes diz is 25m but keep at this length scale as this assumption has
c been used  other hoice calculations above.

        hfrz       = 6.0e-4
        celciusfrz = -1.85
        do mg=1,ln2
        sstocean=tintp(mg)-273.15
          if((cice(mg)).and.(sstocean.lt.celciusfrz))then
           hoice(mg,lg) = hoice(mg,lg) + rhocw1*diz*hfrz*
     &         (sstocean-celciusfrz)/((0.5*diz)**2)
c          write(6,*)'hoice ',hoice(mg,lg),mg,lg
          endif
        enddo

C ........................................................................ end
      endif

c****
c**** Set up values for various surface types
c****

      do 372 mg=1,ln2
      if(.not.sea(mg))go to 372
C**   SEA POINTS
c    correct for temp at elevation
      tstar(mg)=tintp(mg)-0.0065*z4(mg)
      als(mg)=0.05/(cosz(mg)+0.15)
  372 continue

      do mg=1,ln2
      if(.not.land(mg))then
        wg(mg)=-7777777.0
        wg2(mg)=-7777777.0
      endif
      enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  Next for Mk3 nsib start up only
      if(header2(50:50).ne.'M')then
      do 3800 mg=1,ln2
      if(land(mg))wb(mg,1,lg)=wg(mg)
 3800 continue
      endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do 373 mg=1,ln2
C       SET VALUES FOR SNOW FREE LAND
C       WETFAC FROM SURFACE SOIL MOISTURE FACTOR
      if(land(mg).and.(snowd(mg).eq.0.))then 
        select case (lsm_type)
        case ("nsib ")
          wgsat=pssat(mg,lg)
        case default
          wgsat=0.36 ! ??
        end select
        wetfac(mg)=min(1. , wb(mg,1,lg)/wgsat)
        gamms(mg)=3.5e+05
      endif
  373 continue

c CABLE does not require these, but keep cls, wetfac and gamms in case:
      IF (lsm_type .eq. "cable") THEN
       where(land.and.(snowd.gt.0.))
        cls=fcls
        wetfac=wetfac*(0.008*tstar-1.18528)
        gamms=4.857e+04
       endwhere
      ELSE

      notsnalb = header2(44:48).ne.'SNALB'
      if(.not.(lsm_type .eq. "nsib "))notsnalb=.true.

      if ( notsnalb ) then

c      must be commented once snage and ssdn read from the input file
      do 3731 mg=1,ln2
      if(land(mg))then
        snage(mg,lg) = 0.0
        ssdnn(mg,lg) = 140.0
      endif
 3731 continue

      do 3732 mg=1,ln2
      if(land(mg).and.(snowd(mg).gt.0.))then
C**    LAND POINTS,SNOW COVERED. SNOWD IN MM,10CM SNOW=1CM OF WATER.
        snalb = 0.8
        if ( tstar(mg) .ge. 273.09 ) snalb = 0.6
        als(mg)=als(mg) + (snalb-als(mg))*sqrt(snowd(mg)*0.1)
        if(snowd(mg).gt.10.0) als(mg)=snalb
c
        cls(mg)=fcls
        wetfac(mg)=wetfac(mg)*(0.008*tstar(mg)-1.18528)
        gamms(mg)=4.857e+04
      endif
 3732 continue

      else ! notsnalb

      do 3733 mg=1,ln2
      if(land(mg).and.(snowd(mg).gt.0.))then
C**    LAND POINTS,SNOW COVERED. SNOWD IN MM,10CM SNOW=1CM OF WATER.
c
c      new snow albedo (needs osnowd from the previous dt)
c
       saveals = als(mg)
c????  if( saveals .ge. 0.6 ) saveals = 0.19
       if( snowd(mg).le. 1.0 ) snage(mg,lg) = 0.0      
c
c      Snow age is dependent on snow crystal growth, surface melting,
c      accumulation of dirt and amount of a new snow.
c
       ttbg = min(tggsl(mg,1,lg), 273.1)
       if( isflag(mg,lg) .eq. 1) ttbg = tggsn(mg,1,lg) 
       ar1 = 5000.0*( 1./273.1 - 1./ttbg) ! crystal growth
       ar2 = min(0.0,10.0*ar1)   ! surface melting
       ar3 = 0.3                 ! accumulation of dirt
CEAK - suggested change to next line
CEAK   if(isoilm(mg,lg).eq.9) ar3=0.008   
       if(isoilm(mg,lg).eq.9) ar3=0.002

       dtau=1.0e-6*( exp(ar1) + exp(ar2) + ar3 ) * dt
       dnsnow=0.
CEAK - suggested change to next 2 lines (i.e. no conditional test)
CEAK   if(isoilm(mg,lg).ne.9.and.snowd(mg).lt.120.)
CEAK & dnsnow=min(1.,0.1*max(0.,snowd(mg)-osnowd(mg,lg)))!new snow (cm H2O)
       dnsnow=min(1.,0.1*max(0.,snowd(mg)-osnowd(mg,lg)))!new snow (cm H2O)
       if(isoilm(mg,lg).eq.9) dnsnow=max(dnsnow,0.0013)
       snage(mg,lg)=max(0.,(snage(mg,lg) + dtau)*(1.-dnsnow))
       snr = snowd(mg)/max(ssdnn(mg,lg),100.)                      
       snrat=min(1.,snr/(snr+0.02))

c      spsnd=max(25.,snowd(mg))
c      if( he(mg).gt.1200. ) spsnd=snowd(mg)
c      cotp0= sqrt(spsnd/(spsnd+max(0.1,0.1*he(mg))))
c      tz4=max(0.0,z4(mg)-1000.)
c      if( z4(mg).gt.3000. ) tz4=2.*z4(mg)
c      cotp1= sqrt(spsnd/(spsnd+max(0.1,0.3*tz4)))
c      cotp=min(cotp0,cotp1)
c      if( z4(mg).gt.4000. ) cotp=min(0.05,cotp)
c      snrat=min(1.,snrat*cotp)

       if(isoilm(mg,lg).eq.9) snrat=min(1.,snr/(snr+0.001))
c
c      Snow albedo is dependent on zenith angle and  snow age.
       alvo = 0.95                             !alb. for vis. on a new snow
       aliro = 0.65                            !alb. for near-infr. " - "
       fage = 1.0 - 1./( 1.0 + snage(mg,lg)) !age factor

c      albedo zenith dependence
c       cs = 0.2, cn = 0.5, b = 2.0
c       fzen = ( 1.0 + 1.0/b ) / ( 1.0 + 2.*b*cczen ) - 1./b
c       alvd = alvo * (1.0-cs*fage)
c       alird = aliro * (1.-cn*fage)
       cczen = max(0.17365, cosz(mg))
       fzen = ( 1.0 + 1./2. ) / ( 1.0 + 2.*2.*cczen ) - 1./2.
       if( cczen .gt. 0.5 ) fzen = 0.0
       fzenm = max ( fzen, 0. )
       alvd = alvo * (1.0-0.2*fage)
       alv = 0.4 * fzenm * (1.0-alvd) + alvd

       alird = aliro * (1.-0.5*fage)
       alir = 0.4 * fzenm * (1.0-alird) + alird
       talb = 0.5 * ( alv + alir )             ! snow albedo 

       tsigmfx=0.0 
       if( sigmf(mg) .gt. 0.02) then
         tsoil=0.5*(0.4*tggsl(mg,2,lg)+0.6*tggsl(mg,3,lg)
     &         + 0.95*tggsl(mg,4,lg) +  0.05*tggsl(mg,5,lg))
         ftsoil=max(0.0,1.0-0.0016*(298.0-tsoil)**2)
         if( tsoil .ge. 298.0 ) ftsoil=1.0
         tsigmfx=max(0.0,sigmf(mg)-pscveg(mg,lg)*(1.0-ftsoil))
         cc=min(1.,snr/max(snr+2.*z0m(mg),0.02))
         tsigmfx=(1.-cc)*tsigmfx
       endif

       alss = (1.-snrat)*saveals+snrat * talb !canopy free surface albedo 
       als(mg)=min(0.8,(1.-tsigmfx)*alss + tsigmfx*saveals)
c
       cls(mg)=fcls
       wetfac(mg)=wetfac(mg)*(0.008*tstar(mg)-1.18528)
       gamms(mg)=4.857e+04
      endif
 3733 continue

      endif ! notsnalb

      ENDIF ! lsm_type .eq. "cable"

      do 374 mg=1,ln2
      if(.not.cice(mg))go to 374
C**    SEA ICE POINTS
      cls(mg)=fcls
      wetfac(mg)=wetfac(mg)*(0.008*tstar(mg)-1.18528)
c- Melting temp for ice (no snow) = 273.00
        if ( tstar(mg) .ge. 272.95 ) then
          als(mg) = 0.55
        else
          als(mg) = 0.65
        end if
        if(lw.eq.22)als(mg)=als(mg)+0.05 !Slightly higher sea ice albedos for R21 
      ci(mg)=2.04/(siced(mg)+0.001)
      cie(mg)=ci(mg)
      gamms(mg)=3.471e+05
      if(snowd(mg).le.0.0)go to 374
C     CONVERT SNOWD MM TO M
      bs(mg)=0.31e+03/(snowd(mg)+1.0)
      cie(mg)=bs(mg)*ci(mg)/(bs(mg)+ci(mg))
      gamms(mg)=4.857e+04
c- Melting temp for ice with snow = 273.15
        if ( tstar(mg) .ge. 273.10 ) then
          snalb = 0.7
        else
          snalb = 0.8
        end if
      als(mg)=snalb
  374 continue

      if(.not.leads)then
        do 212 mg=1,ln2
        pl(mg)=1.0
        if(cice(mg))pl(mg)=0.0
  212   continue
      end if
      do 213 mg=1,ln2
        if(cice(mg))then
          sicef(mg)=1.0-pl(mg)
        else
          sicef(mg)=0.
        endif
  213 continue

      if(.not.lcouple)then ! SPO change
      if(leads)then
c grey ice albedos, for water logged spring melt.
      do 375 mg=1,ln2
      if(.not.cice(mg))go to 375
      if(tstar(mg).lt.272.90)go to 375
      if(pl(mg).lt.0.5.or.dic(mg).gt.0.50)go to 375
      if(tni(mg).lt.272.65)go to 375
      if(tmix(mg).lt.271.70)go to 375
      if(dsn(mg).gt.0.10)go to 375
      als(mg)=0.3
  375 continue
      endif
      endif

      do 370 mg=1,ln2
      if(.not.mlo(mg))go to 370
C**    MIXED LAYER OCEAN POINTS
      als(mg)=0.05/(cosz(mg)+0.15)
      if(imsl(mg,lg).eq.-2)then !Point has just changed from sea to MLO
        tstar(mg)=tintp(mg)
        imsl(mg,lg)=2
      endif
c    correct for temp at elevation
      tstar(mg)=tstar(mg)-0.0065*z4(mg)
  370 continue
 
      do 371 mg=1,ln2
  371 tg(mg)=tstar(mg)

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a)')'After surfset.'
          write(25,1)'wb 1',wb(mg,1,lg),' wb ',wb(mg,ms,lg),
     &               ' tg ',tg(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,f7.2))

      if((mod(mins+int(mstep, 8),1440_8).eq.0_8).and.gwicm)then
c****
c**** create surface maps at end of day only

      al2=log(2.0)

      do 470 mg=1,ln2

      if(land(mg))then
c**    land points,possibly snow covered.
c**    snowd in mm,10cm snow=1cm of water.
         if(snowd(mg).gt.0.01)then
           snwice(mg,lg)='*'
           kx=max(1,min(12,int(2.0+log(snowd(mg))/al2)))
           snmap(mg,lg)=cint(kx)
         else
c     set values for snow free land
           snwice(mg,lg)='.'
           snmap(mg,lg)='.'
         end if
c      if(.not.nsib)then
      if(lsm_type .eq. "nsib ")then
        kx=1+min(11,nint(wb(mg,1,lg)/0.04))
        chwg(mg,lg)=cint(kx)
        kx=1+min(11,nint(wb(mg,ms,lg)/0.04))
        chwb(mg,lg)=cint(kx)
      endif

      end if

  470 continue

      do 471 mg=1,ln2

      snicel(mg,lg)=snwice(mg,lg)
      sifeet(mg)=0.0

      if(cice(mg))then
c**    sea ice points
      sicft=siced(mg)*3.2808
      kx=1+min(11,nint(sicft))
c     if(sicft.lt.0.5)kx=1
      snwice(mg,lg)=cint(kx)
c.... map of ice depth accounting for leads
      if(leads)then
        sicft=(1.0-pl(mg))*siced(mg)*3.2808
      else
        sicft=siced(mg)*3.2808
      endif
      sifeet(mg)=sicft
      kx=1+min(11,nint(sicft))
c     if(sicft.lt.0.5)kx=1
      snicel(mg,lg)=cint(kx)
        if(snowd(mg).gt.0.01)then
          kx=max(1,min(12,int(2.0+log(snowd(mg))/al2)))
          snmap(mg,lg)=cint(kx)
        else 
          snmap(mg,lg)='i'
        end if
      end if

  471 continue

      do 472 ns=1,2
      ma=(ns-1)*lon
      zisumx=0.0
      jicex=0
      do 472 mg=1+ma,lon+ma
      if(cice(mg))then
        jicex=jicex+1
        zisumx=zisumx+sifeet(mg)
      endif
      zisum(lg,ns)=zisumx
      jice(lg,ns)=jicex
  472 continue

      endif

      return
      end
      subroutine surfice(sgold,rgold,fg,eg,tg,coszro2,z4,
     &                   pwetfac,pcls,pals,ptintp,ptg)

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real sgold(ln2)
      real rgold(ln2)
      real fg(ln2)
      real eg(ln2)
      real tg(ln2)
      real coszro2(ln2)
      real z4(ln2)
      real pwetfac(ln2)
      real pcls(ln2)
      real pals(ln2)
      real ptintp(ln2)
      real ptg(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'

      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)

C Global data blocks

C Local work arrays and variables
      real surfb
      real tfact

      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c Set surface values for ice/leads points
 
c Split the leads into cooling regions where ice is close to maximum
c concentration (95% for both hemispheres) and slow new growth in open water
c by covering it with thin white covering albedo of 0.3 for grease ice nilas
c Set tempearature to be a combination of ocean and ice surface temperatures.

         do 270 mg=1,ln2
         pwetfac(mg)=1.0
         pcls(mg)=1.0
         surfb=sgold(mg)-rgold(mg)-fg(mg)-eg(mg)
         if(surfb.lt.0.and.pl(mg).lt.0.05)then
            pals(mg)=0.3
            tfact=0.3
            ptintp(mg)=tfact*tmix(mg)+(1-tfact)*tg(mg)
         else
c     mixed layer ocean points
            pals(mg)=0.05/(coszro2(mg)+0.15)
            ptintp(mg)=tmix(mg)
         endif
c     correct mlo t* for z4 factor
c still do for both cases as readjusted in surfupl
          ptg(mg)=ptintp(mg)-0.0065*z4(mg)
c set dummy land values for (vector) entire latitude row calculations
          if(land(mg))ptg(mg)=tg(mg)
 270     continue

      return
      end
