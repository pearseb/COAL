c nsib comparison is modified to lsm_type = "nsib " 
c AJA 2009/01/22
c
c Inserted code to detect floating point exceptions before they occur. If
c the square root of a negative number is about to be taken, the run is aborted
c and the values of key variables are written to standard output.
c Note that this is only done within the calculation of the screen temperature,
c as this is the only point at which I have ever encountered FPEs.
c SJP 2004/01/05
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: hsflux.f,v $
c Revision 1.73  2001/06/04 02:27:06  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.72  2001/02/22 05:34:46  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.71  2000/11/14 03:11:40  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.70  2000/08/16 02:59:26  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.69  2000/06/20 02:08:35  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.68  1999/05/20 06:23:58  rot032
c HBG changes to V5-2
c
c Revision 1.67  1998/12/10  01:08:39  ldr
c HBG changes to V5-1-21
c
c Revision 1.66  1998/01/30  04:36:15  ldr
c Fix from SPO for Cd diagnostic.
c
c Revision 1.65  1997/12/23  00:23:39  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.64  1997/12/17  23:23:08  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.63.2.1  1997/12/19  02:03:19  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.63  1996/10/24  01:02:51  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.62  1996/06/13  02:06:41  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.61  1996/04/17  01:45:13  ldr
c Tighten up moisture conservation and improve some related diagnostics.
c
c Revision 1.60  1996/03/21  03:18:47  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.59  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.58  1995/11/10  00:46:21  ldr
c Put back the bugfix to the 10m wind loop which caused crash at T63 - this
c was not included in revs 1.53-1.57.
c
c Revision 1.57  1995/08/31  04:30:42  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.56  1995/08/29  02:43:38  ldr
c Merge of HBG's corrections to V4-7-13h with LDR's changes for run g61.
c
c Revision 1.55  1995/08/18  06:58:36  ldr
c Merge of LDR's and HBG's changes to V4-7-12l.
c
c Revision 1.53.1.1  1995/08/18  06:14:46  ldr
c LDR's changes to V4-7-12l to bring cloud scheme to run g57.
c
c Revision 1.54  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.54.1.1  1995/08/29  01:30:09  ldr
c HBG's corrections to his V4-7-13h, to make results with hybrid=F agree
c with previous version.
c
c Revision 1.53  1995/08/14  07:32:40  ldr
c Updated qcloud to run g39, run for several years on kaos.
c
c Revision 1.52  1995/08/14  05:21:43  ldr
c Fixed a little error in the 10m wind loop that caused crash on kaos.
c
c Revision 1.51  1995/08/08  02:02:15  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.50  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.39.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.49  1995/06/30  02:44:42  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.48.1.1  1995/11/10  05:23:03  ldr
c EAK's changes to soil types (increase to 9 soil types).
c
c Revision 1.48  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.47  1994/12/22  23:51:06  ldr
c Mods to cloud scheme of V4-6-13l to produce faster version for HBG and
c to allow vertically subgrid cloud below inversions in BL. Set slowice=T
c for accurate (but slow) ice fall calculation. Also, Ri now includes hlfusionn
c
c Revision 1.46  1994/11/17  03:07:28  mrd
c Remove a possible divide by zero in 2m RH calculation for neutral conditions.
c
c Revision 1.45  1994/09/13  12:15:43  mrd
c Split 10m wind loop for better vectorization.
c
c Revision 1.44  94/09/13  11:42:06  mrd
c Add daily mean 10m windspeed and 2m relative humidity as history variables.
c 
c Revision 1.43  94/09/12  12:50:08  ldr
c Changes to qcloud scheme of V4-6 to create run E32, except that E32 has
c enhanced collection of cloud water by falling ice.
c 
c Revision 1.42  94/08/08  17:21:27  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.41  94/08/08  13:16:23  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.40  94/08/04  16:55:33  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.39  94/06/17  13:44:45  ldr
c Pass in T, q rather than Tliq, qtot for qcloud scheme.
c 
c Revision 1.38  94/05/13  14:55:45  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c 
c Revision 1.37  94/04/29  15:17:27  ldr
c Minor changes for 18L: start deep and shallow conv from k=3, increase Cd.
c 
c Revision 1.36  94/04/29  15:10:31  ldr
c Changes from HBG to reduce high cirrus, increase convective cirrus to
c previous value (2*cvx), to improve jet again. Also change low cloud and
c Cd over sea to achieve balance.
c 
c Revision 1.35  94/03/30  16:32:43  ldr
c Set min windspeed to 2m/s.
c 
c Revision 1.34  94/03/30  10:23:33  ldr
c Increase aftsea to .0014 and use JMcG's skin temperature.
c
c     INPUT/OUTPUT
c     Input:   from common/cnsta in CNSTA.f
c                  algf - log(sig(k))
c
c              from common/fewflags in FEWFLAGS.f
c                  debug    - flag to control column debugging
c                  lgdebug  - flag to control latitude debugging
c                  mgdebug  - flag to control longitude debugging
c                  insdebug - flag to control hemisphere debugging
c                  hybrid   - if true, use new hybrid vertical coordinate
c                  lsm_type - if "nsib ", use "New SIB" land surface scheme
c                  qcloud   - T if using prognostic cloud scheme
c
c              from common/hybrpr in HYBRPR.f
c                  pdpsk - (pressure/surface pressure)**cappa 
c                  prf   - pressure at full levels 
c
c              from common/levdata in this subroutine
c                  aftsea - neutral drag coefficient over sea
c
c              from common/logimsl in LOGIMSL.f
c                  logical variables for surface type:
c                  cice - seaice     mlo - mixed layer ocean
c                  land, sea
c
c              from common/timex in TIMEX.f
c                  dt - time step in seconds
c
c              from arguments
c                  ns - hemisphere index      ipass  - model pass counter
c                  lg  - latitude index        pg     - surface pressure(mbs) 
c                  tg  - surface temperature   qtg    - mixing ratio
c                  u   - zonal wind            v      - meridional wind
c                  snowd - snow depth          gam    -(L/cp)dqsdt
c                  ttg - physical temperature  wetfac - soil westness
c                  qfg - frozen H2O cloud/ice mixing ratio (kg/kg)
c                  qlg - liquid H2O cloud/ice mixing ratio (kg/kg)
c                  sg  - solar absorbed in ground 
c                  cfrac - cloudy fraction of grid box                 
c
c     Output:  from arguments
c                  degdw  - derivative of eg with respect to soil wetness
c                  taux   - zonal component of wind stress
c                  tauy   - meridional component of wind stress
c                  v10m   - 10m wind speed
c                  zmin   - height of lowest model level
c
c     In/Out:  from common/cdmapx in this subroutine
c                  cdavg - mean drag coefficient 
c                  pcdavg -  drag coefficient  over open water
c
c              from common/surf in this subroutine
c                  aftfhg - reciprocal of air resistance for bare ground
c                  aftfh  - reciprocal of air resistance for ground 
c                           covered with foliage
c
c              from arguments
c                  cduv  - momentum drag coefficient
c                  degdt - derivative of eg with respect to temp
c                  dfgdt - derivative of fg with respect to temp,
c                  fg    - sensible heat flux  eg - latent heat flux
c                  risav - saved value of Richardson number
c                  tscrn - temperature at screen level
c                  vmod  - magnitude of wind speed
c 
c
c If NSIB land surface scheme is used this routine evaluates:
c -  reciprocal of air resistance for bare ground (aftfhg)
c    and ground covered with foliage (aftfh).
c -  roughness length is decreased linearly with an increasing 
c    snow depth. For bare ground z0 (zobg=0.01m) decreases to 0.00024m when
c    snow depth >= 12cm,  for vegetation z0 (z0mx=z0m(mg)) decreases to 0.01m 
c    when snow depth >=10 times the  initial z0m(mg).
c
c    Outputs:
c    aftfhg - reciprocal of air resistance for bare ground
c    aftfh  - reciprocal of air resistance for ground covered with foliage
c

      subroutine hsflux(ipass,lg,wetfac,pg,tg,ttg,qtg,u,v,sg,snowd,!Inputs
     &                  qlg,qfg,cfrac,gam,
     &                  vmod,taux,tauy,cduv,degdt,dfgdt,eg,fg,     !Outputs
     &                  degdw,tscrn,v10m,risav,ustar,als,v10n,zmin)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /VEGDAT/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'CPARAMS.f'
      real zscreen,delta
      parameter (zscreen=2.0)
      parameter (delta=1/epsil-1)

C Argument list
      integer ipass,lg
      real wetfac(ln2)
      real pg(ln2)
      real tg(ln2)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real u(ln2,nl)
      real v(ln2,nl)
      real sg(ln2)
      real snowd(ln2)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real cfrac(ln2,nl)
      real gam(ln2,nl)
      real vmod(ln2)
      real taux(ln2)
      real tauy(ln2)
      real cduv(ln2)
      real degdt(ln2)
      real dfgdt(ln2)
      real eg(ln2)
      real fg(ln2)
      real degdw(ln2)
      real tscrn(ln2)
      real v10m(ln2)
      real v10n(ln2)
      real risav(ln2)
      real ustar(ln2)
      real als(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'LOGIMSL.f'

      real rsmin,z0m,sigmf,vegt
      common/vegdat/ rsmin(ln2),z0m(ln2),sigmf(ln2),vegt(ln2)

C Global data blocks
      include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'TIMEX.f'

      real cdavg,pcdavg
      common/cdmapx/cdavg(ln2,lat),pcdavg(ln2,lat)

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      real zosav(ln2),xsav(ln2),afsav(ln2),aftsav(ln2)
      real rhosav(ln2),rstssav(ln2),drstsav(ln2)
      real fhsav(ln2),fhbgsav(ln2) ! ,fmbgsav(ln2)
      real factch(ln2),tsurf(ln2)
      real zobg(ln2),zoland(ln2),aflandg(ln2),aftlandg(ln2)
      real zmin(ln2)
c     Temporaries for screen level calculations.
      real theta(ln2), tstarx(ln2), windsp(ln2)
      real afland(ln2),aftland(ln2)

      integer it
      integer mg
      integer ns

      real a
      real af
      real afroot
      real aft
      real al
      real alzz
      real b
      real betac
      real betaq
      real betaqt
      real betat
      real betatt
      real c
      real ch
      real cm
      real con
      real con1
      real conh
      real conw
      real consea
      real daf
      real dden
      real deltheta
      real deg
      real den
      real denha
      real denhabg
      real denma
      real dfm
      real drst
      real fh
      real fhbg
      real fm
      real hlrvap
      real qc
      real qtgair
      real qtgnet
      real qtot
      real qt1
      real rho
      real ri
      real rich
      real root
      real rootbg
      real rsts
      real taftfh
      real taftfhg
      real tauxy
      real tliq
      real x
      real y
      real zminavn
      real zminavs
      real zo
      real z1
      real z0mx
      real temp_sjp

C Local data, functions etc
      real zoice
c     data zoice/0.0001/
c see if other more major changes have compensated for this minor adjustment
c SPO 29/1/93
CHBG  data zoice/0.001/
      data zoice/0.01/ ! 9/12/97 HBG
      integer nevap
      real wgmax
      data nevap/1/,wgmax/0.36/
c     can set coefficients for Louis scheme based on Kansas data:
      real bprm,cms,chs,vkar
c     data bprm/4.7/,cms/7.4/,chs/5.3/,vkar/.35/,rdiv/.74/
c     or use the following improved coefficients with vkar=.4:
      data bprm/5./,cms/5./,chs/2.6/,vkar/.4/
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

c     Switches automatically to form  (eg - qairterm) *wg/wgmax for nevap=0
c     Evap has dt* so as to give diagnostic evap per timestep
c
c     N.B. cduv is returned as drag coeffs mult by vmod
c     ttg, u, v, qtg are current values
c     pg is surface pressure(mbs), tg is surface temperature
c     sig gives sigma values
c     fg is sensible heat flux (fg in SURFUPA/B)
c     eg is latent heat flux (eg in SURFUPA/B)

c aftsea has been preset in initax.f

      hlrvap=hl/rvap
c     cp specific heat at constant pressure joule/kgm/deg

      if(hybrid)then
        do mg=1,ln2
          zmin(mg)=-rdry*ttg(mg,1)*log(prf(mg,1)/pg(mg))/grav
        enddo
      else
c     zmin is the approx height of the lowest level in the model
        zminavn=-rdry*ttg(lon/2,1)*algf(1)/grav
        zminavs=-rdry*ttg(lon+lon/2,1)*algf(1)/grav
        do mg=1,lon
          zmin(mg)=zminavn
          zmin(mg+lon)=zminavs
        enddo
      endif
      do 10 mg=1,ln2
      zoland(mg)=0.1682 ! Old (non-NSIB) value
   10 zobg(mg)=0.01
c                   decrease roughness length if snow present
      if(lsm_type .eq. "nsib ")then
        do 12 mg=1,ln2
        if(land(mg))then
        z0mx=z0m(mg)
        if(snowd(mg).gt.0.0)then
          z0mx=max(z0m(mg)-(snowd(mg)/
     &    (10.0*z0m(mg)*100.0))*(z0m(mg)-0.01),0.01)
          zobg(mg)=max(zobg(mg)-(snowd(mg)/12.)*0.00976, 0.00024)
        end if
        zoland(mg)=z0mx
        end if
   12   continue
      end if

      do 14 mg=1,ln2
      alzz=log(zmin(mg)/zobg(mg))
      aflandg(mg)=(vkar/alzz)**2
      aftlandg(mg)=vkar**2/( alzz * (2.+alzz) )
   14 continue
      do 16 mg=1,ln2
      if(land(mg)) then
        alzz=log(zmin(mg)/zoland(mg))
        afland(mg)=(vkar/alzz)**2
        aftland(mg)=vkar**2/( alzz * (2.+alzz) )
      endif
   16 continue

c     precompute x, vmod, ri, rsts, drst

      do 2 mg=1,ln2
      windsp(mg)=sqrt(u(mg,1)**2+v(mg,1)**2)      
      vmod(mg)=max( windsp(mg), 2.)
c     if(sea(mg).or.mlo(mg))    
c    1vmod(mg)=max( sqrt(u(mg,1)**2+v(mg,1)**2) , 5.)

c     following changes represent skin temperature effect
      tsurf(mg)=tg(mg)
C      if(sea(mg).or.mlo(mg))then
C        tsurf(mg)=tsurf(mg)+.008*sg(mg)/(1.+.25*windsp(mg)**2)
Cc       tsurf(mg)=tsurf(mg)+.008*max((sg(mg)-200.0),0.0)
Cc    &               /(1.+.25*windsp(mg)**2)
C      endif

      rsts=qsat(100.0*pg(mg),tsurf(mg))
C Zeng et al, 1998, J Climate, 2628-2644 : 0.98 factor for qs over salt water
      if(sea(mg).or.mlo(mg).or.((ipass.eq.2).and.cice(mg)))
     &                                               rsts=0.98*rsts
      drst=rsts*hlrvap/tsurf(mg)**2
chbg  if(qcloud.and..not.land(mg))then
      if(qcloud)then
chbg    tliq=ttg(mg,1)+hlcp*qlg(mg,1)+hlscp*qfg(mg,1)
        tliq=ttg(mg,1)-hlcp*qlg(mg,1)-hlscp*qfg(mg,1)
        qtot=qtg(mg,1)+qlg(mg,1)+qfg(mg,1)
        qc=qlg(mg,1)+qfg(mg,1)
        betat=1/ttg(mg,1)
	betaq=delta/(1+delta*qtg(mg,1)-qc)
        al=1/(1+gam(mg,1))
        betac=cfrac(mg,1)* al * (hlcp*betat - betaq/(1-epsil) )
        betatt=(betat-(gam(mg,1)/hlcp)*betac) !Beta_t_tilde factor sigk removed
        betaqt=betaq+betac !Beta_q_tilde
        x=grav*zmin(mg) *
     &       (betatt*(tliq-tsurf(mg)*pdpsk(mg,1))+betaqt*(qtot-rsts))
      else
        x=grav*zmin(mg)*(1.-tsurf(mg)*pdpsk(mg,1)/ttg(mg,1))
c     To include virtual temperature term as well use
c       x=grav*zmin(mg)*(1.-tsurf(mg)*pdpsk(mg,1)/ttg(mg,1)
c    &             +.61*wetfac(mg)*(qtg(mg,1)-rsts))
      endif
      xsav(mg)=x
      risav(mg)=x*vmod(mg)**(-2)
      factch(mg)=1.
      rstssav(mg)=rsts
      drstsav(mg)=drst
   2  continue

c Invoke this code over sea-ice if on second pass (i.e. over leads)
      do 63 mg=1,ln2
      if(sea(mg).or.mlo(mg).or.(cice(mg).and.ipass.eq.2))then
        ri=risav(mg)
        x=xsav(mg)
c       This is in-line OCENZO using latest coefficient, i.e. .018
c       derivatives here are with respect to zo
        consea=vmod(mg)*.018/grav
c       consea=vmod(mg)*.032/grav
c    try 032 in hbc 12-4-91
        zo=.001
        if(x.gt.0.)then
          fm=vmod(mg) /(1.+bprm*ri)**2
          con=consea*fm
c         do 3 it=1,3
c         ---- do 3  has been manually unwound to aid the vectoriser
          afroot=vkar/log(zmin(mg)/zo)
          af=afroot**2
          daf=2.*af*afroot/(vkar*zo)
          zo=max(1.5e-5,zo-(zo-con*af)/(1.-con*daf))
c   3     continue
          afroot=vkar/log(zmin(mg)/zo)
          af=afroot**2
          daf=2.*af*afroot/(vkar*zo)
          zo=max(1.5e-5,zo-(zo-con*af)/(1.-con*daf))
          afroot=vkar/log(zmin(mg)/zo)
          af=afroot**2
          daf=2.*af*afroot/(vkar*zo)
          zo=max(1.5e-5,zo-(zo-con*af)/(1.-con*daf))
        else
c         ---- (do 6 does not need to be manually unwound)
          do 6 it=1,3
          afroot=vkar/log(zmin(mg)/zo)
          af=afroot**2
          daf=2.*af*afroot/(vkar*zo)
          con1=cms*2.*bprm*sqrt(-x*zmin(mg)/zo)
          den=vmod(mg)+af*con1
          dden=con1*(daf-.5*af/zo)
          fm=vmod(mg)-(2.*bprm *x)/den
          dfm=2.*bprm*x*dden*den**(-2)
          zo=max(1.5e-5,zo-(zo-consea*af*fm)/
     &                     (1.-consea*(daf*fm+af*dfm)))
   6      continue
        endif
          af=(vkar/log(zmin(mg)/zo))**2
	  afsav(mg)=af
c         following is suggested new af for sea points
c for haw run         aftsav(mg)=.0014
          aftsav(mg)=aftsea
c         aftsav(mg)=.00085
c         aftsav(mg)=.0014
	  zosav(mg)=zo
       endif    
   63 continue

c     Following are sea ice points
      do 65 mg=1,ln2
         if(cice(mg).and.ipass.eq.1)then
           zosav(mg)=zoice
           alzz=log(zmin(mg)/zoice)
           afsav(mg)=(vkar/alzz)**2
           aftsav(mg)=vkar**2/( alzz * (2.+alzz) )
         endif
   65 continue

c     Following are land or snow points
      do 67 mg=1,ln2
	 if(land(mg))then
           zosav(mg)=zoland(mg)
           afsav(mg)=afland(mg)
           aftsav(mg)=aftland(mg)
           factch(mg)=sqrt(7.4)
         endif
   67 continue

c     Having settled on zo (and thus af) now do actual fh and fm calculations
      do 8  mg=1,ln2
      x=xsav(mg)
      ri=risav(mg)
      zo=zosav(mg)
      af=afsav(mg)
      aft=aftsav(mg)
      
      if(x.gt.0.)then
        fm=vmod(mg) /(1.+bprm*ri)**2
        fh=fm
c        fmbg=fm
        fhbg=fh
      else
        root=sqrt(-x*zmin(mg)/zo)
c       First do momentum
        denma=vmod(mg)+cms*2.*bprm*af*root
        fm=vmod(mg)-(2.*bprm *x)/denma
c       n.b. fm denotes ustar**2/(vmod(mg)*af)
c       Now heat ; allow for smaller zo via aft and factch
        denha=vmod(mg)+chs*2.*bprm*factch(mg)*aft*root
        fh=vmod(mg)-(2.*bprm *x)/denha
c--  values for Nsib scheme
        rootbg=sqrt(-x*zmin(mg)/zobg(mg))
c        denmabg=vmod(mg)+cms*2.*bprm*aflandg(mg)*rootbg
c        fmbg=vmod(mg)-(2.*bprm *x)/denmabg
        denhabg=vmod(mg)+chs*2.*bprm*factch(mg)*aftlandg(mg)*rootbg
        fhbg=vmod(mg)-(2.*bprm *x)/denhabg
      endif
       fhsav(mg)=fh
c      fmbgsav(mg)=fmbg
       fhbgsav(mg)=fhbg

c     cduv is now drag coeff *vmod
      cduv(mg) =af*fm
c     cdtq(mg) =aft*fh
c     rho=100.0*pg(mg)*sig(1)/(rdry*ttg(mg,1))
      rho=100.0*pg(mg)/(rdry*tsurf(mg))
      rhosav(mg)=rho
      theta(mg)=ttg(mg,1)/pdpsk(mg,1)
c     Surface stresses taux, tauy: diagnostic only
      taux(mg)=rho*cduv(mg)*u(mg,1)
      tauy(mg)=rho*cduv(mg)*v(mg,1)
c     N.B. potential evaporation is now eg+eg2
      conh=rho*aft*cp
      conw=rho*aft*hl
      fg(mg)=conh*fh*(tsurf(mg)-theta(mg))
      dfgdt(mg)=conh*fh
      x=xsav(mg)
      rsts=rstssav(mg)
      drst=drstsav(mg)

      if(qcloud)then
        qt1=qtg(mg,1)+qlg(mg,1)+qfg(mg,1)
      else
        qt1=qtg(mg,1)
      endif
      if(nevap.eq.0) then
        eg(mg)= wetfac(mg)*conw*fh*(rsts-qt1)
        degdt(mg)= wetfac(mg)*conw*fh*drst
        degdw(mg)= conw*fh*(rsts-qt1)/wgmax
      else
c       following trick reduces -ve evap to 1/10th the value (assumes nevap=1)
        qtgnet=rsts*wetfac(mg) -qt1
        qtgair=rsts*wetfac(mg)-max(qtgnet,.1*qtgnet)
c       qtgair=qt1
c       egnet=eg(mg)*wg(mg)/wgmax  +eg2(mg)
        eg(mg)= conw*fh*(wetfac(mg)*rsts-qtgair)
        deg=fh*drst*conw*wetfac(mg)
c       following reduces degdt by factor of 10 for dew
        degdt(mg)=.55*deg+sign(.45*deg,qtgnet)
        degdw(mg)= conw*fh*rsts/wgmax
      endif
c     Surface drag mapping
      cdavg(mg,lg)=cdavg(mg,lg)+(2-ipass)*cduv(mg)/vmod(mg)
      pcdavg(mg,lg)=pcdavg(mg,lg)+(ipass-1)*cduv(mg)/vmod(mg)
    8 continue

      if(ipass.eq.1)then
       do mg=1,ln2
c Whitecaps is high surface stress regions
c giving max 5% extra surface albedo for large stress.
        if(sea(mg).or.mlo(mg))then
         tauxy=sqrt(taux(mg)**2+tauy(mg)**2) ! Stress magnitude
         als(mg)=als(mg)+0.25*max(0.0,min(0.2,tauxy-0.15))
        endif
       enddo
      endif

c                             compute reciprocal of the air resistance
      if(lsm_type .eq. "nsib ")then
        do 80 mg=1,ln2
        if(land(mg)) then
          aft=aftsav(mg)
          fh=fhsav(mg)
C***      if(nsteps.eq.0)then  !Uncomment to initialize new R42 run only (LDR)
C***        fhbg=fhbgsav(mg)
C***        aftfh(mg,lg) = aft*fh
C***        aftfhg(mg,lg) = aftlandg(mg)*fhbg 
C***      endif
          af=afsav(mg)
          taftfh = aft*fh
          fhbg=fhbgsav(mg)
          taftfhg = aftlandg(mg)*fhbg
          if(sg(mg).gt.0.0.and.sg(mg).lt.100.)then
            taftfh=max(0.5*aftfh(mg,lg),taftfh)
            taftfhg=max(0.5*aftfhg(mg,lg),taftfhg)
          end if
          if(sg(mg).ge.100.)then
            taftfh=min(2.*aftfh(mg,lg),taftfh)
            taftfhg=min(2.*aftfhg(mg,lg),taftfhg)
          end if
c                             prevent air resistance falling below half
c                             of the neutral conditions value
          aftfh(mg,lg)=max(taftfh,0.5*af*vmod(mg))
          aftfhg(mg,lg)=max(taftfhg,0.1*aflandg(mg)*vmod(mg))
        endif
 80     continue
      end if


c**** Calculate screen temperature
c---- Should make this part switchable (it is relatively expensive code)
      do 881 mg=1,ln2
        z1 = zscreen
        zo=zosav(mg)
        if(land(mg))then
          aftsav(mg)=vkar**2/( log(z1/zo) * (2.+log(z1/zo)) )
        else if(sea(mg).or.mlo(mg).or.(cice(mg).and.ipass.eq.2))then
c         Approximate factor to convert the level 1 aft to a 2m value.
          aftsav(mg)=2.*aftsav(mg)
        else if(cice(mg).and.ipass.eq.1)then
          aftsav(mg)=(vkar/log(z1/zo))**2
        end if
  881 continue
      do 88 mg=1,ln2
c     Calculate the fundamental scaling constants
      ustar(mg) = sqrt ( cduv(mg) * vmod(mg) )
      rho=rhosav(mg)
      tstarx(mg)= -fg(mg) / ( rho*cp*ustar(mg))
c     If the roughness length is greater than the screen height the
c     screen temperature is simply the surface temperature.
      if ( zosav(mg) .ge. zscreen ) then
        tscrn(mg) = tsurf(mg)
      else
        zo=zosav(mg)
        ri=risav(mg)
        z1 = zscreen
        afroot=vkar/log(z1/zo)
        af=afroot**2
        aft=aftsav(mg)
        c = -grav*z1/tsurf(mg) * tstarx(mg)/ustar(mg)**2 *
     &       af**1.5 / aft
        if ( tstarx(mg) .ge. 0. ) then 
c         Stable case
c         Solve quadratic for rich, form ax^2 + bx + c
c         because ac < 0 this has two real solutions, only one of which is 
c         positive and so appropriate for the stable case. 
c         There is a slight error here because we use tg rather than theta1
c         to calculate c. However this is a second order effect.
          a = bprm
	  b = 1.
          temp_sjp = b*b - 4.0*a*c
          if (temp_sjp .lt. 0.0) then
            write (*, *)
            write (*, *) "ABORTING: Fatal error in HSFLUX"
            write (*, *)
            write (*, *) "mg = ", mg
            write (*, *) "lg = ", lg
            write (*, *)
            write (*, *) "tsurf(mg)  = ", tsurf(mg)
            write (*, *) "tstarx(mg) = ", tstarx(mg)
            write (*, *)
            stop
          end if
          rich = ( -b + sqrt(temp_sjp) ) /(2.*a)
c         Use the positive solution
          deltheta = tstarx(mg) * afroot / aft * ( 1.+bprm*rich)
        else
c         Unstable case
c         Starting value of the Richardson number.
c         The iteration calculates fm and fh as above but using the Richardson
c         number directly rather than x
          y = sqrt(c)
          root=sqrt(z1/zo)
          cm=2.*bprm*cms*af*root
          ch=2.*bprm*chs*aft*factch(mg)*root
          fm=1.+2.*bprm*y*y/(1.+cm*y)
          fh=1.+2.*bprm*y*y/(1.+ch*y)
          y = sqrt(c) * fm**0.75 / sqrt(fh)
          fm=1.+2.*bprm*y*y/(1.+cm*y)
          fh=1.+2.*bprm*y*y/(1.+ch*y)
          deltheta=tstarx(mg)*sqrt(af*fm)/(aft*fh)
        end if
        tscrn(mg) = tsurf(mg) + deltheta
      end if

   88 continue

c  Calculate 10m winds, as for screen temperature
      do mg=1,ln2
         if(land(mg))then
            aftsav(mg)=vkar**2/( log(z1/zo) * (2.+log(z1/zo)) )
         else if(sea(mg).or.mlo(mg).or.(cice(mg).and.ipass.eq.2))then
c           Approximate factor to convert the level 1 aft to a 10m value.
c           Want 1.4 times the level 1 value. This was already multiplied by
c           2. for the temperature calculation.
c            aftsav(mg)=1.4*aftsav(mg)
            aftsav(mg)=0.7*aftsav(mg)
         else if(cice(mg).and.ipass.eq.1)then
            aftsav(mg)=(vkar/log(z1/zo))**2
         end if
      end do
      do mg=1,ln2
         zo=zosav(mg)
         ri=risav(mg)
         z1 = 10.
         afroot=vkar/log(z1/zo)
         af=afroot**2
         aft=aftsav(mg)
         c = -grav*z1/tsurf(mg) * tstarx(mg)/ustar(mg)**2 *
     &        af**1.5 / aft
         v10n(mg) = windsp(mg)/vmod(mg) * ustar(mg) / sqrt(af) !Neutral 10m wind for DMS flux
         if ( tstarx(mg) .ge. 0. ) then 
            a = bprm
            b = 1.
            rich = ( -b + sqrt(b*b-4.*a*c) ) /(2.*a)
c           Use the positive solution
            v10m(mg) = windsp(mg)/vmod(mg) * ustar(mg) *
     &                 ( 1.+bprm*rich) / afroot
         else
            y = sqrt(c)
            root=sqrt(z1/zo)
            cm=2.*bprm*cms*af*root
            ch=2.*bprm*chs*aft*factch(mg)*root
            fm=1.+2.*bprm*y*y/(1.+cm*y)
            fh=1.+2.*bprm*y*y/(1.+ch*y)
            y = sqrt(c) * fm**0.75 / sqrt(fh)
            fm=1.+2.*bprm*y*y/(1.+cm*y)
            v10m(mg) = windsp(mg)/vmod(mg) * ustar(mg) / sqrt(af*fm)
         end if
      end do

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,i1)')'After hsflux. IPASS = ',ipass
          write(25,1)'eg ',eg(mg),' ev(mm) ',eg(mg)*dt/hl,' fg ',fg(mg)
          write(25,1)'degdt ',degdt(mg),' dfgdt ',dfgdt(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,f9.3))

      return
      end
