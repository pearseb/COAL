c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c the rcondx is set for both nsib and cable
c AJA 2009/04/08
c
c downward shortwave sgdn is passed to surfupa for CABLE
c AJA 2009/03/13
c
c runoff argument is added for surfupa like in surfupb
c AJA 2009/02/04
c
c vmod: surface windspeed  is
c passed into surfupa.f
c AJA 2009/02/02
c
c nsib comparison is modified to lsm_type = "nsib "
c AJA 2009/01/22
c
c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c Removed unnecessary "include 'MACHINE.f'".
c SJP 2001/11/22
c
c $Log: radin.f,v $
c Revision 1.163  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.162  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.161.1.1  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.161  2001/06/04 02:26:57  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.160  2001/03/07 04:29:00  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.159  2001/02/28 04:36:39  rot032
c Further tidy ups from HBG
c
c Revision 1.158  2001/02/22 05:34:41  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.157  2001/02/12 05:39:52  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.156  2000/12/08 03:58:53  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.155  2000/11/21 01:05:37  rot032
c Change ompl to sicef for coupled_aero.
c
c Revision 1.154  2000/11/14 06:43:01  rot032
c Fixes to HBG code for kaos.
c
c Revision 1.152  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.151  2000/06/20 07:12:03  rot032
c Remove some f90 syntax for NEC.
c
c Revision 1.149  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.148  2000/06/20 02:27:10  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.147  2000/06/20 02:17:17  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.146  1999/06/30 05:42:06  rot032
c Changes for SCM.
c
c Revision 1.145  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.144  1998/12/10  04:05:31  ldr
c Minor fixes to HBG code.
c
c Revision 1.143  1998/12/10  00:55:43  ldr
c HBG changes to V5-1-21
c
c Revision 1.142  1997/12/23  00:23:37  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.141  1997/12/17  23:22:53  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.140  1997/11/24  23:25:28  ldr
c Indirect sulfates changes to give version used for GRL paper
c
c Revision 1.139  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.139.2.1  1997/12/19  02:08:07  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.138  1997/07/24  05:42:50  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.137  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.136  1997/06/11  02:21:30  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.135  1997/05/19  07:31:53  ldr
c Implement zmean cfrac diagnostic for dcloud scheme.
c
c Revision 1.134  1996/12/23  03:58:02  mrd
c Add new gravity wave drag scheme as an option.
c
c Revision 1.133  1996/12/19  23:40:50  ldr
c EAK's snow albedo changes merged and debugged.
c
c Revision 1.132  1996/11/22  03:41:41  ldr
c Comments and tidy-ups from LDR; results slightly changed by extra
c calculation of qsg in progcld.
c
c Revision 1.131  1996/10/24  01:03:12  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.130  1996/06/13  02:23:00  ldr
c Merge of TIE and LDR changes.
c
c Revision 1.128.1.1  1996/06/13  02:07:46  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.129  1996/06/13  01:51:03  ldr
c Update qcloud to run 24E (works for 24L and 18L)
c
c Revision 1.128  1996/04/17  02:19:12  ldr
c Replace the endif accidentally deleted after call conv.
c
c Revision 1.127  1996/04/17  01:45:13  ldr
c Tighten up moisture conservation and improve some related diagnostics.
c
c Revision 1.126  1996/03/25  04:31:29  ldr
c Replace clhx by hlcp in conv, rainda and radin (needed for V5-0-7).
c
c Revision 1.125  1996/03/21  03:19:01  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.124  1996/02/19  04:09:59  ldr
c Generalize for 24 levels.
c
c Revision 1.123  1996/01/09  06:18:16  ldr
c This is the version used for the long run H06 and write-up.
c
c Revision 1.122  1995/11/30  04:40:16  ldr
c Final updates to qcloud scheme to bring it to long run h63 used for
c BMRC workshop and write-up.
c
c Revision 1.121  1995/11/23  06:09:40  ldr
c Merge of LDR and EAK changes since V4-7-21l.
c
c Revision 1.120  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.119.1.1  1995/11/23  06:03:31  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.119  1995/08/29  01:57:14  ldr
c Merge of HBG's corrections to V4-7-13h with LDR's changes to bring qcloud
c to run g61.
c
c Revision 1.118  1995/08/18  06:58:36  ldr
c Merge of LDR's and HBG's changes to V4-7-12l.
c
c Revision 1.116.1.2  1995/08/29  01:50:06  ldr
c Update qcloud to run g61.
c
c Revision 1.116.1.1  1995/08/18  06:14:46  ldr
c LDR's changes to V4-7-12l to bring cloud scheme to run g57.
c
c Revision 1.117  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.116  1995/08/08  02:02:18  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.115  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.96.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.114  1995/06/30  02:44:43  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.113.1.1  1995/11/10  05:23:03  ldr
c EAK's changes to soil types (increase to 9 soil types).
c
c Revision 1.113  1995/05/04  04:46:27  ldr
c Just a bit of tidying up.
c
c Revision 1.112  1995/05/02  06:10:36  ldr
c Changes to V4-6-16l from LDR - mainly affecting qcloud scheme.
c This version used for run F35.
c
c Revision 1.111  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.110  1995/02/27  01:40:34  ldr
c Bugfixes received from HBG in response to problems on Fujitsu.
c
c Revision 1.109.1.1  1995/02/27  01:47:08  ldr
c LDR's changes to radin of V4-6-15l - working towards conservative treatment
c of snow in qcloud version. This to be merged with HBG revision 1.110.
c
c Revision 1.109  1994/12/22  23:51:06  ldr
c Mods to cloud scheme of V4-6-13l to produce faster version for HBG and
c to allow vertically subgrid cloud below inversions in BL. Set slowice=T
c for accurate (but slow) ice fall calculation. Also, Ri now includes hlfusionn
c
c Revision 1.108  1994/12/15  06:24:42  ldr
c Further mods to LDR's cloud scheme - this version is as used in E58 run
c and as described in seminar of 27/10/94.
c
c Revision 1.107  94/09/13  09:58:46  mrd
c Add daily mean 10m windspeed and 2m relative humidity as history variables.
c 
c Revision 1.106  94/09/12  16:31:22  ldr
c A merge of MRD's and LDR's changes since V4-6.
c 
c Revision 1.105  94/09/12  10:57:13  mrd
c Added daily average cloud, surface temperature and screen temperature
c as history variables.
c 
c Revision 1.104  94/09/12  10:56:10  mrd
c Added possbility of saving two history files at different frequencies.
c 
c Revision 1.103.2.1  94/09/12  12:50:20  ldr
c Changes to qcloud scheme of V4-6 to create run E32, except that E32 has
c enhanced collection of cloud water by falling ice.
c 
c Revision 1.103  94/08/10  11:40:03  ldr
c Corrections to zonal mean cloud prints.
c 
c Revision 1.102  94/08/09  11:20:44  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.101  94/08/08  16:14:34  ldr
c Move extra rhg calculation from radin into newrain where it belongs.
c 
c Revision 1.100  94/08/08  13:16:25  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.99  94/08/05  16:43:58  ldr
c Remove the spurious argument in call conv - was introduced after V4530
c but is no good with HBG's new conv.
c 
c Revision 1.98  94/08/04  17:07:03  ldr
c Merge of LDR and HBG changes.
c 
c Revision 1.97  94/08/04  10:32:53  ldr
c Split zonal mean cloud water diagnostic into separate liquid and frozen
c components. Also, omit cloud calculation over leads.
c 
c Revision 1.96.1.1  94/08/04  16:56:23  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.96  94/06/17  15:02:12  mrd
c Moved Hal's extrapolation of the top level temperature from radin to radfs
c to make the radiation code more modular.
c 
c Revision 1.95  94/06/17  14:52:47  ldr
c Merge of 1.92.1.1 with 1.94.
c 
      subroutine radin(lg,tdt,exptyp)

      implicit none

!$OMP THREADPRIVATE ( /ENERX/ )
!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /GIANT3/ )
!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /PHYSICAL/ )
!$OMP THREADPRIVATE ( /RELATE/ )
!$OMP THREADPRIVATE ( /RVARSG/ )
!$OMP THREADPRIVATE ( /STATS/ )
!$OMP THREADPRIVATE ( /WORK1/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f' !Coupled aerosol treatment (ECHAM) parameters
      include 'ECTRAC.f' !Coupled aerosol treatment (ECHAM) arrays

C Argument list
      integer lg
      real tdt
      character*50 exptyp 

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real uzon
      common/enerx/uzon(nl,2)

      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)
      real datice(ln2,14)
      equivalence (datice,dsn)

      include 'GIANT3.f'
      include 'HYBRPR.f'
      include 'LOGIMSL.f'
      include 'PHYSICAL.f'
      include 'RELATE.f'
      include 'RVARSG.f'

      real hevap,hrn,hflux,hcld,hcll,hclm,hclh
     & ,hrg,hrt,hsout,hsg,htst,htscrn,hwfg,hwfb,hsnd,hsid
     & ,hrtclr,hsoutclr,htb2,htb3,hvmod,htaux,htauy,hrnc
     & ,hrnf,hintrc,hals,htgg,htgf,hpev
     & ,hsev,hrgclr,hsgclr,hrgdn,hsgdn
     & ,hperc,hclc,hdtm,hqlp,hpwc,hrevap
     & ,hssubl,hpreci,hreffl,hcliq,hdms
     & ,hso2,hso4,hso2dd,hso4dd,hso2wd
     & ,hso4wd,hsem,hso2oh,hso2h2,hso2o3,hqfp,hdms1,hso21,hso41
     & ,hdmsoh,hdmsn3,hscwd
      common/stats/hevap(ln2),hrn(ln2),hflux(ln2)
     & ,hcld(ln2),hcll(ln2),hclm(ln2),hclh(ln2)
     & ,hrg(ln2),hrt(ln2)
     & ,hsout(ln2),hsg(ln2)
     & ,htst(ln2),htscrn(ln2)
     & ,hwfg(ln2),hwfb(ln2),hsnd(ln2),hsid(ln2)
     & ,hrtclr(ln2),hsoutclr(ln2),htb2(ln2),htb3(ln2)
     & ,hvmod(ln2),htaux(ln2),htauy(ln2)
     & ,hrnc(ln2)
     & ,hrnf(ln2),hintrc(ln2),hals(ln2),htgg(ln2),htgf(ln2),hpev(ln2)
     & ,hsev(ln2),hrgclr(ln2),hsgclr(ln2),hrgdn(ln2),hsgdn(ln2)
     & ,hperc(ln2),hclc(ln2),hdtm(ln2),hqlp(ln2),hpwc(ln2),hrevap(ln2)
     & ,hssubl(ln2),hpreci(ln2),hreffl(ln2),hcliq(ln2),hdms(ln2)
     & ,hso2(ln2),hso4(ln2),hso2dd(ln2),hso4dd(ln2),hso2wd(ln2)
     & ,hso4wd(ln2),hsem(ln2),hso2oh(ln2),hso2h2(ln2),hso2o3(ln2)
     & ,hqfp(ln2),hdms1(ln2),hso21(ln2),hso41(ln2),hdmsoh(ln2)
     & ,hdmsn3(ln2),hscwd(ln2)

      include 'WORK1.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'GAUSL.f'
      include 'GWDDATA.f'
      include 'HIST.f'
      include 'NEST.f'
      include 'PRINTT.f'
      include 'RADAV.f'
      include 'SURF1.f'
      include 'TIMEX.f'
      include 'TRACEBLK.f'
      include 'ZMEAN.f'

      real drunof,dsis,dsisb,dsfw,athf,athfa,atsi
      common/aforce/drunof(ln2,lat),dsis(ln2,lat),dsisb(ln2,lat)
     & ,dsfw(ln2,lat),athf(ln2,lat),athfa(ln2,lat),atsi(ln2,lat)
! rjm
       include 'gas.h'

      integer idfrz
      common/icedvz/idfrz(lat,2)

      real hlat_a,clat_a
      common/latheat/hlat_a(ln2,lat,nl),clat_a(ln2,lat,nl)

      real statsice
      common/masiv5/statsice(ln2,14,lat)

      real tgmax,tgmin,tsmax,tsmin,htsmax,htsmin
     & ,extsmax,extsmin,tggmax,tggmin,tgfmax,tgfmin
     & ,htggmax,htggmin,htgfmax,htgfmin,hpmc
      common/maxmin/tgmax(ln2,lat),tgmin(ln2,lat),
     &              tsmax(ln2,lat),tsmin(ln2,lat),
     &              htsmax(ln2,lat),htsmin(ln2,lat),
     &              extsmax(ln2,lat),extsmin(ln2,lat),
     &              tggmax(ln2,lat),tggmin(ln2,lat),
     &              tgfmax(ln2,lat),tgfmin(ln2,lat),
     &              htggmax(ln2,lat),htggmin(ln2,lat),
     &              htgfmax(ln2,lat),htgfmin(ln2,lat)
     &             ,hpmc(ln2,lat)

      real qbneg,qmean 
      common/negqb/qbneg(lat,2,nl),qmean(lat,2,nl)

      real devap,drain,dvgg
      common/prtar/devap(ln2,lat),drain(ln2,lat),dvgg(ln2,lat)

      real sshf,slhf,sswr,slwr
      common/prtar1/sshf(ln2,lat),slhf(ln2,lat),sswr(ln2,lat),
     &slwr(ln2,lat)

      real alat,ha,coszm,taudam,fjd,r1,dlt,slag,dhr
      common/radsav/alat(lat2),ha(lat2),coszm(lat2),taudam(lat2),
     &      fjd,r1,dlt,slag,dhr

      real brain
      common/rblock/brain(ln2,lat)

      integer insdiag,mgdiag,lgdiag,nsdiag
      common/sdiag/insdiag(15),mgdiag(15),lgdiag(15),nsdiag

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

      real zrad
      common/zradc/zrad(lat,2)

C Local work arrays and variables
C Local work arrays with corresponding "leads" arrays {often p---()}
C (Some relate to common blocks as shown below)
c
c-    /giant3/ cl(ln2), cll(ln2), clm(ln2), clh(ln2)
      real    pcl(ln2),pcll(ln2),pclm(ln2),pclh(ln2)

c-    /rvarsg/ pg(ln2), ttg(ln2,nl), tscrn(ln2), als(ln2)
      real             pttg(ln2,nl),ptscrn(ln2),pals(ln2)

      real  u(ln2,nl), v(ln2,nl)
      real pu(ln2,nl),pv(ln2,nl)

      real  sg(ln2), sint(ln2), sout(ln2)
      real psg(ln2),psint(ln2),psout(ln2)

      real  precs(ln2), precc(ln2), tintp(ln2)
      real pprecs(ln2),pprecc(ln2),ptintp(ln2)

      real  cduv(ln2), degdt(ln2), dfgdt(ln2), cls(ln2)
      real pcduv(ln2),pdegdt(ln2),pdfgdt(ln2),pcls(ln2)

      real  wetfac(ln2), eg(ln2), evap(ln2), fg(ln2), condx(ln2)
      real pwetfac(ln2),peg(ln2),pevap(ln2),pfg(ln2),pcondx(ln2)

      real  htk(ln2,nl), degdw(ln2)
      real phtk(ln2,nl),pdegdw(ln2)

      real  tg(ln2), qtg(ln2,nl), qsg(ln2,nl), rhg(ln2,nl), gam(ln2,nl)
      real ptg(ln2),pqtg(ln2,nl),pqsg(ln2,nl),prhg(ln2,nl),pgam(ln2,nl)

      real  rt(ln2), rg(ln2), blg(ln2)
      real prt(ln2),prg(ln2),pblg(ln2)

      real  vmod(ln2), taux(ln2), tauy(ln2)
      real pvmod(ln2),ptaux(ln2),ptauy(ln2)

      real  uten(ln2,nl), vten(ln2,nl), dqgdt(ln2,nl)
      real puten(ln2,nl),pvten(ln2,nl),pdqgdt(ln2,nl)

      real  dkevm(ln2,nl),dkegw(ln2,nl)
      real pdkevm(ln2,nl)

      real  sga(ln2), sgold(ln2), rgold(ln2)
      real psga(ln2),psgold(ln2),prgold(ln2)

      real  rtclr(ln2), soutclr(ln2), rgclr(ln2), sgclr(ln2)
      real prtclr(ln2),psoutclr(ln2),prgclr(ln2),psgclr(ln2)

      real  rgdn(ln2), sgdn(ln2)
      real prgdn(ln2),psgdn(ln2)

      real  qlg(ln2,nl), qfg(ln2,nl), ccov(ln2,nl), cldliq(ln2)
      real pqlg(ln2,nl),pqfg(ln2,nl),pccov(ln2,nl),pcldliq(ln2)

      real  qevc(ln2,nl), refflm(ln2)
      real pqevc(ln2,nl),prefflm(ln2)

      real  preci(ln2),rpreci(ln2) !Tau and Tau-1 values of frozen precip
      real ppreci(ln2)

      real  rhscrn(ln2), v10m(ln2)
      real prhscrn(ln2),pv10m(ln2)
      real v10n(ln2)              !10m wind corrected to neutral stability

c-    Convection
      real  fmp(ln2,nl)
      real pfmp(ln2,nl)
      real fmu(ln2,nl)            ! updraught flux (s**(-1)) +ve
      real fmd(ln2,nl)            ! downdraught flux (s**(-1)) +ve
      integer  kta(ln2), kba(ln2)
      integer pkta(ln2),pkba(ln2)

c-    Atmospheric latent heating and cloud
      real  hlat(ln2,nl), clat(ln2,nl)
      real phlat(ln2,nl),pclat(ln2,nl)

c-    Tracers
      real  cten(ln2,nl,ntrace)
      real pcten(ln2,nl,ntrace)

c-    Coupled aerosol treatment

      real xtg(ln2,nl,ntrac)   !Tracer mixing ratio [kg/kg]
      real xtgsav(ln2,nl,ntrac)!Tracer mixing ratio saved before chemistry [kg/kg]
      real pxtg(ln2,nl,ntrac)  !Tracer mixing ratio over leads [kg/kg] (not used)
      real xso4(ln2,nl)        !SO4 mixing ratio (=xtg(:,:,3)) [kg/kg]
      real xtu(ln2,nl,ntrac)   !Tracer mixing ratio in convective updraught [kg/kg]
      real fscav(ln2,nl,ntrac) !Convective tracer scavenging fraction (b/w 0 and 1)
      real xtem(ln2,ntrac)     !Sfc. flux of tracer passed to vertical 
                               !                          mixing [kg/m2/s]
      real mcmax(ln2)          !Maximum skin reservoir depth passed back
                               !                          from surfa
      real conwd(ln2,ntrac)    !Diagnostic only: SO2 wet deposition
      real so2wd(ln2)          !Diagnostic only: SO2 wet deposition
      real so4wd(ln2)          !Diagnostic only: SO4 wet deposition
      real so2oh(ln2)          !Diagnostic only: SO2 oxidation by OH
      real so2h2(ln2)          !Diagnostic only: Aqueous SO2 oxidation by H202
      real so2o3(ln2)          !Diagnostic only: Aqueous SO2 oxidation by O3
      real dmsoh(ln2)          !Diagnostic only: DMS oxidation by OH
      real dmsn3(ln2)          !Diagnostic only: DMS oxidation by NO3
      real so2dd(ln2)          !Diagnostic only: SO2 dry deposition
      real so4dd(ln2)          !Diagnostic only: SO4 dry deposition
      real sem(ln2)            !Diagnostic only: S emission

c-    Local work arrays without corresponding leads arrays
c-     Radiation
      real coszro2(ln2),taudar2(ln2)
c-     Gravity Wave Drag
      real bvnf(ln2,nl)
c-     Ice model
      real ompl(ln2),qvent(ln2)
      real sicef(ln2) !Sea-ice fraction
      integer il(ln2)
c-     Surface schemes
      real cie(ln2),ci(ln2),bs(ln2),gamms(ln2)
      real wgp(ln2),runoff(ln2)
      real tgx(ln2),ustar(ln2)
      real rich(ln2,nl),risurf(ln2)
      real perc(ln2)    ! Deep soil moisture percolation
      real rcondx(ln2),totpev(ln2),scalev(ln2),hfrm(ln2)
      real sublmi(ln2),snowmi(ln2),plold(ln2),sicold(ln2)
      real tddd(ln2)
      real surfb(ln2),flx(ln2)     !Coupled model
c-     Prognostic cloud
      real qlgsav(ln2,nl),qfgsav(ln2,nl),cfsav(ln2,nl)
      real clcon(ln2,nl),qccon(ln2,nl),cfrac(ln2,nl)
      real cldcon(ln2)
      real qevap(ln2,nl),qsubl(ln2,nl),cdso4(ln2,nl)
      real qlpath(ln2),fluxc(ln2,nl)
      real qfpath(ln2)
      real CCW(ln2,nl)
      real rainx(ln2)
      integer kbase(ln2),ktop(ln2) !Shallow conv.
c-     Work arrays
c     real tmass(ln2)
      real frht(ln2,nl)
      integer ifroz(2)

      logical radstep
      logical solarfit 

      integer ilatn
      integer ilats
      integer ineg
      integer ipass
      integer k
      integer lgn
      integer lil
      integer ma
      integer mg
      integer ns
      integer nt

      real csqr
      real c100g
      real dkedif
      real edifz
      real ezon
      real plx
      real qmk
      real tzon
      real uzox
      real zalb
      real zcfrac
      real zcl
      real zclc
      real zclfor
      real zclh
      real zcll
      real zclm
      real zcolwv
      real zeg
      real zev
c     real zfrac
      real zhf
      real zhtt
      real zsehf
      real zqlpath
      real zrg
      real zrh
      real zrn
      real zrnc
      real zrns
      real zrt
      real zrtclr
      real zsbal
      real zsbl
      real zsbl_l,zsbl_i,zsbl_s
      real zsg
      real zsid
      real zsin
      real zsnd
      real zsou
      real zsouclr
      real zts
      real ztscrn
      real ztss
      real z_tmax
      real z_tmin

      real zmin(ln2)

C Local data, functions etc
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

c**** Set pressures on hybrid vertical levels
      call radhyb

C**
C**   UNITS ARE KGM/M**2 , JOULES/M**2/SEC ETC
C**     BUT P* (SURF PRESS) IS IN MBS=100 KG/M/SEC**2
C**
      if(nrad.eq.0)then
        radstep=.false.
      else
        radstep = mod(mins/int(mstep, 8),int(nrad, 8)).eq.0_8
      endif
      solarfit=.true.
      csqr=sqrt(acsq(lg)*eradsq)

      if (semice) then
c**** Copy sea-ice model variables into work array
        do 10 k=1,14
        do 10 mg=1,ln2
          datice(mg,k)=statsice(mg,k,lg)
 10     continue
      end if

c**** First pass (ipass=1) means no leads
      ipass=1

c**** PRESET ALL UNSPECIFIED VARIABLES TO ZERO
      do 20 mg=1,ln2
         cldcon(mg)=0.0 ! Not set if qcloud=F
         cl(mg)=0.0     ! Not set if not radstep
         clh(mg)=0.0    !  "
         cll(mg)=0.0    !  "
         clm(mg)=0.0    !  "
         hfrm(mg)=0.0
         precc(mg)=0.0
         preci(mg)=0.0
         precs(mg)=0.0
         qvent(mg)=0.0
         rainx(mg)=0.0
         runoff(mg)=0.0
         sg(mg)=0.0
         sint(mg)=0.0
         snowmi(mg)=0.0
         sout(mg)=0.0
         sublmi(mg)=0.0
         tintp(mg)=0.0
 20      continue
      do 30 k=1,nl
         do 30 mg=1,ln2
         clat(mg,k)=0.0
         cfrac(mg,k)=0.0
 30      dkevm(mg,k)=0.0

c**** TEMP AND MOISTURE VARIABLES
      call radvars(csqr,u,v,qtg,dqgdt,ipass,pttg,exptyp)

c**** Copy cloud model variables from big arrays into working arrays
      call setqcld(lg,ipass,qlg,qfg,ccov,gam,pgam)
      
      if(coupled_aero)then
c**** Set up aerosol array xtg from bigxtg + optional diagnostics
        fscav(:,:,:)=0.
        conwd(:,:)=0.
        call c_aero1(lg,tdt,qtg,qlg,qfg,u,v, !Inputs
     &               mcmax,xtg)              !Outputs
      else
        xtg(:,:,:)=0.
      endif

c**** Calculate zenith angle for surface albedos and for the solarfit
c**    calculation.
      if (nrad.ne.0) then
         ilatn=lg
         ilats=lat2p-lg
         call zenith(fjd,r1,dlt,slag,alat(ilatn),ha(ilatn),
     &        float(mstep)/60.,lon,coszro2(1),taudar2(1))
         call zenith(fjd,r1,dlt,slag,alat(ilats),ha(ilats),
     &        float(mstep)/60.,lon,coszro2(1+lon),taudar2(1+lon))
      end if

!      WRITE(43,*) 'Before surfset, lg = ', lg
!      WRITE(43,'(200f5.2)') als
c**** SET UP SURFACE VALUES :LAND,SEA,MLO,CICE,ALBEDOES ETC.
c     CALL TIMER('SURFSET ',3)
      call surfset(lg,z4,tstar1,tstar2,alsf1,alsf2,
     &             wg2,he,snowd,siced,coszro2,tstar,wg,tintp,
     &             tg,als,cls,wetfac,cie,ci,bs,gamms,sicef)
c     CALL TIMER('SURFSET ',4)
!      WRITE(43,*) 'After surfset, lg = ', lg
!      WRITE(43,'(200f5.2)') als

      do mg=1,ln2
         pg(mg)=pn(mg)
      enddo

c**** Initialize surface radiation fields sg and rg
c**    from saved data in /masiv3/
      call radsgrg(ipass,lg,solarfit,nrad,coszro2,taudar2,tg,
     &  sg,rg,blg)

cAJA      if (lsm_type .eq. "nsib ") then
c**** precip from tau-1 used for nsib scheme
         do 130 mg=1,ln2
            rpreci(mg)=opreci(mg,lg) !Ice precip from tau-1
 130        rcondx(mg)=ocondx(mg,lg)
cAJA      end if

c**** SURFACE FLUXES (HEAT,EVAP) FROM HSFLUX
c**     SENSIBLE HEAT FLUX FG (JOULES/M**2/SEC)
c**     EVAPORATION EG (JOULES/M**2/SEC)
c     CALL TIMER('HSFLUX  ',3)
      call hsflux(ipass,lg,wetfac,pg,tg,ttg,qtg,u,v,sg,snowd,    !Inputs
     &            qlg,qfg,ccov,gam,
     &            vmod,taux,tauy,cduv,degdt,dfgdt,eg,fg,
     &            degdw,tscrn,v10m,risurf,ustar,als,v10n,zmin)   !Outputs
c     CALL TIMER('HSFLUX  ',4)
       wind10(:,lg)=v10n
*      print*,'rjmv1',lg,wind10(:,lg)
!      WRITE(43,*) 'After hsflux, lg = ', lg
!      WRITE(43,'(200f5.2)') als

!      IF (lg .eq. 23) WRITE(45,*) 'wb,wbice = ',wb(4,1,lg),wbice(4,1,lg)
c**** FIRST PART OF UPDATE OF TG AND SURFACE FLUXES
!      WRITE(46,*) 'Before surfupa, lg = ', lg
c     CALL TIMER('SURFUPA ',3)
      call surfupa(lg,degdt,dfgdt,degdw,cie,gamms,sg,    !Inputs
     &             tintp,rg,blg,z4,wg,wg2,snowd,rcondx,  !  "
     &             pg,ttg,als,qtg,cls,vmod,zmin,rpreci,  !  "
     &             eg,fg,tg,tb2,tb3,                     !In and Out 
     &    wgp,tddd,totpev,scalev,hfrm,he,mcmax,runoff,perc,tscrn)   !Outputs
                                                    ! tscrn is overwritten by CABLE values
c     CALL TIMER('SURFUPA ',4)
!      WRITE(46,*) 'lg, perc = ', lg, perc
!      WRITE(46,*) 'lg,hperc = ', lg,hperc
!      WRITE(46,*) 'After surfupa, lg = ', lg
!      WRITE(43,'(200f5.2)') als
!      write(38,'(a21,i3,128f9.3)') 'lg, snowd after upa: ', lg, snowd

c**** Compute RH at screen level
      call radrhs(pg,eg,fg,ttg,qtg,tscrn,rhscrn)

c**** Call to aerosol subroutines
      if(coupled_aero)then
        call c_aero2(lg,tdt,mcmax,tintp,v10n,sicef,ustar,eg,fg,rcondx, !Inputs
     &               xtg,                        !In & Out
     &               xtem, so2dd, so4dd, sem)    !Outputs
      endif

c**** VERTICAL MIXING
c     CALL TIMER('HVERTMX ',3)
      call hvertmx(ipass,lg,tdt,cduv,fg,eg,pg,tg,
     &             qlg,qfg,ccov,gam,tscrn,ustar,xtem,    !Inputs
     &             qtg,dqgdt,u,v,ttg,uten,vten,xtg,      !In and Out
     &             dkevm,cten,rich,kbase,ktop)           !Outputs
c     CALL TIMER('HVERTMX ',4)

c**** Call to aerosol decay routines
      if(coupled_aero)call c_aero3(lg,tdt,  !Inputs
     &                             xtg)     !In and out

c**** GRAVITY WAVE DRAG
c     CALL TIMER('GWDRAG  ',3)
      if ( ngwd .eq. 1 ) then
         call gwdrag(lg,tdt,he,tg,ttg,pg,u,v,         !Inputs
     &               ron,son,bvnf,dkegw)              !Outputs
      else
         call ngwdrag(lg,tdt,pg,tg,ttg,u,v,sd_gwd(1,lg),
     &                gamma_gwd(1,lg),theta_gwd(1,lg),
     &                slope_gwd(1,lg),prf,prh,ron,son,dkegw)
         call bvnfcalc(tg,ttg,prf,dprf,pdpsk,bvnf)
      end if
c     CALL TIMER('GWDRAG  ',4)

c     Save the current temperature in the latent heating array
      do k=1,nl
         do mg=1,ln2
            hlat(mg,k) = ttg(mg,k)
         end do
      end do

c**** LARGE SCALE RAIN
      if((.not.qcloud).and.(.not.kuoflag))then
c     CALL TIMER('RAINDA  ',3)
        call rainda(lg,    !Inputs
     &              ttg,qtg,precs, !In and Out
     &              qsg,rhg,gam)           !Outputs
c     CALL TIMER('RAINDA  ',4)
      endif

c**** CONVECTION : UK scheme or Kuo or Modified Arakawa
c     CALL TIMER('CONV    ',3)
      IF(ukconv)THEN
        xso4(:,:)=xtg(:,:,3)
        call convukmo(ipass,lg,tdt,pg,qsg,xso4,         !In
     &                ttg,qtg,rhg,precc,preci,fscav,    !In & Out
     &                qlg,qfg,                          !In & Out
     &                fmp,fmu,fmd,kta,kba,fluxc,rainx,qevc,CCW) !Out

      ELSEIF(kuoflag)THEN
        call hkuo(ipass,lg,tdt,pg,dqgdt, !Inputs
     &            ttg,qtg,precs,precc,   !In and Out
     &            rhg,                   !Output
     &            fmp,kta,kba,fluxc,rainx,qevc)
      ELSE
        call conv(ipass,lg,tdt,pg,qsg,gam,   !Inputs
     &                ttg,qtg,rhg,precc,qfg,   !In and Out
     &                fmp,kta,kba,fluxc,rainx,qevc)   !Outputs
      ENDIF
c     CALL TIMER('CONV    ',4)

c**** MOMENTUM MIXING BY CONVECTION
c     CALL TIMER('CVMIX   ',3)
      if(.not.kuoflag)then
        call cvmix(ipass,lg,fmp,fmu,fmd,kta,kba,u,v,pg,tdt,fscav, !Inputs
     &             uten,vten,dkevm,cten,xtg,conwd,          !In and Out
     &             xtu)                                     !Outputs
      endif
c     CALL TIMER('CVMIX   ',4)

c**** Call aerosol diagnostics
      if(coupled_aero)then
        so2wd(:)=conwd(:,2)
        so4wd(:)=conwd(:,3)
        call c_aero5(lg,tdt,xtg) !Inputs
        xtgsav(:,:,:)=xtg(:,:,:) !Save for stats
      endif

c**** New prognostic cloud scheme.
      if(qcloud)then

c**   Zero the arrays for summing cloud variables for radfs
        if(mod(mins,1440_8).eq.0_8)call rads_zer(lg)

        do mg=1,ln2
          ktop(mg)=max(kta(mg),ktop(mg))
          if(kba(mg).gt.0)kbase(mg)=min(kba(mg),kbase(mg))
        enddo

        call progcld(lg,tdt,radstep,land,fluxc,kbase,ktop,   !Inputs
     &               rainx,CCW,sg,fscav,xtu,                 !Inputs
     &               ttg,qtg,qlg,qfg,precs,cfrac,ccov,xtg,   !In and out
     &               conwd,so2wd,so4wd,                      !In and out
     &               so2oh,so2h2,so2o3,dmsoh,dmsn3,          !Outputs
     &               qlpath,preci,cldcon,qlgsav,qfgsav,cfsav,!Outputs
     &               rhg,qccon,clcon,qevap,qsubl,cdso4,qfpath)

c**   Accumulate quantities for passing into radiation
        call rads_acc(lg,cfsav,qccon,qlgsav,qfgsav)

c**   Call aerosol diagnostics
        if(coupled_aero)call c_aero6(lg,tdt,xtg,so2wd,so4wd) !Inputs

      endif ! qcloud

      do mg=1,ln2
c      Force any snowfall to melt in lowest level over water
        if(sea(mg).or.mlo(mg))then
         ttg(mg,1)=ttg(mg,1)-
     &    hlfcp*grav*preci(mg)/(0.5*dprf(mg,1)*100.0)
         preci(mg)=0.0
        endif
      enddo

c     CONDX=CONDENSATION/TIMESTEP IN KGM/M**2(=MMS)
      do 140 mg=1,ln2
 140     condx(mg)=precs(mg)+precc(mg)

c If running qcloud scheme: condx, precc and precs still refer to total precip
c i.e. liquid + solid. Solid precip is in new variable preci (LDR 3/95).
c The tau-1 value of preci is passed into surfb in variable rpreci, for 
c consistency with the treatment of condx there.

c     Calculate the latent heating from the temperature difference.
c      (Old temperature saved in hlat)
      do k=1,nl
         do mg=1,ln2
            hlat(mg,k) = ttg(mg,k) - hlat(mg,k)
         end do
      end do

c**** SURFUPB - THE MAIN SURFACE UPDATING ROUTINE (2ND PART)
      ifroz(1)=idfrz(lg,1)
      ifroz(2)=idfrz(lg,2)
      do 142 mg=1,ln2
        plold(mg)=pl(mg)
  142   sicold(mg)=siced(mg)

c     CALL TIMER('SURFUPB ',3)
      call surfupb(lg,tintp,ttg,degdt,dfgdt,eg,fg,ci,            !Inputs
     &             gamms,sg,rg,blg,condx,rcondx,preci,rpreci,he, !  "
     &             snowd,siced,tg,cie,bs,wgp,                    !In/Out
     &             tstar,sublmi,snowmi,il,perc,ifroz,wg,wg2,     !Outputs
     &             runoff)                                       !Outputs
c     CALL TIMER('SURFUPB ',4)
!      write(38,'(a21,i3,128f9.3)') 'lg, snowd after upb: ', lg, snowd
c     EVAP=EVAPORATION/TIMESTEP IN KGM/M**2(=MMS)
      do 210 mg=1,ln2
  210    evap(mg)=dt*eg(mg)/hl

c**** SOLAR AND LONG WAVE RADIATION : FELS/SCHWARTZKOPF
c     use sst rather than those corrected for elevation
c     Use NRAD=0 for no radiation (Use for debugging only - LDR 3/92)

      do 145 mg=1,ln2
         rgold(mg)=rg(mg)
  145    sgold(mg)=sg(mg)

      if (radstep) then

c**   Average cloud quantities for passing into radiation
        if(qcloud)call rads_avg(lg,cfsav,qccon,qlgsav,qfgsav)

c     CALL TIMER('RADFS   ',3)
        call radfs(ipass,lg,pg,tg,z4,ttg,qtg,rhg,als,
     &             cfsav,bvnf,qlgsav,qfgsav,qccon,cdso4,         !Inputs
     &             rg,rt,rtclr,rgclr,sg,sint,sout,soutclr,sgclr,
     &             htk,cll,clm,clh,clcon,sga,clat,refflm,cldliq) !Outputs
c     CALL TIMER('RADFS   ',4)
!      WRITE(46,*) 'After radfs, lg = ', lg
!      WRITE(43,'(200f5.2)') als

c**   Zero the arrays for summing cloud variables for radfs
        if(qcloud)call rads_zer(lg)

c**   Copy radiation output fields into global arrays to 
c      save for next NRAD timesteps
        call radh_copy(lg,ipass,sg,sga,tg,rg,htk,
     &    cl,cll,clm,clh,rgdn,blg,sgdn,als)
!      WRITE(46,*) 'After radh_copy, lg = ', lg
!      WRITE(43,'(200f5.2)') als

      else

c**   Copy radiation heating fields from global to working arrays. 
c       call radh_old(lg,ipass,htk,tg,ifroz)
        call radh_old(lg,ipass,htk,rg,ifroz)

      end if ! radstep

c
c---------------------------------------------------------------------
c Start section for leads
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c Calculate total number of points at this latitude which have sea ice
c---------------------------------------------------------------------
      lil=0
      do 220 mg=1,ln2
 220     lil=lil+il(mg)
c---------------------------------------------------------------------
c---- section for leads : only used if leads=.true. and if sea-ice
c---- points do have a fraction of open sea.
c---------------------------------------------------------------------
c 
      if (leads.and.lil.gt.0) then

c**** Second pass through selected physics routines for leads in sea-ice
       ipass=2

c**** Only zero out p.... (leads) arrays at relevant latitudes
         do 230 mg=1,ln2
            pals(mg)=0.0
            pcl(mg)=0.0
            pcll(mg)=0.0
            pclm(mg)=0.0
            pclh(mg)=0.0
            pcondx(mg)=0.0
            peg(mg)=0.0
            pevap(mg)=0.0
            pfg(mg)=0.0
            pprecc(mg)=0.0
            ppreci(mg)=0.0
            pprecs(mg)=0.0
            prg(mg)=0.0
            prgdn(mg)=0.0
            prgold(mg)=0.0
            prt(mg)=0.0
            prtclr(mg)=0.0
            prgclr(mg)=0.0
            psg(mg)=0.0
            psgdn(mg)=0.0
            psgold(mg)=0.0
            psint(mg)=0.0
            psout(mg)=0.0
            psoutclr(mg)=0.0
            psgclr(mg)=0.0
            ptaux(mg)=0.0
            ptauy(mg)=0.0
            ptg(mg)=0.0
            ptintp(mg)=0.0
            ptscrn(mg)=0.0
            prhscrn(mg)=0.0
            pv10m(mg)=0.0
 230        pvmod(mg)=0.0
         do 240 k=1,nl
            do 240 mg=1,ln2
            pclat(mg,k)=0.0
            pdkevm(mg,k)=0.0
            phlat(mg,k)=0.0
            phtk(mg,k)=0.0
            pqevc(mg,k)=0.0
            pqsg(mg,k)=0.0
            pqtg(mg,k)=0.0
            prhg(mg,k)=0.0
            pttg(mg,k)=0.0
            pu(mg,k)=0.0
            puten(mg,k)=0.0
            pv(mg,k)=0.0
            pvten(mg,k)=0.0
 240   continue
         if(sltrace)then
          do 245 nt=1,ntrace
            do 245 k=1,nl
            do 245 mg=1,ln2
 245        pcten(mg,k,nt)=0.0
         endif
c            pxtg(:,:,:)=xtg(:,:,:)
            pxtg(:,:,:)=0.

c**** TEMP AND MOISTURE VARIABLES
      call radvars(csqr,pu,pv,pqtg,pdqgdt,ipass,pttg,exptyp)

c**** Copy cloud model variables from big arrays into working arrays
      call setqcld(lg,ipass,pqlg,pqfg,pccov,gam,pgam)
        
c**** No need to call zenith again.

c**** Set surface values over ice-leads only
c     CALL TIMER('SURFICE ',3)
      call surfice(sgold,rgold,fg,eg,tg,coszro2,z4,
     &             pwetfac,pcls,pals,ptintp,ptg)
c     CALL TIMER('SURFICE ',4)

c**** Initialize surface radiation fields psg and prg
c**    from saved data in /pmasiv3/
      call radsgrg(ipass,lg,solarfit,nrad,coszro2,taudar2,ptg,
     &  psg,prg,pblg)

c**** SURFACE FLUXES (HEAT,EVAP) FROM HSFLUX
c     CALL TIMER('HSFLUXP ',3)
      call hsflux(ipass,lg,pwetfac,pg,ptg,pttg,pqtg,pu,pv,psg,snowd,!Inputs
     &            pqlg,pqfg,pccov,pgam,
     &            pvmod,ptaux,ptauy,pcduv,pdegdt,pdfgdt,peg,pfg, !Inputs
     &            pdegdw,ptscrn,pv10m,risurf,ustar,pals,v10n,zmin) !Outputs
c     CALL TIMER('HSFLUXP ',4)
       wind10(:,lg)=v10n(:)
*      print*,'rjmv2',wind10(:,lg)

c**** SURFUPA NOT CALLED FOR LEADS

c**** Compute RH at screen level
      call radrhs(pg,peg,pfg,pttg,pqtg,ptscrn,prhscrn)

c**** VERTICAL MIXING
c     CALL TIMER('HVERTMXP',3)
      call hvertmx(ipass,lg,tdt,pcduv,pfg,peg,pg,ptg,
     &             pqlg,pqfg,pccov,pgam,ptscrn,ustar,xtem,      !Inputs
     &             pqtg,pdqgdt,pu,pv,pttg,puten,pvten,pxtg,      !In and Out
     &             pdkevm,pcten,rich,kbase,ktop)            !Outputs
c     CALL TIMER('HVERTMXP',4)

c**** GRAVITY WAVE NOT CALLED FOR LEADS

c     Save the current temperature in the latent heating array
      do k=1,nl
         do mg=1,ln2
            phlat(mg,k) = pttg(mg,k)
         end do
      end do

c**** LARGE SCALE RAIN
      if((.not.qcloud).and.(.not.kuoflag))then
c     CALL TIMER('RAINDAP ',3)
        call rainda(lg,               !Inputs
     &              pttg,pqtg,pprecs,     !In and Out
     &              pqsg,prhg,pgam)       !Outputs
c     CALL TIMER('RAINDAP ',4)
      endif

c**** CONVECTION : UK scheme or Kuo or Modified Arakawa
c     CALL TIMER('CONVP   ',3)
      IF(ukconv)THEN
        xso4(:,:)=pxtg(:,:,3)
        call convukmo(ipass,lg,tdt,pg,pqsg,xso4,
     &                pttg,pqtg,prhg,pprecc,ppreci,fscav,
     &                pqlg,pqfg,
     &                pfmp,fmu,fmd,pkta,pkba,fluxc,rainx,pqevc,CCW)
      ELSEIF(kuoflag)THEN
        call hkuo(ipass,lg,tdt,pg,pdqgdt,    !Inputs
     &            pttg,pqtg,pprecs,pprecc,   !In and Out
     &            prhg,                      !Output
     &            pfmp,pkta,pkba,fluxc,rainx,pqevc)
      ELSE
        call conv(ipass,lg,tdt,pg,pqsg,pgam,   !Inputs
     &                pttg,pqtg,prhg,pprecc,pqfg,   !In and Out
     &                pfmp,pkta,pkba,fluxc,rainx,pqevc) !Outputs
      ENDIF
c     CALL TIMER('CONVP   ',4)

c**** MOMENTUM MIXING BY CONVECTION
c     CALL TIMER('CVMIXP  ',3)
      if(.not.kuoflag)then
        call cvmix(ipass,lg,pfmp,fmu,fmd,pkta,pkba,pu,pv,pg,tdt,fscav, !Inputs
     &             puten,pvten,pdkevm,pcten,pxtg,conwd,          !In and Out
     &             xtu)                                           !Outputs
      endif
c     CALL TIMER('CVMIXP  ',4)

c**** Omit call to new cloud scheme for leads.

      do mg=1,ln2
c      Force any snowfall to melt in lowest level over water
        if(il(mg).gt.0)then
         pttg(mg,1)=pttg(mg,1)-
     &    hlfcp*grav*ppreci(mg)/(0.5*dprf(mg,1)*100.0)
         ppreci(mg)=0.0
        endif
      enddo

c     PCONDX=CONDENSATION/TIMESTEP IN KGM/M**2(=MMS)
         do 310 mg=1,ln2
 310        pcondx(mg)=pprecs(mg)+pprecc(mg)

c     Calculate the latent heating from the temperature difference.
c      (Old temperature saved in phlat)
      do k=1,nl
         do mg=1,ln2
            phlat(mg,k) = pttg(mg,k) - phlat(mg,k)
         end do
      end do

c**** Update sea-surface (leads) at sea ice points
c     CALL TIMER('SURFUPL ',3)
      call surfupl(lg,z4,psg,prg,pblg,peg,pfg,
     &                   il,snowd,siced,ptg,qvent,
     &                   tstar,tg,snowmi)
c     CALL TIMER('SURFUPL ',4)
!      write(38,'(a21,i3,128f9.3)') 'lg, snowd after upl: ', lg, snowd
c     PEVAP=EVAPORATION/TIMESTEP IN KGM/M**2(=MMS)
         do 380 mg=1,ln2
  380       pevap(mg)=dt*peg(mg)/hl

c**** SOLAR AND LONG WAVE RADIATION : FELS/SCHWARTZKOPF
c     use sst rather than those corrected for elevation
c     Use NRAD=0 for no radiation (Use for debugging only - LDR 3/92)

      do 315 mg=1,ln2
         prgold(mg)=prg(mg)
  315    psgold(mg)=psg(mg)

      if (radstep) then

c     CALL TIMER('RADFSP  ',3)
        call radfs(ipass,lg,pg,ptg,z4,pttg,pqtg,prhg,pals,
     &           cfsav,bvnf,qlgsav,qfgsav,qccon,cdso4,          !Inputs
     &   prg,prt,prtclr,prgclr,psg,psint,psout,psoutclr,psgclr,
     &   phtk,pcll,pclm,pclh,clcon,psga,pclat,prefflm,pcldliq)  !Outputs
c     CALL TIMER('RADFSP  ',4)

c**   Copy radiation output fields into global arrays to 
c      save for next NRAD timesteps
        call radh_copy(lg,ipass,psg,psga,ptg,prg,phtk,
     &    pcl,pcll,pclm,pclh,prgdn,pblg,psgdn,pals)

      else

c**   Copy radiation heating fields from global to working arrays. 
c       call radh_old(lg,ipass,phtk,ptg,ifroz)
        call radh_old(lg,ipass,phtk,prg,ifroz)

      end if ! radstep

      end if ! (leads.and.lil.gt.0)
c-----------------------------------------------------------------------
c End of leads section
c-----------------------------------------------------------------------
c***********************************************************************


c**** retain tg() in tgx() for case of non-leads
      do 390 mg=1,ln2
 390     tgx(mg)=tg(mg)

c***********************************************************************

      if (.not.leads.or.(leads.and.(lil.eq.0))) then
         do 400 mg=1,ln2
            ompl(mg)=1.0-pl(mg)
            psgold(mg)=0.0
            prgold(mg)=0.0
            pfg(mg)=0.0
 400        peg(mg)=0.0
      end if

c***********************************************************************
c Collect slab ocean model driving data (before combining
c fields for non-leads and leads)

      call radslab(lg,sgold,rgold,fg,eg,hfrm,
     &                psgold,prgold,pfg,peg,
     &                flx,surfb)

c***********************************************************************
c Accumulate atmospheric fluxes used to calculate flux correction for
c coupled model (before combining fields for non-leads and leads

      if (lcouple.or.savefcor.or.statsflag) then
        call radcoupl(lg,lil,il,sgold,rgold,fg,eg,runoff,
     &         condx,evap,flx,sicold,plold,siced,sublmi,
     &         pcondx,pevap,snowmi)
      end if

c***********************************************************************
c Gather statistics on the MLO temperature error. This will be
c used as part of the Delta-T correction for the coupled model.

      if (.not.lcouple.and.(savefcor.or.statsflag)) then
        do 420 mg=1,ln2
  420   if(mlo(mg)) hdtm(mg)=hdtm(mg)+(tg(mg)-tintp(mg))
      end if

c***********************************************************************
c.... Gather ice model stresses and ocean stresses.

      call radstres(lg,lil,il,taux,tauy,ptaux,ptauy)

c***********************************************************************
c Gather solar data used by ocean model

      if (lcouple.or.savefcor) then
        call ocsave1(lg,sgold,psgold,pl,il,lil)
      end if

c***********************************************************************
c-----------------------------------------------------------------------
c Weight stats by proportion of leads over cice points if leads are
c turned on.  pl(mg) is proportion of leads.
c-----------------------------------------------------------------------
c---- Do the weighting only if required
c-----------------------------------------------------------------------
      if (leads) then
         do 450 mg=1,ln2
 450        ompl(mg)=1.0-pl(mg)
      end if

      if (leads.and.lil.gt.0) then
         call comb_1l(ompl,pl,condx ,pcondx )
         call comb_1l(ompl,pl,preci ,ppreci )
c     eg=ENERGY REQUIRED FOR EVAPORATION(+SUBL) IN WATTS/M**2
         call comb_1l(ompl,pl,eg    ,peg    )
c     evap=EVAPORATION/TIMESTEPIN KGM/M**2(=MMS OF WATER)
         call comb_1l(ompl,pl,evap  ,pevap  )
         call comb_1l(ompl,pl,fg    ,pfg    )
         call comb_1l(ompl,pl,precc ,pprecc )
         call comb_1l(ompl,pl,precs ,pprecs )
         call comb_1l(ompl,pl,rgold ,prgold )
         call comb_1l(ompl,pl,sgold ,psgold )
         call comb_1l(ompl,pl,taux  ,ptaux  )
         call comb_1l(ompl,pl,tauy  ,ptauy  )
         call comb_1l(ompl,pl,tg    ,ptg    )
         call comb_1l(ompl,pl,tscrn ,ptscrn )
         call comb_1l(ompl,pl,vmod  ,pvmod  )
         call comb_1l(ompl,pl,als   ,pals   )
         call comb_1l(ompl,pl,rhscrn,prhscrn)
         call comb_1l(ompl,pl,v10m  ,pv10m  )

c        Only accumulate radiation statistics on a radiation step,
c        except for rg and sg which are re-calculated every step.
	 if (radstep) then 
           call comb_1l(ompl,pl,rt    ,prt    )
           call comb_1l(ompl,pl,rtclr ,prtclr )
           call comb_1l(ompl,pl,rgclr ,prgclr )
           call comb_1l(ompl,pl,sint  ,psint  )
           call comb_1l(ompl,pl,sout  ,psout  )
           call comb_1l(ompl,pl,soutclr,psoutclr)
           call comb_1l(ompl,pl,sgclr ,psgclr )
           call comb_1l(ompl,pl,rgdn  ,prgdn  )
           call comb_1l(ompl,pl,sgdn  ,psgdn  )
           call comb_1l(ompl,pl,cl    ,pcl    )
           call comb_1l(ompl,pl,cll   ,pcll   )
           call comb_1l(ompl,pl,clm   ,pclm   )
           call comb_1l(ompl,pl,clh   ,pclh   )
	 end if
         call comb_nl(ompl,pl,u     ,pu     )
         call comb_nl(ompl,pl,v     ,pv     )
         call comb_nl(ompl,pl,htk   ,phtk   )
         call comb_nl(ompl,pl,ttg   ,pttg   )
         call comb_nl(ompl,pl,qtg   ,pqtg   )
         call comb_nl(ompl,pl,uten  ,puten  )
         call comb_nl(ompl,pl,vten  ,pvten  )
         call comb_nl(ompl,pl,dkevm ,pdkevm )
         call comb_nl(ompl,pl,hlat  ,phlat  )
         call comb_nl(ompl,pl,clat  ,pclat  )
         call comb_nl(ompl,pl,qevc  ,pqevc  )
         if(qcloud)then
           call comb_nl(ompl,pl,qlg   ,pqlg   )
           call comb_nl(ompl,pl,qfg   ,pqfg   )
         endif
            
c Combine aerosols for leads and non-leads
c         if(coupled_aero)then
c            call c_aero7(lg,tdt,xtg,pxtg,ompl)
c... Could use :
c...       do nt=1,ntrac
c...         call comb_nl(ompl,pl,xtg(1,1,nt),pxtg(1,1,nt))
c...       enddo
c         endif

c Don't calculate RH here in qcloud version, since the atmosphere can be 
c supersaturated at low levels over the leads at this point, because progcld
c has not been called over the leads (to save computer time).
c The extra moisture will be converted to cloud water at the next timestep;
c This should be tidied in Mark 3 by simplifying the leads stuff (LDR 11/95).

         if(.not.qcloud) call comb_nl(ompl,pl,rhg   ,prhg   )

         if(sltrace)then
           do nt=1,ntrace
             call comb_nl(ompl,pl,cten(1,1,nt),pcten(1,1,nt))
           enddo
         endif

      end if ! (leads.and.lil.gt.0)
c-----------------------------------------------------------------------
c----  end of leads weighting section
c-----------------------------------------------------------------------


c**** CREATE FRICTION TERMS (FOR USE EXPLICITLY IN DYNM)
c**     Add all vertical mixing changes to ron and son

      do 510 k=1,nl
         do 510 mg=1,ln2
         ron(mg,k)=(ron(mg,k)+uten(mg,k))*muf(mg,k)/csqr
  510    son(mg,k)=(son(mg,k)+vten(mg,k))*muf(mg,k)/csqr

c**** Add in tracer tendencies due to convection and vertical mixing
      if (sltrace) then
        lgn=lat2p-lg
        do 525 nt=1,ntrace
          do 520 k=1,nl
          do 520 mg=1,lon
           con(mg,lgn,k,nt)=con(mg,lgn,k,nt)+tdt*cten(mg,k,nt)
           ma=mg+lon
           con(mg,lg ,k,nt)=con(mg,lg ,k,nt)+tdt*cten(ma,k,nt)
  520     continue
  525   continue
      end if


c**** Calculate the energy loss due to momentum diffusion
c**    and vertical stress terms. Add to temp tendency.

      do 677 mg=1,ln2
c     TMASS= MASS OF UNIT COLM OF ATMOS (KGM/M**2)
c        tmass(mg)=pn(mg)*100.0/grav
c**   Add in vent heat from ice to lowest level (if any)
  677 htk(mg,1)=htk(mg,1)+qvent(mg)

c**** TOTAL frictional and RADIATIVE HEATING
      c100g=cp*100.0/grav
      do 660 k=1,nl
      do 660 mg=1,ln2
        dkedif=(u(mg,k)*fun(mg,k)+v(mg,k)*fvn(mg,k))*csqr/muf(mg,k)
        frht(mg,k)=(dkedif-dkevm(mg,k)-dkegw(mg,k))/cp
c**     HTK : RAD TENDANCIES (Watts/M**2) converted to
c**      RAD TEMP TENDANCIES (DEG/SEC) PER LEVEL THICKNESS
        htk(mg,k)=htk(mg,k)/(c100g*dprf(mg,k))
 660    ttg(mg,k)=ttg(mg,k)+tdt*(frht(mg,k)+htk(mg,k))

C**** RESET TEMPS MODIFIED BY LATENT HEAT,CONV ADJUSTMENT ETC
c**   RESET MOISTURE FIELDS (DUE TO RAINFALL ETC)

c**** Transfer from ttg() into ten(), and from qtg() into rmn()
      call radvars2(qtg)

      If (semice) Then
c**** Copy sea-ice model variables back into big storage array

c**   Stats for the seaice scheme
        if (leads.and.lil.gt.0) then
          do 580 mg=1,ln2
            if (pl(mg).gt.0.0) then
              plsum(mg)=plsum(mg)+ompl(mg)
              vol(mg)=vol(mg)+(dic(mg)*ompl(mg))
            end if
 580      continue
        endif

        do 750 k=1,14
        do 750 mg=1,ln2
          statsice(mg,k,lg)=datice(mg,k)
 750    continue

      Endif

c**** Copy cloud model variables back into big arrays
      if (qcloud) call setqcld(lg,3,qlg,qfg,ccov,gam,pgam)


c---------------------------------------------------------------------
c**** Gather Statistics
c---------------------------------------------------------------------

c**********************
c qcloud stats
c**********************

      if(.not.qcloud)then
        do k=1,nl
          do ns=1,2
          ma=(ns-1)*lon
          zcfrac=0.0
          do mg=1+ma,lon+ma
            zcfrac=zcfrac+clat(mg,nl+1-k)
          enddo
c Factor nrad accounts for the fact that clat=0 except on rad steps.
          cfracz(lg,ns,k)=cfracz(lg,ns,k)+nrad*zcfrac/lon ! cloud fraction
         enddo
         enddo
      endif

      if(qcloud)then
        do ns=1,2
        ma=(ns-1)*lon
        zqlpath=0.
        do mg=1+ma,lon+ma
          if (.not. land(mg)) then
            zqlpath=zqlpath+qlpath(mg)
          endif
        enddo
        qlpz(lg,ns)=qlpz(lg,ns)+zqlpath !Cloud liq water path
        enddo

c Calculate monthly mean precipitable water diagnostic 
c and evaporation and sublimation of precip.
      
        do k=1,nl
          do mg=1,ln2
            hpwc(mg)=hpwc(mg)+100.*qtg(mg,k)*dprf(mg,k)/grav
            hrevap(mg)=hrevap(mg) !Includes convective rain
     &              +50.*(qevap(mg,k)+qevc(mg,k))*dprf(mg,k)/grav
            hssubl(mg)=hssubl(mg)+50.*qsubl(mg,k)*dprf(mg,k)/grav
c Factor 50 instead of 100 accounts for leapfrog timestep.
          enddo
        enddo
        do mg=1,ln2
          hqlp(mg)=hqlp(mg)+qlpath(mg)
          hqfp(mg)=hqfp(mg)+qfpath(mg)
          hpreci(mg)=hpreci(mg)+preci(mg)
        enddo
      endif

c Gather aerosol statisticss
        if(coupled_aero)then
          xtgsav(:,:,:)=0.5*(xtgsav(:,:,:)+xtg(:,:,:)) !Avge before/after chemistry
          call c_aero8(lg,tdt,xtgsav,so2wd,so4wd,so2oh,so2h2,so2o3
     &   ,so2dd,so4dd,sem,dmsoh,dmsn3,conwd)

          do nt=1,3
            do k=1,nl
              do mg=1,ln2
                hxtg(mg,lg,k,nt)=hxtg(mg,lg,k,nt)+xtg(mg,k,nt)
     &                           *1.e12*28.96/32. ! in pptv
              enddo
            enddo
          enddo
        endif

c**********************
c Gather general stats
c**********************

      do 490 mg=1,ln2
         hvmod(mg)=hvmod(mg)+vmod(mg)
         htaux(mg)=htaux(mg)+taux(mg)
         htauy(mg)=htauy(mg)+tauy(mg)
         hrg(mg)=hrg(mg)+rgold(mg)
         hsg(mg)=hsg(mg)+sgold(mg)
         hrnc(mg)=hrnc(mg)+precc(mg)
         hrn(mg)=hrn(mg)+condx(mg)
c Convective cloud stats are best done every timestep...
         hclc(mg)=hclc(mg)+cldcon(mg)
c COLLECT SURFACE STATS
         htscrn(mg)=htscrn(mg)+tscrn(mg)
         htst(mg)=htst(mg)+tgx(mg)
         hwfg(mg)=hwfg(mg)+wg(mg)
         hwfb(mg)=hwfb(mg)+wg2(mg)
         htb2(mg)=htb2(mg)+tb2(mg)
         htb3(mg)=htb3(mg)+tb3(mg)
         hsnd(mg)=hsnd(mg)+snowd(mg)
         hsid(mg)=hsid(mg)+siced(mg)
         hrnf(mg)=hrnf(mg)+runoff(mg)   !runoff
         hals(mg)=hals(mg)+als(mg)       !surface albedo
c use this slot for looking at u*3 stats
c        u3=(taux(mg)**2+tauy(mg)**2)**0.75
c        fflux(mg)=fflux(mg)+u3                     
         hflux(mg)=hflux(mg)+fg(mg)
c COLLECT EVAP/SUBLM STATS (EVAP IF NEGATIVE=RAINFALL)
         hevap(mg)=hevap(mg)+0.5*(evap(mg)+abs(evap(mg)))
         hrn(mg)=hrn(mg)-0.5*(evap(mg)-abs(evap(mg)))
c Set potential evap and scaling evap to equal evap over non-land points
         if(.not.land(mg))then
            totpev(mg)=evap(mg)
            scalev(mg)=evap(mg)
         endif
         opreci(mg,lg)=preci(mg) ! tau-1 precip for nsib scheme
         ocondx(mg,lg)=condx(mg) ! tau-1 precip for nsib scheme
         drain(mg,lg)=drain(mg,lg)+condx(mg)
         brain(mg,lg)=brain(mg,lg)+condx(mg)
         devap(mg,lg)=devap(mg,lg)+evap(mg)
         slhf(mg,lg)=slhf(mg,lg)+eg(mg)
         sshf(mg,lg)=sshf(mg,lg)+fg(mg)
         sswr(mg,lg)=sswr(mg,lg)+sgold(mg)
         slwr(mg,lg)=slwr(mg,lg)+rgold(mg)
c Update daily extreme temperatures and global extreme surface temperatures.
         tsmax(mg,lg)=max(tscrn(mg),tsmax(mg,lg))
         tsmin(mg,lg)=min(tscrn(mg),tsmin(mg,lg))
         tgmin(mg,lg)=min(tgmin(mg,lg),tstar(mg))
         tgmax(mg,lg)=max(tgmax(mg,lg),tstar(mg))
 490  continue


c These are from the land surface scheme (NSiB/CABLE)
      do 571 mg=1,ln2
         hintrc(mg)=hintrc(mg)+tddd(mg)  !canopy interception
         htgg(mg)=htgg(mg)+tggsl(mg,1,lg) !bare ground temp
         htgf(mg)=htgf(mg)+tgf(mg,lg) !vegetated ground temp
         hpev(mg)=hpev(mg)+totpev(mg)     !potential evap
         hsev(mg)=hsev(mg)+scalev(mg)     !scaling evap
         hperc(mg)=hperc(mg)+perc(mg)     !moisture percolation
c Update daily extreme temperatures and global extreme surface temperatures.
         tggmin(mg,lg)=min(tggmin(mg,lg),tggsl(mg,1,lg))
         tggmax(mg,lg)=max(tggmax(mg,lg),tggsl(mg,1,lg))
         tgfmin(mg,lg)=min(tgfmin(mg,lg),tgf(mg,lg))
         tgfmax(mg,lg)=max(tgfmax(mg,lg),tgf(mg,lg))
         hpmc(mg,lg)=hpmc(mg,lg)+pmc(mg,lg)
  571 continue
!      WRITE(46,*) 'lg,hperc2 = ', lg,hperc

      if (nestflag) then
         do 500 mg=1,ln2
            nest_alb(mg,lg)=als(mg)
            nest_rain(mg,lg)=nest_rain(mg,lg)+condx(mg)
 500     continue
      end if

      if (radstep) then 

c**************************************************************
c     Accumulate radiation statistics on a radiation step,
c     (except for rg and sg which are re-calculated every step).
c**************************************************************

	do mg=1,ln2
          hrt(mg)=hrt(mg)+rt(mg)
          hrtclr(mg)=hrtclr(mg)+rtclr(mg)
          hrgclr(mg)=hrgclr(mg)+rgclr(mg)
          hsoutclr(mg)=hsoutclr(mg)+soutclr(mg)
          hsgclr(mg)=hsgclr(mg)+sgclr(mg)
          hsout(mg)=hsout(mg)+sout(mg)
          hsgdn(mg)=hsgdn(mg)+sgdn(mg)
          hrgdn(mg)=hrgdn(mg)+rgdn(mg)
          hcld(mg)=hcld(mg)+cl(mg)
          hcll(mg)=hcll(mg)+cll(mg)
          hclm(mg)=hclm(mg)+clm(mg)
          hclh(mg)=hclh(mg)+clh(mg)
	end do
        if(qcloud)then !Reff for liq. clouds, weighted by liq. cloud fraction
          do mg=1,ln2
            hreffl(mg)=hreffl(mg)+1.e6*refflm(mg)*cldliq(mg)
            hcliq(mg)=hcliq(mg)+cldliq(mg)
          enddo
        endif

c**********************
c Zonal means at radstep
c**********************

        do ns=1,2
        ma=(ns-1)*lon

        zclm=0.0
        zclh=0.0
        zcll=0.0
        zcl=0.0
        zrt=0.0
        zrtclr=0.0
        zsin=0.0
        zsou=0.0
        zsouclr=0.0
        zalb=0.0
        do mg=1+ma,lon+ma
          zclm=zclm+clm(mg)
          zclh=zclh+clh(mg)
          zcll=zcll+cll(mg)
          zcl=zcl+cl(mg)
          zrt=zrt+rt(mg)
          zrtclr=zrtclr+rtclr(mg)
          zsin=zsin+sint(mg)
          zsou=zsou+sout(mg)
          zsouclr=zsouclr+soutclr(mg)
          zalb=zalb+sout(mg)/(sint(mg)+0.1e-10)
	end do
        clmz(lg,ns)=clmz(lg,ns)+zclm/lon
        zclmk(lg,ns)=zclm
        clhz(lg,ns)=clhz(lg,ns)+zclh/lon
        zclhk(lg,ns)=zclh
        cllz(lg,ns)=cllz(lg,ns)+zcll/lon
        zcllk(lg,ns)=zcll
        clz(lg,ns)=clz(lg,ns)+zcl/lon
        zclk(lg,ns)=zcl
        rtz(lg,ns)=rtz(lg,ns)+zrt/lon
        rtclrz(lg,ns)=rtclrz(lg,ns)+zrtclr/lon
        zrtk(lg,ns)=zrt
        sinz(lg,ns)=sinz(lg,ns)+zsin/lon
        zsink(lg,ns)=zsin
        souz(lg,ns)=souz(lg,ns)+zsou/lon
        souclrz(lg,ns)=souclrz(lg,ns)+zsouclr/lon
        zsouk(lg,ns)=zsou
        ineg=0
        do 730 mg=1+ma,lon+ma
 730      if (sint(mg).gt.0.0) ineg=ineg+1
        if (ineg.gt.0) zalb=zalb/ineg
        albz(lg,ns)=albz(lg,ns)+zalb
        zalbk(lg,ns)=zalb

        zrad(lg,ns)=zsin-zsou-zrt
        if (clforflag) then
           zclfor=(zrtclr+zsouclr)-(zrt+zsou)
        else
           zclfor=0.
        end if
        zclfork(lg,ns)=zclfor

        enddo ! ns=1,2

      end if

c**************************************************************
c General zonal means
c**************************************************************

      do ns=1,2
      ma=(ns-1)*lon

      zsbl=0.0
      zsbl_l=0.0
      zsbl_i=0.0
      zsbl_s=0.0
      ztss=0.0
c     zfrac=0.0
      zsehf=0.0
      zsnd=0.0
      zsid=0.0
      do 530 mg=1+ma,lon+ma
         if(land(mg))zsbl_l=zsbl_l+surfb(mg)
         if(cice(mg))zsbl_i=zsbl_i+surfb(mg)
         if(sea(mg).or.mlo(mg))zsbl_s=zsbl_s+surfb(mg)
         if(.not.land(mg))then
            zsbl=zsbl+athf(mg,lg)
c**   non land temperature average
            ztss=ztss+tgx(mg)
         endif
         if(sea(mg).or.mlo(mg))then
            zsehf=zsehf+athf(mg,lg)
c           zfrac=zfrac+1.0
         endif
         plx=1.0
         if(leads.and.cice(mg))plx=1.-pl(mg)
         zsnd=zsnd+snowd(mg)*plx
         zsid=zsid+siced(mg)*plx
  530 continue
      sblz(lg,ns)=sblz(lg,ns)+zsbl
      sblz_l(lg,ns)=sblz_l(lg,ns)+zsbl_l/lon
      sblz_i(lg,ns)=sblz_i(lg,ns)+zsbl_i/lon
      sblz_s(lg,ns)=sblz_s(lg,ns)+zsbl_s/lon
      tssz(lg,ns)=tssz(lg,ns)+ztss
c     zfrack(lg,ns)=zfrac
      zsndk(lg,ns)=zsnd
      zsidk(lg,ns)=zsid


      zsbal=0.0
      zrn=0.0
      zrnc=0.0
      zrns=0.0
      zrg=0.0
      zsg=0.0
      zclc=0.0
      ztscrn=0.0
      zts=0.0
      z_tmin=999.
      z_tmax=-999.
      zhf=0.0
      zev=0.0
c     zeg=0.0
      do 540 mg=1+ma,lon+ma
         zsbal=zsbal+surfb(mg)
         zrn=zrn+condx(mg)
         zrns=zrns+precs(mg)
         zrnc=zrnc+precc(mg)
         zrg=zrg+rgold(mg)
         zsg=zsg+sgold(mg)
c Convective cloud stats are best done every timestep...
         zclc=zclc+cldcon(mg)
         ztscrn=ztscrn+tscrn(mg)
         zts=zts+tgx(mg)
         z_tmin=min(z_tmin,tgx(mg))
         z_tmax=max(z_tmax,tgx(mg))
         zhf=zhf+fg(mg)
         zev=zev+evap(mg)
c        zeg=zeg+eg(mg)
 540  continue
      zsbalk(lg,ns)=zsbal
      rnz(lg,ns)=rnz(lg,ns)+zrn/lon
      rnsz(lg,ns)=rnsz(lg,ns)+zrns/lon
      rncz(lg,ns)=rncz(lg,ns)+zrnc/lon
      zrnk(lg,ns)=zrn
      rgz(lg,ns)=rgz(lg,ns)+zrg/lon
      zrgk(lg,ns)=zrg
      sgz(lg,ns)=sgz(lg,ns)+zsg/lon
      zsgk(lg,ns)=zsg
      clcz(lg,ns)=clcz(lg,ns)+zclc/lon
      zclck(lg,ns)=zclc
      taz(lg,ns)=taz(lg,ns)+ztscrn/lon
      tsz(lg,ns)=tsz(lg,ns)+zts/lon
      ztsk(lg,ns)=zts
      z_tmink(lg,ns)=z_tmin
      z_tmaxk(lg,ns)=z_tmax
      hfz(lg,ns)=hfz(lg,ns)+zhf/lon
      zhfk(lg,ns)=zhf
      evz(lg,ns)=evz(lg,ns)+zev/lon
      zevk(lg,ns)=zev
      zeg=zev*hl/dt

c**   TAKE CARE OF NUMBER OF LAND OR NON-LAND POINTS IN SUB PRZAV
c**    FOR GWZ,TSLZ(LAND ONLY) AND SBLZ,TSSZ(NON-LAND ONLY).
      gwz(lg,ns)=gwz(lg,ns)+zgw(ns)
      zgwk(lg,ns)=zgw(ns)
      tslz(lg,ns)=tslz(lg,ns)+ztsl(ns)
      ztslk(lg,ns)=ztsl(ns)
      zrunofk(lg,ns)=zrunof(ns)
      zfrmk(lg,ns)=zfrm(ns)
      zflxik(lg,ns)=zflxi(ns)

c**********************
c Energy balance terms
c**********************
      asbalz(lg,ns)=asbalz(lg,ns)+zrad(lg,ns)/lon
      atbalz(lg,ns)=atbalz(lg,ns)+(zrad(lg,ns)-zsg+zrg+zhf+
     & (zrn*hl/dt))/lon
      sfbalz(lg,ns)=sfbalz(lg,ns)+
     & (zsg-zrg-zhf-zeg)/lon
c    & (zsg-zrg-zhf-zeg+zfrm(ns)+zflxi(ns))/lon
      flxicz(lg,ns)=flxicz(lg,ns)+zflxi(ns)/lon ! Sub-ice flux
      flxmlz(lg,ns)=flxmlz(lg,ns)+zfrm(ns)/lon  ! MLO term
      flxsez(lg,ns)=flxsez(lg,ns)+zsehf/lon ! Sea+mlo heat flux
      rebalz(lg,ns)=rebalz(lg,ns)+((zrn*hl/dt)-zeg)/lon

      enddo ! ns=1,2

c**********************
c Digital cloud maps
c**********************
      if (radstep.and.cldm) then
          if(mod(mins+int(nrad*mstep, 8),1440_8).eq.0_8) call cloudm(lg)
      endif

c**********************
c zonal mean moisture values
c**********************

      do 640 k=1,nl
        do ns=1,2
         ma=(ns-1)*lon
         zcolwv=0.0
         zrh=0.0
         qmk=0.0
         do 630 mg=1+ma,lon+ma
            zcolwv=zcolwv+
     &       (qtg(mg,k)+qlg(mg,k)+qfg(mg,k))*dprf(mg,k)
            zrh=zrh+rhg(mg,k)
 630        qmk=qmk+qtg(mg,k)*pn(mg)
         colwvz(lg,ns)=colwvz(lg,ns)+zcolwv*100.0/(grav*lon)
         rhz(lg,ns,k)=rhz(lg,ns,k)+zrh*100.0/lon
         qmean(lg,ns,k)=qmk/lon
        enddo
  640 continue

c**********************
c level zonal means
c**********************

      do 702 k=1,nl
        do ns=1,2
         ma=(ns-1)*lon
         tzon=0.0  ! zonal mean of temperature
         ezon=0.0  ! zonal mean of KE
         uzox=0.0  ! zonal mean of wind (uzon)
         zhtt=0.0  ! zonal mean of heating
         edifz=0.0 ! zonal mean of frictional heating
         do 700 mg=1+ma,lon+ma
            zhtt=zhtt+htk(mg,k)
            ezon=ezon+(u(mg,k)**2+v(mg,k)**2)*0.5
            uzox=uzox+u(mg,k)/csqr
            tzon=tzon+ttg(mg,k)
  700       edifz=edifz+frht(mg,k)
         htz(lg,ns,k)=htz(lg,ns,k)+zhtt/lon
         edifzk(lg,ns,k)=edifz*86400.0
         tzonk(lg,ns,k)=tzon
         uzon(k,ns)=uzox
         ezonk(lg,ns,k)=ezon
        enddo
  702 continue

c**** Output the detailed surface diagnostics for selected points
      if (sdiagflag) then
         do 740 nt=1,nsdiag
            if (lg.eq.lgdiag(nt)) then
               ns=insdiag(nt)
               mg=mgdiag(nt)
               ma=mg+(ns-1)*lon
           call surfdiag(nt,mins,ns,lg,mg,ma,tg,tb2,tb3,snowd,siced,
     &           precs,precc,eg,fg,sgold,sint,sout,soutclr,rgold,
     &           rt,rtclr,als,cl,cll,clm,clh,pg,ttg,qtg,u,v,cduv,tscrn,
     &           runoff,tddd,totpev,scalev)
            end if
 740     continue
      end if
c BP output 4 points for checking
c      IF (lg .eq. 9) THEN
c        k=11
c        WRITE(42,'(43f9.2)') ttg(k,1), tg(k), tgf(k,lg), tgg(k,lg),
c     &      sg(k), rg(k), mc(k,lg), osnowd(k,lg), snage(k,lg),
c     &      ssdnn(k,lg), ssdn3(k,1,lg), ssdn3(k,2,lg), ssdn3(k,3,lg),
c     &      smass(k,1,lg), smass(k,2,lg), smass(k,3,lg),
c     &      als(k), eg(k), fg(k), gflux(k,lg), sgflux(k,lg), snowd(k),
c     &      tggsl(k,1,lg), tggsl(k,2,lg), tggsl(k,3,lg), tggsl(k,4,lg),
c     &      tggsl(k,5,lg), tggsl(k,6,lg), tggsn(k,1,lg), tggsn(k,2,lg),
c     &      tggsn(k,3,lg), wb(k,1,lg), wb(k,2,lg), wb(k,3,lg), 
c     &      wb(k,4,lg), wb(k,5,lg), wb(k,6,lg), wbice(k,1,lg),
c     &      wbice(k,2,lg), wbice(k,3,lg), wbice(k,4,lg), wbice(k,5,lg),
c     &      wbice(k,6,lg)
c        k=25
c        WRITE(43,'(43f9.2)') ttg(k,1), tg(k), tgf(k,lg), tgg(k,lg),
c     &      sg(k), rg(k), mc(k,lg), osnowd(k,lg), snage(k,lg),
c     &      ssdnn(k,lg), ssdn3(k,1,lg), ssdn3(k,2,lg), ssdn3(k,3,lg), 
c     &      smass(k,1,lg), smass(k,2,lg), smass(k,3,lg),
c     &      als(k), eg(k), fg(k), gflux(k,lg), sgflux(k,lg), snowd(k),
c     &      tggsl(k,1,lg), tggsl(k,2,lg), tggsl(k,3,lg), tggsl(k,4,lg),
c     &      tggsl(k,5,lg), tggsl(k,6,lg), tggsn(k,1,lg), tggsn(k,2,lg),
c     &      tggsn(k,3,lg), wb(k,1,lg), wb(k,2,lg), wb(k,3,lg), 
c     &      wb(k,4,lg), wb(k,5,lg), wb(k,6,lg), wbice(k,1,lg), 
c     &      wbice(k,2,lg), wbice(k,3,lg), wbice(k,4,lg), wbice(k,5,lg),
c     &      wbice(k,6,lg)
c      ENDIF
c      IF (lg .eq. 21) THEN
c        k=87
cc     'tk, trad, tv, tgg1, SWdn, netLW, cansto, osnowd, snage, ssdnn, ssdn(1-3), smass(1-3), albedo, fe, fh, gflux, sgflux, snowd, tggsl(1-6), tggsn(1-3), wb(1-6), wbice(1-6)'
c        WRITE(44,'(43f9.2)') ttg(k,1), tg(k), tgf(k,lg), tgg(k,lg),
c     &      sg(k), rg(k), mc(k,lg), osnowd(k,lg), snage(k,lg),
c     &      ssdnn(k,lg), ssdn3(k,1,lg), ssdn3(k,2,lg), ssdn3(k,3,lg), 
c     &      smass(k,1,lg), smass(k,2,lg), smass(k,3,lg),
c     &      als(k), eg(k), fg(k), gflux(k,lg), sgflux(k,lg), snowd(k), 
c     &      tggsl(k,1,lg), tggsl(k,2,lg), tggsl(k,3,lg), tggsl(k,4,lg),
c     &      tggsl(k,5,lg), tggsl(k,6,lg), tggsn(k,1,lg), tggsn(k,2,lg),
c     &      tggsn(k,3,lg), wb(k,1,lg), wb(k,2,lg), wb(k,3,lg),
c     &      wb(k,4,lg), wb(k,5,lg), wb(k,6,lg), wbice(k,1,lg),
c     &      wbice(k,2,lg), wbice(k,3,lg), wbice(k,4,lg), wbice(k,5,lg),
c     &      wbice(k,6,lg)
c      ENDIF
c      IF (lg .eq. 23) THEN
c        k=4
c        WRITE(45,'(43f9.2)') ttg(k,1), tg(k), tgf(k,lg), tgg(k,lg),
c     &      sg(k), rg(k), mc(k,lg), osnowd(k,lg), snage(k,lg),
c     &      ssdnn(k,lg), ssdn3(k,1,lg), ssdn3(k,2,lg), ssdn3(k,3,lg), 
c     &      smass(k,1,lg), smass(k,2,lg), smass(k,3,lg),
c     &      als(k), eg(k), fg(k), gflux(k,lg), sgflux(k,lg), snowd(k), 
c     &      tggsl(k,1,lg), tggsl(k,2,lg), tggsl(k,3,lg), tggsl(k,4,lg),
c     &      tggsl(k,5,lg), tggsl(k,6,lg), tggsn(k,1,lg), tggsn(k,2,lg),
c     &      tggsn(k,3,lg), wb(k,1,lg), wb(k,2,lg), wb(k,3,lg),
c     &      wb(k,4,lg), wb(k,5,lg), wb(k,6,lg), wbice(k,1,lg),
c     &      wbice(k,2,lg), wbice(k,3,lg), wbice(k,4,lg), wbice(k,5,lg),
c     &      wbice(k,6,lg)
c      ENDIF

c**** Accumulated variables can only be saved in first history file.
      if (savehist(1)) then  
         call hist_acc(mins,lg,radstep,
     &       condx,precc,evap,fg,sgold,sgclr,sgdn,sint,sout,
     &       soutclr,rgold,rgclr,rgdn,rt,rtclr,tscrn,tg,
     &       taux,tauy,runoff,tddd,totpev,scalev,cl,rhscrn,v10m)
      end if

c**** Accumulate latent heating statistics
      do k=1,nl
         do mg=1,ln2
            hlat_a(mg,lg,k) = hlat_a(mg,lg,k) + hlat(mg,k)
         end do
      end do
c**** Accumulate cloud statistics
      do k=1,nl
         do mg=1,ln2
            clat_a(mg,lg,k) = clat_a(mg,lg,k) + clat(mg,nl+1-k)
         end do
      end do
c - put total cloud in at level nl (no cloud at top level)
C***      do mg=1,ln2
C***        clat_a(mg,lg,nl) = clat_a(mg,lg,nl) + cl(mg)
C***      end do

c**** Replace aerosol bigxtg using xtg
      if(coupled_aero)call c_aero9(lg,tdt,xtg)

      return
      end
