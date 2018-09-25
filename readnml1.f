c Integrating Jess Trevena's freshwater hosing source code into the model.
c SJP 2009/08/06
c
c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Enhancing user control over the coupling between the atmosphere and ocean.
c SJP 2009/04/21
c
c nsib variable is taken  out and replaced with a character(5) variable
c lsm_type = "nsib "  
c AJA 2009/01/22
c
c Rewrite the section of code that checks that the values of MSTEP and NRAD are
c valid, thus fixing a bug whereby MSTEP had to divide into 120 minutes, even
c when NRAD was negative. When NRAD is positive and MSTEP*NRAD does not equal
c 120, the model now resets the value of NRAD, rather than simply aborting.
c SJP 2008/02/03
c
c Re-inserted the line "include 'PARAMS.f'" into subroutines READNML1 and
c READNML2, as this line is no longer included via the header file OPARAMS.f.
c SJP 2007/05/31
c
c Transferred COMMON blocks to separate header files, as follows:
c /MDAY/  ->  MDAY.f
c SJP 2007/05/29
c
c Modified to read CRELAX_FLAG and CRELAX_TAU from namelist input. The default
c values are set to .FALSE. and 3650.0 respectively.
c SJP 2006/01/05
c
c Modified for five-digit year numbers.
c SJP 2004/09/29
c
c Modified to read RCRITL, RCRITS, REFAC1 and REFAC2 from namelist input. The
c default values are set to 0.75, 0.85, 0.85 and 0.95 respectively - see
c newcloud.f and cloud2.f
c SJP 2004/09/24
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Added the flag SUBICE, which controls whether or not the sub-ice heat input
c is used.
c SJP 2004/01/05
c
c Modified to add BPYEAR and CSOLAR to the namelist input. The default
c values (0 and 1366.7733 W/m^2 respectively) correspond to the values
c previously hard-coded into the model (0 and 1.96 ly/min respectively).
c SJP 2001/12/23
c
c $Log: readnml1.f,v $
c Revision 1.74  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.73  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.72.1.1  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.72  2001/02/22 05:34:45  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.71  2001/02/12 05:39:56  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.70  2000/11/16 01:42:34  rot032
c Warning message for LBRMAP (obsolete).
c
c Revision 1.69  2000/11/14 06:55:45  rot032
c New PI_emissions flag.
c
c Revision 1.68  2000/11/14 03:11:39  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.67  2000/08/16 02:59:26  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.66  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.65  2000/06/20 02:08:35  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.64.1.1  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.64  1999/06/30 05:29:37  rot032
c Mods for SCM.
c
c Revision 1.63  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.62  1999/06/17 05:09:02  rot032
c More appropriate defaults for control variables.
c
c Revision 1.61  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.60  1998/12/10  00:56:01  ldr
c HBG changes to V5-1-21
c
c Revision 1.59  1998/05/26  01:20:31  ldr
c Change some of the flag defaults.
c
c Revision 1.58  1997/12/23  00:55:47  ldr
c Merge of LDR changes with HBG changes.
c
c Revision 1.57  1997/12/19  02:03:19  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.56.1.1  1997/12/19  06:06:20  ldr
c Flexible treatment of direct and indirect aerosol effects.
c
c Revision 1.56  1997/10/07  03:18:16  ldr
c  Final corrections for V5-1.
c
c Revision 1.55  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.54  1997/06/11  02:21:32  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.53  1996/12/23  03:58:02  mrd
c Add new gravity wave drag scheme as an option.
c
c Revision 1.52  1996/06/13  02:08:07  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.51  1996/04/17  01:32:52  ldr
c Changes from TIE to allow initialization of monthly mean netcdf files
c from within model.
c
c Revision 1.50  1996/03/21  03:19:02  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.49  1996/02/19  04:10:02  ldr
c Generalize for 24 levels.
c
c Revision 1.48  1996/02/08  00:21:01  ldr
c Tidy up reading of qcloud stuff in filerd and remove obsolete slowice flag.
c
c Revision 1.47  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.46  1994/12/22  23:51:06  ldr
c Mods to cloud scheme of V4-6-13l to produce faster version for HBG and
c to allow vertically subgrid cloud below inversions in BL. Set slowice=T
c for accurate (but slow) ice fall calculation. Also, Ri now includes hlfusionn
c
c Revision 1.45  94/09/13  09:51:37  mrd
c Add daily mean 10m windspeed and 2m relative humidity as history variables.
c 
c Revision 1.44  94/09/09  14:53:58  mrd
c Added daily average cloud, surface temperature and screen temperature
c as history variables.
c 
c Revision 1.43  94/09/09  14:14:53  mrd
c Added possbility of saving two history files at different frequencies.
c 
c Revision 1.42  94/08/09  12:33:53  ldr
c Rationalize zonal mean diagnostics: Use conzp flag to control HBG's
c cloud map and introduce savezonu, savezonv for dump of zmean u,v fields.
c Also remove savezcls which is fairly irrelevant.
c 
c Revision 1.41  94/08/08  17:22:29  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.40  94/08/08  14:41:43  ldr
c Add flag c_sflg for cloud stats and make ksc flag obsolete.
c 
c Revision 1.39  94/06/30  10:02:37  mrd
c Add sea-level pressure as a history variable.
c 
c Revision 1.38  94/06/30  09:53:54  ldr
c Add more zonal mean diagnostic plot options.
c 
c Revision 1.37  94/03/30  10:20:07  ldr
c Add new debug diagnostic.
c 
c Revision 1.36  93/11/29  14:48:19  ldr
c Use icu_sflg and icv_sflg for treatment of ice U,V diagnostics, so they
c are consistent with other "stats".
c 
c Revision 1.35  93/11/03  13:07:43  ldr
c Make locean a namelist flag and Replace silly flag zfl_sflg with savezflux.
c 
c Revision 1.34  93/10/11  15:43:37  ldr
c Create new namelist groups histvars and statvars.
c 
c Revision 1.33  93/10/08  11:19:46  ldr
c Corrected hflg data statement and namelist statement.
c 
      subroutine readnml1

c SET DEFAULT VALUES OF NAMELIST PARAMETERS
c AND OTHER PRESET VALUES.
C THEN READ NAMELISTS

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'OPARAMS.f'

C Global data blocks
      include 'AOGCM.f'
      include 'CLOUDPAR.f'
      include 'CRELAX.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'HIST.f'
      include 'HOSING.f'
      include 'NEST.f'
      include 'ORBRAD.f'
      include 'PRINTT.f'
      include 'STFLAGS.f'
      include 'TIMEX.f'
      include 'MDAY.f'

      integer kuocb,kbconv
      real rhcut,rhkuo,rhsat
      common/kuocom/kuocb,kbconv,rhcut,rhkuo,rhsat

      integer nst63,lgt63,mgt63,monthset
      common /scm_spacetime/ nst63,lgt63,mgt63,monthset

      integer insdiag,mgdiag,lgdiag,nsdiag
      common /sdiag/ insdiag(15),mgdiag(15),lgdiag(15),nsdiag

C Local work arrays and variables
      integer kk
      integer lg
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c************************************************************************
c****  Set default values for various parameters in common blocks: ******
c************************************************************************

c..   'AOGCM.f'
c..   integer isync
c..   real ovolume
c..   common /aogcm/ isync, ovolume
      isync = 1
      ovolume = 1.27e18

c..   'CLOUDPAR.f'
c..   real rcritl
c..   real rcrits
c..   real refac1
c..   real refac2
c..   common /cloudpar/ rcritl, rcrits, refac1, refac2

      rcritl = 0.75
      rcrits = 0.85
      refac1 = 0.85
      refac2 = 0.95

c..   'CRELAX.f'
c..   logical crelax_flag
c..   real crelax_tau, crelax_athf, crelax_atsf
c..   common /crelax/ crelax_flag, crelax_tau, crelax_athf(imt, jmt),
c..  &                crelax_atsf(imt, jmt)

      crelax_flag = .false.
      crelax_tau = 3650.0   

c..   'FEWFLAGS.f'
c..   logical ifwds,iener,ispec,cdmap,mlomap,filewrflag
c..  &    ,ltrace,statsflag,rainflag,tempflag,uvtflag,qflux
c..  &    ,lcouple,locean,sdiagflag,saveqflux,savegbrf,saveglmean
c..  &    ,semice,leads,idyn,plotheat,plotclds,plotnetr,plotevrn
c..  &    ,clforflag,sltrace,chenflag,laust,impvor,savefcor
c..  &    ,saveicuv,qcloud,autoname,savezflux,debug,hybrid
c..  &    ,ukconv,kuoflag,ncarpbl,jsfdiff,amipo3,coupled_aero,fluxadj
c..  &    ,ncepagcm,PI_emissions,newriver,subice
c..   integer incd,nsstop,months,nsemilag,lastmonth,nthreads,idayp
c..  &    ,ksc,glmean_interval,ndstop,insdebug,lgdebug,mgdebug,ngwd
c..  &    ,naerosol_d,naerosol_i(2)
c..   character*5 lsm_type
c..   common/fewflags/ifwds,iener,ispec,cdmap,mlomap,filewrflag
c..  &    ,ltrace,statsflag,rainflag,tempflag,uvtflag,qflux
c..  &    ,lcouple,locean,sdiagflag,saveqflux,savegbrf,saveglmean
c..  &    ,semice,leads,idyn,plotheat,plotclds,plotnetr,plotevrn
c..  &    ,clforflag,sltrace,chenflag,laust,impvor,savefcor
c..  &    ,saveicuv,qcloud,autoname,savezflux,debug,hybrid
c..  &    ,ukconv,kuoflag,ncarpbl,jsfdiff,amipo3,coupled_aero,fluxadj
c..  &    ,ncepagcm,PI_emissions,newriver
c..  &    ,incd,nsstop,months,nsemilag,lastmonth,nthreads,idayp
c..  &    ,ksc,glmean_interval,ndstop,insdebug,lgdebug,mgdebug,ngwd
c..  &    ,naerosol_d,naerosol_i,subice,lsm_type

      iener       =.false.
      ispec       =.false.
      cdmap       =.false.
      mlomap      =.false.
      filewrflag  =.false.
      ltrace      =.false.
      statsflag   =.false.
      rainflag    =.false.
      tempflag    =.false.
      uvtflag     =.false.
      qflux       =.false.
      lcouple     =.false.
      locean      =.false.
      sdiagflag   =.false.
      saveqflux   =.false.
      savegbrf    =.false.
      saveglmean  =.false.
      semice      =.true.
      leads       =.true.
      idyn        =.true.
      plotheat    =.false.
      plotclds    =.false.
      plotnetr    =.false.
      plotevrn    =.false.
      clforflag   =.true.
      sltrace     =.false.
      chenflag    =.true.
      laust       =.false.
      impvor      =.true.
      savefcor    =.false.
      saveicuv    =.false. ! Obsolete
      qcloud      =.true.
      autoname    =.false.
      savezflux   =.false.
      debug       =.false.
      hybrid      =.true.
      ukconv      =.true.  ! if true, kuoflag must be false
      kuoflag     =.false. ! if true, ukconv   "   "   "
c                          !  or both false 
      ncarpbl     =.false.
      jsfdiff     =.false.
      amipo3      =.true.
      coupled_aero=.false.
      fluxadj     =.true.
      ncepagcm    =.false.
      PI_emissions=.false.
      newriver    =.false.
      incd        =-1
      nsstop      =0
      months      =1
      nsemilag    =1
      lastmonth   =0
      nthreads    =0
      idayp       =0
      ksc         =-1
      glmean_interval=480
      ndstop      =0
      insdebug    =1
      lgdebug     =1
      mgdebug     =1
      ngwd        =1
      naerosol_d  =2
      naerosol_i(:)=0
      subice = .true.
      lsm_type    = "cable" ! or "nsib "

c..   'FILES.f'
c..   character*32 irfilename,isfilename,orfilename,osfilename
c..   character*32 co2_datafile,o3_datafile      
c..   character*32 so4_direct_file,so4_ind_rad_file,so4_ind_rain_file
c..   common/files32/irfilename,isfilename,orfilename,osfilename
c..  &    ,co2_datafile,o3_datafile
c..  &    ,so4_direct_file,so4_ind_rad_file,so4_ind_rain_file
c..   character*13 str
c..   common/files10/str
c..   character*5 runtype
c..   common/files03/runtype

      runtype='earth'
      irfilename=' '
      isfilename=' '
      orfilename=' '
      osfilename=' '
      co2_datafile='co2_data.hbg18'
      o3_datafile='o3_data.hbg18'
      so4_direct_file=' '
      so4_ind_rad_file=' '
      so4_ind_rain_file=' '

c..   'HIST.f'
c..   logical savehist(2)
c..   integer hist_interval(2), histid(2), histset(2)
c..   integer tid(nl,2),uid(nl,2),vid(nl,2),qid(nl,2),pid(2),plid(2)
c..   integer cldid(2),cllid(2),clmid(2),clhid(2),ichid(2),icmid(2),
c..  &        icbid(2),ictid(2),alsid(2)
c..   logical all_hflg(2),
c..  &        zht_hflg(2), t_hflg(2), u_hflg(2), v_hflg(2), q_hflg(2),
c..  &        psf_hflg(2), imsl_hflg(2),tsu_hflg(2), tb2_hflg(2),
c..  &        tb3_hflg(2), als_hflg(2), wfg_hflg(2), wfb_hflg(2),
c..  &        snd_hflg(2), sid_hflg(2), rnd_hflg(2), rnc_hflg(2),
c..  &        evp_hflg(2), hfl_hflg(2), sgn_hflg(2), sgc_hflg(2),
c..  &        sgd_hflg(2), sit_hflg(2), sot_hflg(2), soc_hflg(2),
c..  &        rgn_hflg(2), rgc_hflg(2), rgd_hflg(2), rtu_hflg(2),
c..  &        rtc_hflg(2), tax_hflg(2), tay_hflg(2), cld_hflg(2),
c..  &        cll_hflg(2), clm_hflg(2), clh_hflg(2), ich_hflg(2),
c..  &        icm_hflg(2), ict_hflg(2), icb_hflg(2), tsc_hflg(2),
c..  &        tsh_hflg(2), tsl_hflg(2), run_hflg(2), int_hflg(2),
c..  &        tgg_hflg(2), tgf_hflg(2), pev_hflg(2), sev_hflg(2),
c..  &        ico_hflg(2), tgl_hflg(2), tgh_hflg(2), psl_hflg(2),
c..  &        tsca_hflg(2), tsua_hflg(2), clda_hflg(2), rhsa_hflg(2),
c..  &        v10ma_hflg(2)
c..   common /hist_control/ savehist,hist_interval,histid,histset,
c..  &       tid,uid,vid,qid,pid,plid,cldid,cllid,clmid,clhid,
c..  &       ichid,icmid,icbid,ictid,alsid,
c..  &       all_hflg, zht_hflg, t_hflg, u_hflg, v_hflg, q_hflg,
c..  &       psf_hflg, imsl_hflg, tsu_hflg, tb2_hflg, tb3_hflg,
c..  &       als_hflg, wfg_hflg, wfb_hflg, snd_hflg, sid_hflg,
c..  &       rnd_hflg, rnc_hflg, evp_hflg, hfl_hflg, sgn_hflg,
c..  &       sgc_hflg, sgd_hflg, sit_hflg, sot_hflg, soc_hflg,
c..  &       rgn_hflg, rgc_hflg, rgd_hflg, rtu_hflg, rtc_hflg,
c..  &       tax_hflg, tay_hflg, cld_hflg, cll_hflg, clm_hflg,
c..  &       clh_hflg, ich_hflg, icm_hflg, ict_hflg, icb_hflg,
c..  &       tsc_hflg, tsh_hflg, tsl_hflg, run_hflg, int_hflg,
c..  &       tgg_hflg, tgf_hflg, pev_hflg, sev_hflg, ico_hflg,
c..  &       tgl_hflg, tgh_hflg, psl_hflg, tsca_hflg, tsua_hflg,
c..  &       clda_hflg, rhsa_hflg, v10ma_hflg
      do kk=1,2
        savehist(kk)=.false.
        hist_interval(kk)=1440

        all_hflg(kk)=.false.
        zht_hflg(kk)=.false.
        t_hflg(kk)=.false.
        u_hflg(kk)=.false.
        v_hflg(kk)=.false.
        q_hflg(kk)=.false.
        psf_hflg(kk)=.false.
        imsl_hflg(kk)=.false.
        tsu_hflg(kk)=.false.
        tb2_hflg(kk)=.false.
        tb3_hflg(kk)=.false.
        als_hflg(kk)=.false.
        wfg_hflg(kk)=.false.
        wfb_hflg(kk)=.false.
        snd_hflg(kk)=.false.
        sid_hflg(kk)=.false.
        rnd_hflg(kk)=.false.
        rnc_hflg(kk)=.false.
        evp_hflg(kk)=.false.
        hfl_hflg(kk)=.false.
        sgn_hflg(kk)=.false.
        sgc_hflg(kk)=.false.
        sgd_hflg(kk)=.false.
        sit_hflg(kk)=.false.
        sot_hflg(kk)=.false.
        soc_hflg(kk)=.false.
        rgn_hflg(kk)=.false.
        rgc_hflg(kk)=.false.
        rgd_hflg(kk)=.false.
        rtu_hflg(kk)=.false.
        rtc_hflg(kk)=.false.
        tax_hflg(kk)=.false.
        tay_hflg(kk)=.false.
        cld_hflg(kk)=.false.
        cll_hflg(kk)=.false.
        clm_hflg(kk)=.false.
        clh_hflg(kk)=.false.
        ich_hflg(kk)=.false.
        icm_hflg(kk)=.false.
        ict_hflg(kk)=.false.
        icb_hflg(kk)=.false.
        tsc_hflg(kk)=.false.
        tsh_hflg(kk)=.false.
        tsl_hflg(kk)=.false.
        run_hflg(kk)=.false.
        int_hflg(kk)=.false.
        tgg_hflg(kk)=.false.
        tgf_hflg(kk)=.false.
        pev_hflg(kk)=.false.
        sev_hflg(kk)=.false.
        ico_hflg(kk)=.false.
        tgl_hflg(kk)=.false.
        tgh_hflg(kk)=.false.
        psl_hflg(kk)=.false.
        tsca_hflg(kk)=.false.
        tsua_hflg(kk)=.false.
        clda_hflg(kk)=.false.
        rhsa_hflg(kk)=.false.
        v10ma_hflg(kk)=.false.
      enddo

c..   'HOSING.f'
c..   logical hosing_flag
c..   real hosing_rate
c..   common /hosing/ hosing_flag, hosing_rate

      hosing_flag = .false.
      hosing_rate = 1.0

c..   'NEST.f'
c..   real nest_alb(ln2,lat)
c..   real nest_rain(ln2,lat)
c..   integer nest_interval
c..   integer nest_start
c..   logical nestflag
c     common/nest/nest_alb,nest_rain,nest_interval,nest_start,nestflag
      do lg=1,lat
      do mg=1,ln2
        nest_rain(mg,lg)=0.0
      enddo
      enddo
      nest_interval=480
      nest_start=240
      nestflag=.false.

c..   ORBRAD.f
c..   common /orbrad/ bpyear, csolar, ec, peril, oblqty
      bpyear = 0
      csolar = 1366.77333333333

c..   PRINTT.f
c..   logical zavgp,dynzp,phz1p,phz2p,conzp,dynfp,rhnmm
c..  & ,gmap1,cvrnm,gwicm,cldm,gmap2,rainm,evapm,pmslm,surfm
c..  & ,savezqcl,savezcfr,savezqlp,savezont,savezonr,savezonw
c..  & ,savezonu,savezonv
c..   common/printt/zavgp,dynzp,phz1p,phz2p,conzp,dynfp,rhnmm
c..  & ,gmap1,cvrnm,gwicm,cldm,gmap2,rainm,evapm,pmslm,surfm
c..  & ,savezqcl,savezcfr,savezqlp,savezont,savezonr,savezonw
c..  & ,savezonu,savezonv
      zavgp=.true.
      dynzp=.true.
      phz1p=.true.
      phz2p=.true.
      conzp=.true.
      dynfp=.false.
      rhnmm=.false.
      gmap1=.false.
      cvrnm=.false.
      gwicm=.false.
      cldm=.false.
      gmap2=.false.
      rainm=.false.
      evapm=.false.
      pmslm=.false.
      surfm=.false.
      savezqcl=.false.
      savezcfr=.false.
      savezqlp=.false.
      savezont=.false.
      savezont=.false.
      savezonr=.false.
      savezonw=.false.
      savezonu=.false.
      savezonv=.false.

c..   STFLAGS.f
C Surface statistics flags
c..   logical 
c..  &  evp_sflg, pev_sflg, sev_sflg ,rnd_sflg, rnc_sflg, hfl_sflg
c..  & ,wfg_sflg, wfb_sflg, run_sflg ,per_sflg, int_sflg, psl_sflg
c..  & ,vmo_sflg, tax_sflg, tay_sflg ,tsu_sflg, tsc_sflg, tb2_sflg
c..  & ,tb3_sflg, tgg_sflg, tgf_sflg ,thd_sflg, tld_sflg, thg_sflg
c..  & ,tlg_sflg, thf_sflg, tlf_sflg ,thm_sflg, tlm_sflg, dtm_sflg
c..  & ,rsv_sflg
C Cloud statistics
c..  &  cld_sflg, cll_sflg, clm_sflg ,clh_sflg
C Radiation statistics
c..  &  rgn_sflg, rgd_sflg, rgc_sflg ,sgn_sflg, sgd_sflg, sgc_sflg
c..  & ,rtu_sflg, rtc_sflg, sot_sflg ,soc_sflg, als_sflg
C Snow & ice statistics
c..  &  snd_sflg, sid_sflg, ico_sflg ,itf_sflg, isf_sflg, icu_sflg
c..  & ,icv_sflg, div_sflg, gro_sflg ,ire_sflg, ich_sflg
C Fresh water flux into ocean statistics
c..  & ,fwf_sflg
C Qcloud statistics
c..  &  sno_sflg, rev_sflg, ssb_sflg ,clc_sflg, lwp_sflg, pwc_sflg
c..  & ,ref_sflg
C Model level statistics
c..  &    u_sflg,   v_sflg,   t_sflg ,  q_sflg,  rh_sflg,   g_sflg
c..  & ,  c_sflg,   l_sflg
c..
c..   common/stflags/
c..  &  evp_sflg, pev_sflg, sev_sflg ,rnd_sflg, rnc_sflg, hfl_sflg
c..  & ,wfg_sflg, wfb_sflg, run_sflg ,per_sflg, int_sflg, psl_sflg
c..  & ,vmo_sflg, tax_sflg, tay_sflg ,tsu_sflg, tsc_sflg, tb2_sflg
c..  & ,tb3_sflg, tgg_sflg, tgf_sflg ,thd_sflg, tld_sflg, thg_sflg
c..  & ,tlg_sflg, thf_sflg, tlf_sflg ,thm_sflg, tlm_sflg, dtm_sflg
c..  & ,rsv_sflg
c..  & ,cld_sflg, cll_sflg, clm_sflg ,clh_sflg
c..  & ,rgn_sflg, rgd_sflg, rgc_sflg ,sgn_sflg, sgd_sflg, sgc_sflg
c..  & ,rtu_sflg, rtc_sflg, sot_sflg ,soc_sflg, als_sflg
c..  & ,snd_sflg, sid_sflg, ico_sflg ,itf_sflg, isf_sflg, icu_sflg
c..  & ,icv_sflg, div_sflg, gro_sflg ,ire_sflg, ich_sflg
c..  & ,fwf_sflg
c..  & ,sno_sflg, rev_sflg, ssb_sflg ,clc_sflg, lwp_sflg, pwc_sflg
c..  & ,ref_sflg
c..  & ,  u_sflg,   v_sflg,   t_sflg ,  q_sflg,  rh_sflg,   g_sflg
c..  & ,  c_sflg,   l_sflg

      evp_sflg=.true.
      pev_sflg=.false.
      sev_sflg=.false.
      rnd_sflg=.true.
      rnc_sflg=.true.
      hfl_sflg=.true.
      wfg_sflg=.true.
      wfb_sflg=.true.
      run_sflg=.true.
      per_sflg=.true.
      int_sflg=.false.
      psl_sflg=.true.
      vmo_sflg=.false.
      tax_sflg=.true.
      tay_sflg=.true.
      tsu_sflg=.true.
      tsc_sflg=.true.
      tb2_sflg=.true.
      tb3_sflg=.true.
      tgg_sflg=.false.
      tgf_sflg=.false.
      thd_sflg=.false.
      tld_sflg=.false.
      thg_sflg=.false.
      tlg_sflg=.false.
      thf_sflg=.false.
      tlf_sflg=.false.
      thm_sflg=.false.
      tlm_sflg=.false.
      dtm_sflg=.false.
      rsv_sflg=.false.

      cld_sflg=.true.
      cll_sflg=.true.
      clm_sflg=.true.
      clh_sflg=.true.

      rgn_sflg=.true.
      rgd_sflg=.false.
      rgc_sflg=.false.
      sgn_sflg=.true.
      sgd_sflg=.false.
      sgc_sflg=.false.
      rtu_sflg=.true.
      rtc_sflg=.false.
      sot_sflg=.true.
      soc_sflg=.false.
      als_sflg=.true.

      snd_sflg=.true.
      sid_sflg=.true.
      ico_sflg=.true.
      itf_sflg=.false.
      isf_sflg=.false.
      icu_sflg=.true.
      icv_sflg=.true.
      div_sflg=.false.
      gro_sflg=.false.
      ire_sflg=.false.
      ich_sflg=.false.

      fwf_sflg=.false.

      sno_sflg=.true.
      rev_sflg=.true.
      ssb_sflg=.true.
      clc_sflg=.true.
      lwp_sflg=.true.
      pwc_sflg=.true.
      ref_sflg=.true.

        u_sflg=.true.
        v_sflg=.true.
        t_sflg=.true.
        q_sflg=.true.
       rh_sflg=.true.
        g_sflg=.false.
        c_sflg=.false.
        l_sflg=.false.

c..   TIMEX.f
c..   real dt,ratlm
c..   integer mstep,minw,kday,kdayp,nstep,ndays
c..  & ,kdays,month,ldays,nsteps,iyear,ncepstrt
c..  & ,nrad
c..   integer*8 mins
c..   logical start,dcflag
c..   common/timex/dt,ratlm
c..  & ,mins,mstep,minw,kday,kdayp,nstep,ndays
c..  & ,kdays,month,ldays,nsteps,iyear,ncepstrt
c..  & ,nrad
c..  & ,start,dcflag
      mstep=30
      minw=480
      nrad=4
      dcflag=.true.

c..   integer kuocb,kbconv
c..   real rhcut,rhkuo,rhsat
c..   common/kuocom/kuocb,kbconv,rhcut,rhkuo,rhsat
      kuocb=2
      kbconv=nl
      rhcut=0.5
      rhkuo=0.8
      rhsat=1.0

c..   common/mday/mdays(12)
      mdays( 1)=31
      mdays( 2)=28
      mdays( 3)=31
      mdays( 4)=30
      mdays( 5)=31
      mdays( 6)=30
      mdays( 7)=31
      mdays( 8)=31
      mdays( 9)=30
      mdays(10)=31
      mdays(11)=30
      mdays(12)=31

c..   common /scm_spacetime/ nst63,lgt63,mgt63,monthset
      nst63=1
      lgt63=29
      mgt63=141
      monthset=7

c..   common /sdiag/ insdiag(15),mgdiag(15),lgdiag(15),nsdiag
      insdiag( 1)=2
      insdiag( 2)=2
      insdiag( 3)=2
      insdiag( 4)=1
      insdiag( 5)=2
      insdiag( 6)=1
      insdiag( 7)=2
      insdiag( 8)=2
      insdiag( 9)=1
      insdiag(10)=1
      insdiag(11)=2
      insdiag(12)=0
      insdiag(13)=0
      insdiag(14)=0
      insdiag(15)=0

      lgdiag( 1)=17
      lgdiag( 2)=21
      lgdiag( 3)=4
      lgdiag( 4)=19
      lgdiag( 5)=27
      lgdiag( 6)=14
      lgdiag( 7)=9
      lgdiag( 8)=28
      lgdiag( 9)=9
      lgdiag(10)=9
      lgdiag(11)=24
      lgdiag(12)=0
      lgdiag(13)=0
      lgdiag(14)=0
      lgdiag(15)=0

      mgdiag( 1)=27
      mgdiag( 2)=25
      mgdiag( 3)=20
      mgdiag( 4)=15
      mgdiag( 5)=54
      mgdiag( 6)=48
      mgdiag( 7)=27
      mgdiag( 8)=38
      mgdiag( 9)=60
      mgdiag(10)=18
      mgdiag(11)=24
      mgdiag(12)=0
      mgdiag(13)=0
      mgdiag(14)=0
      mgdiag(15)=0

      nsdiag=11

c************************************************************************
c****  End of setting default value  ************************************
c************************************************************************

c************************************************************************
c****  Read in namelists (which may reset some of the above values) *****
c************************************************************************
      call readnml2

      return
      end
C---------------------------------------------------------------------
      subroutine readnml2

c Routine to read AGCM namelist input

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'OPARAMS.f'

C Global data blocks
      include 'AOGCM.f'
      include 'CLOUDPAR.f'
      include 'CRELAX.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'HIST.f'
      include 'HOSING.f'
      include 'NEST.f'
      include 'ORBRAD.f'
      include 'PRINTT.f'
      include 'STFLAGS.f'
      include 'TIMEX.f'

      integer kuocb,kbconv
      real rhcut,rhkuo,rhsat
      common/kuocom/kuocb,kbconv,rhcut,rhkuo,rhsat

      integer insdiag,mgdiag,lgdiag,nsdiag
      common /sdiag/ insdiag(15),mgdiag(15),lgdiag(15),nsdiag

      integer nst63,lgt63,mgt63,monthset
      common / scm_spacetime / nst63,lgt63,mgt63,monthset

C Local work arrays and variables
      integer kchk ! for checking convection switches

C Local data, functions etc
      character*4 chtst
      data chtst /'OLD'/
      logical lbrmap
      data lbrmap /.false./ 

C Start code : ----------------------------------------------------------

c NAMELIST INPUT

      namelist /control/ mstep,chtst,incd,nrad,dcflag,filewrflag,
     &     irfilename,isfilename,orfilename,osfilename,nestflag,
     &     nest_interval,nest_start,runtype,qflux,nsstop,months,
     &     semice,leads,idyn,nsemilag,lsm_type,lastmonth,co2_datafile,
     &     o3_datafile,nthreads,impvor,lcouple,locean,ndstop,qcloud,
     &     autoname,hybrid,ngwd,naerosol_d,naerosol_i,
     &     ukconv,ncarpbl,jsfdiff,amipo3,coupled_aero,ncepagcm,
     &     so4_direct_file,so4_ind_rad_file,so4_ind_rain_file,
     &     PI_emissions,newriver,bpyear,csolar,rcritl,rcrits,
     &     refac1,refac2

      namelist /diagnostics/ dynfp,minw,iener,ispec,zavgp,dynzp,phz1p,
     &     phz2p,conzp,gmap1,cvrnm,gwicm,rhnmm,cldm,gmap2,rainm,evapm,
     &     pmslm,surfm,cdmap,mlomap,idayp,glmean_interval,laust,lbrmap,
     &     ltrace,rainflag,tempflag,uvtflag,statsflag,sdiagflag,nsdiag,
     &     insdiag,mgdiag,lgdiag,saveglmean,saveqflux,savegbrf
     &     ,plotheat,plotclds,plotnetr,plotevrn,clforflag,sltrace,
     &     savefcor,saveicuv,hist_interval,savehist,savezflux,
     &     debug,insdebug,lgdebug,mgdebug,savezqcl,savezcfr,savezqlp,
     &     savezont,savezonr,savezonw,savezonu,savezonv

      namelist /histvars/
     &     all_hflg, zht_hflg, t_hflg, u_hflg, v_hflg, q_hflg,
     &     psf_hflg, imsl_hflg, tsu_hflg, tb2_hflg, tb3_hflg,
     &     als_hflg, wfg_hflg, wfb_hflg, snd_hflg, sid_hflg,
     &     rnd_hflg, rnc_hflg, evp_hflg, hfl_hflg, sgn_hflg,
     &     sgc_hflg, sgd_hflg, sit_hflg, sot_hflg, soc_hflg,
     &     rgn_hflg, rgc_hflg, rgd_hflg, rtu_hflg, rtc_hflg,
     &     tax_hflg, tay_hflg, cld_hflg, cll_hflg, clm_hflg,
     &     clh_hflg, ich_hflg, icm_hflg, ict_hflg, icb_hflg,
     &     tsc_hflg, tsh_hflg, tsl_hflg, run_hflg, int_hflg,
     &     tgg_hflg, tgf_hflg, pev_hflg, sev_hflg, ico_hflg,
     &     tgl_hflg, tgh_hflg, psl_hflg, tsca_hflg, tsua_hflg,
     &     clda_hflg, rhsa_hflg, v10ma_hflg

      namelist /statvars/
     &     evp_sflg, pev_sflg, sev_sflg ,rnd_sflg, rnc_sflg, hfl_sflg
     &    ,wfg_sflg, wfb_sflg, run_sflg ,per_sflg, int_sflg, psl_sflg
     &    ,vmo_sflg, tax_sflg, tay_sflg ,tsu_sflg, tsc_sflg, tb2_sflg
     &    ,tb3_sflg, tgg_sflg, tgf_sflg ,thd_sflg, tld_sflg, thg_sflg
     &    ,tlg_sflg, thf_sflg, tlf_sflg ,thm_sflg, tlm_sflg, dtm_sflg
     &    ,rsv_sflg
     &    ,cld_sflg, cll_sflg, clm_sflg ,clh_sflg
     &    ,rgn_sflg, rgd_sflg, rgc_sflg ,sgn_sflg, sgd_sflg, sgc_sflg
     &    ,rtu_sflg, rtc_sflg, sot_sflg ,soc_sflg, als_sflg
     &    ,snd_sflg, sid_sflg, ico_sflg ,itf_sflg, isf_sflg, icu_sflg
     &    ,icv_sflg, div_sflg, gro_sflg ,ire_sflg, ich_sflg
     &    ,fwf_sflg
     &    ,sno_sflg, rev_sflg, ssb_sflg ,clc_sflg, lwp_sflg, pwc_sflg
     &    ,ref_sflg
     &    ,  u_sflg,   v_sflg,   t_sflg ,  q_sflg,  rh_sflg,   g_sflg
     &    ,  c_sflg,   l_sflg

      namelist /params/ kuoflag,ksc,kuocb,rhcut,rhkuo,rhsat

      namelist /scmpars/ nst63, lgt63, mgt63, monthset

      namelist /coupling/ crelax_flag, crelax_tau, fluxadj, isync,
     &                    ovolume, subice, hosing_flag, hosing_rate

c************************************************************************
c****  READ IN NAMELISTS - CHECK FOR INCONSISTENCIES ********************
c************************************************************************

      if(SCM)then
        read(5,scmpars)
        write(6,scmpars)
      endif

      read(5,control)
      write(6,control)

      read(5,diagnostics)
      write(6,diagnostics)

      read(5,statvars)
      write(6,statvars)

      read(5,histvars)
      write(6,histvars)

      read(5,params)
      write(6,params)

      if (lcouple .or. savefcor) then
        read (5, coupling)
        write (6, coupling)
      end if

      if(incd.ne.-1)then
        write(6,*)'Error: Replace obsolete INCD with NDSTOP in namelist'
        write(0,*)'Error: Replace obsolete INCD with NDSTOP in namelist'
        stop
      endif

      if ( saveicuv ) then
        write(6,*)'Error: Namelist flag saveicuv is obsolete - use stats
     &flag with icu_sflg and/or icv_sflg'
        write(0,*)'Error: Namelist flag saveicuv is obsolete - use stats
     &flag with icu_sflg and/or icv_sflg'
        stop
      endif

      if ( ksc.ge.0 ) then
        write(6,*)'Error: Namelist flag ksc is obsolete - It is now hard
     & coded in hvertmx for safety.'
        write(0,*)'Error: Namelist flag ksc is obsolete - It is now hard
     & coded in hvertmx for safety.'
        stop
      endif

      if (lbrmap) then
        write(6,*)'Warning: LBRMAP is obsolete and is no longer used'
      endif

c---- check for only one convection type
      kchk=0
      if(kuoflag)kchk=kchk+1
      if(ukconv)kchk=kchk+1
      if(kchk.gt.1)then
        print *,'Requesting more than 1 of ukconv,'
        print *,'or kuoflag.   Only one or none allowed'
        stop
      endif

c---- For NCEP starts for AGCM runs only only
      ncepstrt=-1
      if(ncepagcm)then
        if(lcouple.or.locean.or.qflux)then
          print *,'ncepagcm=',ncepagcm
          print *,'But 1 or more of following is true:'
          print *,'lcouple=',lcouple
          print *,'locean=',locean
          print *,'qflux=',qflux
          print *,'Are you sure? Model stopped in readnml1.f'
          stop
        endif
        ncepstrt=2
        chtst='old'
      endif

      ifwds=chtst.eq.'new'.or.chtst.eq.'NEW'

c---- J.S.Frederiksen diffusion only coded for T63 at moment
      if(mw.ne.64)jsfdiff=.false.

c...  Check that the values of MSTEP and NRAD are valid.
c...
c...  (1) NRAD = 0. This specifies no radiation, but MSTEP must still divide
c...                into 1440 minutes (= 1 day).
      if (nrad .eq. 0) then
        if(mod(1440,mstep).ne.0)then
          print *,'***  No radiation, but mstep does not divide into'
          print *,'***  the number of minutes per day (1440)'
          print *,'***  mstep=',mstep
          stop
        endif
      end if

c...  (2) NRAD > 0. This specifies the default radiation timestep of 120
c...                minutes (= 2 hours). MSTEP must therefore divide into 120
c...                minutes, and MSTEP*NRAD must be equal to 120.
      if (nrad .gt. 0) then
        if(mod(120,mstep).ne.0)then
          print *,'***  mstep does not allow 2 hr radiation steps'
          print *,'***  mstep=',mstep
          stop
        endif
        if(nrad*mstep.ne.120)then
          print *,'***  nrad not correct for radiation ever 2hrs'
          print *,'***  Changing the value of nrad to ',120/mstep
          nrad = 120/mstep
        endif
      end if

c...  (3) NRAD < 0. This specifies a radiation timestep other than 120 minutes
c...                (= 2 hours). However, the radiation timestep must still
c...                divide into 1440 minutes (= 1 day).
      if (nrad .lt. 0) then
        nrad=-nrad
        if(mod(1440,nrad*mstep).ne.0)then
          print *,'***  nrad (negative) for radiation every ',nrad
          print *,'***  steps of ',mstep,' each. This does not match'
          print *,'***  1440 mins per day.'
          stop
        endif
      end if
        
c---- To use the new river routing scheme, the resolution must be T63
      if(newriver)then
        if(mw.ne.64)then
          print *,'newriver = T, but resolution not T63'
          stop
        endif
      endif

      return
      end
