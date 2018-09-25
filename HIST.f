c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c $Log: HIST.f,v $
c Revision 1.12  2001/02/12 05:39:48  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.11  2000/11/14 03:11:36  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.10  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.9  1995/10/04  06:38:16  mrd
c Changes for more efficient writing of history files.
c
c Revision 1.8  1994/09/13  09:51:17  mrd
c Add daily mean 10m windspeed and 2m relative humidity as history variables.
c
c Revision 1.7  94/09/09  14:53:44  mrd
c Added daily average cloud, surface temperature and screen temperature
c as history variables.
c 
c Revision 1.6  94/09/09  14:10:28  mrd
c Added possbility of saving two history files at different frequencies.
c 
c Revision 1.5  94/06/28  16:23:38  mrd
c Add sea-level pressure as a history variable.
c 
c Revision 1.4  93/10/08  09:14:47  mrd
c Added maximum and minimum daily surface temperature to history.
c 
c Revision 1.3  93/10/07  11:52:18  mrd
c Added flags to control writing of individual variables to daily history.
c 
c Revision 1.2  93/08/18  15:05:32  mrd
c Modified history to use netCDF
c 
c Revision 1.1  93/07/23  15:02:41  mrd
c Initial revision
c 
c This file defines the arrays used for the history archive
c Variables are arrays of size 2 to allow saving 2 history files at
c different frequencies.
      integer hist_interval(2), histid(2), histset(2)
      logical savehist(2)
c  Logical variables to control which fields are saved
      logical all_hflg(2),
     &        zht_hflg(2), t_hflg(2), u_hflg(2), v_hflg(2), q_hflg(2),
     &        psf_hflg(2), imsl_hflg(2),tsu_hflg(2), tb2_hflg(2),
     &        tb3_hflg(2), als_hflg(2), wfg_hflg(2), wfb_hflg(2),
     &        snd_hflg(2), sid_hflg(2), rnd_hflg(2), rnc_hflg(2),
     &        evp_hflg(2), hfl_hflg(2), sgn_hflg(2), sgc_hflg(2),
     &        sgd_hflg(2), sit_hflg(2), sot_hflg(2), soc_hflg(2),
     &        rgn_hflg(2), rgc_hflg(2), rgd_hflg(2), rtu_hflg(2),
     &        rtc_hflg(2), tax_hflg(2), tay_hflg(2), cld_hflg(2),
     &        cll_hflg(2), clm_hflg(2), clh_hflg(2), ich_hflg(2),
     &        icm_hflg(2), ict_hflg(2), icb_hflg(2), tsc_hflg(2),
     &        tsh_hflg(2), tsl_hflg(2), run_hflg(2), int_hflg(2),
     &        tgg_hflg(2), tgf_hflg(2), pev_hflg(2), sev_hflg(2),
     &        ico_hflg(2), tgl_hflg(2), tgh_hflg(2), psl_hflg(2),
     &        tsca_hflg(2), tsua_hflg(2), clda_hflg(2), rhsa_hflg(2),
     &        v10ma_hflg(2)
      integer tid(nl,2),uid(nl,2),vid(nl,2),qid(nl,2),pid(2),plid(2)
      integer cldid(2),cllid(2),clmid(2),clhid(2),ichid(2),icmid(2),
     &        icbid(2),ictid(2),alsid(2)
      common /hist_control/ savehist,hist_interval,histid,histset,
     &       tid,uid,vid,qid,pid,plid,cldid,cllid,clmid,clhid,
     &       ichid,icmid,icbid,ictid,alsid,
     &       all_hflg, zht_hflg, t_hflg, u_hflg, v_hflg, q_hflg,
     &       psf_hflg, imsl_hflg, tsu_hflg, tb2_hflg, tb3_hflg,
     &       als_hflg, wfg_hflg, wfb_hflg, snd_hflg, sid_hflg,
     &       rnd_hflg, rnc_hflg, evp_hflg, hfl_hflg, sgn_hflg,
     &       sgc_hflg, sgd_hflg, sit_hflg, sot_hflg, soc_hflg,
     &       rgn_hflg, rgc_hflg, rgd_hflg, rtu_hflg, rtc_hflg,
     &       tax_hflg, tay_hflg, cld_hflg, cll_hflg, clm_hflg,
     &       clh_hflg, ich_hflg, icm_hflg, ict_hflg, icb_hflg,
     &       tsc_hflg, tsh_hflg, tsl_hflg, run_hflg, int_hflg,
     &       tgg_hflg, tgf_hflg, pev_hflg, sev_hflg, ico_hflg,
     &       tgl_hflg, tgh_hflg, psl_hflg, tsca_hflg, tsua_hflg,
     &       clda_hflg, rhsa_hflg, v10ma_hflg

c     Archive arrays for DARLAM
      integer nhcomp ! Use this to save less data from eg T63 Mk3 runs
c....  nhcomp = 1 ==> All grid ponts saved
c....  nhcomp = 2 ==> Every other grid point saved  (T63 down to sort of T32) 
c....  nhcomp = 3 ==> Every third grid point saved  (T63 down to sort of T21)
c
      parameter (nhcomp=1) ! For T63 resolution data for DARLAM
c     parameter (nhcomp=2) ! For T63 data saved as T32 for DARLAM
c     parameter (nhcomp=3) ! For T63 data saved as RT1 for DARLAM
c....
      integer nhist,lonnh,lat2nh
      parameter (nhist=1)  !Memory saving trick - set nhist=1 to use savehist
      parameter (lonnh=nhist*(lon-1)+1)   ! =  lon if nhist=1
      parameter (lat2nh=nhist*(lat2-1)+1) ! = lat2 if nhist=1
      real rain_a   (lonnh,lat2nh)
      real precc_a  (lonnh,lat2nh)
      real evap_a   (lonnh,lat2nh)
      real fg_a     (lonnh,lat2nh)
      real sg_a     (lonnh,lat2nh)
      real sgclr_a  (lonnh,lat2nh)
      real sgdn_a   (lonnh,lat2nh)
      real sint_a   (lonnh,lat2nh)
      real sout_a   (lonnh,lat2nh)
      real soutclr_a(lonnh,lat2nh)
      real rg_a     (lonnh,lat2nh)
      real rgclr_a  (lonnh,lat2nh)
      real rgdn_a   (lonnh,lat2nh)
      real rt_a     (lonnh,lat2nh)
      real rtclr_a  (lonnh,lat2nh)
      real taux_a   (lonnh,lat2nh)
      real tauy_a   (lonnh,lat2nh)
      real runoff_a (lonnh,lat2nh)
      real cint_a   (lonnh,lat2nh)
      real pev_a    (lonnh,lat2nh)
      real sev_a    (lonnh,lat2nh)
      real tscrn_sav(lonnh,lat2nh)
      real tscrn_max(lonnh,lat2nh)
      real tscrn_min(lonnh,lat2nh)
      real tscrn_a  (lonnh,lat2nh)
      real tg_max   (lonnh,lat2nh)
      real tg_min   (lonnh,lat2nh)
      real tg_a     (lonnh,lat2nh)
      real cld_a    (lonnh,lat2nh)
      real rhscrn_a (lonnh,lat2nh)
      real v10m_a   (lonnh,lat2nh)
      real tgrid    (lonnh,lat2nh,nl)
      real pslgrid  (lonnh,lat2nh)
      real cldgrid  (lonnh,lat2nh)
      real cllgrid  (lonnh,lat2nh)
      real clmgrid  (lonnh,lat2nh)
      real clhgrid  (lonnh,lat2nh)
      real alsgrid  (lonnh,lat2nh)
      integer ichgrid(lonnh,lat2nh)
      integer icmgrid(lonnh,lat2nh)
      integer icbgrid(lonnh,lat2nh)
      integer ictgrid(lonnh,lat2nh)
      common /hist/ rain_a,precc_a,evap_a,fg_a,sg_a,sgclr_a,sgdn_a,
     &              sint_a,sout_a,soutclr_a,rg_a,rgclr_a,rgdn_a,
     &              rt_a,rtclr_a,taux_a,tauy_a,runoff_a,cint_a,
     &              pev_a,sev_a,tscrn_sav,tscrn_max,tscrn_min,
     &              tscrn_a,tg_max,tg_min,tg_a,cld_a,rhscrn_a,v10m_a,
     &              tgrid,pslgrid,cldgrid,cllgrid,clmgrid,clhgrid,
     &              alsgrid,ichgrid,icmgrid,icbgrid,ictgrid
