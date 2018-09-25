c $Log: STFLAGS.f,v $
c Revision 1.11  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.10  2001/02/12 05:39:51  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.9  2000/11/14 03:11:37  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.8  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.7  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.6  1997/12/19  02:03:15  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.5  1997/06/11  02:21:29  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.4  1996/04/17  01:32:52  ldr
c Changes from TIE to allow initialization of monthly mean netcdf files
c from within model.
c
c Revision 1.3  1994/08/08  13:12:49  ldr
c Add flag c_sflg for cloud stats.
c
c Revision 1.2  93/11/29  14:49:21  ldr
c Use icu_sflg and icv_sflg for treatment of ice U,V diagnostics, so they
c are consistent with other "stats".
c 
c Revision 1.1  93/10/08  10:11:29  ldr
c Initial revision
c 

C Surface statistics flags
      logical 
     &  evp_sflg, pev_sflg, sev_sflg ,rnd_sflg, rnc_sflg, hfl_sflg
     & ,wfg_sflg, wfb_sflg, run_sflg ,per_sflg, int_sflg, psl_sflg
     & ,vmo_sflg, tax_sflg, tay_sflg ,tsu_sflg, tsc_sflg, tb2_sflg
     & ,tb3_sflg, tgg_sflg, tgf_sflg ,thd_sflg, tld_sflg, thg_sflg
     & ,tlg_sflg, thf_sflg, tlf_sflg ,thm_sflg, tlm_sflg, dtm_sflg
     & ,rsv_sflg
C Cloud statistics
      logical 
     &  cld_sflg, cll_sflg, clm_sflg ,clh_sflg
C Radiation statistics
      logical 
     &  rgn_sflg, rgd_sflg, rgc_sflg ,sgn_sflg, sgd_sflg, sgc_sflg
     & ,rtu_sflg, rtc_sflg, sot_sflg ,soc_sflg, als_sflg
C Snow & ice statistics
      logical 
     &  snd_sflg, sid_sflg, ico_sflg ,itf_sflg, isf_sflg, icu_sflg
     & ,icv_sflg, div_sflg, gro_sflg ,ire_sflg, ich_sflg
C Fresh water fllux into ocean statistics
      logical fwf_sflg
C Qcloud statistics
      logical 
     &  sno_sflg, rev_sflg, ssb_sflg ,clc_sflg, lwp_sflg, pwc_sflg
     & ,ref_sflg
C Model level statistics
      logical 
     &    u_sflg,   v_sflg,   t_sflg ,  q_sflg,  rh_sflg,   g_sflg
     & ,  c_sflg,   l_sflg

      common/stflags/
     &  evp_sflg, pev_sflg, sev_sflg ,rnd_sflg, rnc_sflg, hfl_sflg
     & ,wfg_sflg, wfb_sflg, run_sflg ,per_sflg, int_sflg, psl_sflg
     & ,vmo_sflg, tax_sflg, tay_sflg ,tsu_sflg, tsc_sflg, tb2_sflg
     & ,tb3_sflg, tgg_sflg, tgf_sflg ,thd_sflg, tld_sflg, thg_sflg
     & ,tlg_sflg, thf_sflg, tlf_sflg ,thm_sflg, tlm_sflg, dtm_sflg
     & ,rsv_sflg
     & ,cld_sflg, cll_sflg, clm_sflg ,clh_sflg
     & ,rgn_sflg, rgd_sflg, rgc_sflg ,sgn_sflg, sgd_sflg, sgc_sflg
     & ,rtu_sflg, rtc_sflg, sot_sflg ,soc_sflg, als_sflg
     & ,snd_sflg, sid_sflg, ico_sflg ,itf_sflg, isf_sflg, icu_sflg
     & ,icv_sflg, div_sflg, gro_sflg ,ire_sflg, ich_sflg
     & ,fwf_sflg
     & ,sno_sflg, rev_sflg, ssb_sflg ,clc_sflg, lwp_sflg, pwc_sflg
     & ,ref_sflg
     & ,  u_sflg,   v_sflg,   t_sflg ,  q_sflg,  rh_sflg,   g_sflg
     & ,  c_sflg,   l_sflg
