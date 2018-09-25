c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c $Log: hist_acc.f,v $
c Revision 1.9  2000/11/14 03:11:39  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.8  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.7  1998/12/10  00:55:52  ldr
c HBG changes to V5-1-21
c
c Revision 1.6  1996/06/13  02:05:35  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.5  1994/09/13  09:51:27  mrd
c Add daily mean 10m windspeed and 2m relative humidity as history variables.
c
c Revision 1.4  94/09/09  14:53:56  mrd
c Added daily average cloud, surface temperature and screen temperature
c as history variables.
c 
c Revision 1.3  94/09/09  14:14:38  mrd
c Added possbility of saving two history files at different frequencies.
c 
c Revision 1.2  93/10/08  09:14:46  mrd
c Added maximum and minimum daily surface temperature to history.
c 
c Revision 1.1  93/08/18  15:06:05  mrd
c Initial revision
c 
c 
      subroutine hist_acc(mins,lg,radstep,
     &           rain,precc,evap,fg,sg,sgclr,sgdn,sint,sout,
     &           soutclr,rg,rgclr,rgdn,rt,rtclr,tscrn,tg,
     &           taux,tauy,runoff,cint,pev,sev,cld,rhscrn,v10m)

c     Accumulate averages for the history.
c     hist_acc is only called for the first history set because accumulated
c     variables may not be saved in the second.

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer*8 mins	    ! Current time
      integer lg            ! Latitude index
      logical radstep
      real rain     (ln2)   ! Total precipitation
      real precc    (ln2)   ! Convective precipitation
      real evap     (ln2)   ! Evaporation
      real fg       (ln2)   ! Sensible heat flux
      real sg       (ln2)   ! Solar absorbed in ground
      real sgclr    (ln2)   ! Clear sky solar absorbed in ground
      real sgdn     (ln2)   ! Solar downwards at surface
      real sint     (ln2)   ! Solar at top
      real sout     (ln2)   ! Solar out at top
      real soutclr  (ln2)   ! Clear sky solar out at top
      real rg       (ln2)   ! Net long wave at ground
      real rgclr    (ln2)   ! Clear sky net LW at ground
      real rgdn     (ln2)   ! LW downwards at surface
      real rt       (ln2)   ! Long wave at top
      real rtclr    (ln2)   ! Clear sky long wave at top
      real tscrn    (ln2)   ! Temperature at screen level
      real tg       (ln2)   ! Surface temperature
      real taux     (ln2)   ! Zonal component of wind stress
      real tauy     (ln2)   ! Meridional component of wind stress
      real runoff   (ln2)   ! Runoff
      real cint     (ln2)   ! Canopy interception
      real pev      (ln2)   ! Potential evaporation
      real sev      (ln2)   ! Scaling evaporation
      real cld      (ln2)   ! Total cloud
      real rhscrn   (ln2)   ! Screen level RH
      real v10m     (ln2)   ! 10m wind speed

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'HIST.f'

C Local work arrays and variables
      integer lgns
      integer ma
      integer mg
      integer ns

C Local data, functions etc

C Start code : ----------------------------------------------------------

c     If a history archive interval has just passed then initialise all
c     the accumulation arrays.

      do ns=1,2

      lgns = lg*(ns-1)+(lat2p-lg)*(2-ns)
      if( mod(mins,int(hist_interval(1), 8)).eq.0_8) then
	 do mg=1,lon
            ma=mg+(ns-1)*lon
	    rain_a(mg,lgns)    = rain(ma) - 0.5*(evap(ma)-abs(evap(ma)))
	    precc_a(mg,lgns)   = precc(ma)
	    evap_a(mg,lgns)    = 0.5*(evap(ma)+abs(evap(ma)))
	    fg_a(mg,lgns)      = fg(ma)
	    sg_a(mg,lgns)      = sg(ma)
	    sgclr_a(mg,lgns)   = sgclr(ma)
	    sgdn_a(mg,lgns)    = sgdn(ma)
	    sint_a(mg,lgns)    = sint(ma)
	    sout_a(mg,lgns)    = sout(ma)
	    soutclr_a(mg,lgns) = soutclr(ma)
	    rg_a(mg,lgns)      = rg(ma)
	    rgclr_a(mg,lgns)   = rgclr(ma)
	    rgdn_a(mg,lgns)    = rgdn(ma)
	    rt_a(mg,lgns)      = rt(ma)
	    rtclr_a(mg,lgns)   = rtclr(ma)
	    taux_a(mg,lgns)    = taux(ma)
	    tauy_a(mg,lgns)    = tauy(ma)
	    runoff_a(mg,lgns)  = runoff(ma)
	    cint_a(mg,lgns)    = cint(ma)
	    pev_a(mg,lgns)     = pev(ma)
	    sev_a(mg,lgns)     = sev(ma)
	    tscrn_sav(mg,lgns) = tscrn(ma)
	    tscrn_max(mg,lgns) = tscrn(ma)
	    tscrn_min(mg,lgns) = tscrn(ma)
            tscrn_a(mg,lgns)   = tscrn(ma)
            tg_max(mg,lgns)    = tg(ma)
            tg_min(mg,lgns)    = tg(ma)
            tg_a(mg,lgns)      = tg(ma)
            cld_a(mg,lgns)     = cld(ma)
            rhscrn_a(mg,lgns)  = rhscrn(ma)
            v10m_a(mg,lgns)    = v10m(ma)
	 end do
c        Note that this assumes that the history interval is a 
c        multiple of the radiation interval. Otherwise the radiative
c        quantities would have to be separately initialised
      else
c        Accumulation
	 do mg=1,lon
            ma=mg+(ns-1)*lon
	    rain_a(mg,lgns)    = rain_a(mg,lgns) + rain(ma) - 
     &                            0.5*(evap(ma)-abs(evap(ma)))
	    precc_a(mg,lgns)   = precc_a(mg,lgns)    + precc(ma)
	    evap_a(mg,lgns)    = evap_a(mg,lgns) + 
     &                            0.5*(evap(ma)+abs(evap(ma)))
	    fg_a(mg,lgns)      = fg_a(mg,lgns)       + fg(ma)
	    sg_a(mg,lgns)      = sg_a(mg,lgns)       + sg(ma)
	    rg_a(mg,lgns)      = rg_a(mg,lgns)       + rg(ma)
	    taux_a(mg,lgns)    = taux_a(mg,lgns)     + taux(ma)
	    tauy_a(mg,lgns)    = tauy_a(mg,lgns)     + tauy(ma)
	    runoff_a(mg,lgns)  = runoff_a(mg,lgns)   + runoff(ma)
	    cint_a(mg,lgns)    = cint_a(mg,lgns)     + cint(ma)
	    pev_a(mg,lgns)     = pev_a(mg,lgns)      + pev(ma)
	    sev_a(mg,lgns)     = sev_a(mg,lgns)      + sev(ma)
	    tscrn_sav(mg,lgns) = tscrn(ma)
	    tscrn_max(mg,lgns) = max(tscrn_max(mg,lgns), tscrn(ma))
	    tscrn_min(mg,lgns) = min(tscrn_min(mg,lgns), tscrn(ma))
	    tscrn_a(mg,lgns)   = tscrn_a(mg,lgns)    + tscrn(ma)
	    tg_max(mg,lgns)    = max(tg_max(mg,lgns), tg(ma))
	    tg_min(mg,lgns)    = min(tg_min(mg,lgns), tg(ma))
	    tg_a(mg,lgns)      = tg_a(mg,lgns)       + tg(ma)
	    rhscrn_a(mg,lgns)  = rhscrn_a(mg,lgns)   + rhscrn(ma)
	    v10m_a(mg,lgns)    = v10m_a(mg,lgns)     + v10m(ma)
	 end do
	 if ( radstep ) then
	   do mg=1,lon
              ma=mg+(ns-1)*lon
	      sgclr_a(mg,lgns)   = sgclr_a(mg,lgns)    + sgclr(ma)
	      sgdn_a(mg,lgns)    = sgdn_a(mg,lgns)     + sgdn(ma)
	      sint_a(mg,lgns)    = sint_a(mg,lgns)     + sint(ma)
	      sout_a(mg,lgns)    = sout_a(mg,lgns)     + sout(ma)
	      soutclr_a(mg,lgns) = soutclr_a(mg,lgns)  + soutclr(ma)
	      rgclr_a(mg,lgns)   = rgclr_a(mg,lgns)    + rgclr(ma)
	      rgdn_a(mg,lgns)    = rgdn_a(mg,lgns)     + rgdn(ma)
	      rt_a(mg,lgns)      = rt_a(mg,lgns)       + rt(ma)
	      rtclr_a(mg,lgns)   = rtclr_a(mg,lgns)    + rtclr(ma)
              cld_a(mg,lgns)     = cld_a(mg,lgns)      + cld(ma)
	   end do
	end if
      end if

      enddo ! ns=1,2

      return
      end
