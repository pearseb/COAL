c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: radfs.f,v $
c Revision 1.60  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.59  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.58.1.1  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.58  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.57  2001/02/12 05:39:48  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.56  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.55  1999/06/30 05:29:37  rot032
c Mods for SCM.
c
c Revision 1.54  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.53  1998/12/10  00:55:37  ldr
c HBG changes to V5-1-21
c
c Revision 1.52  1997/12/17  23:22:47  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.51  1997/11/27  05:35:00  ldr
c Flexible treatment of direct and indirect aerosol effects.
c
c Revision 1.50  1997/06/11  02:21:29  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.49  1996/11/22  05:13:05  ldr
c Comments and tidy-ups from LDR.
c
c Revision 1.48  1996/06/13  02:07:42  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.47  1996/02/19  04:09:57  ldr
c Generalize for 24 levels.
c
c Revision 1.46  1995/11/23  06:03:31  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.45  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.44  1995/08/08  02:02:17  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.43  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.36.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.42  1994/12/08  00:21:56  mrd
c Add declarations to work with implicit none.
c
c Revision 1.41  1994/09/12  12:50:16  ldr
c Changes to qcloud scheme of V4-6 to create run E32, except that E32 has
c enhanced collection of cloud water by falling ice.
c
c Revision 1.40  94/08/09  11:20:52  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.39  94/08/08  16:15:00  ldr
c Add clat as an argument to first cloud2 call also.
c 
c Revision 1.38  94/08/08  13:11:43  ldr
c Pass back clat from cloud2 as total cloud fraction in box.
c 
c Revision 1.37  94/08/04  16:56:20  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.36  94/06/15  17:18:27  mrd
c Moved Hal's extrapolation of the top level temperature from radin to radfs
c to make the radiation code more modular.
c 
c Revision 1.35  94/05/03  11:57:22  ldr
c Simplify conversion of fluxes to heating rates (from MRD).
c 
c Revision 1.34  93/12/17  15:33:34  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.33  93/11/23  16:45:46  mrd
c Removed unused csolar parameter (moved to initfs).
c 
c Revision 1.32  93/11/03  11:32:04  ldr
c Let the radiation scheme see the new progrnostic clouds.
c 
c  This routine gives the interface between the RADIN routine and the
c  Fels - Schwarzkopf radiation code.
c
c  Inputs
c  PG 	     Surface pressure (mb)
c  TG        Surface temperature
c  TTG       Temperature
c  QTG       Mixing Ratio
c  RHG       Relative humidity
c  ALBEDO    Surface albedo
c  SIG       Sigma levels
c
c  Outputs
c  RG        Net long wave heating at the ground (i.e., including
c            sigma Ts^4 term)
c  RT        Long wave at the top of the atmosphere
c  RTCLR     Clear sky long wave at the top of the atmosphere
c  RGCLR     Clear sky long wave at the surface
c  SG        Short wave at the ground.
c  SINT      Solar in at top
c  SOUT      Solar out at top
c  SOUTCLR   Clear sky solar out at top
c  SGCLR     Clear sky solar at surface
c  HTK       Total atmospheric heating
c
      subroutine radfs(ipass,lg,pg,tg,z4,ttg,qtg,rhg,albedo,
     &             cfrac,bvnf,qlg,qfg,qccon,cdso4,
     &             rg,rt,rtclr,rgclr,
     &             sg,sint,sout,soutclr,sgclr,htk,cll,clm,clh,clcon,
     &             sga,clat,refflm,cldliq)

      implicit none

!$OMP THREADPRIVATE ( /CLDCOM/ )
!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /KDACOM/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /LWOUT/ )
!$OMP THREADPRIVATE ( /RADISW/ )
!$OMP THREADPRIVATE ( /RDFLUX/ )
!$OMP THREADPRIVATE ( /CLRFLX/ )
!$OMP THREADPRIVATE ( /SRCCOM/ )
!$OMP THREADPRIVATE ( /SWOCOM/ )
!$OMP THREADPRIVATE ( /TFCOM/ )

C Global parameters
      include 'HCON.f'
      include 'PARAMS.f'
      include 'RDPARM.f'
c     parameters for the aerosol calculation
      real beta_ave, alpha
c      parameter(beta_ave = 0.29, alpha = 8.00)
      parameter(beta_ave = 0.233, alpha = 5.44) !Smaller forcing

C Argument list
      integer ipass,lg
      real pg(ln2)
      real tg(ln2)
      real z4(ln2)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real rhg(ln2,nl)
      real albedo(ln2)
      real cfrac(ln2,nl)
      real bvnf(ln2,nl)
      real qlg(ln2,nl) !Cloud liq water
      real qfg(ln2,nl) !Cloud ice
      real qccon(ln2,nl) !Conv cloud water mixing ratio
      real cdso4(ln2,nl)
      real rg(ln2)
      real rt(ln2)
      real rtclr(ln2)
      real rgclr(ln2)
      real sg(ln2)
      real sint(ln2)
      real sout(ln2)
      real soutclr(ln2)
      real sgclr(ln2)
      real htk(ln2,nl)
      real cll(ln2)
      real clm(ln2)
      real clh(ln2)
      real clcon(ln2,nl) !Conv cloud amount
      real sga(ln2)
      real clat(ln2,nl)  !Cloud amount diagnostic
      real refflm(ln2)
      real cldliq(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'CLDCOM.f'
      include 'HYBRPR.f'
      include 'KDACOM.f'
      include 'LOGIMSL.f'
      include 'LWOUT.f'
      include 'RADISW.f'
      include 'RDFLUX.f' ! includes CLRFLX
      include 'SRCCOM.f'
      include 'SWOCOM.f'
      include 'TFCOM.f'

C Global data blocks
      include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'SULFATE1.f'
      include 'TIMEX.f'

      real alat,ha,coszm,taudam,fjd,r1,dlt,slag,dhr
      common /radsav/ alat(lat*2),ha(lat*2),coszm(lat*2),taudam(lat*2),
     &      fjd,r1,dlt,slag,dhr

      integer nst63,lgt63,mgt63,monthset
      common / scm_spacetime / nst63,lgt63,mgt63,monthset !Needed for SCM version only

C Local work arrays and variables
c  Next 2 variables dimensioned up for SCM version
      real coszen(ilon*2),tauda2(ilon*2)
      real duo3n(nl,2)

      integer ip
      integer ilatn
      integer ilats
      integer k
      integer kr
      integer mg
      integer ns

      real atxn ! Constant for extrapolating top level temperature
      real btxn ! Constant  "
      real cosz
      real delta
      real fjdr

      logical cldoff

C Local data, functions etc
      logical solarfit
      save solarfit
      data solarfit /.true./

C Start code : ----------------------------------------------------------

c     Initialization stuff has been moved to initfs.
c     Other stuff (e.g. call zenith) has been moved to radin.      

      ilatn=lg
      ilats=lat2+1-lg

c     Time (hours) between radiation calculations
      dhr=float(nrad*mstep)/60.
c Due to radfs called "last" in radin, the time indicator (fjdr) must 
c      have float(mstep)/1440 added
      fjdr=fjd+float(mstep)/1440. ! Radiation for next "nrad" period
      call zenith(fjdr,r1,dlt,slag,alat(ilatn),ha(ilatn),dhr,ilon,
     &            coszen(1),tauda2(1))
      call zenith(fjdr,r1,dlt,slag,alat(ilats),ha(ilats),dhr,ilon,
     &            coszen(1+ilon),tauda2(1+ilon))

      if(SCM)then
        if(nst63.eq.1)then
          coszro(1)=coszen(mgt63)
          taudar(1)=tauda2(mgt63)
        else
          coszro(1)=coszen(mgt63+ilon)
          taudar(1)=tauda2(mgt63+ilon)
        endif
      else
        do ip=1,ln2
          coszro(ip)=coszen(ip)
          taudar(ip)=tauda2(ip)
        enddo
      endif

      if(amipo3)then
        call o3set_amip ( alat(ilatn), mins, prh, 1, qo3 )
        call o3set_amip ( alat(ilats), mins, prh, 2, qo3 )
      else
c     Set up ozone for this time and latitude
        call o3set(alat(ilatn),mins,duo3n(1,1),sig)
        call o3set(alat(ilats),mins,duo3n(1,2),sig)
c         Conversion of o3 from units of cm stp to gm/gm
        do 215 k=1,nl
        do 215 ip=1,lon
          qo3(ip,k) = duo3n(k,1)*1.01325e+02/(pg(ip)*1.e3)
          qo3(ip+lon,k) = duo3n(k,2)*1.01325e+02/(pg(ip+lon)*1.e3)
  215   continue
      endif

c     Constants for extrapolating top level temperature
      atxn=sig(lm1)/(sig(lm1)-sig(l))
      btxn=sig(l)/(sig(lm1)-sig(l))

c     Set up basic variables, reversing the order of the vertical levels
      do 196 ip=1,ln2
  196 temp(ip,lp1) = tg(ip)
      do 198 ip=1,ln2
      if(sea(ip).or.mlo(ip).or.((ipass.eq.2).and.cice(ip)))then
        temp(ip,lp1) = tg(ip)-z4(ip)*0.0065
      end if
  198 continue
      do 210 k=1,nl
        do 210 ip=1,ln2
          kr=nl+1-k
          press(ip,kr)=prf(ip,k)       * 1000. ! mb converted to cgs
          temp(ip,kr)=ttg(ip,k)
c         Set min value to avoid numerical problems
c         rh2o(ip,kr)=max(qtg(ip,k) ,1.e-7)
  210     rh2o(ip,kr)=max(qtg(ip,k) ,0.3e-05)

      do 200 ip=1,ln2
c       Hal's scheme for extrapolating the top level temperature
        temp(ip,1)=atxn*temp(ip,1)-btxn*temp(ip,2)+0.0
        
        press(ip,lp1) = pg(ip) * 1000.
ch      temp(ip,lp1) = tg(ip)
        cuvrf(ip,1)=albedo(ip)   ! surface albedo
        cirrf(ip,1)=albedo(ip)
        cirab(ip,1)=zero
c       Calculate half level pressures and temperatures by linear interp
        do 220 k=1,nlp
	  kr=nlp+1-k
          press2(ip,kr)=prh(ip,k)        * 1000. ! mb converted to cgs
  220   continue
        do 230 k=2,nl
	  kr=nlp+1-k
	  temp2(ip,kr) =
     &    (temp(ip,kr)*(prf(ip,k)-prh(ip,k))
     &    +temp(ip,kr-1)*(prh(ip,k)-prf(ip,k-1)))/
     &    (prf(ip,k)-prf(ip,k-1))
c     &    (ttg(ip,k-1)*(prf(ip,k)-prh(ip,k))
c     &    +ttg(ip,k)*(prh(ip,k)-prf(ip,k-1)))/
c     &    (prf(ip,k)-prf(ip,k-1))
  230   continue
        temp2(ip,1) = temp(ip,1)
	temp2(ip,lp1) = temp(ip,lp1)
  200 continue

c Apply sulfate modification to the albedo

      if(naerosol_d.gt.0)then
        do ip=1,ln2
          cosz = max ( coszro(ip), 1.e-10)
          delta =  beta_ave*alpha*so4dir(ip,lg)*
     &                                 (1-albedo(ip))**2/cosz
          cuvrf(ip,1)=min(1., delta+albedo(ip)) ! surface albedo
          cirrf(ip,1)=cuvrf(ip,1)
        enddo
      endif

c  Clear sky calculation
      if (clforflag) then
        cldoff=.true.
c     set up cloud for this time and latitude
        if ( qcloud ) then
          call cloud2(cldoff,lg,ttg,qlg,qfg,cfrac,
     &                clcon,qccon,cdso4,                              !Inputs
     &                clat,cll,clm,clh,refflm,cldliq)                 !Outputs
        else
          call cloud(ipass,lg,alat(ilatn),alat(ilats)
     &       ,cll,clm,clh,cldoff,rhg,clcon,bvnf,clat)
        endif
        call swr89(fsw,hsw,sg,ufsw,dfsw,press,press2,coszro,
     &       taudar,rh2o,rrco2,ssolar,qo3,nclds,
     &       ktopsw,kbtmsw,cirab,cirrf,cuvrf,camt)
        do ip=1,ln2
          soutclr(ip) = ufsw(ip,1)*h1m3 ! solar out top
	  sgclr(ip)   = -fsw(ip,lp1)*h1m3  ! solar absorbed at the surface
        end do
      else
        do ip=1,ln2
          soutclr(ip) = 0.
	  sgclr(ip) = 0.
        end do
      endif

c  Cloudy sky calculation
      cldoff=.false.
      if ( qcloud ) then
        call cloud2(cldoff,lg,ttg,qlg,qfg,cfrac,
     &                clcon,qccon,cdso4,                              !Inputs
     &             clat,cll,clm,clh,refflm,cldliq)                    !Outputs
      else
        call cloud(ipass,lg,alat(ilatn),alat(ilats)
     &     ,cll,clm,clh,cldoff,rhg,clcon,bvnf,clat)
      endif
      call swr89(fsw,hsw,sg,ufsw,dfsw,press,press2,coszro,
     &                   taudar,rh2o,rrco2,ssolar,qo3,nclds,
     &                   ktopsw,kbtmsw,cirab,cirrf,cuvrf,camt)
      do ip=1,ln2
	  sint(ip) = dfsw(ip,1)*h1m3   ! solar in top
	  sout(ip) = ufsw(ip,1)*h1m3   ! solar out top
	  sg(ip)   = sg(ip)*h1m3       ! solar absorbed at the surface
      end do

      call clo89
      call lwr88

c     set up the diagnostics to be passed back
      do 300 ip=1,ln2
        rt(ip) = ( gxcts(ip)+flx1e1(ip) ) * h1m3    ! longwave at top
c       clear sky longwave at the top
        rtclr(ip) = ( gxctsclr(ip)+flx1e1clr(ip) ) * h1m3
        rg(ip) = grnflx(ip)*h1m3   ! longwave at surface
	rgclr(ip) = grnflxclr(ip)*h1m3   ! clear sky longwave at surface
  300 continue
      do 310 k=1,nl
      do 310 ip=1,ln2
c         total heating rate
c---- note : htk now in Watts/M**2 (no pressure at level weighting)
c     Convert from cgs to SI units
          htk(ip,nl+1-k) = 0.001*(hsw(ip,k)+heatra(ip,k))
  310 continue

      if ( solarfit ) then
c       Calculate the amplitude of the diurnal cycle of solar radiation
c       at the surface (using the value for the middle of the radiation
c       step) and use this value to get solar radiation at other times.
c       Use the zenith angle and daylight fraction calculated in zenith
c       to remove these factors.

        do ip=1,ln2
           if ( taudar(ip) .le. 0. ) then
c             The sun isn't up at all over the radiation period so no 
c             fitting need be done.
              sga(ip) = 0.
           else
              sga(ip) = sg(ip) / (coszro(ip)*taudar(ip))
           end if
	end do
      else
	do ip=1,ln2
	  sga(ip) = 0.
	end do
      end if

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,i1)')'After radfs. IPASS = ',ipass
          write(25,91)'clat ',(clat(mg,nl+1-k),k=1,nl)
          write(25,91)'clcon',(clcon(mg,k),k=1,nl)
          write(25,1)'cll ',cll(mg),' clm ',clm(mg),' clh ',clh(mg)
          write(25,1)'sg ',sg(mg),' rg ',rg(mg),' zen ',acos(coszro(mg))
     &                                                  *57.296
          write(25,1)'sint ',sint(mg),' sout ',sout(mg),' soutclr ',
     &         soutclr(mg)
          write(25,91)'htk ',(htk(mg,k),k=1,nl)
          write(25,91)'hsw ',(.001*hsw(mg,nl+1-k),k=1,nl)
          write(25,91)'hlw ',(.001*heatra(mg,nl+1-k),k=1,nl)
          write(25,91)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qccon ',(qccon(mg,k),k=1,nl)
          write(25,*)
        endif
      endif
 1    format(3(a,f10.3))
 9    format(a,30g10.3)
 91   format(a,30f10.3)


      return
      end
