c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c $Log: surfdiag.f,v $
c Revision 1.12  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.11  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.10  1998/12/10  00:55:28  ldr
c HBG changes to V5-1-21
c
c Revision 1.9  1997/12/19  02:03:09  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.8  1993/10/07  12:09:56  ldr
c Move machine parameter statement into new include file MACHINE.f.
c
c Revision 1.7  93/08/19  15:10:57  ldr
c Minor cosmetic changes.
c 
c Revision 1.6  93/08/19  11:22:32  ldr
c Put lock around write statement on SGI.
c 
c Revision 1.5  93/07/27  14:16:47  mrd
c Added extra sib fields to surfdiag diagnostics
c 
c Revision 1.4  92/12/09  14:44:34  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  91/05/17  14:43:35  mrd
c *** empty log message ***
c 
c Revision 1.2  91/03/13  13:00:44  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:38:10  ldr
c Initial release V3-0
c 
      subroutine surfdiag(ip,mins,ns,lg,mg,ma,tg,tb2,tb3,snowd,
     &              siced,
     &              precs,precc,eg,fg,sg,sint,sout,soutclr,rg,rt,
     &              rtclr,alb,cl,cll,clm,clh,pg,ttg,qtg,u,v,cduv,tscrn,
     &              runoff,cint,totpev,scalev)
c  IP      Index of point
c  MINS    Current step
c  INS     Hemisphere index
c  LG      Latitude index
c  MG      Longitude index
c  MA      Longitude index for combined NH & SH arrays
c  TG      Surface temperature
c  TB2     Lower level ground temperature
c  TB3     Lower level ground temperature
c  WG      Upper level ground wetness
c  WB      Lower level ground wetness
c  SNOWD   Snow depth
c  SICED   Sea-ice thickness
c  PRECS   Large scale precipitation
c  PRECC   Convective precipitation
c  EG      Latent heat flux
c  FG      Sensible heat flux
c  SG      Solar absorbed in ground
c  SINT    Solar at top
c  SOUT    Solar out at top
c  SOUTCLR Clear sky solar out at top
c  RG      Net long wave at ground
c  RT      Long wave at top
c  RTCLR   Clear sky long wave at top
c  ALB     Surface albedo
c  CL      Total cloud
c  CLL     Low level cloud
c  CLM     Medium level cloud
c  CLH     High level cloud
c  PG      Surface pressure
c  T1      Level 1 atmospheric temperature
c  Q1      Level 1 mixing ratio
c  U1      Level 1 zonal wind
c  V1      Level 1 meridional wind
c  CDUV    Momentum drag coefficient
c  TSCRN   Temperature at screen level (1.5m)
c  RUNOFF  Runoff
c  CINT    Caonpy interception
c  TGG     Ground surface temperature
c  TGF     Canopy temperature
c  TOTPEV  Potential evaporation
c  SCALEV  Scaling evaporation

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ip
      integer*8 mins
      integer ns
      integer lg
      integer mg
      integer ma
      real tg(ln2)
      real tb2(ln2)
      real tb3(ln2)
      real snowd(ln2)
      real siced(ln2)
      real precs(ln2)
      real precc(ln2)
      real eg(ln2)
      real fg(ln2)
      real sg(ln2)
      real sint(ln2)
      real sout(ln2)
      real soutclr(ln2)
      real rg(ln2)
      real rt(ln2)
      real rtclr(ln2)
      real alb(ln2)
      real cl(ln2)
      real cll(ln2)
      real clm(ln2)
      real clh(ln2)
      real pg(ln2)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real u(ln2,nl)
      real v(ln2,nl)
      real cduv(ln2)
      real tscrn(ln2)
      real runoff(ln2)
      real cint(ln2)
      real totpev(ln2)
      real scalev(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'SURF1.f'

C Local work arrays and variables
      integer k

C Local data, functions etc

C Start code : ----------------------------------------------------------

      write(60,1000) ip,mins,ns,lg,mg,
     &     tg(ma),tb2(ma),tb3(ma),wb(ma,ms,lg),snowd(ma),
     &     siced(ma),
     &     precs(ma),precc(ma),eg(ma),fg(ma),sg(ma),sint(ma),
     &     sout(ma),soutclr(ma),rg(ma),rt(ma),rtclr(ma),alb(ma),
     &     cl(ma),cll(ma),clm(ma),clh(ma),pg(ma),
     &     (ttg(ma,k),k=1,2),(qtg(ma,k),k=1,1),
     &     (u(ma,k),k=1,1),(v(ma,k),k=1,1),cduv(ma),tscrn(ma),
     &     runoff(ma),cint(ma),tggsl(ma,1,lg),tgf(ma,lg),
     &     totpev(ma),scalev(ma)
 1000 format(i3,i8,3i3,/,(6e12.6))
      return
      end
