c Re-implementing the "limited free slip" boundary condition at R21 and R42
c resolutions. The new implementation is simpler, and corrects the bugs that 
c were present in the previous version.
c SJP 2009/08/04
c
c Tidying up the source code for the sea ice model, particularly: fixing the
c bug in the calculation of the latitudes; and removing the buggy and
c experimental code which sought to impose "limited slip at low resolution".
c SJP 2009/04/25
c
c $Log: COMDICE.f,v $
c Revision 1.21  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.20  1998/12/10  00:55:44  ldr
c HBG changes to V5-1-21
c
c Revision 1.19  1997/12/19  02:03:16  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.18  1996/10/24  01:02:26  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.17  1996/03/21  03:18:27  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.16  1995/08/31  04:30:39  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.15  1995/07/26  07:28:38  ldr
c Merge Mickles speedups to ice scheme (V4-5-30mic) into V4-7-3l.
c
c Revision 1.14  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.13  1994/08/08  17:20:44  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.12  94/08/04  16:53:47  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.11.1.1  1995/07/26  07:24:08  ldr
c Speedups from Mickles replacing character mask with integers.
c
c Revision 1.11  94/03/30  12:33:30  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.10  93/12/17  15:31:19  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.9  93/11/29  14:40:20  ldr
c Merge of HBG and SPO changes.
c 
c Revision 1.8  93/11/29  11:38:18  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.7.1.1  93/11/29  12:32:46  ldr
c SPO's changes to rationalize printing of icuv statistics.
c 
c Revision 1.7  93/03/16  16:46:09  ldr
c Removed redundant tstx from common block dicemask.
c 
c Revision 1.6  93/03/15  17:20:59  ldr
c HBG changes to reduce size of /comdice common block.
c 
c Revision 1.5  93/02/03  12:43:19  ldr
c SPO's minor changes to coupled model and sea-ice scheme.
c 
c Revision 1.4  92/12/09  14:42:05  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c
c     Parameters and commons for the dynamical seaice model.
c
c        
c        PLON     = lsx longitude dimension
c        PLAT     = lsx latitude dimension
c        PICELAY  = number of seaice layers in seaice model (= PSOILAY)
c        PSNOLAY  = number of snow layers in snow model (ge 2)
c
c
      INTEGER  PLON,  PLAT,    PICELAY, PSNOLAY 
c
      PARAMETER (
     &            PLON  = lon ,         PLAT  = lat*2,             
     &            PSNOLAY = 1,         picelay = 1)
      real dtf
      common/icetime/dtf

c
c radius of earth (m)
      real tmelt ! freezing point of water (except seaice bottom)
      real rhoi  ! density of ice (all types)
      parameter (rhoi=0.9e3,tmelt=273.05)

c rotation rate of earth (rad s-1)
      real omega
C      parameter (omega=(2.*pi/86400.))
      parameter (omega=7.2921233e-5) ! Standard value for omega chosen by
c                                      HBG and used in PHYSPARAMS.f

c turning angle for water stress (deg)
      real wturn,aturn
c     parameter (wturn=25.,aturn=20.)
      parameter (wturn=25.,aturn=0.)

c drag coefficients for water and air ((Nl m-2)/(m s-1))
      real csubw, csuba
c     parameter (csubw=1.5*0.6524, csuba=0.01256)
      parameter (csubw=3.0*0.6524, csuba=0.01256)

c length scale of ice pressure "transmittal", used only in cavit2 (m)
      real streslen
      parameter (streslen=1.e5)
      integer ioterm
      parameter(ioterm=6)

c flag for which cavitating iteration method to use.
c if true,  solve for free u,v just once and increment u,v (cavit)
c if false, solve for "free+pressure" u,v and increment pressure(cavit2)
      logical flato
      parameter (flato=.true.)

c max, min numbers of cavit iterations,and convergence criterion (m s-1)
      integer nloop, nloopmin
      real convinc
c     parameter (nloop=100, nloopmin=15, convinc=.0001  )
      parameter (nloop=150, nloopmin=15, convinc=.0001  )

c        EPSILON      = small quantity to avoid zero-divides and other
c                       truncation or machine-limit troubles with small
c                       values. Should be slightly greater than O(1) 
c                       machine precision.
      real epsilon
      parameter (epsilon=1.e-12)
c grids 
      real alon(0:plon+1)
      real alonu(0:plon+1)
      real slonu(0:plon+1)
      real clonu(0:plon+1)
      real blat(plat)
      real blon(plon)  
      real alat(0:plat+1)
      real alatu(plat+2)
      real blatu(plat+1)
      real blonu(0:plon+1)
      real slat(0:plat+1)
      real slatu(plat+1)
      real deltx(0:plat+1)
      real deltxu(plat+1)
      real delty(0:plat+1)
      real deltyu(plat+1)
      real srfarea(plat)
      real srfareau(plat+1)
      integer latha, lathb, latua, latub
      common /dicegrid/ alon, alonu, slonu, clonu, blat, blon, alat,
     &                  alatu, blatu, blonu, slat, slatu, deltx,
     &                  deltxu, delty, deltyu, srfarea, srfareau,
     &                  latha, lathb, latua, latub

c water-stess turning angles 
      real coswturn(plat+1), sinwturn(plat+1)
      real cosaturn(plat+1), sinaturn(plat+1)
      common /diceturn/ coswturn, sinwturn, cosaturn, sinaturn

c u-grid land-ocean mask umask, and h-grid "divergence" mask
      real umask(0:plon+1,plat+1)
      real dmask(plon,plat)
      real adgain(plon,plat)
      integer imask(plon,plat), kma(plon,plat)
      logical cmask(0:plon+1,plat+1)
      common /dicemask/ umask, dmask, adgain, imask, kma, cmask

      real u(0:plon+1,plat+1)
      real v(0:plon+1,plat+1)
      real ocu(plon+2,plat+2,2)
      real ocv(plon+2,plat+2,2)
      real uav(plon,plat)
      real vav(plon,plat)
      real fav(plon,plat)
      real uo(0:plon+1,plat+1)
      real vo(0:plon+1,plat+1)    
      real wl(plon,plat)
      real wd(plon,plat)

      common/ ocwind/ u, v, ocu, ocv, uav, vav, fav, uo, vo    
     & ,wl ,wd

c diagnostic output 
      integer loop(2), icolaps(plon,plat)
      real rmsinc(2), rmsvel(2)
      common /dicediag/ loop, rmsinc, rmsvel, icolaps

c work space
      real work(plon,plat)
      real work2(plon,plat)
      real workh(0:plon+1,0:plat+1)
      real worku(0:plon+1,plat+1)
      real worku2(0:plon+1,plat+1)
      common /dicework/ work, work2, workh, worku, worku2

      real workp(plon,plat),work2nd(plon,plat)
      real workl(plon,plat),workdf(plon,plat)
      common /divwork/ workp,work2nd,workl,workdf
