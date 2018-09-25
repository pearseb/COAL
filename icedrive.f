c Correcting an apparent energy conservation error within the final loop. At
c the two other locations within the model where mixed-layer ocean points are
c changed to sea points, the energy change is accounted for. The final loop is
c therefore modified so that the energy change is also accounted for here.
c SJP 2009/05/11
c
c Tidying up the source code for the sea ice model, particularly: fixing the
c bug in the calculation of the latitudes; and removing the buggy and 
c experimental code which sought to impose "limited slip at low resolution".
c SJP 2009/04/25
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Transferred COMMON blocks to separate header files, as follows:
c /MDAY/  ->  MDAY.f
c SJP 2007/05/29
c
c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c $Log: icedrive.f,v $
c Revision 1.41  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.40  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.39  2000/06/20 02:08:34  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.38  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.37  1998/12/10  00:55:50  ldr
c HBG changes to V5-1-21
c
c Revision 1.36  1997/12/23  00:23:37  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.35  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.34.1.1  1997/12/19  02:03:17  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.34  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.33  1996/10/24  01:02:55  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.32  1996/06/13  02:06:50  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.31  1996/03/21  03:18:50  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.30  1995/02/27  01:40:34  ldr
c Bugfixes received from HBG in response to problems on Fujitsu.
c
c Revision 1.29  1994/08/08  17:21:33  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.28  94/08/04  16:55:38  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.27  94/03/30  16:13:12  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.26  93/12/17  15:32:50  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.25  93/12/17  11:51:38  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.24  93/11/29  16:09:17  ldr
c Corrected length of common block hm.
c 
c Revision 1.23  93/11/29  11:38:30  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.22  93/11/03  11:44:16  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.21  93/10/05  13:06:24  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.20  93/08/10  15:27:22  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c
c     INPUT/OUTPUT
c     Input:   from common/aticstr in this subroutine
c                  ttaux,ttauy - ice stress accumulator arrays
c
c              from common/dicemask in COMDICE.f
c                  adgain - advection of ice
c
c              from common/fewflags in FEWFLAGS.f
c                  lcouple - flag to run model in coupled mode, default F
c                  qflux   - set to T for qflux run, default F
c
c              from common/freeze in FREEZE.f
c                  tfi - temperature at bottom of ice
c
c              from common/gausl in GAUSL.f
c                  w - normal guass weights
c
c              from common/ocwind in COMDICE.f
c                  ocu - ocean currents, zonal
c                  ocv - ocean currents, meridional
c
c              from common/timex in TIMEX.f
c                  dt - time step in seconds
c                  month - month counter
c                  ratlm - fraction of month elapsed, runs 0-1
c
c              from arguments
c                  dtice - time step for ice scheme = 1 hr 
c                  id - day identifier (0 during first day etc)
c
c     Output:  from common/dicegrid in COMDICE.f
c                  latua,latub - latitudinal indices bounding the extent
c                                of ice (u-grid)
c
c              from common/icetime in COMDICE.f
c                  dtf - dynamic seaice time step, seconds
c
c     In/Out:  from common/aforce in this subroutine
c                  athf - heat flux
c                  athfa - heat flux with zero sub-ice heat input
c                  dsfw - vol. fresh water on ice
c                  dsis - vol. ice thickness*concentration
c
c              from common/curent in this subroutine
c                  occur - ocean currents
c
c              from common/ficeug in this subroutine
c                  ficeu - ice concentration on u-grid 
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from common/masiv4 in MASIV4.f
c                  savegrid - global surface variable
c
c              from common/masiv5 in this subroutine
c                  statsice - global statistics relating to sea-ice scheme
c
c              from common/ocwind in COMDICE.f
c                  fav - average ice concentration
c                  hav - advective ice
c                  u - eastward  velocity on u-grid
c                  uav - average zonal velocities
c                  v - northward velocity on u-grid 
c                  vav - average meridional velocities
c                  wl cumulative array for residual divergence
c                  wd cumulative array for effect of rheology on divergence
c
c              from common/red
c                  hav - advective ice
c 
      subroutine icedrive(id,dtice,ijdyn)

c Dummy driver for dynamic ice model. All arrays dimensioned here will
c be parts of the agcm or ogcm. They are described in dynice.

      implicit none
C Global parameters
      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

C Argument list
      integer id
      real dtice
      integer ijdyn

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'FREEZE.f' !
      include 'GAUSL.f'
      include 'LSMI.f'
      include 'MASIV4.f'
      include 'TIMEX.f'
      include 'MDAY.f'
      
      include 'gas.h'

      real drunof,dsis,dsisb,dsfw,athf,athfa,atsi
      common/aforce/drunof(ln2,lat),dsis(ln2,lat),dsisb(ln2,lat)
     & ,dsfw(ln2,lat),athf(ln2,lat),athfa(ln2,lat),atsi(ln2,lat)

      real ttaux,ttauy
      common/aticstr/ttaux(ln2,lat),ttauy(ln2,lat)

      real ochf,occur,g50dt
      integer khf
      common/curent/ochf(ln2,lat),occur(ln2,lat),khf(ln2,lat)
     &,g50dt(ln2,lat)

      real ficeu,fice
      common/ficeug/ficeu(0:plon+1,plat+1),fice(plon,plat)

      integer idfrz
      common/icedvz/idfrz(lat,2)

      real hmo,hmo2
      common/hm/hmo(ln2,lat),hmo2(lon,lat2,2)

      real statsice
      common/masiv5/statsice(ln2,14,lat)

      real redice,groice,hav
      common/red/redice(ln2,lat),groice(ln2,lat),hav(ln2,lat)

C Local work arrays and variables
      real                           stor(plon,plat),
     &   hice(plon,plat),            ticel(plon,plat,2),
     &   hsno(plon,plat),            tsno(plon,plat),
     &   fuwindi(plon,plat),         fvwindi(plon,plat),
     &   uocn(0:plon+1,plat+1),      vocn(0:plon+1,plat+1),
     &   ficex(plon,plat),hicex(plon,plat)
      real hmox(plon)
      real hcap1(plon)
      logical mlox(plon),seax(plon)

      integer i
      integer imlgp
      integer j
      integer lg
      integer ma
      integer mg
      integer ns

      real b1
      real b2
      real daysm
      real delt
      real div
      real emj
      real hsn
      real qheat
      real qtot
      real tintp
      real tmlo
      real ts1
      real tst1
      real tst2

C Local data, functions etc
      real shi,shs,qice,hcap
      data shi/1.88837e6/,shs/6.9069e5/,qice/3.36e8/,hcap/2.095e8/

C Start code : ----------------------------------------------------------

      latua=1
      latub=plat+1

      dtf=dtice
      div=float(ijdyn)
c Ice dynamics time loop
      
c.... Get the ocean model velocities :

      if(.not.lcouple)then

c.... ocu,ocv hold ocean model time averaged velocities : (plon+2,plat+2)
c.... These are set for i=1,plon+2 and j=2,plat+1 only with 
c.... j=plat+1 (NP row) being zeroes. (see flatset.f)
c---- Data comming in is at First of each month (not middle)
c---- for ocu,ocv
      daysm=float(mdays(month))
c---- id=ldays is day identifier (0 during first day etc)
      b2=id/daysm

      else

c.... ocu,ocv hold ocean model instantaneous velocities : (plon+2,plat+2)
c.... (see icecurr.f).
      b2=0.0

      endif

      b1=1.0-b2
      do 868 j=2,plat+1
      do 868 i=1,plon+1
      uocn(i,j)=b1*ocu(i,j,1)+b2*ocu(i,j,2)
  868 vocn(i,j)=b1*ocv(i,j,1)+b2*ocv(i,j,2)

c.... wrap around for i=0 element
      do 867 j=2,plat+1
      uocn(0,j)=uocn(plon,j)
  867 vocn(0,j)=vocn(plon,j)

      call polefilt (uocn, vocn, 86, 3)

c.... Extend the atmospheric ice stresses over adjacent non-ice points
c.... Done before the ice model averages these stresses to the ice-model
c.... U/V grid (offset from atmos model).
c.... Note that the U/V offset is different to the ocean model U/V offset
c....
c....      H(i-1,j  )   H(i  ,j  )
c....
c....            U(i  ,j  )          : This is a SW offset
c.... 
c....      H(i-1,j-1)   H(i  ,j-1)
c....
c.... (The ocean model has a NE offset).
      call icetau

      do 845 j=1,plat
      ns=2-(j-1)/lat
      lg=(ns-1)*j+(lat2+1-j)*(2-ns)
      do 845 i=1,plon
      ma=i+(ns-1)*plon
      fuwindi(i,j)=ttaux(ma,lg)/div
      fvwindi(i,j)=ttauy(ma,lg)/div
      hsno(i,j)=statsice(ma,1,lg)
      hice(i,j)=statsice(ma,2,lg)
      tsno(i,j)=statsice(ma,3,lg)
      ticel(i,j,1)=statsice(ma,4,lg)
      ticel(i,j,2)=statsice(ma,5,lg)
      stor(i,j)=statsice(ma,7,lg)
      if(imsl(ma,lg).eq.1)then 
      fice(i,j)=1-statsice(ma,8,lg)
      else
      fice(i,j)=0.0
      endif
 845  continue

c.... Extend the ice values (fice,hice) over adjacent land points.
c.... This creates arrays ficex,hicex.
c.... Done before the ice model uses ltou/ltoh to interpolate these
      call icefhx(fice,ficex,hice,hicex)

c       Step ice and snow through one dynamic time step

        call dynice(fice, hice, ticel, hsno, tsno, fuwindi, fvwindi,
     &              uocn, vocn, ficeu, ficex, hicex, stor)

       do 560 j=1,plat
       do 560 i=1,plon
       if(ficeu(i,j).le.0.005)then
         u(i,j)=0.0
         v(i,j)=0.0
       endif
 560   continue

       do 433 j=1,plat
       do 433 i=1,plon
       fav(i,j)=fav(i,j)+fice(i,j)
       wl(i,j)=wl(i,j)+workl(i,j)
       wd(i,j)=wd(i,j)+workdf(i,j)
       uav(i,j)=uav(i,j)+u(i,j)
 433   vav(i,j)=vav(i,j)+v(i,j)
c rjm
	do 443 j=1,plat/2
	do 443 i=1,plon
	 sice(i+plon,j) = fice(i,j)
 443	 sice(i,j) = fice(i,plat+1-j)
c rjm


c.... Set indicator for "Equatorward ice growth" to zero
c.... (Used in radin.f in same way as ifroz)
      do ns=1,2
      do lg=1,lat
        idfrz(lg,ns)=0
      enddo
      enddo

c.... Loop over all latitudes checking for ice drift onto
c.... MLO and SEA points. May need to reset the IMSL mask.

      do 874 j=1,plat
      ns=2-(j-1)/lat
      lg=(ns-1)*j+(lat2+1-j)*(2-ns)

      if(qflux)then !Use IGW's qflux code with variable MLO depth
       do 52 i=1,plon
       ma=i+(ns-1)*plon
  52   hmox(i)=hmo(ma,lg)
      end if

c.... Set up masks for ice over MLO or SEA points
      do 56 i=1,plon
      ma=i+(ns-1)*plon
      hav(ma,lg)=hav(ma,lg)+adgain(i,j)
      mlox(i)=.false.
      seax(i)=.false.
      if(fice(i,j).gt.0.0)then
        if(iabs(imsl(ma,lg)).eq.2)mlox(i)=.true.
        if(imsl(ma,lg).eq.3)seax(i)=.true.
      end if
   56 continue

c Use IGW's qflux code with variable MLO depth

      if(qflux)then
        do i=1,plon
          hcap1(i)=hcap/50.*hmox(i)
        enddo
      else
        do i=1,plon
          hcap1(i)=hcap*2.
        enddo
      endif

c.... Remove small ice amounts over M or S points
      do 843 i=1,plon
      ma=i+(ns-1)*plon

      if((mlox(i).or.seax(i)).and.(fice(i,j).le.0.005))then

      qheat=(hice(i,j)*0.5*shi*((tfi-ticel(i,j,1))+(tfi-ticel(i,j,2))))
     &+((tfi-tsno(i,j))*hsno(i,j)*shs)
      hsn=hice(i,j)+0.1*hsno(i,j)
      qtot=(qheat+qice*hsn)*fice(i,j)
      delt=qtot/((1-fice(i,j))*hcap1(i))
      if(mlox(i))savegrid(ma,3,lg)=savegrid(ma,3,lg)-delt
      athf(ma,lg)=athf(ma,lg)-qtot/dt
      athfa(ma,lg)=athfa(ma,lg)-qtot/dt
      dsfw(ma,lg)=dsfw(ma,lg)+0.1*hsno(i,j)*fice(i,j)
      dsis(ma,lg)=dsis(ma,lg)+hice(i,j)*fice(i,j)
      fice(i,j)=0.0
      hice(i,j)=0.0
      ticel(i,j,1)=tmelt
      ticel(i,j,2)=tmelt
      tsno(i,j)=tmelt
      hsno(i,j)=0.0  

      endif
 843  continue

c.... Convert to Ice points if enough ice over M or S points
c.... and the temp is below critical, else melt the ice.
      do 844 i=1,plon
      ma=i+(ns-1)*plon

      if((mlox(i).or.seax(i)).and.(fice(i,j).gt.0.005))then

      qheat=(hice(i,j)*0.5*shi*((tfi-ticel(i,j,1))+(tfi-ticel(i,j,2))))
     &+((tfi-tsno(i,j))*hsno(i,j)*shs)
      hsn=hice(i,j)+0.1*hsno(i,j)
      qtot=(qheat+qice*hsn)*fice(i,j)
      ts1=savegrid(ma,3,lg)

      if(ts1.ge.273.15)then
      delt=qtot/((1-fice(i,j))*hcap1(i))
      if(mlox(i))savegrid(ma,3,lg)=savegrid(ma,3,lg)-delt
      athf(ma,lg)=athf(ma,lg)-qtot/dt
      athfa(ma,lg)=athfa(ma,lg)-qtot/dt
      dsfw(ma,lg)=dsfw(ma,lg)+0.1*hsno(i,j)*fice(i,j)
      dsis(ma,lg)=dsis(ma,lg)+hice(i,j)*fice(i,j)
      fice(i,j)=0.0
      hice(i,j)=0.0
      ticel(i,j,1)=tmelt
      ticel(i,j,2)=tmelt
      tsno(i,j)=tmelt
      hsno(i,j)=0.0  

      else

      imsl(ma,lg)=1
      if(imsl(ma,lg+1).eq.3)imsl(ma,lg+1)=-2
      statsice(ma,9,lg)=savegrid(ma,3,lg)
      idfrz(lg,ns)=1

      endif ! 273.15

      endif 
 844  continue

c.... If Ice point and too little ice, then remove + change point to M
      do 842 i=1,plon
      ma=i+(ns-1)*plon

      if(imsl(ma,lg).eq.1.and.fice(i,j).le.0.005)then

      fice(i,j)=0.0
      hice(i,j)=0.0
      ticel(i,j,1)=tmelt
      ticel(i,j,2)=tmelt
      tsno(i,j)=tmelt
      hsno(i,j)=0.0  
      imsl(ma,lg)=-2
      if((imsl(ma,lg+1).eq.2).and.(.not.qflux)) then
            imsl(ma,lg+1)=3
c.... Equatorward MLO point changed to a SEA point. Temperature jump.
c.... Account for energy.
            tmlo=savegrid(ma,3,lg+1)
            tst1=savegrid(ma,5,lg+1)
            tst2=savegrid(ma,6,lg+1)
            tintp=tst1*(ratlm-1.0)-tst2*ratlm
            emj=hcap1(i)*(tintp-tmlo)/dt
c.... Add as term in calculation of Q-flux
            occur(ma,lg+1)=occur(ma,lg+1)+emj
      endif

      endif
 842  continue     

c.... Recheck all points for sub-min ice.
      do  841 i=1,plon

      if(fice(i,j).le.0.005)then

      fice(i,j)=0.0
      hice(i,j)=0.0
      ticel(i,j,1)=tmelt
      ticel(i,j,2)=tmelt
      tsno(i,j)=tmelt
      hsno(i,j)=0.0  

      endif
 841  continue     

      do 848 i=1,plon
      ma=i+(ns-1)*plon

      if(imsl(ma,lg).eq.1)then
      statsice(ma,8,lg)=1-fice(i,j)
      else
      statsice(ma,8,lg)=0.0
      endif

      statsice(ma,1,lg)=hsno(i,j)
      statsice(ma,2,lg)=hice(i,j)
      statsice(ma,3,lg)=tsno(i,j)
      statsice(ma,4,lg)=ticel(i,j,1)
      statsice(ma,5,lg)=ticel(i,j,2)
      statsice(ma,7,lg)=stor(i,j)

      if(imsl(ma,lg).lt.4)then
      savegrid(ma,12,lg)=hsno(i,j)*100.0
      savegrid(ma,13,lg)=hice(i,j)  
      endif

 848  continue

 874  continue ! j=1,plat

c.... Check the whole IMSL mask at this stage. Remove any isolated
c.... MLO points or adjacent MLO points between latitude rows
c.... These may arise due to E-W ice drift.

c.... Modify the excess M points to S points (2 to 3)
c.... (the mask is 1=Ice, 2 or -2 is MLO, 3 is Sea, 4 is Land)
c.... This must not be done if a Q-flux run (all M points!) and is
c.... not needed for a coupled run (all reset to S points in SURFSET).

      if(lcouple.or.qflux)return

      do 41 mg=1,ln2
      do 41 lg=1,lat-1
      if((imsl(mg,lg).eq.1).or.(imsl(mg,lg).eq.4))go to 41
c..   either M or S at (mg,lg)
c..   check for M at (mg,lg+1) and reset to S
        imlgp=iabs(imsl(mg,lg+1))
        if (imlgp .eq. 2) then
          imsl(mg, lg+1) = 3
c.... Equatorward MLO point changed to a SEA point. Temperature jump.
c.... Account for energy.
          tmlo=savegrid(mg,3,lg+1)
          tst1=savegrid(mg,5,lg+1)
          tst2=savegrid(mg,6,lg+1)
          tintp=tst1*(ratlm-1.0)-tst2*ratlm
          emj=hcap*2.0*(tintp-tmlo)/dt
c.... Add as term in calculation of Q-flux
          occur(mg,lg+1)=occur(mg,lg+1)+emj
        endif
   41 continue

      return
      end
