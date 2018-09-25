c $Log: icetau.f,v $
c Revision 1.7  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.6  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.5  1996/10/24  01:03:28  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1994/08/08  16:52:50  ldr
c Fix up RCS header at top of file.
c
c     INPUT/OUTPUT
c     Input:   from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c     In/Out:  from common/aticstr in this subroutine
c                  ttaux, ttauy - ice stress accumulator arrays
c
      subroutine icetau

c Filler routine for ice stresses from atmos model
c Set all non-ice stresses to zero.
c Then average points with stresses to get expansion
c of stresses over the adjacent non-ice points.
c This is done so that the change to the ice U/V-grid
c (which is different to the atmos model U/V-grid)
c done by timefilt,polefilt, gets sufficient surrounding
c values with realistic ice magnitude stresses.

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'LSMI.f'

c   atmos-ice stresses
      real ttaux,ttauy
      common/aticstr/ttaux(ln2,lat),ttauy(ln2,lat)

C Local work arrays and variables
      real staux(0:lon+1,0:lat2+1),stauy(0:lon+1,0:lat2+1)
      integer mski(0:lon+1,0:lat2+1),msksum(lon,lat2)

      integer i
      integer isum
      integer j
      integer lg
      integer ma
      integer mg
      integer ns

C Local data, functions etc

C Start code : ----------------------------------------------------------

c.... Set all non-ice point stresses to zero
      do 10 lg=1,lat
      do 10 mg=1,ln2
      if(imsl(mg,lg).ne.1)then
        ttaux(mg,lg)=0.0
        ttauy(mg,lg)=0.0
      endif
   10 continue

c.... Make global (i,j) arrays from (mg,lg) arrays with
c.... wrap around and extend N & S.
c.... Create a mask of 1 for ice, else 0
      do 20 j=1,lat2
      ns=2-(j-1)/lat
      lg=(ns-1)*j+(lat2+1-j)*(2-ns)
      do 30 mg=1,lon
      ma=mg+(ns-1)*lon
      i=mg
      mski(i,j)=0
      if(imsl(ma,lg).eq.1)mski(i,j)=1
      staux(i,j)=ttaux(ma,lg)
   30 stauy(i,j)=ttauy(ma,lg)
      mski(0,j)=mski(lon,j)
      staux(0,j)=staux(lon,j)
      stauy(0,j)=stauy(lon,j)
      mski(lon+1,j)=mski(1,j)
      staux(lon+1,j)=staux(1,j)
   20 stauy(lon+1,j)=stauy(1,j)
      do 40 i=0,lon+1
      mski(i,0)=mski(i,1)
      staux(i,0)=staux(i,1)
      stauy(i,0)=stauy(i,1)
      mski(i,lat2+1)=mski(i,lat2)
      staux(i,lat2+1)=staux(i,lat2)
   40 stauy(i,lat2+1)=stauy(i,lat2)

      do 45 j=1,lat2
      do 45 i=1,lon
      msksum(i,j)=0
      if(mski(i,j).eq.0)then
      msksum(i,j)=
     &  mski(i+1,j  )+mski(i+1,j+1)+mski(i  ,j+1)+mski(i-1,j+1)
     & +mski(i-1,j  )+mski(i-1,j-1)+mski(i  ,j-1)+mski(i+1,j-1)
      endif
   45 continue

c.... Run through the mask checking all points that have a
c.... zero in the mask (i.e. a non ice point).
c.... If zero, try to create a stress at that point by averaging
c.... the 8 surrounding points, using only values that are 
c.... non-zero (ice values)
      do 50 j=1,lat2
      ns=2-(j-1)/lat
      lg=(ns-1)*j+(lat2+1-j)*(2-ns)
      do 60 i=1,lon
      isum=msksum(i,j)
      if(isum.gt.0)then
      ma=i+(ns-1)*lon
      ttaux(ma,lg)  =( staux(i+1,j  )+staux(i+1,j+1)
     & +staux(i  ,j+1)+staux(i-1,j+1)+staux(i-1,j  )
     & +staux(i-1,j-1)+staux(i  ,j-1)+staux(i+1,j-1) )/isum
      ttauy(ma,lg)  =( stauy(i+1,j  )+stauy(i+1,j+1)
     & +stauy(i  ,j+1)+stauy(i-1,j+1)+stauy(i-1,j  )
     & +stauy(i-1,j-1)+stauy(i  ,j-1)+stauy(i+1,j-1) )/isum
      endif
   60 continue
   50 continue
     
      return
      end
