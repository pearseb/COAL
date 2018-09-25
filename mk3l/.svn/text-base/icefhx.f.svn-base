c $Log: icefhx.f,v $
c Revision 1.7  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.6  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.5  1996/10/24  01:03:28  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1994/08/08  16:51:12  ldr
c Fix up RCS header at top of file.
c
c     INPUT/OUTPUT
c     Input:   from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from arguments
c                  fice - ice concentration on h-grid  
c                  hice - ice thickness on h-grid
c
c     In/Out:  from arguments
c                  ficex, hicex - work arrays
c
      subroutine icefhx(fice,ficex,hice,hicex)

c Filler routine for ice values (fractional cover and depth).
c Set all land values to zero.
c Then average points with values to get expansion
c over the adjacent land points.
c This is done so that the change to the ice u-grid
c which uses area weighting
c gets sufficient surrounding vales with realistic values.

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real fice(lon,lat2)
      real ficex(lon,lat2)
      real hice(lon,lat2)
      real hicex(lon,lat2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'LSMI.f'

C Local work arrays and variables
      real fx(0:lon+1,0:lat2+1),hx(0:lon+1,0:lat2+1)

      integer i
      integer isum
      integer j
      integer lg
      integer ma
      integer mg
      integer ns

C Local data, functions etc
      logical first
      data first/.true./
      integer mskl(0:lon+1,0:lat2+1),msklsum(lon,lat2)
      save first,mskl,msklsum

C Start code : ----------------------------------------------------------

      if(first)then

c.... Create a mask of 0 for land, else 1
      do 11 j=1,lat2
      ns=2-(j-1)/lat
      lg=(ns-1)*j+(lat2+1-j)*(2-ns)
      do 11 mg=1,lon
      i=mg
      ma=mg+(ns-1)*lon
      mskl(i,j)=1
      if(imsl(ma,lg).eq.4)mskl(i,j)=0
   11 continue

      do 21 j=1,lat2
      mskl(0,j)=mskl(lon,j)
   21 mskl(lon+1,j)=mskl(1,j)

      do 41 i=0,lon+1
      mskl(i,0)=mskl(i,1)
   41 mskl(i,lat2+1)=mskl(i,lat2)

      do 45 j=1,lat2
      do 45 i=1,lon
      msklsum(i,j)=0
      if(mskl(i,j).eq.0)then
      msklsum(i,j)=
     &  mskl(i+1,j  )+mskl(i+1,j+1)+mskl(i  ,j+1)+mskl(i-1,j+1)
     & +mskl(i-1,j  )+mskl(i-1,j-1)+mskl(i  ,j-1)+mskl(i+1,j-1)
      endif
   45 continue

        first=.false.
      end if

c.... Set all land point values to zero
      do 10 j=1,lat2
      do 10 i=1,lon
      ficex(i,j)=fice(i,j)*mskl(i,j)
   10 hicex(i,j)=hice(i,j)*mskl(i,j)

c.... Make temporary global (i,j) arrays with
c.... wrap around and extend N & S.
      do 20 j=1,lat2
      do 30 i=1,lon
      fx(i,j)=ficex(i,j)
   30 hx(i,j)=hicex(i,j)
      fx(0,j)=fx(lon,j)
      hx(0,j)=hx(lon,j)
      fx(lon+1,j)=fx(1,j)
   20 hx(lon+1,j)=hx(1,j)

      do 40 i=0,lon+1
      fx(i,0)=fx(i,1)
      hx(i,0)=hx(i,1)
      fx(i,lat2+1)=fx(i,lat2)
   40 hx(i,lat2+1)=hx(i,lat2)

c.... Run through the mask checking all points that have a
c.... zero in the mask (i.e. a land point).
c.... If zero, try to create a value at that point by averaging
c.... the 8 surrounding points, using only values that are 
c.... non-zero (non-land values)
      do 50 j=1,lat2
      do 50 i=1,lon
      isum=msklsum(i,j)
      if(isum.gt.0)then
      ficex(i,j)=(
     &  fx(i+1,j  )+fx(i+1,j+1)+fx(i  ,j+1)+fx(i-1,j+1)
     & +fx(i-1,j  )+fx(i-1,j-1)+fx(i  ,j-1)+fx(i+1,j-1)
     & )/isum
      hicex(i,j)=(
     &  hx(i+1,j  )+hx(i+1,j+1)+hx(i  ,j+1)+hx(i-1,j+1)
     & +hx(i-1,j  )+hx(i-1,j-1)+hx(i  ,j-1)+hx(i+1,j-1)
     & )/isum
      endif
   50 continue
     
      return
      end
