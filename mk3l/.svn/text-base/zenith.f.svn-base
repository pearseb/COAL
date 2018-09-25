c $Log: zenith.f,v $
c Revision 1.7  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.6  1998/12/10  00:55:54  ldr
c HBG changes to V5-1-21
c
c Revision 1.5  1996/10/24  01:03:22  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1994/03/23  10:20:28  mrd
c Input time now ordinary time, converted to be zero at noon inside zenith.
c
c Revision 1.3  93/12/17  15:34:30  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.2  93/08/03  11:55:57  ldr
c  Replace intrinsic functions with generic forms where possible.
c 
c Revision 1.1  91/02/22  16:38:24  ldr
c Initial release V3-0
c 
      subroutine zenith(fjd,r,dlt,slag,xlat,ha,dhr,nlng,coszrs,frac)
C
C     ZENITH COMPUTES EFFECTIVE MEAN COSINE OF ZENITH ANGLE AND DAYLIGHT
C       FRACTION FROM LATITUDE AND PARAMETERS GIVEN BY SUBROUTINE SOLAR
C
C     INPUT ARGUMENTS TO CIRCULAR FUNCTIONS ARE IN RADIANS
C
C Global parameters
      include 'PHYSPARAMS.f'

C Argument list
      real fjd
      real r
      real dlt
      real slag
      real xlat
      real ha
      real dhr
      integer nlng
      real coszrs(nlng)
      real frac(nlng)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      logical rise,set

C Local data, functions etc
      data tpi/6.2831853/,pid24/.13089969/

C Start code : ----------------------------------------------------------

      rs=1.
      cvpr=tpi/float(nlng)
      gha=fjd*tpi+slag+pi    ! Add pi because time is zero at local noon
      arg=dhr*pid24
      sinfac=sin(arg)/arg
      ss=sin(xlat)*sin(dlt)
      cc=cos(xlat)*cos(dlt)
      if(ha.gt.0.0)cons=rs*(ss+cc*sin(ha)/ha)
C#DO MODE(SCALAR)
      do 20 i = 1,nlng
      xlng=cvpr*(float(i)-1)
C     HLOC=GHA+XLNG
      hloc=gha+xlng+arg
C     LOCAL HOUR ANGLE SHIFTED BY HALF OF THE AVERAGING PERIOD
      rise=.false.
      set=.false.
      hloc=mod(hloc,tpi)
      if(hloc.gt.pi) hloc=hloc-tpi
      hlpar=hloc+arg
      armhl=arg-hloc
      if(hlpar.gt.ha) set=.true.
      if(armhl.gt.ha) rise=.true.
      if(rise.and.set) go to 7
      if(hlpar.gt.pi) go to 8
      if(armhl.gt.pi) go to 9
      if(set) go to 3
      if(rise) go to 4
      frac(i)=1.0
      coszrs(i)=rs*(ss+cc*cos(hloc)*sinfac)
      go to 2
  3   delsh=.5*(ha+armhl)
      go to 5
  4   delsh=.5*(ha+hlpar)
  5   if(delsh.le.0.0) go to 6
      frac(i)=delsh/arg
      coszrs(i)=rs*(ss+cc*cos(ha-delsh)*sin(delsh)/delsh)
      go to 2
  7   if(ha.le.0.) go to 6
      frac(i)=ha/arg
      coszrs(i)=cons
      go to 2
  8   dele=.5*max(hlpar+ha-tpi,0.)
      delw=.5*max(ha+armhl,0.)
      go to 10
  9   dele=.5*max(ha+hlpar,0.)
      delw=.5*max(armhl+ha-tpi,0.)
 10   frac(i)=(dele+delw)/arg
      if(frac(i).eq.0.) go to 11
      coszrs(i)=rs*(ss+cc*(cos(ha-dele)*sin(dele)+
     &                     cos(ha-delw)*sin(delw))/(dele+delw))
      go to 2
    6 frac(i)=0.0
   11 coszrs(i)=0.0
  2   continue
      coszrs(i) = min(1.0,coszrs(i))
      frac(i) = min(1.0,frac(i))
   20 continue
      return
C END OF ZENITH
      end
