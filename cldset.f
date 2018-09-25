c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Added IMPLICIT NONE statement, plus variable declarations for PI, TPI, RLAG,
c YEAR, DATE, RANG, RSIN1, RCOS1, RCOS2, THETA, L, K, CCD, CCD2, CCD3, CCD4,
c KKTH and KKBH.
c SJP 2007/05/28
c
c $Log: cldset.f,v $
c Revision 1.7  1998/12/10 00:55:46  ldr
c HBG changes to V5-1-21
c
c Revision 1.6  1996/10/24  01:02:31  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.5  1996/06/13  02:05:41  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.4  1993/08/03  11:55:46  ldr
c  Replace intrinsic functions with generic forms where possible.
c
c Revision 1.3  92/11/30  10:54:41  mrd
c Changed date calculation to use the model year of exactly 365,
c not 365.25 days.
c 
c Revision 1.2  92/04/16  16:38:44  ldr
c Moved initialization stuff to initfs.
c 
c Revision 1.1  91/02/22  16:37:01  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT
c     Input:   from common/clddat in this subroutine
c                  ccd, ccd2, ccd3, ccd4, kkbh, kkth
c
c              from arguments
c                  alat - latitude (from -pi/2 to pi/2)
c
c     Output:  from arguments
c                  ch - high cloud                cm - mid cloud
c                  icb - bottom of low cloud     ich - level of high cloud 
c                  icm - level of medium cloud   ict - top of low cloud
c
c     In/Out:  from arguments
c                  mins - current model time in mins
c
 
      subroutine cldset(alat,mins,ch,cm,cl,ich,icm,ict,icb)
c
c  This routine interpolates in latitude and time to set the cloud
c  amounts.
c  INPUT
c    ALAT    latitude (from -pi/2 to pi/2)
c    MINS    current model time in mins
c  OUTPUT
c    CD      cloud fraction
c    KTH     Level of cloud top
c    KBH     Level of cloud base
C
      implicit none

C Global parameters
      real pi, tpi, rlag, year
      parameter(pi=3.141592653589793,tpi=2.*pi,rlag=14.8125)
      parameter(year=365)

C Argument list
      real alat
      integer*8 mins
      real ch
      real cm
      real cl
      integer ich
      integer icm
      integer ict
      integer icb

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
c          Winter       Spring       Summer       Autumn       (NH)
      integer kkth, kkbh
      real ccd, ccd2, ccd3, ccd4
      common /clddat/ ccd(37,5),ccd2(37,5),ccd3(37,5),ccd4(37,5),
     &                kkth(37,5),kkbh(37,5)

C Local work arrays and variables
      integer l, k
      real cd(5), kth(5), kbh(5), date, rang, rsin1, rcos1, rcos2,
     &     theta

C Local data, functions etc
      external cldblk
c     logical start
c     data start / .true. /
c     save start

C Start code : ----------------------------------------------------------

C
c   Move this to initfs
C***      if ( start ) then
C***c       Rearrange the seasonal mean cloud data to allow interpolation
C***c       Define the amplitudes of the mean, annual and semi-annual cycles
C***	call reset(ccd,ccd2,ccd3,ccd4,37*5)
C***	start = .false.
C***      end if

c     Convert time to day number
      date = mod( real(mins)/1440., year)
      rang = tpi*(date-rlag)/year
      rsin1 = sin(rang)
      rcos1 = cos(rang)
      rcos2 = cos(2.0*rang)
C
      theta=90.-alat*180./pi
      l = (theta +2.5)/5. + 1.0
      do 5 k = 1,5
        cd(k)=ccd(l,k)+rsin1*ccd2(l,k)+rcos1*ccd3(l,k)+rcos2*ccd4(l,k)
        kth(k)=kkth(l,k)
        kbh(k)=kkbh(l,k)
    5 continue
      ch =cd(2)
      cm =cd(3)
      cl =cd(4)
      ich=kth(2)/2-1
      icm=kth(3)/2-1
      ict=kth(4)/2-1
      icb=kbh(4)/2-1
      return
      end
