c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c $Log: o3set.f,v $
c Revision 1.9  1998/12/10 00:55:32  ldr
c HBG changes to V5-1-21
c
c Revision 1.8  1996/10/24  01:03:05  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.7  1992/12/09  14:44:03  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c
c Revision 1.6  92/10/20  16:20:26  mrd
c Fixed error in date calculation (introduced in previous change).
c 
c Revision 1.5  92/09/16  16:02:18  mrd
c Changed date calculation to use the model year of exactly 365, 
c not 365.25 days.
c 
c Revision 1.4  92/05/11  15:13:40  ldr
c Put Include PARAMS.f in source file rather than in RDPARM.f, 
c to avoid nested includes.
c 
c Revision 1.3  92/04/16  16:37:30  ldr
c Moved initialization stuff to iitfs.
c 
c Revision 1.2  92/04/15  12:33:02  mrd
c Restructured radiation code include files and data input
c 
c Revision 1.1  91/02/22  16:37:45  ldr
c Initial release V3-0
c 
      subroutine o3set(alat,mins,duo3n,sigma)
c
c  This routine interpolates in latitude and time to set the ozone 
c  amounts.
c  INPUT
c    ALAT    latitude (from -pi/2 to pi/2)
c    MINS    current model time in mins
c  OUTPUT
c    DUO3N   ozone mixing ratio
c

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'RDPARM.f'
      parameter (tpi=2.*pi,rlag=14.8125)
      parameter (year=365)

C Argument list
      real alat
      integer*8 mins
      real duo3n(nl)
      real sigma(nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
c          winter       spring       summer       autumn       (nh)
      common /o3dat/ dduo3n(37,nl),ddo3n2(37,nl),ddo3n3(37,nl),
     &     ddo3n4(37,nl)

C Local work arrays and variables

C Local data, functions etc
c     logical start
c     data start / .true. /
c     save start

C Start code : ----------------------------------------------------------

c This moved to initfs
c
C***      if ( start ) then
C***c       Rearrange the seasonal mean o3 data to allow interpolation
C***c       Define the amplitudes of the mean, annual and semi-annual cycles.
C***	call o3_read(sigma)
C***	call reset(dduo3n,ddo3n2,ddo3n3,ddo3n4,37*nl)
C***	start = .false.
C***      end if

c     Convert time to day number
c     date = amod( float(mins)/1440., year)
c     Use year of exactly 365 days
      date = real(mod(mins,525600_8))/1440.
      rang = tpi*(date-rlag)/year
      rsin1 = sin(rang)
      rcos1 = cos(rang)
      rcos2 = cos(2.0*rang)
c
      theta=90.-alat*180./pi
      ilat = theta/5.
      angle = 5 * ilat
      than = (theta-angle)/5.
      ilat = ilat+1
      do 1 m = 1,nl
      do3  = dduo3n(ilat,m) + rsin1*ddo3n2(ilat,m) 
     &     + rcos1*ddo3n3(ilat,m) + rcos2*ddo3n4(ilat,m)
      do3p = dduo3n(ilat+1,m) + rsin1*ddo3n2(ilat+1,m) 
     &     + rcos1*ddo3n3(ilat+1,m) + rcos2*ddo3n4(ilat+1,m)
      duo3n(m)=do3+than*(do3p-do3)
    1 continue
      return
      end
