c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c $Log: nestarc.f,v $
c Revision 1.16  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.15  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.14  1996/10/24  01:03:04  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.13  1996/03/21  03:18:55  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.12  1994/08/08  17:21:44  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.11  93/10/05  13:06:39  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.10  93/08/19  15:08:42  ldr
c Minor cosmetic changes.
c 
c Revision 1.9  93/06/28  14:38:53  ldr
c Tidy up data statements for blocks /nest and /brian1 to keep SGI happy.
c 
c Revision 1.8  92/12/09  15:06:58  ldr
c Put masiv4 in include file for ogcm compatibility.
c 
c Revision 1.7  92/12/09  14:44:01  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.6  92/06/04  12:03:15  ldr
c Fixed data statement so this will compile on SGI for R42.
c 
c Revision 1.5  92/01/29  15:40:52  ldr
c Fixed writing of soil temperatures.
c 
c Revision 1.4  91/05/15  11:08:53  ldr
c Put common /nest into include file.
c 
c Revision 1.3  91/03/13  12:59:10  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c
c     INPUT/OUTPUT
c     Input:  from common/masiv4 in MASIV4.f
c                 savegrid - global surface variables
c
c             from common/fldri in FLDRI.f
c                  tei - spectral temp, imaginary part
c                  ter - spectral temp, real part
c
c             from common/nest in NEST.f
c                 nest_alb - albedo
c                 nest_interval - mins between calls to nestarc
c
c             from common/rmgrid in RMGRID.f
c                 rmg - pressure weighted moistures, current timestep
c
c             from common/worka in WORKA.f
c                 tmean - preset isothermal temperature
c
c             from arguments
c                 ldays - month count of days  
c                 mins  - current model time in mins
c
c     In/Out: from common/fldri in FLDRI.f
c                 pi - spectral pressure, imaginary part
c                 pr - spectral pressure, real part
c                 psii - spectral stream function, imaginary part
c                 psir - spectral stream function, real part
c                 xhii - velocity potental, imaginary part
c                 xhir - velocity potental, real part
c
c             from common/nest in NEST.f
c                 nest_rain - rain since last archive in mm/day (gridpt form)
c
c
      subroutine nestarc(mins,ldays)

c Subroutine to archive spectral and gridpoint quantities needed for
c running the nested limited area model. (LDR 9/1990)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer*8 mins
      integer ldays

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FLDRI.f'
      include 'MASIV4.f'
      include 'NEST.f'
      include 'RMGRID.f'
      include 'WORKA.f'

C Local work arrays and variables
      real qout(lon,lat2),tout(lon,lat2),zout(lon,lat2)
      real albout(lon,lat2),wgout(lon,lat2),wbout(lon,lat2)
      real tb2out(lon,lat2),tb3out(lon,lat2),rainout(lon,lat2)
      complex cx(lw,mw)

      integer k
      integer lg
      integer ll
      integer mg
      integer mm

      real avgrn
      real tempi
      real tempr

C Local data, functions etc

C Start code : ----------------------------------------------------------

      write(20)mins,ldays

C**** WRITE P*

      do 70 mm=1,mw
      do 70 ll=1,lw
   70 cx(ll,mm)=cmplx(prr(ll,mm),pri(ll,mm))
      write(20)cx

C     Write PSI(flux)

      do 11 k=1,nl
        do 10 mm=1,mw
          do 10 ll=1,lw
   10       cx(ll,mm)=cmplx(psir(ll,mm,k),psii(ll,mm,k))
   11   write(20) cx

C**** WRITE XHI(flux)

      do 14 k=1,nl
        do 12 mm=1,mw
          do 12 ll=1,lw
   12       cx(ll,mm)=cmplx(xhir(ll,mm,k),xhii(ll,mm,k))
   14   write(20) cx

c WRITE TEMPERATURE (P* weighted)

      do 20 k=1,nl
        do 18 mm=1,mw
          do 18 ll=1,lw
            tempr=ter(ll,mm,k)+tmean(k)*prr(ll,mm)
            tempi=tei(ll,mm,k)+tmean(k)*pri(ll,mm)
   18       cx(ll,mm)=cmplx(tempr,tempi)
   20   write(20) cx

C**** OUTPUT GRID MOISTURE (gridpoint form)

      do 29 k=1,nl
          do 27 lg=1,lat
            do 27 mg=1,lon
C---- NOTE NOTE
C---- This version has NO pressure weighting on qout
C---- NOTE NOTE
              qout(mg,lg)=rmg(mg+lon,k,lg)
   27         qout(mg,lat2p-lg)=rmg(mg,k,lg)
        write(20)qout
   29 continue

c WRITE SURFACE GEOPOTENTIAL (gridpoint form)

      do 30 lg=1,lat
        do 30 mg=1,lon
          zout(mg,lg)=grav*savegrid(mg+lon,1,lg)
   30     zout(mg,lat2p-lg)=grav*savegrid(mg,1,lg)

      write(20) zout

c WRITE RAIN SINCE LAST ARCHIVE IN MM/DAY (gridpoint form)
c ALSO ZERO RAIN ARRAY

      avgrn=1440./nest_interval
      do 32 lg=1,lat
        do 32 mg=1,lon
          rainout(mg,lg)=nest_rain(mg+lon,lg)*avgrn
   32     rainout(mg,lat2p-lg)=nest_rain(mg,lg)*avgrn

      write(20)rainout

      do 34 lg=1,lat
        do 34 mg=1,ln2
   34     nest_rain(mg,lg)=0.0

c WRITE SURF TEMP (gridpoint form)

      do 40 lg=1,lat
        do 40 mg=1,lon
          tout(mg,lg)=savegrid(mg+lon,3,lg)
   40     tout(mg,lat2p-lg)=savegrid(mg,3,lg)
      write(20)tout

c WRITE SOIL TEMPERATURE (UPPER then LOWER)
c NO! This is lower then upper!!

      do 41 lg=1,lat
        do 41 mg=1,lon
          tb2out(mg,lg)=savegrid(mg+lon,9,lg)
          tb2out(mg,lat2p-lg)=savegrid(mg,9,lg)
          tb3out(mg,lg)=savegrid(mg+lon,11,lg)
   41     tb3out(mg,lat2p-lg)=savegrid(mg,11,lg)
      
      write(20) tb2out
      write(20) tb3out

c WRITE SOIL MOISTURE (TOP)

      do 50 lg=1,lat
        do 50 mg=1,lon
          wgout(mg,lg)=savegrid(mg+lon,8,lg)
   50     wgout(mg,lat2p-lg)=savegrid(mg,8,lg)

      write(20) wgout

c WRITE SOIL MOISTURE (LOWER)

      do 60 lg=1,lat
        do 60 mg=1,lon
          wbout(mg,lg)=savegrid(mg+lon,10,lg)
   60     wbout(mg,lat2p-lg)=savegrid(mg,10,lg)

      write(20) wbout

c WRITE ALBEDO

      do 65 lg=1,lat
        do 65 mg=1,lon
          albout(mg,lg)=nest_alb(mg+lon,lg)
   65     albout(mg,lat2p-lg)=nest_alb(mg,lg)

      write(20)albout

      return
      end
