c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Modified to call ORBPAR at the start of each run, in order to obtain the
c Earth's orbital parameters. Also modified to obtain CSOLAR via /orbrad/,
c rather than it being declared here as a parameter. CSOLAR is now specified
c in W/m^2, rather than ly/min.
c SJP 2001/12/23.
c 
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: initfs.f,v $
c Revision 1.22  2001/02/12 05:39:55  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.21  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.20  1998/12/10  00:55:54  ldr
c HBG changes to V5-1-21
c
c Revision 1.19  1997/12/17  23:22:58  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.18  1996/10/24  01:02:58  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.17  1996/06/13  02:07:00  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.16  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.14.3.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.15  1994/08/08  17:21:37  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.14  94/03/23  09:50:47  mrd
c Added new orbital calculation suitable for paleoclimate work.
c 
c Revision 1.13  93/08/19  15:08:31  ldr
c Minor cosmetic changes.
c 
c Revision 1.12  93/08/03  16:31:06  ldr
c Simplify common block radsav to avoid misaligned data on SGI 64 bit.
c 
c Revision 1.11  93/07/22  09:22:00  mrd
c Removed saving of purely diagnostic radiative quanitities. Output 
c diagnostics for these are now accumulated only on radiation steps.
c 
c Revision 1.10  93/06/23  11:42:41  mrd
c Generalised number of levels in o3dat
c 
c Revision 1.9  92/12/09  14:43:40  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.8  92/12/01  14:03:28  ldr
c Tidied up and moved HBG's new clddia initialization to initfs.
c 
c Revision 1.7  92/11/30  14:31:45  mrd
c Cleaned up solarfit code
c
c     INPUT/OUTPUT
c     Input:   from common/clddiab in this subroutine
c                  cldlev - cloud levels
c
c              from common/orbrad in ORBRAD.f
c                  bpyear - epoch, in years before 1950
c                  csolar - solar constant, in W/m^2
c
c              from common/radisw2 in RADISW.f
c                  rrvco2 - CO2 volume mixing ratio 
c
c              from common/radnml in this subroutine
c                  dcflag - if T use diagnostic clouds, if F use
c                           fixed clouds, default F
c
c              from arguments
c                  icall - if zero, then at start of run
c
c     In/Out:  from common/clddiab in this subroutine
c                  icld - k indices for low, mid and high cloud 
c
c              from common/cnsta in CNSTA.f
c                  sig - sigma values
c
c              from common/radsav in this subroutine
c                  alat - latitude (from -pi/2 to pi/2)
c                  dlt  - declination of sun
c                  slag - apparent sun lag angle (west of mean sun is plus)
c                  coszm - cosine(solar zenith angle)
c                  fjd - day number, runs 1-365
c                  ha - hour angle of sun
c                  r1 - sun earth distance relative to standard distance
c                  taudam - fraction of day that is daylight
c
c              from common/timex in TIMEX.f
c                  mins  - current model time in mins
c                  mstep - timestep in minutes 
c                  nrad - no timesteps between calls to the radiation scheme,
c                         set to zero for no radiation, default 4
c
c     Output:  from common/radisw2 in RADISW.f
c                  rrco2  - mass mixing ratio (g/g) of CO2
c                  ssolar - solar constant (at present,in ly/min). 

      subroutine initfs(icall)

      implicit none

!$OMP THREADPRIVATE ( /CLDCOM/ )
!$OMP THREADPRIVATE ( /KDACOM/ )
!$OMP THREADPRIVATE ( /LWOUT/ )
!$OMP THREADPRIVATE ( /RADISW/ )
!$OMP THREADPRIVATE ( /RDFLUX/ )
!$OMP THREADPRIVATE ( /CLRFLX/ )
!$OMP THREADPRIVATE ( /SRCCOM/ )
!$OMP THREADPRIVATE ( /SWOCOM/ )
!$OMP THREADPRIVATE ( /TFCOM/ )

c To initialize for FS radiation scheme. (LDR 4/92)
c This code taken from radfs.f - better here for SGI.
c Called with argument 0 for start of run, argument 1 for once each timestep

C Global parameters
      include 'HCON.f'
      include 'PARAMS.f'
      include 'RDPARM.f'

C Argument list
      integer icall

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'CLDCOM.f'
      include 'KDACOM.f'
      include 'LWOUT.f'
      include 'RADISW.f'
      include 'RDFLUX.f' ! includes CLRFLX
      include 'SRCCOM.f'
      include 'SWOCOM.f'
      include 'TFCOM.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'GAUSL.f'
      include 'ORBRAD.f'
      include 'TIMEX.f'
      include 'CNSTA.f'

      real alat,ha,coszm,taudam,fjd,r1,dlt,slag,dhr
      common /radsav/ alat(lat*2),ha(lat*2),coszm(lat*2),taudam(lat*2),
     &      fjd,r1,dlt,slag,dhr

c Stuff from o3set
      real dduo3n,ddo3n2,ddo3n3,ddo3n4
      common /o3dat/ dduo3n(37,nl),ddo3n2(37,nl),ddo3n3(37,nl),
     &               ddo3n4(37,nl)

c Stuff from cldset
      real ccd,ccd2,ccd3,ccd4
      integer kkth,kkbh
      common /clddat/ ccd(37,5),ccd2(37,5),ccd3(37,5),ccd4(37,5),
     &                kkth(37,5),kkbh(37,5)

c Stuff from clddia
      integer icld
      real cldlev
      common/clddiab/ icld(2,3), cldlev(2,3)  !Communicates with clddia

C Local work arrays and variables
      real sigf(l)

      integer i
      integer i1
      integer j
      integer k
      integer*8 kount

      real alp
      real orbyear
      real lymin_to_wm2
      parameter (lymin_to_wm2 = 41840.0 / 60.0)

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(icall.eq.0)then  !At start of run

c Call ORBPAR to derive the Earth's orbital parameters
        orbyear = 1950.0 - real(bpyear)
        call orbpar(orbyear, ec, oblqty, peril)
        write (*, *) ""
        write (*, *) "Epoch = ", bpyear, " years BP"
        if (bpyear .lt. 1950) then
          write (*, *) "Calendar year = ", 1950-bpyear, " AD"
        else
          write (*, *) "Calendar year = ", bpyear-1949, " BC"
        end if
        write (*, *) ""
        write (*, *) "Solar constant = ", csolar, " W/m^2"
        write (*, *) ""
        write (*, *) "Orbital parameters:"
        write (*, *) "Eccentricity = ", ec
        write (*, *) "Obliquity    = ", oblqty, " degrees"
        write (*, *) "Longitude of perihelion = ", peril, " degrees"
        write (*, *) ""

        call hconst
 	call co2_read(sig)
        call table
        rrco2=rrvco2*ratco2mw
c     Calculate latitude (from PI/2 to -PI/2)
        do 100 j=1,lat
          alat(j) = asin(coa(j))
          alat(2*lat+1-j) = - alat(j)
 100    continue

        if(amipo3)then
c         AMIP2 ozone
          call o3read_amip
          print *,'AMIP2 ozone input'
        else
c Stuff from o3set
c       Rearrange the seasonal mean O3 data to allow interpolation
c       Define the amplitudes of the mean, annual and semi-annual cycles
          call o3_read(sig)
          call reset(dduo3n,ddo3n2,ddo3n3,ddo3n4,37*nl)
        endif

c Stuff from cldset
c       Rearrange the seasonal mean cloud data to allow interpolation
c       Define the amplitudes of the mean, annual and semi-annual cycles
	if(.not.dcflag) call reset(ccd,ccd2,ccd3,ccd4,37*5)

c Stuff from clddia
C
C---------------------------------------------------------------------*
C START COMPUTATION                                                   *
C---------------------------------------------------------------------*
C
C---------------------------------------------------------------------*
C ON FIRST CALCULATION CALCULATE THE HIGH, MIDDLE AND LOW             *
C CLOUD LEVEL LIMITS                                                  *
C---------------------------------------------------------------------*
C
        do k=1,l
          sigf(k) = sig(l+1-k)
        enddo
        i1 = 1
        do 25 j = 1,3
          do 20 i = i1,lm1
            if(sigf(i).lt.cldlev(1,j)) goto 20
            if(icld(1,j).eq.0) icld(1,j) = i
            if(sigf(i).lt.cldlev(2,j)) icld(2,j) = i
 20       continue
          i1 = icld(2,j) + 1
 25     continue
C---------------------------------------------------------------------*
 

      return

c******************************************************************************

      else      !Once each timestep

        kount = mins/mstep
        if((.not.SCM.and.(mod(kount,int(nrad, 8)).ne.0_8)).or.
     &     (SCM.and.(nsteps.gt.0).and.(mod(kount,int(nrad, 8)).ne.0_8)))
     &    then

c          Use stored radiation values
           fjd = fjd + float(mstep)/1440.
        else

c          Calculate solar parameters for radiation calculation
           fjd = real(mod(mins,525600_8))/1440.
           call solargh(fjd,r1,dlt,alp,slag,lat2,alat,ha,coszm,taudam)

c          Modify solar constant according to orbital radius.
           ssolar = csolar / (r1**2)
           ssolar = ssolar / lymin_to_wm2
        end if
      end if

      return
      end
