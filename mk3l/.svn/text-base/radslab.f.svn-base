c HFACA removed from /FICECON/, as this array is never used.
c SJP 2003/06/20
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: radslab.f,v $
c Revision 1.3  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.2  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine radslab(lg,sgold,rgold,fg,eg,hfrm,
     &                      psgold,prgold,pfg,peg,
     &                      flx,surfb)

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg
      real sgold(ln2)
      real rgold(ln2)
      real fg(ln2)
      real eg(ln2)
      real hfrm(ln2)
      real psgold(ln2)
      real prgold(ln2)
      real pfg(ln2)
      real peg(ln2)
      real flx(ln2)
      real surfb(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'

      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)

C Global data blocks
      include 'DATICE1.f'
      include 'FEWFLAGS.f'

      real ochf,occur,g50dt
      integer khf
      common/curent/ochf(ln2,lat),occur(ln2,lat),khf(ln2,lat)
     &,g50dt(ln2,lat)

      real flxia
      common/ficecon/flxia(lat,2)

C Local work arrays and variables
      real qfa(ln2)

      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

      do mg=1,lon
c.... flx() is any sub-ice flux
        flx(mg)=flxia(lg,1)+hoice(mg,lg)
        flx(mg+lon)=flxia(lg,2)+hoice(mg+lon,lg)
      enddo

      do mg=1,ln2
c.... qfa() is the timestep q-flux, surfb() is the timestep surface flux
        qfa(mg)=0.0
        surfb(mg)=sgold(mg)-rgold(mg)-fg(mg)-eg(mg)
      enddo

      do mg=1,ln2
        if (sea(mg))  qfa(mg) = g50dt(mg,lg) - surfb(mg)
        if (mlo(mg))  qfa(mg) = hfrm(mg)
        if (cice(mg)) qfa(mg) = flx(mg)
      enddo

      if(leads)then
        do mg=1,ln2
        if (cice(mg)) surfb(mg) = surfb(mg)*(1.0-pl(mg))
     &     +(psgold(mg)-prgold(mg)-pfg(mg)-peg(mg))*pl(mg)
        enddo
      endif

      do mg=1,ln2
        occur(mg,lg)=occur(mg,lg)+qfa  (mg)
        ochf (mg,lg)=ochf (mg,lg)+surfb(mg)
      enddo

      return
      end
