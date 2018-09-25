c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: radstres.f,v $
c Revision 1.2  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine radstres(lg,lil,il,taux,tauy,ptaux,ptauy)

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg
      integer lil
      integer il(ln2)
      real taux(ln2)
      real tauy(ln2)
      real ptaux(ln2)
      real ptauy(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'

C Global data blocks
      include 'FEWFLAGS.f'

      real ttaux,ttauy
      common/aticstr/ttaux(ln2,lat),ttauy(ln2,lat) !MKS Kgm/m/sec**2

      real otaux,otauy
      common/atocstr/otaux(ln2,lat),otauy(ln2,lat) !    dynes/cm**2

C Local work arrays and variables
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c***********************************************************************
c.... These are the stresses for the ice model.
c.... If ice is present at a point, taux & tauy will be atmos/ice stress
      do mg=1,ln2
         ttaux(mg,lg)=ttaux(mg,lg)+taux(mg)
         ttauy(mg,lg)=ttauy(mg,lg)+tauy(mg)
      enddo


c***********************************************************************
c.... Gather ocean stresses. Land values will be present.
c.... These stresses must be expanded over land areas,
c.... and then transferred to the ocean (U/V) grid which does
c.... not correspond to the atmos stress grid.
c.... This will be done by routine ocntau.f
c.... Ocean model (MOM) needs stresses in dynes/cm**2
c.... 1 dyne = 1 gm cm/sec**2
c.... AGCM stresses in Kgm/m/sec**2 ==> factor of 10

      IF(leads)THEN

         if(lil.gt.0)then
         do 435 mg=1,ln2
            if(il(mg).eq.0)then
            otaux(mg,lg)=10.*taux(mg)
            otauy(mg,lg)=10.*tauy(mg)
            else
            otaux(mg,lg)=10.*ptaux(mg)
            otauy(mg,lg)=10.*ptauy(mg)
c..  note that the remaining ocean stress (from ice/ocean drag) will
c..  be added in ice model (see dynice.f) weighted by (1-pl)
c..  and must allow for 2 atmos steps per 1 ice dynamics step.
            endif
 435     continue
         else
         do 430 mg=1,ln2
            otaux(mg,lg)=10.*taux(mg)
            otauy(mg,lg)=10.*tauy(mg)
 430     continue
         endif

      ELSE

         do 438 mg=1,ln2
            if(cice(mg))then
              otaux(mg,lg)=0.0
              otauy(mg,lg)=0.0
            else
              otaux(mg,lg)=10.*taux(mg)
              otauy(mg,lg)=10.*tauy(mg)
            endif
 438     continue

      ENDIF

      return
      end
