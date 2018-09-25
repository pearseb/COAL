c Modified by HBG in accordance with revised treatment of SNOWMI in surfupl.f.
c SJP 2005/02/24
c
c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: radcoupl.f,v $
c Revision 1.3  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.2  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine radcoupl(lg,lil,il,sgold,rgold,fg,eg,runoff,
     &             condx,evap,flx,sicold,plold,siced,sublmi,
     &             pcondx,pevap,snowmi)

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg
      integer lil
      integer il(ln2)
      real sgold(ln2)
      real rgold(ln2)
      real fg(ln2)
      real eg(ln2)
      real runoff(ln2)
      real condx(ln2)
      real evap(ln2)
      real flx(ln2)
      real sicold(ln2)
      real plold(ln2)
      real siced(ln2)
      real sublmi(ln2)
      real pcondx(ln2)
      real pevap(ln2)
      real snowmi(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'

      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)

C Global data blocks
      include 'FEWFLAGS.f'
      real drunof,dsis,dsisb,dsfw,athf,athfa,atsi
      common/aforce/drunof(ln2,lat),dsis(ln2,lat),dsisb(ln2,lat)
     & ,dsfw(ln2,lat),athf(ln2,lat),athfa(ln2,lat),atsi(ln2,lat)
      real flxia
      common/ficecon/flxia(lat,2)

C Local work arrays and variables
      integer mg, ns

C Local data, functions etc

C Start code : ----------------------------------------------------------

      do 416 mg=1,ln2
        athf  (mg,lg)=sgold(mg)-rgold(mg)-fg(mg)-eg(mg)
        athfa (mg,lg)=athf(mg,lg)
        if(land(mg))then
          athf(mg,lg)=0.0
          athfa(mg,lg)=0.0
        end if
        drunof(mg,lg)=runoff(mg)*0.001
        dsis  (mg,lg)=0.0
        dsisb (mg,lg)=0.0
        dsfw  (mg,lg)=(condx(mg)-evap(mg))*0.001
  416 continue

      If(leads)then

      If(lil.gt.0)Then
      do 417 mg=1,ln2
        ns = (mg-1) / lon + 1
        if(il(mg).gt.0)then
        athf  (mg,lg)=-flx(mg)
        athfa (mg,lg)=athf(mg,lg) + flxia(lg, ns)
        tflux (mg)   =tflux(mg)-flx(mg)
        dsis  (mg,lg)=sicold(mg)*(1.0-plold(mg))
     &                    -(siced(mg)+sublmi(mg))*(1.0-pl(mg))
        sflux (mg)   =sflux(mg)+dsis(mg,lg)
        dsisb (mg,lg)=sublmi(mg)*(1.0-pl(mg))
        dsfw  (mg,lg)=(pcondx(mg)-pevap(mg))*0.001*pl(mg)
     &                +snowmi(mg) ! snowmi already factored by (1-pl)
      end if
  417 continue
      End If

      Else

      do 418 mg=1,ln2
      ns = (mg-1) / lon + 1
      if((sicold(mg).gt.0.0).or.cice(mg))then
        athf  (mg,lg)=-flx(mg)
        athfa (mg,lg)=athf(mg,lg) + flxia(lg, ns)
        dsis  (mg,lg)=sicold(mg)-(siced(mg)+sublmi(mg))
        dsisb (mg,lg)=sublmi(mg)
        dsfw  (mg,lg)=snowmi(mg)
      end if
  418 continue

      End if

      return
      end
