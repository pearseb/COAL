c $Log: radsgrg.f,v $
c Revision 1.2  2001/02/22 05:34:41  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine radsgrg(ipass,lg,solarfit,nrad,coszro2,taudar2,tg,
     &  sg,rg,blg)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ipass
      integer lg
      logical solarfit
      integer nrad
      real coszro2(ln2)
      real taudar2(ln2)
      real tg(ln2)
      real sg(ln2)
      real rg(ln2)
      real blg(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'MASIV3.f'
      include 'PMASIV3.f'

C Local work arrays and variables
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** Initialize surface radiation fields sg and rg

      IF(ipass.eq.1)THEN

      if (solarfit.and.nrad.ne.0) then
c**** Calculate the solar using the saved amplitude.
         do mg=1,ln2
            sg(mg)=sgamp(mg,lg)*coszro2(mg)*taudar2(mg)
         enddo
      else
         do mg=1,ln2
            sg(mg)=sgsav(mg,lg)
         enddo
      end if
      do  mg=1,ln2
         blg(mg)=stefbo*tg(mg)**4
         rg(mg)=rgsav(mg,lg)+blg(mg)
      enddo

      ELSE ! ipass=2

c---- the tg below is taken from ptg in radin.f
c---- the sg,rg,blg below are returned as psg,prg,pblg in radin.f

      if (solarfit.and.nrad.ne.0) then
c**** Calculate the solar using the saved amplitude.
         do mg=1,ln2
            sg(mg)=psgamp(mg,lg)*coszro2(mg)*taudar2(mg)
         enddo
      else
         do mg=1,ln2
            sg(mg)=psgsav(mg,lg)
         enddo
      endif
      do  mg=1,ln2
         blg(mg)=stefbo*tg(mg)**4
         rg(mg)=prgsav(mg,lg)+blg(mg)
      enddo

      END IF

      return
      end
