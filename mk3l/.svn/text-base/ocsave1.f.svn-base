c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c$Log: ocsave1.f,v $
cRevision 1.14  2001/02/22 05:34:45  rot032
cChanges from HBG to complete concatenation of NH/SH latitudes
c
cRevision 1.13  2000/06/20 02:08:35  rot032
cHBG changes to V5-3-3, mainly for coupled model
c
cRevision 1.12  1999/05/20 06:23:57  rot032
cHBG changes to V5-2
c
c Revision 1.11  1998/12/10  00:55:57  ldr
c HBG changes to V5-1-21
c
c Revision 1.10  1997/12/17  23:23:00  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.9  1995/09/07  06:45:24  ldr
c Add task commons.
c
c Revision 1.8  1993/11/03  11:44:24  ldr
c Changes to V4-4-24l from HBG for coupled model
c
c Revision 1.7  93/06/16  14:18:28  ldr
c HBG tidy ups.
c 
c Revision 1.6  93/02/03  12:44:51  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
c Revision 1.5  92/12/10  09:55:15  ldr
c Minor fixes.
c 
c Revision 1.4  92/12/09  15:24:13  ldr
c Put /datice in include file and change name to datice1.
c 
c Revision 1.3  92/12/09  12:20:06  ldr
c Some changes to get coupled common blocks into include files.
c 
c Revision 1.2  92/10/08  17:22:44  ldr
c SPO's changes to sea-ice and MLO scheme for coupled model.
c 
c Revision 1.1  91/05/20  15:06:05  ldr
c Initial revision
c 
c Revision 1.1  91/05/20  14:53:23  ldr
c Initial revision
c 
      subroutine ocsave1(lg,sg,psg,pl,il,lil)

c Routine to save solar data for ogcm

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg
      real sg(ln2)
      real psg(ln2)
      real pl(ln2)
      integer il(ln2)
      integer lil

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'

      real drunof,dsis,dsisb,dsfw,athf,athfa,atsi
      common/aforce/drunof(ln2,lat),dsis(ln2,lat),dsisb(ln2,lat)
     & ,dsfw(ln2,lat),athf(ln2,lat),athfa(ln2,lat),atsi(ln2,lat)

C Local work arrays and variables
      real sgx(ln2)

      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------


       do 10 mg=1,ln2
        sgx(mg)=sg(mg)
 10    continue

      IF(leads)THEN

       if(lil.gt.0)then
        do 30 mg=1,ln2   
         if(il(mg).gt.0)then
c..       sgx(mg)=pl(mg)*psg(mg)
c.. Note - ice model uses surface energy balance in its mixed layer
          sgx(mg)=0.0
         endif
 30     continue
       endif

      END IF

      do 60 mg=1,ln2
        atsi(mg,lg)=sgx(mg)
   60 continue

      return
      end
