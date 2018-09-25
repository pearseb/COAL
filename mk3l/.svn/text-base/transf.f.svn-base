c $Log: transf.f,v $
c Revision 1.7  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.6  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.5  1996/10/24  01:03:21  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1993/10/05  13:07:48  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.3  92/12/09  14:44:46  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.2  91/03/13  13:01:16  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:38:21  ldr
c Initial release V3-0
c 
      subroutine transf

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FLDRI.f'
      include 'FLDMRI.f'
      include 'RMGRID.f'

C Local work arrays and variables
      integer k
      integer lg
      integer ll
      integer mg
      integer mm

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**** TO REPLACE T-1 FIELDS BY T FIELDS FOR FWD STEPS.
      do 10 k=1,nl
      do 10 mm=1,mw
      do 10 ll=1,lw
      temr(ll,mm,k)=ter(ll,mm,k)
      temi(ll,mm,k)=tei(ll,mm,k)
      psimr(ll,mm,k)=psir(ll,mm,k)
      psimi(ll,mm,k)=psii(ll,mm,k)
      xhimr(ll,mm,k)=xhir(ll,mm,k)
   10 xhimi(ll,mm,k)=xhii(ll,mm,k)

      do 11 lg=1,lat
      do 11 k=1,nl
      do 11 mg=1,ln2
   11 rmmg(mg,k,lg)=rmg(mg,k,lg)

      do 12 mm=1,mw
      do 12 ll=1,lw
      pmr(ll,mm)=prr(ll,mm)
   12 pmi(ll,mm)=pri(ll,mm)

      return
      end
