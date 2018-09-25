c nsib comparison is modified to lsm_type = "nsib "
c AJA 2009/01/22
c
c Removed unnecessary "include 'MACHINE.f'"
c SJP 2001/11/22
c
c $Log: zerogi.f,v $
c Revision 1.26  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.25  1997/12/19  02:03:19  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.24  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.23  1996/12/20  01:36:42  ldr
c Tidy-ups to EAK snow stuff.
c
c Revision 1.22  1996/10/24  01:02:25  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.21  1996/03/21  03:18:26  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.20  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.19  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.18  1994/12/15  06:50:42  ldr
c Remove unnecessary divide of edbar, which is not used currently in the code.
c
c Revision 1.17  94/08/08  17:20:42  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.16  93/12/17  15:31:12  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.15  93/12/06  16:54:58  ldr
c Changes to V4-4-45l to get T63 running on SGI (Seca).
c 
c Revision 1.14  93/10/15  14:15:53  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c 
c Revision 1.13  93/10/07  12:10:16  ldr
c Move machine parameter statement into new include file MACHINE.f.
c 
c Revision 1.12  93/10/07  11:33:55  ldr
c Correction to SGI code in Hal's revision 1.10.
c 
c Revision 1.11  93/10/05  13:04:38  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.10  93/08/19  15:11:11  ldr
c Minor cosmetic changes.
c 
c Revision 1.9  93/07/14  14:51:00  ldr
c  ECMWF implicit treatment of vorticity eqn is now an option (HBG).
c 
      subroutine zerogi(ind,exptyp)

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ind
      character*50 exptyp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
c     include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'FLDMRI.f'
      include 'FLDRI.f'
c     include 'WORKRID.f'
c     include 'TIMEX.f'
      include 'WORKA.f'
      character*50 header2
      common/rshead/header2
      common/ubplm/utr(nl,lw),ur(lw,nl)

C Local work arrays and variables

C Local data, functions etc
c     include 'LLMAX.f'  !Statement function

C Start code : ----------------------------------------------------------

C**** MANIPULATION ROUTINE USED FOR PHYS & DYNM (almost redundant)

      if(ind.ne.1)go to 20
C**** IND=1 ONLY : BEFORE PHYS LOOP

      do 114 ll=1,lw 
      do 114 k=1,nl
 114  utr(k,ll)=0.0

      return

C**** IND=2 ONLY : BEFORE DYNM LOOP
   20 continue

      do 152 k=1,nl
      do 152 ll=1,lw
  152 ur(ll,k)=utr(k,ll)

c On raw restart with new CHEN temperature variable, set tau-1 temperature
c to equal tau value.

      if(chenflag.and.exptyp(42:45).ne.'CHEN')then
        do k=1,nl
          do mm=1,mw
            do ll=1,lw
              temr(ll,mm,k)=ter(ll,mm,k)
              temi(ll,mm,k)=tei(ll,mm,k)
            enddo
          enddo
        enddo
        exptyp(42:45)='CHEN'
        print*,'zerogi: Finished initializing CHEN temperature variable'
      endif

c On raw restart with EAK's new 9 soil types, set the header2
c New soil data is initialized in surfa if header2='     '

      if((lsm_type .eq. "nsib ").and.header2(38:42).ne.'SOIL9')then
        header2(38:42)='SOIL9'
        print*,'New NSIB 9 soil types initialized in surfa.'
      endif

c Similar treatment for initialization of snow albedo data

      if((lsm_type .eq. "nsib ").and.header2(44:48).ne.'SNALB')then
        header2(44:48)='SNALB'
        print*,'New NSIB snow albedo data initialized in surfset.'
      endif

c Similar treatment for initialization of snow albedo data

      if((lsm_type .eq. "nsib ").and.header2(50:50).ne.'M')then
        header2(50:50)='M'
        print*,'New NSIB soil moisture, temperature and snow data initia
     &lized in surfa.'
        print*,'NB: Model needs 2 years to equilibrate for this scheme.'
      endif
      
      pbar=prr(1,1)/sq2

      return
      end
