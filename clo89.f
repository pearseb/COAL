c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/08
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c $Log: clo89.f,v $
c Revision 1.10  1998/12/10 00:56:02  ldr
c HBG changes to V5-1-21
c
c Revision 1.9  1997/12/17  23:23:07  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.8  1996/06/13  02:05:45  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.7  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.6  93/08/19  15:07:45  ldr
c Minor cosmetic changes.
c 
c Revision 1.5  93/06/23  14:30:28  mrd
c Made padding of vtemp common block resolution independent
c 
c Revision 1.4  92/05/11  15:13:19  ldr
c Put Include PARAMS.f in source file rather than in RDPARM.f, 
c to avoid nested includes.
c 
c Revision 1.3  92/04/16  16:40:57  ldr
c 
c Reinstated common vtemp (needed for SGI.)
c 
c Revision 1.2  92/04/15  11:59:08  mrd
c Restructured radiation code include files and data input
c 
c Revision 1.1  91/02/22  16:37:01  ldr
c Initial release V3-0
c 
      subroutine clo89
 
!$OMP THREADPRIVATE ( /CLDCOM/ )
!$OMP THREADPRIVATE ( /RADISW/ )
!$OMP THREADPRIVATE ( /VTEMP/ )

c     subroutine clo88 computes cloud transmission functions for the
c  longwave code,using code written by bert katz (301-763-8161).
c  and modified by dan schwarzkopf in december,1988.
c                inputs:          (common block)
c      camt,ktop,kbtm,nclds,emcld   radisw
c                output:  
c      cldfac                       cldcom
c 
c          called by:      radmn or model routine 
c          calls    : 
c 

C Global parameters
      include 'PARAMS.f'
      include 'RDPARM.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'CLDCOM.f'
      include 'RADISW.f'
      common /vtemp/ tempc(lp1,lp1,imax),cldfip(lp1,lp1)
      common /vtemp/ dummy(imax*(l*l+15*l+18)-lp1*lp1)

C Global data blocks

C Local work arrays and variables
      integer i, j

C Local data, functions etc

C Start code : ----------------------------------------------------------

      do 11 ip=1,imax
      if (nclds(ip).eq.0) then
CSJP        do 29 i=1,lp1*lp1
CSJP        tempc(i,1,ip)=1.
        do 29 j = 1, lp1
          do 29 i = 1, lp1
            tempc(i, j, ip) = 1.0
29      continue
      endif
      if (nclds(ip).ge.1) then
          xcld=1.-camt(ip,2)*emcld(ip,2)
           k1=ktop(ip,2)+1
           k2=kbtm(ip,2)
CSJP          do 31 i=1,lp1*lp1
CSJP              cldfip(i,1)=1.
           do 31 j = 1, lp1
             do 31 i = 1, lp1
               cldfip(i, j) = 1.0
31        continue
          do 41 k=k1,lp1
          do 41 kp=1,k2
               cldfip(kp,k)=xcld
41        continue
          do 43 k=1,k2
          do 43 kp=k1,lp1
              cldfip(kp,k)=xcld
43        continue
CSJP            do 61 i=1,lp1*lp1
CSJP          tempc(i,1,ip)=cldfip(i,1)
          do 61 j = 1, lp1
            do 61 i = 1, lp1
              tempc(i, j, ip) = cldfip(i, j)
61        continue
      endif
      if (nclds(ip).ge.2) then
        do 21 nc=2,nclds(ip)
          xcld=1.-camt(ip,nc+1)*emcld(ip,nc+1)
           k1=ktop(ip,nc+1)+1
           k2=kbtm(ip,nc+1)
CSJP          do 32 i=1,lp1*lp1
CSJP              cldfip(i,1)=1.
           do 32 j = 1, lp1
             do 32 i = 1, lp1
               cldfip(i, j) = 1.0
32        continue
          do 42 k=k1,lp1
          do 42 kp=1,k2
               cldfip(kp,k)=xcld
42        continue
          do 44 k=1,k2
          do 44 kp=k1,lp1
              cldfip(kp,k)=xcld
44        continue
CSJP            do 62 i=1,lp1*lp1
CSJP          tempc(i,1,ip)=tempc(i,1,ip)*cldfip(i,1)
          do 62 j = 1, lp1
            do 62 i = 1, lp1
              tempc(i, j, ip) = tempc(i, j, ip) * cldfip(i, j)
62        continue
21        continue
      endif
11    continue
      do 70 ip=1,imax
CSJP      do 70 i=1,lp1*lp1
CSJP         cldfac(ip,i,1)=tempc(i,1,ip)
        do 70 j = 1, lp1
          do 70 i = 1, lp1
            cldfac(ip, i, j) = tempc(i, j, ip)
70      continue
      return
      end 
