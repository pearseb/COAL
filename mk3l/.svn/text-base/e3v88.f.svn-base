c Add a call to NINT to resolve a warning issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/08
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: e3v88.f,v $
c Revision 1.11  1998/12/10 00:55:56  ldr
c HBG changes to V5-1-21
c
c Revision 1.10  1997/12/17  23:22:59  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.9  1996/06/13  02:06:26  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.8  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.7  93/12/17  15:32:20  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.6  93/08/19  15:07:56  ldr
c Minor cosmetic changes.
c 
c Revision 1.5  93/06/23  14:30:31  mrd
c Made padding of vtemp common block resolution independent
c 
c Revision 1.4  92/05/11  15:13:27  ldr
c Put Include PARAMS.f in source file rather than in RDPARM.f, 
c to avoid nested includes.
c 
c Revision 1.3  92/04/16  16:41:39  ldr
c 
c Reinstated common vtemp (needed for SGI.)
c 
c Revision 1.2  92/04/15  12:21:23  mrd
c Restructured radiation code include files and data input
c 
c Revision 1.1  91/02/22  16:37:15  ldr
c Initial release V3-0
c 
c     subroutine e3v88 computes nearby layer transmissivities for 
c  h2o using a table lookup of the pre-computed e3 function 
c ( described in ref. (4)). 
c         inputs:                 (common blocks,args.) 
c       tv,av                      argument list
c       em3                        tabcom 
c          outputs: 
c       emv                        argument list
c 
c       called by  :    fst88 
c       calls      :    alog10h,alog10v 
 
      subroutine e3v88(emv,tv,av) 
 
!$OMP THREADPRIVATE ( /VTEMP/ )

C Global parameters
      include 'HCON.f'
      include 'PARAMS.f'
      include 'RDPARM.f'

C Argument list
      real emv(imax,llp1)
      real tv(imax,llp1)
      real av(imax,llp1)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common / vtemp / it(imax,llp1)
      common / vtemp / fyo(imax,llp1),
     &                 ww1(imax,llp1),
     &                 tval(imax,llp1),dt(imax,llp1),ww2(imax,llp1),
     &                 uval(imax,llp1),du(imax,llp1),
     &                 fxo(imax,llp1),tmp3(imax,llp1)
      common /vtemp/ dummy(imax*(2*l*l-3*l+9))

C Global data blocks
      include 'TABCOM.f'
      dimension em3v(5040)
      equivalence (em3v(1),em3(1,1))

C Local work arrays and variables

C Local data, functions etc
      integer i, j

C Start code : ----------------------------------------------------------

c---the following loop replaces a double loop over i (1-imax) and
c   k (1-llp1)
CSJP      do 203 i=1,imax*llp1
CSJP        fxo(i,1)=aint(tv(i,1)*hp1)
CSJP        tmp3(i,1)=log10(av(i,1))+h16e1
CSJP        dt(i,1)=tv(i,1)-ten*fxo(i,1)
CSJP        fyo(i,1)=aint(tmp3(i,1)*ten)
CSJP        du(i,1)=tmp3(i,1)-hp1*fyo(i,1)
      do 203 j = 1, llp1
        do 203 i = 1, imax
          fxo(i, j) = aint(tv(i, j)*hp1)
          tmp3(i, j) = log10(av(i, j)) + h16e1
          dt(i, j) = tv(i, j) - ten * fxo(i, j)
          fyo(i, j) = aint(tmp3(i, j)*ten)
          du(i, j) = tmp3(i, j) - hp1 * fyo(i, j)
c---obtain index for table lookup; this value will have to be
c   decremented by 9 to account for table temps starting at 100K.
CSJP        it(i,1)=fxo(i,1)+fyo(i,1)*h28e1
CSJP        ww1(i,1)=ten-dt(i,1)
CSJP        ww2(i,1)=hp1-du(i,1)
CSJP        emv(i,1)=ww1(i,1)*ww2(i,1)*em3v(it(i,1)-9)+
CSJP     &           ww2(i,1)*dt(i,1)*em3v(it(i,1)-8)+
CSJP     &           ww1(i,1)*du(i,1)*em3v(it(i,1)+19)+
CSJP     &           dt(i,1)*du(i,1)*em3v(it(i,1)+20)
          it(i, j) = nint(fxo(i, j) + fyo(i, j) * h28e1)
          ww1(i, j) = ten - dt(i, j)
          ww2(i, j) = hp1 - du(i, j)
          emv(i, j) = ww1(i, j) * ww2(i, j) * em3v(it(i, j)-9) +
     &                ww2(i, j) * dt(i, j) * em3v(it(i, j)-8) +
     &                ww1(i, j) * du(i, j) * em3v(it(i, j)+19) +
     &                dt(i, j) * du(i, j) * em3v(it(i, j)+20)
203   continue
      return
      end 
