c $Log: trim.f,v $
c Revision 1.11  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.10  1998/12/10  00:55:53  ldr
c HBG changes to V5-1-21
c
c Revision 1.9  1997/12/17  23:22:57  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.8  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.7  93/12/17  15:34:22  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.6  93/08/19  15:11:09  ldr
c Minor cosmetic changes.
c 
c Revision 1.5  92/12/09  14:44:46  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.4  92/10/20  16:54:19  ldr
c Changes to V4-0 for SGI (mainly to get nsib and seaice stuff running.)
c 
c Revision 1.3  92/06/16  12:22:16  ldr
c Replaced working common /trimw with local arrays.
c 
c Revision 1.2  91/03/13  13:01:18  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:38:22  ldr
c Initial release V3-0
c 
      subroutine trim(rhs,a,c,u,it,e,temp,b)

C Global parameters
      include 'PARAMS.f'

C Argument list
      real rhs(ln2,nl)
      real a(ln2,nl)
      real c(ln2,nl)
      real u(ln2,nl)
      integer it
      real e(ln2,nl)
      real temp(ln2,nl)
      real b(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real g(ln2,nl)

C Local data, functions etc

C Start code : ----------------------------------------------------------

C     N.B. WE NOW ALWAYS ASSUME B = 1-A-C
C     RHS SHARED BETWEEN TRIM AND VERTMIX

C
C     THIS ROUTINE SOLVES THE SYSTEM
C       A(K)*U(K-1)+B(K)*U(K)+C(K)*U(K+1)=RHS(K)    FOR K=2,Nl-1
C       WITH  B(K)*U(K)+C(K)*U(K+1)=RHS(K)          FOR K=1
C       AND   A(K)*U(K-1)+B(K)*U(K)=RHS(K)          FOR K=Nl
C
C     THE THOMAS ALGORITHM IS USED
C
      if(it.ne.0)go to 13
      do 1002 k=1,nl
        do 1002 mg=1,ln2
          b(mg,k)=1-a(mg,k)-c(mg,k)
 1002 continue

      do 1004 mg=1,ln2
1004  e(mg,1)=c(mg,1)/b(mg,1)
      do 10 k=2,nl-1
      do 10 mg=1,ln2
      temp(mg,k)= 1/(b(mg,k)-a(mg,k)*e(mg,k-1))
10    e(mg,k)=c(mg,k)*temp(mg,k)
C     PRINT *,'TEMP ',TEMP
C     USE PRECOMPUTED VALUES WHEN AVAILABLE
C12   PRINT *,'E ',IT,E
C     PRINT *,'TEMP ',IT,TEMP
C
13    do 14 mg=1,ln2
14    g(mg,1)=rhs(mg,1)/b(mg,1)
      do 16 k=2,nl-1
      do 16 mg=1,ln2
16    g(mg,k)=(rhs(mg,k)-a(mg,k)*g(mg,k-1))*temp(mg,k)
      do 1602 mg=1,ln2
1602  u(mg,nl)=(rhs(mg,nl)-a(mg,nl)*g(mg,nl-1))/
     &        (b(mg,nl)-a(mg,nl)*e(mg,nl-1))
      do 20 kk=1,nl-1
      k=nl-kk
      do 20 mg=1,ln2
20    u(mg,k)=g(mg,k)-e(mg,k)*u(mg,k+1)
      return
      end
