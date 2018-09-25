c $log$
      subroutine trim3x(a,b,c,rhs,y,ns,npts)

C Global parameters
      include 'PARAMS.f'
      parameter (ms3=ms+3)

C Argument list
      real a(ln2,1)
      real b(ln2,1)
      real c(ln2,1)
      real rhs(ln2,1)
      real y(ln2,1)
      integer ns
      integer npts

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real e(ln2,ms3),g(ln2,ms3),temp(ln2,ms3)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c     this routine solves the system
c       a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,ns-1
c       with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
c       and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=ns

c     the Thomas algorithm is used

      do 400 mg=1,npts
  400   e(mg,1)=c(mg,1)/b(mg,1)
      do 12 k=2,ns-1
        do 410 mg=1,npts
        temp(mg,k)= 1./(b(mg,k)-a(mg,k)*e(mg,k-1))
  410   e(mg,k)=c(mg,k)*temp(mg,k)
   12 continue

      do 420 mg=1,npts
  420 g(mg,1)=rhs(mg,1)/b(mg,1)
      do 16 k=2,ns-1
        do 430 mg=1,npts
  430   g(mg,k)=(rhs(mg,k)-a(mg,k)*g(mg,k-1))*temp(mg,k)
   16 continue

c     do back substitution to give answer now
      do 440 mg=1,npts
      y(mg,ns)=(rhs(mg,ns)-a(mg,ns)*g(mg,ns-1))/
     .        (b(mg,ns)-a(mg,ns)*e(mg,ns-1))
  440 continue
      do 20 k=ns-1,1,-1
        do 450 mg=1,npts
  450   y(mg,k)=g(mg,k)-e(mg,k)*y(mg,k+1)
   20 continue

      return
      end
