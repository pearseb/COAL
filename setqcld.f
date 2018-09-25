c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: setqcld.f,v $
c Revision 1.2  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine setqcld(lg,index,qlg,qfg,ccov,gam,pgam)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /RVARSG/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      integer index
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real ccov(ln2,nl)
      real gam(ln2,nl)
      real pgam(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'RVARSG.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'QCLOUD1.f'

C Local work arrays and variables
      integer k
      integer mg

      real dqsdt
      real pk
      real qs

C Local data, functions etc
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

      if(index.eq.3)go to 100

c Copy cloud model variables from big arrays into working arrays
      if(qcloud)then
        do k=1,nl
          do mg=1,ln2
            qlg(mg,k)=qlb(mg,k,lg)
            qfg(mg,k)=qfb(mg,k,lg)
            ccov(mg,k)=cfb(mg,k,lg) !Estimate initial ccov from t-1 value
          enddo
        enddo
        if(index.eq.1)then
          do k=1,nl
          do mg=1,ln2
c Work out gam=(L/cp)*dqsdt which is needed for hsflux, hvertmx
            pk=100.0*prf(mg,k)
            qs=qsat(pk,ttg(mg,k))
            dqsdt=qs*hl/(rvap*ttg(mg,k)**2)
            gam(mg,k)=(hl/cp)*dqsdt
            pgam(mg,k)=gam(mg,k)  !Set leads value now
          enddo
          enddo
        endif
      else
        do k=1,nl
          do mg=1,ln2
            qlg(mg,k)=0.
            qfg(mg,k)=0.
          enddo
        enddo
      endif
 
      return

  100 continue

c Copy cloud model variables back into big arrays

      do k=1,nl
        do mg=1,ln2
          qlb(mg,k,lg)=qlg(mg,k)
          qfb(mg,k,lg)=qfg(mg,k)
          cfb(mg,k,lg)=ccov(mg,k)
        enddo
      enddo

      return
      end
