c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: rainda.f,v $
c Revision 1.34  2001/02/22 05:34:41  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.33  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.32  1998/12/10  00:55:39  ldr
c HBG changes to V5-1-21
c
c Revision 1.31  1997/12/17  23:22:50  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.30  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.29  1996/10/24  01:03:13  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.28  1996/06/13  02:08:05  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.27  1996/03/25  04:31:29  ldr
c Replace clhx by hlcp in conv, rainda and radin (needed for V5-0-7).
c
c Revision 1.26  1996/03/21  03:19:02  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.25  1995/08/31  04:30:47  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.24  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.23  1995/06/30  02:44:43  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.22  1994/08/08  17:22:28  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.21  94/08/08  13:16:27  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.20  94/08/04  16:56:31  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.19  93/12/23  15:35:08  ldr
c Put establ back into an include file.
c 
c Revision 1.18  93/12/17  15:33:39  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.17  93/09/07  10:26:44  ldr
c Move dry adjustment from rainda to conv.
c 
c Revision 1.16  93/09/06  12:28:11  ldr
c Use correct formula for qs and dqsdt (i.e. don't subtract es from p).
c 
c Revision 1.15  93/08/19  15:09:09  ldr
c Minor cosmetic changes.
c 
c Revision 1.14  93/08/17  16:37:22  ldr
c Slight change to lock around the floating exception trap.
c 
c Revision 1.13  93/07/26  12:57:27  ldr
c Put locks around infinite do loop trap in SGI code.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  debug    - flag to control column debugging
c                  lgdebug  - flag to control latitude debugging
c                  mgdebug  - flag to control longitude debugging
c                  insdebug - flag to control hemisphere debugging
c
c              from common/hybrpr in HYBRPR.f
c                  dprf - pressure thickness at each sigma level
c                  prf  - pressure at full levels 
c
c              from common/levdata in this subroutine
c                  klowrn -  low rainfall (up to sig=0.9)
c                  kmidrn - middle rain (up to sig=0.45) 
c
c              from common/printt in PRINTT.f
c                  cvrnm - convection/rain mapping indicator
c
c              from common/timex in TIMEX.f
c                  mins  - current model time in mins
c                  mstep - timestep in minutes 
c
c     Output:  from common/cloud1 in this subroutine
c                  rainp - rain map 
c
c     In/Out:  from arguments
c                  gam - (L/cp)dqsdt 
c                  precs - array for accumulated large scale precipitation
c                  rhg - relative humidity      qsg - ice value
c                  ttg - physical temperature   qtg - mixing ratio
c 
      subroutine rainda (lg,    !Inputs
     &                   ttg,qtg,precs, !In and Out
     &                   qsg,rhg,gam)           !Outputs                  

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real precs(ln2)
      real qsg(ln2,nl)
      real rhg(ln2,nl)
      real gam(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'PRINTT.f'
      include 'TIMEX.f'

      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

C Local work arrays and variables
      integer irain(ln2,nl)
      real pk(ln2)
      logical reevap

      integer irn
      integer jh
      integer k
      integer mg
      integer ns

      real dqrx
      real hlrvap
      real rex
      real rhmax

C Local data, functions etc
      character*1 cont(8)
      data cont/' ','L','M','4','H','5','6','7'/
      data rhmax/1.00/
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

      hlrvap=hl/rvap
c**** check for rainfall
      do 10 k=1,nl
         do 10 mg=1,ln2
 10      irain(mg,k)=0
c---- Either rain from all levels, no evaporating
c---- Or rain from top down, re-evaporating
      reevap=.true.
c      reevap=.false.
      If(.not.reevap)Then

      do 15 k=1,nl
      do mg=1,ln2
        pk(mg)=100.0*prf(mg,k)
      enddo
      if(ukconv)then
        call QSATU(qsg(1,k),ttg(1,k),pk,ln2)
      else
        do mg=1,ln2
           qsg(mg,k)=qsat(pk(mg),ttg(mg,k))
        enddo
      endif
      do 152 mg=1,ln2
      rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
      gam(mg,k)=hlcp*qsg(mg,k)*hlrvap/ttg(mg,k)**2
      if(rhg(mg,k).gt.rhmax)then
cx      cloud(mg,k)=1.0
        irain(mg,k)=1
        dqrx=(qtg(mg,k)-rhmax*qsg(mg,k))/(1.0+rhmax
     &    *gam(mg,k))
        ttg(mg,k)=ttg(mg,k)+hlcp*dqrx
        qtg(mg,k)=qtg(mg,k)-dqrx
        qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dqrx
        rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
        precs(mg)=precs(mg)+0.5*dqrx*dprf(mg,k)*100.0/grav
      end if
  152 continue
   15 continue

      else

      k=nl
      do mg=1,ln2
        pk(mg)=100.0*prf(mg,k)
      enddo
      if(ukconv)then
        call QSATU(qsg(1,k),ttg(1,k),pk,ln2)
      else
        do mg=1,ln2
           qsg(mg,k)=qsat(pk(mg),ttg(mg,k))
        enddo
      endif
      do 20 mg=1,ln2
         rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
         gam(mg,k)=hlcp*qsg(mg,k)*hlrvap/ttg(mg,k)**2
         if (rhg(mg,k).gt.rhmax) then
c       cloud(mg,k)=1.0
            irain(mg,k)=1
            dqrx=(qtg(mg,k)-rhmax*qsg(mg,k))/(1.0+rhmax*gam(mg,k))
            ttg(mg,k)=ttg(mg,k)+hlcp*dqrx
            qtg(mg,k)=qtg(mg,k)-dqrx
            rex=dqrx*dprf(mg,k)/dprf(mg,k-1)
            qtg(mg,k-1)=qtg(mg,k-1)+rex
            ttg(mg,k-1)=ttg(mg,k-1)-hlcp*rex
         end if
 20   continue
      do 30 k=nl-1,2,-1
      do mg=1,ln2
        pk(mg)=100.0*prf(mg,k)
      enddo
      if(ukconv)then
        call QSATU(qsg(1,k),ttg(1,k),pk,ln2)
      else
        do mg=1,ln2
           qsg(mg,k)=qsat(pk(mg),ttg(mg,k))
        enddo
      endif
         do 30 mg=1,ln2
         rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
         gam(mg,k)=hlcp*qsg(mg,k)*hlrvap/ttg(mg,k)**2
         if (rhg(mg,k).gt.rhmax) then
c       cloud(mg,k)=1.0
            irain(mg,k)=1
            dqrx=(qtg(mg,k)-rhmax*qsg(mg,k))/(1.0+rhmax*gam(mg,k))
            ttg(mg,k)=ttg(mg,k)+hlcp*dqrx
            qtg(mg,k)=qtg(mg,k)-dqrx
            rex=dqrx*dprf(mg,k)/dprf(mg,k-1)
            qtg(mg,k-1)=qtg(mg,k-1)+rex
            ttg(mg,k-1)=ttg(mg,k-1)-hlcp*rex
            qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dqrx
            rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
         end if
 30   continue
c---- source of large scale rainfall (precs) is now only
c---- from level 1
      k=1
      do mg=1,ln2
        pk(mg)=100.0*prf(mg,k)
      enddo
      if(ukconv)then
        call QSATU(qsg(1,k),ttg(1,k),pk,ln2)
      else
        do mg=1,ln2
           qsg(mg,k)=qsat(pk(mg),ttg(mg,k))
        enddo
      endif
      do 40 mg=1,ln2
         rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
         gam(mg,k)=hlcp*qsg(mg,k)*hlrvap/ttg(mg,k)**2
         if (rhg(mg,k).gt.rhmax) then
c       cloud(mg,k)=1.0
            irain(mg,k)=1
            dqrx=(qtg(mg,k)-rhmax*qsg(mg,k))/(1.0+rhmax*gam(mg,k))
            ttg(mg,k)=ttg(mg,k)+hlcp*dqrx
            qtg(mg,k)=qtg(mg,k)-dqrx
            qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dqrx
            rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
            precs(mg)=precs(mg)+0.5*dqrx*dprf(mg,k)*100.0/grav
         end if
 40   continue

      Endif

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,i1)')'After rainda.'
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,99)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,99)'qsg ',(qsg(mg,k),k=1,nl)
          write(25,91)'rhg ',(rhg(mg,k),k=1,nl)
          write(25,1)'precs ',precs(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,g10.3))
 91   format(a,18f10.3)
 99   format(a,18g10.3)

      if ((mod(mins+int(mstep, 8),1440_8).eq.0_8).and.cvrnm) then
c**** update the rainfall map at end of day only
c**** Low rainfall (up to sig=0.9), Middle Rain (up to sig=0.45), High above
c**** (see initax.f for klowrn,kmidrn)
         do 50 mg=1,ln2
            irn=1
            jh=0
            do 52 k=kmidrn+1,nl
   52       jh=jh+irain(mg,k)
            if (jh.ge.1) irn=5
            jh=0
            do 54 k=klowrn+1,kmidrn
   54       jh=jh+irain(mg,k)
            if (jh.ge.1) irn=irn+2
            jh=0
            do 56 k=1,klowrn
   56       jh=jh+irain(mg,k)
            if (jh.ge.1) irn=irn+1
 50         rainp(mg,lg)=cont(irn)
      end if


      return
      end
