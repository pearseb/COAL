c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Form two different versions of one of the parallel regions, one for qcloud=T
c and one for qcloud=F. This enables us to retain the previous loop fusion,
c while ensuring correct results when qcloud=F.
c SJP 2004/05/25
c
c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c REDUCTION operations commented out, for reproducibility.
c SJP 2004/01/05
c
c Manual parallelisation of JMCGSLT - but only those loops that are executed
c with the way that we configure the model. Loops that are only executed if
c SLTRACE, COUPLED_AERO or DEBUG are T are not parallelised. Outer two loops
c in second parallelised loop swapped, so that the outermost loop is now of
c size LAT, and some prognostic cloud scheme work incorporated into one of the
c other loops for improved performance.
c SJP 2003/04/25
c
c Manual parallelisation of ENFORCE_CONQ - note that as we use a REDUCTION
c operation on DELQPOS and DELQNEG, the results will vary slightly from run to
c run on multiple CPUs. The results will not necessarily be less accurate than
c those obtained on one CPU, however.
c SJP 2003/04/02
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/03/29
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: jmcgslt.f,v $
c Revision 1.46  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.45  2001/06/04 02:26:57  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.44  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.43  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.42  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.41  2000/06/20 02:08:33  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.40.1.1  2000/06/20 02:27:10  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.40  1999/06/16 06:21:54  rot032
c HBG changes to V5-3
c
c Revision 1.39  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.38  1998/12/10  00:55:46  ldr
c HBG changes to V5-1-21
c
c Revision 1.37  1998/01/30  05:05:37  ldr
c Further parallelization stuff for NEC from MRD.
c
c Revision 1.36  1996/10/24  01:02:59  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.35  1996/06/13  02:07:01  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.34  1996/03/21  03:18:53  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.33  1996/02/19  23:23:07  ldr
c Do vertical SLT properly for 24 level version.
c
c Revision 1.32  1995/11/23  05:53:46  ldr
c Cleaner SGI implementation.
c
c Revision 1.31  1995/08/29  03:46:58  ldr
c Added "shared" directive for muf for kaos.
c
c Revision 1.30  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.29  1995/08/08  02:02:16  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.28  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.23.2.2  1995/07/04  06:23:13  ldr
c Corrections to V4-5-30mic to make it work on cherax.
c
c Revision 1.23.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.27  1994/08/08  17:21:39  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.26  94/08/05  16:43:37  ldr
c Tidy up IGW's tracer conservation code so that it works on SGI.
c 
c Revision 1.25  94/08/04  17:07:01  ldr
c Merge of LDR and HBG changes.
c 
c Revision 1.24  94/07/11  17:03:36  ldr
c Go back to old clumsy method for delqpos, delqneg, so that results do
c not vary with no. of processors on SGI.
c 
c Revision 1.23.1.1  94/08/04  16:55:50  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.23  94/05/13  14:55:52  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c 
c Revision 1.22  94/04/28  17:20:41  ldr
c Remove JMcG's max/min based code which is slow on Seca, and tidy up
c calculation of delqpos, delqneg by declaring them REDUCTION.
c 
c Revision 1.21  94/03/22  15:37:07  ldr
c Make 2 pass conservation fix the standard and do no conservation on qc.
c 
c Revision 1.20  93/10/05  13:06:33  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c     INPUT/OUTPUT
c     Input:   from common/cnsta in CNSTA.f
c                  dsk - sigma thicknesses (1 to nl)
c
c              from common/fewflags in FEWFLAGS.f
c                  qcloud - T if using prognostic cloud scheme
c                  sltrace - if T compute tracer transport by semi-
c                            Langrangian advection and by mixing and
c                            convection in the vertical
c
c              from common/gausl in GAUSL.f
c                  w - gaussian weights*cos(lat)
c
c              from common/hybarr in HYBARR.f
c                  dadnf,dbdnf - hybrid vertical coordinate arrays, full level
c
c              from common/timex in TIMEX.f
c                  mins  - current model time in mins
c                  mstep - timestep in minutes 
c
c              from common/uvpgd in UVPGD.f
c                  dpsdt - omega/sigma at full levels (formerly dpstar/dt)
c
c     Output:  from common/rmgrid in RMGRID.f
c                  rgt - save moisture tendency for possible use in Kuo scheme
c
c     In/Out:  from common/cnsta in CNSTA.f
c                  sig - sigma values for levels 1 to nl
c
c              from common/fewflags in FEWFLAGS.f
c                  glmean_interval - interval in mins between prints of global
c                                    means to standard output
c
c              from common/qcloud1 in QCLOUD1.f
c                  qfb - cloud ice at t  qfbm - cloud ice at t-1
c
c              from common/rmgrid in RMGRID.f
c                  rmg   - pressure weighted moistures, current timestep
c                  rmmg  -    "        "        "     , previous timestep
c
c              from common/traceblk in TRACEBLK.f
c                  con - concentration of atmospheric tracers at t
c                  conm - concentration of atmospheric tracers at t-1 
c                  conmnth - monthly mean tracer concentration, tracer no. nt
c
c              from arguments
c                  exptyp - restart file header
c                  tdt -  2*timestep (for 3 time level scheme)
c                  tot_water - global integral of moisture
c 
c 
      subroutine jmcgslt(tdt,tot_water,exptyp)

c Interface to JMcG's Semi-Lagrangian moisture advection scheme (LDR 7/92)
c Asselin filter is applied to the old value (RMMG) after the time-step.
c In this routine, the grid runs from 1 to lat2 and 1 to lon, so that
c e.g index -1 in this routine corresponds to 1 in horadv.

c INPUTS:
c tdt   - 2*timestep (for 3 time level scheme)
c sig   - Sigma values for levels 1 to nl
c dsk   - Sigma thicknesses (1 to nl).
c rmg   - Pressure weighted moistures at current timestep
c rmmg  -    "        "        "      "  previous timestep
c sdot  - sdot at full levels 1 to nl
c ureal - Zonal wind in m/s
c vreal - Meridional wind in m/s
c pstar - Surface pressure
c dpsdt - omega/sigma at full levels (formerly dpstar/dt)
c w     - Gaussian weights

c OUTPUTS:
c rmg   - New moisture field after timestep
c rmmg  - New value of moisture at t-1 - has had Asselin filter applied.
c rgt   - Save moisture tendency for possible use in Kuo scheme.
c tot_water - Global integral of moisture      

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f'
      integer il,jl,kl
      parameter(il=lon,jl=lat2,kl=nl)
      integer il4,jl4
      parameter(il4=il+4,jl4=jl+4)
      real asf,asfx
      parameter(asf=0.05, asfx=1.0-2*asf) !Asselin filter coefficients

C Argument list
      real tdt
      real tot_water
      character*50 exptyp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'    !Input sig, sigh, dsk
      include 'ECTRAC.f'   !Input and output bigxt, bigxtm
      include 'FEWFLAGS.f' !Input sltrace
      include 'GAUSL.f'    !Input w
      include 'HYBARR.f'
      include 'QCLOUD1.f'  !Input and output qlb, qlbm, qfb, qfbm
      include 'RMGRID.f'   !Input and output rmg, rmmg. Output rgt.
      include 'TIMEX.f'    !Input mins, mstep
      include 'TRACEBLK.f' !Input and output con, conm
      include 'UVPGD.f'    !Input ureal,vreal,pstar,dpsdt,sdot

C Local work arrays and variables
      real qq(-1:lon+2,-1:lat2+2,nl),delq(ln2,lat,nl,3)
      real ql(-1:lon+2,-1:lat2+2,nl),qf(-1:lon+2,-1:lat2+2,nl)
      real muf(lon,nl,lat2),wdsk(lon,nl,lat2)
      real st(il,jl,kl)
      real conmsav(-1:lonx,-1:lat2x,nl,ntrace) !NX is memory saving trick
c      real tottrace(ntrace),altrace(ntrace)
      real xtm(-1:lone+2,-1:lat2e+2,nl,ntrac),
     &     xt(-1:lone+2,-1:lat2e+2,nl,ntrac)
      real xtmsav(-1:lone+2,-1:lat2e+2,nl,ntrac)
      real delx(lone,lat2e,nl,ntrac)
      integer kdel(il,jl,kl)

      logical aspline ! For Akima spline
      logical ncepx

      integer i
      integer j
      integer k
      integer kq
      integer kql
      integer lg
      integer lgns
      integer ma
      integer mg
      integer nfield
      integer ns
      integer nt
      integer ntrace1

      real alpha
      real alphax
      real anhk
      real bnhk
      real contemp
      real delt
      real deltneg
      real deltnegs
      real deltpos
      real deltposs
      real dp
      real one_on_alpha
      real rhf
      real rhodz
      real totq
      real tot_trace
      real wdskmuf
      real xint

C Local data, functions etc
      logical firstcall,secondcall
      save firstcall,secondcall
      data firstcall,secondcall / 2*.true. /
      integer ncons
      data ncons/2/   ! for double-pass conservation

C Start code : ----------------------------------------------------------

C**** COMPUTE layer pressure WEIGHTS (if .not.hybrid then muf=pstar)

!$OMP  PARALLEL DO SCHEDULE (DYNAMIC)
!$OMP& PRIVATE (anhk, bnhk, k, lg, lgns, mg, ns)
!$OMP& SHARED  (anh, bnh, dadnf, dbdnf, dsk, hybrid, muf, pstar, w,
!$OMP&          wdsk)

      do lgns=1,lat2

        ns=2-(lgns-1)/lat
        lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)

c---- compute muf  = dp/d(coord) for hybrid coords (=P* if sigma)
c---- compute wdsk = w(lg)*dsk(k) if not hybrid (w=Gaussian weight)
c----              = w(lg)*d(p)/P* if hybrid :
c---- half level pressures = anh(k)+bnh(k)*pstar(mg,lgns)
c---- pressure thickness d(p) at level k = prh(mg,k)-prh(mg,k+1)
        if(hybrid)then
          do k=1,nl
            anhk=(anh(k)-anh(k+1))
            bnhk=(bnh(k)-bnh(k+1))
            do mg=1,lon
              muf(mg,k,lgns)=dadnf(k)+dbdnf(k)*pstar(mg,lgns)
              wdsk(mg,k,lgns)=w(lg)*(anhk/pstar(mg,lgns)+bnhk)
            enddo
          enddo
        else
          do k=1,nl
            do mg=1,lon
              muf(mg,k,lgns)=pstar(mg,lgns)
              wdsk(mg,k,lgns)=w(lg)*dsk(k)
            enddo
          enddo
        endif

      enddo ! lgns=1,lat2

!$OMP END PARALLEL DO

c     Diagnostic of total water:

      tot_water=0.

      if(mod(mins+int(mstep, 8),int(glmean_interval, 8)).eq.0_8)then

!!$OMP  PARALLEL DO SCHEDULE (DYNAMIC)
!!$OMP& PRIVATE   (k, lg, lgns, ma, mg, ns)
!!$OMP& REDUCTION (+ : tot_water)
!!$OMP& SHARED    (muf, rmg, wdsk)

        do lg=1,lat
          do ns=1,2
            lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
            do k=1,nl
              do mg=1,lon
                ma=mg+(ns-1)*lon
                tot_water=tot_water+rmg(ma,k,lg)*wdsk(mg,k,lgns)
     &                              *muf(mg,k,lgns)
              enddo
            enddo
          enddo
        enddo

!!$OMP END PARALLEL DO

        tot_water=0.5*tot_water/float(lon)*100./grav
      endif

c Tracers (if any)
      if(sltrace) then
        ntrace1=2
c        tottrace(1)=0.
c        altrace(1)=0.
        do nt=ntrace1,ntrace
          tot_trace=0.
          do lgns=1,lat2
            do k=1,nl
              do mg=1,lon
                tot_trace=tot_trace+
     &          con(mg,lgns,k,nt)*wdsk(mg,k,lgns)*muf(mg,k,lgns)
                conmnth(mg,lgns,k,nt)=conmnth(mg,lgns,k,nt)
     &                                        +con(mg,lgns,k,nt)
              enddo
            enddo
          enddo
          tot_trace=tot_trace*0.5/float(lon)*100./grav 
c          tottrace(nt)=tot_trace
        enddo
      endif


c     Set up moisture array qq (runs from south to north pole)

      if (qcloud) then

!$OMP  PARALLEL DO SCHEDULE (DYNAMIC)
!$OMP& PRIVATE (k, lg, mg)
!$OMP& SHARED  (qf, qfbm, ql, qlbm, qq, rmmg)

      do lg=1,lat
        do k=1,nl
          do mg=1,lon
            qq(mg,lg,k)=rmmg(mg+lon,k,lg)
            qq(mg,lat2p-lg,k)=rmmg(mg,k,lg)
            ql(mg,lg,k)=qlbm(mg+lon,k,lg)
            ql(mg,lat2p-lg,k)=qlbm(mg,k,lg)
            qf(mg,lg,k)=qfbm(mg+lon,k,lg)
            qf(mg,lat2p-lg,k)=qfbm(mg,k,lg)

          enddo
        enddo
      enddo

!$OMP END PARALLEL DO

      else

!$OMP  PARALLEL DO SCHEDULE (DYNAMIC)
!$OMP& PRIVATE (k, lg, mg)
!$OMP& SHARED  (qq, rmmg)

      do lg=1,lat
        do k=1,nl
          do mg=1,lon
            qq(mg,lg,k)=rmmg(mg+lon,k,lg)
            qq(mg,lat2p-lg,k)=rmmg(mg,k,lg)
          enddo
        enddo
      enddo

!$OMP END PARALLEL DO

      end if

c Set min top level moisture to help with strat SLT boundary condition

!$OMP  PARALLEL DO SCHEDULE (DYNAMIC)
!$OMP& PRIVATE (lg, mg)
!$OMP& SHARED  (qq)

      do lg=1,lat2
        do mg=1,lon
          qq(mg,lg,nl)=max(qq(mg,lg,nl),0.1e-05)
        enddo
      enddo

!$OMP END PARALLEL DO

      if(exptyp(17:20).ne.'SLAG')then
        if(firstcall.or.secondcall)then
          print*,'Warning: truncating top level q'
          if(.not.firstcall)then
            secondcall=.false.
            exptyp(17:20)='SLAG'
          endif
          firstcall=.false.
          do 385 j=1,jl
          do 385 i=1,il
 385      qq(i,j,kl)=min(qq(i,j,kl), 4.0 e-8) !Truncate top level q
        endif
      endif

c Vertical advection:
c First calculate departure point.

c If using Akima spline in vertical advection, set aspline to true
      aspline=.true.
c      aspline=.false.

      if(aspline)then
        call vadvect3(tdt,sdot,                !Inputs
     &                st,kdel)                 !Output (departure points)
      else
        call vadvect(tdt,sdot,                 !Inputs
     &                st)                      !Output (departure points)
      endif

c Interpolate to find field value qq at departure point st

      if(kl.gt.18)then
        nfield=2       !Don't use sigma^3 treatment for 24 level model
      else
        nfield=1
      endif

      if(aspline)then
        call vinterp3(st,kdel,nfield, !Input
     &                qq)             !In and out
      else
        call vinterp(st,nfield,       !Input
     &               qq)              !In and out
      endif
      
c Vertical tracer advection: interpolate to find field value con(;,;,;,nt)
c at departure point st.

      if ( sltrace ) then
        do nt=ntrace1,ntrace           !loop over no. of tracers
          do k=1,nl
            do lgns=1,lat2
              do mg=1,lon
                conmsav(mg,lgns,k,nt)=conm(mg,lgns,k,nt) !Save for Asselin filt
              enddo
            enddo
          enddo
          nfield=2
          if(aspline)then
            call vinterp3(st,kdel,nfield,   !Input
     &                    conm(-1,-1,1,nt)) !In and out
          else
            call vinterp(st,nfield,         !Input
     &                   conm(-1,-1,1,nt))  !In and out
          endif
        enddo
      endif

c Vertical advection of ECHAM tracers

      if ( coupled_aero ) then

        do nt=1,ntrac           !loop over no. of tracers

c Reorder ECHAM arrays for SLT.

          do k=1,nl
            do lg=1,lat
              do mg=1,lon
                xtm(mg,lat2p-lg,k,nt)=bigxtm(mg,k,nt,lg) !NH
                xtm(mg,lg,k,nt)=bigxtm(mg+lon,k,nt,lg)   !SH
                xtmsav(mg,lat2p-lg,k,nt)=bigxtm(mg,k,nt,lg) !NH
                xtmsav(mg,lg,k,nt)=bigxtm(mg+lon,k,nt,lg)   !SH
                xt(mg,lat2p-lg,k,nt)=bigxt(mg,k,nt,lg) !NH
                xt(mg,lg,k,nt)=bigxt(mg+lon,k,nt,lg)   !SH
              enddo
            enddo
          enddo

          if(debug)then
            lgns=lgdebug*(insdebug-1)+(lat2p-lgdebug)*(2-insdebug)
            mg=mgdebug
            xint=0.
            do k=1,nl
              dp=100.*wdsk(mg,k,lgns)*pstar(mg,lgns)/w(lgdebug)
              rhodz=dp/grav
              xint=xint+(xtm(mg,lgns,k,nt))*rhodz
            enddo
            write(25,*)
            write(25,*)'Before V advection, nt = ',nt
            write(25,1)' xint ',xint
            write(25,9)'xtm ',(xtm(mg,lgns,k,nt),k=1,nl)
            write(25,9)'sdot ',(sdot(mg,k,lgns),k=1,nl)

            totq=0.
            do ns=1,2
              do lg=1,lat
                lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
                do k=1,nl
                  do mg=1,lon
                   ma=mg+(ns-1)*lon
                   totq=totq+(rmmg(ma,k,lg)+qlbm(ma,k,lg)+qfbm(ma,k,lg))
     &                   *wdsk(mg,k,lgns)*muf(mg,k,lgns)
                  enddo
                enddo
              enddo
            enddo
            totq=0.5*totq/float(lon)*100./grav
            write(25,*)'Before adv, totq = ',totq

          endif

 1        format(3(a,g12.5))
 9        format(a,30g10.3)


c Use the TVD scheme for vertical tracer advection...

          call vadvtvd (tdt,sdot,        !Inputs
     &                  xtm(-1,-1,1,nt)) !In & out
          
          if(debug)then
            lgns=lgdebug*(insdebug-1)+(lat2p-lgdebug)*(2-insdebug)
            mg=mgdebug
            xint=0.
            do k=1,nl
              dp=100.*wdsk(mg,k,lgns)*pstar(mg,lgns)/w(lgdebug)
              rhodz=dp/grav
              xint=xint+(xtm(mg,lgns,k,nt))*rhodz
            enddo
            write(25,*)'After V advection, nt = ',nt
            write(25,9)'kdel ',(float(kdel(mg,lgns,k)),k=1,nl)
            write(25,9)'st ',(st(mg,lgns,k),k=1,nl)
            write(25,1)' xint ',xint
            write(25,9)'xtm ',(xtm(mg,lgns,k,nt),k=1,nl)
          endif
        enddo
      endif

c Vertical advection of cloud water and cloud ice

      if ( qcloud ) then
        nfield=3
        if(aspline)then
          call vinterp3(st,kdel,3,
     &                  ql)
          call vinterp3(st,kdel,3,
     &                  qf)
        else
          call vinterp(st,3,
     &                 ql)
          call vinterp(st,3,
     &                 qf)
        endif
      endif

      if(qcloud)then
        kql=3 ! number of qcloud water variables
      else
        kql=1
      endif
      if(coupled_aero)kql=kql+ntrac

C***      if (debug ) then
C***        nt=2
C***        k=2
C***        lgns=lgdebug*(insdebug-1)+(lat2p-lgdebug)*(2-insdebug)          
C***        mg=mgdebug
C***        write(26,*)'lgns, mg ',lgns,mg
C***        write(26,*)'Before hadv, nt = ',nt, ' xtm 3x3 '
C*** 93     format(3g13.4)
C***        write(26,93)xtm(mg-1,lgns+1,k,nt),xtm(mg,lgns+1,k,nt),
C***     &                     xtm(mg+1,lgns+1,k,nt)
C***        write(26,93)xtm(mg-1,lgns,k,nt),xtm(mg,lgns,k,nt),
C***     &                     xtm(mg+1,lgns,k,nt)
C***        write(26,93)xtm(mg-1,lgns-1,k,nt),xtm(mg,lgns-1,k,nt),
C***     &                     xtm(mg+1,lgns-1,k,nt)
C***        write(26,1)' ureal ',ureal(mg,lgns,k),' vreal ',vreal(mg,lgns,k)
C***      endif

      
c For SGI, unscoped variables default to "share", while on CRAY, all
c must be explicitly scoped "private" or "shared".

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (k)
!$OMP& SHARED  (tdt,qq,ql,qf,delq,muf,wdsk,aspline,xtm,xtmsav,delx)

      DO k=1,nl ! Multi-tasked over k levels

c Horizontal advection.
        call jmcghor(k,tdt,qq,ql,qf,delq,xtm,xtmsav,delx,
     &           muf,aspline)

      ENDDO ! Multi-tasked over k levels

!$OMP END PARALLEL DO

      if (debug ) then
        do nt=1,ntrac
          lgns=lgdebug*(insdebug-1)+(lat2p-lgdebug)*(2-insdebug)
          mg=mgdebug
          xint=0.
          do k=1,nl
            dp=100.*wdsk(mg,k,lgns)*pstar(mg,lgns)/w(lgdebug)
            rhodz=dp/grav
            xint=xint+(xtm(mg,lgns,k,nt))*rhodz
          enddo
          write(25,*)'After hadv, nt = ',nt
          write(25,1)' xint ',xint
          write(25,9)'xtm ',(xtm(mg,lgns,k,nt),k=1,nl)
        enddo
      endif

c Enforce conservation for water vapour and cloud water variables.
c Then do time-step, applying Asselin filter.
c Save the vapour tendency (rgt) for possible use in Kuo convection.
c NB: Not saving rgt any more, since Kuo is no longer used (LDR 4/2000)
c Time-step and filter cloud water and cloud ice if qcloud (kq=2,3)
c Remove extreme small traces left after SLT (set values < 1.0e-20 to zero)

*odir cncall
*odir concur
      do kq=1,kql
        if(kq.eq.1)then
          call enforce_conq(2,muf,wdsk,delq(1,1,1,kq),      !In
     &                      rmmg,rmg)                       !In and out
        elseif(kq.eq.2)then
          call enforce_conq(2,muf,wdsk,delq(1,1,1,kq),      !In
     &                      qlbm,qlb)                       !In and out
        elseif(kq.eq.3)then
          call enforce_conq(2,muf,wdsk,delq(1,1,1,kq),      !In
     &                      qfbm,qfb)                       !In and out
        elseif(kq.ge.4)then
          nt=kq-3
          call enforce_cont(2,muf,wdsk,delx(1,1,1,nt),        !In
     &                      xtmsav(-1,-1,1,nt),xt(-1,-1,1,nt))!In and out

        endif
      enddo

c Calculate global moisture after advection and fixup

      if(debug)then
        totq=0.
        do ns=1,2
          do lg=1,lat
            lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
            do k=1,nl
              do mg=1,lon
                ma=mg+(ns-1)*lon
                totq=totq+(rmg(ma,k,lg)+qlb(ma,k,lg)+qfb(ma,k,lg))
     &               *wdsk(mg,k,lgns)*muf(mg,k,lgns)
              enddo
            enddo
          enddo
        enddo
        totq=0.5*totq/float(lon)*100./grav
        write(25,*)'After fixup, totq = ',totq
      endif

c Reorganise arrays for aerosol tracers

      if(coupled_aero)then
        do nt=1,ntrac

          if(debug)then
            lgns=lgdebug*(insdebug-1)+(lat2p-lgdebug)*(2-insdebug)
            mg=mgdebug
            xint=0.
            do k=1,nl
              dp=100.*wdsk(mg,k,lgns)*pstar(mg,lgns)/w(lgdebug)
              rhodz=dp/grav
              xint=xint+(xt(mg,lgns,k,nt))*rhodz
            enddo
            write(25,*)'After fixup, nt = ',nt
            write(25,1)' xint ',xint
            write(25,9)'xt ',(xt(mg,lgns,k,nt),k=1,nl)
          endif

          do k=1,nl
            do lg=1,lat
              do mg=1,lon
                bigxt(mg,k,nt,lg)=xt(mg,lat2p-lg,k,nt) !NH
                bigxt(mg+lon,k,nt,lg)=xt(mg,lg,k,nt)   !SH
                bigxtm(mg,k,nt,lg)=xtmsav(mg,lat2p-lg,k,nt) !NH
                bigxtm(mg+lon,k,nt,lg)=xtmsav(mg,lg,k,nt)   !SH
              enddo
            enddo
          enddo
        
        enddo
      endif


c Tracers (if any)
      if ( sltrace ) then

c Time-step and filter tracer arrays con, conm
c add for conservation igw 8-8-93

        do nt=ntrace1,ntrace
          deltpos=0.
          deltneg=0.
          do k=1,nl
            do lgns=1,lat2
              deltposs=0.
              deltnegs=0.
              do mg=1,lon
                delt=conm(mg,lgns,k,nt)-conmsav(mg,lgns,k,nt)
                rhf=max(0.0,delt)
                wdskmuf=muf(mg,k,lgns)*wdsk(mg,k,lgns)
                deltposs=deltposs+      rhf *wdskmuf
                deltnegs=deltnegs+(delt-rhf)*wdskmuf
              enddo
              deltpos=deltpos+deltposs
              deltneg=deltneg+deltnegs
            enddo
          enddo
          alphax=1.
          if(deltpos.eq.0.) go to 998
          alphax = -deltneg/deltpos
c          altrace(nt)=alphax
          
          if(alphax.gt.1)then   !2-sided alpha if > 1
            alpha = sqrt (-deltneg/deltpos)
            one_on_alpha = 1./(alpha)
          else
            alpha = -deltneg/deltpos
            one_on_alpha = 1.
          endif
          
          do k=1,nl
            do lgns=1,lat2
              do mg=1,lon
                delt=conm(mg,lgns,k,nt)-conmsav(mg,lgns,k,nt) 
                if(delt.gt.0.) then
                  delt=alpha*delt
                else
                  delt=one_on_alpha*delt
                endif
                conm(mg,lgns,k,nt)=conmsav(mg,lgns,k,nt)+delt 
              enddo
            enddo
          enddo
          
 998      continue
        enddo
c     end conservation method
        
        ncepx=(ncepstrt.eq.2).or.(ncepstrt.eq.0)
        do nt=ntrace1,ntrace
          do k=1,nl
            do lgns=1,lat2
              do mg=1,lon
                contemp=conm(mg,lgns,k,nt) !t+1 value
                if(.not.ncepx)
     &          conm(mg,lgns,k,nt)=asfx*con(mg,lgns,k,nt)+
     &               asf*conmsav(mg,lgns,k,nt)
                con(mg,lgns,k,nt)=contemp
              enddo
            enddo
          enddo
        enddo
        
c        write(89,'(i10,3e13.6,3f10.3)') mins+mstep,(tottrace(jj),jj=1,3)
c     &       ,(altrace(jj),jj=1,3)
        
      endif

      return
      end
c******************************************************************************
      subroutine jmcghor(k,tdt,qq,ql,qf,delq,xtm,xtmsav,delx,
     &           muf,aspline)

      implicit none
C Global parameters
      include 'PARAMS.f'
c     include 'PHYSPARAMS.f'
      include 'ECPARM.f'
      integer il,jl,kl
      parameter(il=lon,jl=lat2,kl=nl)
      integer il4,jl4
      parameter(il4=il+4,jl4=jl+4)

C Argument list
      integer k ! level
      real tdt
      real qq(-1:lon+2,-1:lat2+2,nl)
      real ql(-1:lon+2,-1:lat2+2,nl)
      real qf(-1:lon+2,-1:lat2+2,nl)
      real delq(ln2,lat,nl,3)
      real xtm(-1:lone+2,-1:lat2e+2,nl,ntrac)
      real xtmsav(-1:lone+2,-1:lat2e+2,nl,ntrac)
      real delx(lone,lat2e,nl,ntrac)
      real muf(lon,nl,lat2)
      logical aspline ! For Akima spline

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'    !Input sig
      include 'FEWFLAGS.f' !Input sltrace, coupled_aero
      include 'QCLOUD1.f'  !Input and output qlb, qlbm, qfb, qfbm
      include 'RMGRID.f'   !Input and output rmg, rmmg. Output rgt.
      include 'TRACEBLK.f' !Input and output con, conm
      include 'UVPGD.f'    !Input ureal,vreal,pstar,dpsdt,sdot

C Local work arrays and variables
      real xg(il4,jl4),yg(il4,jl4)
      integer idel(il4,jl4),jdel(il4,jl4)

      integer lg
      integer lgns
      integer ma
      integer mg
      integer ns
      integer nt
      integer ntrace1

      real bigq
      real qqx

C Local data, functions etc

C Start code : ----------------------------------------------------------

        if(.not.aspline)then
c Cunning sigma^3 method: (!) (when not aspline)
c First set up RHS of (#17) to be advected
        do lgns=1,lat2
          do mg=1,lon
            bigq=qq(mg,lgns,k)/(pstar(mg,lgns)*sig(k))**3
            qq(mg,lgns,k)=bigq
     &                   -1.5*tdt*dpsdt(mg,lgns,k)*bigq/pstar(mg,lgns)
          enddo
        enddo
        endif

c Horizontal advection.
c Arrays are 3d here and 2d in hadvect, so pass them in by the address
c of the first required element.
c First calculate coordinates of departure point

        call hadvect(ureal(1,1,k),vreal(1,1,k),tdt,   !Inputs
     &               xg,yg,idel,jdel)                         !Outputs

c Interpolate to find value of water vapour qq at departure point

        call hinterp(0,.false.,xg,yg,idel,jdel, !Inputs
     &               qq(-1,-1,k))     !In and out

c Horizontal advection of cloud water and cloud ice 

        if ( qcloud ) then
          call hinterp (1,.true.,xg,yg,idel,jdel, !Inputs
     &                  ql(-1,-1,k))             !In and out

          call hinterp (1,.true.,xg,yg,idel,jdel, !Inputs
     &                  qf(-1,-1,k))             !In and out
        endif

c Horizontal tracer advection.
        
        if ( sltrace ) then
          ntrace1=2
          do nt=ntrace1,ntrace
            call hinterp(1,.true.,xg,yg,idel,jdel, !Inputs
     &                   conm(-1,-1,k,nt)) !In and out
          enddo
        endif

c Horizontal advection for aerosol tracers
        
        if ( coupled_aero ) then
          do nt=1,ntrac
            call hinterp(1,.true.,xg,yg,idel,jdel, !Inputs
     &                   xtm(-1,-1,k,nt)) !In and out
            do lgns=1,lat2
              do mg=1,lon
                xtm(mg,lgns,k,nt)=max(0.,xtm(mg,lgns,k,nt))
                delx(mg,lgns,k,nt)=(xtm(mg,lgns,k,nt)
     &               -xtmsav(mg,lgns,k,nt))*muf(mg,k,lgns)
              enddo
            enddo
          enddo
        endif

        if(k.eq.2.and.debug ) then
          nt=2
          lgns=lgdebug*(insdebug-1)+(lat2p-lgdebug)*(2-insdebug)
          mg=mgdebug
          write(26,*)'After hadv, nt = ',nt
 1        format(3(a,g13.4))
C*** 93       format(3g13.4)
          write(26,1)'xtm ',xtm(mg,lgns,k,nt)
          write(26,1)'xg ',xg(mg+2,lgns+2),' yg ',yg(mg+2,lgns+2)
          write(26,1)'idel ',idel(mg+2,lgns+2),' jdel',jdel(mg+2,lgns+2)
        endif

c Compute pressure weighted time change in water variables (delq)
c Sum the positive and negative increments

        do lgns=1,lat2
          ns=2-(lgns-1)/lat
          lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)

c-- Water vapour
          if(.not.aspline)then
c           End of cunning sigma^3 method: apply JMcG's equation 17.
            do mg=1,lon
              qq(mg,lgns,k)=qq(mg,lgns,k)*(pstar(mg,lgns)*sig(k))**3 /
     &             (1.+1.5*tdt*dpsdt(mg,lgns,k)/pstar(mg,lgns))
            enddo
          endif
          do mg=1,lon
            ma=mg+(ns-1)*lon
            qq(mg,lgns,k)=max(qq(mg,lgns,k),0.)
            qqx=(qq(mg,lgns,k)-rmmg(ma,k,lg))*muf(mg,k,lgns)
            delq(ma,lg,k,1)=qqx
          enddo

          if(qcloud)then

            do mg=1,lon
              ma=mg+(ns-1)*lon
c-- Cloud water (kq=2)
              ql(mg,lgns,k)=max(ql(mg,lgns,k),0.)
              qqx=(ql(mg,lgns,k)-qlbm(ma,k,lg))*muf(mg,k,lgns)
              delq(ma,lg,k,2)=qqx
c-- Cloud ice   (kq=3)
              qf(mg,lgns,k)=max(qf(mg,lgns,k),0.)
              qqx=(qf(mg,lgns,k)-qfbm(ma,k,lg))*muf(mg,k,lgns)
              delq(ma,lg,k,3)=qqx
            enddo
          endif

        enddo

      return
      end
c******************************************************************************
      subroutine enforce_conq ( ncons,muf,wdsk,delq,    !In     
     &                          qqm,qq)                 !In and out

      implicit none
C Global parameters
      include 'PARAMS.f'
      real asf,asfx
      parameter(asf=0.05, asfx=1.0-2*asf) !Asselin filter coefficients

C Argument list
      integer ncons
      real muf(lon,nl,lat2)
      real wdsk(lon,nl,lat2)
      real delq(ln2,lat,nl)
      real qqm(ln2,nl,lat)
      real qq(ln2,nl,lat)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'TIMEX.f'    !Input ncepstrt

C Local work arrays and variables
      real fact(2)
      logical ncepx

      integer k
      integer lg
      integer lgns
      integer ma
      integer mg
      integer ns

      real alphax
      real delqneg
      real delqpos
      real delqq
      real facm
      real qmuf
      real qqx
      real rgtx
      real rhf

C Local data, functions etc

C Start code : ----------------------------------------------------------

c Adjust moistures to ensure conservation
c delqpos is the sum of all positive changes of moisture over globe
c delqneg is the sum of all negative changes of moisture over globe  
c Alpha is chosen to satisfy alpha*delqpos + delqneg/alpha = 0
c NB alpha=sqrt(alphax) - alpha is not used in coding due to JMcG's 
c tricky fact(2), i.e. one_on_alpha=fact(2), alpha=alphax*fact(2).


c Sum to get the global delqpos, delqneg

      delqpos=0.
      delqneg=0.

!!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!!$OMP& PRIVATE   (k, lg, lgns, ma, mg, ns, rhf)
!!$OMP& REDUCTION (+ : delqneg, delqpos)
!!$OMP& SHARED    (delq, wdsk)

c      do k=1,nl
c        do ns=1,2
c          do lg=1,lat
c            lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
      do lg=1,lat
        do ns=1,2
          lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
          do k=1,nl
            do mg=1,lon
              ma=mg+(ns-1)*lon
              rhf=max(0.0,delq(ma,lg,k))
              delqpos=delqpos+     rhf *wdsk(mg,k,lgns)
              delqneg=delqneg+(delq(ma,lg,k)-rhf)*wdsk(mg,k,lgns)
            enddo
          enddo
        enddo
      enddo

!!$OMP END PARALLEL DO

      alphax = -delqneg/delqpos
      fact(2)=1./sqrt(alphax)   ! standard 2-sided
      fact(1)=min(1.,fact(2))

      delqpos=0.
      delqneg=0.

c Now multiply through delq for all points by alpha or 1/alpha,
c depending on whether delq is positive or negative at that point.
c Double pass conservation

      if(ncons.eq.2)then        ! for double-pass conservation

        if(alphax.lt.1.)Then !i.e., if global q change is positive
c       N.B. second pass not needed if alphax.ge.1

!!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!!$OMP& PRIVATE   (delqq, facm, k, lg, lgns, ma, mg, ns, qmuf, qqx, rhf)
!!$OMP& REDUCTION (+ : delqneg, delqpos)
!!$OMP& SHARED    (alphax, delq, fact, muf, ncons, qqm, wdsk)

c          do k=1,nl
c            do ns=1,2
c              do lg=1,lat
c                lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
          do lg=1,lat
            do ns=1,2
              lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
              do k=1,nl
                do mg=1,lon
                  ma=mg+(ns-1)*lon
                  facm=1.0
                  if(delq(ma,lg,k).gt.0)facm=alphax
                  delqq=facm*delq(ma,lg,k)*fact(ncons)
                  qmuf=qqm(ma,k,lg)*muf(mg,k,lgns)
                  qqx=max(0.0,(qmuf+delqq))-qmuf
                  delq(ma,lg,k)=qqx
                  rhf=max(0.0,qqx)
                  delqpos=delqpos+     rhf *wdsk(mg,k,lgns)
                  delqneg=delqneg+(qqx-rhf)*wdsk(mg,k,lgns)
                enddo
              enddo
            enddo
          enddo

!!$OMP END PARALLEL DO

          alphax = -delqneg/delqpos
          fact(2)=1.            ! 2 pass just use simple 1-sided 2nd time

        endif                   ! (alphax.lt.1.)
      endif                     ! ncons.eq.2



c Now do time-step, applying Asselin filter.
c NB: Not saving rgt any more, since Kuo is no longer used (LDR 4/2000)
c Time-step and filter cloud water and cloud ice if qcloud (kq=2,3)
c Remove extreme small traces left after SLT (set values < 1.0e-20 to zero)

      ncepx=(ncepstrt.eq.2).or.(ncepstrt.eq.0)

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (delqq, facm, k, lg, lgns, ma, mg, ns, rgtx)
!$OMP& SHARED  (alphax, delq, fact, muf, ncepx, ncons, qq, qqm)

      do lg=1,lat
        do ns=1,2
          lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
          do k=1,nl
            do mg=1,lon
              ma=mg+(ns-1)*lon
              facm=1.0
              if(delq(ma,lg,k).gt.0)facm=alphax
              delqq=facm*delq(ma,lg,k)*fact(ncons)
              rgtx=qqm(ma,k,lg)+delqq/muf(mg,k,lgns)
              if(rgtx.lt.1.0e-20)rgtx=0.0
              if(.not.ncepx)
     &        qqm(ma,k,lg)=asfx*qq(ma,k,lg)+
     &             asf*qqm(ma,k,lg)
              qq(ma,k,lg)=rgtx
            enddo
          enddo
        enddo
      enddo

!$OMP END PARALLEL DO

      return
      end
c******************************************************************************
      subroutine enforce_cont ( ncons,muf,wdsk,delq,    !In     
     &                          qqm,qq)                 !In and out

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'ECPARM.f'
      real asf,asfx
      parameter(asf=0.05, asfx=1.0-2*asf) !Asselin filter coefficients

C Argument list
      integer ncons
      real muf(lon,nl,lat2)
      real wdsk(lon,nl,lat2)
      real delq(lone,lat2e,nl)
      real qqm(-1:lone+2,-1:lat2e+2,nl)
      real qq(-1:lone+2,-1:lat2e+2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'TIMEX.f'    !Input ncepstrt

C Local work arrays and variables
      real fact(2)
      logical ncepx

      integer k
      integer lgns
      integer mg

      real alphax
      real delqneg
      real delqpos
      real delqq
      real facm
      real qmuf
      real qqx
      real rgtx
      real rhf

C Local data, functions etc

C Start code : ----------------------------------------------------------

c This routine is essentially the same as enforce_conq
c Only real diff is that arrays are organised differently

c Adjust tracers to ensure conservation
c delqpos is the sum of all positive changes of tracer over globe
c delqneg is the sum of all negative changes of tracer over globe  
c Alpha is chosen to satisfy alpha*delqpos + delqneg/alpha = 0
c NB alpha=sqrt(alphax) - alpha is not used in coding due to JMcG's 
c tricky fact(2), i.e. one_on_alpha=fact(2), alpha=alphax*fact(2).


c Sum to get the global delqpos, delqneg

      delqpos=0.
      delqneg=0.
      
      do k=1,nl
        do lgns=1,lat2
          do mg=1,lon
            rhf=max(0.0,delq(mg,lgns,k))
            delqpos=delqpos+     rhf *wdsk(mg,k,lgns)
            delqneg=delqneg+(delq(mg,lgns,k)-rhf)*wdsk(mg,k,lgns)
          enddo
        enddo
      enddo


      alphax = -delqneg/delqpos
      fact(2)=1./sqrt(alphax)   ! standard 2-sided
      fact(1)=min(1.,fact(2))

      delqpos=0.
      delqneg=0.

c Now multiply through delq for all points by alpha or 1/alpha,
c depending on whether delq is positive or negative at that point.
c Double pass conservation

      if(ncons.eq.2)then        ! for double-pass conservation

        if(alphax.lt.1.)Then !i.e., if global q change is positive
c       N.B. second pass not needed if alphax.ge.1

          do k=1,nl
            do lgns=1,lat2
              do mg=1,lon
                facm=1.0
                if(delq(mg,lgns,k).gt.0)facm=alphax
                delqq=facm*delq(mg,lgns,k)*fact(ncons)
                qmuf=qqm(mg,lgns,k)*muf(mg,k,lgns)
                qqx=max(0.0,(qmuf+delqq))-qmuf
                delq(mg,lgns,k)=qqx
                rhf=max(0.0,qqx)
                delqpos=delqpos+     rhf *wdsk(mg,k,lgns)
                delqneg=delqneg+(qqx-rhf)*wdsk(mg,k,lgns)
              enddo
            enddo
          enddo

          alphax = -delqneg/delqpos
          fact(2)=1.            ! 2 pass just use simple 1-sided 2nd time

        endif                   ! (alphax.lt.1.)
      endif                     ! ncons.eq.2



c Now do time-step, applying Asselin filter.
c NB: Not saving rgt any more, since Kuo is no longer used (LDR 4/2000)
c Time-step and filter cloud water and cloud ice if qcloud (kq=2,3)
c Remove extreme small traces left after SLT (set values < 1.0e-20 to zero)

      ncepx=(ncepstrt.eq.2).or.(ncepstrt.eq.0)
      do k=1,nl
        do lgns=1,lat2
          do mg=1,lon
            facm=1.0
            if(delq(mg,lgns,k).gt.0)facm=alphax
            delqq=facm*delq(mg,lgns,k)*fact(ncons)
            rgtx=qqm(mg,lgns,k)+delqq/muf(mg,k,lgns)
            if(.not.ncepx)
     &      qqm(mg,lgns,k)=asfx*qq(mg,lgns,k)+
     &           asf*qqm(mg,lgns,k)
            qq(mg,lgns,k)=rgtx
          enddo
        enddo
      enddo

      return
      end
