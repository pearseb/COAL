c Form two different versions of one of the parallel regions, one for qcloud=T
c and one for qcloud=F. This enables us to retain the previous loop fusion,
c while ensuring correct results when qcloud=F.
c SJP 2004/05/25
c
c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Manually parallelised. As we always set qcloud=T, we merge two of the loops
c for optimal performance.
c SJP 2003/04/24
c
c $Log: assel.f,v $
c Revision 1.15  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.14  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.13  2000/06/20 02:27:10  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.12  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.11  1996/10/24  01:02:28  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.10  1994/05/13  14:54:58  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c
c Revision 1.9  93/12/06  16:55:20  ldr
c Changes to V4-4-45l to get T63 running on SGI (Seca).
c 
c Revision 1.8  93/10/05  13:05:22  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.7  93/09/15  17:09:40  ldr
c Changes to get inital version of cloud water scheme going.
c 
c Revision 1.6  93/01/26  16:26:38  ldr
c Changes to implement Semi-Lagrangian transport of tracers.
c 
c Revision 1.5  92/12/09  14:42:49  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.4  92/04/16  16:32:58  ldr
c Corrected declaration of pm in block fieldm. (No change to answers.)
c 
c Revision 1.3  92/03/19  11:37:54  ldr
c Corrected order of subscripts in /field/ (Error didn't affect results). 
c 
c Revision 1.2  91/03/13  12:56:08  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:36:46  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT:
c     Input:   from common/fewflags in FEWFLAGS.f
c                  ifwds   - default F, only used when doing forward start
c                            from new restart file
c                  qcloud  - T if using prognostic cloud scheme
c                  sltrace - if T compute tracer transport by semi-
c                            Langrangian advection and by mixing and
c                            convection in the vertical
c
c              from common/fldri in FLDRI.f
c                  tei - spectral temp, imaginary part
c                  ter - spectral temp, real part
c
c              from common/qcloud1 in QCLOUD1.f
c                  qfb - cloud ice at t 
c                  qlb - cloud liquid water at t
c
c              from common/rmgrid in RMGRID.f
c                  rmg - pressure weighted moistures at current timestep
c
c              from com mon/traceblk in TRACEBLK.f
c                  conm - concentration of atmospheric tracers at t-1
c
c     In/Out:  from common/fldmri in FLDMRI.f
c                  temi - spectral temp, imaginary part at t-1
c                  temr - spectral temp, real part at t-1
c
c              from common/qcloud1 in 
c                  qfbm - cloud ice at t-1
c                  qlbm - cloud liquid water at t-1
c
c              from common/rmgrid in RMGRID.f
c                  rmmg - pressure weighted mositures at previous timestep
c
c              from common/traceblk in TRACEBLK.f
c                  con - concentration of atmospheric tracers at t
c 
C$MP_SCHEDTYPE=INTERLEAVE

      subroutine assel

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'ECPARM.f'

C Argument list

C Global data blocks
      include 'ECTRAC.f'
      include 'FEWFLAGS.f'  !Input ifwds, sltrace, qcloud
      include 'FLDRI.f'
      include 'FLDMRI.f'
      include 'QCLOUD1.f'   !Qcloud arrays
      include 'RMGRID.f'
      include 'TIMEX.f'
      include 'TRACEBLK.f'  !Tracer arrays con, conm

C Local work arrays and variables
      integer k
      integer lg
      integer lgns
      integer ll
      integer mg
c     integer mm ! defined in LLMAX
      integer nt

      real asf

C Local data, functions etc
      include 'LLMAX.f'     !Statement function

C Start code : ----------------------------------------------------------

C**** TO ADD THE PHYSICALLY ADJUSTED TEMP AND
C**** MOISTURE COMPONENTS OF THE ASSELIN FILTER.
C**** T+1 ASSELIN FILTER TERM FOR TEMP ADDED( NOW AT T )
      asf=0.05
      if((ncepstrt.eq.2).or.(ncepstrt.eq.0))asf=0.0

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (k, ll, mm)
!$OMP& SHARED  (asf, tei, temi, temr, ter)

      do 10 k=1,nl
      do 10 mm=1,mw
      do 10 ll=1,llmax(mm)
      temr(ll,mm,k)=temr(ll,mm,k)+asf*ter(ll,mm,k)
   10 temi(ll,mm,k)=temi(ll,mm,k)+asf*tei(ll,mm,k)

!$OMP END PARALLEL DO

C**** ADD GRID POINT MOISTURE ASSELIN FILTER
C**** Prognostic cloud scheme stuff moved to here. SJP 2003/04/24

      if (qcloud) then

c  Apply Asselin filter to cloud liquid water and cloud frozen water

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (k, lg, mg)
!$OMP& SHARED  (asf, qfb, qfbm, qlb, qlbm, rmg, rmmg)

      do 12 lg=1,lat
      do 12 k=1,nl
      do 12 mg=1,ln2
      qlbm(mg,k,lg)=qlbm(mg,k,lg)+asf*qlb(mg,k,lg)
      qfbm(mg,k,lg)=qfbm(mg,k,lg)+asf*qfb(mg,k,lg)
   12 rmmg(mg,k,lg)=rmmg(mg,k,lg)+asf*rmg(mg,k,lg)
 
!$OMP END PARALLEL DO

      else

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (k, lg, mg)
!$OMP& SHARED  (asf, rmg, rmmg)

      do 13 lg=1,lat
      do 13 k=1,nl
      do 13 mg=1,ln2
   13 rmmg(mg,k,lg)=rmmg(mg,k,lg)+asf*rmg(mg,k,lg)

!$OMP END PARALLEL DO

      end if

c Apply Asselin filter to tracer array

      if ( sltrace ) then

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (k, lgns, mg, nt)
!$OMP& SHARED  (asf, con, conm)

        do nt=1,ntrace
          do k=1,nl
            do lgns=1,lat2
              do mg=1,lon
                conm(mg,lgns,k,nt)=conm(mg,lgns,k,nt)
     &                            +asf*con(mg,lgns,k,nt)
              enddo
            enddo
          enddo
        enddo

!$OMP END PARALLEL DO

      endif

c  Apply Asselin filter to ECHAM tracer fields

      if (coupled_aero) then
        bigxtm(:,:,:,:)=bigxtm(:,:,:,:)+asf*bigxt(:,:,:,:)
      endif

      return
      end
