c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c REDUCTION operation commented out, for reproducibility.
c SJP 2004/01/05
c
c Manually parallelised.
c SJP 2003/04/24
c
c $Log: uvreal.f,v $
c Revision 1.6  2001/11/07 07:02:00  rot032
c Do max wind diagnostic at all levels.
c
c Revision 1.5  1996/10/24 01:03:22  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1996/03/21  03:19:12  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.3  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.2  1993/07/06  16:39:31  ldr
c      Tidied up global_umax calculation (HBG).
c
c Revision 1.1  93/03/15  15:45:19  ldr
c Initial revision
c 
      subroutine uvreal(ind,global_umax)
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'UVPGD.f'
      include 'GAUSL.f'
      include 'HYBARR.f'

      if(ind.eq.1)then
c---- Replace the u,v weighted values by their non-weighted equivalents
c---- (to be done before SLT scheme, and before calculation of Umx)

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (csqr, csqr2, ins, k, lg, lgns, mg, pgdx)
!$OMP$ SHARED  (acsq, dadnf, dbdnf, pgd, ugd, ureal, vgd, vreal)

        do 20 lg=1,lat
          csqr2=acsq(lg)*eradsq
          csqr=sqrt(csqr2)
        do 20 ins=1,2
          lgns=lg*(ins-1)+(lat2p-lg)*(2-ins)
        do 40 k=1,nl
           do 40 mg=1,lon
           pgdx=dadnf(k)+dbdnf(k)*pgd(mg,lgns)
           ureal(mg,lgns,k)=ugd(mg,lgns,k)*csqr/pgdx
   40      vreal(mg,lgns,k)=vgd(mg,lgns,k)*csqr/pgdx
   20   continue

!$OMP END PARALLEL DO

      end if

      if(ind.eq.2)then
c.... find maximum wind
        global_umax=0.

!!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!!$OMP& DEFAULT   (NONE)
!!$OMP& PRIVATE   (k, lgns, mg)
!!$OMP& REDUCTION (max : global_umax)
!!$OMP$ SHARED    (ureal, vreal)

        do 50 k=1,nl
        do 50 lgns=1,lat2
        do 50 mg=1,lon
        global_umax=max(global_umax,
     &     sqrt(ureal(mg,lgns,k)**2+vreal(mg,lgns,k)**2)  )
   50   continue

!!$OMP END PARALLEL DO

      end if

      return
      end
