c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/03/29
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: just_fm.f,v $
c Revision 1.3  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.2  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.1  2001/02/12 05:36:25  rot032
c Initial revision
c
      subroutine just_fm(ind)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ind

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'LSMI.f'
      include 'MASIV4.f'
      include 'RADAV.f'
      include 'TIMEX.f'

      real ochf,occur,g50dt
      integer khf
      common/curent/ochf(ln2,lat),occur(ln2,lat),khf(ln2,lat)
     &,g50dt(ln2,lat)

C Local work arrays and variables
      real emja(ln2)

      integer lg
      integer ma
      integer mg
      integer ns

      real emj
      real sumemj
      real tintp
      real tmlo
      real tst1
      real tst2

C Local data, functions etc
      real hcap
      data hcap/2.095e8/

C Start code : ----------------------------------------------------------

C This routine deals with sea points that have just frozen
C   or ice points that have just melted. Called by phys.f or physseca.f

      if(ind.eq.2)then

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (emj, emja, lg, ma, mg, ns, sumemj, tintp, tmlo, tst1,
!$OMP&          tst2)
!$OMP& SHARED  (dt, flxmjz, hcap, imsl, just_frozen, just_melted,
!$OMP&          occur, ratlm, savegrid)

        do lg=1,lat-1

         do mg=1,ln2
c Just frozen point - if equatorward point is a SEA point
c   then change to a MLO point
          if ( just_frozen(mg,lg) ) then
           if(imsl(mg,lg+1).eq.3) imsl(mg,lg+1)=-2
           just_frozen(mg,lg)=.false.
           just_melted(mg,lg)=.false.
          endif
         enddo

         do mg=1,ln2
          emja(mg)=0.0
c Just melted point - if equatorward point is a MLO point
c   then change to a SEA point
          if ( just_melted(mg,lg) ) then
           if(imsl(mg,lg+1).eq.2) then
            imsl(mg,lg+1)=3
c.... Equatorward MLO point changed to a SEA point. Temperature jump.
c.... Account for energy.
            tmlo=savegrid(mg,3,lg+1)
            tst1=savegrid(mg,5,lg+1)
            tst2=savegrid(mg,6,lg+1)
            tintp=tst1*(ratlm-1.0)-tst2*ratlm
c.... mixed layer depth now 100m (was 50m) - multiply hcap by 2
            emj=hcap*2.*(tintp-tmlo)/dt
c.... Add as term in calculation of Q-flux
            occur(mg,lg+1)=occur(mg,lg+1)+emj
            emja(mg)=emj
           endif
           just_melted(mg,lg)=.false.
          endif
         enddo

         do ns=1,2
          sumemj=0.0
          ma=(ns-1)*lon
          do mg=1+ma,lon+ma
           sumemj=sumemj+emja(mg)
          enddo
          flxmjz(lg,ns)=flxmjz(lg,ns)+sumemj/lon
         enddo

        enddo

!$OMP END PARALLEL DO

        return

      endif

      if(ind.eq.1)then

C Set logical arrays to false (done at start of model run only)

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (lg, mg)
!$OMP& SHARED  (just_frozen, just_melted)

        do lg=1,lat
         do mg=1,ln2
          just_frozen(mg,lg)=.false.
          just_melted(mg,lg)=.false.
         enddo
        enddo

!$OMP END PARALLEL DO

      endif

      return
      end
