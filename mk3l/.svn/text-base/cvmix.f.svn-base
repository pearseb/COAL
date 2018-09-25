c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: cvmix.f,v $
c Revision 1.21  2001/06/04 02:27:04  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.20  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.19  2000/06/20 02:08:35  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.18.1.1  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.18  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.17  1998/12/10  00:55:58  ldr
c HBG changes to V5-1-21
c
c Revision 1.16  1997/12/17  23:23:01  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.15  1996/10/24  01:02:35  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.14  1996/06/13  02:06:07  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.13  1996/02/19  04:09:44  ldr
c Generalize for 24 levels.
c
c Revision 1.12  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.11  1994/08/08  17:20:57  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.10  94/08/08  13:16:21  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.9  93/06/16  14:09:59  ldr
c Now solved implicitly (using trim) - HBG.
c 
c Revision 1.8  93/01/26  17:26:57  ldr
c  Merge of tracer changes to V4-1 with other changes since V4-1.
c 
c Revision 1.7  93/01/26  16:52:58  ldr
c Corrected cvmix (with simple coding) and compatible conv.
c 
c Revision 1.5  92/11/20  15:07:02  ldr
c Vertical mixing of momentum by convection now done by adjacent level mixing
c (i.e. none from cloud base to cloud top.) Merged with SPO's change from
c ron, son variables to radd, sadd (additions to ron. son) by LDR.
c 
c Revision 1.3  92/06/16  12:02:32  ldr
c Rewritten to pass physical variables as arguments rather than in common,
c in order to implement sea-ice model.
c 
c Revision 1.2  91/03/13  12:57:30  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c
c     INPUT/OUTPUT:
c     Input:   from common/fewflags in FEWFLAGS.f
c                  sltrace  - if T compute tracer transport by semi-
c                             Langrangian advection and by mixing and
c                             convection in the vertical
c
c              from common/hybrpr in HYBRPR.f
c                  dprf - pressure thickness at each sigma level
c
c              from common/traceblk in TRACEBLK.f
c                  con - concentration of atmospheric tracers at time, t
c
c              from arguments
c                  fmp - convective mass flux
c                  kba - level of cloud base   kta - level of cloud top
c                  lg  - latitude index        pg  - surface pressure(mbs)
c                  u   - zonal wind velocity m/s
c                  v   - meridional wind velocity m/s
c                  tdt - 2 X timestep (dt)
c
c     In/Out:  from arguments
c                  cten - rate of change of tracer concentration due to 
c                         convection
c                  uten, vten  - effective time derivative of u, v
c                  dkevm - change in KE  
c                  xtg - ECHAM model tracer mixing ratio (kg/kg)
c
      subroutine cvmix(ipass,lg,fmp,fmu,fmd,kta,kba,u,v,pg,tdt,fscav, !Inputs
     &                 uten,vten,dkevm,cten,xtg,conwd,      !In and out
     &                 xtu)                                       !Output

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'ECPARM.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real fmp(ln2,nl)
      real fmu(ln2,nl)            ! updraught flux (hPa/s) +ve
      real fmd(ln2,nl)            ! downdraught flux (hPa/s) +ve
      integer kta(ln2)
      integer kba(ln2)
      real u(ln2,nl)
      real v(ln2,nl)
      real pg(ln2)
      real tdt
      real fscav(ln2,nl,ntrac) !Convective tracer scavenging fraction (b/w 0 and 1)
      real uten(ln2,nl)
      real vten(ln2,nl)
c     real cten(ln2,nl,ntrace) ! Must be after TRACEBLK.f (below)
      real dkevm(ln2,nl)
      real xtg(ln2,nl,ntrac)
      real xtu(ln2,nl,ntrac)
      real conwd(ln2,ntrac)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Global data blocks
      include 'FEWFLAGS.f' !Input sltrace
      include 'TRACEBLK.f'
      real cten(ln2,nl,ntrace)
c      real omegac(ln2,nl,lat)
c      common/omblk/omegac


C Local work arrays and variables
      real au(ln2,nl),cu(ln2,nl),rhs(ln2,nl)
      real wke(ln2,nl),wkt(ln2,nl),wkb(ln2,nl)
      real fout(ln2,nl),fthru(ln2,nl),dxd(ln2,nl),fu(ln2,nl)
      real fmav(ln2,0:nl) !Updraft mass flux at half levels
      real xscav(ln2,nl,ntrac),dxs(ln2,nl,ntrac)

C Local data, functions etc

C Start code : ----------------------------------------------------------

C****
C**** CALCULATE ADDITIONS TO TENDENCIES UTEN,VTEN DUE TO CONVECTIVE MIXING.
c**** All additions to RON, SON are made later in radin.
C****
C**** VERTICAL MIXING OF MOMENTUM BY CONVECTION
C****   MASS FLUX FMP=0 FOR NON CONVECTIVE LEVELS, >0 FOR Kb TO Kt
C****

c---- new method from Mar 93
c---- now only mixes adjacent levels, no mixing kb to kt
c---- note that fmp has been factored by 0.5 by half due to adjacent
c---- level (double) mixing.
c---- Implicit method now in use.

c---- set up tri-diag matrices.
         fact=0.5
      do 55 mg=1,ln2
      au(mg,1)=0.0
   55 cu(mg,nl)=0.0
      do 60 k=2,nl
      do 60 mg=1,ln2
      fmpav=fact*0.5*(fmp(mg,k)+fmp(mg,k-1))*tdt
      dskx=dprf(mg,k)/pg(mg)
      au(mg,k)=-fmpav/dskx
      dskmx=dprf(mg,k-1)/pg(mg)
   60 cu(mg,k-1)=-fmpav/dskmx
      do 62 mg=1,ln2
      kb=kba(mg)
      if(kb.eq.0)go to 62
      kt=kta(mg)
      au(mg,kb)=0.0
      cu(mg,kb-1)=0.0
      cu(mg,kt)=0.0
      au(mg,kt+1)=0.0
   62 continue

c---- Do the U field
      do 65 k=1,nl
      do 65 mg=1,ln2
   65 rhs(mg,k)=u(mg,k)
      call trim(rhs,au,cu,rhs,0,wke,wkt,wkb)
      do 70 k=1,nl
      do 70 mg=1,ln2
      uten(mg,k)=uten(mg,k)+(rhs(mg,k)-u(mg,k))/tdt
   70 u(mg,k)=rhs(mg,k)
c---- Do the V field
      do 75 k=1,nl
      do 75 mg=1,ln2
   75 rhs(mg,k)=v(mg,k)
      call trim(rhs,au,cu,rhs,1,wke,wkt,wkb)
      do 80 k=1,nl
      do 80 mg=1,ln2
      vten(mg,k)=vten(mg,k)+(rhs(mg,k)-v(mg,k))/tdt
   80 v(mg,k)=rhs(mg,k)

c---- Compute the KE change per pairs of levels for calculating heating
c---- caused by vertical mixing
      do 90 k=1,nl-1
      do 90 mg=1,ln2
      adkedt=(dprf(mg,k)*cu(mg,k)/tdt)*
     & ((u(mg,k)-u(mg,k+1))**2+(v(mg,k)-v(mg,k+1))**2)
     &  /(dprf(mg,k)+dprf(mg,k+1))
      dkevm(mg,k  )=dkevm(mg,k  )+adkedt
   90 dkevm(mg,k+1)=dkevm(mg,k+1)+adkedt

c Vertical mixing of aerosols by cumulus convection (rewritten by LDR 5/01).

      if(coupled_aero.and.ipass.eq.1) then

c       omegac(:,:,lg)=fmu(:,:) !updraft massflux in hPa/s

c Calculate the updraft mass flux at the half levels
c fmav(:,k) is the flux at top of layer k

c Whole column version...

        fmav(:,0)=0.
        fmav(:,nl)=0.
        do k=1,nl-1
          fmav(:,k)=0.5*(fmu(:,k)+fmu(:,k+1))
        enddo

c Inside convective cloud version...

C***        fmav(:,:)=0.
C***        do mg=1,ln2
C***          if(kba(mg).gt.0)then
C***            do k=kba(mg)-1,kta(mg)-1
C***              fmav(mg,k)=0.5*(fmu(mg,k)+fmu(mg,k+1))
C***            enddo
C***          endif
C***        enddo

c Do the updraft mixing, including effect of compensating subsidence.
c First, precalculate the fraction leaving (fout) and fraction going through (fthru)
c at each level.

        fout(:,:)=0.
        fthru(:,:)=1.


c Work out the fractions entrained and detrained in each layer
c Coding assumes mass flux fmu is zero in top layer, or a lack of conservation will occur.

        do k=1,nl
          do mg=1,ln2
            if(fmav(mg,k).gt.fmav(mg,k-1))then !Entrainment from layer k
              cdt=(fmav(mg,k)-fmav(mg,k-1))*tdt/dprf(mg,k) !fmav in mb/s
c              fout(mg,k)=1.-exp(-cdt)
              fout(mg,k)=min(cdt/(1.+0.5*cdt), 1.)
c              fout(mg,k)=min(cdt, 1.)
              fthru(mg,k)=1.
            elseif(fmav(mg,k-1).gt.0.)then    !Detrainment into layer k
              fout(mg,k)=0.
              fin=max(0., (fmav(mg,k-1)-fmav(mg,k))/fmav(mg,k-1))
              fthru(mg,k)=1.-fin
            endif
          enddo
        enddo

c Now, apply the tendencies, including the effect of scavenging.
c Criterion kba.gt.0 means that it is only applied when moist convection is diagnosed.

        xscav(:,:,:)=0.
        do nt=1,ntrac
          dxd(:,:)=0.
          fu(:,:)=0.
          do k=1,nl-1
            do mg=1,ln2
              if(kba(mg).gt.0)then
                dp=dprf(mg,k)
                dpp=dprf(mg,k+1)
                xscav(mg,k,nt)=fu(mg,k)*fscav(mg,k,nt)/dp  !Mixing ratio scavenged in layer k
                fu(mg,k)=fu(mg,k)*(1.-fscav(mg,k,nt))      !Flux going thru layer k
                dxe=xtg(mg,k,nt)*fout(mg,k)                !Entrainment in layer k
                fu(mg,k+1)=(fu(mg,k)+dxe*dp)*fthru(mg,k+1) !Flux going thru layer k+1
                dxd(mg,k+1)=(fu(mg,k)+dxe*dp)*(1-fthru(mg,k+1))/dpp !Detrainment
                xtg(mg,k,nt)=xtg(mg,k,nt)-dxe+dxd(mg,k)
              endif
            enddo
          enddo
        enddo

c Save estimate of mixing ratio in convective updraft for xtchemie.

        do nt=1,ntrac
          do k=1,nl
            do mg=1,ln2
              xtu(mg,nlp-k,nt)=xtg(mg,k,nt)
            enddo
          enddo
        enddo


c Work out the tendencies due to the compensating subsidence term.

        dxs(:,:,:)=0.
        do nt=1,ntrac
          do k=nl-1,1,-1
            do mg=1,ln2
              if(kba(mg).gt.0)then
                dp=dprf(mg,k)
                dpp=dprf(mg,k+1)
                cdt=fmav(mg,k)*tdt/dprf(mg,k+1)
                fsub=min(cdt/(1.+0.5*cdt), 1.)
                dxs(mg,k+1,nt)=dxs(mg,k+1,nt)-fsub*xtg(mg,k+1,nt)
                dxs(mg,k,nt)=dxs(mg,k,nt)+fsub*xtg(mg,k+1,nt)*dpp/dp
              endif
            enddo
          enddo
        enddo
        
        xtg(:,:,:)=xtg(:,:,:)+dxs(:,:,:)


c Calculate the total convective scavenging

        do nt=1,ntrac
          if(lwetdep(nt))then
            do k=1,nl
              do mg=1,ln2
                conwd(mg,nt)=conwd(mg,nt)
     &                      +xscav(mg,k,nt)*100.*dprf(mg,k)/(grav*tdt)
              enddo
            enddo
          endif
        enddo


        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            write(25,'(3(a,i2))')'After cvmix updraught... '
            write(25,91)'fout ',(fout(mg,k),k=1,nl)
            write(25,91)'fthru ',(fthru(mg,k),k=1,nl)
            write(25,91)'fs so2 ',(fscav(mg,k,2),k=1,nl)
            write(25,91)'fs so4 ',(fscav(mg,k,3),k=1,nl)
            write(25,1)'SO2 convwetdep (kg/m2) = ',conwd(mg,2)*tdt
            write(25,1)'SO4 convwetdep (kg/m2) = ',conwd(mg,3)*tdt
            write(25,*)
          endif
        endif
 1      format(3(a,g10.3))
 91     format(a,30f10.3)


C***        if(debug)then
C***          if(lg.eq.lgdebug)then
C***            ns=insdebug
C***            mg=mgdebug+(ns-1)*lon
C***            write(25,'(3(a,i2))')'After cvmix downdraught... '
C***            write(25,91)'fout ',(fout(mg,k),k=1,nl)
C***            write(25,91)'fthru ',(fthru(mg,k),k=1,nl)
C***            write(25,*)
C***          endif
C***        endif

      endif

c Vertical mixing of tracers by cumulus convection (the old tracer scheme used by IGW)

      if ( sltrace ) then
        lgn=lat2p-lg

        do nt=1,ntrace
          do k=1,nl
          do mg=1,lon
            rhs(mg,k)=con(mg,lgn,k,nt)
            rhs(mg+lon,k)=con(mg,lg,k,nt)
          enddo
          enddo
          call trim(rhs,au,cu,rhs,1,wke,wkt,wkb)
          do k=1,nl
          do mg=1,lon
            cten(mg,k,nt)=cten(mg,k,nt)+
     &            (rhs(mg,k)-con(mg,lgn,k,nt))/tdt
            cten(mg+lon,k,nt)=cten(mg+lon,k,nt)+
     &            (rhs(mg+lon,k)-con(mg,lg,k,nt))/tdt
          enddo
          enddo
        enddo

      endif


      return
      end
