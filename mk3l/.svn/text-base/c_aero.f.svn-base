c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c $Log: c_aero.f,v $
c Revision 1.10  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.9  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.8.1.1  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.8  2001/06/04 02:26:54  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.7  2001/03/07 04:28:58  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.6  2001/02/22 06:46:47  rot032
c Merge LDR and HBG changes.
c
      subroutine c_aero1(lg,tdt,qtg,qlg,qfg,u,v,  !Inputs
     &                   mcmax,xtg)               !Outputs

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /PHYSICAL/ )
!$OMP THREADPRIVATE ( /RVARSG/ )
!$OMP THREADPRIVATE ( /STATS/ )
!$OMP THREADPRIVATE ( /SURFBC/ )
!$OMP THREADPRIVATE ( /WORK1/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f' !Coupled aerosol treatment (ECHAM) parameters

C Argument list ! May be used for some entry points only
      integer lg
      real tdt
      real xtg(ln2,nl,ntrac)   !Tracer mixing ratio [kg/kg]
      real pxtg(ln2,nl,ntrac)  !Tracer mixing ratio over leads [kg/kg]
      real xtem(ln2,ntrac)     !Sfc. flux of tracer passed to vertical
                               !                          mixing [kg/m2/s]
      real mcmax(ln2)          !Maximum skin reservoir depth passed back
                               !                          from surfa
      real so2wd(ln2)          !Diagnostic only: SO2 wet deposition
      real so4wd(ln2)          !Diagnostic only: SO4 wet deposition
      real conwd(ln2,ntrac)    !Diagnostic only: Convective wet deposition
      real so2oh(ln2)          !Diagnostic only: SO2 oxidation by OH
      real so2h2(ln2)          !Diagnostic only: Aqueous SO2 oxidation by H202
      real so2o3(ln2)          !Diagnostic only: Aqueous SO2 oxidation by O3
      real dmsoh(ln2)          !Diagnostic only: DMS oxidation by OH
      real dmsn3(ln2)          !Diagnostic only: DMS oxidation by NO3
      real so2dd(ln2)          !Diagnostic only: SO2 dry deposition
      real so4dd(ln2)          !Diagnostic only: SO4 dry deposition
      real sem(ln2)            !Diagnostic only: S emission
      real tg(ln2)
      real v10m(ln2)
      real eg(ln2)
      real fg(ln2)
      real rcondx(ln2)      
c      real Veff(ln2)  
      real qtg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real u(ln2,nl)
      real v(ln2,nl)
      real ompl(ln2)
      real sicef(ln2)
      real ustar(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)

      include 'HYBRPR.f'
      include 'LOGIMSL.f'
      include 'PHYSICAL.f'
      include 'RVARSG.f'

      real hevap,hrn,hflux,hcld,hcll,hclm,hclh
     & ,hrg,hrt,hsout,hsg,htst,htscrn,hwfg,hwfb,hsnd,hsid
     & ,hrtclr,hsoutclr,htb2,htb3,hvmod,htaux,htauy,hrnc
     & ,hrnf,hintrc,hals,htgg,htgf,hpev
     & ,hsev,hrgclr,hsgclr,hrgdn,hsgdn
     & ,hperc,hclc,hdtm,hqlp,hpwc,hrevap
     & ,hssubl,hpreci,hreffl,hcliq,hdms
     & ,hso2,hso4,hso2dd,hso4dd,hso2wd
     & ,hso4wd,hsem,hso2oh,hso2h2,hso2o3,hqfp,hdms1,hso21,hso41
     & ,hdmsoh,hdmsn3,hscwd
      common/stats/hevap(ln2),hrn(ln2),hflux(ln2)
     & ,hcld(ln2),hcll(ln2),hclm(ln2),hclh(ln2)
     & ,hrg(ln2),hrt(ln2)
     & ,hsout(ln2),hsg(ln2)
     & ,htst(ln2),htscrn(ln2)
     & ,hwfg(ln2),hwfb(ln2),hsnd(ln2),hsid(ln2)
     & ,hrtclr(ln2),hsoutclr(ln2),htb2(ln2),htb3(ln2)
     & ,hvmod(ln2),htaux(ln2),htauy(ln2)
     & ,hrnc(ln2)
     & ,hrnf(ln2),hintrc(ln2),hals(ln2),htgg(ln2),htgf(ln2),hpev(ln2)
     & ,hsev(ln2),hrgclr(ln2),hsgclr(ln2),hrgdn(ln2),hsgdn(ln2)
     & ,hperc(ln2),hclc(ln2),hdtm(ln2),hqlp(ln2),hpwc(ln2),hrevap(ln2)
     & ,hssubl(ln2),hpreci(ln2),hreffl(ln2),hcliq(ln2),hdms(ln2)
     & ,hso2(ln2),hso4(ln2),hso2dd(ln2),hso4dd(ln2),hso2wd(ln2)
     & ,hso4wd(ln2),hsem(ln2),hso2oh(ln2),hso2h2(ln2),hso2o3(ln2)
     & ,hqfp(ln2),hdms1(ln2),hso21(ln2),hso41(ln2),hdmsoh(ln2)
     & ,hdmsn3(ln2),hscwd(ln2)

      real egg,condxg,tsigmf,evapxf,ewwp,pmcmax
      common/surfbc/egg(ln2),condxg(ln2),tsigmf(ln2),evapxf(ln2)
     &             ,ewwp(ln2),pmcmax(ln2)

      include 'WORK1.f'

C Global data blocks
      include 'ECTRAC.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'LSMI.f'
      include 'SURF1.f'
      include 'TIMEX.f'

      real prhos,pcss,pcnsd,pbch,pssat,pswilt,phsbh,psfc,psucs
     & ,phyds,pcdr3,prlaim,prlais,pslveg,pscveg
      integer ip2bp3,ipbp2
      common/psoilvegpr/prhos(ln2,lat),pcss(ln2,lat),
     &pcnsd(ln2,lat),pbch(ln2,lat),pssat(ln2,lat),
     &pswilt(ln2,lat),phsbh(ln2,lat),ip2bp3(ln2,lat),
     &ipbp2(ln2,lat),psfc(ln2,lat),psucs(ln2,lat),
     &phyds(ln2,lat),pcdr3(ln2,lat),prlaim(ln2,lat),
     &prlais(ln2,lat),pslveg(ln2,lat),pscveg(ln2,lat)

C Local work arrays and variables

c-    Coupled aerosol treatment

      real xtm1(ln2,nl,ntrac)  !Ditto, vertically inverted for ECHAM code
      real xte(ln2,nl,ntrac)   !Tracer tendencies [kg/kg/s]
      real rho1(ln2)           !Density of air in surface layer
      real mc1(ln2)            !Skin reservoir depth passed back from surfa
      real wgmax(ln2)          !field capacity of soil passed back from surfa
      real aphp1(ln2,nlp)      !P on half levels, inverted k index and SI units

      integer k
      integer mg
      integer nt
      integer ns

      real qint
      real xint1
      real xint2
      real xint3
      real rhodz
      real rhodz1
      real rhodz2

C Local data, functions etc

C Start code : ----------------------------------------------------------

c Note that the tracer arrays bigxt, xtg are stored
c with the usual CSIRO usage (k=1 at surface).

        do nt=1,ntrac
          do k=1,nl
            do mg=1,ln2
              xtg(mg,k,nt)=max(0., bigxt(mg,k,nt,lg))
            enddo
          enddo
        enddo

        do mg=1,ln2
          mcmax(mg)=0.04 !Do this to stop xtemiss getting confused over bare-land points
        enddo

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            open(25,file='debug.'//runtype)
            write(25,'(a,i4)')'NSTEPS = ',nsteps
            write(25,'(a,3i3)')'IPASS=1, before physics. NS,LG,MG= ',
     &           ns,lg,mgdebug
            write(25,'(a,i2)')'imsl ',imsl(mg,lg)
            write(25,1)'snowd ',snowd(mg),' siced ',siced(mg),
     &           ' pg-1000 ',pn(mg)-1000
            write(25,1)'wg ',wg(mg),' wb ',wb(mg,ms,lg),' tstar ',
     &      tstar(mg)
            write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
            write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
            write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
            write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
            write(25,91)'u ',(u(mg,k),k=1,nl)
            write(25,91)'v ',(v(mg,k),k=1,nl)
            write(25,*)
            if(coupled_aero)then
              write(25,9)'DMS ',(xtg(mg,k,1),k=1,nl)
              write(25,9)'SO2 ',(xtg(mg,k,2),k=1,nl)
              write(25,9)'SO4 ',(xtg(mg,k,3),k=1,nl)
            endif
            qint=0.
            xint1=0.
            xint2=0.
            xint3=0.
            do k=1,nl
              rhodz=100.*dprf(mg,k)/grav
              qint=qint+(qtg(mg,k)+qlg(mg,k)+qfg(mg,k))*rhodz
              xint1=xint1+(xtg(mg,k,1))*rhodz
              xint2=xint2+(xtg(mg,k,2))*rhodz
              xint3=xint3+(xtg(mg,k,3))*rhodz
            enddo
            write(25,1)'integrated q (mm) = ',qint
            if(coupled_aero)then
              write(25,1)'integrated DMS (kg/m2) = ',xint1
              write(25,1)'integrated SO2 (kg/m2) = ',xint2
              write(25,1)'integrated SO4 (kg/m2) = ',xint3
            endif
            write(25,*)

 1          format(4(a,g12.5))
 9          format(a,30g10.3)
 91         format(a,30f10.3)
          endif

        endif

        return

c**************************************************************

       entry c_aero2(lg,tdt,mcmax,tg,v10m,sicef,ustar,eg,fg,rcondx,  !Inputs
     &               xtg,                 !In & Out
     &               xtem, so2dd, so4dd, sem)   !Outputs

        rho1(:)=100.*prf(:,1)/(rdry*ttg(:,1)) !density of air

        do mg=1,ln2
          sem(mg)=0.
          so2dd(mg)=0.
          so4dd(mg)=0.
          mc1(mg)=mc(mg,lg)  !skin reservoir depth passed back from surfa        
          wgmax(mg)=psfc(mg,lg) !field capacity of soil passed back from surfa 
        enddo

c Calculate sub-grid Vgust
c Mesoscale enhancement follows Redelsperger et al. (2000), J. Climate 13, 402-421.
c Equation numbers follow Fairall et al. 1996, JGR 101, 3747-3764.

C***        do mg=1,ln2
C***          rrate = 8640.*rcondx(mg)/tdt !Rainfall rate in cm/day (really should be conv. only)
C***
C***c Calculate convective scaling velocity (Eq.17) and gustiness velocity (Eq.16)
C***
C***          zi=600.
C***          beta=0.65
C***          Wstar3 = max ( 0., (( grav*zi/tg(mg) )
C***     &    * ( fg(mg)/(rho1(mg)*cp) + 0.61*tg(mg)*eg(mg)/(rho1(mg)*hl))))
C***          Vgust_free = beta*Wstar3**(1./3)
C***
C***c Calculate the Redelsperger-based Vgust_deep if deep convection is present, and then take
C***c the maximum of these. Note that Redelspreger gives two other parameterizations, based on
C***c the updraft or downdraught mass fluxes respectively.
C***
C***          Vgust_deep = (19.8*rrate**2/(1.5+rrate+rrate**2))**0.4
C***          Vgust = max(Vgust_free, Vgust_deep)
C***
C***c Calculate effective 10m wind (Eq. 15)
C***
C***          veff(mg) = sqrt ( v10m(mg)**2 + Vgust**2 )
C***
C***        enddo


c Need to invert vertical levels for ECHAM code... Don't you hate that?

        do k=1,nlp
          aphp1(:,nlp+1-k)=100.*prh(:,k)
        enddo

        do k=1,nl
          xtm1(:,nlp-k,:)=xtg(:,k,:)
        enddo

c Emission and dry deposition

        call xtemiss
     &   (tdt, lg, xtm1, rho1, tg, sicef, v10m, aphp1,              !Inputs
     &    land, tsigmf, 1.e-3*snowd, mcmax, mc1, wgmax, wg, ustar,  !Inputs
     &    xte, xtem, so2dd, so4dd, sem)                             !Outputs

c Add tendency to tracer mixing ratio (invert vertical levels)

        do k=1,nl
          xtg(:,k,:)=xtg(:,k,:)+xte(:,nlp-k,:)*tdt
        enddo

c Add in contribution to sulfur emission at higher levels
        
C***        do nt=1,3  !Three sulfur species
C***          do k=1,nl
C***            sem(:)=sem(:)+xte(:,nlp-k,nt)*100.*dprf(:,k)/grav
C***          enddo
C***        enddo

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon

            xint1=0.
            xint2=0.
            xint3=0.
            k=2
            rhodz=100.*dprf(mg,k)/grav
            rhodz1=100.*dprf(mg,1)/grav
            rhodz2=100.*dprf(mg,2)/grav
            write(25,*)'After xtemiss, before hvermtx'
            write(25,1)'sem*tdt ',tdt*sem(mg)
            write(25,9)'DMS ',(xtg(mg,k,1),k=1,nl)
            write(25,9)'SO2 ',(xtg(mg,k,2),k=1,nl)
            write(25,9)'SO4 ',(xtg(mg,k,3),k=1,nl)
            write(25,*)

            xint1=0.
            xint2=0.
            xint3=0.
            do k=1,nl
              rhodz=100.*dprf(mg,k)/grav
              xint1=xint1+(xtg(mg,k,1))*rhodz
              xint2=xint2+(xtg(mg,k,2))*rhodz
              xint3=xint3+(xtg(mg,k,3))*rhodz
            enddo
            write(25,1)'integrated DMS (kg/m2) = ',xint1
            write(25,1)'integrated SO2 (kg/m2) = ',xint2
            write(25,1)'integrated SO4 (kg/m2) = ',xint3
            write(25,*)

            xint1=xtem(mg,1)*tdt
            xint2=xtem(mg,2)*tdt
            xint3=xtem(mg,3)*tdt
            do k=1,nl
              rhodz=100.*dprf(mg,k)/grav
              xint1=xint1+xte(mg,nlp-k,1)*rhodz*tdt
              xint2=xint2+xte(mg,nlp-k,2)*rhodz*tdt
              xint3=xint3+xte(mg,nlp-k,3)*rhodz*tdt
            enddo
            write(25,1)'integrated DMS emission+dep (kg/m2) = ',xint1,
     &           ' surf emission (kg/m2) ',xtem(mg,1)*tdt,
     &           ' lev 1 emission (kg/m2)', xte(mg,nl,1)*rhodz1*tdt
            write(25,1)'integrated SO2 emission+dep (kg/m2) = ',xint2,
     &           ' surf emission (kg/m2) ',xtem(mg,2)*tdt,
     &           ' lev 2 emission (kg/m2) ', xte(mg,nl-1,2)*rhodz2*tdt,
     &           ' lev 1 emission+dep (kg/m2)', xte(mg,nl,2)*rhodz1*tdt
            write(25,1)'integrated SO4 emission+dep (kg/m2) = ',xint3,
     &           ' surf emission (kg/m2) ',xtem(mg,3)*tdt,
     &           ' lev 2 emission (kg/m2) ', xte(mg,nl-1,3)*rhodz2*tdt,
     &           ' lev 1 emission+dep (kg/m2)', xte(mg,nl,3)*rhodz1*tdt
          endif
        endif

        return

c**************************************************************
      entry c_aero3(lg,tdt, !Inputs
     &              xtg)  !In and out

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            xint1=0.
            xint2=0.
            xint3=0.
            do k=1,nl
              rhodz=100.*dprf(mg,k)/grav
              xint1=xint1+(xtg(mg,k,1))*rhodz
              xint2=xint2+(xtg(mg,k,2))*rhodz
              xint3=xint3+(xtg(mg,k,3))*rhodz
            enddo
            write(25,1)'integrated DMS (kg/m2) = ',xint1
            write(25,1)'integrated SO2 (kg/m2) = ',xint2
            write(25,1)'integrated SO4 (kg/m2) = ',xint3
            write(25,*)
          endif
        endif

c Decay of hydrophobic black and organic carbon into hydrophilic forms

        do k=1,nl
          xtm1(:,nlp-k,:)=xtg(:,k,:)
        enddo
        call xtsink (tdt, xtm1, !Inputs
     &               xte)       !Output
        do k=1,nl
          xtg(:,k,:)=xtg(:,k,:)+xte(:,nlp-k,:)*tdt
        enddo

        xte(:,:,:)=0.

        return

c**************************************************************
      entry c_aero4(lg,tdt,xtg,pxtg,  !Inputs
     &              so2wd,so4wd)      !Outputs

      so2wd(:)=0.
      so4wd(:)=0.

      xtg(:,:,:)=pxtg(:,:,:)

      do k=1,nl
        do mg=1,ln2
          so2wd(mg)=so2wd(mg)
     &         +(pxtg(mg,k,2)-xtg(mg,k,2))*100.*dprf(mg,k)
     &         /(grav*tdt)
          so4wd(mg)=so4wd(mg)
     &         +(pxtg(mg,k,3)-xtg(mg,k,3))*100.*dprf(mg,k)
     &         /(grav*tdt)
        enddo
      enddo

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          qint=0.
          xint1=0.
          xint2=0.
          xint3=0.
          do k=1,nl
            rhodz=100.*dprf(mg,k)/grav
            qint=qint+(qtg(mg,k)+qlg(mg,k)+qfg(mg,k))*rhodz
            xint1=xint1+(xtg(mg,k,1))*rhodz
            xint2=xint2+(xtg(mg,k,2))*rhodz
            xint3=xint3+(xtg(mg,k,3))*rhodz
          enddo
          write(25,1)'integrated DMS (kg/m2) = ',xint1
          write(25,1)'integrated SO2 (kg/m2) = ',xint2
          write(25,1)'integrated SO4 (kg/m2) = ',xint3
          write(25,*)
          write(25,1)'SO2 convwetdep (kg/m2) = ',so2wd(mg)
          write(25,1)'SO4 convwetdep (kg/m2) = ',so4wd(mg)
          write(25,*)
        endif
      endif

      return

c**************************************************************
      entry c_aero5(lg,tdt,xtg)

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            xint1=0.
            xint2=0.
            xint3=0.
            do k=1,nl
              rhodz=100.*dprf(mg,k)/grav
              xint1=xint1+(xtg(mg,k,1))*rhodz
              xint2=xint2+(xtg(mg,k,2))*rhodz
              xint3=xint3+(xtg(mg,k,3))*rhodz
            enddo
            write(25,1)'After cvmix...'
            write(25,1)'integrated DMS (kg/m2) = ',xint1
            write(25,1)'integrated SO2 (kg/m2) = ',xint2
            write(25,1)'integrated SO4 (kg/m2) = ',xint3
            write(25,*)
          endif
        endif

        return

c**************************************************************
      entry c_aero6(lg,tdt,xtg,so2wd,so4wd)

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            xint1=0.
            xint2=0.
            xint3=0.
            do k=1,nl
              rhodz=100.*dprf(mg,k)/grav
              xint1=xint1+(xtg(mg,k,1))*rhodz
              xint2=xint2+(xtg(mg,k,2))*rhodz
              xint3=xint3+(xtg(mg,k,3))*rhodz
            enddo
            write(25,1)'SO2 tot wetdep (kg/m2) = ',so2wd(mg)*tdt
            write(25,1)'SO4 tot wetdep (kg/m2) = ',so4wd(mg)*tdt
            write(25,1)'integrated DMS (kg/m2) = ',xint1
            write(25,1)'integrated SO2 (kg/m2) = ',xint2
            write(25,1)'integrated SO4 (kg/m2) = ',xint3
            write(25,*)
          endif
        endif

        return

c**************************************************************
      entry c_aero7(lg,tdt,xtg,pxtg,ompl)

        do nt=1,ntrac
          do k=1,nl
            do mg=1,ln2
              xtg(mg,k,nt)=xtg(mg,k,nt)*ompl(mg)+pxtg(mg,k,nt)*pl(mg)
            enddo
          enddo
        enddo

        return

c**************************************************************
      entry c_aero8(lg,tdt,xtg,so2wd,so4wd,so2oh,so2h2,so2o3
     &   ,so2dd,so4dd,sem,dmsoh,dmsn3,conwd)

        rho1(:)=100.*prf(:,1)/(rdry*ttg(:,1)) !density of air

        do k=1,nl
          do mg=1,ln2
c Factor 1.e6 for sulfur fields to convert to mg/m2, times factor 100 for P
            hdms(mg)=hdms(mg)+1.e8*xtg(mg,k,1)*dprf(mg,k)/grav
            hso2(mg)=hso2(mg)+1.e8*xtg(mg,k,2)*dprf(mg,k)/grav
            hso4(mg)=hso4(mg)+1.e8*xtg(mg,k,3)*dprf(mg,k)/grav
          enddo
        enddo
        do mg=1,ln2
          hso2dd(mg)=hso2dd(mg)+so2dd(mg)*1.e12  !In 10^-12 kg/m2/s
          hso4dd(mg)=hso4dd(mg)+so4dd(mg)*1.e12  ! ditto
          hso2wd(mg)=hso2wd(mg)+so2wd(mg)*1.e12  !In 10^-12 kg/m2/s
          hso4wd(mg)=hso4wd(mg)+so4wd(mg)*1.e12  ! ditto
          hscwd(mg)=hscwd(mg)+(conwd(mg,2)+conwd(mg,3))*1.e12  ! ditto
          hsem(mg)  =hsem(mg)  +sem(mg)  *1.e12  ! ditto
          hso2oh(mg)=hso2oh(mg)+so2oh(mg)*1.e12  !In 10^-12 kg/m2/s
          hso2h2(mg)=hso2h2(mg)+so2h2(mg)*1.e12  !In 10^-12 kg/m2/s
          hso2o3(mg)=hso2o3(mg)+so2o3(mg)*1.e12  !In 10^-12 kg/m2/s
          hdmsoh(mg)=hdmsoh(mg)+dmsoh(mg)*1.e12  !In 10^-12 kg/m2/s
          hdmsn3(mg)=hdmsn3(mg)+dmsn3(mg)*1.e12  !In 10^-12 kg/m2/s
          hdms1(mg)=hdms1(mg)+1.e9*xtg(mg,1,1)*rho1(mg) !In ugS/m3
c          hdms1(mg)=hdms1(mg)+xtg(mg,1,1)*1.e12*28.96/32. !In pptv
          hso21(mg)=hso21(mg)+1.e9*xtg(mg,1,2)*rho1(mg) !In ugS/m3
c          hso21(mg)=hso21(mg)+xtg(mg,1,2)*1.e12*28.96/32. !In pptv
          hso41(mg)=hso41(mg)+1.e9*xtg(mg,1,3)*rho1(mg) !In ugS/m3
        enddo

        return

c**************************************************************
      entry c_aero9(lg,tdt,xtg)

        do nt=1,ntrac
          do k=1,nl
            do mg=1,ln2
              bigxt(mg,k,nt,lg)=xtg(mg,k,nt)
            enddo
          enddo
        enddo

        return

      end
