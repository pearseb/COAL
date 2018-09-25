c Modified by HBG for improved treatment of snowmelt, and of the case where
c the updated value of PL(MG) exceeds 1.0.
c 2005/02/24
c
c HFACA removed from /FICECON/, as this array is never used.
c SJP 2003/06/20
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: surfupl.f,v $
c Revision 1.44  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.43  2001/02/28 04:36:37  rot032
c Further tidy ups from HBG
c
c Revision 1.42  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.41  2000/11/21 01:05:37  rot032
c Correct setting of pl for full ice cover.
c
c Revision 1.40  2000/11/14 03:11:36  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.39  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.38  1998/12/10  00:55:36  ldr
c HBG changes to V5-1-21
c
c Revision 1.37  1997/12/23  00:23:35  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.36  1997/12/17  23:22:46  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.35  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.34.1.1  1997/12/19  02:03:12  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.34  1996/10/24  01:03:19  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.33  1996/03/21  03:19:10  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.32  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.28.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.31  1994/08/08  17:22:50  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.30  94/08/08  13:16:31  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.29  94/08/04  16:56:44  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.28  94/03/30  12:35:24  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.27  93/12/17  15:34:02  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.26  93/12/17  12:06:02  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.25  93/12/06  17:32:38  ldr
c Get rid of some irritating warnings on the VP.
c 
c Revision 1.24.1.1  93/12/17  11:51:51  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.24  93/11/29  16:09:44  ldr
c Corrected length of common block hm.
c 
c Revision 1.23  93/11/29  11:38:42  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
      subroutine surfupl(lg,z4,psg,prg,pblg,peg,pfg,
     &                   il,snowd,siced,ptg,qvent,
     &                   tstar,tg,snowmi)

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /RELATE/ )

c Physical variables prefixed by "p" are values over leads at cice points
c Those not prefixed by "p" are values over the ice at cice points.
     
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real z4(ln2)
      real psg(ln2)
      real prg(ln2)
      real pblg(ln2)
      real peg(ln2)
      real pfg(ln2)
      integer il(ln2)
      real snowd(ln2)
      real siced(ln2)
      real ptg(ln2)
      real qvent(ln2)
      real tstar(ln2)
      real tg(ln2)
      real snowmi(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)
      include 'RELATE.f'

C Global data blocks
      include 'DATICE1.f'
      include 'FEWFLAGS.f'
      include 'FREEZE.f'
      include 'QFLXDAT.f'
      include 'TIMEX.f'

      real flxia
      common/ficecon/flxia(lat,2)

      real hmo,hmo2
      common/hm/hmo(ln2,lat),hmo2(lon,lat2,2)

      real redice,groice,hav
      common/red/redice(ln2,lat),groice(ln2,lat),hav(ln2,lat)

C Local work arrays and variables
      real flx(ln2)
      real hcap1(ln2)
      real plmin(ln2)
      real smelt(ln2)

      real dicold
      real dz1
      real fint
      real hbot
      real hc
      real hi
      real hiinc
      real hinew
      real htop
      real htopx
      real pc
      real pdirad
      real plold
      real qheat
      real qhtice
      real sumflx
      real tfim
      real tmixo

      integer ma
      integer mg
      integer ns

C Local data, functions etc
      real qice,hcap,tfrz
      data qice,hcap,tfrz/3.35e8,2.095e8,273.15/    ! BP aug2010
      real qsnow
      data qsnow/3.35e8/
      real shi,shs
      data shi/1.8837e6/,shs/6.9069e5/

C Start code : ----------------------------------------------------------

c     routine to update the lead/ice ratios
c
c set up heat capacity of leads using leads mixed layer depth.
      if(.not.qflux)then
        dz1=100.0
        do 290 mg=1,ln2
  290   hcap1(mg)=hcap/50.*dz1
      elseif(qflux)then
        do 291 mg=1,ln2
  291   hcap1(mg)=hcap/50.*hmo(mg,lg)
      endif

c**** Sub-ice flux of heat from ocean 
      do 292 mg=1,lon
        flx(mg)=flxia(lg,1)+hoice(mg,lg)
        flx(mg+lon)=flxia(lg,2)+hoice(mg+lon,lg)
c set minimum values of pl for Arctic and Antarctic
        plmin(mg)=0.005
        plmin(mg+lon)=0.02
  292 continue

c compute zonal heat flux into the ice.
      do 2951 ns=1,2
        ma=(ns-1)*lon
        sumflx=0.0
        do 295 mg=1+ma,lon+ma
  295     if(il(mg).gt.0)sumflx=sumflx+flx(mg)
 2951   zflxi(ns)=sumflx

      do mg=1,ln2
       smelt(mg)=0.0
      enddo
 
c compute heating of leads => ice change and/or temp change
      do 300 mg=1,ln2
      if(il(mg).eq.0)goto 300
      ptg(mg)=ptg(mg)+z4(mg)*0.0065
c put heat flux into mixed layer in leads from ocean below (should increase
c rate of ice melt and slow down ice advance
c      flxi has been preset in ICECON
c---- note that Flxi is applied to both the ice covered part
c---- and leads part of the dynamic ice model grid point.

      dicold=dic(mg)
      plold=pl(mg)
      tmixo=tmix(mg)
      htop=psg(mg)-prg(mg)-peg(mg)-pfg(mg)+flx(mg)
      pdirad=4.0*pblg(mg)/ptg(mg)
      hbot=(hcap1(mg)/dt)+ pdirad
c     tnew=ptg(mg)+htop/hbot
c hiinc is the equivalent increment in ice depth from  the energy balance     
      hiinc=-hcap1(mg)/qice*(htop/hbot)      
c.... if hiinc < 0, then heating
c.... if hiinc > 0, then cooling.

      IF(hiinc.le.0.0)THEN
c....
c.... Heating is implied. A decrease in ice depth may occur,
c.... depending upon how warm the water is.
        if(tmix(mg).gt.272.15)then
c..     The leads water is warmer than 272.15. The heating is used to
c..     melt laterally, and also warm the leads water.
c..  50% goes into lateral melting 50% to mixed layer warming
c..  except if the temperature is > 273.15
          tmix(mg)=min(tmixo+0.5*htop/hbot,tfrz)    ! BP aug2010
c fint=50% except when tmix exceeds 273.15 then fint < 50 % and
c hiinc has to be more than 50 % (ie. 1-fint)
c (note :htop >= 0.0 here since hiinc <= 0)
          htopx=max(htop,0.1e-20)
          fint=min(1.0,max(0.0,(tmix(mg)-tmixo)*hbot/htopx))
c fint=0.5 except when tmix has been reset 273.15
          hiinc=hiinc*(1-fint)
        else  
c..     The leads water is cooler than 272.15. The heating is used to
c..     warm the leads water only.
          hiinc=0.0
          tmix(mg)=tmix(mg)+htop/hbot
        endif 

      ELSE
c.... (hiinc > 0)
c....
c.... Cooling of the leads water. An increase in ice depth may occur,
c.... depending upon how cold the water is.
        if(tmix(mg).gt.tfi)then
c..     The leads water is warmer than 271.3. The cooling reduces the 
c..     leads temp to tfi, then produces ice.
          tmix(mg)=max(tmixo+htop/hbot,tfi)
c (note : htop <=0 since hiinc >= 0)
          htopx=min(htop,-0.1e-20)
          fint=min(1.0,max(0.0,(tmix(mg)-tmixo)*hbot/htopx))
c fint will be 1.0 except when tmix(mg) tried to become colder
c than tfi when fint -> 0 
          hiinc=hiinc*(1.0-fint)
        endif
c..   (For leads temp < tfi, ice production only)

      END IF

c adjust the computed ice change to account for internal ice temps.
c qheat is amount of heat required to adjust newly frozen ice (at temp tfi)
c to the temperature of the ice (1 or 2 levels),
c or the amount of heat required to adjust the computed ice melt (at temp tmix)
c to the temperature of the ice (1 or 2 levels) and overlying snow.
c.... check for 1 or 2 levels of ice
      tfim=tfi
      if(hiinc.lt.0.0)tfim=tmix(mg)
      qhtice=tfim-tn1(mg)
      if(dic(mg).gt.0.2)qhtice=0.5*( (tfim-tn1(mg)) + (tfim-tn2(mg)) )
      qheat=dic(mg)*shi*qhtice

      hi=dic(mg)
      if(hiinc.gt.0.0)then
c.... creating only ice (no change to snow)
c.... at feezing temp tfi (tfim=tfi when hiinc>0)
c add in storage term if it exists
        qheat=qheat+sto(mg)
        hinew=hiinc*(qice*hi-qheat)/(hi*qice)
      else
c.... melting ice at temp of leads water (tfim=tmix)
c.... melting overlying snow cover at temp of formation tfi.
        if(dic(mg).gt.0.25)then
          qheat=qheat+ dsn(mg)*shs*(tfi-tn0(mg))
          hi=dic(mg)+0.1*dsn(mg)
        end if
        hinew=hiinc*(qice*hi-qheat)/(hi*qice)
      end if
c.... hi now contains (depth of ice) when creating ice by adding to thin ice
c....                         from below or by adding to thick ice laterally.
c....     or          (depth of ice) when melting thin ice from below.
c....     or          (depth of ice+snow) when melting thick ice laterally.

c     ptf(mg)=htop            
c     psf(mg)=hinew*1000./dt
 

c.... change either ice depth or ice extent
      if(dic(mg).lt.0.25)then
c if the the ice is thin then use hinew to change ice depth
        hc=hinew*pl(mg)/(1-pl(mg))
        pc=0.0
      else
c if the ice is thicker use hinew to change leads area
        hc=0.0
        pc=pl(mg)*hinew/hi
        if(pc.gt.0.0)then
c if pc>0 (ice cover has grown) then
c transfer volume of current snow across new ice
          dsn(mg)=dsn(mg)*(1.0-plold)/(1.0-(pl(mg)-pc))
        else
c if pc=0 (heating into leads), no ice change
c   and hiinc=hinew=pc=smelt=0.0
c if pc<0 (ice cover has melted) then
c melt the snow on the fraction (pc) of melted ice
c check pc does not make pl exceed max allowable value:
          if((pl(mg)-pc).gt.1.0)pc=pl(mg)-1.0
          smelt(mg)=-0.1*dsn(mg)*pc  ! > 0 (m of water)
          qvent(mg)=-sto(mg)*pc/dt 
        endif
      endif
      dic(mg)=dic(mg)+hc  
      groice(mg,lg)=groice(mg,lg)+hc
      pl(mg)=pl(mg)-pc

c set minimum values of pl for Arctic and Antarctic.
      if(pl(mg).lt.plmin(mg))then
        dic(mg)=dic(mg)*(1-pl(mg))/(1-plmin(mg))
        dsn(mg)=dsn(mg)*(1-pl(mg))/(1-plmin(mg))
        redice(mg,lg)=redice(mg,lg)+(dic(mg)-dicold)
        pl(mg)=plmin(mg)
      endif

c Check limits
      if((pl(mg).ge.1.0).or.(dic(mg).le.0.0))then
        smelt(mg)=smelt(mg)+
     &            (1.0-plold+pc)*max(dsn(mg),0.0)*0.1 ! m of water
        pl(mg)=1.0
        dsn(mg)=0.0
        dic(mg)=0.0
        tg(mg)=tfi
        tstar(mg)=tfi
c imsl and il switching occurs at next sea ice step
      endif

      snowd(mg)=dsn(mg)*100.
      siced(mg)=dic(mg)

c mix the 'leads' and 'under ice' water
      tmix(mg)=pl(mg)*tmix(mg)+(1-pl(mg))*tmixo

 300  continue

c Form snow melt on ice (and multiply by fractional area)
c snowmi() comes from seaice.f (>=0).
c Multiply it by the seaice fraction and add any snow melt from
c  leads expansion in this routine (already a fractional change)
      do mg=1,ln2
       if(il(mg).ne.0)then
        snowmi(mg)=(1.0-pl(mg))*snowmi(mg)+smelt(mg)
       endif
      enddo

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,3i3)')'After surfupl. IPASS=2.'
          write(25,1)'snowd ',snowd(mg),' siced ',siced(mg)
          write(25,1)'tg ',tg(mg),' pl ',pl(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,f7.2))

      return
      end
