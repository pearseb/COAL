c $Log: advect.f,v $
c Revision 1.11  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.10  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.9  1998/12/10  00:55:44  ldr
c HBG changes to V5-1-21
c
c Revision 1.8  1997/12/19  02:03:16  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.7  1996/10/24  01:02:27  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/06/13  02:05:32  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.5  1996/03/21  03:18:30  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.4  1993/12/17  15:31:39  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.3  93/03/16  16:44:02  ldr
c Removed redundant tstx from common block dicemask.
c 
c Revision 1.2  92/07/16  12:52:19  ldr
c Minor changes to dynamical seaice model to correct runs not starting in
c January, make general for R42 and remove need for 'mask' file.
c 
c Revision 1.1  92/06/16  11:54:42  ldr
c Initial revision
c 
c     INPUT/OUTPUT
c     Input:   from common/dicegrid in COMDICE.f
c                  deltxu - east-west widths of u-grid boxes
c                  delty  - nl-south widths of h-grid boxes
c                  latha,lathb - h-grid)latitudinal indices bounding
c                  srfarea - surface area of grid boxes
c
c              from common/icetime in COMDICE.f
c                  dtf - dynamic seaice time step, seconds
c
c              from common/ocwind in COMDICE.f
c                  u - eastward  velocity on u-grid
c                  v - northward velocity on u-grid 
c
c              from arguments
c                  factor - post-advection qty to divide advected field to 
c                           get result 
c                  icode  - see below
c
c     In/Out:  from common/dicework in COMDICE.f
c                  workh  - h grid array
c
c              from arguments
c                  field - h-grid field
c                  result - advected field/factor
c
c 
      subroutine advect (result,field,factor,icode,ihem,pres,streng)

c Advects an h-grid field "field" by u-grid velocities (u,v) through
c time dtf. Divides advected field by "factor" to give "result".

c The flux-form used here guarantees conservation of the advected 
c fields. Note that "upstream" values are advected - which guarantees
c that recovered thicknesses and temperatures will have reasonable
c values (as long as the CFL condition is satisfied so no fields become
c negative). Given the CFL condition, ice concentration (and other 
c fields) will approach but never reach zero in a continuously
c divergent grid box - the thermodynamic ice model can kill off small
c ice quantities if it wants.

c All arguments are supplied except "result", and on h-grid except u,v.
c result = final advected result (not necessarily the advected field)
c field  = pre-advection field to be advected
c factor = post-advection qty to divide advected field to get result
c u      = eastward  velocity (u-grid)
c v      = northward velocity (u-grid)
c dtf     = time step
c icode  =  0 for ice concentration,
c           1 for ice mass,
c           2 for ice heat,
c           3 for ice brine resevoir,
c          10 for snow concentration,
c          11 for snow mass,
c          12 for snow heat.

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'COMDICE.f' !

C Argument list
      real result(plon,plat)
      real field(plon,plat)
      real factor(plon,plat)
      integer icode
      integer ihem
      real pres(plon,plat)
      real streng(plon,plat)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real div(plon,plat)
      logical ncaradv

      integer i
      integer icnt
      integer j

      real areacom
      real areain
      real areainc
      real areanorm
      real avgcom
      real avgin
      real avginc
      real dudx
      real fb
      real fl
      real fr
      real ft
      real renorm
      real ul
      real ur
      real vb
      real vt

C Local data, functions etc

C Start code : ----------------------------------------------------------



c Transfer h-grid field to be advected to a wrapped h-grid array

      call ltoh (field, workh)

c for NCAR advect method use ncaradv=.true.
c     ncaradv=.false.
      ncaradv=.true.

      IF(ncaradv)THEN

*PDIR SERIAL
c Same as NCAR with div set to zero when pres < ice strength
      avgin = 0.
      areain  = 0.
      do j=latha,lathb
      do i=1,plon
         if (workh(i,j).gt.0.0) then
            avgin = avgin + workh(i,j)*srfarea(j)
            areain  = areain  + srfarea(j)
         endif
      enddo
      enddo
      if (areain.gt.0.0) then
          avgin = avgin/areain
      else
          avgin = 0.
      endif
c    Compute the divergence field on the h-grid for use in the
c     revised transport equations.
c ** Correction for converging ice:  if pressure is >0, < ice strength,
c      then zero out the divergence at that point
c
        do j=latha,lathb
        do i=1,plon
          ul = 0.5 * (u(i,j)   + u(i,j+1))
          ur = 0.5 * (u(i+1,j) + u(i+1,j+1))
          vb = 0.5 * (v(i,j)   + v(i+1,j))
          vt = 0.5 * (v(i,j+1) + v(i+1,j+1))

          div(i,j)=((ur-ul)*delty(j) +(vt*deltxu(j+1)-vb*deltxu(j)) )
     &              /srfarea(j)

          if(flato)then
          if (pres(i,j).gt.0.and.pres(i,j).lt.streng(i,j))
c         if (pres(i,j).gt.0.and.pres(i,j).eq.streng(i,j))
     &     div(i,j) =-1.0E-11

          if(kma(i,j).eq.7)div(i,j)=-1.e-11
          endif

        enddo
        enddo

c Do advection
c *****  Uses the advection equation of B. Briegleb that is reformulated
c     in terms of the divergence field, which is pre-computed in icedia.f
c   Very different from Pollard's flux-form advection, but gives the
c    identical answer.  For the 'incompressible' grid points, div=0, set
c    in icedia.f  
c
cjcc                                             10/26/94
      do 400 j=latha,lathb
c     do 400 j=latha,lathb-1
cjcc
      do 400 i=1,plon
        ul = 0.5 * (u(i,j)   + u(i,j+1))
        ur = 0.5 * (u(i+1,j) + u(i+1,j+1))
        vb = 0.5 * (v(i,j)   + v(i+1,j))
        vt = 0.5 * (v(i,j+1) + v(i+1,j+1))

c       fl = max(ul,0.)*workh(i-1,j) + min(ul,0.)*workh(i,j)
c       fr = max(ur,0.)*workh(i,j)   + min(ur,0.)*workh(i+1,j)
c       fb = max(vb,0.)*workh(i,j-1) + min(vb,0.)*workh(i,j)
c       ft = max(vt,0.)*workh(i,j)   + min(vt,0.)*workh(i,j+1)
c
c jww 
        fl = 0.
        if (ul .gt. 0. ) fl = workh(i-1,j)
        if (ul .lt. 0. ) fl = workh(i,j)

        fr = 0.
        if (ur .gt. 0. ) fr = workh(i,j)
        if (ur .lt. 0. ) fr = workh(i+1,j)

        fb = 0.
        if (vb .gt. 0. ) fb = workh(i,j-1)
        if (vb .lt. 0. ) fb = workh(i,j)

        ft = 0.
        if (vt .gt. 0. ) ft = workh(i,j)
        if (vt .lt. 0. ) ft = workh(i,j+1)

        dudx = (ur-ul)*delty(j)/srfarea(j)

        
        field(i,j) = workh(i,j) + dtf 
     &              *( -(2.*div(i,j)*((fl+fr+fb+ft)*.25) )
     &              +((div(i,j)-dudx)*((fl+fr)*.50))
     &              + (dudx*((fb+ft)*.50))
     &              +((ul+ur)*.50)*((fl-fr)*delty(j)/srfarea(j))
     &      +((vb*deltxu(j)+vt*deltxu(j+1))*.50)*((fb-ft)/srfarea(j)) )
c
c       field(i,j) = workh(i,j) + dt
c    &             * (  fl*delty(j)  - fr*delty(j)
c    &                + fb*deltxu(j) - ft*deltxu(j+1) ) / srfarea(j)
c jww
  400 continue

c     add pole point treatment
cjww -- use gradient form of advection equation, assuming div =0, and
c         du/dx = dv/dy = 0
c
c     if(ihem.eq.2) then
c     j=plat
c     sum=0.0
c     do i=1,plon
c       vb = 0.5 * (v(i,j)   + v(i+1,j))
c       fb = max(vb,0.)*workh(i,j-1) + min(vb,0.)*workh(i,j)
c       sum = sum  + fb*deltxu(j)
c     enddo
c
c     field(1,j) = workh(1,j) + dtf * sum /(srfarea(j)*plon)
c     do i=2,plon
c       field(i,j) = field(1,j)             
c     enddo
c     endif

cjcc
cjww
c

c    Calculate area average after advection
c
      if(flato)then
      avginc = 0.
      avgcom = 0.
      areainc = 0.
      areacom = 0.
      icnt=0
      do j=latha,lathb
      do i=1,plon
         if (field(i,j).gt.0.0) then
           if (div(i,j).eq.-1.0E-11) then
             avginc = avginc + field(i,j)*srfarea(j)
             areainc=areainc + srfarea(j)
            icnt=icnt+1
           else
             avgcom = avgcom + field(i,j)*srfarea(j)
             areacom=areacom + srfarea(j)
           endif
         endif
      enddo
      enddo
      if (areainc.gt.0.) then
           avginc = avginc/areainc
      else
           avginc = 0.0
      endif
      if (areacom.gt.0.0) then
          avgcom = avgcom/areacom
      else
          avgcom = 0.0
      endif
c     write(6,*)avginc,avgcom,areainc,areacom,avgin,areain
c
c Calculate renormalization factor
 
      if (avgcom.ne.0.0) then
        renorm =  (avgin*areain - avginc*areainc)/(avgcom*areacom)
      else
         renorm = 1.0
      endif
c Put upper limit on renorm
      if (areacom.gt.0.) then
          areanorm = areainc/areacom
      else
          areanorm = 100. 
      endif
      if (areanorm.ge.15.0) then
          if (avginc.ne.0.0) then
               renorm = avgin/avginc
          else
               renorm = 1.0
          endif
      endif
c
c  Renormalize field after advection
c
c special treatment for icode 3 as for many of the grid points this is zero 
c anyway as no heat is stored internally.
      if(icode.eq.3)renorm=1.0
      do j=latha,lathb
      do i=1,plon
           if (div(i,j).ne.-1.0E-11) then
             field(i,j) = renorm * field(i,j)
           endif
      enddo
      enddo

cjww
      endif
*PDIR ENDSERIAL

      ELSE ! IF(.not.ncaradv) 

c Do advection old way

      do 500 j=latha,lathb
 
      do 500 i=1,plon
        ul = 0.5 * (u(i,j)   + u(i,j+1))
        ur = 0.5 * (u(i+1,j) + u(i+1,j+1))
        vb = 0.5 * (v(i,j)   + v(i+1,j))
        vt = 0.5 * (v(i,j+1) + v(i+1,j+1))
   
        fl = max(ul,0.)*workh(i-1,j) + min(ul,0.)*workh(i,j)
        fr = max(ur,0.)*workh(i,j)   + min(ur,0.)*workh(i+1,j)
        fb = max(vb,0.)*workh(i,j-1) + min(vb,0.)*workh(i,j)
        ft = max(vt,0.)*workh(i,j)   + min(vt,0.)*workh(i,j+1)

        field(i,j) = workh(i,j) + dtf
     &             * (  fl*delty(j)  - fr*delty(j)
     &                + fb*deltxu(j) - ft*deltxu(j+1) ) / srfarea(j)
  500 continue

      ENDIF ! IF(ncaradv)


c Optionally smear latitude bands nearest poles (on h-grid, j=1 or plat)

c     do 450 j=1,plat,plat-1
c       if (latha.le.j .and. j.le.lathb) then
c         zf = 0.
c         do 452 i=1,plon
c           zf = zf + field(i,j)
c 452     continue
c         zf = zf/plon
c         do 454 i=1,plon
c           field(i,j) = zf
c 454     continue
c       endif
c 450 continue


c For the primary (first-advected) field of ice concentration, only
c set "result" (ie, allow creation in a previously ice-free box) if
c the new concentration is above a small limit, since arbitrarily small
c concentrations might cause machine-limit problems.

      if (icode.eq.0) then
        do 501 j=latha,lathb
        do 501 i=1,plon
          if (field(i,j).gt.epsilon) result(i,j) = field(i,j)
  501   continue

c For all other fields, divide by "factor" if it is non-zero.
c "factor" always involves the advected ice concentration, so will be
c zero if advection was cut off by the small limit in loop 500. Also
c it is always zero for land points, where we mustn't change "result"
c since that would corrupt soil thicknesses and land snow amounts.

      else
        do 510 j=latha,lathb
        do 510 i=1,plon
          if (factor(i,j).gt.0.) result(i,j) = field(i,j) / factor(i,j)
  510   continue
      endif


c For ice and snow concentrations, impose maximum of 1.

      if (icode.eq.0 .or. icode.eq.10) then
        do 600 j=latha,lathb
        do 600 i=1,plon
          result(i,j) = min (result(i,j), 1.)
  600   continue
      endif

       if(icode.eq.3)then
       do 700 j=latha,lathb
       do 700 i=1,plon
       if(result(i,j).lt.epsilon)result(i,j)=0.0
 700   continue
       endif

c Check - advection should never create negative values

c     do 900 j=latha,lathb                                    !masscomp
c     do 900 i=1,plon                                         !masscomp
c       if (field(i,j).lt.0.) write(ioterm,902) i,j,icode,    !masscomp
c    &  field(i,j),workh(i,j),factor(i,j),result(i,j)         !masscomp
c 902   format(' *** Warning. advect. i,j,icode,',            !masscomp
c    &         'field,workh,factor,result:',3i4,4e13.4)       !masscomp
c 900 continue                                                !masscomp

      return
      end
