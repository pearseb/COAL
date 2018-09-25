c Added IMPLICIT NONE statement, plus variable declarations for I, J, DIV,
c PRES, PRESINC, STRENG and ZP.
c SJP 2007/05/28
c
c $Log: cavit2.f,v $
c Revision 1.4  1996/10/24 01:02:30  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.3  1996/03/21  03:18:31  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.2  1993/12/17  15:31:45  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  92/06/16  11:54:43  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:  from common/dicegrid in COMDICE.f
c                 deltxu - east-west widths of u-grid boxes
c                 delty  - nl-south widths of h-grid boxes
c                 latha, lathb - h-grid)latitudinal indices bounding the
c                                extent of ice
c                 srfarea - surface area of grid boxes
c
c             from common/dicemask in COMDICE.f
c                 dmask - h-grid "divergence" mask
c
c     In/Out: from common/ocwind COMDICE.f
c                 u - eastward  velocity on u-grid
c                 v - northward velocity on u-grid 
c
c             from arguments
c                 div    - divergence of u,v field
c                 pres   - internal ice "pressure" on h-grid
c                 streng - internal ice strength on h-grid
c 
      subroutine cavit2 (div, pres, streng)

c Accounts for ice rheology using an alternate cavitating-fluid method
c (flato=.false.). In each iteration, the pressure field is accumulated
c by an amount proportional to the current convergence, but only up
c to a maximum determined by the local ice strength. No increment is
c applied where there is divergence. The amount of pressure increment
c is the value that would produce the corresponding amount of
c divergence in icefree with no other forcing.

c u      = eastward  velocity on u-grid (modified)
c v      = northward velocity on u-grid (modified)
c div    = divergence of u,v field (returned)
c pres   = internal ice "pressure" on h-grid (modified)
c streng = internal ice strength on h-grid (supplied)

      implicit none
 
      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension
     &  div(plon,plat),         pres(plon,plat),
     &  streng(plon,plat)

      dimension presinc(plon)

      integer i, j
      real div, pres, presinc, streng, zp

c Set the ring of velocities at 90N to pure slab motion

      call polefilt (u, v, 90, 4)


c Calculate the divergence, paralleling calculations in advect,
c i.e, using 2-point differences and using spherical metrics delt*

      do 100 j=latha,lathb
      do 100 i=1,plon
        div(i,j) = ( - 0.5 * (u(i,j)   + u(i,j+1))   * delty(j)
     &               + 0.5 * (u(i+1,j) + u(i+1,j+1)) * delty(j)
     &               - 0.5 * (v(i,j)   + v(i+1,j))   * deltxu(j)
     &               + 0.5 * (v(i,j+1) + v(i+1,j+1)) * deltxu(j+1)
     &             ) / srfarea(j)
  100 continue


c Increment the pressure where the flow is convergent, limited by
c the ice strength.

      do 200 j=latha,lathb
        do 202 i=1,plon
          if (dmask(i,j).gt.0.) then                           !masscomp
            presinc(i) = -min(div(i,j),0.) / dmask(i,j)        !masscomp
          else                                                 !masscomp
            presinc(i) = 0.                                    !masscomp
          endif                                                !masscomp
c         don't use epsilon here...dmask is O(1/srfarea), = O(epsilon)
cray      presinc(i)= cvmgt ( -min(div(i,j),0.)/max(dmask(i,j),1.e-50),
cray *                        0., dmask(i,j).gt.1.e-50 )

          pres(i,j) = min ( pres(i,j)+presinc(i), streng(i,j) )
  202   continue
  200 continue


c Smear pres at latitude nearest pole, to avoid instability

      do 300 j=1,plat,plat-1
        if (latha.le.j .and. j.le.lathb) then
          zp = 0.
          do 310 i=1,plon
            zp = zp + pres(i,j)
  310     continue
          zp = zp/plon
          do 312 i=1,plon
            pres(i,j) = zp
  312     continue
        endif
  300 continue

      return
      end
