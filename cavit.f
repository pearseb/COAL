c Re-implementing the "limited free slip" boundary condition at R21 and R42
c resolutions. The new implementation is simpler, and corrects the bugs that 
c were present in the previous version.
c SJP 2009/08/04
c
c Tidying up the source code for the sea ice model, particularly: fixing the
c bug in the calculation of the latitudes; and removing the buggy and 
c experimental code which sought to impose "limited slip at low resolution".
c SJP 2009/04/25
c
c Added IMPLICIT NONE statement, plus variable declarations for I, J, A,
c ASRAT, DIV, DU, DV, FIDIV, PRES, PRESINC, STRENG, UOLD, UWRAPSAV, VOLD and
c VWRAPSAV.
c SJP 2007/05/28
c
c Loops 100 and 200 merged for improved performance, as are the final three
c loops (limited slip).
c SJP 2003/04/24
c
c $Log: cavit.f,v $
c Revision 1.9  1998/12/10 00:55:31  ldr
c HBG changes to V5-1-21
c
c Revision 1.8  1997/12/19  02:03:10  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.7  1996/10/24  01:02:29  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/03/21  03:18:31  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.5  1995/07/26  07:28:38  ldr
c Merge Mickles speedups to ice scheme (V4-5-30mic) into V4-7-3l.
c
c Revision 1.4  1994/09/14  12:00:33  ldr
c Split do 300 loop to give Seca a fair go.
c
c Revision 1.3.1.1  1995/07/26  07:24:08  ldr
c Speedups from Mickles replacing character mask with integers.
c
c Revision 1.3  94/03/30  12:33:49  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.2  93/12/17  15:31:44  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.1  92/06/16  11:54:43  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:  from common/dicegrid in COMDICE.f
c                 deltxu - east-west widths of u-grid boxes
c                 delty  - nl-south widths of h-grid boxes
c                 latha, lathb - h-grid)latitudinal indices bounding
c                 latua, latub - u-grid)the extent of ice
c                 srfarea - surface area of grid boxes
c
c             from common/dicemask in COMDICE.f
c                 dmask - h-grid "divergence" mask
c                 umask - u-grid land-ocean mask 
c
c             from common/diceturn in COMDICE.f
c                 coswturn - water stress turning angle
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
      subroutine cavit (div, pres, streng)

c Accounts for ice rheology using Flato's cavitating-fluid method
c (flato=.true.). In each iteration, nothing is done for divergent
c h-grid points, and for convergent h-grid points a cancelling amount 
c of divergence is added by incrementing u,v at its corners. "Pressure"
c is accumulated each time this is done by an amount proportional to 
c the increment of divergence, up to a maximum pressure that is a 
c function of ice thickness and concentration. If the maximum pressure
c is exceeded, the flow is allowed to converge ("collapse").

c Since the increments in u and v affect convergence at previously
c adjusted h-grid points, the process must be iterated globally.
c Each iteration reduces kinetic energy but conserves momentum (at
c least for uniform hice and fice) and so convergence is guaranteed.

c u      = eastward  velocity on u-grid (modified)
c v      = northward velocity on u-grid (modified)
c div    = divergence of u,v field (returned)
c pres   = internal ice "pressure" on h-grid (modified)
c streng = internal ice strength on h-grid (supplied)

      implicit none

      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'COMDICE.f' !

      dimension
     &  div(plon,plat),         pres(plon,plat),
     &  streng(plon,plat)

      dimension
     &  presinc(plon),
     &  uwrapsav(plat+1),       vwrapsav(plat+1)

      dimension uold(0:plon+1,plat+1), vold(0:plon+1,plat+1)

      integer i, j
      real a, asrat, div, du, dv, pres, presinc, streng,
     &     uold, uwrapsav, vold, vwrapsav

c Calculate the divergence, paralleling calculations in advect,
c i.e, using 2-point differences and using spherical metrics delt*
c
c Increment the pressure where the flow is convergent. If the 
c accumulated pressure would exceed the maximum, reduce the increment
c of divergence so that the pressure becomes the maximum, and
c allow some convergence (ie, collapse).

      do 200 j=latha,lathb
        a = csubw * 0.5*(coswturn(j)+coswturn(j+1))
        do 202 i=1,plon
          div(i,j) = ( - 0.5 * (u(i,j)   + u(i,j+1))   * delty(j)
     &                 + 0.5 * (u(i+1,j) + u(i+1,j+1)) * delty(j)
     &                 - 0.5 * (v(i,j)   + v(i+1,j))   * deltxu(j)
     &                 + 0.5 * (v(i,j+1) + v(i+1,j+1)) * deltxu(j+1)
     &               ) / srfarea(j)

c         presinc #1: constant * mean incoming normal velocity,
c                 #2: constant * convergence.
c         Small differences in results with obstacle near pole, but
c         it's unclear which is more realistic.

c         presinc(i) = 2.*a*streslen * (-min(div(i,j),0.)) * srfarea(j)
c    *                 / (.5*delty(j)+.25*deltxu(j)+.25*deltxu(j+1))
          presinc(i) = 2.*a*(streslen**2) * (-min(div(i,j),0.))

          div(i,j) = min(div(i,j),0.) * (streng(i,j)-pres(i,j))
     &               / max (streng(i,j)-pres(i,j), presinc(i), epsilon)
          workp(i,j) = -min(div(i,j),0.)
          pres(i,j) = min (pres(i,j)+presinc(i), streng(i,j))
  202   continue
  200 continue

c Save u and v values at plon+1, for wraparound in loop 400 below

      do 250 j=latua,latub
        uwrapsav(j) = u(plon+1,j)
        vwrapsav(j) = v(plon+1,j)
  250 continue

c At this point div is zero at divergent points, and is equal either
c to the (-) convergence at convergent points, or to a reduced
c (-) convergence at those points where collapse has occured in this
c iteration. (It is zero at previously collapsed points). So for
c non-zero (-ve) div points, add in a cancelling amount of divergence
c by incrementing the u-grid velocities at the corners.

c This calculation of the u,v increments is an unfolding of the
c divergence calculation. It makes use of dmask, the divergence value
c corresponding to unit outgoing velocities at all non-masked corners,
c except with the u-velocities weighted by asrat, the mean aspect ratio
c of the grid box. (This weighting allows for non-square box dimensions
c assuming the divergence is applied isotropically, and significantly
c improves convergence with obstacle near pole.) Note dmask is zero
c only for land or for one-grid-cell-wide straits.

      do 300 j=latha,lathb
        asrat = 0.5*(deltxu(j)+deltxu(j+1)) / delty(j)
        do 302 i=1,plon
          dv = 0.                                              !masscomp
          if (dmask(i,j).gt.0.) dv = -div(i,j)/dmask(i,j)      !masscomp
          du = dv*asrat
          u(i+1,j)   = u(i+1,j)   + du*umask(i+1,j)
          u(i,j)     = u(i,j)     - du*umask(i,j)
          u(i+1,j+1) = u(i+1,j+1) + du*umask(i+1,j+1)
          u(i,j+1)   = u(i,j+1)   - du*umask(i,j+1)
          v(i+1,j)   = v(i+1,j)   - dv*umask(i+1,j)
          v(i,j)     = v(i,j)     - dv*umask(i,j)
          v(i+1,j+1) = v(i+1,j+1) + dv*umask(i+1,j+1)
          v(i,j+1)   = v(i,j+1)   + dv*umask(i,j+1)
  302   continue
  300 continue

c Allow for the east-west boundary, since the westernmost u-grid points
c (1,j) have not yet been influenced by the easternmost h-grid boxes 
c (plon,j). [u,v]wrapsav hold the previous iteration's u,v at plon+1.

      do 400 j=latua,latub
        u(1,j) = u(1,j) + u(plon+1,j) - uwrapsav(j)
        v(1,j) = v(1,j) + v(plon+1,j) - vwrapsav(j)

c       also set wraparounds

        u(plon+1,j) = u(1,j)
        v(plon+1,j) = v(1,j)

        u(0,j) = u(plon,j)
        v(0,j) = v(plon,j)
  400 continue

c...  Apply a limited free slip boundary condition at low resolution. The
c...  velocities at coastal gridpoints, which would otherwise be equal to zero,
c...  are set equal to the average of the velocities at the surrounding ocean
c...  gridpoints.

      if (mw .lt. 64) then ! limited slip at low resolution

c...  Save a copy of the velocity field
        do j = 1, plat+1
          do i = 0, plon+1
            uold(i, j) = u(i, j)
            vold(i, j) = v(i, j)
          end do
        end do

c...  Set the velocities at coastal gridpoints equal to the average of the
c...  velocities at the surrounding ocean gridpoints. The average velocities
c...  are calculated by applying a two-dimensional boxcar smoother of width
c...  three in both the x- and y-directions. This is implemented by calculating
c...  the sum of the velocities at the gridpoints (i-1:i+1, j-1:j+1) and then
c...  dividing by 9. The velocity at gridpoint (i, j) does not need to be
c...  included when calculating the sum, as this gridpoint must be land on the
c...  u-grid if CMASK is true.
        do j = latua, latub
          do i = 1, plon
            if (cmask(i, j)) then
              u(i, j) = (umask(i, j+1) * uold(i, j+1) +
     &                   umask(i+1, j+1) * uold(i+1, j+1) +
     &                   umask(i+1, j) * uold(i+1, j) +
     &                   umask(i+1, j-1) * uold(i+1, j-1) +
     &                   umask(i, j-1) * uold(i, j-1) +
     &                   umask(i-1, j-1) * uold(i-1, j-1) +
     &                   umask(i-1, j) * uold(i-1, j) +
     &                   umask(i-1, j+1) * uold(i-1, j+1)) / 9.0
              v(i, j) = (umask(i, j+1) * vold(i, j+1) + 
     &                   umask(i+1, j+1) * vold(i+1, j+1) +
     &                   umask(i+1, j) * vold(i+1, j) + 
     &                   umask(i+1, j-1) * vold(i+1, j-1) +
     &                   umask(i, j-1) * vold(i, j-1) + 
     &                   umask(i-1, j-1) * vold(i-1, j-1) +
     &                   umask(i-1, j) * vold(i-1, j) + 
     &                   umask(i-1, j+1) * vold(i-1, j+1)) / 9.0
            end if
          end do
          u(plon+1, j) = u(1, j)
          v(plon+1, j) = v(1, j)
          u(0, j) = u(plon, j)
          v(0, j) = v(plon, j)
        end do

      end if ! limited slip at low resolution

      return
      end
