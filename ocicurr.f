# 1 "ocicurr.F"
c Minor enhancements to the coupling of the new, high-resolution version of the
c ocean model to the atmosphere.
c SJP 2007/12/21
c
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c Replaced the line "include 'OCEANPARS.f'" with "include 'OPARAMS.f'",
c enabling the header file OCEANPARS.f to be removed from the model source
c code.
c SJP 2007/05/31
c
c Transferred COMMON blocks to separate header files, as follows:
c /OCUR/  ->  OCUR.f
c SJP 2007/05/29
c
c     INPUT    from common/ocur in this subroutine
c                   uco, vco - ocean model instantaneous velocities
c
c     OUTPUT   from common/ocwind in COMDICE.f
c                  ocu - ocean currents, zonal
c                  ocv - ocean currents, meridional
c
      subroutine ocicurr
C
C Get ocean model surface velocities for ice model IF COUPLED
C
C If not coupled, then velocities are interpolated from
C  monthly mean data sets - see icedrive.f
C
      include 'PARAMS.f'
      include 'OPARAMS.f'
      include 'COMDICE.f'
      include 'FEWFLAGS.f'
      include 'ONEDIM.f'
C
C Need oceanic and ice models to have matching U grids for this code
C
CMOM1 R21 AGCM, R21 OGCM
      include 'O2A.f'
C
CMOM2 T63 AGCM, T63_2 OGCM
CMOM2  (ima=lon=plon,jma=lat2=plat)
CMOM2 common/CSIROocur/Cuco(ima+2,jma+2),Cvco(ima+2,jma+2)
C
      common/CSIROocur/Cuco(plon+2,plat+2),Cvco(plon+2,plat+2)

c.... Get the coupled ocean model velocities for the ice model:
c.... Put them into the arrays holding Atmos/Ice model ocean currents
c.... (i.e. From uco into ocu; and from vco into ocv).

c---- For MOM1 (R21 model):
c.... uco,vco hold ocean model instantaneous velocities
c.....  See step.f for MOM1 code.
c---- For MOM2 (T63 model):
c.... Cuco,Cvco hold ocean model instantaneous velocitie
c....   See gasbc.F for MOM2 code.

c---- Notes for MOM1 code
c.... The ocean array dimensions are (plon+2,plat+2)
c.... These are already set for i=1,plon+2 and j=2,plat+1 only with
c.... j=plat+1 (NP row) being zeroes.
c.... Because of a NE shift in the ocean UV grid relative to the ocean T
c.... grid and a SW shift in the ice UV grid relative to the ice T grid
c.... AND the fact that the expanded ocean grid is from (1:plon+2,1:plat+2)
c.... whereas the expanded ice grid is from (0:plon+1,1:plat+1)
c.... THEN the (i,j) index for the ocean U/V data is exactly the same
c.... spatial position as the same (i,j) index for the ice U/V point!
c.... (I think !#?#)
c.... See notes on the CSIRO model grids.

      if(lcouple)then

      if(plon.eq.64)then 
c---- MOM1 (Mk2) R21 currents

# 89


c...  When using the high-resolution version of the ocean model, the surface
c...  velocities for each ice model gridbox are obtained by calculating the
c...  area-weighted average of the surface velocities for the corresponding
c...  region of the OGCM grid.
c...
c...  The OGCM does not calculate velocities for the northernmost two rows of
c...  the grid, so we assume that the velocities at these latitudes are equal
c...  to zero.

        do j = 2, plat+1
          jo = 2 * j - 1
          jom1 = jo - 1
          jop1 = jo + 1

          do i = 1, plon+1
            io = 2 * i - 1
            iom1 = io - 1
            iop1 = io + 1
            if (iom1 .eq. 0) iom1 = imt - 2
          
            if (j .eq. plat+1) then

              ocu(i, j, 1) = (2.0 * cs(jom1) * uco(io, jom1) +
     &                        cs(jom1) * uco(iop1, jom1) +
     &                        cs(jom1) * uco(iom1, jom1)) /
     &                       (4.0 * cs(jom1) + 8.0 * cs(jo) +
     &                        4.0 * cs(jop1))
              ocv(i, j, 1) = (2.0 * cs(jom1) * vco(io, jom1) +
     &                        cs(jom1) * vco(iop1, jom1) +
     &                        cs(jom1) * vco(iom1, jom1)) /
     &                       (4.0 * cs(jom1) + 8.0 * cs(jo) +
     &                        4.0 * cs(jop1))

            else

              ocu(i, j, 1) = (4.0 * cs(jo) * uco(io, jo) +
     &                        2.0 * cs(jop1) * uco(io, jop1) +
     &                        2.0 * cs(jo) * uco(iop1, jo) +
     &                        2.0 * cs(jom1) * uco(io, jom1) +
     &                        2.0 * cs(jo) * uco(iom1, jo) +
     &                        cs(jop1) * uco(iop1, jop1) +
     &                        cs(jom1) * uco(iop1, jom1) +
     &                        cs(jom1) * uco(iom1, jom1) +
     &                        cs(jop1) * uco(iom1, jop1)) /
     &                       (4.0 * cs(jom1) + 8.0 * cs(jo) +
     &                        4.0 * cs(jop1))
              ocv(i, j, 1) = (4.0 * cs(jo) * vco(io, jo) +
     &                        2.0 * cs(jop1) * vco(io, jop1) +
     &                        2.0 * cs(jo) * vco(iop1, jo) +
     &                        2.0 * cs(jom1) * vco(io, jom1) +
     &                        2.0 * cs(jo) * vco(iom1, jo) +
     &                        cs(jop1) * vco(iop1, jop1) +
     &                        cs(jom1) * vco(iop1, jom1) +
     &                        cs(jom1) * vco(iom1, jom1) +
     &                        cs(jop1) * vco(iom1, jop1)) /
     &                       (4.0 * cs(jom1) + 8.0 * cs(jo) +
     &                        4.0 * cs(jop1))
            end if

            ocu(i, j, 1) = ocu(i, j, 1) * 0.01
            ocv(i, j, 1) = ocv(i, j, 1) * 0.01
            ocu(i, j, 2) = 0.0
            ocv(i, j, 2) = 0.0

          end do
        end do



      endif

      if(plon.eq.192)then 
c---- MOM2 (Mk3) T63 currents
        do j=2,plat+1
        do i=1,plon+1
          ocu(i,j,1)=Cuco(i,j)*.01
          ocu(i,j,2)=0.0
          ocv(i,j,1)=Cvco(i,j)*.01
          ocv(i,j,2)=0.0
        enddo
        enddo
      endif

      endif ! lcouple

      return
      end
