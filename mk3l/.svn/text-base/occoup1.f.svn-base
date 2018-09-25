c (1) Modified for the conversion of the coupled model restart file to netCDF.
c (2) Modified so that flux adjustments are only applied when FLUXADJ=T.
c SJP 2008/03/08
c
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Split off from ocean.f. LCOUPLE removed as an argument, as this is always
c true when this subroutine is called.
c SJP 2007/06/20

      subroutine occoup1
C
C=======================================================================
C The coupled model needs ocean/atmos coupling installed.            ===
C=======================================================================
C
      implicit none

      include 'PARAMS.f'
      include 'OPARAMS.f'
      include 'OCEAN_NML.f'
      include 'COEFFS.f'
      include 'FCOR.f'
      include 'DATICE1.f'   !agcm
      include 'ORESTART.f'
      include 'FEWFLAGS.f'

      integer i, j, mt, itt_rest

c......define some contants for the atmos and ocean......
      rhocw1=rhocw

c Andy set arck 0.01e-4 antk 0.05e-4 changed to current values to preven
c excessive ice thicknesses in early stages of coupled model on 29/1/93
      arck=0.03e-4
      antk=0.10e-4

c     zero out flux corrections (only needed if non-coupled)
      if (fluxadj) then
        do mt=1,2
          do j=1,jmt
            do i=1,imt
              hfcor(i,j,mt)=0.0
              sfcor(i,j,mt)=0.0
              txcor(i,j,mt)=0.0
              tycor(i,j,mt)=0.0
            enddo
          enddo
        enddo
      end if

c...  Read the surface fluxes from the coupled model restart file, unless the
c...  ocean model is being initialised from scratch
      if (nfirst .eq. 0) then
        write (*, *)
        write (*, *) "Reading data from coupled model restart file..."
        write (*, *)
        call oflux_read(itt_rest)
        if (itt_rest .ne. itt) then
          write (*, *)
          write (*, *)
     &      "***  WARNING: Timestep counter in coupled model restart"
          write (*, *)
     &      "***           file does not match ocean model timestep"
          write (*, *) "***"
          write (*, *) "***  Restart file:   ITT = ", itt_rest
          write (*, *) "***  Ocean model:    ITT = ", itt
          write (*, *)
        end if
      end if

      return
      end
