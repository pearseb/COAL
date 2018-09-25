c $Log: SWOCOM.f,v $
c Revision 1.3  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/12/17  15:31:32  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  92/04/15  11:13:39  mrd
c Initial revision
c 
c     common block swocom contains quantities outputted by the 
c     shortwave radiation code to the external module:
      real fsw(imax,lp1)  ! Net radiation (up-down) in CGS units at all
                          ! pressure levels
      real dfsw(imax,lp1) ! Downward radiation at all pressure levels
      real ufsw(imax,lp1) ! Upward radiation at all pressure levels
      real hsw(imax,l)    ! Shortwave heating rates in K/day for pressure
                          ! layers. 
      common / swocom / fsw, dfsw, ufsw, hsw
