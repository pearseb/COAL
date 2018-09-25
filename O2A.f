c Defines the COMMON block /O2A/, which contains the coupling fields passed
c from the ocean to the atmosphere.
c SJP 2007/12/19

      real osst, osal, uco, vco
      common /o2a/ osst(imt, jmt, 2), osal(imt, jmt, 2),
     &             uco(imt, jmt), vco(imt, jmt)

