c Defines the new COMMON block /GIANT5/. This contains the array PHNF,
c which was previously contained within the COMMON block /GIANT4/.
c SJP 2009/03/12

      real phnf
      common/giant5/phnf(ln2,nl)
