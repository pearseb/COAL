c Purpose
c -------
c Defines the COMMON block /HOSING/, which contains parameters relevant to
c freshwater hosing.
c
c Contains
c --------
c HOSING_FLAG	If true, apply freshwater hosing
c HOSING_RATE	Freshwater hosing rate [Sv]
c
c History
c -------
c 2009 Aug 6	Steven Phipps	Original version

      logical hosing_flag
      real hosing_rate
      common /hosing/ hosing_flag, hosing_rate
