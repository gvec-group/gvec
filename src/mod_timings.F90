!===================================================================================================================================
! Copyright (c) 2014 Lorenz Hüdepohl
! License: LGPL-3.0-only
!
! simple wrappers for mapping perflib calls to the ftimings library
! ftimings is available at https://gitlab.mpcdf.mpg.de/loh/ftimings
!===================================================================================================================================

module timings
  use ftimings
  type(timer_t) :: timer
end module
