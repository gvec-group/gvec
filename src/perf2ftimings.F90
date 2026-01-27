!===================================================================================================================================
! Copyright (c) 2014 Lorenz Hüdepohl
! License: LGPL-3.0-only
!
! simple wrappers for mapping perflib calls to the ftimings library
! ftimings is available at https://gitlab.mpcdf.mpg.de/loh/ftimings
!===================================================================================================================================

subroutine perfinit
  use timings

  call timer%enable()
end subroutine perfinit

subroutine perfon(label)
  use timings
  character(*), intent(in) :: label

  call timer%start(trim(adjustl(label)))
end subroutine perfon

subroutine perfoff(label)
  use timings
  character(*), intent(in) :: label

  call timer%stop(trim(adjustl(label)))
end subroutine perfoff

subroutine perfout(label)
  use timings
  character(*), intent(in) :: label

  call timer%print()
end subroutine perfout

subroutine perf_context_start(label)
  character(*), intent(in) :: label

end subroutine perf_context_start

subroutine perf_context_end()

end subroutine perf_context_end
