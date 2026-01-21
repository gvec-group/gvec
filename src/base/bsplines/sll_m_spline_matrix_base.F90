! Copyright (c) INRIA
! License: CECILL-B
!
!> @ingroup splines
!> @brief
!> Abstract class for small matrix library with basic operations:
!> set matrix element, factorize, solve and write to output
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_matrix_base

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_working_precision, only: f64

  implicit none

  public :: sll_c_spline_matrix

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  type, abstract :: sll_c_spline_matrix
  contains
    procedure(i_sub_mat_copy     ), deferred :: mat_copy
    procedure(i_sub_mat_add      ), deferred :: mat_add
    procedure(i_sub_set_element  ), deferred :: set_element
    procedure(i_sub_set_element  ), deferred :: add_element
    procedure(i_sub_get_element  ), deferred :: get_element
    procedure(i_sub_matvec_prod  ), deferred :: matvec_prod
    procedure(i_sub_factorize    ), deferred :: factorize
    procedure(i_sub_solve_inplace), deferred :: solve_inplace
    procedure(i_sub_write        ), deferred :: write
    procedure(i_sub_free         ), deferred :: free
  end type sll_c_spline_matrix

  abstract interface

  !-----------------------------------------------------------------------------
  subroutine i_sub_mat_copy( self, tocopy )
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(inout) :: self
    class(sll_c_spline_matrix), intent(in   ) :: tocopy
  end subroutine i_sub_mat_copy

  !-----------------------------------------------------------------------------
  subroutine i_sub_mat_add( self, a,amat,b,bmat )
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(inout) :: self
    real(wp)                  , intent(in   ) :: a
    class(sll_c_spline_matrix), intent(in   ) :: amat
    real(wp)                  , intent(in   ) :: b
    class(sll_c_spline_matrix), intent(in   ) :: bmat
  end subroutine i_sub_mat_add

  !-----------------------------------------------------------------------------
  subroutine i_sub_set_element( self, i, j, a_ij )
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(inout) :: self
    integer                   , intent(in   ) :: i
    integer                   , intent(in   ) :: j
    real(wp)                  , intent(in   ) :: a_ij
  end subroutine i_sub_set_element

  !-----------------------------------------------------------------------------
  function i_sub_get_element( self, i, j) result( a_ij )
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(inout) :: self
    integer                   , intent(in   ) :: i
    integer                   , intent(in   ) :: j
    real(wp)                                  :: a_ij
  end function i_sub_get_element

  !-----------------------------------------------------------------------------
  function i_sub_matvec_prod( self, v_in) result(v_out)
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(in) :: self
    real(wp)                  , intent(in) :: v_in(:)
    real(wp)                               :: v_out(size(v_in))
  end function i_sub_matvec_prod

  !-----------------------------------------------------------------------------
  subroutine i_sub_factorize( self )
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(inout) :: self
  end subroutine i_sub_factorize

  !-----------------------------------------------------------------------------
  subroutine i_sub_solve_inplace( self, nrhs, bx )
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(in   ) :: self
    integer                   , intent(in   ) :: nrhs
    real(wp),dimension(*)     , intent(inout) :: bx
  end subroutine i_sub_solve_inplace

  !-----------------------------------------------------------------------------
  subroutine i_sub_write( self, unit, fmt )
   import sll_c_spline_matrix, wp
    class(sll_c_spline_matrix), intent(in) :: self
    integer         , optional, intent(in) :: unit
    character(len=*), optional, intent(in) :: fmt
  end subroutine i_sub_write

  !-----------------------------------------------------------------------------
  subroutine i_sub_free( self )
   import sll_c_spline_matrix
    class(sll_c_spline_matrix), intent(inout) :: self
  end subroutine i_sub_free

  end interface

end module sll_m_spline_matrix_base
