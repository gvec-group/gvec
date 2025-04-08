! Copyright (c) INRIA
! License: CECILL-B
!
!> @ingroup splines
!> @brief   Derived type for dense matrices
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_matrix_dense

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_matrix_base, only: &
    sll_c_spline_matrix

  use iso_fortran_env, only: output_unit

  implicit none

  public :: sll_t_spline_matrix_dense

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  type, extends( sll_c_spline_matrix ) :: sll_t_spline_matrix_dense

    integer :: n
    logical :: factorized
    integer , allocatable :: ipiv(:)
    real(wp), allocatable :: a(:,:)

  contains

    procedure :: init          => s_spline_matrix_dense__init
    procedure :: mat_copy      => s_spline_matrix_dense__mat_copy
    procedure :: mat_add       => s_spline_matrix_dense__mat_add
    procedure :: set_element   => s_spline_matrix_dense__set_element
    procedure :: add_element   => s_spline_matrix_dense__add_element
    procedure :: get_element   => s_spline_matrix_dense__get_element
    procedure :: matvec_prod   => s_spline_matrix_dense__matvec_prod
    procedure :: factorize     => s_spline_matrix_dense__factorize
    procedure :: solve_inplace => s_spline_matrix_dense__solve_inplace
    procedure :: write         => s_spline_matrix_dense__write
    procedure :: free          => s_spline_matrix_dense__free

  end type sll_t_spline_matrix_dense

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__init( self, n )
    class(sll_t_spline_matrix_dense), intent(  out) :: self
    integer                         , intent(in   ) :: n

    SLL_ASSERT( n > 0 )

    self%n = n
    allocate( self%ipiv(n) )
    allocate( self%a(n,n) )
    self%a(:,:) = 0.0_wp
    self%factorized=.FALSE.

  end subroutine s_spline_matrix_dense__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__mat_copy( self,tocopy)
    class(sll_t_spline_matrix_dense), intent(inout) :: self
    class(sll_c_spline_matrix      ), intent(in   ) :: tocopy

    select type(tocopy); type is(sll_t_spline_matrix_dense)
    SLL_ASSERT( tocopy%n  == self%n  )

    self%n          = tocopy%n
    self%a(:,:)     = tocopy%a(:,:)
    self%ipiv(:)    = tocopy%ipiv(:)
    self%factorized = tocopy%factorized
    end select
  end subroutine s_spline_matrix_dense__mat_copy

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__mat_add( self,a,amat,b,bmat) !self=a*amat+b*bmat
    class(sll_t_spline_matrix_dense), intent(inout) :: self
    real(wp)                        , intent(in   ) :: a
    class(sll_c_spline_matrix      ), intent(in   ) :: amat
    real(wp)                        , intent(in   ) :: b
    class(sll_c_spline_matrix      ), intent(in   ) :: bmat

    select type(amat); type is(sll_t_spline_matrix_dense)
    select type(bmat); type is(sll_t_spline_matrix_dense)
    SLL_ASSERT( amat%n == self%n )
    SLL_ASSERT( bmat%n == self%n )
    SLL_ASSERT( .NOT.amat%factorized )
    SLL_ASSERT( .NOT.bmat%factorized )

    self%a(:,:) = a*amat%a(:,:)+b*bmat%a(:,:)
    self%ipiv(:)=0
    self%factorized=.FALSE.
    end select
    end select
  end subroutine s_spline_matrix_dense__mat_add

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__set_element( self, i, j, a_ij )
    class(sll_t_spline_matrix_dense), intent(inout) :: self
    integer                         , intent(in   ) :: i
    integer                         , intent(in   ) :: j
    real(wp)                        , intent(in   ) :: a_ij

    self%a(i,j) = a_ij

  end subroutine s_spline_matrix_dense__set_element

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__add_element( self, i, j, a_ij )
    class(sll_t_spline_matrix_dense), intent(inout) :: self
    integer                         , intent(in   ) :: i
    integer                         , intent(in   ) :: j
    real(wp)                        , intent(in   ) :: a_ij

    self%a(i,j) = self%a(i,j)+a_ij

  end subroutine s_spline_matrix_dense__add_element

  !-----------------------------------------------------------------------------
  function s_spline_matrix_dense__get_element( self, i, j ) result( a_ij )
    class(sll_t_spline_matrix_dense), intent(inout) :: self
    integer                         , intent(in   ) :: i
    integer                         , intent(in   ) :: j
    real(wp)                                        :: a_ij

    a_ij=self%a(i,j)

  end function s_spline_matrix_dense__get_element

  !-----------------------------------------------------------------------------
  function s_spline_matrix_dense__matvec_prod( self, v_in) result(v_out )
    class(sll_t_spline_matrix_dense), intent(in) :: self
    real(wp)                        , intent(in) :: v_in(:)
    real(wp)                                     :: v_out(size(v_in))
    integer                                      :: j

    SLL_ASSERT( size(v_in,1) == self%n )

    DO j=1,self%n
      v_out(j)=DOT_PRODUCT(self%a(:,j),v_in(:))
    END DO

  end function s_spline_matrix_dense__matvec_prod

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__factorize( self )
    class(sll_t_spline_matrix_dense), intent(inout) :: self

    integer :: info

    character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_dense % factorize"
    character(len=256) :: err_msg

    SLL_ASSERT( size(self%a,1) == self%n )
    SLL_ASSERT( size(self%a,2) == self%n )
    SLL_ASSERT( .not.self%factorized )

    ! Perform LU decomposition using Lapack (A=PLU)
    call dgetrf( self%n, self%n, self%a, self%n, self%ipiv, info )

    if ( info < 0 ) then
      write( err_msg, '(i0,a)' ) abs(info), "-th argument had an illegal value"
      SLL_ERROR(this_sub_name,err_msg)
    else if ( info > 0 ) then
      write( err_msg, "('U(',i0,',',i0,')',a)" ) info, info, &
           " is exactly zero. The factorization has been completed, but the factor" &
           //" U is exactly singular, and division by zero will occur if it is used to" &
           //" solve a system of equations."
      SLL_ERROR(this_sub_name,err_msg)
    end if
    self%factorized=.TRUE.

  end subroutine s_spline_matrix_dense__factorize

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__solve_inplace( self, nrhs,bx )
    class(sll_t_spline_matrix_dense), intent(in   ) :: self
    integer                         , intent(in   ) :: nrhs
    real(wp),dimension(*)           , intent(inout) :: bx

    integer :: info

    character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_dense % solve_inplace"
    character(len=256) :: err_msg

    SLL_ASSERT( size(self%a,1) == self%n )
    SLL_ASSERT( size(self%a,2) == self%n )
!    SLL_ASSERT( size(bx)  == self%n*nrhs )
    SLL_ASSERT( self%factorized )

    ! Solve linear system PLU*x=b using Lapack
    call dgetrs( 'N', self%n, nrhs, self%a, self%n, self%ipiv, bx, self%n, info )

    if ( info < 0 ) then
      write( err_msg, '(i0,a)' ) abs(info), "-th argument had an illegal value"
      SLL_ERROR(this_sub_name,err_msg)
    end if

    end subroutine s_spline_matrix_dense__solve_inplace

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__write( self, unit, fmt )
    class(sll_t_spline_matrix_dense), intent(in) :: self
    integer         , optional      , intent(in) :: unit
    character(len=*), optional      , intent(in) :: fmt

    integer :: i
    integer :: unit_loc
    character(len=32) :: fmt_loc

    if ( present( unit ) ) then
      unit_loc = unit
    else
      unit_loc = output_unit
    end if

    if ( present( fmt  ) ) then
      fmt_loc = fmt
    else
      fmt_loc = 'es9.1'
    end if

    write(fmt_loc,'(a)') "('(',i0,'" // trim(fmt_loc) // ")')"
    write(fmt_loc,fmt_loc) size(self%a,2)

    write(unit_loc,*) 'factorized?=',self%factorized
    do i = 1, size(self%a,1)
      write(unit_loc,fmt_loc) self%a(i,:)
    end do

  end subroutine s_spline_matrix_dense__write

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__free( self )
    class(sll_t_spline_matrix_dense), intent(inout) :: self

    self%n = -1
    self%factorized = .FALSE.
    if ( allocated( self%ipiv ) ) deallocate( self%ipiv )
    if ( allocated( self%a    ) ) deallocate( self%a    )

  end subroutine s_spline_matrix_dense__free

end module sll_m_spline_matrix_dense
