! Copyright (c) INRIA
! License: CECILL-B
!
!> @ingroup splines
!> @brief   Derived type for banded matrices
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_matrix_banded

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_matrix_base, only: &
    sll_c_spline_matrix

  use iso_fortran_env, only: output_unit

  implicit none

  public :: sll_t_spline_matrix_banded

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  type, extends( sll_c_spline_matrix ) :: sll_t_spline_matrix_banded

    integer :: n
    integer :: kl
    integer :: ku
    logical :: factorized
    integer, allocatable :: ipiv(:)
    real(wp), allocatable :: q(:,:)

  contains

    procedure :: init          => s_spline_matrix_banded__init
    procedure :: reset         => s_spline_matrix_banded__reset
    procedure :: mat_copy      => s_spline_matrix_banded__mat_copy
    procedure :: mat_add       => s_spline_matrix_banded__mat_add
    procedure :: set_element   => s_spline_matrix_banded__set_element
    procedure :: add_element   => s_spline_matrix_banded__add_element
    procedure :: get_element   => s_spline_matrix_banded__get_element
    procedure :: matvec_prod   => s_spline_matrix_banded__matvec_prod
    procedure :: factorize     => s_spline_matrix_banded__factorize
    procedure :: solve_inplace => s_spline_matrix_banded__solve_inplace
    procedure :: write         => s_spline_matrix_banded__write
    procedure :: free          => s_spline_matrix_banded__free

  end type sll_t_spline_matrix_banded

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__init( self, n, kl, ku )
    class(sll_t_spline_matrix_banded), intent(  out) :: self
    integer                          , intent(in   ) :: n
    integer                          , intent(in   ) :: kl
    integer                          , intent(in   ) :: ku

    SLL_ASSERT( n  >  0 )
    SLL_ASSERT( kl >= 0 )
    SLL_ASSERT( ku >= 0 )
    SLL_ASSERT( kl < n  )
    SLL_ASSERT( ku < n  )

    self%n  = n
    self%kl = kl
    self%ku = ku

    ! Given the linear system A*x=b, we assume that A is a square (n by n)
    ! with ku super-diagonals and kl sub-diagonals.
    ! All non-zero elements of A are stored in the rectangular matrix q, using
    ! the format required by DGBTRF (LAPACK): diagonals of A are rows of q.
    ! q has 2*kl rows for the subdiagonals, 1 row for the diagonal, and ku rows
    ! for the superdiagonals. (The kl additional rows are needed for pivoting.)
    ! The term A(i,j) of the full matrix is stored in q(i-j+2*kl+1,j).

    allocate( self%ipiv(n) )
    allocate( self%q(2*kl+ku+1,n) )
    self%q(:,:) = 0.0_wp
    self%ipiv(:) = 0
    self%factorized=.FALSE.

  end subroutine s_spline_matrix_banded__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__reset( self)
    class(sll_t_spline_matrix_banded), intent(inout) :: self

    self%q(:,:) = 0.0_wp
    self%ipiv(:) = 0
    self%factorized=.FALSE.

  end subroutine s_spline_matrix_banded__reset

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__mat_copy( self,tocopy)
    class(sll_t_spline_matrix_banded), intent(inout) :: self
    class(sll_c_spline_matrix       ), intent(in   ) :: tocopy

    select type(tocopy); type is(sll_t_spline_matrix_banded)
    SLL_ASSERT( tocopy%n  == self%n  )
    SLL_ASSERT( tocopy%kl == self%kl )
    SLL_ASSERT( tocopy%ku == self%ku )

    self%n         = tocopy%n
    self%kl        = tocopy%kl
    self%ku        = tocopy%ku
    self%q(:,:)    = tocopy%q(:,:)
    self%ipiv(:)   = tocopy%ipiv(:)
    self%factorized= tocopy%factorized
    end select

  end subroutine s_spline_matrix_banded__mat_copy

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__mat_add( self,a,amat,b,bmat) !self=a*amat+b*bmat
    class(sll_t_spline_matrix_banded), intent(inout) :: self
    real(wp)                         , intent(in   ) :: a
    class(sll_c_spline_matrix       ), intent(in   ) :: amat
    real(wp)                         , intent(in   ) :: b
    class(sll_c_spline_matrix       ), intent(in   ) :: bmat

    select type(amat); type is(sll_t_spline_matrix_banded)
    select type(bmat); type is(sll_t_spline_matrix_banded)
    SLL_ASSERT( amat%n  == self%n  )
    SLL_ASSERT( amat%kl == self%kl )
    SLL_ASSERT( amat%ku == self%ku )
    SLL_ASSERT( bmat%n  == self%n  )
    SLL_ASSERT( bmat%kl == self%kl )
    SLL_ASSERT( bmat%ku == self%ku )
    SLL_ASSERT( .not.amat%factorized )
    SLL_ASSERT( .not.bmat%factorized )

    self%q(:,:) = a*amat%q(:,:)+b*bmat%q(:,:)
    self%ipiv(:)=0
    self%factorized=.FALSE.
    end select
    end select
  end subroutine s_spline_matrix_banded__mat_add

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__set_element( self, i, j, a_ij )
    class(sll_t_spline_matrix_banded), intent(inout) :: self
    integer                          , intent(in   ) :: i
    integer                          , intent(in   ) :: j
    real(wp)                         , intent(in   ) :: a_ij

    SLL_ASSERT( max( 1, j-self%ku ) <= i )
    SLL_ASSERT( i <= min( self%n, j+self%kl ) )

    self%q(self%kl+self%ku+1+i-j,j) = a_ij

  end subroutine s_spline_matrix_banded__set_element

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__add_element( self, i, j, a_ij )
    class(sll_t_spline_matrix_banded), intent(inout) :: self
    integer                          , intent(in   ) :: i
    integer                          , intent(in   ) :: j
    real(wp)                         , intent(in   ) :: a_ij

    SLL_ASSERT( max( 1, j-self%ku ) <= i )
    SLL_ASSERT( i <= min( self%n, j+self%kl ) )

    self%q(self%kl+self%ku+1+i-j,j) = self%q(self%kl+self%ku+1+i-j,j)+ a_ij

  end subroutine s_spline_matrix_banded__add_element

  !-----------------------------------------------------------------------------
  function s_spline_matrix_banded__get_element( self, i, j ) result( a_ij )
    class(sll_t_spline_matrix_banded), intent(inout) :: self
    integer                          , intent(in   ) :: i
    integer                          , intent(in   ) :: j
    real(wp)                                         :: a_ij

    SLL_ASSERT( max( 1, j-self%ku ) <= i )
    SLL_ASSERT( i <= min( self%n, j+self%kl ) )

    a_ij = self%q(self%kl+self%ku+1+i-j,j)

  end function s_spline_matrix_banded__get_element

  !-----------------------------------------------------------------------------
  function s_spline_matrix_banded__matvec_prod( self, v_in) result(v_out )
    class(sll_t_spline_matrix_banded), intent(in) :: self
    real(wp)                         , intent(in) :: v_in(:)
    real(wp)                                      :: v_out(size(v_in))
    integer                                       :: j,imin,imax

    SLL_ASSERT( size(v_in,1) == self%n )
    SLL_ASSERT( .not.self%factorized   )

    DO j=1,self%n
      imin=max(1,j-self%ku)
      imax=min(self%n,j+self%kl)
      v_out(j)=DOT_PRODUCT(self%q(self%kl+self%ku+1+imin-j:self%kl+self%ku+1+imax-j,j),v_in(imin:imax))
    END DO

  end function s_spline_matrix_banded__matvec_prod

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__factorize( self )
    class(sll_t_spline_matrix_banded), intent(inout) :: self

    integer :: info

    character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_banded % factorize"
    character(len=256) :: err_msg

    SLL_ASSERT( .not.self%factorized   )
    ! Perform LU decomposition of matrix q with Lapack
    call dgbtrf( self%n, self%n, self%kl, self%ku, self%q, 2*self%kl+self%ku+1, &
                 self%ipiv, info )

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

  end subroutine s_spline_matrix_banded__factorize

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__solve_inplace( self, nrhs, bx )
    class(sll_t_spline_matrix_banded), intent(in   ) :: self
    integer                          , intent(in   ) :: nrhs
    real(wp),dimension(*)            , intent(inout) :: bx

    integer :: info

    character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_banded % solve_inplace"
    character(len=256) :: err_msg

!    SLL_ASSERT( size(bx)  == self%n*nrhs )
    SLL_ASSERT( self%factorized   )

    call dgbtrs( 'N', self%n, self%kl, self%ku, nrhs, self%q, 2*self%kl+self%ku+1, &
                 self%ipiv, bx, self%n, info )

    if ( info < 0 ) then
      write( err_msg, '(i0,a)' ) abs(info), "-th argument had an illegal value"
      SLL_ERROR(this_sub_name,err_msg)
    end if

  end subroutine s_spline_matrix_banded__solve_inplace

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__write( self, unit, fmt )
    class(sll_t_spline_matrix_banded), intent(in) :: self
    integer         , optional       , intent(in) :: unit
    character(len=*), optional       , intent(in) :: fmt

    integer :: i, j
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
      fmt_loc = 'es10.2'
    end if

    write(unit_loc,'(a6,a6,1X,A2,1X,a6,2X,a)') &
          "i","jmin","..","jmax","values..."
    write(fmt_loc,'(a)') "(" // trim(fmt_loc) // ")"

    write(unit_loc,*) 'factorized?=',self%factorized
    do i = 1, self%n
      write(unit_loc,'(i6,i6,1X,A2,1X,i6,2X)',advance='no') &
           i,max( 1, i-self%ku ),'..',min( self%n, i+self%kl )
      do j = max( 1, i-self%ku ), min( self%n, i+self%kl )
        write(unit_loc,fmt_loc,advance='no') self%q(self%kl+self%ku+1+i-j,j)
      end do
      write(unit_loc,*)
    end do

  end subroutine s_spline_matrix_banded__write

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_banded__free( self )
    class(sll_t_spline_matrix_banded), intent(inout) :: self

    self%n  = -1
    self%kl = -1
    self%ku = -1
    self%factorized = .FALSE.
    if ( allocated( self%ipiv ) ) deallocate( self%ipiv )
    if ( allocated( self%q    ) ) deallocate( self%q    )

  end subroutine s_spline_matrix_banded__free

end module sll_m_spline_matrix_banded
