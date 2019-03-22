!===================================================================================
!>
!!# Module **gvec_to_gene_c_bind**
!!
!!
!!
!===================================================================================

module modgvec_gvec_to_gene_c_bind

  use modgvec_gvec_to_gene, only: &
       init_gvec_to_gene, &
       gvec_to_gene_coords, &
       finalize_gvec_to_gene

  use modgvec_globals, only: WP

  use, intrinsic :: iso_c_binding, only : C_INT, C_FLOAT, C_DOUBLE, C_LONG_DOUBLE, &
       C_CHAR, C_NULL_CHAR

  implicit none
  private

  ! select a proper C working precision (CWP), default C_FLOAT
  ! adjust this parameters according to modgvec_globals
  integer, parameter :: CWP = &
       merge(C_LONG_DOUBLE, &
       merge(C_DOUBLE, C_FLOAT, &
       WP .eq. selected_real_kind(15,307)), &
       WP .eq. selected_real_kind(33,307))

  public init_gvec_to_gene_c, gvec_to_gene_coords_c, finalize_gvec_to_gene_c, &
       test_print_file_name_c, test_pass_arrays_shift_c

contains

  subroutine init_gvec_to_gene_c(fileName) bind(c,name='init_gvec_to_gene')
    character(kind=C_CHAR, len=1), dimension(*) :: fileName
    character(len=:), allocatable :: fileName_f

    call c_to_f_string(fileName, fileName_f)
    call init_gvec_to_gene(fileName_f)
  end subroutine init_gvec_to_gene_c

  subroutine gvec_to_gene_coords_c( &
       nthet, nzeta, spos_in, theta_star_in, &
       zeta_in, theta_out, cart_coords) &
       bind(c,name='gvec_to_gene_coords')
    integer(kind=C_INT), value :: nthet, nzeta
    real(kind=CWP), value :: spos_in
    real(kind=CWP), dimension(nthet,nzeta), intent(in) :: theta_star_in, zeta_in
    real(kind=CWP), intent(out) :: theta_out(nthet,nzeta)
    real(kind=CWP), intent(out) :: cart_coords(3,nthet,nzeta)

    call gvec_to_gene_coords( &
         nthet, nzeta, &
         spos_in, theta_star_in, zeta_in, &
         theta_out, cart_coords)
  end subroutine gvec_to_gene_coords_c

  subroutine finalize_gvec_to_gene_c() bind(c,name='finalize_gvec_to_gene')
    call finalize_gvec_to_gene()
  end subroutine finalize_gvec_to_gene_c

  !===================================================================================

  subroutine test_print_file_name_c(fileName) bind(c,name='test_print_file_name')
    character(kind=C_CHAR, len=1), dimension(*) :: fileName

    character(len=:), allocatable :: fileName_f

    call c_to_f_string(fileName, fileName_f)
    write(*,*) "Tests: ", fileName_f
  end subroutine test_print_file_name_c

  subroutine test_pass_arrays_shift_c( &
       nthet, nzeta, arr_in, arr_out) &
       bind(c, name = "test_pass_arrays_shift")
    integer(kind=C_INT), value :: nthet, nzeta
    real(kind=CWP),intent(in) :: arr_in(nthet,nzeta)
    real(kind=CWP),intent(out) :: arr_out(nthet,nzeta)

    call simple_shift(arr_in, arr_out, nthet*nzeta)

  contains
    subroutine simple_shift(a_in, a_out, size)
      integer(kind=C_INT), intent(in) :: size
      real(kind=CWP),intent(in) :: a_in(size)
      real(kind=CWP),intent(out) :: a_out(size)
      integer :: i
      forall(i=1:size) a_out(i) = a_in(i) + i
    end subroutine simple_shift
  end subroutine test_pass_arrays_shift_c

  !===================================================================================

  subroutine c_to_f_string(s,f)
    character(kind=C_CHAR, len=1), intent(in) :: s(*)
    character(len=:), allocatable, intent(out) :: f
    integer :: nchars

    nchars = 0
    do while( s(nchars+1).ne.C_NULL_CHAR )
       nchars = nchars + 1
    end do

    allocate(character(len=nchars) :: f)
    if( storage_size(f).eq.storage_size(s)*nchars) then
       f = transfer(s(1:nchars), f)
    else
       stop "can't transfer C_CHAR array to fortran character, do explicit copy!"
    end if
  end subroutine c_to_f_string

end module modgvec_gvec_to_gene_c_bind
