!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================

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

  use modgvec_output_vtk, only: WriteDataToVTK

  use, intrinsic :: iso_c_binding, only : C_INT, C_FLOAT, C_DOUBLE, C_LONG_DOUBLE, &
       C_CHAR, C_NULL_CHAR, C_PTR, c_f_pointer

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
       test_print_file_name_c, test_pass_arrays_shift_c, test_int_array_c

contains

  subroutine write_data_to_vtk_c( &
       dim1,vecDim,nVal,NPlot,nElems,strlen,VarNames_c,Coord,Values,FileString_c) &
       bind(c,name='write_data_to_vtk')
    ! INPUT VARIABLES
    ! dimension of the data (either 2=quads or 3=hexas)
    INTEGER(kind=C_INT), value :: dim1
    INTEGER(kind=C_INT), value :: vecdim ! dimension of coordinates
    INTEGER(kind=C_INT), value :: nVal ! Number of nodal output variables
    ! Number of output points per element: (nPlot+1)**dim1
    INTEGER(kind=C_INT),INTENT(IN) :: NPlot(dim1)
    INTEGER(kind=C_INT), value :: nElems ! Number of output elements
    ! CoordinatesVector
    REAL(kind=CWP),INTENT(IN) :: Coord(vecdim,1:PRODUCT(Nplot+1),nElems)
    integer(kind=C_INT), value :: strlen ! lenght of strings in VarNames_c
    ! Names of all variables that will be written out
    character(kind=C_CHAR,len=1), intent(in) :: VarNames_c(strlen,nVal)
    REAL(kind=CWP),INTENT(IN) :: Values(nVal,1:PRODUCT(Nplot+1),nElems) ! Statevector
    ! Output file name
    character(kind=C_CHAR, len=1), dimension(*) :: FileString_c

    character(len=strlen) :: VarNames_f(nVal)
    character(len=strlen) :: VarName
    character(len=:), allocatable :: FileString_f
    integer :: val, nchar1, nchar2

    nchar1 = get_c_string_length(FileString_c)

    if(nchar1.gt.256) stop "c string length is restricted to 256"

    ! retrieve file name into Fortran string
    call c_to_f_string(FileString_c, FileString_f)

    ! save all var names into VarNames_f
    do val = 1, nVal
      nchar2 = get_c_string_length(VarNames_c(:,val))
      if(nchar2.gt.strlen) stop "c string length is restricted to strlen"
      VarName = get_fortran_string(VarNames_c(:,val),nchar2)
      VarNames_f(val) = trim(VarName(1:nchar2))
    end do

    call WriteDataToVTK( &
         dim1,vecDim,nVal,NPlot,nElems,VarNames_f,Coord,Values,FileString_f)

  end subroutine write_data_to_vtk_c

  subroutine init_gvec_to_gene_c(fileName) bind(c,name='init_gvec_to_gene')
    character(kind=C_CHAR, len=1), dimension(*) :: fileName
    character(len=:), allocatable :: fileName_f
    integer :: nchar

    nchar = get_c_string_length(fileName)
    if(nchar.gt.256) stop "c string length is restricted to 256"
    call c_to_f_string(fileName, fileName_f)
    call init_gvec_to_gene(fileName_f(1:nchar))

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

  subroutine test_int_array_c(dim1,Nplot) &
       bind(c,name='test_int_array')

    INTEGER(kind=C_INT), value :: dim1
    INTEGER(kind=C_INT),INTENT(IN) :: NPlot(dim1)

    write(*,*) Nplot, "Nplot"
    write(*,*) "product", PRODUCT(Nplot+1)

  end subroutine test_int_array_c

subroutine test_print_char_rank2_array_c(strlen,nval,varnames_c) &
       bind(c,name='test_print_char_rank2_array')
    integer(kind=C_INT), value, intent(in) :: nval, strlen
    character(kind=C_CHAR,len=1), intent(in) :: varnames_c(strlen,nval)

    character(len=strlen) :: VarNames_f(nval)
    character(len=strlen) :: varname
    integer(kind=C_INT) :: val
    integer :: nchar

    write(*,*) "test_print_char_rank2_array nval: ", nval
    do val = 1, nval
      nchar = get_c_string_length(varnames_c(:,val))
      if(nchar.gt.strlen) stop "c string length is restricted to strlen"
      varname = get_fortran_string(varnames_c(:,val),nchar)
      VarNames_f(val) = trim(varname(1:nchar))
      write(*,*), val, ": VarNames_f: ", VarNames_f(val)
    end do

  end subroutine test_print_char_rank2_array_c

  subroutine test_print_file_name_c(fileName) bind(c,name='test_print_file_name')
    character(kind=C_CHAR, len=1), dimension(*) :: fileName

    character(len=256) :: fileName_f
    integer :: nchar

    nchar = get_c_string_length(fileName)

    if(nchar.gt.256) stop "c string length is restricted to 256"

    fileName_f = get_fortran_string(fileName,nchar)

    write(*,*) "Tests: ", fileName_f(1:nchar)
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

  function get_c_string_length(s) result(nchars)
    character(kind=C_CHAR, len=1), intent(in) :: s(*)
    integer :: nchars

    nchars = 0
    do while( s(nchars+1).ne.C_NULL_CHAR )
      nchars = nchars + 1
    end do
  end function get_c_string_length

  function get_fortran_string(s,nchars) result(f)
    character(kind=C_CHAR, len=1), intent(in) :: s(*)
    integer, intent(in) :: nchars
    character(len=nchars) :: f

    if( storage_size(f).eq.storage_size(s)*nchars) then
      f = transfer(s(1:nchars), f)
    else
      stop "can't transfer C_CHAR array to fortran character, do explicit copy!"
    end if
  end function get_fortran_string

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
