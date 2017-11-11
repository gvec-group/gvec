!===================================================================================================================================
! Copyright (c) 2017-2018 Florian Hindenlang <hindenlang@gmail.com>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!===================================================================================================================================

! Abbrevations
#ifndef __FILENAME__ 
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#if MPI
#  define SWRITE IF(MPIRoot) WRITE
#else
#  define SWRITE WRITE
#endif
#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)

!boundary condition for zero,odd and even m modes

!boundary types
#define NBC_TYPES    6

#define BC_TYPE_OPEN      1 
#define BC_TYPE_NEUMANN   2
#define BC_TYPE_DIRICHLET 3
#define BC_TYPE_SYMM      4
#define BC_TYPE_SYMMZERO  5
#define BC_TYPE_ANTISYMM  6

!grid types
#define GRID_TYPE_UNIFORM 1 
#define GRID_TYPE_SQRT_S  2 
#define GRID_TYPE_S2      3 
#define GRID_TYPE_BUMP    4 
