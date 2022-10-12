!===================================================================================================================================
! Copyright (c) 2017 - 2018 Florian Hindenlang <hindenlang@gmail.com>
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/hopr)

! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
 
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# MODULE **READ IN TOOLS"
!! 
!! Readin routines for the input file 
!!
!===================================================================================================================================
MODULE MODgvec_ReadInTools
! MODULES
USE MODgvec_Globals
USE ISO_VARYING_STRING
IMPLICIT NONE
PRIVATE

PUBLIC::FillStrings
PUBLIC::TRYREAD
PUBLIC::GETSTR
PUBLIC::CNTSTR
PUBLIC::GETINT
PUBLIC::GETREAL
PUBLIC::GETLOGICAL
PUBLIC::GETINTARRAY
PUBLIC::GETREALARRAY
PUBLIC::GETINTALLOCARRAY
PUBLIC::GETREALALLOCARRAY

PUBLIC::IgnoredStrings
!===================================================================================================================================

INTERFACE TRYREAD
  MODULE PROCEDURE TRYREAD
END INTERFACE

INTERFACE GETSTR
  MODULE PROCEDURE GETSTR
END INTERFACE

INTERFACE CNTSTR
  MODULE PROCEDURE CNTSTR
END INTERFACE

INTERFACE GETINT
  MODULE PROCEDURE GETINT
END INTERFACE

INTERFACE GETREAL
  MODULE PROCEDURE GETREAL
END INTERFACE

INTERFACE GETLOGICAL
  MODULE PROCEDURE GETLOGICAL
END INTERFACE

INTERFACE GETINTARRAY
  MODULE PROCEDURE GETINTARRAY
END INTERFACE

INTERFACE GETREALARRAY
  MODULE PROCEDURE GETREALARRAY
END INTERFACE

INTERFACE IgnoredStrings
  MODULE PROCEDURE IgnoredStrings
END INTERFACE

INTERFACE FillStrings
  MODULE PROCEDURE FillStrings
END INTERFACE

INTERFACE FindStr
  MODULE PROCEDURE FindStr
END INTERFACE

INTERFACE LowCase
  MODULE PROCEDURE LowCase
END INTERFACE

INTERFACE GetNewString
  MODULE PROCEDURE GetNewString
END INTERFACE

INTERFACE DeleteString
  MODULE PROCEDURE DeleteString
END INTERFACE

TYPE tString
  TYPE(Varying_String)::Str
  TYPE(tString),POINTER::NextStr,PrevStr
END TYPE tString

TYPE(tString),POINTER::FirstString

CONTAINS

!===================================================================================================================================
!> Read string from specified unit
!!
!===================================================================================================================================
FUNCTION TRYREAD(UnitLoc,Key,abortOpt)
! MODULES
USE MODgvec_Globals,ONLY:abort
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      !! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: unitLoc
LOGICAL,INTENT(IN),OPTIONAL          :: abortOpt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: TRYREAD
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                              :: stat
CHARACTER(LEN=255)                   :: tmp
LOGICAL                              :: abortLoc=.TRUE.
!===================================================================================================================================
IF(PRESENT(abortOpt)) abortLoc=abortOpt
TRYREAD=.TRUE.
READ(unitLoc,*,IOSTAT=stat) tmp
IF(stat.NE.0)              TRYREAD=.FALSE.
IF(TRIM(Key).NE.TRIM(tmp)) TRYREAD=.FALSE.

IF(.NOT.TRYREAD.AND.abortLoc)&
  CALL abort(__STAMP__,&
             'Keyword '//TRIM(Key)//' not found in file.')
END FUNCTION TRYREAD


!===================================================================================================================================
!> Read string named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
!! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
!! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!!
!===================================================================================================================================
FUNCTION GETSTR(Key,Proposal)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      !! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal !! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255)                   :: GetStr   !! String read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=8)                     :: DefMsg  
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,GetStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,GetStr,DefMsg)
END IF
SWRITE(UNIT_StdOut,'(a3,a30,a3,a33,a3,a7,a3)')' | ',TRIM(Key),' | ', TRIM(GetStr),' | ',TRIM(DefMsg),' | '
END FUNCTION GETSTR

!===================================================================================================================================
!> Counts all occurances of string named "key" from inifile and store in "CNTSTR". If keyword "Key" is not found in ini file,
!! the default value "Proposal" is used for "CNTSTR" (error if "Proposal" not given).
!! Inifile was read in before and is stored as list of character strings starting with "FirstString".
!!
!===================================================================================================================================
FUNCTION CNTSTR(Key,Proposal)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      !! Search for this keyword in ini file
INTEGER         ,OPTIONAL,INTENT(IN) :: Proposal !! Default values as integer 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                              :: CntStr   !! Number of parameters named "Key" in inifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=LEN(Key))              :: TmpKey  
TYPE(tString),POINTER                :: Str1  
!===================================================================================================================================

CntStr=0
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)

! Search
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) CntStr=CntStr+1
  ! Next string in list
  Str1=>Str1%NextStr
END DO
IF (CntStr.EQ.0) THEN
  IF (PRESENT(Proposal)) THEN
    CntStr=Proposal
  ELSE
    SWRITE(UNIT_StdOut,*) 'Inifile missing necessary keyword item : ',TRIM(TmpKey)
    CALL abort(__STAMP__, &
         'Code stopped during inifile parsing!')
  END IF
END IF
END FUNCTION CNTSTR

!===================================================================================================================================
!> Read integer named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
!! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
!! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!!
!===================================================================================================================================
FUNCTION GETINT(Key,Proposal,quiet_def_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: Key          !! Search for this keyword in ini file
INTEGER,OPTIONAL,INTENT(IN) :: Proposal     !! Default values as integer scalar 
LOGICAL,OPTIONAL,INTENT(IN) :: quiet_def_in !! flag to be quiet if DEFAULT is taken 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                     :: GetInt  !! Integer read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)          :: HelpStr,ProposalStr
CHARACTER(LEN=8)            :: DefMsg  
INTEGER                     :: ioerr
LOGICAL                     :: quiet_def
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,intScalar=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*,IOSTAT=ioerr)GetInt
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (integer):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
quiet_def=.FALSE.
IF(PRESENT(quiet_def_in))THEN
  IF(INDEX(TRIM(DefMsg),"DEFAULT").NE.0)THEN
    quiet_def=quiet_def_in
  END IF
END IF
IF(.NOT.quiet_def) THEN
  SWRITE(UNIT_StdOut,'(a3,a30,a3,i33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetInt,' | ',TRIM(DefMsg),' | '
END IF
END FUNCTION GETINT


!===================================================================================================================================
!> Read real named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
!! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
!! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!!
!===================================================================================================================================
FUNCTION GETREAL(Key,Proposal,quiet_def_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key          !! Search for this keyword in ini file
REAL(wp)        ,OPTIONAL,INTENT(IN) :: Proposal     !! Default values as real scalar
LOGICAL         ,OPTIONAL,INTENT(IN) :: quiet_def_in !! flag to be quiet if DEFAULT is taken 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES                                 
REAL(wp)                             :: GetReal  !! Real read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=500)                   :: HelpStr,ProposalStr  
CHARACTER(LEN=8)                     :: DefMsg  
INTEGER                              :: ioerr
LOGICAL                              :: quiet_def
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,realScalar=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
! Find values of pi in the string
READ(HelpStr,*,IOSTAT=ioerr)GetReal
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (real):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
quiet_def=.FALSE.
IF(PRESENT(quiet_def_in))THEN
  IF(INDEX(TRIM(DefMsg),"DEFAULT").NE.0)THEN
    quiet_def=quiet_def_in
  END IF
END IF
IF(.NOT.quiet_def) THEN
  SWRITE(UNIT_StdOut,'(a3,a30,a3,e33.5,a3,a7,a3)')' | ',TRIM(Key),' | ', GetReal,' | ',TRIM(DefMsg),' | '
END IF
END FUNCTION GETREAL


!===================================================================================================================================
!> Read logical named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
!! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
!! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!!
!===================================================================================================================================
FUNCTION GETLOGICAL(Key,Proposal,quiet_def_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key          !! Search for this keyword in ini file
LOGICAL         ,OPTIONAL,INTENT(IN) :: Proposal     !! Default values as logical 
LOGICAL         ,OPTIONAL,INTENT(IN) :: quiet_def_in !! flag to be quiet if DEFAULT is taken 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES                                  
LOGICAL                              :: GetLogical !! Logical read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr,ProposalStr 
CHARACTER(LEN=8)                     :: DefMsg  
INTEGER                              :: ioerr
LOGICAL                              :: quiet_def
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,logScalar=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*,IOSTAT=ioerr)GetLogical
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (logical):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
quiet_def=.FALSE.
IF(PRESENT(quiet_def_in))THEN
  IF(INDEX(TRIM(DefMsg),"DEFAULT").NE.0)THEN
    quiet_def=quiet_def_in
  END IF
END IF
IF(.NOT.quiet_def) THEN
  SWRITE(UNIT_StdOut,'(a3,a30,a3,l33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetLogical,' | ',TRIM(DefMsg),' | '
END IF
END FUNCTION GETLOGICAL


!===================================================================================================================================
!> Read array of "nIntegers" integer values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
!! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
!! list of character strings starting with "FirstString".
!!
!===================================================================================================================================
FUNCTION GETINTARRAY(Key,nIntegers,Proposal,quiet_def_in)
! MODULES
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key              !! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nIntegers        !! Number of values in array
INTEGER         ,OPTIONAL,INTENT(IN) :: Proposal(:)      !! Default values as integer array 
LOGICAL         ,OPTIONAL,INTENT(IN) :: quiet_def_in     !! flag to be quiet if DEFAULT is taken 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                   :: GetIntArray(nIntegers)      !! Integer array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255+nIntegers*50) :: HelpStr,ProposalStr 
CHARACTER(LEN=8)                :: DefMsg  
INTEGER                         :: iInteger  
INTEGER                         :: ioerr
TYPE(varying_string)            :: separator,astr,bstr
LOGICAL                         :: quiet_def
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,intarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
astr=var_str(TRIM(helpstr))
iInteger=0
Separator="X"
DO WHILE(LEN(CHAR(separator)) .NE. 0)
  iInteger=iInteger+1
  CALL split(astr,bstr," ",separator,back=.false.) !bStr is string in front of @
END DO
IF(iInteger.NE.nIntegers)THEN
  WRITE(*,'(A,I4,A,I4)')'PROBLEM IN READIN OF LINE (integer Array), number of elements : ', iInteger, ' .NE. ',nIntegers
  WRITE(*,*) '"',TRIM(key),' = ',TRIM(helpStr),'"'
  STOP     
END IF
READ(HelpStr,*,IOSTAT=ioerr)GetIntArray
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (integer array):'
  WRITE(*,*) '"',TRIM(key),' = ',TRIM(helpStr),'"'
  STOP     
END IF
quiet_def=.FALSE.
IF(PRESENT(quiet_def_in))THEN
  IF(INDEX(TRIM(DefMsg),"DEFAULT").NE.0)THEN
    quiet_def=quiet_def_in
  END IF
END IF
IF(.NOT.quiet_def) THEN
  SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                                 'Integer array of size (',nIntegers,') | ',TRIM(DefMsg),' | '
  DO iInteger=0,nIntegers-1
    IF ((iInteger.GT.0) .AND. (MOD(iInteger,8).EQ.0)) THEN
      SWRITE(UNIT_stdOut,*)
      SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
    END IF
    SWRITE(UNIT_stdOut,'(i5)',ADVANCE='NO')GetIntArray(iInteger+1)
  END DO
  SWRITE(UNIT_stdOut,*)
END IF !quiet_def

END FUNCTION GETINTARRAY


!===================================================================================================================================
!> Allocate and read integer array of unknown length "nIntegers" integer values named "Key" from ini file. 
!! If keyword "Key" is not found in setup file, the default
!! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
!! list of character strings starting with "FirstString".
!!
!===================================================================================================================================
SUBROUTINE GETINTALLOCARRAY(Key,GetIntArray,nIntegers,Proposal,quiet_def_in)
! MODULES
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key          !! Search for this keyword in ini file
INTEGER         ,OPTIONAL,INTENT(IN) :: Proposal(:)  !! Default values as integer array 
LOGICAL         ,OPTIONAL,INTENT(IN) :: quiet_def_in !! flag to be quiet if DEFAULT is taken 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)       :: nIntegers        !! Number of values in array
INTEGER,ALLOCATABLE       :: GetIntArray(:)   !! Integer array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=50*255)     :: HelpStr,ProposalStr 
CHARACTER(LEN=8)          :: DefMsg  
INTEGER                   :: iInteger  
INTEGER                   :: ioerr
TYPE(varying_string)      :: separator,astr,bstr
LOGICAL                   :: quiet_def
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,intarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
astr=var_str(TRIM(helpstr))
nIntegers=0
Separator="X"
DO WHILE(LEN(CHAR(separator)) .NE. 0)
  nIntegers=nIntegers+1
  CALL split(astr,bstr," ",separator,back=.false.) !bStr is string in front of @
END DO
IF(ALLOCATED(GetIntArray)) DEALLOCATE(GetIntArray)
ALLOCATE(GetIntArray(nIntegers))
READ(HelpStr,*,IOSTAT=ioerr)GetIntArray
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (integer array):'
  WRITE(*,*) '"',TRIM(key),' = ',TRIM(helpStr),'"'
  STOP     
END IF
quiet_def=.FALSE.
IF(PRESENT(quiet_def_in))THEN
  IF(INDEX(TRIM(DefMsg),"DEFAULT").NE.0)THEN
    quiet_def=quiet_def_in
  END IF
END IF
IF(.NOT.quiet_def) THEN
  SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                                 'Integer array of size (',nIntegers,') | ',TRIM(DefMsg),' | '
  DO iInteger=0,nIntegers-1
    IF ((iInteger.GT.0) .AND. (MOD(iInteger,8).EQ.0)) THEN
      SWRITE(UNIT_stdOut,*)
      SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
    END IF
    SWRITE(UNIT_stdOut,'(i5)',ADVANCE='NO')GetIntArray(iInteger+1)
  END DO
  SWRITE(UNIT_stdOut,*)
END IF !quiet_def
END SUBROUTINE GETINTALLOCARRAY


!===================================================================================================================================
!> Read array of "nReals" real values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
!! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
!! list of character strings starting with "FirstString".
!!
!===================================================================================================================================
FUNCTION GETREALARRAY(Key,nReals,Proposal,quiet_def_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key          !! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nReals       !! Number of values in array
REAL(wp)        ,OPTIONAL,INTENT(IN) :: Proposal(:)  !! Default values as real array
LOGICAL         ,OPTIONAL,INTENT(IN) :: quiet_def_in !! flag to be quiet if DEFAULT is taken 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)                  :: GetRealArray(nReals)        !! Real array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255+nReals*50) :: HelpStr,ProposalStr 
CHARACTER(LEN=8)             :: DefMsg  
INTEGER                      :: iReal  
INTEGER                      :: ioerr  
TYPE(varying_string)         :: separator,astr,bstr
LOGICAL                      :: quiet_def
!===================================================================================================================================


IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,realarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
astr=var_str(TRIM(helpstr))
iReal=0
Separator="X"
DO WHILE(LEN(CHAR(separator)) .NE. 0)
  iReal=iReal+1
  CALL split(astr,bstr," ",separator,back=.false.) !bStr is string in front of @
END DO
IF(iReal.NE.nReals)THEN
  WRITE(*,'(A,I4,A,I4)')'PROBLEM IN READIN OF LINE (RealArray), number of elements : ', iReal, ' .NE. ',nReals
  WRITE(*,*) '"',TRIM(key),' = ',TRIM(helpStr),'"'
  STOP     
END IF

READ(HelpStr,*,IOSTAT=ioerr)GetRealArray
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (RealArray):'
  WRITE(*,*) '"',TRIM(key),' = ',TRIM(helpStr),'"'
  STOP     
END IF
quiet_def=.FALSE.
IF(PRESENT(quiet_def_in))THEN
  IF(INDEX(TRIM(DefMsg),"DEFAULT").NE.0)THEN
    quiet_def=quiet_def_in
  END IF
END IF
IF(.NOT.quiet_def) THEN
  SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                                 'Real array of size (',nReals,') | ',TRIM(DefMsg),' | '
  DO iReal=0,nReals-1
    IF ((iReal.GT.0) .AND. (MOD(iReal,8).EQ.0)) THEN
      SWRITE(UNIT_stdOut,*)
      SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
    END IF
    SWRITE(UNIT_stdOut,'(f5.2)',ADVANCE='NO')GetRealArray(iReal+1)
  END DO
  SWRITE(UNIT_stdOut,*)
END IF !quiet_def

END FUNCTION GETREALARRAY


!===================================================================================================================================
!> Read array of "nReals" real values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
!! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
!! list of character strings starting with "FirstString".
!!
!===================================================================================================================================
SUBROUTINE GETREALALLOCARRAY(Key,GetRealArray,nReals,Proposal,quiet_def_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key          !! Search for this keyword in ini file
REAL(wp)        ,OPTIONAL,INTENT(IN) :: Proposal(:)  !! Default values as real array
LOGICAL         ,OPTIONAL,INTENT(IN) :: quiet_def_in !! flag to be quiet if DEFAULT is taken 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)       :: nReals           !! Number of values in array
REAL(wp),ALLOCATABLE      :: GetRealArray(:)  !! Real array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=50*255)     :: HelpStr,ProposalStr 
CHARACTER(LEN=8)          :: DefMsg  
INTEGER                   :: iReal  
INTEGER                   :: ioerr  
TYPE(varying_string)      :: separator,astr,bstr
LOGICAL                   :: quiet_def
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,realarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
astr=var_str(TRIM(helpstr))
nReals=0
Separator="X"
DO WHILE(LEN(CHAR(separator)) .NE. 0)
  nReals=nReals+1
  CALL split(astr,bstr," ",separator,back=.false.) !bStr is string in front of @
END DO
IF(ALLOCATED(GetRealarray)) DEALLOCATE(GetRealArray)
ALLOCATE(GetRealArray(nReals))

READ(HelpStr,*,IOSTAT=ioerr)GetRealArray
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (RealArray):'
  WRITE(*,*) '"',TRIM(key),' = ',TRIM(helpStr),'"'
  STOP     
END IF
quiet_def=.FALSE.
IF(PRESENT(quiet_def_in))THEN
  IF(INDEX(TRIM(DefMsg),"DEFAULT").NE.0)THEN
    quiet_def=quiet_def_in
  END IF
END IF
IF(.NOT.quiet_def) THEN
  SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                                 'Real array of size (',nReals,') | ',TRIM(DefMsg),' | '
  DO iReal=0,nReals-1
    IF ((iReal.GT.0) .AND. (MOD(iReal,8).EQ.0)) THEN
      SWRITE(UNIT_stdOut,*)
      SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
    END IF
    SWRITE(UNIT_stdOut,'(f5.2)',ADVANCE='NO')GetRealArray(iReal+1)
  END DO
  SWRITE(UNIT_stdOut,*)
END IF !quiet_def

END SUBROUTINE GETREALALLOCARRAY


!===================================================================================================================================
!> Prints out remaining strings in list after read-in is complete
!!
!===================================================================================================================================
SUBROUTINE IgnoredStrings()
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tString),POINTER                  :: Str1  
!===================================================================================================================================
IF(MPIroot)THEN !<<<<
  Str1=>FirstString
  IF(ASSOCIATED(str1))THEN
    WRITE(UNIT_stdOut,'(132("-"))')
    WRITE(UNIT_stdOut,'(A)')" THE FOLLOWING INI-FILE PARAMETERS WERE IGNORED:"
    DO WHILE(ASSOCIATED(Str1))
      WRITE(UNIT_stdOut,'(A4,A)')" |- ",TRIM(CHAR(Str1%Str))
      Str1=>Str1%NextStr
    END DO
    WRITE(UNIT_stdOut,'(132("-"))')
  END IF
END IF !MPIroot !<<<<
END SUBROUTINE IgnoredStrings


!===================================================================================================================================
!> Read ini file and put each line in a string object. All string objects are connected to a list of string objects starting
!! with "firstString". MUST BE CALLED IN THE VERY BEGINNING OF THE PROGRAM!
!===================================================================================================================================
SUBROUTINE FillStrings(IniFile)
! MODULES
USE_MPI !<<<<
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)            :: IniFile                    !! Name of ini file to be read in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tString),POINTER          :: Str1=>NULL(),Str2=>NULL()  
CHARACTER(LEN=255)             :: HelpStr,Str  
CHARACTER(LEN=300)             :: Filename  
TYPE(Varying_String)           :: aStr,bStr,Separator  
INTEGER                        :: stat,iniUnit,nLines,i,iError !<<<<
LOGICAL                        :: file_exists !<<<<
CHARACTER(LEN=255),ALLOCATABLE :: FileContent(:) !<<<<
CHARACTER(LEN=1)               :: tmpChar='' !<<<<
!===================================================================================================================================
!READ FROM FILE ONLY ON MPIroot
IF(MPIroot)THEN !<<<<
  FileName = TRIM(IniFile)
  ! Get name of ini file
  WRITE(UNIT_StdOut,*)'| Reading from file "',TRIM(filename),'":'
  INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
  IF (.NOT.file_exists) THEN
    CALL Abort(__STAMP__,&
        "Ini file does not exist.")
  END IF

  OPEN(NEWUNIT= iniUnit,        &
       FILE   = TRIM(filename), &
       STATUS = 'OLD',          &
       ACTION = 'READ',         &
       ACCESS = 'SEQUENTIAL',   &
       IOSTAT = stat)
  IF(stat.NE.0)THEN
    CALL abort(__STAMP__,&
      "Could not open ini file.")
  END IF

  ! parallel IO: ROOT reads file and sends it to all other procs
  nLines=0
  stat=0
  DO
    READ(iniunit,"(A)",IOSTAT=stat)tmpChar
    IF(stat.NE.0)EXIT
    nLines=nLines+1
  END DO
END IF !MPIroot !<<<<

!broadcast number of lines, read and broadcast file content
#if MPI
CALL MPI_BCAST(nLines,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError) !<<<<
#endif
ALLOCATE(FileContent(nLines))

IF(MPIroot)THEN !<<<<
  !read file
  REWIND(iniUnit)
  READ(iniUnit,'(A)') FileContent
  CLOSE(iniUnit)
END IF !MPIroot !<<<<
!BROADCAST FileContent
#if MPI
CALL MPI_BCAST(FileContent,LEN(FileContent)*nLines,MPI_CHARACTER,0,MPI_COMM_WORLD,iError) !<<<<
#endif

NULLIFY(Str1,Str2)
DO i=1,nLines !<<<<
  IF(.NOT.ASSOCIATED(Str1)) CALL GetNewString(Str1)
  ! Read line from file
  aStr=FileContent(i)
  Str=aStr
  ! Remove comments with "!"
  CALL Split(aStr,Str1%Str,"!")
  ! Remove comments with "#"
  CALL Split(Str1%Str,bStr,"#")
  Str1%Str=bStr
  ! Remove "%" sign from old ini files, i.e. mesh% disc% etc.
  CALL Split(Str1%Str,bStr,"%",Separator,Back=.false.)
  ! If we have a newtype ini file, take the other part
  IF(LEN(CHAR(Separator)).EQ.0) Str1%Str=bStr
  ! Remove blanks
  Str1%Str=Replace(Str1%Str," ","",Every=.true.)
  ! Replace brackets
  Str1%Str=Replace(Str1%Str,"(/","",Every=.true.)
  Str1%Str=Replace(Str1%Str,"/)","",Every=.true.)
  ! Replace commas
  Str1%Str=Replace(Str1%Str,","," ",Every=.true.)
  ! Lower case
  CALL LowCase(CHAR(Str1%Str),HelpStr)
  ! If we have a remainder (no comment only)
  IF(LEN_TRIM(HelpStr).GT.2) THEN
    Str1%Str=Var_Str(HelpStr)
    IF(.NOT.ASSOCIATED(Str2)) THEN
      FirstString=>Str1
    ELSE
      Str2%NextStr=>Str1
      Str1%PrevStr=>Str2
    END IF
    Str2=>Str1
    CALL GetNewString(Str1)
  END IF
END DO

!find line continuation "&" and merge strings (can be multiple lines)
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF(INDEX(CHAR(Str1%str),'&').NE.0)THEN !found "&"
    CALL Split(Str1%Str,aStr,"&") !take part in front of "&"
    Str2=>Str1%nextStr
    Str1%Str=Var_str(CHAR(aStr)//CHAR(Str2%Str))
    CALL deleteString(Str2) 
    !do not go to next  string as long as there are "&" in the string  
  ELSE
    Str1=>Str1%NextStr !nothing to be done
  END IF
END DO

CALL UserDefinedVars()


END SUBROUTINE FillStrings

!===================================================================================================================================
!> Get the user defined variables 
!!
!===================================================================================================================================
SUBROUTINE UserDefinedVars()
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                          :: i,j,nDefVars
TYPE(Varying_String),ALLOCATABLE :: DefVar(:,:)
TYPE(tString),POINTER            :: Str1  
TYPE(Varying_String)             :: vStr,vStr1,vStr2
LOGICAL                          :: found
!===================================================================================================================================
nDefVars=CNTSTR('defvar',0)
IF(nDefVars.EQ.0) RETURN
SWRITE(UNIT_StdOut,'(A,I4,A)')' | Found ',nDefVars,' UserDefined variables: '
ALLOCATE(DefVar(2,nDefVars))
DO i=1,nDefVars
  CALL GetDefVar(DefVar(:,i))
  !check if part of the variable name was used before
  DO j=1,i-1
    IF (INDEX(TRIM(CHAR(DefVar(1,i))),TRIM(CHAR(DefVar(1,j)))).NE.0) THEN
      SWRITE(UNIT_StdOut,*) '!! WARNING !! Problem with DEFVAR ', TRIM(CHAR(DefVar(1,i)))
      SWRITE(UNIT_StdOut,*) '  a part of this variable name was already used in DEFVAR ' ,TRIM(CHAR(DefVar(1,j)))
      CALL abort(__STAMP__, &
         'DEFVAR: do not reuse same strings for variable names! Code stopped during inifile parsing!')
    END IF
  END DO
  Str1=>FirstString
  DO WHILE(ASSOCIATED(Str1))
    vStr=Str1%Str
    vStr2=Str1%Str
    CALL Split(vStr2,vStr1,"=",back=.FALSE.) 
    found=.FALSE.
    IF (INDEX(TRIM(CHAR(vStr2)),TRIM(CHAR(DefVar(1,i)))).NE.0) THEN
      found=.TRUE.
      vStr2=replace(vStr2,TRIM(CHAR(DefVar(1,i))),TRIM(CHAR(DefVar(2,i))),Every=.TRUE.)
    END IF
    IF(Found)THEN
      !SWRITE(UNIT_StdOut,*)'DEBUG, ',TRIM(CHAR(Str1%str))
      Str1%Str=CHAR(vStr1)//'='//CHAR(vStr2)
      !SWRITE(UNIT_StdOut,*)' >>>>>>',TRIM(CHAR(Str1%str))
    END IF !found
    ! Next string in list
    Str1=>Str1%NextStr
  END DO !WHILE Str1 associated
END DO !i=1,nDefVars

END SUBROUTINE UserDefinedVars 


!===================================================================================================================================
!> Get the user defined variables 
!!
!===================================================================================================================================
SUBROUTINE GetDefVar(DefVar)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(Varying_String):: DefVar(2)   !! Name, Value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)  :: HelpStr  
CHARACTER(LEN=255)  :: aStr  
CHARACTER(LEN=8)    :: DefMsg  
TYPE(Varying_String):: vStr,vStr1,vStr2,vStrTmp,vStr_narr
INTEGER             :: DefVarInt
REAL(wp)            :: DefVarReal
INTEGER,ALLOCATABLE :: DefVarIntArr(:)
REAL(wp),ALLOCATABLE:: DefVarRealArr(:)
INTEGER             :: nArr
LOGICAL             :: DefVarIsInt
LOGICAL             :: DefVarIsIntarray
LOGICAL             :: DefVarIsReal
LOGICAL             :: DefVarIsRealarray
!===================================================================================================================================
CALL FindStr('DEFVAR',HelpStr,DefMsg)
CALL LowCase(HelpStr,aStr)
vStr=aStr

CALL Split(vStr,vStr1,":",back=.FALSE.) !split after first occurence in vStr, and put first part in vStr1, second part in vStr
CALL Split(vStr,vStr2,":",back=.TRUE.)  !split after last occurence in Vstr and put second part in vStr1 and first part in Vstr

vStrTmp=vStr1
CALL Split(vStrtmp,vStr1,"~",back=.FALSE.) !first part in vStr1, second part in VStrtmp
CALL Split(vStrtmp,vStr_narr,"~",back=.TRUE.) !second part in VStr_narr, first part in VStrtmp
CALL Split(vStr_narr,vStrTmp,")",back=.TRUE.) !first part VStr_narr

DefVarIsInt      =(CHAR(vStr1).EQ.'(int)') 
DefVarIsIntarray =(CHAR(vStr1).EQ.'(int') !array  must be of format (int~n)
DefVarIsReal     =(CHAR(vStr1).EQ.'(real)') 
DefVarIsRealarray=(CHAR(vStr1).EQ.'(real') 


IF(.NOT.((DefVarIsInt).OR.(DefVarIsIntArray).OR.(DefVarIsReal).OR.(defVarIsRealarray) ))THEN
  SWRITE(UNIT_StdOut,*) 'DEFVAR not correctly defined: ',TRIM(HelpStr)
    CALL abort(__STAMP__, &
         'Code stopped during inifile parsing!')
END IF

IF(DefVarIsIntArray.OR.DefVarIsRealArray)THEN
  aStr=CHAR(vstr_narr)
  READ(aStr,*) nArr
END IF

!now take the second part of the definition ( nvar = xxx)
vStr=VStr2
CALL Split(vStr,vStr1,"=",back=.FALSE.) 
CALL Split(vStr,vStr2,"=",back=.TRUE.) 
DefVar(1)=vStr1
DefVar(2)=vStr2
aStr=CHAR(DefVar(2))
IF(DefVarIsInt) THEN
  READ(aStr,*)DefVarInt
  SWRITE(UNIT_StdOut,'(A3,A30,A3,I33,A3,A7,A3)')  ' | ',TRIM(CHAR(DefVar(1))),' | ', DefVarInt      ,' | ','=>INT  ',' | '
ELSE IF(DefVarIsReal)THEN
  READ(aStr,*)DefVarReal
  SWRITE(UNIT_StdOut,'(A3,A30,A3,E33.5,A3,A7,A3)')' | ',TRIM(CHAR(DefVar(1))),' | ', DefVarReal     ,' | ','=>REAL ',' | '
ELSE IF(DefVarIsIntArray)THEN
  ALLOCATE(DefVarIntArr(nArr))
  READ(aStr,*)DefVarIntArr(:)
  SWRITE(UNIT_StdOut,'(A3,A30,A31,I4,A4,A7,A3,'//TRIM(CHAR(vStr_narr))//'(X,I4))')' | ',TRIM(CHAR(DefVar(1))), &
        ' |      Integer array of size (', nArr,') | ','=>INT  ',' | ',DefVarIntArr(1:narr)
  WRITE(aStr,'('//TRIM(CHAR(vStr_narr))//'(X,I8))') DefVarIntArr !overwrite
  DEALLOCATE(DefVarIntArr)
  DefVar(2)=aStr
ELSE IF(DefVarIsRealArray)THEN
  ALLOCATE(DefVarRealArr(nArr))
  READ(aStr,*)DefVarRealArr(:)
  SWRITE(UNIT_StdOut,'(A3,A30,A31,I4,A4,A7,A3,'//TRIM(CHAR(vStr_narr))//'(X,E8.1))')' | ',TRIM(CHAR(DefVar(1))), &
        ' |         Real array of size (', nArr,') | ','=>REAL ',' | ',DefVarRealArr(1:narr)
  WRITE(aStr,'('//TRIM(CHAR(vStr_narr))//'(X,E23.15))') DefVarRealArr !overwrite
  DEALLOCATE(DefVarRealArr)
  DefVar(2)=aStr
END IF

END SUBROUTINE GetDefVar 


!===================================================================================================================================
!> Create and initialize new string object.
!!
!===================================================================================================================================
SUBROUTINE GetNewString(Str)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tString),POINTER,INTENT(INOUT) :: Str !! New string
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
NULLIFY(Str)
ALLOCATE(Str)
NULLIFY(Str%NextStr,Str%PrevStr)
END SUBROUTINE GetNewString


!===================================================================================================================================
!> Remove string "Str" from list of strings witFirstString,h first element "DirstString" and delete string.
!!
!===================================================================================================================================
SUBROUTINE DeleteString(Str)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tString),POINTER,INTENT(INOUT) :: Str         !! String to delete
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF (ASSOCIATED(Str%NextStr)) Str%NextStr%PrevStr=>Str%PrevStr
IF (ASSOCIATED(Str,FirstString)) THEN
  FirstString=>Str%NextStr
ELSE
  Str%PrevStr%NextStr=>Str%NextStr
END IF
DEALLOCATE(Str)
NULLIFY(Str)
END SUBROUTINE DeleteString


!===================================================================================================================================
!> Find parameter string containing keyword "Key" in list of strings starting with "FirstString" and return string "Str" without
!! keyword. If keyword is not found in list of strings, return default values "Proposal" (error if not given).
!! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!!
!===================================================================================================================================
SUBROUTINE FindStr(Key,Str,DefMsg,Proposal)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key         !! Search for this keyword in ini file
CHARACTER(LEN=8),INTENT(INOUT)       :: DefMsg      !! Default message = keyword not found, return default parameters (if available)
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal    !! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT)         :: Str         !! Parameter string without keyword
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=LEN(Key))              :: TmpKey   
TYPE(tString),POINTER                :: Str1  
LOGICAL                              :: Found  
!===================================================================================================================================
DefMsg='*CUSTOM'
! Convert to lower case
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)
Found=.FALSE.
Str1=>FirstString
DO WHILE(.NOT.Found)
  IF (.NOT.ASSOCIATED(Str1)) THEN
    IF (.NOT.PRESENT(Proposal)) THEN
      SWRITE(UNIT_StdOut,*) 'Inifile missing necessary keyword item : ',TRIM(TmpKey)
      CALL abort(__STAMP__, &
           'Code stopped during inifile parsing!')
    ELSE ! Return default value
!      CALL LowCase(TRIM(Proposal),Str)
      Str=TRIM(Proposal)
      IF (Str(1:1).NE.'@') THEN
        DefMsg='DEFAULT'
      END IF
      RETURN
    END IF ! (.NOT.PRESENT(Proposal))
  END IF ! (.NOT.ASSOCIATED(Str1))

  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) THEN
    Found=.TRUE.
    Str1%Str=replace(Str1%Str,TRIM(TmpKey)//'=',"",Every=.TRUE.)
    Str=TRIM(CHAR(Str1%Str))
    ! Remove string from list
    CALL DeleteString(Str1)
  ELSE
    ! Next string in list
    Str1=>Str1%NextStr
  END IF

END DO
END SUBROUTINE FindStr


!===================================================================================================================================
!> Transform upper case letters in "Str1" into lower case letters, result is "Str2"
!!
!==================================================================================================================================
SUBROUTINE LowCase(Str1,Str2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: Str1 !! Input string 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT) :: Str2 !! Output string, lower case letters only
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: iLen,nLen,Upper  
CHARACTER(LEN=*),PARAMETER   :: lc='abcdefghijklmnopqrstuvwxyz'  
CHARACTER(LEN=*),PARAMETER   :: UC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'  
LOGICAL                      :: HasEq  
!===================================================================================================================================
HasEq=.FALSE.
Str2=Str1
nLen=LEN_TRIM(Str1)
DO iLen=1,nLen
  ! Transformation stops at "="
  IF(Str1(iLen:iLen).EQ.'=') HasEq=.TRUE.
  Upper=INDEX(UC,Str1(iLen:iLen))
  IF ((Upper > 0).AND. .NOT. HasEq) Str2(iLen:iLen) = lc(Upper:Upper)
END DO
END SUBROUTINE LowCase

!===================================================================================================================================
!> Get logical, integer, real, integer array or real array and transform it to string in the proposal format
!!
!===================================================================================================================================
SUBROUTINE ConvertToProposalStr(ProposalStr,LogScalar,IntScalar,realScalar,Intarr,realarr)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT        VARIABLES
LOGICAL ,INTENT(IN),OPTIONAL   :: LogScalar
INTEGER ,INTENT(IN),OPTIONAL   :: intScalar
REAL(wp),INTENT(IN),OPTIONAL   :: realScalar
INTEGER ,INTENT(IN),OPTIONAL   :: intarr(:)
REAL(wp),INTENT(IN),OPTIONAL   :: realarr(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: ProposalStr
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=LEN(ProposalStr)) :: str_tmp
TYPE(VARYING_STRING) :: tmpstr
!===================================================================================================================================
  IF(PRESENT(logScalar))THEN
    IF(logScalar)THEn
      str_tmp='T'
    ELSE
      str_tmp='F'
    END IF
  ELSEIF(PRESENT(intscalar))THEN
    WRITE(str_tmp,'(I10)')intScalar
  ELSEIF(PRESENT(realScalar))THEN
    WRITE(str_tmp,'(E23.15)')realScalar
  ELSEIF(PRESENT(intarr))THEN
    WRITE(str_tmp,'(*(I8,:,","))')intarr(:)
  ELSEIF(PRESENT(realarr))THEN
    WRITE(str_tmp,'(*(E21.11,:,","))')realarr(:)
  ELSE 
    ProposalStr=" "
    RETURN
  END IF
  tmpstr=VAR_STR(str_tmp)
  tmpstr=Replace(tmpstr," ","",Every=.true.)
  tmpstr=Replace(tmpstr,","," ",Every=.true.)
  ProposalStr=TRIM(CHAR(tmpstr))
END SUBROUTINE ConvertToProposalStr

END MODULE MODgvec_ReadInTools
