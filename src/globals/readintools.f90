!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
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
USE MODgvec_Globals, ONLY:wp,UNIT_stdout,MPIroot,abort,MAXLEN,enter_subregion,exit_subregion
!USE ISO_VARYING_STRING
IMPLICIT NONE
PRIVATE


PUBLIC::GETSTR
PUBLIC::CNTSTR
PUBLIC::GETINT
PUBLIC::GETREAL
PUBLIC::GETLOGICAL
PUBLIC::GETINTARRAY
PUBLIC::GETREALARRAY
PUBLIC::GETINTALLOCARRAY
PUBLIC::GETREALALLOCARRAY
PUBLIC::remove_blanks,replace,split

PUBLIC::IgnoredStrings
PUBLIC::FillStrings
PUBLIC::FinalizeReadIn
!===================================================================================================================================

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
#if defined(NVHPC)
  CHARACTER(LEN=MAXLEN) :: Str  !! ONLY NVHPC COMPILER DOES NOT SEEM TO WORK WITH ALLOCATABLE CHARACTERS (SIGSEV!)
#else
  CHARACTER(LEN=:),ALLOCATABLE::Str
#endif
  TYPE(tString),POINTER::NextStr,PrevStr
END TYPE tString


LOGICAL,PUBLIC::ReadInDone=.FALSE.
TYPE(tString),POINTER::FirstString


CONTAINS

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
CHARACTER(LEN=512)                   :: GetStr   !! String read from setup file or initialized with default value
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
TmpKey=remove_blanks(TmpKey)

! Search
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF (INDEX(Str1%Str,TRIM(TmpKey)//'=').EQ.1) CntStr=CntStr+1
  ! Next string in list
  Str1=>Str1%NextStr
END DO
IF (CntStr.EQ.0) THEN
  IF (PRESENT(Proposal)) THEN
    CntStr=Proposal
  ELSE
    CALL abort(__STAMP__, &
         "missing necessary parameter '"//TRIM(TmpKey)//"'", &
         TypeInfo="MissingParameterError")
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
CHARACTER(LEN=MAXLEN)       :: HelpStr,ProposalStr
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
READ(HelpStr,'(I8)',IOSTAT=ioerr)GetInt
IF(ioerr.NE.0)THEN
  CALL abort(__STAMP__, &
       "Problem reading parameter '"//TRIM(key)//"', expected integer, got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
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
!> Read real named "key" from setup file and store in "GETREAL". If keyword "Key" is not found in ini file,
!! the default value "Proposal" is used for "GETREAL" (error if "Proposal" not given).
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
CHARACTER(LEN=MAXLEN)                :: HelpStr,ProposalStr
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
  CALL abort(__STAMP__,&
       "Problem reading parameter '"//TRIM(key)//"', expected real, got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
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
!> Read logical named "key" from setup file and store in "GETLOGICAL". If keyword "Key" is not found in ini file,
!! the default value "Proposal" is used for "GETLOGICAL" (error if "Proposal" not given).
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
CHARACTER(LEN=MAXLEN)                :: HelpStr,ProposalStr
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
  WRITE(UNIT_stdout,*)'PROBLEM IN READIN OF LINE (logical):'
  WRITE(UNIT_stdout,*) TRIM(key),' = ',TRIM(helpStr)
  CALL abort(__STAMP__, &
       "Problem reading parameter '"//TRIM(key)//"', expected logical, got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
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
CHARACTER(LEN=MAXLEN)           :: HelpStr,ProposalStr
CHARACTER(LEN=8)                :: tmpstrarray(nIntegers)
CHARACTER(LEN=8)                :: DefMsg,tmpstr
INTEGER                         :: iInteger
INTEGER                         :: ioerr
LOGICAL                         :: quiet_def
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,intarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
iInteger=1+count_sep(Key,helpstr,",")

IF(iInteger.NE.nIntegers)THEN
  WRITE(tmpstr,'(I8)')nIntegers
  CALL abort(__STAMP__,&
       "Problem reading parameter '"//TRIM(key)//"', expected integer array of size "//TRIM(tmpstr)//", got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
END IF

READ(HelpStr,*,IOSTAT=ioerr)tmpstrarray
IF(ioerr.NE.0)THEN
  CALL abort(__STAMP__,&
       "Problem reading parameter '"//TRIM(key)//"', expected integer array, got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
END IF
DO iInteger=1,nIntegers
   READ(tmpstrarray(iInteger),"(I8)",IOSTAT=ioerr) getIntArray(iInteger)
   IF(ioerr.NE.0)THEN
     CALL abort(__STAMP__,&
          "Problem reading parameter '"//TRIM(key)//"', expected integer array, got '= "//TRIM(helpStr)//"'", &
          TypeInfo="InvalidParameterError")
   END IF
END DO
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
CHARACTER(LEN=MAXLEN)     :: HelpStr,ProposalStr
CHARACTER(LEN=8)          :: DefMsg,tmpstr
INTEGER                   :: iInteger
INTEGER                   :: ioerr
LOGICAL                   :: quiet_def
CHARACTER(LEN=8),ALLOCATABLE :: tmpstrarray(:)
!===================================================================================================================================

IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,intarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
nIntegers=1+count_sep(Key,helpstr,",")

IF(ALLOCATED(GetIntArray)) DEALLOCATE(GetIntArray)
ALLOCATE(GetIntArray(nIntegers),tmpstrarray(nIntegers))
READ(HelpStr,*,IOSTAT=ioerr)tmpstrarray
IF(ioerr.NE.0)THEN
  WRITE(tmpstr,'(I8)')nIntegers
  CALL abort(__STAMP__,&
       "Problem reading parameter '"//TRIM(key)//"', expected integer array of size "//TRIM(tmpstr)//", got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
END IF
DO iInteger=1,nIntegers
  READ(tmpstrarray(iInteger),"(I8)",IOSTAT=ioerr) GetIntArray(iInteger)
  IF(ioerr.NE.0)THEN
    WRITE(tmpstr,'(I8)')nIntegers
    CALL abort(__STAMP__,&
         "Problem reading parameter '"//TRIM(key)//"', expected integer array of size "//TRIM(tmpstr)//", got '= "//TRIM(helpStr)//"'", &
         TypeInfo="InvalidParameterError")
  END IF
END DO
DEALLOCATE(tmpstrarray)
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
CHARACTER(LEN=MAXLEN)        :: HelpStr,ProposalStr
CHARACTER(LEN=8)             :: DefMsg,tmpstr
INTEGER                      :: iReal
INTEGER                      :: ioerr
LOGICAL                      :: quiet_def
!===================================================================================================================================


IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,realarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
iReal=1+count_sep(Key,helpstr,",")
IF(iReal.NE.nReals)THEN
  WRITE(tmpstr,'(I8)')nReals
  CALL abort(__STAMP__,&
       "Problem reading parameter '"//TRIM(key)//"', expected real array of size "//TRIM(tmpstr)//", got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
END IF

READ(HelpStr,*,IOSTAT=ioerr)GetRealArray
IF(ioerr.NE.0)THEN
  CALL abort(__STAMP__, &
       "Problem reading parameter '"//TRIM(key)//"', expected real array, got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
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
CHARACTER(LEN=MAXLEN)     :: HelpStr,ProposalStr
CHARACTER(LEN=8)          :: DefMsg,tmpstr
INTEGER                   :: iReal
INTEGER                   :: ioerr
LOGICAL                   :: quiet_def
!===================================================================================================================================
IF (PRESENT(Proposal)) THEN
  CALL ConvertToProposalStr(ProposalStr,realarr=Proposal)
  CALL FindStr(Key,HelpStr,DefMsg,ProposalStr)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
!count number of components
nReals=1+count_sep(Key,helpstr,",")

IF(ALLOCATED(GetRealarray)) DEALLOCATE(GetRealArray)
ALLOCATE(GetRealArray(nReals))

READ(HelpStr,*,IOSTAT=ioerr)GetRealArray
IF(ioerr.NE.0)THEN
  WRITE(tmpstr,'(I8)')nReals
  CALL abort(__STAMP__,&
       "Problem reading parameter '"//TRIM(key)//"', expected real array of size "//TRIM(tmpstr)//", got '= "//TRIM(helpStr)//"'", &
       TypeInfo="InvalidParameterError")
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
TYPE(tString),POINTER                  :: Str1, Str2
!===================================================================================================================================
IF(MPIroot)THEN !<<<<
  Str1=>FirstString
  IF(ASSOCIATED(str1))THEN
    WRITE(UNIT_stdOut,'(132("-"))')
    WRITE(UNIT_stdOut,'(A)')" THE FOLLOWING INI-FILE PARAMETERS WERE IGNORED:"
    DO WHILE(ASSOCIATED(Str1))
      WRITE(UNIT_stdOut,'(A4,A)')" |- ",TRIM(Str1%Str)
      Str2=>Str1%NextStr
      CALL DeleteString(Str1) ! remove string from the list -> no strings should be left
      Str1=>Str2
    END DO
    WRITE(UNIT_stdOut,'(132("-"))')
  END IF
END IF !MPIroot !<<<<
END SUBROUTINE IgnoredStrings


!===================================================================================================================================
!> Reset global variables
!!
!===================================================================================================================================
SUBROUTINE FinalizeReadIn()
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (ASSOCIATED(FirstString)) THEN
  CALL IgnoredStrings()
END IF
ReadInDone=.FALSE.
END SUBROUTINE FinalizeReadIn


!===================================================================================================================================
!> Read ini file and put each line in a string object. All string objects are connected to a list of string objects starting
!! with "firstString". MUST BE CALLED IN THE VERY BEGINNING OF THE PROGRAM!
!===================================================================================================================================
SUBROUTINE FillStrings(IniFile)
! MODULES
USE MODgvec_MPI, ONLY: par_BCast
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)            :: IniFile                    !! Name of ini file to be read in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tString),POINTER          :: Str1=>NULL(),Str2=>NULL()
CHARACTER(LEN=MAXLEN)          :: HelpStr,Str
CHARACTER(LEN=300)             :: Filename
INTEGER                        :: stat,iniUnit,nLines,i !<<<<
LOGICAL                        :: file_exists !<<<<
CHARACTER(LEN=MAXLEN),ALLOCATABLE :: FileContent(:) !<<<<
CHARACTER(LEN=1)               :: tmpChar='' !<<<<
!===================================================================================================================================
! do nothing if FillStrings was already called
IF (ReadInDone) RETURN
CALL enter_subregion("read-parameterfile")
!READ FROM FILE ONLY ON MPIroot
IF(MPIroot)THEN !<<<<
  FileName = TRIM(IniFile)
  ! Get name of ini file
  WRITE(UNIT_StdOut,*)'| Reading from file "',TRIM(filename),'":'
  INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
  IF (.NOT.file_exists) THEN
    CALL Abort(__STAMP__,&
        "parameter file '"//TRIM(filename)//"' file does not exist.",TypeInfo="FileNotFoundError")
  END IF

  OPEN(NEWUNIT= iniUnit,        &
       FILE   = TRIM(filename), &
       STATUS = 'OLD',          &
       ACTION = 'READ',         &
       ACCESS = 'SEQUENTIAL',   &
       IOSTAT = stat)
  IF(stat.NE.0)THEN
    CALL abort(__STAMP__,&
      "Could not open parameter file '"//TRIM(filename)//"'.")
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
CALL par_BCast(nLines,0)
ALLOCATE(FileContent(nLines))

IF(MPIroot)THEN !<<<<
  !read file
  REWIND(iniUnit)
  READ(iniUnit,'(A)') FileContent
  CLOSE(iniUnit)
END IF !MPIroot !<<<<
!BROADCAST FileContent
CALL par_BCast(FileContent,0)
!#if MPI
!CALL MPI_BCAST(FileContent,LEN(FileContent)*nLines,MPI_CHARACTER,0,worldComm,iError) !<<<<
!#endif

NULLIFY(Str1,Str2)
DO i=1,nLines !<<<<
  IF(.NOT.ASSOCIATED(Str1)) CALL GetNewString(Str1)
  ! Read line from file
  Str=FileContent(i)
  ! Remove comments with "!"
  CALL Split(Str,Str,"!")
  ! Remove comments with "#"
  CALL Split(Str,Str,"#")
  Str=remove_blanks(Str)
  Str=Replace(Str,"(/","")
  Str=Replace(Str,"/)","")
  ! Replace brackets
  ! DO NOT Replace commas, used for array dimensions!
  !Str1%Str=Replace(Str1%Str,","," ")
  ! Lower case
  CALL LowCase(TRIM(Str),HelpStr)
  ! If we have a remainder (no comment only)
  IF(LEN_TRIM(HelpStr).GT.2) THEN
    Str1%Str=TRIM(HelpStr)
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
DEALLOCATE(FileContent)

!find line continuation "&" and merge strings (can be multiple lines)
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF(INDEX((Str1%str),'&').NE.0)THEN !found "&"
    CALL Split(Str1%Str,HelpStr,"&") !take part in front of "&"
    Str2=>Str1%nextStr
#if(!defined(NVHPC))
    DEALLOCATE(Str1%Str)
#endif /* ONLY NVHPC COMPILER DOES NOT SEEM TO WORK WITH ALLOCATABLE CHARACTERS (SIGSEV!) */
    Str1%Str=TRIM(HelpStr)//TRIM(Str2%Str)
    CALL deleteString(Str2)
    !do not go to next  string as long as there are "&" in the string
  ELSE
    Str1=>Str1%NextStr !nothing to be done
  END IF
END DO
!check len_trim<MAXLEN
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF(LEN_TRIM(Str1%Str).EQ.MAXLEN)THEN
    CALL abort(__STAMP__,&
      "parameter file readin: Line of input file might be longer than MAXLEN.",Intinfo=MAXLEN)
  END IF
  Str1=>Str1%NextStr !nothing to be done
END DO

ReadInDone = .TRUE.
CALL exit_subregion("read-parameterfile")
END SUBROUTINE FillStrings

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
#if (!defined(NVHPC))
DEALLOCATE(Str%Str)
#endif /* ONLY NVHPC COMPILER DOES NOT SEEM TO WORK WITH ALLOCATABLE CHARACTERS (SIGSEV!) */
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
TmpKey=remove_blanks(TmpKey)
Found=.FALSE.
Str1=>FirstString
DO WHILE(.NOT.Found)
  IF (.NOT.ASSOCIATED(Str1)) THEN
    IF (.NOT.PRESENT(Proposal)) THEN
      CALL abort(__STAMP__, &
           "missing necessary parameter '"//TRIM(TmpKey)//"'", &
           TypeInfo="MissingParameterError")
    ELSE ! Return default value
!      CALL LowCase(TRIM(Proposal),Str)
      IF(LEN_TRIM(Proposal).LE.LEN(Str))THEN
        Str=TRIM(Proposal)
      ELSE
        CALL abort(__STAMP__,&
          'parameter readin: proposal string of parameter '//TRIM(Key)//' does not fit into output string!')
      END IF


      IF (Str(1:1).NE.'@') THEN
        DefMsg='DEFAULT'
      END IF
      RETURN
    END IF ! (.NOT.PRESENT(Proposal))
  END IF ! (.NOT.ASSOCIATED(Str1))

  IF (INDEX(Str1%Str,TRIM(TmpKey)//'=').EQ.1) THEN
    Found=.TRUE.
    Str1%Str=replace(Str1%Str,TRIM(TmpKey)//'=',"")
    IF(LEN_TRIM(Str1%str).LE.LEN(Str))THEN
      Str=TRIM(Str1%Str)
    ELSE
      CALL abort(__STAMP__,&
        'parameter readin: string of parameter '//TRIM(Key)//' does not fit into output string!')
    END IF
    ! Remove string from list
    CALL DeleteString(Str1)
  ELSE
    ! Next string in list
    Str1=>Str1%NextStr
  END IF

END DO
END SUBROUTINE FindStr


!===================================================================================================================================
!> Transform upper case letters in "Str1" into lower case letters, result is "Str2", but only up the the equal sign.
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
!===================================================================================================================================
Str2=Str1
nLen=LEN_TRIM(Str1)
DO iLen=1,nLen
  IF(Str1(iLen:iLen).EQ.'=') EXIT ! Transformation stops at "="
  Upper=INDEX(UC,Str1(iLen:iLen))
  IF (Upper > 0) Str2(iLen:iLen) = lc(Upper:Upper)
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
  ProposalStr=TRIM(remove_blanks(str_tmp))
END SUBROUTINE ConvertToProposalStr

PURE FUNCTION remove_blanks(str_in) RESULT(str_out)
  IMPLICIT NONE
  !-------------------------------------------
  !input
  CHARACTER(LEN=*),INTENT(IN) :: str_in
  !output
  CHARACTER(LEN=LEN(str_in))  :: str_out
  !-------------------------------------------
  ! LOCAL VARIABLES
  INTEGER :: len_in,i,j
  !==============================================================================
  len_in=LEN_TRIM(str_in)
  str_out=""
  j=1
  DO i=1,len_in
    IF (str_in(i:i).NE.' ') THEN
      str_out(j:j)=str_in(i:i)
      j=j+1
    END If
  END DO
END FUNCTION remove_blanks

PURE FUNCTION replace(str_in,find,rep) RESULT(str_out)
  IMPLICIT NONE
  !-------------------------------------------
  ! input
  CHARACTER(LEN=*),INTENT(IN) :: str_in
  CHARACTER(LEN=*),INTENT(IN) :: find
  CHARACTER(LEN=*),INTENT(IN) :: rep
  ! output
  CHARACTER(LEN=LEN(str_in)) :: str_out
  !-------------------------------------------
  ! LOCAL VARIABLES
  CHARACTER(LEN=LEN(str_in)) :: str_tmp
  INTEGER :: i_find,lfind,lrep
  !=============================================================================
  str_out=""
  str_tmp=TRIM(str_in)
  i_find=INDEX(str_tmp,TRIM(find))
  lfind=LEN_TRIM(find)
  lrep=LEN_TRIM(rep)
  DO WHILE (i_find > 0)
    str_out=TRIM(str_out)//str_tmp(1:i_find-1)//TRIM(rep)
    str_tmp=str_tmp(i_find+lfind:)
    i_find=INDEX(str_tmp,TRIM(find))
  END DO
  str_out=TRIM(str_out)//TRIM(str_tmp)
END FUNCTION replace

SUBROUTINE split(str_in,bStr,separator)
  IMPLICIT NONE
  !-------------------------------------------
  ! input
  CHARACTER(LEN=*),INTENT(IN) :: str_in
  CHARACTER(LEN=1),INTENT(IN) :: separator
  ! output
  CHARACTER(LEN=*),INTENT(OUT) :: bStr
  !-------------------------------------------
  ! LOCAL VARIABLES
  INTEGER :: i_sep
  !==============================================================================
  bstr=TRIM(str_in)
  i_sep = INDEX(bstr,separator)
  IF (i_sep > 0) THEN
      bstr=bstr(1:i_sep-1)
  END IF
END SUBROUTINE split

FUNCTION count_sep(Key,str_in,separator) RESULT(n_sep)
  IMPLICIT NONE
  !-------------------------------------------
  ! input
  CHARACTER(LEN=*),INTENT(IN) :: Key
  CHARACTER(LEN=*),INTENT(IN) :: str_in
  CHARACTER(LEN=1),INTENT(IN) :: separator
  ! output
  INTEGER :: n_sep
  !-------------------------------------------
  ! LOCAL VARIABLES
  CHARACTER(LEN=LEN(str_in)) :: str_tmp
  INTEGER :: len_in,i
  !==============================================================================
  n_sep=0
  len_in=LEN_TRIM(str_in)
  str_tmp=TRIM(str_in)
  IF(str_tmp(1:1).EQ.separator) THEN
    CALL abort(__STAMP__,&
         "parameter '"//TRIM(Key)//"', problem with count separator:  first character should not be a separator!", &
         TypeInfo="InvalidParameterError")
  END IF
  DO i=2,len_in-1
    IF (str_tmp(i:i).EQ.separator) THEN
      n_sep=n_sep+1
    END IF
  END DO
  IF(str_tmp(len_in:len_in).EQ.separator) THEN
    CALL abort(__STAMP__,&
         "parameter '"//TRIM(Key)//"', problem with count separator: last character should not be a separator!", &
         TypeInfo="InvalidParameterError")
  END IF
END FUNCTION count_sep

END MODULE MODgvec_ReadInTools
