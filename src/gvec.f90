!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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


!===================================================================================================================================
!> 
!!# **GVEC** Driver program 
!!
!===================================================================================================================================
PROGRAM GVEC
USE MOD_Globals    ,ONLY: fmt_sep
USE MOD_MHDEQ      ,ONLY: InitMHDEQ,FinalizeMHDEQ
USE MOD_ReadInTools,ONLY: IgnoredStrings 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: nArgs
CHARACTER(LEN=255)      :: Parameterfile
!===================================================================================================================================
nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs.EQ.1)THEN
  CALL GET_COMMAND_ARGUMENT(1,Parameterfile)
ELSE
  STOP 'parameterfile not given, usage: "executable parameter.ini"'
END IF


WRITE(*,'(132("*"))')
WRITE(*,'(("*"),130X,("*"))')
WRITE(*,'(("*"),13X,A104,13X,("*"))')'         GGGGGGGGGGGGG    VVVVVVVV           VVVVVVVV    EEEEEEEEEEEEEEEEEEEEEE           CCCCCCCCCCCCC '  
WRITE(*,'(("*"),13X,A104,13X,("*"))')'      GGG::::::::::::G    V::::::V           V::::::V    E::::::::::::::::::::E        CCC::::::::::::C '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'    GG:::::::::::::::G    V::::::V           V::::::V    E::::::::::::::::::::E      CC:::::::::::::::C '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'   G:::::GGGGGGGG::::G    V::::::V           V::::::V    EE::::::EEEEEEEEE::::E     C:::::CCCCCCCC::::C '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'  G:::::G       GGGGGG     V:::::V           V:::::V       E:::::E       EEEEEE    C:::::C       CCCCCC '
WRITE(*,'(("*"),13X,A104,13X,("*"))')' G:::::G                    V:::::V         V:::::V        E:::::E                C:::::C               '
WRITE(*,'(("*"),13X,A104,13X,("*"))')' G:::::G                     V:::::V       V:::::V         E::::::EEEEEEEEEE      C:::::C               '
WRITE(*,'(("*"),13X,A104,13X,("*"))')' G:::::G    GGGGGGGGGG        V:::::V     V:::::V          E:::::::::::::::E      C:::::C               '
WRITE(*,'(("*"),13X,A104,13X,("*"))')' G:::::G    G::::::::G         V:::::V   V:::::V           E:::::::::::::::E      C:::::C               '
WRITE(*,'(("*"),13X,A104,13X,("*"))')' G:::::G    GGGGG::::G          V:::::V V:::::V            E::::::EEEEEEEEEE      C:::::C               '
WRITE(*,'(("*"),13X,A104,13X,("*"))')' G:::::G        G::::G           V:::::V:::::V             E:::::E                C:::::C               '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'  G:::::G       G::::G            V:::::::::V              E:::::E       EEEEEE    C:::::C       CCCCCC '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'   G:::::GGGGGGGG::::G             V:::::::V             EE::::::EEEEEEEE:::::E     C:::::CCCCCCCC::::C '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'    GG:::::::::::::::G              V:::::V              E::::::::::::::::::::E      CC:::::::::::::::C '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'      GGG::::::GGG:::G               V:::V               E::::::::::::::::::::E        CCC::::::::::::C '
WRITE(*,'(("*"),13X,A104,13X,("*"))')'         GGGGGG   GGGG                VVV                EEEEEEEEEEEEEEEEEEEEEE           CCCCCCCCCCCCC '
WRITE(*,'(("*"),130X,("*"))')
WRITE(*,'(132("*"))')
WRITE(*,*)


!initialization phase
CALL InitMHDEQ()

CALL IgnoredStrings()
! do something

!finalization phase

CALL FinalizeMHDEQ()

WRITE(*,fmt_sep)
WRITE(*,'(A)') ' GVEC SUCESSFULLY FINISHED'
WRITE(*,fmt_sep)

END PROGRAM GVEC


