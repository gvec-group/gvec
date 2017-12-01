!===================================================================================================================================
! Copyright (c) 2017 - 2018 Florian Hindenlang <hindenlang@gmail.com>
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/hopr)
!
! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
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
#include "defines.h"
#define SIZEOF_F(x) STORAGE_SIZE(x)/8

!===================================================================================================================================
!>
!!# Module **Output VTK**
!!
!! Write to unstructured VTK file 
!!
!===================================================================================================================================
MODULE MOD_Output_VTK
! MODULES
USE MOD_Globals, ONLY: wp
IMPLICIT NONE
PRIVATE

!INTERFACE WriteDataToVTK
!  MODULE PROCEDURE WriteDataToVTK
!END INTERFACE

PUBLIC::WriteDataToVTK
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write 3D point data to VTK format
!!
!===================================================================================================================================
SUBROUTINE WriteDataToVTK(dim1,vecDim,nVal,NPlot,nElems,VarNames,Coord,Values,FileString)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim1                    !! dimension of the data (either 2=quads or 3=hexas)
INTEGER,INTENT(IN)            :: vecdim                  !! dimension of coordinates 
INTEGER,INTENT(IN)            :: nVal                    !! Number of nodal output variables
INTEGER,INTENT(IN)            :: NPlot(dim1)             !! Number of output points per element : (nPlot+1)**dim1
INTEGER,INTENT(IN)            :: nElems                  !! Number of output elements
REAL,INTENT(IN)               :: Coord(vecdim,1:PRODUCT(Nplot+1),nElems)      ! CoordinatesVector 
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          !! Names of all variables that will be written out
REAL,INTENT(IN)               :: Values(nVal,1:PRODUCT(Nplot+1),nElems)   !! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString              !! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iVal,iElem,Offset,nBytes,nVTKElems,nVTKCells,ivtk=44
INTEGER            :: INTdummy
REAL(KIND=8)       :: FLOAT64dummy
INTEGER            :: Vertex(2**dim1,PRODUCT(Nplot+1)*nElems)  ! ?
INTEGER            :: ProdNplot,ProdNplot_p1,NPlot_p1(dim1),NodeID,NodeIDElem,ElemType  ! ?
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2  ! ?
CHARACTER(LEN=300) :: Buffer
CHARACTER(LEN=255) :: tmpVarName,tmpVarNameY,tmpVarNameZ
INTEGER            :: StrLen,iValVec,nValVec,VecOffset(0:nVal)
LOGICAL            :: isVector,maybeVector
CHARACTER(LEN=1)   :: strvecdim
CHARACTER(LEN=1)   :: lf
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITE DATA TO VTX XML BINARY (VTU) FILE "'//TRIM(FileString)//'" ...'
NPlot_p1  =(Nplot(:)+1)
ProdNPlot  =PRODUCT(Nplot(:))
ProdNPlot_p1  =PRODUCT(Nplot_p1(:))

IF(vecdim.LT.dim1) THEN
  WRITE(*,*)'WARNING:data dimension should be <= vecdim! dim1= ',dim1,' vecdim= ',vecdim
  STOP
END IF
! Line feed character
lf = char(10)
WRITE(strvecdim,'(I1)') vecdim

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify file type
nVTKElems=ProdNPlot_p1*nElems
nVTKCells=ProdNPlot*nElems
Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKElems
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
       &NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=0
WRITE(StrOffset,'(I16)')Offset
!accout for vectors: 
! if Variable Name ends with an X and the following have the same name with Y and Z 
! then it forms a vector variable (X is omitted for the name) 

iVal=0 !scalars
iValVec=0 !scalars & vectors
VecOffset(0)=0
DO WHILE(iVal.LT.nVal)
  iVal=iVal+1
  iValVec=iValVec+1
  tmpVarName=TRIM(VarNames(iVal)) 
  StrLen=LEN(TRIM(tmpVarName))
  maybeVector=(iVal+vecdim-1.LE.nVal)
  isVector=.FALSE.
  IF(maybeVector)THEN
    SELECT CASE(vecdim)
    CASE(2)
      tmpVarNameY=TRIM(VarNames(iVal+1))
      isVector=((iVal+2.LE.nVal).AND.(INDEX(tmpVarName( StrLen:StrLen),"X").NE.0) &
                                .AND.(INDEX(tmpVarNameY(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Y").NE.0))
    CASE(3)
      tmpVarNameY=TRIM(VarNames(iVal+1))
      tmpVarNameZ=TRIM(VarNames(iVal+2)) 
      isVector=((iVal+2.LE.nVal).AND.(INDEX(tmpVarName( StrLen:StrLen),"X").NE.0) &
                                .AND.(INDEX(tmpVarNameY(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Y").NE.0) &
                                .AND.(INDEX(tmpVarNameZ(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Z").NE.0))

    END SELECT
  END IF !maybevector

  IF(isvector)THEN !variable is a vector!
    tmpVarName=tmpVarName(:StrLen-1)
    Buffer='        <DataArray type="Float64" Name="'//TRIM(tmpVarName)//'" NumberOfComponents="'//strvecdim// &
           &'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    Offset=Offset+SIZEOF_F(INTdummy)+vecdim*nVTKElems*SIZEOF_F(FLOAT64dummy)
    WRITE(StrOffset,'(I16)')Offset
    VecOffset(iValVec)=VecOffset(iValVec-1)+vecdim
    iVal=iVal+vecdim-1 !skip the Y (& Z) components
  ELSE
    Buffer='        <DataArray type="Float64" Name="'//TRIM(tmpVarName)// &
           &'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    Offset=Offset+SIZEOF_F(INTdummy)+nVTKElems*SIZEOF_F(FLOAT64dummy)
    WRITE(StrOffset,'(I16)')Offset
    VecOffset(iValVec)=VecOffset(iValVec-1)+1
  END IF !isvector
END DO !iVal <=nVal
nValVec=iValVec

Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify cell data
Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="'//strvecdim// &
'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF_F(INTdummy)+vecdim*nVTKElems*SIZEOF_F(FLOAT64dummy)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF_F(INTdummy)+2**dim1*nVTKElems*SIZEOF_F(INTdummy)
WRITE(StrOffset,'(I16)')Offset
! Offsets
Buffer='        <DataArray type="Int32" Name="offsets" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF_F(INTdummy)+nVTKElems*SIZEOF_F(INTdummy)
WRITE(StrOffset,'(I16)')Offset
! Elem types
Buffer='        <DataArray type="Int32" Name="types" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
! Prepare append section
Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
! Write leading data underscore
Buffer='_';WRITE(ivtk) TRIM(Buffer)

! Write binary raw data into append section
! Point data
nBytes = nVTKElems*SIZEOF_F(FLOAT64dummy)
DO iValVec=1,nValVec
  WRITE(ivtk) (vecOffset(iValVec)-vecOffset(iValVec-1))*nBytes, &
              REAL(Values(VecOffSet(iValVec-1)+1:VecOffset(iValVec),:,:),8)
END DO !iValVec
! Points
nBytes = nBytes * vecdim
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coord(:,:,:),8)
! Connectivity
SELECT CASE(dim1)
CASE(2)
  NodeID = 0
  NodeIDElem = 0
  DO iElem=1,nElems
    DO j=1,NPlot(2)
      DO i=1,NPlot(1)
        NodeID = NodeID+1
        !visuQuadElem
        Vertex(:,NodeID) = (/                  &
          NodeIDElem+(i-1)+  j   * NPlot_p1(1) ,    & !P4
          NodeIDElem+(i-1)+ (j-1)* NPlot_p1(1) ,    & !P1(CGNS=tecplot standard)
          NodeIDElem+ i   + (j-1)* NPlot_p1(1) ,    & !P2
          NodeIDElem+ i   +  j   * NPlot_p1(1)     /) !P3
      END DO
    END DO
    NodeIDElem=NodeIDElem+ProdNPlot_p1
  END DO
CASE(3)
  NodeID=0
  NodeIDElem=0
  DO iElem=1,nElems
    DO k=1,NPlot(3)
      DO j=1,NPlot(2)
        DO i=1,NPlot(1)
          NodeID=NodeID+1
          !
          Vertex(:,NodeID)=(/                                       &
            NodeIDElem+(i-1)+( j   +(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P4(CGNS=tecplot standard)
            NodeIDElem+(i-1)+((j-1)+(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P1
            NodeIDElem+ i   +((j-1)+(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P2
            NodeIDElem+ i   +( j   +(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P3
            NodeIDElem+(i-1)+( j   + k   *NPlot_p1(2))*NPlot_p1(1),      & !P8
            NodeIDElem+(i-1)+((j-1)+ k   *NPlot_p1(2))*NPlot_p1(1),      & !P5
            NodeIDElem+ i   +((j-1)+ k   *NPlot_p1(2))*NPlot_p1(1),      & !P6
            NodeIDElem+ i   +( j   + k   *NPlot_p1(2))*NPlot_p1(1)      /) !P7
        END DO
      END DO
    END DO
    !
    NodeIDElem=NodeIDElem+ProdNPlot_p1
  END DO
END SELECT
nBytes = 2**dim1*nVTKElems*SIZEOF_F(INTdummy)
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
! Offset
nBytes = nVTKElems*SIZEOF_F(INTdummy)
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=2**dim1,2**dim1*nVTKElems,2**dim1)
! Elem type
ElemType =3+3*dim1 !9 VTK_QUAD 12  VTK_HEXAHEDRON
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKElems)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToVTK
 
 

END MODULE MOD_Output_VTK
