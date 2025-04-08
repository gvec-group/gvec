!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
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
MODULE MODgvec_Output_VTK
! MODULES
USE MODgvec_Globals, ONLY: wp
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
USE MODgvec_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim1                    !! dimension of the data (either 1:lines,2=quads or 3=hexas)
INTEGER,INTENT(IN)            :: vecdim                  !! dimension of coordinates
INTEGER,INTENT(IN)            :: nVal                    !! Number of nodal output variables
INTEGER,INTENT(IN)            :: NPlot(dim1)             !! Number of output points per element : (nPlot+1)**dim1
INTEGER,INTENT(IN)            :: nElems                  !! Number of output elements
REAL(wp),INTENT(IN)           :: Coord(vecdim,1:PRODUCT(Nplot+1),nElems)      ! CoordinatesVector
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          !! Names of all variables that will be written out
REAL(wp),INTENT(IN)           :: Values(nVal,1:PRODUCT(Nplot+1),nElems)   !! Statevector
CHARACTER(LEN=*),INTENT(IN)   :: FileString              !! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER     :: kindFloat=8  !set floating point accuracy single (4) double (8), should be equal or lower than input data!
REAL(KIND=kindFloat)  :: FLOATdummy
CHARACTER(LEN=7)      :: strfloat
INTEGER               :: INTdummy
INTEGER               :: sizefloat,sizeInt
INTEGER            :: i,j,k,iVal,iElem,Offset,nBytes,nVTKPoints,nVTKCells,ivtk
INTEGER            :: Vertex(2**dim1,PRODUCT(Nplot)*nElems)  ! ?
INTEGER            :: ProdNplot,ProdNplot_p1,NPlot_p1(dim1),CellID,PointID,ElemType  ! ?
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2  ! ?
CHARACTER(LEN=300) :: Buffer
CHARACTER(LEN=255) :: tmpVarName,tmpVarNameY,tmpVarNameZ
INTEGER            :: StrLen,iValVec,nValVec,VecOffset(0:nVal)
LOGICAL            :: isVector,maybeVector
CHARACTER(LEN=1)   :: strvecdim
CHARACTER(LEN=1)   :: lf
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITE DATA TO VTX XML BINARY (VTU) FILE "'//TRIM(FileString)//'" ...'
ivtk=GETFREEUNIT()
NPlot_p1  =(Nplot(:)+1)
ProdNPlot  =PRODUCT(Nplot(:))
ProdNPlot_p1  =PRODUCT(Nplot_p1(:))

IF(kindFloat.EQ.4) THEN
  strfloat='Float32'
ELSEIF(kindFloat.EQ.8)THEN
  strfloat='Float64'
ELSE
  STOP 'kindFloat not implemented in output vtk'
END IF
sizefloat=SIZEOF_F(FLOATdummy)
sizeInt  =SIZEOF_F(INTdummy)


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
nVTKPoints= ProdNPlot_p1*nElems
nVTKCells = ProdNPlot   *nElems
Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKPoints
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
    Buffer='        <DataArray type="'//strfloat//'" Name="'//TRIM(tmpVarName)//'" NumberOfComponents="'//strvecdim// &
           &'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    Offset=Offset+sizeInt+vecdim*nVTKPoints*sizefloat
    WRITE(StrOffset,'(I16)')Offset
    VecOffset(iValVec)=VecOffset(iValVec-1)+vecdim
    iVal=iVal+vecdim-1 !skip the Y (& Z) components
  ELSE
    Buffer='        <DataArray type="'//strfloat//'" Name="'//TRIM(tmpVarName)// &
           &'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    Offset=Offset+sizeInt+nVTKPoints*sizeFloat
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
Buffer='        <DataArray type="'//strfloat//'" Name="Coordinates" NumberOfComponents="'//strvecdim// &
'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+sizeInt+vecdim*nVTKPoints*sizeFloat
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+sizeInt+2**dim1*nVTKCells*sizeInt
WRITE(StrOffset,'(I16)')Offset
! Offset in connectivity data
Buffer='        <DataArray type="Int32" Name="offsets" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+sizeInt+nVTKCells*sizeInt
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
nBytes = nVTKPoints*sizeFloat
DO iValVec=1,nValVec
  WRITE(ivtk) (vecOffset(iValVec)-vecOffset(iValVec-1))*nBytes
  WRITE(ivtk) REAL(Values(VecOffSet(iValVec-1)+1:VecOffset(iValVec),:,:),kindFloat)
END DO !iValVec
! Point coordinates
nBytes = nVTKPoints * vecdim*sizeFloat
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coord(:,:,:),kindFloat)
! Connectivity
SELECT CASE(dim1)
CASE(1)
  CellID = 0
  PointID= 0
  DO iElem=1,nElems
    DO i=1,NPlot(1)
      CellID = CellID+1
      !visuLineElem
      Vertex(:,CellID) = (/ PointID+(i-1), PointID+ i /)
    END DO
    PointID=PointID+NPlot(1)
  END DO
CASE(2)
  CellID = 0
  PointID= 0
  DO iElem=1,nElems
    DO j=1,NPlot(2)
      DO i=1,NPlot(1)
        CellID = CellID+1
        !visuQuadElem
        Vertex(:,CellID) = (/                  &
          PointID+(i-1)+  j   * NPlot_p1(1) ,    & !P4
          PointID+(i-1)+ (j-1)* NPlot_p1(1) ,    & !P1(CGNS=tecplot standard)
          PointID+ i   + (j-1)* NPlot_p1(1) ,    & !P2
          PointID+ i   +  j   * NPlot_p1(1)     /) !P3
      END DO
    END DO
    PointID=PointID+ProdNPlot_p1
  END DO
CASE(3)
  CellID=0
  PointID=0
  DO iElem=1,nElems
    DO k=1,NPlot(3)
      DO j=1,NPlot(2)
        DO i=1,NPlot(1)
          CellID=CellID+1
          !
          Vertex(:,CellID)=(/                                       &
            PointID+(i-1)+( j   +(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P4(CGNS=tecplot standard)
            PointID+(i-1)+((j-1)+(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P1
            PointID+ i   +((j-1)+(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P2
            PointID+ i   +( j   +(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P3
            PointID+(i-1)+( j   + k   *NPlot_p1(2))*NPlot_p1(1),      & !P8
            PointID+(i-1)+((j-1)+ k   *NPlot_p1(2))*NPlot_p1(1),      & !P5
            PointID+ i   +((j-1)+ k   *NPlot_p1(2))*NPlot_p1(1),      & !P6
            PointID+ i   +( j   + k   *NPlot_p1(2))*NPlot_p1(1)      /) !P7
        END DO
      END DO
    END DO
    !
    PointID=PointID+ProdNPlot_p1
  END DO
END SELECT
nBytes = 2**dim1*nVTKCells*sizeInt
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
! Offset in connectivity
nBytes = nVTKCells*sizeInt
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=2**dim1,2**dim1*nVTKCells,2**dim1)
! Elem type
SELECT CASE(dim1)
CASE(1)
  ElemType =3 !VTK_LINE
CASE(2)
  ElemType =9 !VTK_QUAD
CASE(3)
  ElemType =12  !VTK_HEXAHEDRON
END SELECT
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToVTK



END MODULE MODgvec_Output_VTK
