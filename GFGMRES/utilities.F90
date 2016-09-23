!Copyright (c) 2016 Mikhail Kruglyakov 
!This file is part of GFGMRES.

!GFGMRES is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 2 of the License.

!GFGMRES is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with GFMRES.  If not, see <http://www.gnu.org/licenses/>.

MODULE UTILITIES
	USE, INTRINSIC :: iso_c_binding
        IMPLICIT NONE
	INTEGER, PARAMETER ::REALPARM = SELECTED_REAL_KIND(15)

	INTEGER,PARAMETER::ZERO=0
	INTEGER,PARAMETER::ONE=1

	COMPLEX(REALPARM),PARAMETER ::	C_ZERO=(0e0_REALPARM,0e0_REALPARM)
	COMPLEX(REALPARM),PARAMETER ::	C_ONE=(1e0_REALPARM,0e0_REALPARM)
	REAL(REALPARM)   ,PARAMETER ::	R_ZERO=0e0_REALPARM
	REAL(REALPARM)   ,PARAMETER ::	R_ONE=1e0_REALPARM

        INTERFACE FREE
                MODULE PROCEDURE FreeMatrix,FreeRealVector,FreeComplexVector

        ENDINTERFACE


      CONTAINS
        

        SUBROUTINE FreeComplexVector(v)
                COMPLEX(REALPARM),POINTER,INTENT(INOUT)::v(:)
                IF (ASSOCIATED(v)) THEN
                        DEALLOCATE(v)
                        NULLIFY(v)
                ENDIF
        ENDSUBROUTINE
        SUBROUTINE FreeRealVector(v)
                REAL(REALPARM),POINTER,INTENT(INOUT)::v(:)
                IF (ASSOCIATED(v)) THEN
                        DEALLOCATE(v)
                        NULLIFY(v)
                ENDIF
        ENDSUBROUTINE
        SUBROUTINE FreeMatrix(m)
                COMPLEX(REALPARM),POINTER,INTENT(INOUT)::m(:,:)
                IF (ASSOCIATED(m)) THEN
                        DEALLOCATE(m)
                        NULLIFY(m)
                ENDIF
        ENDSUBROUTINE


      END MODULE

