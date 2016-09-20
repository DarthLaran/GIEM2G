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
                MODULE PROCEDURE::FreeMatrix,FreeRealVector,FreeComplexVector

        ENDINTERFACE
        INTERFACE
                FUNCTION DZNRM2(N,x,INCX) RESULT(RES)
                        IMPORT REALPARM
                        INTEGER,INTENT(IN)::N,INCX
                        COMPLEX(REALPARM)::x(:)
                        REAL(REALPARM)::RES
                ENDFUNCTION

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

