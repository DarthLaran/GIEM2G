MODULE  FGMRES_INTERFACES
     USE UTILITIES, ONLY : REALPARM 
     USE, INTRINSIC :: iso_c_binding

     IMPLICIT NONE
      LOGICAL, PARAMETER::FLEXIBLE=.TRUE.
      LOGICAL, PARAMETER::TRUE_RESIDUAL_AT_RESTART=.TRUE.
      LOGICAL, PARAMETER::ITERATIVE_RESIDUAL_AT_RESTART=.FALSE.

      INTEGER,PARAMETER::CONVERGED=1
      INTEGER,PARAMETER::NON_TRUE_CONVERGED=2

      INTEGER,PARAMETER::NON_CONVERGED=-1
      INTEGER,PARAMETER::INPROCESS=-2
      INTEGER,PARAMETER::INTERUPTED=-3

      INTEGER,PARAMETER::INTERUPT_SOLVER=1  
      INTEGER,PARAMETER::NOINTERUPT=-1
      TYPE RESULT_INFO
	      INTEGER:: Iterations
	      REAL(REALPARM)::bea
	      REAL(REALPARM)::be
	      INTEGER:: stat
      ENDTYPE

      TYPE GMRES_PARAMETERS
	      INTEGER:: MaxIt
	      REAL(REALPARM)::tol
	      LOGICAL::RESTART_RESIDUAL
      ENDTYPE
        INTERFACE 
                SUBROUTINE MatrixVectorMult(A,v_in,v_out)
                        USE UTILITIES,ONLY: REALPARM
                        USE, INTRINSIC :: iso_c_binding
                        TYPE(C_PTR),INTENT(IN)::A
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::v_in(:)
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::v_out(:)
                ENDSUBROUTINE

                SUBROUTINE MANYDOTPRODUCTS(M,v,N,res,ptr)
                        USE UTILITIES,ONLY: REALPARM
                        USE, INTRINSIC :: iso_c_binding
                        COMPLEX(REALPARM), POINTER, INTENT(IN):: M(:,:)
                        COMPLEX(REALPARM), POINTER, INTENT(IN):: v(:)
                        INTEGER                   , INTENT(IN):: N
                        COMPLEX(REALPARM), POINTER, INTENT(IN):: res(:)
                        TYPE(C_PTR)::ptr
                ENDSUBROUTINE

                FUNCTION InformAboutIteration(info) RESULT(interp)
                        IMPORT RESULT_INFO
                        TYPE(RESULT_INFO),INTENT(IN)::info
                        INTEGER::interp
                        !Possibility to interupt solver after analyze of info
                ENDFUNCTION

        END INTERFACE
        PRIVATE:: REALPARM
ENDMODULE
