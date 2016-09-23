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
                        COMPLEX(REALPARM),INTENT(IN)::v_in(:)
                        COMPLEX(REALPARM),INTENT(INOUT)::v_out(:)
                ENDSUBROUTINE

                SUBROUTINE MANYDOTPRODUCTS(M,v,N,res,ptr)
                        USE UTILITIES,ONLY: REALPARM
                        USE, INTRINSIC :: iso_c_binding
                        COMPLEX(REALPARM), INTENT(IN):: M(:,:)
                        COMPLEX(REALPARM), INTENT(IN):: v(:)
                        INTEGER          , INTENT(IN):: N
                        COMPLEX(REALPARM),  INTENT(INOUT):: res(:)
                        TYPE(C_PTR),INTENT(IN)::ptr
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
