MODULE SEQ_FGMRES
        USE UTILITIES, ONLY: REALPARM
        USE FGMRES
        USE FGMRES_INTERFACES
        IMPLICIT NONE        
        PRIVATE
        TYPE SEQ_FGMRES_DATA
                INTEGER::Depth
                TYPE (FGMRES_DATA),POINTER::solvers(:)
        ENDTYPE

	PROCEDURE (MatrixVectorMult),POINTER::null_operator=>NULL()
	PROCEDURE (InformAboutIteration),POINTER::null_inform=>NULL()

        PUBLIC:: SEQ_FGMRES_DATA
        PUBLIC:: INIT_SEQ_FGMRES,SET_OPERATORS,SEQ_FGMRES_SOLVE
        PUBLIC:: SET_DP_INSTANCE
        PUBLIC:: DELETE_SEQ_FGMRES
CONTAINS
        SUBROUTINE INIT_SEQ_FGMRES(solver,N,M,depth,params,l)
                TYPE(SEQ_FGMRES_DATA),INTENT(INOUT)::solver
                INTEGER,INTENT(IN)::N,depth,M(depth)
                TYPE(GMRES_PARAMETERS),INTENT(IN)::params
                INTEGER(C_INTPTR_T),INTENT(OUT),OPTIONAL::l
                TYPE(FGMRES_DATA),POINTER::gsolver
                TYPE(C_PTR)::cptr
                INTEGER::I
                INTEGER(C_INTPTR_T)::length,full_length
                solver%Depth=depth
                ALLOCATE(solver%solvers(Depth))
                gsolver=>solver%solvers(Depth)
                CALL AllocateFGMRES(gsolver,N,M(depth),.FALSE.,full_length)
                gsolver%params=params
                gsolver%params%Maxit=1
                
      		NULLIFY(solver%solvers(Depth)%PrecondLeft)
	        solver%solvers(Depth)%MLeft=C_NULL_PTR
                DO I=Depth-1,1,-1
                        CALL  AllocateFGMRES(solver%solvers(I),N,&
                                &M(I),.TRUE.,length)
                        full_length=full_length+length
                        cptr=C_LOC(solver%solvers(I+1))
        		solver%solvers(I)%PrecondRight=>FGMRES_BASED_PRECONDITIONER
		        solver%solvers(I)%MRight=cptr
		        solver%solvers(I)%MLeft=C_NULL_PTR
        		NULLIFY(solver%solvers(I)%PrecondLeft)
		        solver%solvers(I)%params=params
		        solver%solvers(I)%params%MaxIt=1
                ENDDO
		solver%solvers(1)%params%MaxIt=params%MaxIt
                IF (PRESENT(l)) l=full_length
        ENDSUBROUTINE
        SUBROUTINE SET_OPERATORS(solver,apply,dp,inform,A)
                TYPE(SEQ_FGMRES_DATA),INTENT(IN)::solver
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::apply
		PROCEDURE(MANYDOTPRODUCTS),POINTER,INTENT(IN)::dp
		PROCEDURE (InformAboutIteration),POINTER,INTENT(IN)::Inform
		TYPE(C_PTR),INTENT(IN)::A
                TYPE(FGMRES_DATA),POINTER::gsolver
                INTEGER::I,Depth
                Depth=solver%Depth
                gsolver=>solver%solvers(Depth)
                CALL Set_FGMRES_Operators(gsolver,apply,null_operator,&
                        &null_operator,dp,null_inform,&
                        &A,C_NULL_PTR,C_NULL_PTR)

                DO I=solver%Depth-1,1,-1
                        solver%solvers(I)%ApplyOperator=>apply
                        NULLIFY(solver%solvers(I)%PrecondLeft)
                        solver%solvers(I)%A=A
                        solver%solvers(I)%MLeft=C_NULL_PTR
                        solver%solvers(I)%MANYDP=>dp
                        solver%solvers(I)%Inform=>null_inform
                ENDDO
                solver%solvers(1)%Inform=>inform
        ENDSUBROUTINE
        SUBROUTINE SET_RIGHT_PRECONDTIONER(solver,precond_right,M2)
                TYPE(SEQ_FGMRES_DATA),INTENT(IN)::solver
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::precond_right
		TYPE(C_PTR),INTENT(IN)::M2
                TYPE(FGMRES_DATA),POINTER::gsolver
                INTEGER::I,Depth
                Depth=solver%Depth
                gsolver=>solver%solvers(Depth)
                IF (.NOT. ASSOCIATED(precond_right)) NULLIFY(gsolver%PrecondRight)
                gsolver%Mright=M2
        ENDSUBROUTINE
        SUBROUTINE SET_LEFT_PRECONDITIONER(solver,precond_left,M1)
                TYPE(SEQ_FGMRES_DATA),INTENT(IN)::solver
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::precond_left
		TYPE(C_PTR),INTENT(IN)::M1
                TYPE(FGMRES_DATA),POINTER::gsolver
                INTEGER::I,Depth

                DO I=solver%Depth,1,-1
                        solver%solvers(I)%PrecondLeft=>precond_left
                        IF (.NOT. ASSOCIATED(precond_left)) NULLIFY(solver%solvers(I)%PrecondLeft)
                        solver%solvers(I)%MLeft=M1
                ENDDO
        ENDSUBROUTINE
        SUBROUTINE SET_DP_INSTANCE(solver,dpinstance)
                TYPE(SEQ_FGMRES_DATA),INTENT(IN)::solver
		TYPE(C_PTR),INTENT(IN)::dpinstance
                INTEGER::I,Depth
                Depth=solver%Depth
                DO I=1,Depth
                        solver%solvers(I)%dpinstance=dpinstance
                ENDDO
        ENDSUBROUTINE
        SUBROUTINE SEQ_FGMRES_SOLVE(solver,x,x0,b,info)
                TYPE(SEQ_FGMRES_DATA),INTENT(INOUT)::solver
		COMPLEX(REALPARM),INTENT(INOUT)::x(:)
		COMPLEX(REALPARM),INTENT(IN)::x0(:),b(:)
		TYPE(RESULT_INFO),INTENT(OUT)::info
                CALL GMRES_SOLVE(solver%solvers(1),x,x0,b,info)
        ENDSUBROUTINE
        SUBROUTINE DELETE_SEQ_FGMRES(solver)
                TYPE(SEQ_FGMRES_DATA),INTENT(INOUT)::solver
                INTEGER::I
                DO I=solver%Depth,1,-1
                        CALL DeallocateFGMRES(solver%solvers(I))
                ENDDO
                DEALLOCATE(solver%solvers)
                solver%Depth=-1
        ENDSUBROUTINE

        RECURSIVE    SUBROUTINE FGMRES_BASED_PRECONDITIONER(psolver,v_in,v_out)
		TYPE(C_PTR),INTENT(IN)::psolver
		COMPLEX(REALPARM),INTENT(IN)::v_in(:)
		COMPLEX(REALPARM),INTENT(INOUT)::v_out(:)
		TYPE(FGMRES_DATA),POINTER::solver
                COMPLEX(REALPARM),TARGET::v0(SIZE(v_out))
                TYPE(RESULT_INFO)::info
                INTEGER::I,J
                CALL C_F_POINTER(psolver,solver)
                CALL GMRES_SOLVE(solver,v_out,v_in,v_in,info)
        ENDSUBROUTINE
ENDMODULE
