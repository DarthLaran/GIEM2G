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
        PUBLIC:: DELETE_SEQ_FGMRES
CONTAINS
        SUBROUTINE INIT_SEQ_FGMRES(solver,N,M,depth,params,l)
                TYPE(SEQ_FGMRES_DATA),INTENT(INOUT)::solver
                INTEGER,INTENT(IN)::N,M(depth),depth
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
                gsolver%params%Maxit=1
                
                DO I=Depth-1,1,-1
                        CALL  AllocateFGMRES(solver%solvers(I),N,&
                                &M(I),.TRUE.,length)
                        full_length=full_length+length
                        cptr=C_LOC(solver%solvers(I+1))
        		solver%solvers(I)%PrecondRight=>FGMRES_BASED_PRECONDITIONER
		        solver%solvers(I)%MRight=cptr
		        solver%solvers(I)%params=params
		        solver%solvers(I)%params%MaxIt=1
                ENDDO
		solver%solvers(1)%params%MaxIt=params%MaxIt
                IF (PRESENT(l)) l=full_length
        ENDSUBROUTINE
        SUBROUTINE SET_OPERATORS(solver,apply,precond_right,precond_left,dp,inform,A,M1,M2)
                TYPE(SEQ_FGMRES_DATA),INTENT(IN)::solver
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::apply
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::precond_right
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::precond_left
		PROCEDURE(MANYDOTPRODUCTS),POINTER,INTENT(IN)::dp
		PROCEDURE (InformAboutIteration),POINTER,INTENT(IN)::Inform
		TYPE(C_PTR),INTENT(IN)::A,M1,M2
                TYPE(FGMRES_DATA),POINTER::gsolver
                INTEGER::I,Depth
                Depth=solver%Depth
                gsolver=>solver%solvers(Depth)
                CALL Set_FGMRES_Operators(gsolver,apply,precond_right,precond_left,dp,null_inform,A,M1,M2)

                DO I=solver%Depth-1,1,-1
                        solver%solvers(I)%ApplyOperator=>apply
                        solver%solvers(I)%PrecondLeft=>precond_left
                        IF (.NOT. ASSOCIATED(precond_left)) NULLIFY(solver%solvers(I)%PrecondLeft)
                        solver%solvers(I)%A=A
                        solver%solvers(I)%MLeft=M1
                        solver%solvers(I)%MANYDP=>dp
                        solver%solvers(I)%Inform=>null_inform
                ENDDO
                        solver%solvers(1)%Inform=>inform
        ENDSUBROUTINE

        SUBROUTINE SEQ_FGMRES_SOLVE(solver,x,x0,b,info)
                TYPE(SEQ_FGMRES_DATA),INTENT(IN)::solver
		COMPLEX(REALPARM),POINTER,INTENT(IN)::x(:),x0(:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::b(:)
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
		COMPLEX(REALPARM),POINTER,INTENT(IN)::v_in(:),v_out(:)
		TYPE(FGMRES_DATA),POINTER::solver
                COMPLEX(REALPARM),TARGET::v0(SIZE(v_out))
                TYPE(RESULT_INFO)::info
                INTEGER::I,J
                CALL C_F_POINTER(psolver,solver)
                v0=(0e0_REALPARM,0e0_REALPARM)
                CALL GMRES_SOLVE(solver,v_out,v0,v_in,info)
        ENDSUBROUTINE
ENDMODULE
