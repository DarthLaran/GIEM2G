MODULE IE_SOLVER_MODULE
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE
	USE DATA_TYPES_MODULE
	USE INTEGRAL_EQUATION_MODULE
	USE Timer_Module 
        USE LOGGER_MODULE
        USE APPLY_IE_OPERATOR_MODULE 

        USE FGMRES_INTERFACES
        USE SEQ_FGMRES

	IMPLICIT NONE
	PRIVATE

#define no_compile
!For integral equation kernel!

	PUBLIC:: SolveEquation,SetAnomalySigma
	
	CONTAINS

	SUBROUTINE SolveEquation(ie_op,fgmres_ctl,E_bkg,Esol,E0)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		TYPE (FGMRES_CTL_TYPE),INTENT(IN)::fgmres_ctl
		COMPLEX(REALPARM),POINTER,INTENT(IN)::E_bkg(:,:,:,:),Esol(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN),OPTIONAL::E0(:,:,:,:) !
		COMPLEX(REALPARM),POINTER::guess(:,:,:,:)
		INTEGER::anom_shape(3)
		INTEGER::tmp_shape(3)
		INTEGER::Iz,Nz,Nx,Ny_loc
		INTEGER::Ix,Iy,Ic
		INTEGER(MPI_CTL_KIND)::IERROR
		COMPLEX(REALPARM)::asiga
		REAL(DOUBLEPARM)::time1,time2
                CALL PRINT_BORDER
	        CALL LOGGER('IE solving starts')
		time1=GetTime()
		Nx=ie_op%Nx
		Nz=ie_op%Nz
		Ny_loc=ie_op%Ny_loc
		IF (ie_op%real_space) THEN
                        ALLOCATE(guess(ie_op%Nx,ie_op%Ny_loc,ie_op%Nz,3))
			DO Iz=1,Nz
				Esol(:,:,Iz,:)=ie_op%sqsigb(Iz)*E_bkg(:,:,Iz,:)
			ENDDO
			IF (PRESENT(E0)) THEN
			    DO Iy=1,Ny_loc
				    DO Ix=1,ie_op%Nx
					    DO Ic=1,3
						    DO Iz=1,Nz
							    asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
							    guess(Ix,Iy,Iz,Ic)=E0(Ix,Iy,Iz,Ic)/asiga
						    ENDDO
					    ENDDO
				    ENDDO
			    ENDDO
			ELSE
			    guess=Esol
			ENDIF
		ENDIF
		CALL GIEM2G_FGMRES(ie_op,fgmres_ctl,Esol,guess)


		IF (ie_op%real_space) THEN
			DO Iy=1,Ny_loc
				DO Ix=1,ie_op%Nx
					DO Ic=1,3
						DO Iz=1,Nz
							asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
							Esol(Ix,Iy,Iz,Ic)=asiga*Esol(Ix,Iy,Iz,Ic)
						ENDDO
					ENDDO
				ENDDO
			ENDDO
			    DEALLOCATE(guess)
		ENDIF
		time2=GetTime()
	        CALL LOGGER('Solving has been finished')
		CALL PRINT_CALC_TIME('Total time:                                               ',time2-time1)
		CALL PRINT_CALC_TIME('Time for multiplications:					', ie_op%counter%apply)
		CALL PRINT_CALC_TIME('Time for dot products:					', ie_op%counter%dotprod)
		CALL PRINT_CALC_TIME('Average dot product:					', ie_op%counter%dotprod/ie_op%counter%dotprod_num)
		CALL PRINT_CALC_TIME('Average fftw forward:					', ie_op%counter%mult_fftw/ie_op%counter%mult_num)
		CALL PRINT_CALC_TIME('Average fftw backward:					', ie_op%counter%mult_fftw_b/ie_op%counter%mult_num)
		CALL PRINT_CALC_TIME('Average zgemv:						', ie_op%counter%mult_zgemv/ie_op%counter%mult_num)
		CALL PRINT_CALC_TIME('Average mult:						', ie_op%counter%apply/ie_op%counter%mult_num)
		CALL PRINT_CALC_TIME('Total mult/dp:						', ie_op%counter%apply/ie_op%counter%dotprod)
		CALL PRINT_CALC_TIME('Average mult/dp:						',&
                                & ie_op%counter%apply/ie_op%counter%dotprod*ie_op%counter%dotprod_num/ie_op%counter%mult_num)
		CALL PRINT_CALC_NUMBER('Number of matrix-vector multiplications:                  ', ie_op%counter%mult_num)
		CALL PRINT_CALC_NUMBER('Number of dotproducts:					  ', ie_op%counter%dotprod_num)
                CALL PRINT_BORDER
	END SUBROUTINE

	SUBROUTINE SetAnomalySigma(ie_op,siga,freq)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		REAL(REALPARM),INTENT(IN),POINTER::siga(:,:,:)
		REAL(REALPARM),INTENT(IN)::freq
		REAL(REALPARM)::w
		INTEGER::Ix,Iy,Iz
		IF (ie_op%real_space) THEN
			IF (ASSOCIATED(ie_op%csiga)) DEALLOCATE(ie_op%csiga)
			ALLOCATE(ie_op%csiga(ie_op%Nx,ie_op%Ny_loc,ie_op%Nz))
			w=freq*PI*2
			DO Iz=1,ie_op%Nz
			    DO Iy=1,ie_op%Ny_loc
				DO Ix=1,ie_op%Nx
#ifndef NO_DISPLACEMENT_CURRENTS
					ie_op%csiga(Ix,Iy,Iz)=siga(Iz,Ix,Iy)-C_IONE*w*EPS0	
#else
					ie_op%csiga(Ix,Iy,Iz)=siga(Iz,Ix,Iy)	
#endif
				   ENDDO
			ENDDO
		    ENDDO
		ENDIF

	ENDSUBROUTINE
!-------------------------------------PRIVATE---------------------------------------------------------------------!
	SUBROUTINE GIEM2G_FGMRES(ie_op,fgmres_ctl,Esol, guess)
		TYPE(IntegralEquationOperator),TARGET,INTENT(INOUT)::ie_op
		TYPE (FGMRES_CTL_TYPE),INTENT(IN)::fgmres_ctl
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Esol(:,:,:,:),guess(:,:,:,:)

                TYPE(SEQ_FGMRES_DATA)::seq_solver
                PROCEDURE (MatrixVectorMult),POINTER::Apply=>FULL_APPLY
                PROCEDURE (MANYDOTPRODUCTS),POINTER::MANYDP=>DISTRIBUTED_DOT_PRODUCT
                PROCEDURE (InformAboutIteration),POINTER::InformMe=>Information
		COMPLEX(REALPARM),ALLOCATABLE::solution(:),initial_guess(:),rhs(:)
		COMPLEX(REALPARM),POINTER::pguess(:),pEsol(:)
                TYPE(C_PTR)::ptmp
		INTEGER(8)::Length
		INTEGER(MPI_CTL_KIND),TARGET::comm
                TYPE(C_PTR)::pA,pcomm
		INTEGER(MPI_CTL_KIND)::	IERROR
		REAL(8)::time1,time2
		REAL(8)::full_time
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A, I6, A, ES10.2E2)"
		CHARACTER(LEN=2048,KIND=C_CHAR)::message 
                INTEGER::Nloc
                TYPE(GMRES_PARAMETERS)::params
		INTEGER::MM(2)
                TYPE(RESULT_INFO)::info
		full_time=GetTime()
                CALL DROP_IE_COUNTER(ie_op)
                params%MaxIt=fgmres_ctl%fgmres_maxit
                params%Tol=fgmres_ctl%misfit
                params%RESTART_RESIDUAL=.FALSE.

		MM(1)=fgmres_ctl%fgmres_buf
		MM(2)=fgmres_ctl%gmres_buf
		Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
                CALL INIT_SEQ_FGMRES(seq_solver,Nloc,MM,2,params,Length)
                pA=C_LOC(ie_op)

                CALL SET_OPERATORS(seq_solver,Apply,MANYDP,informMe,pA)
                comm=ie_op%ie_comm
                pcomm=C_LOC(comm)

                CALL SET_DP_INSTANCE(seq_solver,pA)

		WRITE (message,*) 'Solver needs',Length*16.0/1024/1024/1024, 'Gb for workspace'
		
		CALL LOGGER(message)

                ALLOCATE(solution(1:Nloc),initial_guess(1:Nloc),rhs(1:Nloc))
                solution=C_ONE
                initial_guess=C_ONE
                rhs=C_ONE
		IF (ie_op%real_space) THEN
        		Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2

                        ptmp=C_LOC(Esol(1,1,1,1))
                        CALL C_F_POINTER(ptmp,pEsol,(/2*Nloc/))

                        ptmp=C_LOC(guess(1,1,1,1))
                        CALL C_F_POINTER(ptmp,pguess,(/2*Nloc/))

			CALL MPI_SEND(pguess(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX,&
                                & ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)


        		Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2

			CALL MPI_SEND(pEsol(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX,&
                                & ie_op%partner, ie_op%me+ie_op%comm_size,ie_op%ie_comm, IERROR)

                        CALL ZCOPY(Nloc,Esol,ONE,rhs,ONE)
                        CALL ZCOPY(Nloc,pguess,ONE,initial_guess,ONE)

		ELSE
        		Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
			CALL MPI_RECV(initial_guess,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner,&
                                & ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)

	        	Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
			CALL MPI_RECV(rhs,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner,&
                                & ie_op%partner+ie_op%comm_size,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
		ENDIF

                CALL SEQ_FGMRES_SOLVE(seq_solver,solution,initial_guess,rhs,info)
		IF (ie_op%real_space) THEN
	        	Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
			CALL MPI_RECV(pEsol(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX,&
                                & ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
	        	Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
			CALL ZCOPY (Nloc, solution, ONE, pEsol, ONE)
		ELSE
	        	Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
			CALL MPI_SEND(solution,Nloc, MPI_DOUBLE_COMPLEX,&
                                & ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
		ENDIF

                DEALLOCATE(solution,initial_guess,rhs)
                CALL DELETE_SEQ_FGMRES(seq_solver)


		full_time=GetTime()-full_time
		ie_op%counter%solving=full_time

		IF (info%stat==CONVERGED ) THEN
                        WRITE(message,info_fmt) 'GFGMRES converged in', info%Iterations, &
				&' iterations with misfit:',info%be	
		ELSEIF (info%stat==NON_TRUE_CONVERGED) THEN
                        WRITE(message,info_fmt) 'GFGMRES  more or less converged in', info%Iterations, &
				&' iterations with Arnoldy misfit',info%bea, &
                                &' and true misfit', info%be	
                ELSEIF (info%stat==NON_CONVERGED) THEN
                        WRITE(message,info_fmt) 'GFGMRES has NOT converged in', info%Iterations, &
				&' iterations with Arnoldy  misfit:',info%bea	
                ELSE
			WRITE(message,'(A)') 'Unknown behaviour of GFGMRES'
		ENDIF
		CALL LOGGER(message)
	END SUBROUTINE

        SUBROUTINE FULL_APPLY (pA,v_in,v_out)
                        TYPE(C_PTR),INTENT(IN)::pA
                        COMPLEX(REALPARM),INTENT(IN)::v_in(:)
                        COMPLEX(REALPARM),INTENT(INOUT)::v_out(:)
        		TYPE(IntegralEquationOperator),POINTER::ie_op
                        COMPLEX(REALPARM),TARGET::tmp_in(2*SIZE(v_in))
                        COMPLEX(REALPARM),TARGET::tmp_out(2*SIZE(v_out))
                        INTEGER::Nloc
        		INTEGER(MPI_CTL_KIND)::	IERROR
                        CALL C_F_POINTER(pA,ie_op)
		        Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
                        IF (ie_op%real_space) THEN
                                CALL ZCOPY (Nloc, v_in, ONE, tmp_in, ONE)
		                Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
                                CALL MPI_RECV(tmp_in(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX,&
                                        & ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)

                                CALL APPLY_PRECOND(ie_op,tmp_in,tmp_out)
        		        Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2

                                CALL MPI_SEND(tmp_out(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX,&
                                        & ie_op%partner, ie_op%me+ie_op%comm_size,ie_op%ie_comm, IERROR)
                                CALL ZCOPY (Nloc, tmp_out, ONE, v_out, ONE)

                        ELSE
                                CALL MPI_SEND(v_in,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner,&
                                        & ie_op%me,ie_op%ie_comm, IERROR)

                                CALL APPLY_IE_ZEROS(ie_op)
	        	        Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2

                                CALL MPI_RECV(v_out,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner,&
                                        & ie_op%partner+ie_op%comm_size,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
                        ENDIF
        ENDSUBROUTINE




        SUBROUTINE DISTRIBUTED_DOT_PRODUCT(Matrix,v,K,res,ptr)
                COMPLEX(REALPARM),  INTENT(IN):: Matrix(:,:)
                COMPLEX(REALPARM),  INTENT(IN):: v(:)
                INTEGER                   , INTENT(IN):: K
                COMPLEX(REALPARM),  INTENT(INOUT):: res(:)
                TYPE(C_PTR),INTENT(IN)::ptr
       		TYPE(IntegralEquationOperator),POINTER::ie_op
                COMPLEX(REALPARM),TARGET::tmp(K)
		INTEGER(MPI_CTL_KIND)::	IERROR,comm,KK
                REAL(DOUBLEPARM)::time1,time2
                INTEGER::N
		time1=GetTime()
                N=SIZE(V)
                CALL C_F_POINTER(ptr,ie_op)
                CALL ZGEMV('C',N,K,C_ONE,Matrix,N,v,ONE,C_ZERO,res,ONE)

                comm=ie_op%ie_comm
                KK=K
   !             IF( (ie_op%counter%dotprod_num >230).AND. (ie_op%counter%dotprod_num <250))THEN
    !                    PRINT*,'DP $$',ie_op%me, ie_op%counter%dotprod_num,res(1)
    !                    CALL MPI_BARRIER(comm,IERROR)
    !            ENDIF
!                CALL MPI_ALLREDUCE(tmp,res(1:KK),KK,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,IERROR)
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,res(1:KK),KK,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,IERROR)
                IF (IERROR/=0) THEN
                        PRINT*,'DP ##',ie_op%me, ie_op%counter%dotprod_num,res(1),IERROR
                ENDIF
!                IF( (ie_op%counter%dotprod_num >230).AND. (ie_op%counter%dotprod_num <250))THEN
 !                       PRINT*,'DP ##',ie_op%me, ie_op%counter%dotprod_num,res(1)
 !                       CALL MPI_BARRIER(comm,IERROR)
  !              ENDIF
		time2=GetTime()
                ie_op%counter%dotprod_num=ie_op%counter%dotprod_num+K
                ie_op%counter%dotprod=ie_op%counter%dotprod+time2-time1
        ENDSUBROUTINE

		FUNCTION Information(info) RESULT(interp)
			TYPE(RESULT_INFO),INTENT(IN)::info
			INTEGER::interp
                        CHARACTER(LEN=2048,KIND=C_CHAR)::message 
                        WRITE (message,'(A, I5, A, ES10.3E3)') "FGMRES outer  iteration:", info%Iterations,&
                                & " Arnoldy misfit: ", info%bea
                        CALL Logger(message)
                        interp=NOINTERUPT
                ENDFUNCTION

END MODULE  
