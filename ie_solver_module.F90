MODULE IE_SOLVER_MODULE
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE
	USE DATA_TYPES_MODULE
	USE INTEGRAL_EQUATION_MODULE
	USE Timer_Module 
        USE LOGGER_MODULE
	USE APPLY_IE_OPERATOR_MODULE 

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
	        CALL LOGGER('Solver started')
		time1=GetTime()
		Nx=ie_op%Nx
		Nz=ie_op%Nz
		Ny_loc=ie_op%Ny_loc
		ALLOCATE(guess(ie_op%Nx,ie_op%Ny_loc,ie_op%Nz,3))
		IF (ie_op%real_space) THEN
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
	!		    DEALLOCATE(guess)
		ENDIF
                DEALLOCATE(guess)
		time2=GetTime()
	        CALL LOGGER('Solver finished')
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
		TYPE(IntegralEquationOperator),TARGET,INTENT(IN)::ie_op
		TYPE (FGMRES_CTL_TYPE),INTENT(IN)::fgmres_ctl
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Esol(:,:,:,:),guess(:,:,:,:)

                TYPE(SEQ_FGMRES_DATA)::seq_solver
                PROCEDURE (MatrixVectorMult),POINTER::Apply=>FULL_APPLY
                PROCEDURE (MANYDOTPRODUCTS),POINTER::MANYDP=>DISTRIBUTED_DOT_PRODUCT
                PROCEDURE (InformAboutIteration),POINTER::InformMe=>Information
		COMPLEX(REALPARM),POINTER::x(:),x0(:),b(:)
		INTEGER(8)::Length
		INTEGER(MPI_CTL_KIND),TARGET::comm
                TYPE(C_PTR)::pA,pcomm
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
                params%RESTART_RESIDUAL=ITERATIVE_RESIDUAL_AT_RESTART

		MM(1)=fgmres_buf
		MM(2)=gmres_buf

                CALL INIT_SEQ_FGMRES(seq_solver,Nloc,MM,2,params,Length)

                CALL SET_OPERATORS(seq_solver,Apply,MANYDP,informMe,pA)
                comm=ie_op%comm
                pcomm=C_LOC(comm)
                CALL SET_DP_INSTANCE(seq_solver,pcomm)

		WRITE (message,*) 'Solver needs',Length*16.0/1024/1024/1024, 'Gb for workspace'
		
		CALL LOGGER(message)
		Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
		IF (ie_op%real_space) THEN
                        b(1:2*Nloc)=>Esol
                        x0(1:2*Nloc)=>guess
                        x(1:Nloc)=>x0
                ELSE
                        ALLOCATE(b(Nloc))
                        x0(1:Nloc)=>guess
                        x(1:Nloc)=>x0
                ENDIF
		IF (ie_op%real_space) THEN
			CALL MPI_SEND(x0(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
			CALL MPI_SEND(b(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me+ie_op%comm_size,ie_op%ie_comm, IERROR)
		ELSE
			CALL MPI_RECV(x0,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
			CALL MPI_RECV(b,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner+ie_op%comm_size,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
		ENDIF

                CALL SEQ_FGMRES_SOLVE(seq_solver,x,x0(1:Nloc),b(1:Nloc),info)

		IF (ie_op%real_space) THEN
			CALL MPI_RECV(x0(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
			CALL ZCOPY (2*Nloc, x0, ONE, Esol, ONE)
		ELSE
			CALL MPI_SEND(x,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
                        DEALLOCATE(b)
		ENDIF

		full_time=GetTime()-full_time
		ie_op%counter%solving=full_time

		IF (info%stat==CONVERGED ) THEN
                        WRITE(message,info_fmt) 'GFGMRES converged in', info%Iterations, &
				&' iterations with misfit:',info%be	
		ELSEIF (info%stat==NON_TRUE_CONVERGED) THEN
                        WRITE(message,info_fmt) 'GFGMRES  more or less converged in', info%Iterations, &
				&' iterations with Arnoldy misfit',info%bea, &
                                &' and true misfit' info%be	
                ELSEIF (info%stat==NON_CONVERGED) THEN
                        WRITE(message,info_fmt) 'GFGMRES has NOT converged in', info%Iterations, &
				&' iterations with Arnoldy  misfit:',info%bea	
                ELSE
			WRITE(message,'(A)') 'Unknown behaviour of GFGMRES'
		ENDIF
		CALL LOGGER(message)
	END SUBROUTINE

	SUBROUTINE GIEM2G_FGMRES_2(ie_op,fgmres_ctl,Esol, guess)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		TYPE (FGMRES_CTL_TYPE),INTENT(IN)::fgmres_ctl
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Esol(:,:,:,:),guess(:,:,:,:)
		INTEGER :: buf_length,maxit,maxit_precond
		REAL(REALPARM)::misfit
		COMPLEX(REALPARM),POINTER :: work_fgmres(:)
		COMPLEX(REALPARM),POINTER :: work_gmres(:)
		COMPLEX(REALPARM),POINTER ::v_in(:),v_out(:)
		COMPLEX(REALPARM),POINTER ::v_in2(:),v_out2(:)
		COMPLEX(REALPARM),POINTER::tmp_in(:),tmp_out(:);
		COMPLEX(REALPARM),POINTER :: aux(:)
		COMPLEX(REALPARM),POINTER::tmp(:)
		INTEGER :: l_fgmres
		INTEGER :: l_gmres
		INTEGER::m
		integer ::  m2,lwork2
		INTEGER :: revcom, colx, coly, colz, nbscal
		INTEGER :: revcom2, colx2, coly2, colz2, nbscal2
		INTEGER :: irc(7), icntl(7), info(3)
		INTEGER:: irc2(5), icntl2(8), info2(3)
		INTEGER(MPI_CTL_KIND)::	IERROR,comm_inner,comm_outer


		INTEGER::   MATVEC, PRECONDLEFT,PRECONDRIGHT,DOTPROD
		parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)
		INTEGER::Nloc
		INTEGER(8)::N_unknowns
		integer nout
		INTEGER::Nfgmres,Ngmres
		real(8)::  cntl(3), rinfo
		real(8)::  cntl2(5), rinfo2(2)
		real(8):: rn, rx, rc
		REAL(8)::time1,time2
		REAL(8)::full_time
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A, I6, A, ES10.2E2)"
		CHARACTER(LEN=2048,KIND=C_CHAR)::message 

		full_time=GetTime()
                CALL DROP_IE_COUNTER(ie_op)
		Nfgmres=0

		!CALL MPI_COMM_DUP(ie_op%ie_comm, comm_inner, IERROR)
		!CALL MPI_COMM_DUP(ie_op%ie_comm, comm_outer, IERROR)
		comm_inner=ie_op%ie_comm
		comm_outer=ie_op%ie_comm
		CALL init_zfgmres(icntl,cntl)
		CALL init_zgmres(icntl2,cntl2)
		misfit=fgmres_ctl%misfit
		buf_length=fgmres_ctl%fgmres_buf
		maxit=fgmres_ctl%fgmres_maxit
		maxit_precond=fgmres_ctl%gmres_buf

		Nloc=3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/2
		m=buf_length
		l_fgmres = m*m + m*(2*Nloc+6) + 6*Nloc + 1
		l_gmres=maxit_precond*maxit_precond+maxit_precond*(Nloc+6)+6*Nloc+1
		WRITE (message,*) 'Solver needs',(l_fgmres+l_gmres+Nloc+2+4*Nloc)*16.0/1024/1024/1024, 'Gb for workspace'
		
		CALL LOGGER(message)

		ALLOCATE(work_fgmres(l_fgmres),work_gmres(l_gmres))
		ALLOCATE(aux(Nloc+2)) ! possible it is greater than necessary
		ALLOCATE(tmp_in(2*Nloc))
		ALLOCATE(tmp_out(2*Nloc))
		cntl(1) = misfit
		cntl2(1) = misfit

		icntl2(1) = 0
		icntl2(2) = 0
		icntl(6) = maxit 
		IF (ie_op%me==0) THEN
			icntl(1) = 0
			icntl(2)=0
			icntl(3) = 0
			icntl2(3) = 0
		ELSE
			icntl(1) = 0
			icntl(2) = 0
			icntl(3) = 0
			icntl2(3)=0
		ENDIF

		icntl(5)=1
		icntl(4)=fgmres_ctl%ort_type
		icntl2(4) = 0
		icntl2(5) =fgmres_ctl%ort_type
		icntl2(6)=1
		icntl2(7) = maxit_precond
		IF (ie_op%real_space) THEN
			CALL ZCOPY (2*Nloc, guess, ONE, work_fgmres, ONE)
			CALL MPI_SEND(work_fgmres(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
			CALL ZCOPY (2*Nloc, Esol, ONE, work_fgmres(Nloc+1:), ONE)
			CALL MPI_SEND(work_fgmres(2*Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me+ie_op%comm_size,ie_op%ie_comm, IERROR)
		ELSE
			CALL MPI_RECV(work_fgmres,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
			CALL MPI_RECV(work_fgmres(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner+ie_op%comm_size,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
		ENDIF
		N_unknowns=3*ie_op%Nz*ie_op%Nx*ie_op%Ny
		CALL PRINT_CALC_NUMBER('Number of unknowns:', N_unknowns)

		DO
			CALL drive_zfgmres(N_unknowns,Nloc,m,l_fgmres,work_fgmres,irc,icntl,cntl,info,rinfo)
			revcom = irc(1)
			colx   = irc(2)
			coly   = irc(3)
			colz   = irc(4)
			nbscal = irc(5)
			v_in=>work_fgmres(colx:colx+Nloc-1)
			v_out=>work_fgmres(colz:colz+Nloc-1)
			IF ((info(1)==77).OR.(info(1)==0).OR.(info(1)==-4).OR.(info(1)==79)) THEN
				Nfgmres=Nfgmres+1
 				WRITE(message,info_fmt) 'FGMRES interation', Nfgmres,' Arnoldi b.e.:' , rinfo
				CALL LOGGER(message)
			ENDIF
			SELECT CASE(REVCOM)
			CASE (MATVEC)
				IF (ie_op%real_space) THEN
					CALL ZCOPY (Nloc, v_in, ONE, tmp_in, ONE)
					CALL MPI_RECV(tmp_in(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
					CALL APPLY_IE_OP(ie_op,tmp_in,tmp_out)
					CALL MPI_SEND(tmp_out(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me+ie_op%comm_size,ie_op%ie_comm, IERROR)
					CALL ZCOPY (Nloc, tmp_out, ONE, v_out, ONE)

				ELSE
					CALL MPI_SEND(v_in,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
					CALL APPLY_IE_ZEROS(ie_op)
					CALL MPI_RECV(v_out,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner+ie_op%comm_size,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
				ENDIF

			CASE(PRECONDRIGHT)
				work_gmres(1:Nloc)=v_in
				work_gmres(Nloc+1:2*Nloc)=v_in
				Ngmres=0
				DO
					icntl2(1) = 0
					icntl2(2) = 0
					icntl2(3) = 0
					icntl2(4) = 0
					icntl2(5) =fgmres_ctl%ort_type
					icntl2(6)=1
					icntl2(7) = maxit_precond
					CALL drive_zgmres(N_unknowns,Nloc,maxit_precond,l_gmres,&
						& work_gmres(1:),irc2,icntl2,cntl2,info2,rinfo2) 
					revcom2 = irc2(1)
					colx2	= irc2(2) 
					coly2	= irc2(3) 
					colz2	= irc2(4) 
					nbscal2 = irc2(5)
					IF ((info2(1)==77).OR.(info2(1)==0).OR.(info2(1)==-4)) THEN
						Ngmres=Ngmres+1
#ifdef GMRES_VERBOSE
		 				WRITE(message,info_fmt),'GMRES interation', Ngmres,' Arnoldi b.e.:' , rinfo2(1)
						CALL LOGGER(message) 						
#endif
					ENDIF
					SELECT CASE (REVCOM2)
						CASE(MATVEC)
							v_in2=>work_gmres(colx2:colx2+Nloc-1)
							v_out2=>work_gmres(colz2:colz2+Nloc-1)
							IF (ie_op%real_space) THEN
								CALL ZCOPY (Nloc, v_in2, ONE, tmp_in, ONE)
								CALL MPI_RECV(tmp_in(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
								CALL APPLY_IE_OP(ie_op,tmp_in,tmp_out)
								CALL MPI_SEND(tmp_out(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me+ie_op%comm_size,ie_op%ie_comm, IERROR)
								CALL ZCOPY (Nloc, tmp_out, ONE, v_out2, ONE)

							ELSE
								CALL MPI_SEND(v_in2,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
								CALL APPLY_IE_ZEROS(ie_op)
								CALL MPI_RECV(v_out2,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner+ie_op%comm_size,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
							ENDIF
						CASE (DOTPROD)
							time1=GetTime()
							CALL zgemv('C',Nloc,nbscal2,C_ONE,&
								&work_gmres(colx2:),Nloc,&
								&work_gmres(coly2:),1,C_ZERO,aux,1)
							CALL MPI_ALLREDUCE(aux,work_gmres(colz2:),nbscal2,&
								&MPI_DOUBLE_COMPLEX,&
							&MPI_SUM,comm_inner,IERROR)
							time2=GetTime()
							ie_op%counter%dotprod_num=ie_op%counter%dotprod_num+nbscal2
							ie_op%counter%dotprod=ie_op%counter%dotprod+time2-time1
						CASE DEFAULT
							EXIT
					ENDSELECT
				ENDDO
				v_out=work_gmres(1:Nloc)
			CASE (DOTPROD)
				time1=GetTime()
				CALL ZGEMV('C',Nloc,nbscal,C_ONE,work_fgmres(colx:),&
					&Nloc,work_fgmres(coly:),ONE,C_ZERO,aux,ONE)
				CALL MPI_ALLREDUCE(aux,work_fgmres(colz:),nbscal,&
					&MPI_DOUBLE_COMPLEX,MPI_SUM,comm_outer,IERROR)
				time2=GetTime()
				ie_op%counter%dotprod_num=ie_op%counter%dotprod_num+nbscal
				ie_op%counter%dotprod=ie_op%counter%dotprod+time2-time1
			CASE DEFAULT
				EXIT
			END SELECT
		ENDDO
		IF (ie_op%real_space) THEN
			CALL MPI_RECV(work_fgmres(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
			CALL ZCOPY (2*Nloc, work_fgmres, ONE, Esol, ONE)
		ELSE
			CALL MPI_SEND(work_fgmres,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
		ENDIF
		DEALLOCATE(aux,work_fgmres,work_gmres)
		DEALLOCATE(tmp_in,tmp_out)
		!CALL MPI_COMM_FREE(comm_inner,IERROR)
		!CALL MPI_COMM_FREE(comm_outer,IERROR)
		full_time=GetTime()-full_time
		ie_op%counter%solving=full_time
		IF (info(1)==0) THEN
			WRITE(message,info_fmt) 'FGMRES converged in', Nfgmres, &
				&' iterations with b.e:',rinfo	
		ELSEIF (info(1)==-4) THEN
			WRITE(message,info_fmt)'FGMRES not converged in', Nfgmres, &
				&' iterations. Arnoldy b.e. is:',rinfo	
		ELSE
			WRITE(message,'(A)') 'Unknown behaviour of FGMRES'
		ENDIF
		CALL LOGGER(message)
	END SUBROUTINE


        SUBROUTINE FULL_APPLY (pA,v_in,v_out)
                        TYPE(C_PTR),INTENT(IN)::pA
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::v_in(:)
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::v_out(:)
        		TYPE(IntegralEquationOperator),POINTER::ie_op
                        COMPLEX(REALPARM)::tmp_in(SIZE(v_in))
                        COMPLEX(REALPARM)::tmp_out(SIZE(v_out))
                        INTEGER::Nloc
                        CALL C_F_POINTER(pA,ie_op)
                        Nloc=SIZE(v_in)
                        IF (ie_op%real_space) THEN
                                CALL ZCOPY (Nloc, v_in, ONE, tmp_in, ONE)
                                CALL MPI_RECV(tmp_in(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
                                CALL APPLY_IE_OP(ie_op,tmp_in,tmp_out)
                                CALL MPI_SEND(tmp_out(Nloc+1:),Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me+ie_op%comm_size,ie_op%ie_comm, IERROR)
                                CALL ZCOPY (Nloc, tmp_out, ONE, v_out, ONE)

                        ELSE
                                CALL MPI_SEND(v_in,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%me,ie_op%ie_comm, IERROR)
                                CALL APPLY_IE_ZEROS(ie_op)
                                CALL MPI_RECV(v_out,Nloc, MPI_DOUBLE_COMPLEX, ie_op%partner, ie_op%partner+ie_op%comm_size,ie_op%ie_comm, IERROR,MPI_STATUS_IGNORE)
                        ENDIF
        ENDSUBROUTINE




        SUBROUTINE DISTRIBUTED_DOT_PRODUCT(Matrix,v,K,res,ptr)
                COMPLEX(REALPARM), POINTER, INTENT(IN):: Matrix(:,:)
                COMPLEX(REALPARM), POINTER, INTENT(IN):: v(:)
                INTEGER                   , INTENT(IN):: K
                COMPLEX(REALPARM), POINTER, INTENT(IN):: res(:)
                TYPE(C_PTR),INTENT(IN)::ptr
                COMPLEX(REALPARM)::tmp(K)
		INTEGER(MPI_CTL_KIND)::	IERROR,comm
		INTEGER(MPI_CTL_KIND),POINTER::	pcomm
                INTEGER::N
                N=SIZE(V)
                CALL ZGEMV('C',N,K,C_ONE,Matrix,N,v,ONE,C_ZERO,tmp,ONE)
                CALL C_F_POINTER(ptr,pcomm)
                comm=pcomm
                CALL MPI_ALLREDUCE(res,tmp,K,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,IERROR)
        ENDSUBROUTINE

		FUNCTION Information(info) RESULT(interp)
			TYPE(RESULT_INFO),INTENT(IN)::info
			INTEGER::interp
                        PRINT*, "Iteration: " , info%Iterations ,&
                                &"Arnoldy misfit:"  ,info%bea
                ENDFUNCTION

END MODULE  
