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

	SUBROUTINE SolveEquation(int_eq,fgmres_ctl,E_bkg,Esol,E0)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::int_eq
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
		Nx=int_eq%Nx
		Nz=int_eq%Nz
		Ny_loc=int_eq%Ny_loc
		IF (int_eq%real_space) THEN
			DO Iz=1,Nz
				Esol(Iz,:,:,:)=int_eq%sqsigb(Iz)*E_bkg(Iz,:,:,:)
			ENDDO
			ALLOCATE(guess(Nz,3,int_eq%Nx,int_eq%Ny_loc))
			IF (PRESENT(E0)) THEN
				DO Iz=1,Nz
					guess(Iz,:,:,:)=int_eq%sqsigb(Iz)*E0(Iz,:,:,:)
				ENDDO
			ELSE
				guess=Esol
			ENDIF			
		ENDIF
		CALL GIEM2G_FGMRES(int_eq,fgmres_ctl,Esol,guess)


		IF (int_eq%real_space) THEN
			DO Iy=1,Ny_loc
				DO Ix=1,int_eq%Nx
					DO Ic=1,3
						DO Iz=1,Nz
							asiga=C_TWO*int_eq%sqsigb(Iz)/(int_eq%csiga(Iz,Ix,Iy)+CONJG(int_eq%csigb(Iz)))
							Esol(Iz,Ic,Ix,Iy)=asiga*Esol(Iz,Ic,Ix,Iy)
						ENDDO
					ENDDO
				ENDDO
			ENDDO
			    DEALLOCATE(guess)
		ENDIF
		time2=GetTime()
	        CALL LOGGER('Solver finished')
		CALL PRINT_CALC_TIME('Total time:                                               ',time2-time1)
		CALL PRINT_CALC_TIME('Time for multiplications:					', int_eq%counter%apply)
		CALL PRINT_CALC_TIME('Time for dot products:					', int_eq%counter%dotprod)
		CALL PRINT_CALC_TIME('Average dot product:					', int_eq%counter%dotprod/int_eq%counter%dotprod_num)
		CALL PRINT_CALC_TIME('Average fftw forward:					', int_eq%counter%mult_fftw/int_eq%counter%mult_num)
		CALL PRINT_CALC_TIME('Average fftw backward:					', int_eq%counter%mult_fftw_b/int_eq%counter%mult_num)
		CALL PRINT_CALC_TIME('Average zgemv:						', int_eq%counter%mult_zgemv/int_eq%counter%mult_num)
		CALL PRINT_CALC_TIME('Average mult:						', int_eq%counter%apply/int_eq%counter%mult_num)
		CALL PRINT_CALC_TIME('Total mult/dp:						', int_eq%counter%apply/int_eq%counter%dotprod)
		CALL PRINT_CALC_TIME('Average mult/dp:						',&
                                & int_eq%counter%apply/int_eq%counter%dotprod*int_eq%counter%dotprod_num/int_eq%counter%mult_num)
		CALL PRINT_CALC_NUMBER('Number of matrix-vector multiplications:                  ', int_eq%counter%mult_num)
		CALL PRINT_CALC_NUMBER('Number of dotproducts:					  ', int_eq%counter%dotprod_num)
                CALL PRINT_BORDER
	END SUBROUTINE
	SUBROUTINE SetAnomalySigma(int_eq,siga,freq)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::int_eq
		REAL(REALPARM),INTENT(IN),POINTER::siga(:,:,:)
		REAL(REALPARM),INTENT(IN)::freq
		REAL(REALPARM)::w
		IF (int_eq%real_space) THEN
			IF (ASSOCIATED(int_eq%csiga)) DEALLOCATE(int_eq%csiga)
			ALLOCATE(int_eq%csiga(int_eq%Nz,int_eq%Nx,int_eq%Ny_loc))
			w=freq*PI*2
#ifndef NO_DISPLACEMENT_CURRENTS
			int_eq%csiga=siga-C_IONE*w*EPS0	
#else
			int_eq%csiga=siga	
#endif
		ENDIF

	ENDSUBROUTINE
!-------------------------------------PRIVATE---------------------------------------------------------------------!
	SUBROUTINE GIEM2G_FGMRES(int_eq,fgmres_ctl,Esol, guess)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::int_eq
		TYPE (FGMRES_CTL_TYPE),INTENT(IN)::fgmres_ctl
		COMPLEX(REALPARM),POINTER,INTENT(INOUT)::Esol(:,:,:,:),guess(:,:,:,:)
		INTEGER :: buf_length,maxit,maxit_precond
		REAL(REALPARM)::misfit
		COMPLEX(REALPARM),POINTER :: work_fgmres(:)
		COMPLEX(REALPARM),POINTER :: work_gmres(:)
		COMPLEX(REALPARM),POINTER ::v_in(:),v_out(:)
		COMPLEX(REALPARM),POINTER ::v_in2(:),v_out2(:)
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
		LOGICAL :: NextMult
		INTEGER(MPI_CTL_KIND):: IERROR,comm_inner,comm_outer
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
		int_eq%counter%mult_num=0
		int_eq%counter%dotprod_num=0
		int_eq%counter%apply=0d0
		int_eq%counter%mult_fftw=0d0
		int_eq%counter%mult_fftw_b=0d0
		int_eq%counter%mult_zgemv=0d0
		int_eq%counter%dotprod=0d0
		NextMult=.FALSE.
		Nfgmres=0
		IF (.NOT. int_eq%real_space) THEN
			DO
				CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, int_eq%master_proc,int_eq%ie_comm, IERROR)
				IF (NextMult) THEN
					CALL APPLY_IE_OP(int_eq)
				ELSE
					RETURN
				ENDIF
			END DO
		ENDIF
		CALL MPI_COMM_DUP(int_eq%fgmres_comm, comm_inner, IERROR)
		CALL MPI_COMM_DUP(int_eq%fgmres_comm, comm_outer, IERROR)
		
		
		CALL init_zfgmres(icntl,cntl)
		CALL init_zgmres(icntl2,cntl2)
		misfit=fgmres_ctl%misfit
		buf_length=fgmres_ctl%fgmres_buf
		maxit=fgmres_ctl%fgmres_maxit
		maxit_precond=fgmres_ctl%gmres_buf

		Nloc=3*int_eq%Nz*int_eq%Nx*int_eq%Ny_loc
		m=buf_length
		l_fgmres = m*m + m*(2*Nloc+6) + 6*Nloc + 1
		l_gmres=maxit_precond*maxit_precond+maxit_precond*(Nloc+6)+6*Nloc+1
		WRITE (message,*) 'Solver needs',(l_fgmres+l_gmres+Nloc+2)*16.0/1024/1024/1024, 'Gb for workspace'

		CALL LOGGER(message)

		ALLOCATE(work_fgmres(l_fgmres),work_gmres(l_gmres))
		ALLOCATE(aux(Nloc+2)) ! possible it is greater than necessary

		cntl(1) = misfit
		cntl2(1) = misfit

		icntl2(2) = 0
		icntl(6) = maxit 
		IF (int_eq%fgmres_me==0) THEN
			icntl(3) = 220
			icntl2(3) = 320
			icntl(2)=1
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
		CALL ZCOPY (Nloc, guess, ONE, work_fgmres, ONE)
		CALL ZCOPY (Nloc, Esol, ONE, work_fgmres(Nloc+1:), ONE)
		N_unknowns=3*int_eq%Nz*int_eq%Nx*int_eq%Ny
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
				IF (.NOT. int_eq%master) THEN
					CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
					&int_eq%master_proc, int_eq%ie_comm, IERROR)
				ELSE
					NextMult=.TRUE.
					CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
					&int_eq%master_proc, int_eq%ie_comm, IERROR)
				ENDIF
				CALL APPLY_IE_OP(int_eq,v_in,v_out)
			CASE(PRECONDRIGHT)
				work_gmres(1:Nloc)=v_in
				work_gmres(Nloc+1:2*Nloc)=v_in
				Ngmres=0
				DO
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
							IF (.NOT. int_eq%master) THEN
								CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
								&int_eq%master_proc, int_eq%ie_comm, IERROR)
							ELSE
								NextMult=.TRUE.
								CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
								&int_eq%master_proc, int_eq%ie_comm, IERROR)
							ENDIF
							CALL APPLY_IE_OP(int_eq,v_in2,v_out2)
						CASE (DOTPROD)
							time1=GetTime()
							CALL zgemv('C',Nloc,nbscal2,C_ONE,&
								&work_gmres(colx2:),Nloc,&
								&work_gmres(coly2:),1,C_ZERO,aux,1)
							CALL MPI_ALLREDUCE(aux,work_gmres(colz2:),nbscal2,&
								&MPI_DOUBLE_COMPLEX,&
							&MPI_SUM,comm_inner,IERROR)
							time2=GetTime()
							int_eq%counter%dotprod_num=int_eq%counter%dotprod_num+nbscal2
							int_eq%counter%dotprod=int_eq%counter%dotprod+time2-time1
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
				int_eq%counter%dotprod_num=int_eq%counter%dotprod_num+nbscal
				int_eq%counter%dotprod=int_eq%counter%dotprod+time2-time1
			CASE DEFAULT
				EXIT
			END SELECT
		ENDDO
		IF (.NOT. int_eq%master) THEN
			CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
				&int_eq%master_proc, int_eq%ie_comm, IERROR)
		ELSE
			NextMult=.FALSE.
			CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
				&int_eq%master_proc, int_eq%ie_comm, IERROR)
		ENDIF
		CALL ZCOPY (Nloc, work_fgmres, ONE, Esol, ONE)
		DEALLOCATE(aux,work_fgmres,work_gmres)
		CALL MPI_COMM_FREE(comm_inner,IERROR)
		CALL MPI_COMM_FREE(comm_outer,IERROR)
		full_time=GetTime()-full_time
		int_eq%counter%solving=full_time
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

END MODULE  
