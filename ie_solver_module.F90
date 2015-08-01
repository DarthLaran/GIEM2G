MODULE IE_SOLVER_MODULE
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE
	USE DATA_TYPES_MODULE
	USE INTEGRAL_EQUATION_MODULE
	IMPLICIT NONE
	PRIVATE

#define no_compile
!For integral equation kernel!

	PUBLIC:: SolveEquation,SetSigb
	
	CONTAINS

	SUBROUTINE SolveEquation(int_eq,fgmres_ctl,guess)
		TYPE(IntegralEquation),INTENT(INOUT)::int_eq
		TYPE (FGMRES_CTL_TYPE),INTENT(IN)::fgmres_ctl
		COMPLEX(REALPARM),POINTER,INTENT(IN),OPTIONAL::guess(:,:,:,:) !
		COMPLEX(REALPARM),POINTER::tmp(:,:,:,:)
		INTEGER::anom_shape(3)
		INTEGER::tmp_shape(3)
		INTEGER::Iz,Nz,Nx,Ny_loc
		INTEGER(MPI_CTL_KIND)::IERROR
		REAL(8)::time1,time2
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A, ES10.2E3)"
		CALL MPI_BARRIER(int_eq%matrix_comm,IERROR)
		IF (VERBOSE) THEN
			IF (int_eq%master) THEN
				PRINT'(A80)','***********************************************************************************'
				PRINT*, 'Solver started'
			ENDIF
		ENDIF
		time1=MPI_WTIME()
		Nx=int_eq%Nx
		Nz=int_eq%Nz
		Ny_loc=int_eq%Ny_loc
		IF (int_eq%real_space) THEN
			!$OMP PARALLEL	DEFAULT(SHARED) PRIVATE(Iz)
			!$OMP WORKSHARE
				int_eq%sqsigb=SQRT(int_eq%sigb)
				int_eq%asiga=2d0*int_eq%sqsigb/(int_eq%siga+int_eq%sigb)
				int_eq%dsig=(int_eq%siga-int_eq%sigb)
				int_eq%gsig=int_eq%dsig*int_eq%asiga
				int_eq%Esol(:,EX,:,:)=int_eq%sqsigb*int_eq%E_n(:,EX,:,:)
				int_eq%Esol(:,EY,:,:)=int_eq%sqsigb*int_eq%E_n(:,EY,:,:)
				int_eq%Esol(:,EZ,:,:)=int_eq%sqsigb*int_eq%E_n(:,EZ,:,:)
			!$OMP END WORKSHARE
			!$OMP END PARALLEL
			tmp(1:Nz,1:3,1:Nx,1:Ny_loc)=>int_eq%initial_guess
			IF (PRESENT(guess)) THEN
				!$OMP PARALLEL	DEFAULT(SHARED) PRIVATE(Iz)
				!$OMP WORKSHARE
					tmp(:,EX,:,:)=int_eq%sqsigb*guess(:,EX,:,:)
					tmp(:,EY,:,:)=int_eq%sqsigb*guess(:,EY,:,:)
					tmp(:,EZ,:,:)=int_eq%sqsigb*guess(:,EZ,:,:)
				!$OMP END WORKSHARE
				!$OMP END PARALLEL
			ELSE
				!$OMP PARALLEL	DEFAULT(SHARED) PRIVATE(Iz)
				!$OMP WORKSHARE
					tmp(:,EX,:,:)=int_eq%sqsigb*int_eq%E_n(:,EX,:,:)
					tmp(:,EY,:,:)=int_eq%sqsigb*int_eq%E_n(:,EY,:,:)
					tmp(:,EZ,:,:)=int_eq%sqsigb*int_eq%E_n(:,EZ,:,:)
				!$OMP END WORKSHARE
				!$OMP END PARALLEL
			ENDIF			
		ENDIF
		CALL GIEM2G_FGMRES(int_eq,fgmres_ctl)


		IF (int_eq%real_space) THEN
			!$OMP PARALLEL DEFAULT(SHARED)
			!$OMP WORKSHARE
				int_eq%Esol(:,EX,:,:)=int_eq%asiga*int_eq%Esol(:,EX,:,:)
				int_eq%Esol(:,EY,:,:)=int_eq%asiga*int_eq%Esol(:,EY,:,:)
				int_eq%Esol(:,EZ,:,:)=int_eq%asiga*int_eq%Esol(:,EZ,:,:)
			!$OMP END WORKSHARE
			!$OMP END PARALLEL
		ENDIF
		CALL MPI_BARRIER(int_eq%matrix_comm,IERROR)
		time2=MPI_WTIME()
		IF (VERBOSE) THEN
			IF (int_eq%master) THEN
				PRINT*,'Solver finished'
				PRINT info_fmt, 'Total time:							',time2-time1
				PRINT info_fmt, 'Solving time:							',int_eq%counter%solving
				PRINT info_fmt, 'Time for multiplications:				', int_eq%counter%apply
				PRINT info_fmt, 'Time for dot products:					', int_eq%counter%dotprod
				PRINT info_fmt, 'Average dot product:					', int_eq%counter%dotprod/int_eq%counter%dotprod_num
				PRINT info_fmt, 'Average fftw forward:					', int_eq%counter%mult_fftw/int_eq%counter%mult_num
				PRINT info_fmt, 'Average fftw backward:					', int_eq%counter%mult_fftw_b/int_eq%counter%mult_num
				PRINT info_fmt, 'Average zgemv:							', int_eq%counter%mult_zgemv/int_eq%counter%mult_num
				PRINT info_fmt, 'Average mult:							', int_eq%counter%apply/int_eq%counter%mult_num
				PRINT info_fmt, 'Total mult/dp:							', int_eq%counter%apply/int_eq%counter%dotprod
				PRINT info_fmt, 'Average mult/dp:						', int_eq%counter%apply/int_eq%counter%dotprod*int_eq%counter%dotprod_num/int_eq%counter%mult_num
				PRINT*,			'Number of matrix-vector multiplications: ', int_eq%counter%mult_num
				PRINT*,			'Number of dotproducts:					  ', int_eq%counter%dotprod_num
				PRINT'(A80)','***********************************************************************************'
			ENDIF
		ENDIF
	END SUBROUTINE
	SUBROUTINE SetSigb(int_eq,anomaly,bkg)
		TYPE(IntegralEquation),INTENT(INOUT)::int_eq
        TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		INTEGER::anom_shape(3),tmp_shape(3)
		INTEGER::Iz,Nx,Nz,Ny_loc
		Nx=int_eq%Nx
		Nz=int_eq%Nz
		Ny_loc=int_eq%Ny_loc
		IF (int_eq%real_space) THEN
			IF (ASSOCIATED(int_eq%sigb))	THEN
				anom_shape=(/Nz,Nx,Ny_loc/)
				tmp_shape=SHAPE(int_eq%sigb)
				IF ((tmp_shape(1)/=anom_shape(1)).OR.&
					&(tmp_shape(2)/=anom_shape(2)).OR.&
					&(tmp_shape(3)/=anom_shape(3))) THEN
					DEALLOCATE(int_eq%sigb)
					IF(ASSOCIATED(int_eq%sqsigb))DEALLOCATE(int_eq%sqsigb)
					IF(ASSOCIATED(int_eq%asiga))DEALLOCATE(int_eq%asiga)
					IF(ASSOCIATED(int_eq%gsig))DEALLOCATE(int_eq%gsig)
					IF(ASSOCIATED(int_eq%dsig))DEALLOCATE(int_eq%dsig)
					ALLOCATE(int_eq%sigb(Nz,Nx,Ny_loc))
					ALLOCATE(int_eq%sqsigb(Nz,Nx,Ny_loc))
					ALLOCATE(int_eq%asiga(Nz,Nx,Ny_loc))
					ALLOCATE(int_eq%gsig(Nz,Nx,Ny_loc))
					ALLOCATE(int_eq%dsig(Nz,Nx,Ny_loc))
				ENDIF
			ELSE
				IF(ASSOCIATED(int_eq%sqsigb))DEALLOCATE(int_eq%sqsigb)
				IF(ASSOCIATED(int_eq%asiga))DEALLOCATE(int_eq%asiga)
				IF(ASSOCIATED(int_eq%gsig))DEALLOCATE(int_eq%gsig)
				IF(ASSOCIATED(int_eq%dsig))DEALLOCATE(int_eq%dsig)
				ALLOCATE(int_eq%sigb(Nz,Nx,Ny_loc))
				ALLOCATE(int_eq%sqsigb(Nz,Nx,Ny_loc))
				ALLOCATE(int_eq%asiga(Nz,Nx,Ny_loc))
				ALLOCATE(int_eq%gsig(Nz,Nx,Ny_loc))
				ALLOCATE(int_eq%dsig(Nz,Nx,Ny_loc))
			ENDIF
			!$OMP PARALLEL	DEFAULT(SHARED) PRIVATE(Iz)
			!$OMP DO SCHEDULE(GUIDED)
			DO Iz=1,int_eq%Nz
					int_eq%sigb(Iz,:,:)=bkg%sigma(anomaly%Lnumber(Iz))
			ENDDO
			!$OMP ENDDO
			!$OMP END PARALLEL
		ENDIF
	ENDSUBROUTINE

!-------------------------------------PRIVATE---------------------------------------------------------------------!
	SUBROUTINE GIEM2G_FGMRES(int_eq,fgmres_ctl)
		TYPE(IntegralEquation),INTENT(INOUT)::int_eq
		TYPE (FGMRES_CTL_TYPE),INTENT(IN)::fgmres_ctl
		INTEGER :: buf_length,maxit,maxit_precond
		REAL(REALPARM)::misfit
		COMPLEX(REALPARM),POINTER :: work_fgmres(:)
		COMPLEX(REALPARM),POINTER :: work_gmres(:)
		COMPLEX(REALPARM),POINTER ::v_in(:),v_out(:)
		COMPLEX(REALPARM),POINTER ::v_in2(:),v_out2(:)
		COMPLEX(REALPARM),POINTER :: aux(:)
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

		integer nout
		INTEGER::Nfgmres,Ngmres
		real(8)::  cntl(3), rinfo
		real(8)::  cntl2(5), rinfo2(2)
		real(8):: rn, rx, rc
		REAL(8)::time1,time2
		REAL(8)::full_time
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A, I6, A, ES10.2E2)"
		full_time=MPI_Wtime()
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
				CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, int_eq%master_proc,int_eq%matrix_comm, IERROR)
				IF (NextMult) THEN
					CALL APPLY_EQ_OP2(int_eq)
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

		m=buf_length
		l_fgmres = m*m + m*(2*int_eq%Nloc+6) + 6*int_eq%Nloc + 1
		l_gmres=maxit_precond*maxit_precond+maxit_precond*(int_eq%Nloc+6)+6*int_eq%Nloc+1
		IF (int_eq%master) PRINT*,'Solver needs',(l_fgmres+l_gmres+int_eq%Nloc+2)*16.0/1024/1024/1024, 'Gb for workspace'
		ALLOCATE(work_fgmres(l_fgmres),work_gmres(l_gmres))
		ALLOCATE(aux(int_eq%Nloc+2)) ! possible it is greater than necessary

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
		icntl(4)=3
		icntl2(4) = 0
		icntl2(5) =3
		icntl2(6)=1
		icntl2(7) = maxit_precond

		work_fgmres(1:int_eq%Nloc)=int_eq%initial_guess
		work_fgmres(int_eq%Nloc+1:2*int_eq%Nloc)=int_eq%solution
		N64=int_eq%N
		DO
			CALL drive_zfgmres(int_eq%N,int_eq%Nloc,m,l_fgmres,work_fgmres,irc,icntl,cntl,info,rinfo)
			revcom = irc(1)
			colx   = irc(2)
			coly   = irc(3)
			colz   = irc(4)
			nbscal = irc(5)
			v_in=>work_fgmres(colx:colx+int_eq%Nloc-1)
			v_out=>work_fgmres(colz:colz+int_eq%Nloc-1)
			IF ((info(1)==77).OR.(info(1)==0).OR.(info(1)==-4)) THEN
				Nfgmres=Nfgmres+1
				IF ((int_eq%master) .AND.(FGMRES_VERBOSE)) THEN
				       PRINT info_fmt,'FGMRES interation', Nfgmres,' Arnoldi b.e.:' , rinfo
				ENDIF
			ENDIF
			SELECT CASE(REVCOM)
			CASE (MATVEC)
				IF (.NOT. int_eq%master) THEN
					CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
					&int_eq%master_proc, int_eq%matrix_comm, IERROR)
				ELSE
					NextMult=.TRUE.
					CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
					&int_eq%master_proc, int_eq%matrix_comm, IERROR)
				ENDIF
				CALL APPLY_EQ_OP(int_eq,v_in,v_out)
			CASE(PRECONDRIGHT)
				work_gmres(1:int_eq%Nloc)=v_in
				work_gmres(int_eq%Nloc+1:2*int_eq%Nloc)=v_in
				Ngmres=0
				DO
					CALL drive_zgmres(int_eq%N,int_eq%Nloc,maxit_precond,l_gmres,&
						& work_gmres(1:),irc2,icntl2,cntl2,info2,rinfo2) 
					revcom2 = irc2(1)
					colx2	= irc2(2) 
					coly2	= irc2(3) 
					colz2	= irc2(4) 
					nbscal2 = irc2(5)
					IF ((info2(1)==77).OR.(info2(1)==0).OR.(info2(1)==-4)) THEN
						Ngmres=Ngmres+1
						IF ((int_eq%master) .AND.(GMRES_VERBOSE)) THEN
						       PRINT info_fmt,'GMRES interation', Ngmres,' Arnoldi b.e.:' , rinfo2(1)
						ENDIF
					ENDIF
					SELECT CASE (REVCOM2)
						CASE(MATVEC)
							v_in2=>work_gmres(colx2:colx2+int_eq%Nloc-1)
							v_out2=>work_gmres(colz2:colz2+int_eq%Nloc-1)
							IF (.NOT. int_eq%master) THEN
								CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
								&int_eq%master_proc, int_eq%matrix_comm, IERROR)
							ELSE
								NextMult=.TRUE.
								CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
								&int_eq%master_proc, int_eq%matrix_comm, IERROR)
							ENDIF
							CALL APPLY_EQ_OP(int_eq,v_in2,v_out2)
						CASE (DOTPROD)
							time1=MPI_WTIME()
!							PRINT*,'dp1',nbscal2,int_eq%Nloc,colx2,coly2
							CALL zgemv('C',int_eq%Nloc,nbscal2,C_ONE,&
								&work_gmres(colx2:),int_eq%Nloc,&
								&work_gmres(coly2:),1,C_ZERO,aux,1)
!							PRINT*,'dp3'
							CALL MPI_ALLREDUCE(aux,work_gmres(colz2:),nbscal2,&
								&MPI_DOUBLE_COMPLEX,&
							&MPI_SUM,comm_inner,IERROR)
							time2=MPI_WTIME()
							int_eq%counter%dotprod_num=int_eq%counter%dotprod_num+nbscal2
							int_eq%counter%dotprod=int_eq%counter%dotprod+time2-time1
!							PRINT*,'dp2'
						CASE DEFAULT
							EXIT
					ENDSELECT
				ENDDO
				v_out=work_gmres(1:int_eq%Nloc)
			CASE (DOTPROD)
				time1=MPI_WTIME()
				CALL ZGEMV('C',int_eq%Nloc,nbscal,C_ONE,work_fgmres(colx:),&
					&int_eq%Nloc,work_fgmres(coly:),ONE,C_ZERO,aux,ONE)
				CALL MPI_ALLREDUCE(aux,work_fgmres(colz:),nbscal,&
					&MPI_DOUBLE_COMPLEX,MPI_SUM,comm_outer,IERROR)
				time2=MPI_WTIME()
				int_eq%counter%dotprod_num=int_eq%counter%dotprod_num+nbscal
				int_eq%counter%dotprod=int_eq%counter%dotprod+time2-time1
			CASE DEFAULT
				EXIT
			END SELECT
		ENDDO
		IF (.NOT. int_eq%master) THEN
			CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
				&int_eq%master_proc, int_eq%matrix_comm, IERROR)
		ELSE
			NextMult=.FALSE.
			CALL MPI_BCAST(NextMult,1,MPI_LOGICAL, &
				&int_eq%master_proc, int_eq%matrix_comm, IERROR)
		ENDIF
		int_eq%solution=work_fgmres(1:int_eq%Nloc)
		DEALLOCATE(aux,work_fgmres,work_gmres)
		CALL MPI_COMM_FREE(comm_inner,IERROR)
		CALL MPI_COMM_FREE(comm_outer,IERROR)
		full_time=MPI_WTIME()-full_time
		int_eq%counter%solving=full_time
		IF (info(1)==0) THEN
			IF(int_eq%master) PRINT info_fmt,'FGMRES converged in', Nfgmres, &
				&' iterations with b.e:',rinfo	
		ELSEIF (info(1)==-4) THEN
			IF(int_eq%master) PRINT info_fmt,'FGMRES not converged in', Nfgmres, &
				&' iterations. Arnoldy b.e. is:',rinfo	
		ELSE
			IF(int_eq%master) PRINT* ,'Unknown behaviour of FGMRES'
		ENDIF
	END SUBROUTINE

	SUBROUTINE Apply_EQ_OP(int_eq,v_in,v_out)
		TYPE(IntegralEquation),INTENT(INOUT)::int_eq
		COMPLEX(REALPARM),POINTER,INTENT(IN)::v_in(:)
		COMPLEX(REALPARM),POINTER,INTENT(OUT)::v_out(:)
		COMPLEX(REALPARM),POINTER::field_in(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_out(:,:,:,:)
		COMPLEX(REALPARM)::d1,d2
		REAL(8)::time1,time2
		INTEGER::IERROR
		INTEGER ::Ix,Iy,Iz,Ic,Ixy
		time1=MPI_Wtime()
		field_in(1:int_eq%Nz,1:3,1:int_eq%Nx,1:int_eq%Ny_loc)=>v_in
		field_out(1:int_eq%Nz,1:3,1:int_eq%Nx,1:int_eq%Ny_loc)=>v_out
		DO Iy=1,int_eq%Ny_loc
		!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(Ix,Ic,Iz)
			!$OMP DO SCHEDULE(GUIDED)
			DO Ix=1,int_eq%Nx
				DO Ic=1,3
					DO Iz=1,int_eq%Nz
						int_eq%field_in4(Iz,Ic,Ix,Iy)=&
						&field_in(Iz,Ic,Ix,Iy)*int_eq%gsig(Iz,Ix,Iy)
					ENDDO
				ENDDO
			ENDDO
			!$OMP END DO
			!$OMP DO SCHEDULE(GUIDED)
			DO Ix=int_eq%Nx+1,int_eq%Nx2
				DO Ic=1,3
					DO Iz=1,int_eq%Nz
						int_eq%field_in4(Iz,Ic,Ix,Iy)=C_ZERO
					ENDDO
				ENDDO
			ENDDO
			!$OMP END DO
		!$OMP END PARALLEL
		ENDDO
		CALL	Apply_IE_OP(int_eq)
		DO Iy=1,int_eq%Ny_loc
			!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ix,Ic,Iz,d1,d2)
			!$OMP DO SCHEDULE(GUIDED)
			DO Ix=1,int_eq%Nx
				DO Ic=1,3
					DO Iz=1,int_eq%Nz
						d1=field_in(Iz,Ic,Ix,Iy)*int_eq%asiga(Iz,Ix,Iy)
						d2=d1-int_eq%field_out4(Iz,Ic,Ix,Iy)/(int_eq%dz(Iz))
						field_out(Iz,Ic,Ix,Iy)=d2*int_eq%sqsigb(Iz,Ix,Iy)
					ENDDO
				ENDDO
			ENDDO
			!$OMP END DO
			!$OMP END PARALLEL
		ENDDO
!		PRINT*,int_eq%me,'out',v_out(1:5)
		time2=MPI_WTIME()
		int_eq%counter%mult_num=int_eq%counter%mult_num+1
		int_eq%counter%apply=int_eq%counter%apply+time2-time1
	END SUBROUTINE
	SUBROUTINE APPLY_EQ_OP2(int_eq)
		TYPE(IntegralEquation),INTENT(INOUT)::int_eq

		!$OMP PARALLEL	 DEFAULT(SHARED)
		!$OMP WORKSHARE
		int_eq%field_in4=C_ZERO
		!$OMP END WORKSHARE
		!$OMP ENDPARALLEL
		CALL APPLY_IE_OP(int_eq)
	ENDSUBROUTINE
	SUBROUTINE APPLY_IE_OP(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		INTEGER::Ixy,Nz
		REAL(8)::time1,time2,time3,time4
		time1=MPI_WTIME()
		CALL IE_OP_FFTW_FWD(ie_op)
		time2=MPI_WTIME()
		ie_op%counter%mult_fftw=ie_op%counter%mult_fftw+time2-time1
		!$OMP PARALLEL	PRIVATE(Ixy) DEFAULT(SHARED)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ixy=1,ie_op%Nx2Ny_loc
			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EXX,EX,EX,Ixy,C_ONE,C_ZERO)
			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EXY,EY,EX,Ixy,C_ONE,C_ONE)
			CALL IE_OP_ZGEMV_ASYM(ie_op,'N',A_EXZ,EZ,EX,Ixy,C_ONE,C_ONE)

			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EYX,EX,EY,Ixy,C_ONE,C_ZERO)
			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EYY,EY,EY,Ixy,C_ONE,C_ONE)
			CALL IE_OP_ZGEMV_ASYM(ie_op,'N',A_EYZ,EZ,EY,Ixy,C_ONE,C_ONE)

			CALL IE_OP_ZGEMV_ASYM(ie_op,'T',A_EZX,EX,EZ,Ixy,-C_ONE,C_ZERO)
			CALL IE_OP_ZGEMV_ASYM(ie_op,'T',A_EZY,EY,EZ,Ixy,-C_ONE,C_ONE)
			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EZZ,EZ,EZ,Ixy,C_ONE,C_ONE)
		ENDDO
		!$OMP END DO
		!$OMP END  PARALLEL
		time3=MPI_WTIME()
		CALL IE_OP_FFTW_BWD(ie_op)
		time4=MPI_WTIME()
!		PRINT*,ie_op%me,'Tj',ie_op%field_out4(:,EZ,:,:)
		ie_op%counter%mult_fftw_b=ie_op%counter%mult_fftw_b+time4-time3
		ie_op%counter%mult_zgemv=ie_op%counter%mult_zgemv+time3-time2
	ENDSUBROUTINE
	SUBROUTINE IE_OP_ZGEMV_SYMM(ie_op,Tc,c_in,c_out,I,ALPHA,BETA)
			TYPE(IntegralEquation),INTENT(INOUT)::ie_op
			INTEGER,INTENT(IN)::Tc,c_in,c_out,I
			COMPLEX(REALPARM),INTENT(IN)::ALPHA,BETA
			CALL ZSPMV('U',ie_op%Nz,ALPHA,ie_op%G_symm4(:,Tc,I),&
			&ie_op%field_in3(:,c_in,I),ONE,BETA,ie_op%field_out3(:,c_out,I),ONE)
	END SUBROUTINE
	SUBROUTINE IE_OP_ZGEMV_ASYM(ie_op,TRANS,Tc,c_in,c_out,I,ALPHA,BETA)
			TYPE(IntegralEquation),INTENT(INOUT)::ie_op
			CHARACTER,INTENT(IN)::TRANS(*)
			INTEGER,INTENT(IN)::Tc,c_in,c_out,I
			COMPLEX(REALPARM),INTENT(IN)::ALPHA,BETA
			CALL ZGEMV(TRANS,ie_op%Nz,ie_op%Nz,ALPHA,ie_op%G_asym4(:,:,Tc,I),ie_op%Nz,&
			&ie_op%field_in3(:,c_in,I),ONE,BETA,ie_op%field_out3(:,c_out,I),ONE)
	END SUBROUTINE
END MODULE  
