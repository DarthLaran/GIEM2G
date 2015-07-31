MODULE CONTINUATION_FUNCTION_MODULE
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE DATA_TYPES_MODULE
	IMPLICIT NONE
	INTEGER,PARAMETER::COUNTER_ALL=1
	INTEGER,PARAMETER::COUNTER_WT=2
	INTEGER,PARAMETER::COUNTER_LIN=3
	INTEGER,PARAMETER::COUNTER_SQR=4
	INTEGER,PARAMETER::COUNTER_OTHER=5
	TYPE TypeCounter
		REAL(8)::plans
		REAL(8)::tensor_calc(5)
		REAL(8)::tensor_fft
		REAL(8)::mult_fftw
		REAL(8)::mult_zgemv
		REAL(8)::apply
		REAL(8)::dotprod
		REAL(8)::solving
		INTEGER::mult_num
		INTEGER::dotprod_num
	ENDTYPE
	TYPE RC_OPERATOR
		INTEGER::Nx,Ny,Nz,Nr
		INTEGER::Nx2,Ny_loc,Ny_offset
		INTEGER::Nx2Ny2
		INTEGER::Nx2Ny_loc
		COMPLEX(REALPARM),POINTER:: G_E(:,:,:,:,:)!G(Nz,Nr,REXX:REZZ,2*Nx,Ny_loc)
		COMPLEX(REALPARM),POINTER:: G_E_fftw(:,:,:,:,:)! for fft(G)
		COMPLEX(REALPARM),POINTER:: G_E4(:,:,:,:) !for distribudted calc G(Nz,Nr,REXX:REZZ,2*Nx*Ny_loc)

		COMPLEX(REALPARM),POINTER:: G_H(:,:,:,:,:)!G(Nz,Nr,RHXX:RHYZ,2*Nx,Ny_loc)
		COMPLEX(REALPARM),POINTER:: G_H_fftw(:,:,:,:,:)! for fft(G)
		COMPLEX(REALPARM),POINTER:: G_H4(:,:,:,:) !for distribudted calc G(Nz,Nr,RHXX:RHYZ,2*Nx*Ny_loc)


		TYPE(C_PTR)::pG_E
		TYPE(C_PTR)::pG_H
		INTEGER:: matrix_comm
		INTEGER::me
		INTEGER::master_proc
		LOGICAL::master
		TYPE(TypeCounter)::counter
		TYPE (RECEIVER_TYPE),POINTER::recvs(:)
		INTEGER (C_INTPTR_T) ::localsize_in
		INTEGER (C_INTPTR_T) ::localsize_out
		TYPE(C_PTR)::p_in,p_out
		COMPLEX(REALPARM),POINTER::field_in4(:,:,:,:),field_out4(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_in3(:,:,:),field_out3(:,:,:)
		TYPE(C_PTR)::planFwd,planBwd
		REAL(REALPARM),POINTER ::siga(:,:,:)
		REAL(REALPARM),POINTER ::sigb(:,:,:)
		REAL(REALPARM),POINTER ::dsig(:,:,:)
		LOGICAL::real_space
		LOGICAL::fftw_threads_ok
	  ENDTYPE
CONTAINS
	SUBROUTINE PrepareContinuationOperator(rc_op,anomaly,recvs,mcomm,fftw_threads_ok)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		TYPE (RECEIVER_TYPE),POINTER,INTENT(IN)::recvs(:)
		INTEGER,INTENT(IN)::mcomm 
		LOGICAL,OPTIONAL,INTENT(IN)::fftw_threads_ok
		INTEGER::Nx,Ny,Nz
		INTEGER::comm_size,IERROR,me
		Nx=anomaly%Nx
		Ny=anomaly%Ny
		Nz=anomaly%Nz
		rc_op%Nx=Nx
		rc_op%Ny=Ny
		rc_op%Nz=Nz
		rc_op%Nr=SIZE(recvs)
		rc_op%recvs=>recvs
		rc_op%matrix_comm=mcomm
		CALL MPI_COMM_RANK(mcomm, me, IERROR)
		CALL MPI_COMM_SIZE(mcomm,comm_size,IERROR) 
		rc_op%me=me
		if (me==0) THEN
			rc_op%master=.TRUE.
		ELSE
			rc_op%master=.FALSE.
		ENDIF
		IF (mod(2*Ny,comm_size)/=0) THEN
			IF (rc_op%master) THEN
				PRINT*,'Wrong number of processes. GIEM2G forced to halt!'
			ENDIF
			CALL MPI_FINALIZE(IERROR)
			STOP
		ENDIF
		IF (PRESENT(fftw_threads_ok)) THEN
			rc_op%fftw_threads_ok=fftw_threads_ok
		ELSE
			rc_op%fftw_threads_ok=.FALSE.
		ENDIF
		CALL CalcSizesForRC_OP(rc_op)
		CALL AllocateRC_OP(rc_op)
		CALL CalcFFTWPlansRC_OP(rc_op)
		IF (rc_op%Ny_offset >= Ny) THEN	  
			rc_op%real_space=.FALSE.
		ELSE
			rc_op%real_space=.TRUE.
		ENDIF
                rc_op%master_proc=0
		IF ((rc_op%master).AND.(rc_op%Ny_offset/=0)) THEN
			PRINT*,  'Main process has no zero offset.&
				   & Now it is problem. It will be solved... sometime'
				CALL MPI_FINALIZE(IERROR)
				STOP
		ENDIF
		rc_op%siga=>NULL()
		rc_op%sigb=>NULL()
		rc_op%dsig=>NULL()
		anomaly%Ny_loc=rc_op%Ny_loc
	ENDSUBROUTINE
#define no_compile	
	SUBROUTINE ReCalculation(rc_op,Eint,Ea,Ha)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Eint(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(OUT)::Ea(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(OUT)::Ha(:,:,:,:)
		INTEGER::IERROR,Nx
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A, ES10.2E3)"
		REAL(8)::time1,time2
		CALL MPI_BARRIER(rc_op%matrix_comm,IERROR)
		IF (VERBOSE) THEN
			IF (rc_op%master) THEN
				PRINT'(A80)','***********************************************************************************'
				PRINT*, 'Recalculation started'
			ENDIF
		ENDIF
		time1=MPI_WTIME()
		Nx=rc_op%Nx
		IF (rc_op%real_space) THEN
			!$OMP PARALLEL	DEFAULT(SHARED)
			!$OMP WORKSHARE
					rc_op%dsig=(rc_op%siga-rc_op%sigb)
					rc_op%field_in4(:,EX,1:Nx,:)=rc_op%dsig*Eint(:,EX,:,:)
					rc_op%field_in4(:,EY,1:Nx,:)=rc_op%dsig*Eint(:,EY,:,:)
					rc_op%field_in4(:,EZ,1:Nx,:)=rc_op%dsig*Eint(:,EZ,:,:)
					rc_op%field_in4(:,:,Nx+1:,:)=C_ZERO
			!$OMP END WORKSHARE
			!$OMP END PARALLEL
		ELSE
					rc_op%field_in4=C_ZERO
		ENDIF
		CALL APPLY_RC_E_OP(rc_op)
		IF (rc_op%real_space) THEN

			!$OMP PARALLEL	DEFAULT(SHARED)
			!$OMP WORKSHARE
					Ea=rc_op%field_out4(:,:,1:Nx,:)
					rc_op%field_in4(:,EX,1:Nx,:)=rc_op%dsig*Eint(:,EX,:,:)
					rc_op%field_in4(:,EY,1:Nx,:)=rc_op%dsig*Eint(:,EY,:,:)
					rc_op%field_in4(:,EZ,1:Nx,:)=rc_op%dsig*Eint(:,EZ,:,:)
					rc_op%field_in4(:,:,Nx+1:,:)=C_ZERO
				!$OMP END WORKSHARE
			!$OMP END PARALLEL
		ELSE
					rc_op%field_in4=C_ZERO
		ENDIF

		CALL APPLY_RC_H_OP(rc_op)
		IF (rc_op%real_space) THEN
			!$OMP PARALLEL	DEFAULT(SHARED)
			!$OMP WORKSHARE
					Ha=rc_op%field_out4(:,:,1:Nx,:)
				!$OMP END WORKSHARE
			!$OMP END PARALLEL
		ELSE
					Ea=>NULL()
					Ha=>NULL()
		ENDIF
		time2=MPI_WTIME()
		IF (VERBOSE) THEN
			IF (rc_op%master) THEN
				PRINT*,'Recalculation finished'
				PRINT info_fmt, 'Total time:							',time2-time1
				PRINT'(A80)','***********************************************************************************'
			ENDIF
		ENDIF
				
	END SUBROUTINE

	SUBROUTINE CalcSizesForRC_OP(rc_op) 
		TYPE(RC_Operator),INTENT(INOUT)::rc_op
		INTEGER(C_INTPTR_T)::tsize8(2)
		INTEGER(C_INTPTR_T)::Nz3,Nr3
		INTEGER(C_INTPTR_T)::block
		INTEGER(C_INTPTR_T):: CNy,CNy_offset !size and offset for electrical current in Y direction at this process 
		INTEGER(FFTW_COMM_SIZE)::COMM
		REAL(REALPARM)::size_g
		REAL(REALPARM)::tensor_size
		tsize8=(/rc_op%Ny*2,rc_op%Nx*2/)
		
		COMM=rc_op%matrix_comm
		block=FFTW_MPI_DEFAULT_BLOCK!(0)
!----------------------------------------------------------------------------------------------------------!
		Nz3=rc_op%Nz*3
		rc_op%localsize_in = fftw_mpi_local_size_many(FFTW_TWO,tsize8,Nz3,block,COMM, &
			& CNy, CNy_offset)
		rc_op%Ny_loc=INT(CNy,KIND(rc_op%Ny_loc))
		rc_op%Ny_offset=INT(CNy_offset,KIND(rc_op%Ny_offset))
		rc_op%Nx2=rc_op%Nx*2
		rc_op%Nx2Ny_loc=rc_op%Nx2*rc_op%Ny_loc
		rc_op%Nx2Ny2=rc_op%Nx*rc_op%Ny*4

		Nr3=rc_op%Nr*3
		rc_op%localsize_out = fftw_mpi_local_size_many(FFTW_TWO,tsize8,Nr3,block,COMM, &
			& CNy, CNy_offset)


		size_g=rc_op%Nx2Ny_loc*rc_op%Nz*rc_op%Nr*15.0*REALPARM/1024/1024/1024
		IF (VERBOSE) THEN
			IF (rc_op%master) THEN
				PRINT*,'RC matricies needs', size_g, 'Gb per process'
			 ENDIF
		 ENDIF
	ENDSUBROUTINE
	SUBROUTINE AllocateRCMatrix(matrix)
		TYPE(RC_OPERATOR),INTENT(INOUT)::matrix
		INTEGER(C_INTPTR_T)::length
		INTEGER::shape1(1),Nz,Nx2,Ny_loc,N1,N2,Nr
		COMPLEX(REALPARM),POINTER::tmp(:)
		length=matrix%Nx2Ny_loc*matrix%Nz*matrix%Nr*8
		matrix%pG_E=fftw_alloc_complex(length)
		shape1=(/length/)
		CALL c_f_pointer(matrix%pG_E,tmp , shape1)
		Nz=matrix%Nz
		Nr=matrix%Nr
		Nx2=2*matrix%Nx-1
		N1=matrix%Ny_offset
		N2=matrix%Ny_offset+matrix%Ny_loc-1
		
		matrix%G_E(1:Nz,1:Nr,REXX:REZZ,0:Nx2,N1:N2)=>tmp
		matrix%G_E_fftw(1:Nz,1:8,1:Nr,0:Nx2,N1:N2)=>tmp
		matrix%G_E4(1:Nz,1:Nr,REXX:REZZ,1:matrix%Nx2Ny_loc)=>tmp

		matrix%pG_H=fftw_alloc_complex(length)
		CALL c_f_pointer(matrix%pG_H,tmp , shape1)
		
		matrix%G_H(1:Nz,1:Nr,RHXX:RHZY,0:Nx2,N1:N2)=>tmp
		matrix%G_H_fftw(1:Nz,1:7,1:Nr,0:Nx2,N1:N2)=>tmp
		matrix%G_H4(1:Nz,1:Nr,RHXX:RHZY,1:matrix%Nx2Ny_loc)=>tmp
	END SUBROUTINE

	SUBROUTINE AllocateRC_OP(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		INTEGER::shape4(4),shape3(3)
		CALL AllocateRCMatrix(rc_op)
		rc_op%p_in=fftw_alloc_complex(rc_op%localsize_in)
		rc_op%p_out=fftw_alloc_complex(rc_op%localsize_out)
		shape4=(/rc_op%Nz,3,rc_op%Nx2,rc_op%Ny_loc/)
		shape3=(/rc_op%Nz,3,rc_op%Nx2Ny_loc/)
		CALL c_f_pointer(rc_op%p_in,rc_op%field_in4, shape4)
		CALL c_f_pointer(rc_op%p_in,rc_op%field_in3, shape3)
		shape4=(/rc_op%Nr,3,rc_op%Nx2,rc_op%Ny_loc/)
		shape3=(/rc_op%Nr,3,rc_op%Nx2Ny_loc/)
		CALL c_f_pointer(rc_op%p_out,rc_op%field_out4, shape4)
		CALL c_f_pointer(rc_op%p_out,rc_op%field_out3, shape3)
	ENDSUBROUTINE
	SUBROUTINE CalcFFTWPlansRC_OP(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		INTEGER(C_INTPTR_T)::fftwsize(2)
		INTEGER(C_INTPTR_T)::Nz3,Nr3
		INTEGER(C_INTPTR_T)::block
		INTEGER::IERROR
		INTEGER::omp_get_max_threads,nt
		INTEGER(FFTW_COMM_SIZE)::COMM, FFTW_NT
		REAL(8)::time1,time2
		fftwsize=(/rc_op%Ny*2,rc_op%Nx*2/)
		block=FFTW_MPI_DEFAULT_BLOCK
		Nz3=3*rc_op%Nz
		Nr3=3*rc_op%Nr
		CALL MPI_BARRIER(rc_op%matrix_comm,IERROR)
		time1=MPI_WTIME()
		IF (rc_op%fftw_threads_ok) THEN
			NT=OMP_GET_MAX_THREADS()
			FFTW_NT=NT
			CALL FFTW_PLAN_WITH_NTHREADS(FFTW_NT)
		ENDIF
		COMM=rc_op%matrix_comm
!		 CALL fftw_set_timelimit(5d0)
				 
		rc_op%planFWD=fftw_mpi_plan_many_dft(FFTW_TWO,fftwsize,Nz3,block,block,&
		&rc_op%field_in4,rc_op%field_in4,COMM, FFTW_FORWARD ,FFTW_MEASURE)!FFTW_PATIENT);

		rc_op%planBWD=fftw_mpi_plan_many_dft(FFTW_TWO,fftwsize,Nr3,block,block,&
		&rc_op%field_out4,rc_op%field_out4,COMM, FFTW_BACKWARD ,FFTW_MEASURE)!FFTW_PATIENT);
		time2=MPI_WTIME()
		rc_op%counter%plans=time2-time1
		IF (VERBOSE) THEN
			IF (rc_op%master) THEN
				PRINT*,'FFTW3 plan calculations:', time2-time1,'s'
			ENDIF
		ENDIF
	ENDSUBROUTINE
	SUBROUTINE RC_OP_FFTW_FWD(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		CALL fftw_mpi_execute_dft(rc_op%planFWD,rc_op%field_in4,rc_op%field_in4)
	ENDSUBROUTINE
	SUBROUTINE RC_OP_FFTW_BWD(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		CALL fftw_mpi_execute_dft(rc_op%planBWD,rc_op%field_out4,rc_op%field_out4)
		rc_op%field_out4=rc_op%field_out4/rc_op%Nx2Ny2
	ENDSUBROUTINE
	SUBROUTINE CalcFFTofRCTensor(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		INTEGER::Ir
		INTEGER::IERROR
		REAL(8)::time1,time2
		CALL MPI_BARRIER(rc_op%matrix_comm,IERROR)
		time1=MPI_WTIME()
		DO Ir=1,rc_op%Nr
			rc_op%field_in4=rc_op%G_E_fftw(:,1:3,Ir,:,:)
			CALL rc_op_FFTW_FWD(rc_op)
			rc_op%G_E_fftw(:,1:3,Ir,:,:)=rc_op%field_in4

			rc_op%field_in4=rc_op%G_E_fftw(:,4:6,Ir,:,:)
			CALL rc_op_FFTW_FWD(rc_op)
			rc_op%G_E_fftw(:,4:6,Ir,:,:)=rc_op%field_in4

			rc_op%field_in4(:,1:2,:,:)=rc_op%G_E_fftw(:,7:8,Ir,:,:)
			rc_op%field_in4(:,3,:,:)=rc_op%G_H_fftw(:,1,Ir,:,:)
			CALL rc_op_FFTW_FWD(rc_op)
			rc_op%G_E_fftw(:,7:8,Ir,:,:)=rc_op%field_in4(:,1:2,:,:)
			rc_op%G_H_fftw(:,1,Ir,:,:)=rc_op%field_in4(:,3,:,:)

			rc_op%field_in4=rc_op%G_H_fftw(:,2:4,Ir,:,:)
			CALL rc_op_FFTW_FWD(rc_op)
			rc_op%G_H_fftw(:,2:4,Ir,:,:)=rc_op%field_in4

			rc_op%field_in4=rc_op%G_H_fftw(:,5:7,Ir,:,:)
			CALL rc_op_FFTW_FWD(rc_op)
			rc_op%G_H_fftw(:,5:7,Ir,:,:)=rc_op%field_in4
		ENDDO
		time2=MPI_WTIME()
		IF (VERBOSE) THEN
			IF (rc_op%master) THEN
				PRINT*, 'FFT of tensor has been computed: ',time2-time1 ,' s'
			ENDIF
		ENDIF
		rc_op%counter%tensor_fft=time2-time1
	ENDSUBROUTINE

	SUBROUTINE DeleteMatrix(matrix)
		TYPE(RC_OPERATOR),INTENT(INOUT)::matrix
		CALL fftw_free(matrix%pG_E)
		CALL fftw_free(matrix%pG_H)
		matrix%G_E=>NULL()
		matrix%G_E_fftw=>NULL()
		matrix%G_E4=>NULL()
		matrix%G_H=>NULL()
		matrix%G_H_fftw=>NULL()
		matrix%G_H4=>NULL()
	ENDSUBROUTINE
	SUBROUTINE DeleteRC_OP(rc_op)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		CALL fftw_free(rc_op%p_in)
		CALL fftw_free(rc_op%p_out)
		rc_op%field_in4=>NULL()
		rc_op%field_in3=>NULL()
		rc_op%field_out4=>NULL()
		rc_op%field_out3=>NULL()
		CALL fftw_destroy_plan(rc_op%planFwd)
		CALL fftw_destroy_plan(rc_op%planBwd)
		CALL DeleteMatrix(rc_op)
	ENDSUBROUTINE

	SUBROUTINE SetSigbRC(rc_op,anomaly,bkg)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		INTEGER::anom_shape(3),tmp_shape(3)
		INTEGER::Iz,Nx,Nz,Ny_loc
		Nx=rc_op%Nx
		Nz=rc_op%Nz
		Ny_loc=rc_op%Ny_loc
		IF (rc_op%real_space) THEN
			IF (ASSOCIATED(rc_op%sigb))	THEN
				anom_shape=(/Nz,Nx,Ny_loc/)
				tmp_shape=SHAPE(rc_op%sigb)
				IF ((tmp_shape(1)/=anom_shape(1)).OR.&
					&(tmp_shape(2)/=anom_shape(2)).OR.&
					&(tmp_shape(3)/=anom_shape(3))) THEN
					DEALLOCATE(rc_op%sigb)
					IF(ASSOCIATED(rc_op%dsig))DEALLOCATE(rc_op%dsig)
					ALLOCATE(rc_op%sigb(Nz,Nx,Ny_loc))
					ALLOCATE(rc_op%dsig(Nz,Nx,Ny_loc))
				ENDIF
			ELSE
				IF(ASSOCIATED(rc_op%dsig))DEALLOCATE(rc_op%dsig)
				ALLOCATE(rc_op%sigb(Nz,Nx,Ny_loc))
				ALLOCATE(rc_op%dsig(Nz,Nx,Ny_loc))
			ENDIF
			!!$OMP PARALLEL	DEFAULT(SHARED) PRIVATE(Iz)
			!!$OMP DO SCHEDULE(GUIDED)
			DO Iz=1,rc_op%Nz
					rc_op%sigb(Iz,:,:)=bkg%sigma(anomaly%Lnumber(Iz))
			ENDDO
			!!$OMP ENDDO
			!!$OMP END PARALLEL
		ENDIF
	ENDSUBROUTINE

	SUBROUTINE APPLY_RC_E_OP(rc_op)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		INTEGER::Ixy
		REAL(8)::time1,time2,time3,time4
		COMPLEX(REALPARM),POINTER::G(:,:,:,:)
		time1=MPI_WTIME()
		CALL RC_OP_FFTW_FWD(rc_op)
		time2=MPI_WTIME()
		rc_op%counter%mult_fftw=rc_op%counter%mult_fftw+time2-time1
		G=>rc_op%G_E4
		!$OMP PARALLEL	PRIVATE(Ixy) DEFAULT(SHARED)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ixy=1,rc_op%Nx2Ny_loc

			CALL RC_OP_ZGEMV(rc_op,G,REXX,EX,EX,Ixy,C_ONE,C_ZERO)
			CALL RC_OP_ZGEMV(rc_op,G,REXY,EY,EX,Ixy,C_ONE,C_ONE)
			CALL RC_OP_ZGEMV(rc_op,G,REXZ,EZ,EX,Ixy,C_ONE,C_ONE)

			CALL RC_OP_ZGEMV(rc_op,G,REYX,EX,EY,Ixy,C_ONE,C_ZERO)
			CALL RC_OP_ZGEMV(rc_op,G,REYY,EY,EY,Ixy,C_ONE,C_ONE)
			CALL RC_OP_ZGEMV(rc_op,G,REYZ,EZ,EY,Ixy,C_ONE,C_ONE)


			CALL RC_OP_ZGEMV(rc_op,G,REZX,EX,EZ,Ixy,C_ONE,C_ZERO)
			CALL RC_OP_ZGEMV(rc_op,G,REZY,EY,EZ,Ixy,C_ONE,C_ONE)
			CALL RC_OP_ZGEMV(rc_op,G,REZZ,EZ,EZ,Ixy,C_ONE,C_ONE)
		ENDDO
		!$OMP END DO
		!$OMP END  PARALLEL
		time3=MPI_WTIME()
		CALL RC_OP_FFTW_BWD(rc_op)
		time4=MPI_WTIME()
		rc_op%counter%mult_fftw=rc_op%counter%mult_fftw+time4-time3
		rc_op%counter%mult_zgemv=rc_op%counter%mult_zgemv+time3-time2
	ENDSUBROUTINE

	SUBROUTINE APPLY_RC_H_OP(rc_op)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		INTEGER::Ixy
		REAL(8)::time1,time2,time3,time4
		COMPLEX(REALPARM),POINTER::G(:,:,:,:)
		!------ ONLY IN THIS SUBROUTINE !!!!
                INTEGER, PARAMETER :: HX=1
                INTEGER, PARAMETER :: HY=2
                INTEGER, PARAMETER :: HZ=3
		!------ ONLY IN THIS SUBROUTINE !!!!
		time1=MPI_WTIME()
		CALL RC_OP_FFTW_FWD(rc_op)
		time2=MPI_WTIME()
		rc_op%counter%mult_fftw=rc_op%counter%mult_fftw+time2-time1
		G=>rc_op%G_H4
		!$OMP PARALLEL	PRIVATE(Ixy) DEFAULT(SHARED)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ixy=1,rc_op%Nx2Ny_loc

			CALL RC_OP_ZGEMV(rc_op,G,RHXX,EX,HX,Ixy,C_ONE,C_ZERO)
			CALL RC_OP_ZGEMV(rc_op,G,RHXY,EY,HX,Ixy,C_ONE,C_ONE)
			CALL RC_OP_ZGEMV(rc_op,G,RHXZ,EZ,HX,Ixy,C_ONE,C_ONE)

			CALL RC_OP_ZGEMV(rc_op,G,RHYX,EX,HY,Ixy,C_ONE,C_ZERO)
			CALL RC_OP_ZGEMV(rc_op,G,RHYY,EY,HY,Ixy,-C_ONE,C_ONE)
			CALL RC_OP_ZGEMV(rc_op,G,RHYZ,EZ,HY,Ixy,C_ONE,C_ONE)

			CALL RC_OP_ZGEMV(rc_op,G,RHZX,EX,HZ,Ixy,C_ONE,C_ZERO)
			CALL RC_OP_ZGEMV(rc_op,G,RHZY,EY,HZ,Ixy,C_ONE,C_ONE)
		ENDDO
		!$OMP END DO
		!$OMP END  PARALLEL
		time3=MPI_WTIME()
		CALL RC_OP_FFTW_BWD(rc_op)
		time4=MPI_WTIME()
		rc_op%counter%mult_fftw=rc_op%counter%mult_fftw+time4-time3
		rc_op%counter%mult_zgemv=rc_op%counter%mult_zgemv+time3-time2
	ENDSUBROUTINE
	SUBROUTINE RC_OP_ZGEMV(rc_op,G,Tc,c_in,c_out,I,ALPHA,BETA)
			TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G(:,:,:,:)
			INTEGER,INTENT(IN)::c_in,c_out,I,TC
			COMPLEX(REALPARM),INTENT(IN)::ALPHA,BETA
			CALL ZGEMV('T',rc_op%Nz,rc_op%Nr,ALPHA,G(:,:,Tc,I),rc_op%Nz,&
			&rc_op%field_in3(:,c_in,I),ONE,BETA,rc_op%field_out3(:,c_out,I),ONE)
	END SUBROUTINE
ENDMODULE
