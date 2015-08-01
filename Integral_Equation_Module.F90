MODULE INTEGRAL_EQUATION_MODULE
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
		REAL(8)::mult_fftw_b
		REAL(8)::mult_zgemv
		REAL(8)::apply
		REAL(8)::dotprod
		REAL(8)::solving
		INTEGER::mult_num
		INTEGER::dotprod_num
	ENDTYPE
	TYPE IntegralEquation
		INTEGER::Nx,Ny,Nz
		INTEGER::Nx2,Ny_loc,Ny_offset
		INTEGER::Nx2Ny2
		INTEGER::Nx2Ny_loc
		REAL(REALPARM),POINTER::dz(:)
		COMPLEX(REALPARM),POINTER:: G_asym(:,:,:,:,:)!G(Nz,Nz,EXZ:EYZ,2*Nx,Ny_loc)
		COMPLEX(REALPARM),POINTER:: G_asym_fftw(:,:,:,:,:)! for fft(G)
		COMPLEX(REALPARM),POINTER:: G_asym4(:,:,:,:) !for distribudted calc

		COMPLEX(REALPARM),POINTER:: G_symm(:,:,:,:)!G(Nz*(Nz+1)/2,EXX:EZZ,2*Nx,Ny_loc)
		COMPLEX(REALPARM),POINTER:: G_symm_fftw(:,:,:,:,:)! for fft(G)
		COMPLEX(REALPARM),POINTER:: G_symm4(:,:,:) !for distribudted calc

		TYPE(C_PTR)::pG_asym
		TYPE(C_PTR)::pG_symm
		INTEGER(MPI_CTL_KIND)::matrix_comm
		INTEGER(MPI_CTL_KIND)::me
		INTEGER(MPI_CTL_KIND)::master_proc
		LOGICAL::master
		TYPE(TypeCounter)::counter

		INTEGER (C_INTPTR_T) ::localsize
		TYPE(C_PTR)::p_in,p_out
		COMPLEX(REALPARM),POINTER::field_in4(:,:,:,:),field_out4(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_in3(:,:,:),field_out3(:,:,:)
		TYPE(C_PTR)::planFwd,planBwd
		LOGICAL::fftw_threads_ok

		INTEGER::Nloc
		INTEGER(8)::N
		REAL(REALPARM),POINTER ::siga(:,:,:)
		REAL(REALPARM),POINTER ::sigb(:,:,:)
		REAL(REALPARM),POINTER ::dsig(:,:,:)


		REAL(REALPARM),POINTER ::asiga(:,:,:)
		REAL(REALPARM),POINTER ::sqsigb(:,:,:)
		REAL(REALPARM),POINTER ::gsig(:,:,:)

		COMPLEX(REALPARM),POINTER::Esol(:,:,:,:)
		COMPLEX(REALPARM),POINTER::E_n(:,:,:,:)

		COMPLEX(REALPARM),POINTER::solution(:)
		COMPLEX(REALPARM),POINTER::initial_guess(:)
		COMPLEX(REALPARM),POINTER::rhs(:)
		INTEGER(MPI_CTL_KIND)::fgmres_comm
		INTEGER(MPI_CTL_KIND)::fgmres_me
		LOGICAL::real_space
	  ENDTYPE
	PUBLIC::  IntegralEquation,TypeCounter
	PUBLIC::PrepareIntegralEquation
CONTAINS
	SUBROUTINE PrepareIntegralEquation(int_eq,anomaly,mcomm,fftw_threads_ok)
		TYPE(IntegralEquation),INTENT(INOUT)::int_eq
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		INTEGER(MPI_CTL_KIND),INTENT(IN)::mcomm 
		LOGICAL,OPTIONAL,INTENT(IN)::fftw_threads_ok
		INTEGER::Nx,Ny,Nz
		INTEGER(MPI_CTL_KIND)::comm_size,IERROR,me
		INTEGER(MPI_CTL_KIND),PARAMETER::MPI_TWO=2
		INTEGER(MPI_CTL_KIND),PARAMETER::MPI_ONE=1
		Nx=anomaly%Nx
		Ny=anomaly%Ny
		Nz=anomaly%Nz
		int_eq%Nx=Nx
		int_eq%Ny=Ny
		int_eq%Nz=Nz
		int_eq%matrix_comm=mcomm
		CALL MPI_COMM_RANK(mcomm, me, IERROR)
		CALL MPI_COMM_SIZE(mcomm,comm_size,IERROR) 
		int_eq%me=me
		if (me==0) THEN
			int_eq%master=.TRUE.
		ELSE
			int_eq%master=.FALSE.
		ENDIF
		int_eq%master_proc=0
		IF (mod(2*Ny,comm_size)/=0) THEN
			IF (int_eq%master) THEN
				PRINT*,'Wrong number of processes. GIEM2G forced to halt!'
			ENDIF
			CALL MPI_FINALIZE(IERROR)
			STOP
		ENDIF
		IF (PRESENT(fftw_threads_ok)) THEN
			int_eq%fftw_threads_ok=fftw_threads_ok
		ELSE
			int_eq%fftw_threads_ok=.FALSE.
		ENDIF
		CALL PrepareOperatorIE_OP(int_eq)
		int_eq%N=Nx*Ny*Nz*3
		int_eq%Nloc=Nx*3*Nz*int_eq%Ny_loc
		IF (int_eq%Ny_offset >= Ny) THEN	  
			int_eq%real_space=.FALSE.
			CALL MPI_COMM_SPLIT(int_eq%matrix_comm, MPI_TWO, int_eq%me, int_eq%fgmres_comm, IERROR)
		ELSE
			int_eq%real_space=.TRUE.
			CALL MPI_COMM_SPLIT(int_eq%matrix_comm, MPI_ONE, int_eq%me, int_eq%fgmres_comm, IERROR)
			CALL MPI_COMM_RANK(int_eq%fgmres_comm, int_eq%fgmres_me, IERROR)
		ENDIF
		CALL MPI_COMM_RANK(int_eq%fgmres_comm, int_eq%fgmres_me, IERROR)
		IF ((int_eq%master).AND.(int_eq%Ny_offset/=0)) THEN
			PRINT*,  'Main process has no zero offset.&
				   & Now it is problem. It will be solved... sometime'
				CALL MPI_FINALIZE(IERROR)
				STOP
		ENDIF
		IF (int_eq%real_space) THEN
			ALLOCATE(int_eq%solution(int_eq%Nloc))
			ALLOCATE(int_eq%initial_guess(int_eq%Nloc))
			ALLOCATE(int_eq%rhs(int_eq%Nloc))
			int_eq%Esol(1:Nz,1:3,1:Nx,1:int_eq%Ny_loc)=>int_eq%solution
			int_eq%E_n(1:Nz,1:3,1:Nx,1:int_eq%Ny_loc)=>int_eq%rhs
		ELSE
			
			int_eq%Esol=>NULL()
			int_eq%E_n=>NULL()
			int_eq%solution=>NULL()
			int_eq%initial_guess=>NULL()
			int_eq%rhs=>NULL()
		ENDIF
		int_eq%siga=>NULL()
		int_eq%sigb=>NULL()
		int_eq%asiga=>NULL()
		int_eq%sqsigb=>NULL()
		int_eq%gsig=>NULL()
		int_eq%dsig=>NULL()
		anomaly%Ny_loc=int_eq%Ny_loc
		int_eq%dz=anomaly%z(1:Nz)-anomaly%z(0:Nz-1)
	ENDSUBROUTINE
!--------------------------------------------------------------------------------------------------------------------!	
	SUBROUTINE PrepareOperatorIE_OP(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		CALL CalcSizesForIE_OP(ie_op)
		CALL AllocateIE_OP(ie_op)
		CALL CalcFFTWPlansIE_OP(ie_op)
	ENDSUBROUTINE

	SUBROUTINE CalcSizesForIE_OP(ie_op) 
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		INTEGER(C_INTPTR_T)::tsize8(2)
		INTEGER(C_INTPTR_T)::Nz3
		INTEGER(C_INTPTR_T)::block
		INTEGER(C_INTPTR_T):: CNy,CNy_offset !size and offset for electrical current in Y direction at this process 		  INTEGE
		INTEGER(FFTW_COMM_SIZE)::comm
		REAL(REALPARM)::size_symm,size_asym
		REAL(REALPARM)::tensor_size

		tsize8=(/ie_op%Ny*2,ie_op%Nx*2/)
		block=FFTW_MPI_DEFAULT_BLOCK!(0)
		Nz3=ie_op%Nz*3
		comm=ie_op%matrix_comm
		ie_op%localsize = fftw_mpi_local_size_many(FFTW_TWO,tsize8,Nz3,block,comm, &
			& CNy, CNy_offset)
		ie_op%Ny_loc=INT(CNy,KIND(ie_op%Ny_loc))
		ie_op%Ny_offset=INT(CNy_offset,KIND(ie_op%Ny_offset))
		ie_op%Nx2=ie_op%Nx*2
		ie_op%Nx2Ny_loc=ie_op%Nx2*ie_op%Ny_loc
		ie_op%Nx2Ny2=ie_op%Nx*ie_op%Ny*4
		size_asym=ie_op%Nx2Ny_loc*ie_op%Nz*ie_op%Nz*2.0*2*REALPARM/1024/1024/1024
		size_symm=ie_op%Nx2Ny_loc*ie_op%Nz*(ie_op%Nz+1)*2.0*2*REALPARM/1024/1024/1024
		tensor_size=size_asym+size_symm
		IF (VERBOSE) THEN
			IF (ie_op%master) THEN
				PRINT*,'IE matrix needs', tensor_size, 'Gb per process'
			 ENDIF
		 ENDIF
	ENDSUBROUTINE
	SUBROUTINE AllocateIEMatrix(matrix)
		TYPE(IntegralEquation),INTENT(INOUT)::matrix
		INTEGER(C_INTPTR_T)::length
		INTEGER::shape1(1),Nz,Nx2,Ny_loc,N1,N2
		COMPLEX(REALPARM),POINTER::tmp(:)
		length=matrix%Nx2Ny_loc*matrix%Nz*matrix%Nz*2
		matrix%pG_asym=fftw_alloc_complex(length)
		shape1=(/length/)
		CALL c_f_pointer(matrix%pG_asym,tmp , shape1)
		Nz=matrix%Nz
		Nx2=2*matrix%Nx-1
		N1=matrix%Ny_offset
		N2=matrix%Ny_offset+matrix%Ny_loc-1
		
		matrix%G_asym(1:Nz,1:Nz,A_EXZ:A_EYZ,0:Nx2,N1:N2)=>tmp
		matrix%G_asym_fftw(1:Nz,1:2,1:Nz,0:Nx2,N1:N2)=>tmp
		matrix%G_asym4(1:Nz,1:Nz,A_EXZ:A_EYZ,1:matrix%Nx2Ny_loc)=>tmp

		length=matrix%Nx2Ny_loc*matrix%Nz*(matrix%Nz+1)*2
		matrix%pG_symm=fftw_alloc_complex(length)
		shape1=(/length/)
		CALL c_f_pointer(matrix%pG_symm,tmp , shape1)

		matrix%G_symm(1:Nz*(Nz+1)/2,S_EXX:S_EZZ,0:Nx2,N1:N2)=>tmp
		matrix%G_symm_fftw(1:Nz,1:2,1:Nz+1,0:Nx2,N1:N2)=>tmp
		matrix%G_symm4(1:Nz*(Nz+1)/2,S_EXX:S_EZZ,1:matrix%Nx2Ny_loc)=>tmp

		ALLOCATE(matrix%dz(1:Nz))
	END SUBROUTINE
	SUBROUTINE AllocateIE_OP(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		INTEGER::shape4(4),shape3(3)
		CALL AllocateIEMatrix(ie_op)
		ie_op%p_in=fftw_alloc_complex(ie_op%localsize)
		ie_op%p_out=fftw_alloc_complex(ie_op%localsize)
		shape4=(/ie_op%Nz,3,ie_op%Nx2,ie_op%Ny_loc/)
		shape3=(/ie_op%Nz,3,ie_op%Nx2Ny_loc/)
		CALL c_f_pointer(ie_op%p_in,ie_op%field_in4, shape4)
		CALL c_f_pointer(ie_op%p_in,ie_op%field_in3, shape3)
		CALL c_f_pointer(ie_op%p_out,ie_op%field_out4, shape4)
		CALL c_f_pointer(ie_op%p_out,ie_op%field_out3, shape3)
	ENDSUBROUTINE
	SUBROUTINE CalcFFTWPlansIE_OP(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		INTEGER(C_INTPTR_T)::fftwsize(2)
		INTEGER(C_INTPTR_T)::Nz3
		INTEGER(C_INTPTR_T)::block
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER::omp_get_max_threads,nt
		INTEGER(FFTW_COMM_SIZE)::COMM, FFTW_NT
		REAL(8)::time1,time2
		fftwsize=(/ie_op%Ny*2,ie_op%Nx*2/)
		block=FFTW_MPI_DEFAULT_BLOCK
		Nz3=3*ie_op%Nz
		CALL MPI_BARRIER(ie_op%matrix_comm,IERROR)
		time1=MPI_WTIME()
		IF (ie_op%fftw_threads_ok) THEN
			NT=OMP_GET_MAX_THREADS()
			FFTW_NT=NT
			CALL FFTW_PLAN_WITH_NTHREADS(FFTW_NT)
		ENDIF
		COMM=ie_op%matrix_comm
!		 CALL fftw_set_timelimit(5d0)
				 
		ie_op%planFWD=fftw_mpi_plan_many_dft(FFTW_TWO,fftwsize,Nz3,block,block,&
		&ie_op%field_in4,ie_op%field_in4,COMM, FFTW_FORWARD ,FFTW_MEASURE)!FFTW_PATIENT);
		ie_op%planBWD=fftw_mpi_plan_many_dft(FFTW_TWO,fftwsize,Nz3,block,block,&
		&ie_op%field_out4,ie_op%field_out4,COMM, FFTW_BACKWARD ,FFTW_MEASURE)!FFTW_PATIENT);
		time2=MPI_WTIME()
		ie_op%counter%plans=time2-time1
		IF (VERBOSE) THEN
			IF (ie_op%master) THEN
				PRINT*,'FFTW3 plan calculations:', time2-time1,'s'
			ENDIF
		ENDIF
	ENDSUBROUTINE
	SUBROUTINE IE_OP_FFTW_FWD(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		CALL fftw_mpi_execute_dft(ie_op%planFWD,ie_op%field_in4,ie_op%field_in4)
	ENDSUBROUTINE
	SUBROUTINE IE_OP_FFTW_BWD(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		CALL fftw_mpi_execute_dft(ie_op%planBWD,ie_op%field_out4,ie_op%field_out4)
		ie_op%field_out4=ie_op%field_out4/ie_op%Nx2Ny2
	ENDSUBROUTINE
	SUBROUTINE CalcFFTofIETensor(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		INTEGER::Iz,Is,Ia
		INTEGER::Is1,Ia1
		INTEGER(MPI_CTL_KIND)::IERROR
		   REAL(8)::time1,time2
		   CALL MPI_BARRIER(ie_op%matrix_comm,IERROR)
		   time1=MPI_WTIME()
		Is=1
		Ia=1
		ie_op%field_in4=C_ZERO
		DO Iz=1,ie_op%Nz/3 
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_in4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,1,Ia,:,:)=ie_op%field_in4(:,3,:,:)
			Is=Is+1

			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,1,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia+1,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,1,Is,:,:)=ie_op%field_in4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_in4(:,2,:,:)
			ie_op%G_asym_fftw(:,1,Ia+1,:,:)=ie_op%field_in4(:,3,:,:)
			Ia=Ia+1

			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,2,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia+1,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,2,Is,:,:)=ie_op%field_in4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_in4(:,2,:,:)
			ie_op%G_asym_fftw(:,1,Ia+1,:,:)=ie_op%field_in4(:,3,:,:)
			Ia=Ia+1
			Is=Is+1
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_in4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_in4(:,3,:,:)
			Is=Is+1
			Ia=Ia+1
		ENDDO
		IF (Ia==ie_op%Nz-1) THEN
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)	!Nz-1 Nz-1
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_in4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,1,Ia,:,:)=ie_op%field_in4(:,3,:,:)
			Is=Is+1

			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,1,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia+1,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz Nz
			ie_op%G_symm_fftw(:,1,Is,:,:)=ie_op%field_in4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_in4(:,2,:,:)
			ie_op%G_asym_fftw(:,1,Ia+1,:,:)=ie_op%field_in4(:,3,:,:)
			Ia=Ia+1

			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,2,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz Nz
			ie_op%G_symm_fftw(:,2,Is,:,:)=ie_op%field_in4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_in4(:,2,:,:)
			Is=Is+1
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,1:2,Is,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz+1 
			ie_op%G_symm_fftw(:,1:2,Is,:,:)=ie_op%field_in4(:,1:2,:,:)
		ELSEIF (Ia==ie_op%Nz) THEN
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)	!Nz Nz
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_in4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,1,Ia,:,:)=ie_op%field_in4(:,3,:,:)
			Is=Is+1

			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,1:2,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz+1 Nz
			ie_op%G_symm_fftw(:,1:2,Is,:,:)=ie_op%field_in4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_in4(:,3,:,:)
		ELSE
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,1:2,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz+1 Nz
			ie_op%G_symm_fftw(:,1:2,Is,:,:)=ie_op%field_in4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_in4(:,3,:,:)
		ENDIF
		time2=MPI_WTIME()
		IF (VERBOSE) THEN
			IF (ie_op%master) THEN
				PRINT*, 'FFT of tensor has been computed: ',time2-time1 ,' s'
			ENDIF
		ENDIF
		ie_op%counter%tensor_fft=time2-time1
	ENDSUBROUTINE

	SUBROUTINE DeleteMatrix(matrix)
		TYPE(IntegralEquation),INTENT(INOUT)::matrix
		CALL fftw_free(matrix%pG_symm)
		CALL fftw_free(matrix%pG_asym)
		matrix%G_symm=>NULL()
		matrix%G_symm_fftw=>NULL()
		matrix%G_symm4=>NULL()
		matrix%G_asym=>NULL()
		matrix%G_asym_fftw=>NULL()
		matrix%G_asym4=>NULL()
		DEALLOCATE(matrix%dz)
	ENDSUBROUTINE
	SUBROUTINE DeleteIE_OP(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		CALL fftw_free(ie_op%p_in)
		CALL fftw_free(ie_op%p_out)
		ie_op%field_in4=>NULL()
		ie_op%field_in3=>NULL()
		ie_op%field_out4=>NULL()
		ie_op%field_out3=>NULL()
		CALL fftw_destroy_plan(ie_op%planFwd)
		CALL fftw_destroy_plan(ie_op%planBwd)
		CALL DeleteMatrix(ie_op)
	ENDSUBROUTINE

ENDMODULE
