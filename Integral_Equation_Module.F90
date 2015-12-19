MODULE INTEGRAL_EQUATION_MODULE
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE Timer_Module 
	USE DATA_TYPES_MODULE
	USE DISTRIBUTED_FFT_MODULE

	IMPLICIT NONE
	INTEGER,PARAMETER::COUNTER_ALL=1
	INTEGER,PARAMETER::COUNTER_WT=2
	INTEGER,PARAMETER::COUNTER_LIN=3
	INTEGER,PARAMETER::COUNTER_SQR=4
	INTEGER,PARAMETER::COUNTER_OTHER=5
	INTEGER(8),PARAMETER::THREE_64=3
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
		TYPE(DistributedFourierData)::DFD
		COMPLEX(REALPARM),POINTER::field_in4(:,:,:,:),field_out4(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_in3(:,:,:),field_out3(:,:,:)
		LOGICAL::fftw_threads_ok

		INTEGER(8)::N
		INTEGER::Nloc
		COMPLEX(REALPARM),POINTER ::csiga(:,:,:)
		COMPLEX(REALPARM),POINTER ::csigb(:)
		REAL(REALPARM),POINTER ::sqsigb(:)

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
	SUBROUTINE PrepareIntegralEquation(ie_op,anomaly,mcomm,fftw_threads_ok,Nb1)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		INTEGER(MPI_CTL_KIND),INTENT(IN)::mcomm 
		LOGICAL,OPTIONAL,INTENT(IN)::fftw_threads_ok
		INTEGER,OPTIONAL,INTENT(IN)::Nb1
		INTEGER::Nx,Ny,Nz,Nb,Nopt
		INTEGER::Nx2,Ny2,Nc
		INTEGER(8)::Nx64,Ny64,Nc64
		INTEGER(MPI_CTL_KIND)::comm_size,IERROR,me,comm
		INTEGER(MPI_CTL_KIND),PARAMETER::MPI_TWO=2
		INTEGER(MPI_CTL_KIND),PARAMETER::MPI_ONE=1
		INTEGER::NT,OMP_GET_MAX_THREADS
		Nx=anomaly%Nx
		Ny=anomaly%Ny
		Nz=anomaly%Nz
		Nx64=Nx
		Ny64=Ny
		Nc64=3*Nz
		Nc=3*Nz
		Nx2=2*Nx
		Ny2=2*Ny
		IF (PRESENT(Nb1)) THEN
			Nb=Nb1
		ELSE
			Nb=3
		ENDIF
!		CALL MPI_DUP_COMM(mcomm,comm,IERROR)
		comm=mcomm
		ie_op%Nx=Nx
		ie_op%Ny=Ny
		ie_op%Nz=Nz
		ie_op%matrix_comm=comm
		CALL MPI_COMM_RANK(comm, me, IERROR)
		CALL MPI_COMM_SIZE(comm,comm_size,IERROR) 
		ie_op%me=me
		if (me==0) THEN
			ie_op%master=.TRUE.
		ELSE
			ie_op%master=.FALSE.
		ENDIF
		
		ie_op%master_proc=0
		IF (mod(2*Ny,comm_size)/=0) THEN
			IF (ie_op%master) THEN
				PRINT*,'Wrong number of processes. GIEM2G forced to halt!'
			ENDIF
			CALL MPI_FINALIZE(IERROR)
			STOP
		ENDIF
		IF (PRESENT(fftw_threads_ok)) THEN
			ie_op%fftw_threads_ok=fftw_threads_ok
		ELSE
			ie_op%fftw_threads_ok=.FALSE.
		ENDIF
		

		IF (ie_op%fftw_threads_ok) THEN
			NT=OMP_GET_MAX_THREADS()
			CALL FFTW_PLAN_WITH_NTHREADS(NT)
			IF (ie_op%master) PRINT* ,'MULTITHREADED FFTW'
		ENDIF
#ifdef LEGACY_MPI
		Nopt=1
#else
		IF (Nb/=1) THEN
!			CALL TestBlocksNumber(ie_op%DFD,Nx2,Ny2,Nc,comm,Nb,Nopt)
			Nopt=Nb
		ELSE
			Nopt=1
		ENDIF
#endif
		CALL PrepareDistributedFourierData(ie_op%DFD,Nx2,Ny2,Nc,comm,Nopt)
		IF (VERBOSE) THEN
			IF (ie_op%master) THEN
			ENDIF
		ENDIF
		CALL PrepareOperatorIE_OP(ie_op)
		ie_op%N=Nx64*Ny64*Nc64
		ie_op%Nloc=Nx*3*Nz*ie_op%Ny_loc
		IF (ie_op%Ny_offset >= Ny) THEN	  
			ie_op%real_space=.FALSE.
			CALL MPI_COMM_SPLIT(ie_op%matrix_comm, MPI_TWO, ie_op%me, ie_op%fgmres_comm, IERROR)
		ELSE
			ie_op%real_space=.TRUE.
			CALL MPI_COMM_SPLIT(ie_op%matrix_comm, MPI_ONE, ie_op%me, ie_op%fgmres_comm, IERROR)
			CALL MPI_COMM_RANK(ie_op%fgmres_comm, ie_op%fgmres_me, IERROR)
		ENDIF
		CALL MPI_COMM_RANK(ie_op%fgmres_comm, ie_op%fgmres_me, IERROR)
		IF ((ie_op%master).AND.(ie_op%Ny_offset/=0)) THEN
			PRINT*,  'Main process has no zero offset.&
				   & Now it is problem. It will be solved... sometime'
				CALL MPI_FINALIZE(IERROR)
				STOP
		ENDIF
		IF (ie_op%real_space) THEN
			ALLOCATE(ie_op%solution(ie_op%Nloc))
			ALLOCATE(ie_op%initial_guess(ie_op%Nloc))
			ALLOCATE(ie_op%rhs(ie_op%Nloc))
			ALLOCATE(ie_op%csigb(Nz))
			ALLOCATE(ie_op%sqsigb(Nz))
			ie_op%Esol(1:Nz,1:3,1:Nx,1:ie_op%Ny_loc)=>ie_op%solution
			ie_op%E_n(1:Nz,1:3,1:Nx,1:ie_op%Ny_loc)=>ie_op%rhs
		ELSE
			
			ie_op%Esol=>NULL()
			ie_op%E_n=>NULL()
			ie_op%solution=>NULL()
			ie_op%initial_guess=>NULL()
			ie_op%rhs=>NULL()
		ENDIF
		ie_op%csiga=>NULL()
!		ie_op%sigb=>NULL()
!		ie_op%asiga=>NULL()
!		ie_op%sqsigb=>NULL()
!		ie_op%gsig=>NULL()
!		ie_op%dsig=>NULL()
		anomaly%Ny_loc=ie_op%Ny_loc
		ie_op%dz=anomaly%z(1:Nz)-anomaly%z(0:Nz-1)
	ENDSUBROUTINE
!--------------------------------------------------------------------------------------------------------------------!	
	SUBROUTINE PrepareOperatorIE_OP(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		CALL CalcSizesForIE_OP(ie_op)
		CALL AllocateIE_OP(ie_op)
!		CALL CalcFFTWPlansIE_OP(ie_op)
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

		ie_op%Ny_loc=ie_op%DFD%Ny_loc
		ie_op%Ny_offset=ie_op%me*ie_op%Ny_loc
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
		INTEGER::Nx2,Nz,Ny_loc,Nx2Ny_loc
		CALL AllocateIEMatrix(ie_op)
		Nz=ie_op%Nz
		Nx2=ie_op%Nx2
		Ny_loc=ie_op%Ny_loc
		Nx2Ny_loc=ie_op%Nx2Ny_loc
		ie_op%field_in4(1:Nz,1:3,1:Nx2,1:Ny_loc)=>ie_op%DFD%field_out
		ie_op%field_out4(1:Nz,1:3,1:Nx2,1:Ny_loc)=>ie_op%DFD%field_in

		ie_op%field_in3(1:Nz,1:3,1:Nx2Ny_loc)=>ie_op%DFD%field_in
		ie_op%field_out3(1:Nz,1:3,1:Nx2Ny_loc)=>ie_op%DFD%field_out

	ENDSUBROUTINE

	SUBROUTINE IE_OP_FFTW_FWD(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		CALL	CalcDistributedFourier(ie_op%DFD,FFT_FWD)
	ENDSUBROUTINE
	SUBROUTINE IE_OP_FFTW_BWD(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		CALL	CalcDistributedFourier(ie_op%DFD,FFT_BWD)
	ENDSUBROUTINE
	SUBROUTINE CalcFFTofIETensor(ie_op)
		TYPE(IntegralEquation),INTENT(INOUT)::ie_op
		INTEGER::Iz,Is,Ia
		INTEGER::Is1,Ia1
		INTEGER(MPI_CTL_KIND)::IERROR
		   REAL(8)::time1,time2
!		   CALL MPI_BARRIER(ie_op%matrix_comm,IERROR)
		   time1=GetTime()
		Is=1
		Ia=1
		ie_op%field_in4=C_ZERO
		DO Iz=1,ie_op%Nz/3 
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_out4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,1,Ia,:,:)=ie_op%field_out4(:,3,:,:)
			Is=Is+1
			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,1,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia+1,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,1,Is,:,:)=ie_op%field_out4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_out4(:,2,:,:)
			ie_op%G_asym_fftw(:,1,Ia+1,:,:)=ie_op%field_out4(:,3,:,:)
			Ia=Ia+1

			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,2,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia+1,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,2,Is,:,:)=ie_op%field_out4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_out4(:,2,:,:)
			ie_op%G_asym_fftw(:,1,Ia+1,:,:)=ie_op%field_out4(:,3,:,:)
			Ia=Ia+1
			Is=Is+1
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_out4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_out4(:,3,:,:)
			Is=Is+1
			Ia=Ia+1
		ENDDO
		IF (Ia==ie_op%Nz-1) THEN
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)	!Nz-1 Nz-1
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_out4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,1,Ia,:,:)=ie_op%field_out4(:,3,:,:)
			Is=Is+1

			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,1,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia+1,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz Nz
			ie_op%G_symm_fftw(:,1,Is,:,:)=ie_op%field_out4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_out4(:,2,:,:)
			ie_op%G_asym_fftw(:,1,Ia+1,:,:)=ie_op%field_out4(:,3,:,:)
			Ia=Ia+1

			ie_op%field_in4(:,1,:,:)=ie_op%G_symm_fftw(:,2,Is,:,:)
			ie_op%field_in4(:,2,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz Nz
			ie_op%G_symm_fftw(:,2,Is,:,:)=ie_op%field_out4(:,1,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_out4(:,2,:,:)
			Is=Is+1
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,1:2,Is,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz+1 
			ie_op%G_symm_fftw(:,1:2,Is,:,:)=ie_op%field_out4(:,1:2,:,:)
		ELSEIF (Ia==ie_op%Nz) THEN
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,:,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,1,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op)	!Nz Nz
			ie_op%G_symm_fftw(:,:,Is,:,:)=ie_op%field_out4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,1,Ia,:,:)=ie_op%field_out4(:,3,:,:)
			Is=Is+1
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,1:2,Is,:,:)
			ie_op%field_in4(:,3,:,:)=ie_op%G_asym_fftw(:,2,Ia,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) !Nz+1 Nz
			ie_op%G_symm_fftw(:,1:2,Is,:,:)=ie_op%field_out4(:,1:2,:,:)
			ie_op%G_asym_fftw(:,2,Ia,:,:)=ie_op%field_out4(:,3,:,:)
		ELSE
			ie_op%field_in4(:,1:2,:,:)=ie_op%G_symm_fftw(:,1:2,Is,:,:)
			CALL IE_OP_FFTW_FWD(ie_op) 
			ie_op%G_symm_fftw(:,1:2,Is,:,:)=ie_op%field_out4(:,1:2,:,:)
		ENDIF
		time2=GetTime()
		IF (VERBOSE) THEN
			IF (ie_op%master) THEN
				PRINT*, 'FFT of tensor has been computed: ',time2-time1 ,' s'
			ENDIF
		ENDIF
		ie_op%counter%tensor_fft=time2-time1
		CALL PrintTimings(ie_op%dfd)
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
		ie_op%field_in4=>NULL()
		ie_op%field_in3=>NULL()
		ie_op%field_out4=>NULL()
		ie_op%field_out3=>NULL()
		CALL DeleteMatrix(ie_op)
		CALL DeleteDistributedFourierData(ie_op%DFD)
		IF (ie_op%real_space )	DEALLOCATE(ie_op%csigb,ie_op%sqsigb)
	ENDSUBROUTINE

ENDMODULE
