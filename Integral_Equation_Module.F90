MODULE INTEGRAL_EQUATION_MODULE
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE Timer_Module 
	USE DATA_TYPES_MODULE
	USE DISTRIBUTED_FFT_MODULE

	USE LOGGER_MODULE
	USE LOCAL_OMP_FFT_MODULE
	IMPLICIT NONE

	PUBLIC

	INTEGER,PARAMETER::GENERAL_MATRIX=1
	INTEGER,PARAMETER::UNIFORM_MATRIX=2

	INTEGER,PARAMETER::COUNTER_ALL=1
	INTEGER,PARAMETER::COUNTER_WT=2
	INTEGER,PARAMETER::COUNTER_LIN=3
	INTEGER,PARAMETER::COUNTER_SQR=4
	INTEGER,PARAMETER::COUNTER_OTHER=5
	INTEGER(LONG_INT),PARAMETER::THREE_64=3
	
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
	

	TYPE IntegralEquationOperator
		INTEGER::Nx,Ny,Nz

		INTEGER::Nx_loc,Ny_loc,Ny_offset
		INTEGER::NxNy_loc
!---------------------------------------------------------------------------------------------------!
		INTEGER(LONG_INT)::local_length
		TYPE(C_PTR)::buff_ptr

		COMPLEX(REALPARM),POINTER:: G_asym(:,:,:,:,:)!G(Nz,Nz,EXZ:EYZ,Nx_loc,Ny_loc)
		COMPLEX(REALPARM),POINTER:: G_asymT(:,:,:,:,:) !for distribudted calc

		COMPLEX(REALPARM),POINTER:: G_symm(:,:,:,:)!G(Nz*(Nz+1)/2,EXX:EZZ,Nx_loc,Ny_loc)
		COMPLEX(REALPARM),POINTER:: G_symmT(:,:,:,:)!for distribudted calc


		INTEGER::matrix_kind
		TYPE(LOCAL_OMP_FFT_DATA),POINTER::LFFT
		COMPLEX(REALPARM),POINTER:: G_symm5(:,:,:,:)!
!---------------------------------------------------------------------------------------------------!

		INTEGER(MPI_CTL_KIND)::ie_comm
		INTEGER(MPI_CTL_KIND)::me
		INTEGER(MPI_CTL_KIND)::master_proc
		LOGICAL::master
!---------------------------------------------------------------------------------------------------!		     
		LOGICAL::fftw_threads_ok
		TYPE(DistributedFourierData)::DFD

		COMPLEX(REALPARM),POINTER::field_in(:,:,:,:),field_out(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_inT(:,:,:,:,:),field_outT(:,:,:,:,:)

!---------------------------------------------------------------------------------------------------!
		REAL(REALPARM),POINTER::dz(:)
		REAL(REALPARM),POINTER ::sqsigb(:)
		COMPLEX(REALPARM),POINTER::csigb(:)
		COMPLEX(REALPARM),POINTER::csiga(:,:,:)
!---------------------------------------------------------------------------------------------------!

		INTEGER(MPI_CTL_KIND)::fgmres_comm
		INTEGER(MPI_CTL_KIND)::fgmres_me

		LOGICAL::real_space
		TYPE(TypeCounter)::counter
	ENDTYPE
CONTAINS
	FUNCTION CalcLocalKernelSize(anomaly,Np) RESULT(local_length)!Quantity of  complex numbers
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		INTEGER(MPI_CTL_KIND),INTENT(IN)::Np 
		INTEGER(LONG_INT)::local_length


		INTEGER(LONG_INT)::symm_length
		INTEGER(LONG_INT)::asym_length
		INTEGER(LONG_INT)::fft_length
		INTEGER(LONG_INT)::NxNyloc
		INTEGER::Nx,Ny,Nz

		Nz=anomaly%Nz
		Nx=anomaly%Nx
		Ny=anomaly%Ny
		IF ( (MOD(2*Ny,Np)/=0).OR.(MOD(Np,2)/=0) ) THEN
			local_length=-1
			RETURN
		ENDIF

		NxNyloc=(2*Nx*Ny)/Np
		symm_length=2*Nz*(Nz+1)*NxNyloc
		asym_length=2*Nz*Nz*NxNyloc
		local_length=symm_length+asym_length
	    ENDFUNCTION


	SUBROUTINE  PrepareIntegralEquation(ie_op,anomaly,comm,kernel_buff,fft_buff_in,fft_buff_out,kernel_size,fft_size,fftw_threads_ok)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm 
		TYPE(C_PTR),INTENT(IN)::kernel_buff,fft_buff_in,fft_buff_out
		INTEGER(C_INTPTR_T) ,INTENT(IN)::kernel_size,fft_size
		LOGICAL,INTENT(IN)::fftw_threads_ok
		INTEGER::Nz
		CALL PrepareIntegralEquationKernel(ie_op,anomaly,comm,kernel_buff,kernel_size)
		CALL PrepareIntegralEquationOperator(ie_op,comm,fftw_threads_ok,fft_buff_in,fft_buff_out,fft_size) 
		Nz=ie_op%Nz
		IF (ie_op%real_space) THEN
			ALLOCATE(ie_op%csigb(Nz))
			ALLOCATE(ie_op%sqsigb(Nz))
		ENDIF
		ALLOCATE(ie_op%dz(Nz))
		ie_op%csiga=>NULL()
		ie_op%dz=anomaly%dz
	END SUBROUTINE

	SUBROUTINE	PrepareIntegralEquationKernel(ie_op,anomaly,comm,buff_ptr,local_length)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm 
		TYPE(C_PTR),INTENT(IN)::buff_ptr
		INTEGER(C_INTPTR_T)::local_length
		INTEGER::Nx,Ny,Nz
		INTEGER::Nx2,Ny2,Nc
		INTEGER(MPI_CTL_KIND)::comm_size,IERROR,me
		INTEGER(MPI_CTL_KIND),PARAMETER::MPI_TWO=2
		INTEGER(MPI_CTL_KIND),PARAMETER::MPI_ONE=1
		Nx=anomaly%Nx
		Ny=anomaly%Ny
		Nz=anomaly%Nz

		ie_op%Nx=anomaly%Nx
		ie_op%Ny=anomaly%Ny
		ie_op%Nz=anomaly%Nz
		ie_op%ie_comm=comm

		CALL MPI_BARRIER(comm,IERROR);
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
			CALL LOGGER('Wrong number of processes. GIEM2G forced to halt!')
			CALL MPI_FINALIZE(IERROR)
			STOP
		ENDIF
!---------------------------------------------------------------------------------------------!
		ie_op%Nx_loc=Nx
		ie_op%Ny_loc=2*Ny/comm_size
		ie_op%Ny_offset=me*ie_op%Ny_loc
		ie_op%NxNy_loc=ie_op%Nx_loc*ie_op%Ny_loc
		ie_op%local_length=local_length
!---------------------------------------------------------------------------------------------!
!		ie_op%matrix_kind=UNIFORM_MATRIX
		ie_op%matrix_kind=GENERAL_MATRIX
		CALL SetTensorPointers(ie_op,buff_ptr)
		IF (ie_op%Ny_offset >= Ny) THEN   
			ie_op%real_space=.FALSE.
			CALL MPI_COMM_SPLIT(ie_op%ie_comm, MPI_TWO, ie_op%me, ie_op%fgmres_comm, IERROR)
		ELSE
			ie_op%real_space=.TRUE.
			CALL MPI_COMM_SPLIT(ie_op%ie_comm, MPI_ONE, ie_op%me, ie_op%fgmres_comm, IERROR)
			CALL MPI_COMM_RANK(ie_op%fgmres_comm, ie_op%fgmres_me, IERROR)
		ENDIF
		CALL MPI_COMM_RANK(ie_op%fgmres_comm, ie_op%fgmres_me, IERROR)

	ENDSUBROUTINE
!--------------------------------------------------------------------------------------------------------------------!	
	  SUBROUTINE  SetTensorPointers(ie_op,buff_ptr)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		TYPE(C_PTR),INTENT(IN)::buff_ptr
		COMPLEX(REALPARM),POINTER:: ptr(:)
		INTEGER::Nx,Ny,Nz,Nx2
		INTEGER::NxNy_loc,N1,N2
		INTEGER::symm_length,asym_length,loc_length
                INTEGER(C_INTPTR_T)::bl
                TYPE(C_PTR)::pin,pout
		CALL C_F_POINTER(buff_ptr,ptr,(/ie_op%local_length/))
                ptr=C_ZERO
		Nz=ie_op%Nz
		Nx=ie_op%Nx
		Ny=ie_op%Ny
		NxNy_loc=ie_op%NxNy_loc
		Nx2=ie_op%Nx_loc-1
		N1=ie_op%Ny_offset
		N2=ie_op%Ny_offset+ie_op%Ny_loc-1

                IF (ie_op%matrix_kind==GENERAL_MATRIX) THEN
                        ie_op%G_asym(1:Nz,1:Nz,A_EXZ:A_EYZ,0:Nx2,N1:N2)=>ptr
                        ie_op%G_asymT(1:Nz,1:Nz,A_EXZ:A_EYZ,1:ie_op%Ny_loc,0:Nx2)=>ptr

                        ptr=>ptr(SIZE(ie_op%G_asym)+1:)
                        ie_op%G_symm(1:Nz*(Nz+1)/2,S_EXX:S_EZZ,0:Nx2,N1:N2)=>ptr
                        ie_op%G_symmT(1:Nz*(Nz+1)/2,S_EXX:S_EZZ,1:ie_op%Ny_loc,0:Nx2)=>ptr
                        ie_op%G_symm5=>NULL()
                ELSEIF (ie_op%matrix_kind==UNIFORM_MATRIX) THEN

!                        ie_op%G_asym(1:2*Nz,1:2,A_EXZ:A_EYZ,0:Nx2,N1:N2)=>ptr
!                        ie_op%G_asym4(1:2*Nz,1:2,A_EXZ:A_EYZ,1:NxNy_loc)=>ptr

                        ptr=>ptr(SIZE(ie_op%G_asym)+1:)

 !                       ie_op%G_symm(1:4*Nz,S_EXX:S_EZZ,0:Nx2,N1:N2)=>ptr
 !                       ie_op%G_symm4(1:4*Nz,S_EXX:S_EZZ,1:NxNy_loc)=>ptr

 !                       ie_op%G_symm5(1:2*Nz,1:2,S_EXX:S_EZZ,1:NxNy_loc)=>ptr
                       ALLOCATE(ie_op%LFFT)
                       bl=CALC_LOCAL_OMP_FFT_SIZE(2*Nz)
                      pin=fftw_alloc_complex(bl) 
                      pout=fftw_alloc_complex(bl) 
                        CALL PREPARE_LOCAL_OMP_FFT(ie_op%LFFT,2*Nz,pin,pout)

                ENDIF
	ENDSUBROUTINE
!--------------------------------------------------------------------------------------------------------------------!
	 SUBROUTINE PrepareIntegralEquationOperator(ie_op,comm,fftw_threads_ok,fft_buff_in,fft_buff_out,buff_len) 
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm 
		LOGICAL, INTENT(IN)::fftw_threads_ok
		TYPE(C_PTR),INTENT(IN)::fft_buff_in
		TYPE(C_PTR),INTENT(IN)::fft_buff_out
		INTEGER(C_INTPTR_T),INTENT(IN)::buff_len
		INTEGER::Nx2,Nz,Ny_loc,Nx2Ny_loc
		INTEGER::NT!,OMP_GET_MAX_THREADS

		IF (fftw_threads_ok) THEN
			NT=OMP_GET_MAX_THREADS()
			CALL FFTW_PLAN_WITH_NTHREADS(NT)
		ENDIF

		Nz=ie_op%Nz
		Nx2=2*ie_op%Nx
		Ny_loc=ie_op%Ny_loc
		Nx2Ny_loc=Nx2*Ny_loc
		CALL	PrepareDistributedFourierData(ie_op%DFD,Nx2,2*ie_op%Ny,3*Nz,comm,&
				&fft_buff_in,fft_buff_out,buff_len)

		ie_op%field_in(1:Nz,1:3,1:Nx2,1:Ny_loc)=>ie_op%DFD%field_out
		ie_op%field_out(1:Nz,1:3,1:Nx2,1:Ny_loc)=>ie_op%DFD%field_in

		ie_op%field_inT(1:Nz,1:2,1:3,1:Ny_loc,0:Nx2/2-1)=>ie_op%DFD%field_out
		ie_op%field_outT(1:Nz,1:2,1:3,1:Ny_loc,0:Nx2/2-1)=>ie_op%DFD%field_in

	ENDSUBROUTINE

	SUBROUTINE IE_OP_FFTW_FWD(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		CALL	CalcDistributedFourier(ie_op%DFD,FFT_FWD)
	ENDSUBROUTINE
	SUBROUTINE IE_OP_FFTW_BWD(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		CALL	CalcDistributedFourier(ie_op%DFD,FFT_BWD)
	ENDSUBROUTINE

	SUBROUTINE IE_OP_FFTW_FWD_PREPARED(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                CALL CalcPreparedForwardFFT(ie_op%DFD)

	ENDSUBROUTINE
	SUBROUTINE IE_OP_FFTW_FWD_FOR_TENSOR(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                CALL CalcForwardIETensorFFT(ie_op%DFD)

	ENDSUBROUTINE
	SUBROUTINE IE_OP_FFTW_BWD_PREPARED(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                CALL CalcPreparedBackwardFFT(ie_op%DFD)
	ENDSUBROUTINE
#define fft_new

#ifdef fft_new
	SUBROUTINE CalcFFTofIETensor(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		IF (ie_op%matrix_kind==GENERAL_MATRIX) THEN
                        CALL CalcFFTofIETensorGeneral(ie_op)
                ELSEIF (ie_op%matrix_kind==UNIFORM_MATRIX) THEN
                        CALL  CalcFFTofIETensorUniform(ie_op)
                ELSE
                ENDIF
        END SUBROUTINE
	SUBROUTINE CalcFFTofIETensorGeneral(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		INTEGER::N,Nz
		REAL(DOUBLEPARM)::time1,time2
		COMPLEX(REALPARM),POINTER::G(:,:,:,:)
		COMPLEX(REALPARM),POINTER::Gout(:,:,:,:)
		TYPE(C_PTR)::ptr
		time1=GetTime()
		ie_op%field_in=C_ZERO
		G=>ie_op%G_symm
		Gout=>ie_op%G_symmT
		Nz=ie_op%Nz
		N=(Nz*(Nz+1))/2

		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EXX,.FALSE.,Gout)
		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EXY,.TRUE.,Gout)
		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EYY,.FALSE.,Gout)
		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EZZ,.FALSE.,Gout)

		ptr=C_LOC(ie_op%G_asymT(1,1,1,1,0))
		N=Nz*Nz
		CALL C_F_POINTER(ptr,G,(/N,2,ie_op%Nx_loc,ie_op%Ny_loc/))
		CALL C_F_POINTER(ptr,Gout,(/N,2,ie_op%Ny_loc,ie_op%Nx_loc/))


		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,A_EXZ,.TRUE.,Gout)
		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,A_EYZ,.FALSE.,Gout)

		time2=GetTime()
		CALL PRINT_CALC_TIME('FFT of  IE tensor has been computed: ',time2-time1)
		ie_op%counter%tensor_fft=time2-time1
	ENDSUBROUTINE

	SUBROUTINE CalcFFTofIETensorUniform(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		INTEGER::N,Nz
		REAL(DOUBLEPARM)::time1,time2
		COMPLEX(REALPARM),POINTER::G(:,:,:,:)
		TYPE(C_PTR)::ptr
		time1=GetTime()
		ie_op%field_in=C_ZERO
		Nz=ie_op%Nz
		N=4*Nz

		ptr=C_LOC(ie_op%G_symm5(1,1,1,1))

		CALL C_F_POINTER(ptr,G,(/N,4,ie_op%Nx_loc,ie_op%Ny_loc/))

!		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EXX,.FALSE.)
!		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EXY,.TRUE.)
!		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EYY,.FALSE.)
!		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,S_EZZ,.FALSE.)

!		ptr=C_LOC(ie_op%G_asym4(1,1,1,1))

!		CALL C_F_POINTER(ptr,G,(/N,2,ie_op%Nx_loc,ie_op%Ny_loc/))


!		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,A_EXZ,.TRUE.)
!		CALL CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,A_EYZ,.FALSE.)

		time2=GetTime()
		CALL PRINT_CALC_TIME('FFT of  IE tensor has been computed: ',time2-time1)
		ie_op%counter%tensor_fft=time2-time1
	ENDSUBROUTINE
#else
	SUBROUTINE CalcFFTofIETensor(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		INTEGER::Iz,Is,Ia
		INTEGER::Is1,Ia1
		INTEGER(MPI_CTL_KIND)::IERROR
		REAL(DOUBLEPARM)::time1,time2
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

		CALL PRINT_CALC_TIME('FFT of  IE tensor has been computed: ',time2-time1)
		ie_op%counter%tensor_fft=time2-time1
	ENDSUBROUTINE
#endif
	SUBROUTINE CALC_FFT_OF_TENSOR_COMPONENT(ie_op,G,N,comp,negative,Gout)
		TYPE(IntegralEquationOperator),TARGET,INTENT(INOUT)::ie_op
		COMPLEX(REALPARM),INTENT(INOUT)::G(1:,1:,1:,1:)
		COMPLEX(REALPARM),INTENT(INOUT)::Gout(1:,1:,1:,1:)
		COMPLEX(REALPARM),POINTER::p_in(:,:,:)
		COMPLEX(REALPARM),POINTER::p_out(:,:,:)
		INTEGER,INTENT(IN)::comp
		LOGICAL,INTENT(IN)::negative
		COMPLEX(REALPARM)::tmp(3*ie_op%Nz)
		INTEGER::N,Nz,M,Nz3,Nx,l
		INTEGER::Iz,Ix,Iy,s
		Nx=ie_op%Nx
		Nz=ie_op%Nz
		Nz3=3*Nz
		M=1
		l=MIN(Nz3,N-M+1)
                p_in(1:Nz3,1:2*Nx,1:ie_op%Ny_loc)=>ie_op%DFD%field_out
                p_out(1:Nz3,1:ie_op%Ny_loc,1:2*Nx)=>ie_op%DFD%field_out
		DO
			DO Iy=1,ie_op%Ny_loc
				DO Ix=1,Nx
				      CALL ZCOPY(l,G(M:,comp,Ix,Iy),ONE,p_in(:,Ix,Iy),ONE)  
				ENDDO
				p_in(:,Nx+1,Iy)=C_ZERO
				DO Ix=1,Nx-1
				      CALL ZCOPY(l,G(M:,comp,Ix+1,Iy),ONE,p_in(:,2*Nx-Ix+1,Iy),ONE)  
				ENDDO
			ENDDO
			IF (negative) THEN
				p_in(:,Nx+2:,:)=-p_in(:,Nx+2:,:)
			ENDIF
			CALL IE_OP_FFTW_FWD_FOR_TENSOR(ie_op) 
                        !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ix,Iy)
                        !$OMP DO SCHEDULE(GUIDED) 
			DO Ix=1,Nx
                                DO Iy=1,ie_op%Ny
        			      CALL ZCOPY(l,p_out(:,Iy,Ix),ONE,Gout(M:,comp,Iy,Ix),ONE) 
                                ENDDO
			ENDDO
                        !$OMP END DO
        		!$OMP END PARALLEL
			DO Iy=1,ie_op%Ny_loc
				CALL ZCOPY(l,p_out(:,Iy,Nx+1),ONE,tmp,ONE) 
				s=1
				DO Ix=1,Nx
				      Gout(M:M+l-1,comp,Iy,Ix)=Gout(M:M+l-1,comp,Iy,Ix)-tmp(1:l)*s
				      s=-s
				ENDDO
			ENDDO
			M=M+l
			IF (M==N+1) EXIT
			l=MIN(Nz3,N-M+1)
		ENDDO
	END SUBROUTINE


	SUBROUTINE DeleteIE_OP(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		ie_op%G_symm=>NULL()
		ie_op%G_symmT=>NULL()
		ie_op%G_asym=>NULL()
		ie_op%G_asymT=>NULL()
		ie_op%field_in=>NULL()
		ie_op%field_inT=>NULL()
		ie_op%field_out=>NULL()
		ie_op%field_outT=>NULL()
		DEALLOCATE(ie_op%dz)
		CALL DeleteDistributedFourierData(ie_op%DFD)
		IF (ie_op%real_space )	DEALLOCATE(ie_op%csigb,ie_op%sqsigb)
	ENDSUBROUTINE

ENDMODULE
