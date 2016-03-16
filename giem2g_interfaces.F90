MODULE GIEM2G_INTERFACES
	USE CONST_MODULE
	USE DISTRIBUTED_FFT_MODULE
	USE Data_Types_Module
        USE LOGGER_MODULE
	USE Timer_Module 
	USE Calc_IE_Tensor_Module
	USE INTEGRAL_EQUATION_MODULE
        USE APPLY_IE_OPERATOR_MODULE 

        IMPLICIT NONE

        TYPE, BIND(C)::GIEM2G_BACKGROUND_TYPE
		INTEGER(C_INT) :: Nl
                TYPE(C_PTR)::csigb
                TYPE(C_PTR)::thickness
        END TYPE

        TYPE, BIND(C)::GIEM2G_ANOMALY_TYPE
                INTEGER(C_INT) :: Nx,Ny,Nz
		REAL(C_DOUBLE)::dx,dy
                TYPE(C_PTR)::z
                TYPE(C_PTR)::dz
        ENDTYPE

        TYPE,  BIND(C)::GIEM2G_DATA_TYPE
                INTEGER(C_LONG_LONG)::tensor_size
                INTEGER(C_LONG_LONG)::ie_kernel_buffer_length
                INTEGER(C_LONG_LONG)::fft_buffers_length
                INTEGER(C_LONG_LONG)::comm64
                TYPE(C_PTR)::giem2g_tensor
                TYPE(C_PTR)::dz
                TYPE(C_PTR)::sqsigb
                TYPE(C_PTR)::csigb
                TYPE(C_PTR)::kernel_buffer
                TYPE(C_PTR)::fft_buffer_in
                TYPE(C_PTR)::fft_buffer_out
        ENDTYPE


	INTERFACE GIEM2G_CALC_LOCAL_KERNEL_SIZE
		MODULE PROCEDURE CalcLocalKernelSize
	END INTERFACE
	INTERFACE GIEM2G_CALC_LOCAL_FFT_SIZE
		MODULE PROCEDURE CalcLocalFFTSize
	END INTERFACE


        TYPE(IntegralEquationOperator),POINTER::tmp_ptr
	
!----------------------------		BINDINGS ------------------------------------------------!	

!------------------------------------------------------------------------------------------------!
CONTAINS

	SUBROUTINE GIEM2G_CALC_DATA_SIZES(giem2g_anomaly,giem2g_data) BIND(C,NAME='giem2g_calc_data_sizes')
                        TYPE(GIEM2G_ANOMALY_TYPE),VALUE,INTENT(IN)::giem2g_anomaly
                        TYPE(GIEM2G_DATA_TYPE),INTENT(INOUT)::giem2g_data
        		INTEGER(MPI_CTL_KIND)::comm,comm_size,IERROR,me
        		INTEGER(C_LONG_LONG)::tensor_size,buff_len,fft_len
        		TYPE (ANOMALY_TYPE)::anomaly
                        TYPE(IntegralEquationOperator)::tmp
                        CALL CREATE_ANOMALY(giem2g_anomaly,anomaly);

                        comm=INT(giem2g_data%comm64,KIND=MPI_CTL_KIND)

                	CALL MPI_COMM_SIZE(comm,comm_size,IERROR) 
                        giem2g_data%tensor_size=STORAGE_SIZE(tmp)/8

                	buff_len=GIEM2G_CALC_LOCAL_Kernel_Size(anomaly,comm_size) 
                	fft_len=GIEM2G_CALC_LOCAL_FFT_SIZE(2*anomaly%Nx,2*anomaly%Ny,3*anomaly%Nz,comm_size) 
                        giem2g_data%ie_kernel_buffer_length=buff_len
                        giem2g_data%fft_buffers_length=fft_len

                	CALL MPI_COMM_RANK(comm,me,IERROR) 
                	LOGGER_MASTER=me
                	CALL InitTimer
			ALLOCATE(tmp_ptr)
			giem2g_data%giem2g_tensor=C_LOC(tmp_ptr)
		ENDSUBROUTINE

        SUBROUTINE GIEM2G_PREPARE_IE_KERNEL(giem2g_anomaly,giem2g_data) BIND(C,NAME='giem2g_prepare_ie_kernel')  
                        TYPE(GIEM2G_ANOMALY_TYPE),VALUE,INTENT(IN)::giem2g_anomaly
                        TYPE(GIEM2G_DATA_TYPE),VALUE,INTENT(IN)::giem2g_data
	        	TYPE(IntegralEquationOperator),POINTER::ie_op
        		TYPE (ANOMALY_TYPE)::anomaly
                        INTEGER(C_INTPTR_T)::kernel_size,fft_size
               		INTEGER(MPI_CTL_KIND)::comm,comm_size,IERROR
			TYPE(C_PTR)::p1,p2,p3
			INTEGER::Nz
                        comm=INT(giem2g_data%comm64,KIND=MPI_CTL_KIND)

                	CALL MPI_COMM_SIZE(comm,comm_size,IERROR) 
                        CALL C_F_POINTER(giem2g_data%giem2g_tensor,ie_op);
                        CALL CREATE_ANOMALY(giem2g_anomaly,anomaly);
                        kernel_size=giem2g_data%ie_kernel_buffer_length
                        fft_size=giem2g_data%fft_buffers_length

			p1=giem2g_data%kernel_buffer!AllocateBuff(kernel_size)
			p2=giem2g_data%fft_buffer_in!AllocateBuff(fft_size)
			p3=giem2g_data%fft_buffer_out!AllocateBuff(fft_size)

        		CALL PrepareIntegralEquationKernel(ie_op,anomaly,comm,p1,kernel_size)
	        	CALL PrepareIntegralEquationOperator(ie_op,comm,.TRUE.,p2,p3,fft_size) 
			Nz=anomaly%Nz
        		IF (ie_op%real_space) THEN
				p1=giem2g_data%csigb
				p3=giem2g_data%sqsigb
        !                        CALL C_F_POINTER(p1,ie_op%csigb,(/Nz/));
  !                              CALL C_F_POINTER(p3,ie_op%sqsigb,(/Nz/));
	 			ALLOCATE(ie_op%csigb(Nz),ie_op%sqsigb(Nz))	
        		ENDIF
			p3=giem2g_data%dz
!                        CALL C_F_POINTER(p3,ie_op%dz,(/ie_op%Nz/));
 			ALLOCATE(ie_op%dz(Nz))	
	        	ie_op%csiga=>NULL()
        		ie_op%dz=anomaly%dz
!                        PRINT*,ie_op%me,ASSOCIATED(ie_op%csigb)
        ENDSUBROUTINE
       

        SUBROUTINE  GIEM2G_CALC_INTEGRAL_GREEN_TENSOR(ie_ptr,giem2g_bkg,giem2g_anomaly,omega)&
                                        & BIND(C,NAME='giem2g_calc_ie_kernel')  
                        TYPE(C_PTR),VALUE,INTENT(IN)::ie_ptr
                        TYPE(GIEM2G_BACKGROUND_TYPE),VALUE,INTENT(IN)::giem2g_bkg
                        TYPE(GIEM2G_ANOMALY_TYPE),VALUE,INTENT(IN)::giem2g_anomaly
                        REAL(C_DOUBLE),VALUE,INTENT(IN)::omega
                        TYPE(BKG_DATA_TYPE)::bkg
                        TYPE(ANOMALY_TYPE)::anomaly
                	TYPE(IntegralEquationOperator),POINTER::int_eq
                        INTEGER::Nz
                        
                        CALL C_F_POINTER(ie_ptr,int_eq)

                        CALL CREATE_BACKGROUND(giem2g_bkg,bkg,omega);


                        CALL CREATE_ANOMALY(giem2g_anomaly,anomaly);
                        Nz=anomaly%Nz
                        ALLOCATE(anomaly%Lnumber(0:Nz));
			anomaly%Lnumber(1:Nz)=GetLayer((anomaly%z(0:Nz-1)+anomaly%z(1:Nz))/2d0,bkg)
			anomaly%Lnumber(0)=GetLayer(anomaly%z(0)*(1d0-1d-7),bkg)
                        IF (ASSOCIATED(int_eq%csigb)) int_eq%csigb=C_ZERO
                        CALL CalcIntegralGreenTensor(int_eq,bkg,anomaly)
                        DEALLOCATE(anomaly%Lnumber)
                        CALL DELETE_BACKGROUND(bkg)
        ENDSUBROUTINE

        SUBROUTINE  GIEM2G_CALC_FFT_OF_GREEN_TENSOR(ie_ptr)&
                                        & BIND(C,NAME='giem2g_calc_fft_of_ie_kernel')  
                        TYPE(C_PTR),VALUE,INTENT(IN)::ie_ptr
                	TYPE(IntegralEquationOperator),POINTER::int_eq
                        
                        CALL C_F_POINTER(ie_ptr,int_eq)
                        CALL CalcFFTofIETensor(int_eq)
!                        CALL PrintTimings(int_eq%DFD)
        ENDSUBROUTINE

        SUBROUTINE GIEM2G_SET_ANOMALY_CONDUCTIVITY(ie_op_ptr,siga_ptr) &
                                        &   BIND(C,NAME='giem2g_set_anomaly_conductivity')
                        TYPE(C_PTR),VALUE,INTENT(IN)::ie_op_ptr,siga_ptr
        		TYPE(IntegralEquationOperator),POINTER::ie_op
                        INTEGER::localshape(3);
                        CALL C_F_POINTER(ie_op_ptr,ie_op)
                        localshape=(/ie_op%Nx,ie_op%Ny_loc,ie_op%Nz/)
        		IF (ie_op%real_space) THEN
                                CALL C_F_POINTER(siga_ptr,ie_op%csiga,localshape)
                        ENDIF
        ENDSUBROUTINE

        SUBROUTINE GIEM2G_APPLY_INTEGRAL_OPERATOR(ie_op_ptr,in_ptr,out_ptr) &
                                        &	BIND(C,NAME='giem2g_apply_ie_operator')   
                        TYPE(C_PTR),VALUE,INTENT(IN)::ie_op_ptr,in_ptr,out_ptr

        		TYPE(IntegralEquationOperator),POINTER::ie_op
	        	COMPLEX(REALPARM),POINTER ::v_in(:)
        		COMPLEX(REALPARM),POINTER ::v_out(:)
        		INTEGER::localsize(1);
                        CALL C_F_POINTER(ie_op_ptr,ie_op)
                        localsize=(/3*ie_op%Nz*ie_op%Nx*ie_op%Ny_loc/)
        		IF (ie_op%real_space) THEN
                                CALL C_F_POINTER(in_ptr,v_in,localsize)
                                CALL C_F_POINTER(out_ptr,v_out,localsize)
				CALL APPLY_PRECOND(ie_op,v_in,v_out)
                        ELSE
				CALL APPLY_IE_ZEROS(ie_op)
                        ENDIF
        ENDSUBROUTINE

        SUBROUTINE GIEM2G_SET_LOGGER(FUNPTR) BIND(C, NAME='giem2g_set_logger')
                        TYPE(C_FUNPTR),VALUE,INTENT(IN)::FUNPTR
                        CALL C_F_PROCPOINTER(FUNPTR, C_LOGGER)
                        LOGGER=>EXTREME_LOGGER
        ENDSUBROUTINE

        SUBROUTINE GIEM2G_DELETE_OPERATOR(ie_op_ptr) BIND(C, NAME='giem2g_delete_ie_operator')
                        TYPE(C_PTR),VALUE,INTENT(IN)::ie_op_ptr
                	TYPE(IntegralEquationOperator),POINTER::ie_op
                        CALL C_F_POINTER(ie_op_ptr,ie_op)
			CALL DeleteIE_OP(ie_op)
!			DEALLOCATE(ie_op)
        ENDSUBROUTINE
        SUBROUTINE CREATE_BACKGROUND(giem2g_bkg,bkg,w)
                        TYPE(GIEM2G_BACKGROUND_TYPE),INTENT(IN)::giem2g_bkg
                        TYPE(BKG_DATA_TYPE),INTENT(INOUT)::bkg
                        REAL(C_DOUBLE),INTENT(IN)::w
                        COMPLEX(C_DOUBLE),POINTER::ptr(:);
                        INTEGER::N,I
                        N=giem2g_bkg%nl;
                        bkg%Nl=N;
                        CALL C_F_POINTER(giem2g_bkg%thickness,bkg%thick,(/N-1/))
			ALLOCATE (bkg%depth(1:N-1))
			bkg%depth(1)=bkg%thick(1)
			DO I=2,N-1
				bkg%depth(I)=bkg%depth(I-1)+bkg%thick(I)
			END DO

                        bkg%sigma=>NULL();
                        bkg%freq=w/PI/2;
                        bkg%omega=w;
                        CALL C_F_POINTER(giem2g_bkg%csigb,ptr,(/N+1/))

                        bkg%csigma(0:N)=>ptr

        		ALLOCATE (bkg%k(0:N),bkg%k2(0:N))
        		bkg%k2=C_IONE*bkg%csigma*MU0*bkg%omega
        		bkg%k=SQRT(bkg%k2)
	        	bkg%iwm=C_IONE*MU0*bkg%omega
        ENDSUBROUTINE





        SUBROUTINE CREATE_ANOMALY(giem2g_anomaly,anomaly)
                        TYPE(GIEM2G_ANOMALY_TYPE),INTENT(IN)::giem2g_anomaly
                        TYPE(ANOMALY_TYPE),INTENT(INOUT)::anomaly
                        REAL(REALPARM),POINTER::ptr(:)
                        anomaly%Nx=giem2g_anomaly%Nx;
                        anomaly%Ny=giem2g_anomaly%Ny;
                        anomaly%Nz=giem2g_anomaly%Nz
                        anomaly%dx=giem2g_anomaly%dx;
                        anomaly%dy=giem2g_anomaly%dy;
                        CALL  C_F_POINTER(giem2g_anomaly%z,ptr,(/anomaly%Nz+1/));
                        anomaly%z(0:anomaly%nz)=>ptr
                        CALL  C_F_POINTER(giem2g_anomaly%dz,anomaly%dz,(/anomaly%Nz/));
                        anomaly%siga=>NULL();
                        anomaly%epsa=>NULL();
        ENDSUBROUTINE
        SUBROUTINE DELETE_BACKGROUND(bkg)
                        TYPE(BKG_DATA_TYPE),INTENT(INOUT)::bkg
                        DEALLOCATE(bkg%k,bkg%k2,bkg%depth)
        ENDSUBROUTINE
END MODULE
