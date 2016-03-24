MODULE DISTRIBUTED_FFT_MODULE
!VERY VERY C-LIKE MODULE
!THERE ARE A LOT OF POINTERS FOR THE SAME MEMORY
!YES IT'S VERY DIRTY AND UNCOMFORTABLE

!ATTENTION this  Fourier transform pracically is OUT PLACE and DESTROY INPUT

	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE
	USE Timer_Module 
	USE CHECK_MEMORY
	USE LOGGER_MODULE
	IMPLICIT NONE
	PRIVATE
	INTEGER, PARAMETER, PUBLIC :: FFT_FWD=1!FFTW_FORWARD
	INTEGER, PARAMETER, PUBLIC :: FFT_BWD=2!FFTW_BACKWARD
#define print_time
	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_PAT=IOR(FFTW_PATIENT,FFTW_DESTROY_INPUT )
	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_MEA=IOR(FFTW_MEASURE,FFTW_DESTROY_INPUT )

	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_CTL=FFT_MEA

	TYPE PlanXY_TYPE
		TYPE(C_PTR)::planX,planY
	ENDTYPE


	TYPE FFTW_TRANSPOSE_DATA
		TYPE(C_PTR)::plan
		REAL(REALPARM),POINTER::p_in(:),p_out(:)
		REAL(DOUBLEPARM)::plan_time
	ENDTYPE

	TYPE LocalTransposeData
		PROCEDURE(DFD_TRANSPOSE),POINTER,NOPASS::LocalTranspose
		INTEGER::Nopt
	ENDTYPE


	TYPE FFT_COUNTER
		INTEGER::N
		REAL(DOUBLEPARM)::fftx,ffty, x2y_transpose, y2x_transpose
		REAL(DOUBLEPARM)::kernel_total 
		REAL(DOUBLEPARM)::local_transpose(2)
	ENDTYPE

	TYPE DistributedFourierData
		TYPE(C_PTR)::p_in,p_out
		COMPLEX(REALPARM),POINTER::field_in(:),field_out(:)
		COMPLEX(REALPARM),POINTER::field_load_in(:,:,:),field_load_out(:,:,:)

		COMPLEX(REALPARM),POINTER:: field_fft_x_in(:,:,:) !(Nx,Ny_loc,Nm)
		COMPLEX(REALPARM),POINTER:: field_fft_x_out(:,:,:)!(Nx,Ny_loc,Nm)

		COMPLEX(REALPARM),POINTER:: field_fft_y_in(:,:,:) !(K,Ny_loc,Np)
		COMPLEX(REALPARM),POINTER:: field_fft_y_out(:,:,:) !(K,Ny_loc,Np)

		! ATTENTION field_load*(Nc,Nx,Ny) giem2g dimensions order

		INTEGER::Nc,Nx,Ny !Nx,Ny are sizes of 2D FFT transform
		INTEGER::Np !Np is number of processes 
		INTEGER::Ny_loc
		INTEGER::K

		INTEGER(MPI_CTL_KIND)::comm,me
	
		TYPE(FFTW_TRANSPOSE_DATA)::FFTW_TRANSPOSE
		TYPE(LocalTransposeData)::FullLocalTranspose

		TYPE(PlanXY_TYPE)::plan(FFT_FWD:FFT_BWD)

		REAL(DOUBLEPARM)::plans_time(2,FFT_FWD:FFT_BWD)
		REAL(DOUBLEPARM)::local_transpose_time(2)
		TYPE (FFT_COUNTER)::timer(FFT_FWD:FFT_BWD)
	ENDTYPE
	INTERFACE 
		SUBROUTINE DFD_TRANSPOSE (p_in,p_out,Nc,Nxy,Nt)
			USE CONST_MODULE
			COMPLEX(REALPARM),POINTER,INTENT(INOUT)::p_in(:,:),p_out(:,:)
			INTEGER,INTENT(IN)::Nxy,Nc,Nt
		ENDSUBROUTINE
	ENDINTERFACE

	PUBLIC:: DistributedFourierData! CalcDistributedFourier
	PUBLIC:: PrepareDistributedFourierData,DeleteDistributedFourierData
        PUBLIC:: CalcPreparedForwardFFT
        PUBLIC:: CalcPreparedBackwardFFT
        PUBLIC:: CalcForwardIETensorFFT
	PUBLIC:: CalcLocalFFTSize,PrintTimings,DROP_COUNTER
	CONTAINS
	SUBROUTINE  PrepareDistributedFourierData(DFD,Nx,Ny,Nc,comm,fft_buff_in,fft_buff_out,buff_length)
		TYPE (DistributedFourierData),TARGET,INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::Nx,Ny,Nc
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		TYPE(C_PTR),INTENT(IN)::fft_buff_in
		TYPE(C_PTR),INTENT(IN)::fft_buff_out
		INTEGER(C_INTPTR_T),OPTIONAL,INTENT(IN)::buff_length
		INTEGER::Ny_loc
		INTEGER::M,Lblock,My
		INTEGER(MPI_CTL_KIND)::Np
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER(C_INTPTR_T)::Nbuff
		INTEGER::Nbuff32
		INTEGER::Ix,Iy,Ib,Ic
		INTEGER::Lxm,Lkp,Ix0,Iy0
!-----------------------------------------------------------------------------!
		CALL MPI_COMM_SIZE(comm,Np,IERROR) 
		CALL MPI_COMM_RANK(comm,DFD%me,IERROR) 
		M=MODULO(Ny,Np)
		IF (M/=0) THEN
			CALL LOGGER('NOT IMPLEMENTED!! Wrong  Ny and number of processes ')
			RETURN
		ENDIF
		IF (Nc*Nx<Np) THEN
			CALL LOGGER('6*Nz*Nx <Np unsupported anomaly shape, GIEM2G FORCED HALTED')
			CALL MPI_FINALIZE(IERROR)
			STOP
			RETURN
		ENDIF
		IF (PRESENT(buff_length)) THEN
			CALL PrepareDFDLocal(DFD,Nx,Ny,Nc,Np,fft_buff_in, fft_buff_out,buff_length)
		ELSE
			CALL PrepareDFDLocal(DFD,Nx,Ny,Nc,Np,fft_buff_in, fft_buff_out,buff_length)
		ENDIF
		DFD%comm=comm
		CALL CreateAll2AllPlan(DFD)
#ifdef print_time
		CALL PRINT_BORDER
		CALL LOGGER('Block distributed FFT plan calculations at 0 process:')
		CALL PRINT_CALC_TIME('Forward along X:', DFD%plans_time(1,FFT_FWD))
		CALL PRINT_CALC_TIME('Forward along Y:', DFD%plans_time(2,FFT_FWD))
		CALL PRINT_CALC_TIME('Backward along X:', DFD%plans_time(1,FFT_BWD))
		CALL PRINT_CALC_TIME('Backward along Y:', DFD%plans_time(2,FFT_BWD))
		CALL PRINT_CALC_TIME('Optimal local transpose:',DFD%local_transpose_time(1))
		CALL PRINT_CALC_TIME('Local transpose plan:',DFD%local_transpose_time(2))
		CALL PRINT_CALC_TIME('FFTW Transpose plan:',DFD%FFTW_TRANSPOSE%plan_time)
		CALL PRINT_BORDER
#endif
	ENDSUBROUTINE


	SUBROUTINE CreateBlockPlans(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::pNx(1),pNy(1)
		INTEGER::dx,dy,M
		INTEGER::x_shape(3),y_shape(3)
		REAL(DOUBLEPARM)::time1,time2

		pNx(1)=DFD%Nx
		dx=DFD%Nx
		M=DFD%Nc*DFD%Ny_loc

		time1=GetTime()
		DFD%plan(FFT_FWD)%planX=fftw_plan_many_dft(ONE,pNx, M,DFD%field_fft_x_in,pNx,ONE,dx,&
								DFD%field_fft_x_out,pNx,ONE,dx,&
								&FFTW_FORWARD, FFT_CTL)
		time2=GetTime()
		DFD%plans_time(1,FFT_FWD)=time2-time1
		DFD%plan(FFT_BWD)%planX=fftw_plan_many_dft(ONE,pNx, M,DFD%field_fft_x_in,pNx,ONE,dx,&
								DFD%field_fft_x_out,pNx,ONE,dx,&
								&FFTW_BACKWARD, FFT_CTL)
		time1=GetTime()
		DFD%plans_time(1,FFT_BWD)=time1-time2

		pNy(1)=DFD%Ny_loc*DFD%Np
		M=DFD%K

		DFD%plan(FFT_FWD)%planY=fftw_plan_many_dft(ONE,pNy, M,DFD%field_fft_y_in,pNy,M,ONE,&
								DFD%field_fft_y_out,pNy,M,ONE,&
								&FFTW_FORWARD, FFT_CTL)

		time2=GetTime()
		DFD%plans_time(2,FFT_FWD)=time2-time1

		DFD%plan(FFT_BWD)%planY=fftw_plan_many_dft(ONE,pNy, M,DFD%field_fft_y_in,pNy,M,ONE,&
								DFD%field_fft_y_out,pNy,M,ONE,&
								&FFTW_BACKWARD, FFT_CTL)
		time1=GetTime()
		DFD%plans_time(2,FFT_BWD)=time1-time2
	ENDSUBROUTINE
	
	SUBROUTINE CreateAll2AllPlan(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		REAL(REALPARM),POINTER::r_in(:),r_out(:)
		INTEGER::M,I
		INTEGER(C_INTPTR_T)::N,Np
		INTEGER(C_INTPTR_T),PARAMETER::ONE=1
		REAL(REALPARM)::s1,s2
		REAL(DOUBLEPARM)::t1,t2
                INTEGER::NT,OMP_GET_MAX_THREADS
		N=TWO*DFD%K*DFD%Ny_loc
		Np=DFD%Np

		M=N*Np

		CALL c_f_pointer(DFD%p_out,r_in,(/M/))
		CALL c_f_pointer(DFD%p_in,r_out,(/M/))
		t1=GetTime()	
		DFD%FFTW_TRANSPOSE%plan = fftw_mpi_plan_many_transpose(Np, Np, N, ONE, ONE, r_in, r_out, DFD%comm,&
		& FFT_PAT);
		t2=GetTime();
		
		DFD%FFTW_TRANSPOSE%plan_time=t2-t1
		DFD%FFTW_TRANSPOSE%p_in=>r_in
		DFD%FFTW_TRANSPOSE%p_out=>r_out
	ENDSUBROUTINE

	SUBROUTINE CalcPreparedForwardFFT(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		CALL ProcessDistributedFourierKernelSync(DFD,FFT_FWD)
		CALL FinalTransposeRepack(DFD)
	END SUBROUTINE
	
	SUBROUTINE CalcForwardIETensorFFT(DFD,sy)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		COMPLEX(REALPARM),INTENT(IN)::sy
		CALL InitialTransposeIE(DFD)
		CALL ProcessDistributedFourierKernelSyncMirror(DFD,FFT_FWD,sy)
		CALL FinalTransposeIE(DFD)
	END SUBROUTINE

	SUBROUTINE CalcPreparedBackwardFFT(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		CALL InitialTransposeRepack(DFD)
		CALL ProcessDistributedFourierKernelSync(DFD,FFT_BWD)
	END SUBROUTINE
 
!----------------------------------------------------------------------------------------------------------------------------------!

	SUBROUTINE ProcessDistributedFourierKernelSync(DFD,FFT_DIR)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
		INTEGER (MPI_CTL_KIND)::IERROR,comm
		REAL(REALPARM),POINTER::p_send(:),p_recv(:)
		REAL(DOUBLEPARM)::time1,time2
		TYPE(C_PTR)::cp
#ifndef performance_test
		time1=GetTime()
#endif
		DFD%timer(FFT_DIR)%N=DFD%timer(FFT_DIR)%N+1
		p_send=>DFD%FFTW_TRANSPOSE%p_in
		p_recv=>DFD%FFTW_TRANSPOSE%p_out
		CALL DistributedFourierX(DFD,FFT_DIR)

		CALL BlockTransposeXToY(DFD)

		CALL fftw_mpi_execute_r2r(DFD%FFTW_TRANSPOSE%plan,p_send,p_recv)!All2All

		CALL DistributedFourierY(DFD,FFT_DIR)

		CALL  fftw_mpi_execute_r2r(DFD%FFTW_TRANSPOSE%plan,p_send,p_recv)!All2All

		CALL BlockTransposeYToX(DFD)

#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_DIR)%kernel_total=DFD%timer(FFT_DIR)%kernel_total+time2-time1
#endif
	END SUBROUTINE

	
	SUBROUTINE ProcessDistributedFourierKernelSyncMirror(DFD,FFT_DIR,s)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
                COMPLEX(REALPARM),INTENT(IN)::s
		INTEGER (MPI_CTL_KIND)::IERROR,comm
		REAL(REALPARM),POINTER::p_send(:),p_recv(:)
		REAL(DOUBLEPARM)::time1,time2
		TYPE(C_PTR)::cp
#ifndef performance_test
		time1=GetTime()
#endif
		DFD%timer(FFT_DIR)%N=DFD%timer(FFT_DIR)%N+1
		p_send=>DFD%FFTW_TRANSPOSE%p_in
		p_recv=>DFD%FFTW_TRANSPOSE%p_out
		CALL DistributedFourierX(DFD,FFT_DIR)

		CALL BlockTransposeXToY(DFD)

		CALL fftw_mpi_execute_r2r(DFD%FFTW_TRANSPOSE%plan,p_send,p_recv)!All2All

		CALL DistributedFourierYMirror(DFD,FFT_DIR,s)

		CALL  fftw_mpi_execute_r2r(DFD%FFTW_TRANSPOSE%plan,p_send,p_recv)!All2All

		CALL BlockTransposeYToX(DFD)

#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_DIR)%kernel_total=DFD%timer(FFT_DIR)%kernel_total+time2-time1
#endif
	END SUBROUTINE

	SUBROUTINE DistributedFourierY(DFD,FFT_DIR)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
		COMPLEX(REALPARM),POINTER::pin(:,:,:),pout(:,:,:)
		COMPLEX(REALPARM),POINTER::p1(:,:),p2(:,:)
		TYPE(C_PTR)::plan
		REAL(DOUBLEPARM)::time1,time2
                INTEGER::K,Ny,Ny2,Iy
#ifndef performance_test
		time1=GetTime()
#endif
		pin=>DFD%field_fft_y_in
		pout=>DFD%field_fft_y_out
		plan=DFD%plan(FFT_DIR)%planY

                K=DFD%K
                Ny=DFD%Ny
                Ny2=Ny/2
                IF (FFT_DIR==FFT_BWD) THEN
                        p1(1:K,0:Ny-1)=>DFD%field_in
                        p2(1:K,0:Ny-1)=>DFD%field_out

                        DO Iy=1,Ny2-1
                                p1(:,Ny-Iy)=p2(:,Iy+Ny2)
                        ENDDO
                        p2(:,Ny2+1:)=p1(:,Ny2+1:)
                        p2(:,Ny2)=C_ZERO
                ENDIF
		CALL fftw_execute_dft(plan,pin,pout)
                IF (FFT_DIR==FFT_FWD) THEN
                        p1(1:K,0:Ny-1)=>DFD%field_in
                        p2(1:K,0:Ny-1)=>DFD%field_out
                        p2(:,Ny2)=p1(:,0)
                        DO Iy=1,Ny2-1
                                p2(:,Iy+Ny2)=p1(:,Ny-Iy)
                        ENDDO
                        p1(:,Ny2:)=p2(:,Ny2:)
                ENDIF
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_DIR)%ffty=DFD%timer(FFT_DIR)%ffty+time2-time1
#endif
	END SUBROUTINE

	SUBROUTINE DistributedFourierYMirror(DFD,FFT_DIR,s)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
                COMPLEX(REALPARM),INTENT(IN)::s
		COMPLEX(REALPARM),POINTER::pin(:,:,:),pout(:,:,:)
		COMPLEX(REALPARM),POINTER::ptr(:,:)
		TYPE(C_PTR)::plan
                COMPLEX(REALPARM)::t,tmp(DFD%K)
		REAL(DOUBLEPARM)::time1,time2
                INTEGER::K,Ny,Ny2,Iy

#ifndef performance_test
		time1=GetTime()
#endif
                K=DFD%K
                Ny=DFD%Ny
                Ny2=Ny/2
                ptr(1:K,0:Ny-1)=>DFD%field_out
		pin=>DFD%field_fft_y_in
		pout=>DFD%field_fft_y_out
		plan=DFD%plan(FFT_DIR)%planY
                DO Iy=1,Ny2-1
                        ptr(:,Ny-Iy)=s*ptr(:,Iy)
                ENDDO
                ptr(:,Ny2)=C_ZERO
		CALL fftw_execute_dft(plan,pin,pout)
                ptr(1:K,1:Ny)=>DFD%field_in
		t=C_ONE
                tmp=ptr(:,Ny2+1)
		DO Iy=1,Ny
			ptr(:,Iy)=ptr(:,Iy)-t*tmp
			t=-t
		ENDDO
                ptr(:,Ny2+1:)=ptr(:,1:Ny2)
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_DIR)%ffty=DFD%timer(FFT_DIR)%ffty+time2-time1
#endif
	END SUBROUTINE
	SUBROUTINE DistributedFourierX(DFD,FFT_DIR)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
		COMPLEX(REALPARM),POINTER::pin(:,:,:),pout(:,:,:)
		TYPE(C_PTR)::plan
		REAL(DOUBLEPARM)::time1,time2
#ifndef performance_test
		time1=GetTime()
#endif
		pin=>DFD%field_fft_x_in
		pout=>DFD%field_fft_x_out
		plan=DFD%plan(FFT_DIR)%planX
		CALL fftw_execute_dft(plan,pin,pout)
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_DIR)%fftx=DFD%timer(FFT_DIR)%fftx+time2-time1
#endif
	END SUBROUTINE

	SUBROUTINE InitialTransposeIE(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Ix,Iy,Ic,l
		INTEGER::Nx,Ny,Nc,Nt,M
		COMPLEX(REALPARM),POINTER::p_in(:,:,:),p_out(:,:,:)
		REAL(DOUBLEPARM)::time1,time2
#ifndef performance_test
		time1=GetTime()
#endif
		Nx=DFD%Nx
		Ny=DFD%Ny_loc
		Nc=DFD%Nc
		p_in(1:Nc,1:Ny,1:Nx)=>DFD%field_out
		p_out(1:Nx,1:Ny,1:Nc)=>DFD%field_in

		!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ix,Iy,Ic)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ic=1,Nc
			DO Iy=1,Ny
                                DO Ix=1,Nx
					p_out(Ix,Iy,Ic)=p_in(Ic,Iy,Ix)
                                ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_FWD)%local_transpose(1)=time2-time1+DFD%timer(FFT_FWD)%local_transpose(1)
#endif
	ENDSUBROUTINE

	SUBROUTINE FinalTransposeIE(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Ix,Iy,Ic,l
		INTEGER::Nx,Ny,Nc,Nt,M
		COMPLEX(REALPARM),POINTER::p_in(:,:,:),p_out(:,:,:)
		REAL(DOUBLEPARM)::time1,time2
#ifndef performance_test
		time1=GetTime()
#endif
		Nx=DFD%Nx
		Ny=DFD%Ny_loc
		Nc=DFD%Nc
		p_in(1:Nx,1:Ny,1:Nc)=>DFD%field_in
		p_out(1:Nc,1:Ny,1:Nx)=>DFD%field_out

		!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ix,Iy,Ic)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ix=1,Nx
			DO Iy=1,Ny
                                DO Ic=1,Nc
					p_out(Ic,Iy,Ix)=p_in(Ix,Iy,Ic)
                                ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_FWD)%local_transpose(2)=time2-time1+DFD%timer(FFT_FWD)%local_transpose(2)
#endif
	ENDSUBROUTINE

	SUBROUTINE FinalTransposeRepack(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Ix,Iy,Ic,l,Iz
		INTEGER::Nx,Ny,Nc,Nt,M,Nz,Nx2
		COMPLEX(REALPARM),POINTER::p_in(:,:,:,:),p_out(:,:,:,:,:)
		REAL(DOUBLEPARM)::time1,time2
#ifndef performance_test
		time1=GetTime()
#endif
		Nx=DFD%Nx
		Ny=DFD%Ny_loc
		Nc=DFD%Nc
                Nz=Nc/3
                Nx2=Nx/2
		p_in(0:Nx-1,1:Ny,1:Nz,1:3)=>DFD%field_in
		p_out(1:Nz,1:2,1:3,1:Ny,0:Nx2-1)=>DFD%field_out

		DO Iy=1,Ny
                         DO Ic=1,3
                                DO Iz=1,Nz
        				p_out(Iz,1,Ic,Iy,0)=p_in(0,Iy,Iz,Ic)
        				p_out(Iz,2,Ic,Iy,0)=p_in(Nx2,Iy,Iz,Ic)
                                ENDDO
                         ENDDO
        	ENDDO
		!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ix,Iy,Iz,Ic)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ix=1,Nx2-1
			DO Iy=1,Ny
                                DO Ic=1,3
                                        DO Iz=1,Nz
                                                p_out(Iz,1,Ic,Iy,Ix)=p_in(Ix,Iy,Iz,Ic)
                                        ENDDO
                                        DO Iz=1,Nz
                                                p_out(Iz,2,Ic,Iy,Ix)=p_in(Nx-Ix,Iy,Iz,Ic)
                                        ENDDO
                                ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_FWD)%local_transpose(2)=time2-time1+DFD%timer(FFT_FWD)%local_transpose(2)
#endif
	ENDSUBROUTINE

	SUBROUTINE InitialTransposeRepack(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Ix,Iy,Ic,l,Iz
		INTEGER::Nx,Ny,Nc,Nt,M,Nz,Nx2
		COMPLEX(REALPARM),POINTER::p_in(:,:,:,:,:),p_out(:,:,:,:)
		REAL(DOUBLEPARM)::time1,time2
#ifndef performance_test
		time1=GetTime()
#endif
		Nx=DFD%Nx
		Ny=DFD%Ny_loc
		Nc=DFD%Nc
                Nz=Nc/3
                Nx2=Nx/2
		p_in(1:Nz,1:2,1:3,1:Ny,0:Nx2-1)=>DFD%field_in
		p_out(0:Nx-1,1:Ny,1:Nz,1:3)=>DFD%field_out

                DO Ic=1,3
                        DO Iz=1,Nz
                                DO Iy=1,Ny
                                        p_out(0,Iy,Iz,Ic)=p_in(Iz,1,Ic,Iy,0)
                                        p_out(Nx2,Iy,Iz,Ic)=C_ZERO
                                ENDDO
                         ENDDO
        	ENDDO
		!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ix,Iy,Iz,Ic)
		!$OMP DO SCHEDULE(GUIDED) 
                DO Iz=1,Nz
                        DO Iy=1,Ny
                                DO Ix=1,Nx2-1
                                        p_out(Ix,Iy,Iz,1) =  p_in(Iz,1,1,Iy,Ix)
                                ENDDO
                                DO Ix=Nx2+1,Nx-1
                                        p_out(Ix,Iy,Iz,1) =  p_in(Iz,2,1,Iy,Nx-Ix)
                                ENDDO
                        ENDDO
                ENDDO
		!$OMP ENDDO
		!$OMP DO SCHEDULE(GUIDED) 
                DO Iz=1,Nz
                        DO Iy=1,Ny
                                DO Ix=1,Nx2-1
                                        p_out(Ix,Iy,Iz,2) =  p_in(Iz,1,2,Iy,Ix)
                                ENDDO
                                DO Ix=Nx2+1,Nx-1
                                        p_out(Ix,Iy,Iz,2) =  p_in(Iz,2,2,Iy,Nx-Ix)
                                ENDDO
                        ENDDO
                ENDDO
		!$OMP ENDDO
		!$OMP DO SCHEDULE(GUIDED) 
                DO Iz=1,Nz
                        DO Iy=1,Ny
                                DO Ix=1,Nx2-1
                                        p_out(Ix,Iy,Iz,3) =  p_in(Iz,1,3,Iy,Ix)
                                ENDDO
                                DO Ix=Nx2+1,Nx-1
                                        p_out(Ix,Iy,Iz,3) =  p_in(Iz,2,3,Iy,Nx-Ix)
                                ENDDO
                        ENDDO
                ENDDO
		!$OMP ENDDO
                !$OMP WORKSHARE
                DFD%field_in=DFD%field_out
                !$OMP END WORKSHARE
		!$OMP END PARALLEL
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_BWD)%local_transpose(1)=time2-time1+DFD%timer(FFT_BWD)%local_transpose(1)
#endif
	ENDSUBROUTINE
	SUBROUTINE BlockTransposeXToY(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Ip,Iy,Ik,K
		INTEGER::l,Ix,Im,l0,Iy0,m0,Ny_loc
		INTEGER::Nx,Lkp,Np
		INTEGER::x_shape(3),y_shape(3)
		INTEGER::M
		COMPLEX(REALPARM),POINTER::p_in(:),p_out(:,:,:)
		REAL(DOUBLEPARM)::time1,time2
#ifndef performance_test
		time1=GetTime()
#endif
		K=DFD%K
		Ny_loc=DFD%Ny_loc
		Nx=DFD%Nx
		Np=DFD%Np
		p_out=>DFD%field_fft_y_out
		p_in=>DFD%field_out
		M=SIZE(p_in) !length of data
		!Linear writting, nonlinear reading
		!$OMP PARALLEL	DEFAULT(SHARED),PRIVATE(Ip,Iy,Ik,l0,m0,l)
		!$OMP DO SCHEDULE(GUIDED)	
		DO Ip=1,Np
			DO Iy=1,Ny_loc
				DO Ik=1,K
					l0=Ik+(Ip-1)*K !place in XZ
					m0=MODULO(l0-1,Nx)+1	   !Ix=m0
					m0=l0-m0
					l=l0+(Iy-1)*Nx+m0*(Ny_loc-1)! Ix+(Iy-1)Nx+(Ic-1)NxNyLoc
					l=(l+M-ABS(M-l))/2 !if l>M we are in "addition" zone
					p_out(Ik,Iy,Ip)=p_in(l)
				ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL

#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_FWD)%x2y_transpose=time2-time1+DFD%timer(FFT_FWD)%x2y_transpose
#endif
	END SUBROUTINE

	SUBROUTINE BlockTransposeYToX(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Im,Iy,Ix
		INTEGER::l0,K,I,Il0,m0,l
		INTEGER::Ny_loc,Nx,Np,M,Nm
		INTEGER::x_shape(3),y_shape(3)
		INTEGER::sK,dK,sI,dI,sL,dL,I2	
		COMPLEX(REALPARM),POINTER::p_in(:),p_out(:,:,:)
		REAL(DOUBLEPARM)::time1,time2
#ifndef performance_test
		time1=GetTime()
#endif

		Nm=DFD%Nc
		Ny_loc=DFD%Ny_loc
		Nx=DFD%Nx
		Np=DFD%Np
		K=DFD%K		
		p_out=>DFD%field_fft_x_in
		p_in=>DFD%field_out
		!Linear writting, nonlinear reading
		!$OMP PARALLEL	DEFAULT(SHARED),PRIVATE(Im,Iy,Ix,l0,m0,l)
		!$OMP DO SCHEDULE(GUIDED)	
		DO Im=1,Nm
			DO Iy=1,Ny_loc
				DO Ix=1,Nx
					l0=Ix+(Im-1)*Nx !place in XZ
					m0=MODULO(l0-1,K)+1 !Ik=m0
					m0=l0-m0
					l=l0+(Iy-1)*K+m0*(Ny_loc-1)! Ix+(Iy-1)Nx+(Ic-1)NxNyLoc
					p_out(Ix,Iy,Im)=p_in(l)
				ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
#ifndef performance_test
		time2=GetTime()
		DFD%timer(FFT_FWD)%y2x_transpose=time2-time1+DFD%timer(FFT_FWD)%y2x_transpose
#endif
	END SUBROUTINE
	SUBROUTINE DeleteDistributedFourierData(DFD)
		TYPE (DistributedFourierData),TARGET,INTENT(INOUT)::DFD
		INTEGER::Ib
		INTEGER(MPI_CTL_KIND)::IERROR
		CALL fftw_destroy_plan(DFD%plan(FFT_FWD)%planX);	
		CALL fftw_destroy_plan(DFD%plan(FFT_FWD)%planY);	
		CALL fftw_destroy_plan(DFD%plan(FFT_BWD)%planX);	
		CALL fftw_destroy_plan(DFD%plan(FFT_BWD)%planY);
		DFD%field_in=>NULL()
		DFD%field_out=>NULL()
		DFD%field_load_in=>NULL()
		DFD%field_load_out=>NULL()
	ENDSUBROUTINE
	
	FUNCTION CalcLocalFFTSize(Nx,Ny,Nc,Np) RESULT(buff_len)
		INTEGER,INTENT(IN)::Nx,Ny,Nc,Np
		INTEGER::Ny_loc
		INTEGER::M,My,K
		INTEGER(C_INTPTR_T)::buff_len

		Ny_loc=Ny/Np
		My=Nx*Nc
		M=MODULO(My,Np)
			
		IF (M==0) THEN
			K=My/Np
		ELSE
			K=(My-M)/Np+1
		ENDIF

		buff_len=K*Np*Ny_loc
	    ENDFUNCTION

	SUBROUTINE PrepareDFDLocal(DFD,Nx,Ny,Nc,Np,fft_buff_in,fft_buff_out,buff_len)
		TYPE (DistributedFourierData),TARGET,INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::Nx,Ny,Nc,Np
		TYPE(C_PTR),INTENT(IN)::fft_buff_in
		TYPE(C_PTR),INTENT(IN)::fft_buff_out
		INTEGER(C_INTPTR_T),INTENT(IN),OPTIONAL::buff_len
		INTEGER(C_INTPTR_T)::buff_size
		COMPLEX(REALPARM),POINTER::ptr(:)
		INTEGER::Ny_loc
		INTEGER::M,My,K
		INTEGER::Ix,Iy,Ib,Ic
		INTEGER(MPI_CTL_KIND)::IERROR


		Ny_loc=Ny/Np
		My=Nx*Nc
		M=MODULO(My,Np)
			
		IF (M==0) THEN
			K=My/Np
		ELSE
			K=(My-M)/Np+1
		ENDIF
		buff_size=K*Np*Ny_loc
		IF (PRESENT(buff_len)) THEN
			IF (buff_size/=buff_len) THEN
				CALL LOGGER("Wrong buffer size for Distributed FFT")
				CALL LOGGER('GIEM2G FORCED HALTED')
				CALL MPI_FINALIZE(IERROR)
				STOP
			ENDIF
		ENDIF
		CALL c_f_pointer(fft_buff_in,ptr,(/buff_size/))
		DFD%field_in(1:buff_len)=>ptr
		CALL c_f_pointer(fft_buff_out,ptr,(/buff_size/))
		DFD%field_out(1:buff_len)=>ptr

		DFD%p_out=C_LOC(DFD%field_in(1))
		DFD%p_in=C_LOC(DFD%field_out(1))
		DFD%field_load_in(1:Nc,1:Nx,1:Ny_loc)=>DFD%field_out
		DFD%field_load_out(1:Nc,1:Nx,1:Ny_loc)=>DFD%field_in

		DFD%field_fft_x_in(1:Nx,1:Ny_loc,1:Nc)=>DFD%field_in
		DFD%field_fft_x_out(1:Nx,1:Ny_loc,1:Nc)=>DFD%field_out

		DFD%field_fft_y_in(1:K,1:Ny_loc,1:Np)=>DFD%field_out
		DFD%field_fft_y_out(1:K,1:Ny_loc,1:Np)=>DFD%field_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DFD%plans_time=0d0
		

		DFD%Nx=Nx
		DFD%Ny=Ny
		DFD%Nc=Nc
		DFD%Np=Np
		DFD%Ny_loc=Ny_loc
		DFD%K=K
		CALL CreateBlockPlans(DFD)
                CALL DROP_COUNTER(DFD)

		
		
	ENDSUBROUTINE
	SUBROUTINE PrintTimings(DFD)
		TYPE (DistributedFourierData),INTENT(IN)::DFD
		REAL(DOUBLEPARM)::fftx_av(2),ffty_av(2),x2y_av,y2x_av
		REAL(DOUBLEPARM)::kernel_av(2)
		REAL(DOUBLEPARM)::total_av(2)
		REAL(REALPARM)::av_times(16)
		REAL(REALPARM)::max_times(16)
		REAL(REALPARM)::min_times(16)
		REAL(REALPARM)::min_overhead(2)
		REAL(REALPARM)::max_overhead(2)
		REAL(REALPARM)::variation(8)
		REAL(REALPARM)::local_transpose_av(2,2)
		REAL(REALPARM)::max_total(2)
		REAL(REALPARM)::min_total(2)
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER::M,Nfwd,Nbwd
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A36, ES10.2E3, A3, F10.3 A3)"
		Nfwd=DFD%timer(FFT_FWD)%N
		Nbwd=DFD%timer(FFT_BWD)%N
		M=Nfwd+Nbwd
                IF (M==0) RETURN
		x2y_av=DFD%timer(FFT_FWD)%x2y_transpose/M
		y2x_av=DFD%timer(FFT_FWD)%y2x_transpose/M
		IF (Nfwd/=0) THEN 
			fftx_av(1)=DFD%timer(FFT_FWD)%fftx/Nfwd
			ffty_av(1)=DFD%timer(FFT_FWD)%ffty/Nfwd
			kernel_av(1)=DFD%timer(FFT_FWD)%kernel_total/Nfwd
			local_transpose_av(:,1)=DFD%timer(FFT_FWD)%local_transpose/Nfwd
			total_av(1)=kernel_av(1)+local_transpose_av(1,1)+local_transpose_av(2,1)
			
		ELSE
			fftx_av(1)=0d0
			ffty_av(1)=0d0
			kernel_av(1)=0d0
		ENDIF
		IF (Nbwd/=0) THEN 
			fftx_av(2)=DFD%timer(FFT_BWD)%fftx/Nbwd
			ffty_av(2)=DFD%timer(FFT_BWD)%ffty/Nbwd
			kernel_av(2)=DFD%timer(FFT_BWD)%kernel_total/Nbwd
			local_transpose_av(:,2)=DFD%timer(FFT_BWD)%local_transpose/Nbwd
			total_av(2)=kernel_av(2)+local_transpose_av(1,2)+local_transpose_av(2,2)
		ELSE
			fftx_av(2)=0d0
			ffty_av(2)=0d0
			kernel_av(2)=0d0
		ENDIF
		av_times(1)=x2y_av
		av_times(2)=y2x_av
		av_times(3:4)=fftx_av
		av_times(5:6)=ffty_av
		av_times(7:8)=kernel_av
		av_times(9:10)=total_av
		av_times(11:12)=local_transpose_av(:,1)
		av_times(13:14)=local_transpose_av(:,2)
		av_times(15:16)=kernel_av-av_times(1)-av_times(2)-av_times(3:4)-av_times(5:6)
		CALL MPI_REDUCE(av_times,min_times,16,MPI_DOUBLE,MPI_MIN,0,DFD%COMM,IERROR)
		CALL MPI_REDUCE(av_times,max_times,16,MPI_DOUBLE,MPI_MAX,0,DFD%COMM,IERROR)



		CALL LOGGER('===================================================================================')
		IF (Nfwd/=0) THEN 
			CALL PRINT_BORDER
			CALL PRINT_CALC_NUMBER("Number of forward Fourier transforms:", Nfwd)
			CALL PRINT_CALC_TIME_RANGE('Forward total time from:                 ', min_times(9)*Nfwd, max_times(9)*Nfwd)
			CALL PRINT_CALC_TIME_RANGE('Average total forward transform  from:  ', min_times(9), max_times(9))
			CALL PRINT_CALC_TIME_RANGE('Transfer  from:                          ', min_times(15),max_times(15))
			CALL PRINT_CALC_TIME_RANGE('Average forward transform along X from:  ',min_times(3),max_times(3))
			CALL PRINT_CALC_TIME_RANGE('Average forward transform along Y from:  ',min_times(5), max_times(5)) 

			CALL PRINT_CALC_TIME_RANGE('Average X to Y transpose from:           ',min_times(1), max_times(1))
			CALL PRINT_CALC_TIME_RANGE('Average Y to X transpose from:           ',min_times(2),max_times(2))
			CALL PRINT_CALC_TIME_RANGE('Initial transpose from:                  ',min_times(11),max_times(11))
			CALL PRINT_CALC_TIME_RANGE('Final transpose from:                    ',min_times(12),max_times(12))
		ENDIF
	
		IF (Nbwd/=0) THEN 
			CALL PRINT_BORDER
			CALL PRINT_CALC_NUMBER("Number of backward Fourier transforms:", Nbwd)
			CALL PRINT_CALC_TIME_RANGE('Backward total time from:               ', min_times(10)*Nbwd, max_times(10)*Nbwd)
			CALL PRINT_CALC_TIME_RANGE('Average total backward transform from: ', min_times(10), max_times(10))
			CALL PRINT_CALC_TIME_RANGE('Transfer  from:                         ', min_times(16),max_times(16))
			CALL PRINT_CALC_TIME_RANGE('Average forward transform along X from: ',min_times(4),max_times(4))
			CALL PRINT_CALC_TIME_RANGE('Average forward transform along Y from: ',min_times(6), max_times(6)) 

			CALL PRINT_CALC_TIME_RANGE('Average X to Y transpose from:          ',min_times(1), max_times(1))
			CALL PRINT_CALC_TIME_RANGE('Average Y to X transpose from:          ',min_times(2),max_times(2))
			CALL PRINT_CALC_TIME_RANGE('Initial transpose from:                 ',min_times(13),max_times(13))
			CALL PRINT_CALC_TIME_RANGE('Final transpose from:                   ',min_times(14),max_times(14))
		ENDIF
		CALL LOGGER('===================================================================================')
	END SUBROUTINE

        SUBROUTINE DROP_COUNTER(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		DFD%timer(FFT_FWD)%N=0
		DFD%timer(FFT_FWD)%fftx=0d0
		DFD%timer(FFT_FWD)%ffty=0d0
		DFD%timer(FFT_FWD)%x2y_transpose=0d0
		DFD%timer(FFT_FWD)%y2x_transpose=0d0
		DFD%timer(FFT_FWD)%kernel_total=0d0
		DFD%timer(FFT_FWD)%local_transpose(:)=0d0

		DFD%timer(FFT_BWD)%N=0
		DFD%timer(FFT_BWD)%fftx=0d0
		DFD%timer(FFT_BWD)%ffty=0d0
		DFD%timer(FFT_BWD)%x2y_transpose=0d0
		DFD%timer(FFT_BWD)%y2x_transpose=0d0
		
		DFD%timer(FFT_BWD)%kernel_total=0d0
		DFD%timer(FFT_BWD)%local_transpose(:)=0d0
        END SUBROUTINE
END MODULE
