MODULE DISTRIBUTED_FFT_MODULE
!VERY VERY C-LIKE MODULE
!THERE ARE A LOT OF POINTERS FOR THE SAME MEMORY
!YES IT'S DIRTY AND UNCOMFORTABLE

!ATTENTION this  Fourier transform pracically is OUT PLACE and DESTROY INPUT

	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE
	USE Timer_Module 
	USE CHECK_MEMORY
	IMPLICIT NONE
	PRIVATE
	INTEGER, PARAMETER, PUBLIC :: FFT_FWD=1!FFTW_FORWARD
	INTEGER, PARAMETER, PUBLIC :: FFT_BWD=2!FFTW_BACKWARD

!	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_CTL=IOR(FFTW_PATIENT,FFTW_DESTROY_INPUT )
	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_CTL=IOR(FFTW_MEASURE,FFTW_DESTROY_INPUT )
!	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_CTL=IOR(FFTW_EXHAUSTIVE,FFTW_DESTROY_INPUT )
!	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_CTL=IOR(FFTW_ESTIMATE,FFTW_DESTROY_INPUT )

	LOGICAL, PARAMETER, PUBLIC :: NEW_DFD_ALLOC=.TRUE.
	LOGICAL, PARAMETER, PUBLIC :: OLD_DFD_ALLOC=.FALSE.
    TYPE PlanXY_TYPE
                TYPE(C_PTR)::planX,planY
    ENDTYPE

	TYPE FFTW_DATA_TYPE
		TYPE(C_PTR)::plan_fwd,plan_bwd
	ENDTYPE

	TYPE FFTW_TRANSPOSE_DATA
		TYPE(C_PTR)::plan
		REAL(REALPARM),POINTER::p_in(:),p_out(:)
	ENDTYPE

	TYPE BLOCK_TYPE
		INTEGER::Nm !(Nx,Ny_loc,Nm)
                INTEGER::Lxm
		INTEGER::size_x 
		INTEGER::Ix0
		INTEGER::len_x		

		INTEGER::K  !(K,Ny_loc,Np)
        INTEGER::Lkp
        INTEGER::size_y
		INTEGER::len_y	
		INTEGER::Iy0
		INTEGER::chunk_len !K*Ny_loc)
		
		INTEGER::sK,dK


		COMPLEX(REALPARM),POINTER:: field_fft_x_in(:,:,:) !(Nx,Ny_loc,Nm)
		COMPLEX(REALPARM),POINTER:: field_fft_x_out(:,:,:)!(Nx,Ny_loc,Nm)

		COMPLEX(REALPARM),POINTER:: field_fft_y_in(:,:,:) !(K,Ny_loc,Np)
		COMPLEX(REALPARM),POINTER:: field_fft_y_out(:,:,:) !(K,Ny_loc,Np)

		COMPLEX(REALPARM),POINTER:: field_in(:) 
		COMPLEX(REALPARM),POINTER:: field_out(:) 

		INTEGER(MPI_CTL_KIND)::comm
        TYPE(PlanXY_TYPE)::plan(FFT_FWD:FFT_BWD)
		REAL(DOUBLEPARM)::plans_time(2,FFT_FWD:FFT_BWD)
		TYPE(FFT_COUNTER),POINTER::timer(:)
	ENDTYPE

	TYPE FFT_COUNTER
		INTEGER::N
		REAL(DOUBLEPARM)::fftx,ffty, x2y_transpose, y2x_transpose
		REAL(DOUBLEPARM)::kernel_total 
	ENDTYPE

	TYPE DistributedFourierData
		TYPE(C_PTR)::p_in,p_out
		COMPLEX(REALPARM),POINTER::field_in(:),field_out(:)
		COMPLEX(REALPARM),POINTER::field_load_in(:,:,:),field_load_out(:,:,:)
		! ATTENTION field_load*(Nc,Nx,Ny) giem2g dimensions order
		TYPE(FFTW_DATA_TYPE),POINTER::FFTW_TRANSFORM
		TYPE(BLOCK_TYPE),POINTER::block(:)
		INTEGER::Nb
		INTEGER::Nc,Nx,Ny !Nx,Ny are sizes of 2D FFT transform
		INTEGER::Np !Np is number of processes 
		INTEGER::Ny_loc
		TYPE (FFT_COUNTER)::timer(FFT_FWD:FFT_BWD)
		REAL(DOUBLEPARM)::plans_time(2,FFT_FWD:FFT_BWD)
		INTEGER(MPI_CTL_KIND)::comm,me
        
		TYPE(FFTW_TRANSPOSE_DATA)::FFTW_TRANSPOSE

	ENDTYPE


	PUBLIC:: DistributedFourierData
	PUBLIC:: PrepareDistributedFourierData,DeleteDistributedFourierData
	PUBLIC::PrepareDFDLocal
	PUBLIC:: PrintTimings,TestBlocksNumber
	PUBLIC:: CalcDistributedFourier
	PUBLIC:: InitialTranspose, ProcessDistributedFourierKernel, FinalTranspose
	CONTAINS
	SUBROUTINE PrepareDistributedFourierData(DFD,Nx,Ny,Nc,comm,Nb,alloc_mem)
		TYPE (DistributedFourierData),TARGET,INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::Nx,Ny,Nc,Nb
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		LOGICAL,OPTIONAL,INTENT(IN)::alloc_mem
		INTEGER::Ny_loc
		INTEGER::M,Lblock,My
		INTEGER(MPI_CTL_KIND)::Np
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER(C_INTPTR_T)::Nbuff
		INTEGER::Nbuff32
		TYPE(BLOCK_TYPE),POINTER::block
		INTEGER::Ix,Iy,Ib,Ic
		INTEGER::Lxm,Lkp,Ix0,Iy0
!-----------------------------------------------------------------------------!
		CALL MPI_COMM_SIZE(comm,Np,IERROR) 
		CALL MPI_COMM_RANK(comm,DFD%me,IERROR) 
		M=MODULO(Ny,Np)
		IF (M/=0) THEN
			PRINT*, 'NOT IMPLEMENTED!! Wrong  Ny and number of processes '
			RETURN
		ENDIF
                IF (Nc/Nb*Nx<Np) THEN
                        PRINT*, 'VERY SMALL BLOCKS FOR ASYNC FOURIER TRANSFORM, GIEM2G FORCED HALTED'
                        CALL MPI_FINALIZE(IERROR)
                        STOP
                        RETURN
                ENDIF
		CALL PrepareDFDLocal(DFD,Nx,Ny,Nc,Np,Nb,alloc_mem)

		DO Ib=1,Nb
			block=>DFD%block(Ib)
			CALL MPI_COMM_DUP(comm,block%comm, IERROR)
		ENDDO
		IF (Nb==1) THEN
			CALL CreateAll2AllPlan(DFD)
		ENDIF
		CALL CHECK_MEM(dfd%me,0,comm)

		IF (DFD%me==0) THEN
			PRINT'(A)' ,'Block distributed FFT plan calculations at 0 process:'
			PRINT'(A, ES10.2E3 )' ,'Forward along X', DFD%plans_time(1,FFT_FWD)
			PRINT'(A, ES10.2E3 )' ,'Forward along Y', DFD%plans_time(2,FFT_FWD)
			PRINT'(A, ES10.2E3 )' ,'Backward along X', DFD%plans_time(1,FFT_BWD)
			PRINT'(A, ES10.2E3 )' ,'Backward along Y', DFD%plans_time(2,FFT_BWD)
		ENDIF
		DFD%comm=comm
	ENDSUBROUTINE


	SUBROUTINE CreateBlockPlans(block)
		TYPE(BLOCK_TYPE),POINTER,INTENT(INOUT)::block
		INTEGER::pNx(1),pNy(1)
		INTEGER::dx,dy,M
		INTEGER::x_shape(3),y_shape(3)
		REAL(DOUBLEPARM)::time1,time2
		x_shape=SHAPE(block%field_fft_x_in)
		pNx(1)=x_shape(1)
		dx=x_shape(1)
		M=x_shape(2)*x_shape(3)
		time1=GetTime()
		block%plan(FFT_FWD)%planX=fftw_plan_many_dft(ONE,pNx, M,block%field_fft_x_in,pNx,ONE,dx,&
								block%field_fft_x_out,pNx,ONE,dx,&
								&FFTW_FORWARD, FFT_CTL)
		time2=GetTime()
		block%plans_time(1,FFT_FWD)=time2-time1
		block%plan(FFT_BWD)%planX=fftw_plan_many_dft(ONE,pNx, M,block%field_fft_x_in,pNx,ONE,dx,&
								block%field_fft_x_out,pNx,ONE,dx,&
								&FFTW_BACKWARD, FFT_CTL)
		time1=GetTime()
		block%plans_time(1,FFT_BWD)=time1-time2


		y_shape=SHAPE(block%field_fft_y_in)
		pNy(1)=y_shape(2)*y_shape(3)
		M=y_shape(1)
		block%plan(FFT_FWD)%planY=fftw_plan_many_dft(ONE,pNy, M,block%field_fft_y_in,pNy,M,ONE,&
								block%field_fft_y_out,pNy,M,ONE,&
								&FFTW_FORWARD, FFT_CTL)

		time2=GetTime()
		block%plans_time(2,FFT_FWD)=time2-time1

		block%plan(FFT_BWD)%planY=fftw_plan_many_dft(ONE,pNy, M,block%field_fft_y_in,pNy,M,ONE,&
								block%field_fft_y_out,pNy,M,ONE,&
								&FFTW_BACKWARD, FFT_CTL)
		time1=GetTime()
		block%plans_time(2,FFT_BWD)=time1-time2
	ENDSUBROUTINE
	
	SUBROUTINE CreateAll2AllPlan(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		REAL(REALPARM),POINTER::r_in(:),r_out(:)
		INTEGER::M
		INTEGER(C_INTPTR_T)::N,Np
		INTEGER(C_INTPTR_T),PARAMETER::ONE=1
		N=2*DFD%block(1)%chunk_len
		Np=DFD%Np
		M=N*Np
		CALL c_f_pointer(DFD%p_out,r_in,(/M/))
		CALL c_f_pointer(DFD%p_in,r_out,(/M/))
		DFD%FFTW_TRANSPOSE%plan = fftw_mpi_plan_many_transpose(Np, Np, N, ONE, ONE, r_in, r_out, DFD%block(1)%comm,&
		& FFT_CTL);
		DFD%FFTW_TRANSPOSE%p_in=>r_in
		DFD%FFTW_TRANSPOSE%p_out=>r_out
	ENDSUBROUTINE

	SUBROUTINE CalcDistributedFourier(DFD,FFT_DIR)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
                CALL InitialTranspose(DFD)
	            CALL ProcessDistributedFourierKernel(DFD,FFT_DIR)
                CALL FinalTranspose(DFD)
                IF (FFT_DIR/=FFT_BWD) THEN
                        DFD%field_load_out=DFD%field_load_in
                ELSE
                        DFD%field_load_out=DFD%field_load_in/DFD%Nx/DFD%Ny
				ENDIF
	END SUBROUTINE

	
	SUBROUTINE ProcessDistributedFourierKernel(DFD,FFT_DIR)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
		INTEGER (MPI_CTL_KIND)::IERROR,comm
		INTEGER (MPI_CTL_KIND)::x2y_requests(DFD%Nb)
		INTEGER (MPI_CTL_KIND)::y2x_requests(DFD%Nb)
		INTEGER (MPI_CTL_KIND)::indices(DFD%Nb)
		INTEGER::Nrecv,Ib,Jb,chunk
		TYPE(BLOCK_TYPE),POINTER::block
		LOGICAL::recv1(1:DFD%Nb),recv2(1:DFD%Nb+1)
		COMPLEX(REALPARM),POINTER::p_send(:,:,:),p_recv(:,:,:)
		REAL(DOUBLEPARM)::time1,time2
#ifdef LEGACY_MPI
		CALL ProcessDistributedFourierKernelSync(DFD,FFT_DIR)
#else
		IF (DFD%Nb==1) THEN
			CALL ProcessDistributedFourierKernelSync(DFD,FFT_DIR)
			RETURN
		ENDIF
		time1=GetTime()
		recv1=.FALSE.
		recv2=.FALSE.
		recv2(DFD%Nb+1)=.TRUE.
		indices=1
		DFD%timer(FFT_DIR)%N=DFD%timer(FFT_DIR)%N+1
		DO Ib=1,DFD%Nb
			block=>DFD%block(Ib)
			CALL DistributedFourierX(block,FFT_DIR)
			CALL BlockTransposeXToY(block)
			p_send=>block%field_fft_y_out
			p_recv=>block%field_fft_y_in
			chunk=block%chunk_len
			comm=block%comm
			CALL MPI_IALLTOALL(p_send,chunk , MPI_DOUBLE_COMPLEX,&
						p_recv, chunk, MPI_DOUBLE_COMPLEX,&
						comm, x2y_requests(Ib), IERROR)

		ENDDO

		CALL	MPI_WAITSOME(DFD%Nb, x2y_requests, Nrecv,indices,MPI_STATUSES_IGNORE,IERROR)
		DO WHILE (Nrecv/=MPI_UNDEFINED)
			DO Jb=1,Nrecv 
				Ib=indices(Jb)
       			block=>DFD%block(Ib)
				CALL DistributedFourierY(block,FFT_DIR)
				p_send=>block%field_fft_y_out
				p_recv=>block%field_fft_y_in
				chunk=block%chunk_len
				comm=block%comm
				CALL MPI_IALLTOALL(p_send,chunk , MPI_DOUBLE_COMPLEX,&
						p_recv, chunk, MPI_DOUBLE_COMPLEX,&
						comm, y2x_requests(Ib), IERROR)
			ENDDO
			CALL	MPI_WAITSOME(DFD%Nb, x2y_requests, Nrecv,indices,MPI_STATUSES_IGNORE,IERROR)

		ENDDO
		CALL	MPI_WAITSOME(DFD%Nb, y2x_requests, Nrecv,indices,MPI_STATUSES_IGNORE,IERROR)
		DO WHILE (Nrecv/=MPI_UNDEFINED)
			DO Jb=1,Nrecv 
				Ib=indices(Jb)
				recv1(Ib)=.TRUE.
				recv2(Ib)=.TRUE.
			ENDDO
			DO Ib=1,DFD%Nb
        			block=>DFD%block(Ib)
				IF (recv1(Ib).AND.recv2(Ib+1)) THEN
					CALL BlockTransposeYToX(block)
					recv1(Ib)=.FALSE.
				ENDIF
			ENDDO
			CALL	MPI_WAITSOME(DFD%Nb, y2x_requests, Nrecv,indices,MPI_STATUSES_IGNORE,IERROR)
		ENDDO
		time2=GetTime()
		DFD%timer(FFT_DIR)%kernel_total=DFD%timer(FFT_DIR)%kernel_total+time2-time1
#endif
	END SUBROUTINE
 
!----------------------------------------------------------------------------------------------------------------------------------!

	SUBROUTINE ProcessDistributedFourierKernelSync(DFD,FFT_DIR)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::FFT_DIR
		INTEGER (MPI_CTL_KIND)::IERROR,comm
		INTEGER::Nrecv,Ib,Jb,chunk
		TYPE(BLOCK_TYPE),POINTER::block
		REAL(REALPARM),POINTER::p_send(:),p_recv(:)
		REAL(DOUBLEPARM)::time1,time2
		TYPE(C_PTR)::cp
		time1=GetTime()
		DFD%timer(FFT_DIR)%N=DFD%timer(FFT_DIR)%N+1
		block=>DFD%block(1)
		p_send=>DFD%FFTW_TRANSPOSE%p_in
		p_recv=>DFD%FFTW_TRANSPOSE%p_out
		CALL DistributedFourierX(block,FFT_DIR)
		CALL BlockTransposeXToY(block)
!		p_send=>block%field_fft_y_out
!		p_recv=>block%field_fft_y_in
		chunk=block%chunk_len
		comm=block%comm

!		CALL MPI_ALLTOALL(p_send,chunk , MPI_DOUBLE_COMPLEX,&
!					p_recv, chunk, MPI_DOUBLE_COMPLEX,&
!					comm,  IERROR)
!		CALL MPI_ALLTOALL(p_send,2*chunk, MPI_DOUBLE,&
!				p_recv, 2*chunk, MPI_DOUBLE,&
!				comm, IERROR)
		CALL fftw_mpi_execute_r2r(DFD%FFTW_TRANSPOSE%plan,p_send,p_recv)

		CALL DistributedFourierY(block,FFT_DIR)
!		p_send=>block%field_fft_y_out
!		p_recv=>block%field_fft_y_in
!		CALL MPI_ALLTOALL(p_send,chunk , MPI_DOUBLE_COMPLEX,&
!				p_recv, chunk, MPI_DOUBLE_COMPLEX,&
!				comm, IERROR)

!		CALL MPI_ALLTOALL(p_send,2*chunk, MPI_DOUBLE,&
!				p_recv, 2*chunk, MPI_DOUBLE,&
!				comm, IERROR)
		CALL fftw_mpi_execute_r2r(DFD%FFTW_TRANSPOSE%plan,p_send,p_recv)
		CALL BlockTransposeYToX(block)
		time2=GetTime()
		DFD%timer(FFT_DIR)%kernel_total=DFD%timer(FFT_DIR)%kernel_total+time2-time1
	END SUBROUTINE

	

	SUBROUTINE DistributedFourierY(block,FFT_DIR)
		TYPE(BLOCK_TYPE),POINTER,INTENT(INOUT)::block
		INTEGER,INTENT(IN)::FFT_DIR
		COMPLEX(REALPARM),POINTER::pin(:,:,:),pout(:,:,:)
		TYPE(C_PTR)::plan
		REAL(DOUBLEPARM)::time1,time2
		time1=GetTime()
		pin=>block%field_fft_y_in
		pout=>block%field_fft_y_out
		plan=block%plan(FFT_DIR)%planY
		CALL fftw_execute_dft(plan,pin,pout)
		time2=GetTime()
		block%timer(FFT_DIR)%ffty=block%timer(FFT_DIR)%ffty+time2-time1
!		pout=pin
	END SUBROUTINE

	SUBROUTINE DistributedFourierX(block,FFT_DIR)
		TYPE(BLOCK_TYPE),POINTER,INTENT(INOUT)::block
		INTEGER,INTENT(IN)::FFT_DIR
		COMPLEX(REALPARM),POINTER::pin(:,:,:),pout(:,:,:)
		TYPE(C_PTR)::plan
		REAL(DOUBLEPARM)::time1,time2
		time1=GetTime()
		pin=>block%field_fft_x_in
		pout=>block%field_fft_x_out
		plan=block%plan(FFT_DIR)%planX
		CALL fftw_execute_dft(plan,pin,pout)
		time2=GetTime()
		block%timer(FFT_DIR)%fftx=block%timer(FFT_DIR)%fftx+time2-time1
	END SUBROUTINE

        SUBROUTINE InitialTranspose(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Ix,Iy,Ic,l
		INTEGER::Nx,Ny,Nc
		COMPLEX(REALPARM),POINTER::p_in(:,:,:),p_out(:)
		p_in=>DFD%field_load_in
		p_out=>DFD%field_in
		Nx=DFD%Nx
		Ny=DFD%Ny_loc
		Nc=DFD%Nc
		!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ic,Iy,Ix,l)
		!$OMP DO SCHEDULE(GUIDED)
		DO Ic=1,Nc
			DO Iy=1,Ny
				DO Ix=1,Nx
					l=Ix+(Iy-1)*Nx+(Ic-1)*Ny*Nx
					p_out(l)=p_in(Ic,Ix,Iy)
				ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
        ENDSUBROUTINE

        SUBROUTINE FinalTranspose(DFD)
		TYPE (DistributedFourierData),INTENT(INOUT)::DFD
		INTEGER::Ix,Iy,Ic,l
		INTEGER::Nx,Ny,Nc
		COMPLEX(REALPARM),POINTER::p_in(:),p_out(:,:,:)
		p_in=>DFD%field_in
		p_out=>DFD%field_load_in
		Nx=DFD%Nx
		Ny=DFD%Ny_loc
		Nc=DFD%Nc
		!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ic,Iy,Ix,l)
		!$OMP DO SCHEDULE(GUIDED)! COLLAPSE(2)
		DO Iy=1,Ny
			DO Ix=1,Nx
				DO Ic=1,Nc
					l=Ix+(Iy-1)*Nx+(Ic-1)*Ny*Nx
					p_out(Ic,Ix,Iy)=p_in(l)
				ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
        ENDSUBROUTINE

	SUBROUTINE BlockTransposeXToY(block)
		TYPE(BLOCK_TYPE),POINTER,INTENT(INOUT)::block
		INTEGER::Ip,Iy,Ik,K
		INTEGER::l,Ix,Im,l0,Iy0,m0,Ny_loc
		INTEGER::Nx,Lkp,Np
		INTEGER::x_shape(3),y_shape(3)
		INTEGER::M
		COMPLEX(REALPARM),POINTER::p_in(:),p_out(:,:,:)
		REAL(DOUBLEPARM)::time1,time2
		time1=GetTime()
		x_shape=SHAPE(block%field_fft_x_in)
		y_shape=SHAPE(block%field_fft_y_in)
		K=block%K
		Lkp=block%Lkp	
		Ny_loc=x_shape(2)
		Nx=x_shape(1)
		Np=y_shape(3)
		p_out=>block%field_fft_y_out
		p_in=>block%field_out
		M=SIZE(p_in) !length of data
		!Linear writting, nonlinear reading
		!$OMP PARALLEL  DEFAULT(SHARED),PRIVATE(Ip,Iy,Ik,l0,m0,l)
		!$OMP DO SCHEDULE(GUIDED)	
		DO Ip=1,Np
			DO Iy=1,Ny_loc
				DO Ik=1,K
					l0=Ik+(Ip-1)*K+Lkp-1 !place in XZ
					m0=MODULO(l0-1,Nx)+1       !Ix=m0
					m0=l0-m0
					l=l0+(Iy-1)*Nx+m0*(Ny_loc-1)! Ix+(Iy-1)Nx+(Ic-1)NxNyLoc
					l=(l+M-ABS(M-l))/2 !if l>M we are in "addition" zone
					p_out(Ik,Iy,Ip)=p_in(l)
				ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
		time2=GetTime()
		block%timer(FFT_FWD)%x2y_transpose=time2-time1+block%timer(FFT_FWD)%x2y_transpose
	END SUBROUTINE

	SUBROUTINE BlockTransposeYToX(block)
		TYPE(BLOCK_TYPE),POINTER,INTENT(INOUT)::block
		INTEGER::Im,Iy,Ix
		INTEGER::l0,K,I,Il0,m0,l
		INTEGER::Ny_loc,Nx,Np,M,Nm
		INTEGER::x_shape(3),y_shape(3)
		INTEGER::sK,dK,sI,dI,sL,dL,I2	
		COMPLEX(REALPARM),POINTER::p_in(:),p_out(:,:,:)
		REAL(DOUBLEPARM)::time1,time2
		time1=GetTime()
		x_shape=SHAPE(block%field_fft_x_in)
		y_shape=SHAPE(block%field_fft_y_in)

		Nm=block%Nm
		Ny_loc=x_shape(2)
		Nx=x_shape(1)
		Np=y_shape(3)
		

		sK=block%sK
		dK=block%dK


		I2=block%Lkp+block%size_y
		sI=block%Lkp+I2
		dI=block%size_y
        
		sL=2*block%Iy0+block%len_y
		dL=block%len_y
		M=block%Lxm				

		p_out=>block%field_fft_x_in
		p_in=>block%field_out
		!Linear writting, nonlinear reading
		!$OMP PARALLEL  DEFAULT(SHARED),PRIVATE(Im,Iy,Ix,l0,K,I,Il0,m0,l)
		!$OMP DO SCHEDULE(GUIDED)	
		DO Im=1,Nm
			DO Iy=1,Ny_loc
				DO Ix=1,Nx
					l0=Ix+(Im-1)*Nx+M-1 !place in XZ
					K=(sK+dK*SIGN(ONE,l0-I2))/2
					I=(sI+SIGN(dI,l0-I2))/2
					Il0=(sL+SIGN(dL,l0-I2))/2-1
					l0=l0-I+1
					m0=MODULO(l0-1,K)+1          !Ik=m0
					l=Il0+l0+(Iy-1)*K+(l0-m0)*(Ny_loc-1)! Ix+(Iy-1)Nx+(Ic-1)NxNyLoc
					p_out(Ix,Iy,Im)=p_in(l)
				ENDDO
			ENDDO
		ENDDO
		!$OMP ENDDO
		!$OMP END PARALLEL
		time2=GetTime()
		block%timer(FFT_FWD)%y2x_transpose=time2-time1+block%timer(FFT_FWD)%y2x_transpose
	END SUBROUTINE
	SUBROUTINE PrintTimings(DFD,kernel_max)
		TYPE (DistributedFourierData),INTENT(IN)::DFD
		REAL(REALPARM),OPTIONAL,INTENT(OUT)::kernel_max
		REAL(DOUBLEPARM)::fftx_av(2),ffty_av(2),x2y_av,y2x_av
		REAL(DOUBLEPARM)::kernel_av(2)
		REAL(REALPARM)::av_times(8)
		REAL(REALPARM)::max_times(8)
		REAL(REALPARM)::min_times(8)
		REAL(REALPARM)::min_overhead(2)
		REAL(REALPARM)::max_overhead(2)
		REAL(REALPARM)::variation(8)
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER::M,Nfwd,Nbwd
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A36, ES10.2E3, A3, F10.3 A3)"
		Nfwd=DFD%timer(FFT_FWD)%N
		Nbwd=DFD%timer(FFT_BWD)%N
		M=Nfwd+Nbwd
		IF (M/=0) THEN
			x2y_av=DFD%timer(FFT_FWD)%x2y_transpose/M
			y2x_av=DFD%timer(FFT_FWD)%y2x_transpose/M
		ELSE
		IF (DFD%me==0) PRINT*, 'There were NO transforms with this DFD'
			RETURN
		ENDIF
		
		IF (Nfwd/=0) THEN 
			fftx_av(1)=DFD%timer(FFT_FWD)%fftx/Nfwd
			ffty_av(1)=DFD%timer(FFT_FWD)%ffty/Nfwd
			kernel_av(1)=DFD%timer(FFT_FWD)%kernel_total/Nfwd
		ELSE
			fftx_av(1)=0d0
			ffty_av(1)=0d0
			kernel_av(1)=0d0
		ENDIF
		IF (Nbwd/=0) THEN 
			fftx_av(2)=DFD%timer(FFT_BWD)%fftx/Nbwd
			ffty_av(2)=DFD%timer(FFT_BWD)%ffty/Nbwd
			kernel_av(2)=DFD%timer(FFT_BWD)%kernel_total/Nbwd
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
		CALL MPI_REDUCE(av_times,min_times,8,MPI_DOUBLE,MPI_MIN,0,dfd%comm,IERROR)
		CALL MPI_REDUCE(av_times,max_times,8,MPI_DOUBLE,MPI_MAX,0,dfd%comm,IERROR)
		variation=(max_times/min_times-1)*100
		min_overhead(1)=(min_times(7)-min_times(1)-min_times(2)-min_times(3)-min_times(5))/min_times(7)
		max_overhead(1)=(max_times(7)-max_times(1)-max_times(2)-max_times(3)-max_times(5))/max_times(7)

		min_overhead(2)=(min_times(8)-min_times(1)-min_times(2)-min_times(4)-min_times(6))/min_times(8)
		max_overhead(2)=(max_times(8)-max_times(1)-max_times(2)-max_times(4)-max_times(6))/max_times(8)
		IF (DFD%me==0) THEN
			PRINT'(A80)','==================================================================================='
			PRINT'(A, I6, A)',"Block distributed FFT with ", DFD%Nb, " blocks"


			IF (Nfwd/=0) THEN 
				PRINT'(A80)','***********************************************************************************'
				PRINT'(A, I6)',"Number of forward Fourier transforms:", Nfwd
				PRINT info_fmt,"Forward total",max_times(7)*Nfwd," s ",variation(7), " % "
				PRINT info_fmt,"Average kernel forward transform   ",max_times(7), " s ",variation(7), " % "
				PRINT '(A36, F9.3, A3)',"Approximate overhead for forward ",max_overhead(1)*100, " % "

				PRINT info_fmt,"Average forward transform along X  ",max_times(3), " s ",variation(3), " % "
				PRINT info_fmt,"Average forward transform along Y  ",max_times(5), " s ",variation(5), " % "

				PRINT info_fmt,"Average X to Y transpose",max_times(1), " s ",variation(1), " % "
				PRINT info_fmt,"Average Y to X transpose",min_times(2), " s ",variation(2), " % "
			ENDIF
 
			IF (Nbwd/=0) THEN 
				PRINT'(A80)','***********************************************************************************'
				PRINT'(A, I6)',"Number of backward Fourier transforms:", Nbwd

				PRINT info_fmt,"Backward total",max_times(8)*Nbwd," s ",variation(8), " % "
				PRINT info_fmt,"Average kernel backward transform  ",min_times(8), " s ",variation(8), " % "
				PRINT '(A36, F9.3, A3)',"Approximate overhead for backward ",max_overhead(2)*100, " % "

				PRINT info_fmt,"Average backward transform along X ",min_times(4), " s ",variation(4), " % " 
				PRINT info_fmt,"Average backward transform along Y ",min_times(6), " s ",variation(6), " % "
 
				PRINT info_fmt,"Average X to Y transpose",max_times(1), " s ",variation(1), " % "
				PRINT info_fmt,"Average Y to X transpose",min_times(2), " s ",variation(2), " % "
			ENDIF
			PRINT'(A80)','==================================================================================='
		ENDIF
		IF (PRESENT(kernel_max)) THEN
			kernel_max=MAX(max_times(7),max_times(8))
			CALL MPI_BCAST(kernel_max,1,MPI_DOUBLE,0,DFD%comm, IERROR)
		ENDIF
	END SUBROUTINE
	SUBROUTINE DeleteDistributedFourierData(DFD)
		TYPE (DistributedFourierData),TARGET,INTENT(INOUT)::DFD
		INTEGER::Ib
		INTEGER(MPI_CTL_KIND)::IERROR
		IF ((DFD%Nb/=0).AND.ASSOCIATED(DFD%Block)) THEN
			DO Ib=1,DFD%Nb
				CALL	fftw_destroy_plan(DFD%block(Ib)%plan(FFT_FWD)%planX);	
				CALL	fftw_destroy_plan(DFD%block(Ib)%plan(FFT_FWD)%planY);	
				CALL	fftw_destroy_plan(DFD%block(Ib)%plan(FFT_BWD)%planX);	
				CALL	fftw_destroy_plan(DFD%block(Ib)%plan(FFT_BWD)%planY);
				CALL MPI_COMM_FREE(DFD%block(Ib)%comm,IERROR)	
			ENDDO
			DEALLOCATE(DFD%block)
		ELSEIF (ASSOCIATED(DFD%FFTW_TRANSFORM)) THEN
			CALL	fftw_destroy_plan(DFD%FFTW_TRANSFORM%plan_Fwd);	
			CALL	fftw_destroy_plan(DFD%FFTW_TRANSFORM%plan_Bwd);	
		ENDIF
		CALL fftw_free(DFD%p_in)
		CALL fftw_free(DFD%p_out)
		DFD%block=>NULL()
		DFD%field_in=>NULL()
		DFD%field_out=>NULL()
		DFD%field_load_in=>NULL()
		DFD%field_load_out=>NULL()
		DFD%FFTW_TRANSFORM=>NULL()
	ENDSUBROUTINE
	
	SUBROUTINE TestBlocksNumber(DFD,Nx,Ny,Nc,comm,Nbmax,Nopt)
		TYPE (DistributedFourierData),TARGET,INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::Nx,Ny,Nc,Nbmax
		INTEGER,OPTIONAL,INTENT(OUT)::Nopt
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER::Nb,Ic,Iopt
		REAL(DOUBLEPARM)::kernel_time,kernel_min
		REAL(DOUBLEPARM)::time1,time2
		Iopt=1
		kernel_min=HUGE(kernel_time)	
		kernel_time=0d0	
		DO Nb=1,Nbmax
			CALL PrepareDistributedFourierData(DFD,Nx,Ny,Nc,comm,Nb,NEW_DFD_ALLOC)
			DFD%field_out=1d0
			DFD%field_load_in=(1d0,-246d0)
			DO Ic=1,Nc/3
				CALL ProcessDistributedFourierKernel(DFD,FFT_FWD)
				CALL ProcessDistributedFourierKernel(DFD,FFT_BWD)
			ENDDO
			CALL PrintTimings(DFD,kernel_time)
			IF (kernel_time<kernel_min) THEN
				kernel_min=kernel_time
				Iopt=Nb
			ENDIF
			CALL DeleteDistributedFourierData(DFD)
		ENDDO
		IF (PRESENT(Nopt)) THEN
			Nopt=Iopt
			CALL MPI_BCAST(Nopt,1,MPI_INTEGER,0,comm, IERROR)
		ENDIF
	END SUBROUTINE
	SUBROUTINE PrepareDFDLocal(DFD,Nx,Ny,Nc,Np,Nb,alloc_mem)
		TYPE (DistributedFourierData),TARGET,INTENT(INOUT)::DFD
		INTEGER,INTENT(IN)::Nx,Ny,Nc,Nb,Np
		LOGICAL,OPTIONAL,INTENT(IN)::alloc_mem
		INTEGER::Ny_loc
		INTEGER::M,Lblock,My
		INTEGER(C_INTPTR_T)::Nbuff
		INTEGER::Nbuff32
		TYPE(BLOCK_TYPE),POINTER::block
		INTEGER::Ix,Iy,Ib,Ic
		INTEGER::Lxm,Lkp,Ix0,Iy0
		ALLOCATE(DFD%block(1:Nb))
		DFD%FFTW_TRANSFORM=>NULL()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		M=MODULO(Nc,Nb)
		IF (M==0) THEN
			Lblock=Nc/Nb
			DO Ib=1,Nb
				block=>DFD%block(Ib)
				block%Nm=Lblock
			ENDDO
		ELSE
			Lblock=(Nc-M)/Nb
			DO Ib=1,Nb-1
				DFD%block(Ib)%Nm=Lblock
			ENDDO
			DFD%block(Nb)%Nm=(Lblock+M)
		ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Ny_loc=Ny/Np
		Lxm=1
		Ix0=1
		DO Ib=1,Nb
			block=>DFD%block(Ib)
			block%size_x=Nx*block%Nm
			block%Lxm=Lxm
			block%Ix0=Ix0
			block%len_x=block%size_x*Ny_loc
			Lxm=Lxm+block%size_x
			Ix0=Ix0+block%len_x
		ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		My=0
		Iy0=1
		Lkp=1
		DO Ib=1,Nb
			block=>DFD%block(Ib)
			My=My+block%size_x
			M=MODULO(My,Np)
			
			IF (M==0) THEN
				block%K=My/Np
			ELSEIF(Ib/=Nb) THEN
				block%K=(My-M)/Np
			ELSE
				block%K=(My-M)/Np+1
			ENDIF
			My=M
			block%size_y=block%K*Np
			block%chunk_len=block%K*Ny_loc
			block%Iy0=Iy0
			block%Lkp=Lkp
			block%len_y=block%size_y*Ny_loc

			Lkp=Lkp+block%size_y
			Iy0=Iy0+block%len_y
		ENDDO
		Nbuff=Iy0-1

		IF (PRESENT(alloc_mem)) THEN
			IF (alloc_mem.EQV.NEW_DFD_ALLOC) THEN
				DFD%p_in=fftw_alloc_complex(Nbuff) 
				DFD%p_out=fftw_alloc_complex(Nbuff)
			ENDIF
		ELSE 
			DFD%p_in=fftw_alloc_complex(Nbuff) 
			DFD%p_out=fftw_alloc_complex(Nbuff)
		ENDIF

		Nbuff32=Nbuff
		CALL c_f_pointer(DFD%p_out,DFD%field_in,(/Nbuff32/))
		CALL c_f_pointer(DFD%p_in,DFD%field_out,(/Nbuff32/))
		DFD%field_load_in(1:Nc,1:Nx,1:Ny_loc)=>DFD%field_out
		DFD%field_load_out(1:Nc,1:Nx,1:Ny_loc)=>DFD%field_in

		DFD%plans_time=0d0
		DO Ib=1,Nb
			block=>DFD%block(Ib)

			block%field_in=>DFD%field_in
			block%field_out=>DFD%field_out

			block%field_fft_x_in(1:Nx,1:Ny_loc,1:block%Nm)=>DFD%field_in(block%Ix0:)

			block%field_fft_x_out(1:Nx,1:Ny_loc,1:block%Nm)=>DFD%field_out(block%Ix0:)


			block%field_fft_y_in(1:block%K,1:Ny_loc,1:Np)=>DFD%field_out(block%Iy0:)

			block%field_fft_y_out(1:block%K,1:Ny_loc,1:Np)=>DFD%field_in(block%Iy0:)

			block%timer=>DFD%timer
			block%plans_time=0d0

			CALL CreateBlockPlans(block)
			DFD%plans_time=DFD%plans_time+block%plans_time
			IF (Ib/=Nb) THEN
				block%sK=block%K+DFD%block(Ib+1)%K
				block%dK=DFD%block(Ib+1)%K-block%K
			ELSE
				block%sK=2*block%K
				block%dK=0
			ENDIF
		ENDDO
		

		DFD%Nx=Nx
		DFD%Ny=Ny
		DFD%Nc=Nc
		DFD%Np=Np
		DFD%Ny_loc=Ny_loc
		DFD%Nb=Nb
		DFD%timer(FFT_FWD)%N=0
		DFD%timer(FFT_FWD)%fftx=0d0
		DFD%timer(FFT_FWD)%ffty=0d0
		DFD%timer(FFT_FWD)%x2y_transpose=0d0
		DFD%timer(FFT_FWD)%y2x_transpose=0d0
		DFD%timer(FFT_FWD)%kernel_total=0d0

		DFD%timer(FFT_BWD)%N=0
		DFD%timer(FFT_BWD)%fftx=0d0
		DFD%timer(FFT_BWD)%ffty=0d0
		DFD%timer(FFT_BWD)%x2y_transpose=0d0
		DFD%timer(FFT_BWD)%y2x_transpose=0d0
		
		DFD%timer(FFT_BWD)%kernel_total=0d0
!-----------------------------------------------------------------------------!

		
	ENDSUBROUTINE
END MODULE
