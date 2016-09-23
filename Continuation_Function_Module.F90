!Copyright (c) 2016 Mikhail Kruglyakov 
!This file is part of GIEM2G.
!
!GIEM2G is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 2 of the License.
!
!GIEM2G is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with GFMRES.  If not, see <http://www.gnu.org/licenses/>.
!
!
!

MODULE CONTINUATION_FUNCTION_MODULE

	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE
	USE Timer_Module 
	USE DATA_TYPES_MODULE
	USE DISTRIBUTED_FFT_MODULE
	USE LOGGER_MODULE
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
		INTEGER(MPI_CTL_KIND):: matrix_comm
		INTEGER(MPI_CTL_KIND)::me
		INTEGER(MPI_CTL_KIND)::master_proc
		LOGICAL::master
		TYPE(TypeCounter)::counter
		TYPE (RECEIVER_TYPE),POINTER::recvs(:)
		INTEGER (C_INTPTR_T) ::localsize_in
		INTEGER (C_INTPTR_T) ::localsize_out
		TYPE(C_PTR)::p_in,p_out
		COMPLEX(REALPARM),POINTER::field_in4(:,:,:,:),field_out4(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_in3(:,:,:),field_out3(:,:,:)
		TYPE(C_PTR)::planFwd,planBwd
		COMPLEX(REALPARM),POINTER ::csiga(:,:,:)
		COMPLEX(REALPARM),POINTER ::csigb(:)
		LOGICAL::real_space
		LOGICAL::fftw_threads_ok
		TYPE(DistributedFourierData)::DFD_Current
		TYPE(DistributedFourierData)::DFD_Result
	  ENDTYPE
CONTAINS
	SUBROUTINE PrepareContinuationOperator(rc_op,anomaly,recvs,mcomm,fftw_threads_ok,DFD)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		TYPE (RECEIVER_TYPE),POINTER,INTENT(IN)::recvs(:)
		INTEGER(MPI_CTL_KIND),INTENT(IN)::mcomm 
		LOGICAL,OPTIONAL,INTENT(IN)::fftw_threads_ok
		TYPE (DistributedFourierData),OPTIONAL,INTENT(INOUT)::DFD
		INTEGER::Nx,Ny,Nz
		INTEGER::Nx2,Ny2,Nc,Nr3,Nopt
		INTEGER(MPI_CTL_KIND)::comm_size,IERROR,me
       		INTEGER(C_INTPTR_T)::tensor_size,buff_len,fft_len
		TYPE(C_PTR)::p1,p2
		INTEGER::NT
		Nx=anomaly%Nx
		Ny=anomaly%Ny
		Nz=anomaly%Nz
		rc_op%Nx=Nx
		rc_op%Ny=Ny
		rc_op%Nz=Nz
		rc_op%Nr=SIZE(recvs)
		rc_op%recvs=>recvs
		Nc=3*Nz
		Nr3=3*rc_op%Nr
		Nx2=2*Nx
		Ny2=2*Ny
		rc_op%Nx2Ny2=Nx2*Ny2
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
		IF (rc_op%fftw_threads_ok) THEN
			NT=OMP_GET_MAX_THREADS()
			CALL FFTW_PLAN_WITH_NTHREADS(NT)
		ENDIF
		IF (PRESENT(DFD)) THEN
			CALL PrepareDistributedFourierData(rc_op%DFD_Current,Nx2,Ny2,Nc,mcomm,&
                                &DFD%p_in,DFD%p_out)
		ELSE
                	fft_len=CalcLocalFFTSize(Nx2,Ny2,Nc,comm_size) 
			p1=ALLOCATE_BUFF(fft_len)
			p2=ALLOCATE_BUFF(fft_len)
			CALL PrepareDistributedFourierData(rc_op%DFD_Current,Nx2,Ny2,Nc,mcomm,&
                                &p1,p2,fft_len)
		ENDIF

                fft_len=CalcLocalFFTSize(Nx2,Ny2,Nr3,comm_size) 
		p1=ALLOCATE_BUFF(fft_len)
		p2=ALLOCATE_BUFF(fft_len)
		CALL PrepareDistributedFourierData(rc_op%DFD_Result,Nx2,Ny2,Nr3,mcomm,&
                               &p1,p2,fft_len)
		CALL CalcSizesForRC_OP(rc_op)
		CALL AllocateRC_OP(rc_op)

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
		rc_op%csiga=>NULL()
		rc_op%csigb=>NULL()
		anomaly%Ny_loc=rc_op%Ny_loc
	ENDSUBROUTINE
	SUBROUTINE FinalizeRCOperator(rc_op,bkg,anomaly,freq)
		TYPE(RC_Operator),INTENT(INOUT)::rc_op
		TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		REAL(REALPARM),INTENT(IN)::freq
		REAL(REALPARM)::w
		INTEGER::Iz,Ix,Iy
		IF (rc_op%real_space) THEN
			IF (ASSOCIATED(rc_op%csigb)) DEALLOCATE(rc_op%csigb)
			ALLOCATE(rc_op%csigb(rc_op%Nz))
			DO Iz=1,rc_op%Nz
					rc_op%csigb(Iz)=bkg%csigma(anomaly%Lnumber(Iz))
			ENDDO
			IF (ASSOCIATED(rc_op%csiga)) DEALLOCATE(rc_op%csiga)
			ALLOCATE(rc_op%csiga(rc_op%Nx,rc_op%Ny_loc,rc_op%Nz))
			w=freq*PI*2
			DO Iz=1,rc_op%Nz
			    DO Iy=1,rc_op%Ny_loc
				DO Ix=1,rc_op%Nx
#ifndef NO_DISPLACEMENT_CURRENTS
					rc_op%csiga(Ix,Iy,Iz)=anomaly%siga(Iz,Ix,Iy)-C_IONE*w*EPS0	
#else
					rc_op%csiga(Ix,Iy,Iz)=anomaly%siga(Iz,Ix,Iy)	
#endif
				   ENDDO
			ENDDO
		    ENDDO
		ENDIF
	END SUBROUTINE
#define no_compile	
	SUBROUTINE ReCalculation(rc_op,Eint,Ea,Ha)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Eint(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(INOUT)::Ea(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(INOUT)::Ha(:,:,:,:)
		INTEGER(MPI_CTL_KIND)::IERROR
		COMPLEX(REALPARM)::dsig
		INTEGER::Nx
		CHARACTER(LEN=*), PARAMETER  :: info_fmt = "(A, ES10.2E3)"
		REAL(8)::time1,time2
!		CALL MPI_BARRIER(rc_op%matrix_comm,IERROR)
		CALL PRINT_BORDER
		CALL LOGGER('Recalculation started')
		time1=GetTime()
		Nx=rc_op%Nx
		CALL LoadEint(rc_op,Eint)
		CALL APPLY_RC_E_OP(rc_op)
		IF (rc_op%real_space) THEN
			!$OMP PARALLEL	DEFAULT(SHARED)
			!$OMP WORKSHARE
				Ea=rc_op%field_out4(:,:,1:Nx,:)
			!$OMP END WORKSHARE
			!$OMP END PARALLEL
		ENDIF

		CALL LoadEint(rc_op,Eint)
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
		time2=GetTime()
		CALL LOGGER('Recalculation finished')
		CALL PRINT_CALC_TIME('Total tine: ',time2-time1)
		CALL PRINT_BORDER
				
	END SUBROUTINE
	SUBROUTINE LoadEint(rc_op,Eint)
		TYPE(RC_Operator),INTENT(INOUT)::rc_op
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Eint(:,:,:,:)
		INTEGER::Nz,Nx,Ny_loc
		INTEGER::Iz,Ix,Iy,Ic
		COMPLEX(REALPARM)::dsig
		Nx=rc_op%Nx
		Ny_loc=rc_op%Ny_loc
		Nz=rc_op%Nz
		IF (rc_op%real_space) THEN
			DO Iy=1,Ny_loc
				DO Ix=1,Nx
					DO Ic=1,3
						DO Iz=1,Nz
							dsig=(rc_op%csiga(Ix,Iy,Iz)-rc_op%csigb(Iz))
							rc_op%field_in4(Iz,Ic,Ix,Iy)=dsig*Eint(Ix,Iy,Iz,Ic)
						ENDDO
					ENDDO
				ENDDO
			ENDDO
			rc_op%field_in4(:,:,Nx+1:,:)=C_ZERO
		ELSE
			rc_op%field_in4=C_ZERO
		ENDIF
	ENDSUBROUTINE
	SUBROUTINE CalcSizesForRC_OP(rc_op) 
		TYPE(RC_Operator),INTENT(INOUT)::rc_op
		INTEGER(C_INTPTR_T)::tsize8(2)
		INTEGER(C_INTPTR_T)::Nz3,Nr3
		INTEGER(C_INTPTR_T)::block
		INTEGER(C_INTPTR_T):: CNy,CNy_offset !size and offset for electrical current in Y direction at this process 
		INTEGER(FFTW_COMM_SIZE)::COMM
		REAL(REALPARM)::size_g
		REAL(REALPARM)::tensor_size
		
		COMM=rc_op%matrix_comm
!----------------------------------------------------------------------------------------------------------!
		Nz3=rc_op%Nz*3
		rc_op%Ny_loc=rc_op%DFD_Current%Ny_loc
		rc_op%Ny_offset=rc_op%me*rc_op%Ny_loc

		rc_op%Nx2=rc_op%Nx*2
		rc_op%Nx2Ny_loc=rc_op%Nx2*rc_op%Ny_loc
		rc_op%Nx2Ny2=rc_op%Nx*rc_op%Ny*4

		Nr3=rc_op%Nr*3

		size_g=rc_op%Nx2Ny_loc*rc_op%Nz*rc_op%Nr*15.0*2*REALPARM/1024/1024/1024
		CALL PRINT_STORAGE_SIZE("RC matricies needs", size_g)
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

	FUNCTION ALLOCATE_BUFF(length) RESULT(buff_ptr)
		INTEGER(C_INTPTR_T)::length
		TYPE(C_PTR)::buff_ptr
		buff_ptr=fftw_alloc_complex(length)
	    END FUNCTION
	SUBROUTINE AllocateRC_OP(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		INTEGER::shape4(4),shape3(3)
		INTEGER::Nx2,Nz,Ny_loc,Nx2Ny_loc,Nr
		CALL AllocateRCMatrix(rc_op)
		Nz=rc_op%Nz
		Nx2=rc_op%Nx2
		Ny_loc=rc_op%Ny_loc
		Nx2Ny_loc=rc_op%Nx2Ny_loc
		Nr=rc_op%Nr
		
		rc_op%field_in4(1:Nz,1:3,1:Nx2,1:Ny_loc)=>rc_op%DFD_Current%field_out
		rc_op%field_in3(1:Nz,1:3,1:Nx2Ny_loc)=>rc_op%DFD_Current%field_out

		rc_op%field_out4(1:Nr,1:3,1:Nx2,1:Ny_loc)=>rc_op%DFD_Result%field_out
		rc_op%field_out3(1:Nr,1:3,1:Nx2Ny_loc)=>rc_op%DFD_Result%field_out
	ENDSUBROUTINE
	SUBROUTINE RC_OP_FFTW_FWD(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		CALL	CalcDistributedFourier(rc_op%DFD_Current,FFT_FWD)
	ENDSUBROUTINE
	SUBROUTINE RC_OP_FFTW_BWD(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		CALL	CalcDistributedFourier(rc_op%DFD_Result,FFT_BWD)
		rc_op%field_out4=rc_op%field_out4/rc_op%Nx2Ny2
	ENDSUBROUTINE
	SUBROUTINE CalcFFTofRCTensor(rc_op)
		TYPE(rc_operator),INTENT(INOUT)::rc_op
		INTEGER::Ir
		
		INTEGER(MPI_CTL_KIND)::IERROR
		REAL(8)::time1,time2
!!!! ----------------------ATTENTION DIRTY TRICK with output of CalcDistributedFourier(rc_op%DFD_*,FFT_FWD) -----------!!!!
		time1=GetTime()
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
		time2=GetTime()
                CALL PRINT_CALC_TIME('FFT of RC tensor has been computed: ',time2-time1)
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
		CALL DeleteMatrix(rc_op)
		CALL DeleteDistributedFourierData(rc_op%DFD_Current)
		CALL DeleteDistributedFourierData(rc_op%DFD_Result)
	ENDSUBROUTINE


	SUBROUTINE APPLY_RC_E_OP(rc_op)
		TYPE(RC_OPERATOR),INTENT(INOUT)::rc_op
		INTEGER::Ixy
		REAL(8)::time1,time2,time3,time4
		COMPLEX(REALPARM),POINTER::G(:,:,:,:)
		time1=GetTime()
		CALL RC_OP_FFTW_FWD(rc_op)
		time2=GetTime()
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
		time3=GetTime()
		CALL RC_OP_FFTW_BWD(rc_op)
		time4=GetTime()
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
		time1=GetTime()
		CALL RC_OP_FFTW_FWD(rc_op)
		time2=GetTime()
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
		time3=GetTime()
		CALL RC_OP_FFTW_BWD(rc_op)
		time4=GetTime()
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
