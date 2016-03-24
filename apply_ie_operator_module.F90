MODULE APPLY_IE_OPERATOR_MODULE 
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE
	USE DATA_TYPES_MODULE
	USE INTEGRAL_EQUATION_MODULE
	USE Timer_Module 
        USE LOGGER_MODULE
	IMPLICIT NONE
	PRIVATE

!For integral equation kernel!

	INTERFACE APPLY_IE_OP
		MODULE PROCEDURE APPLY_PRECOND, APPLY_IE_ZEROS
	END INTERFACE
        PUBLIC:: APPLY_IE_OP,APPLY_PRECOND, APPLY_IE_ZEROS
	CONTAINS


	SUBROUTINE APPLY_PRECOND(ie_op,v_in,v_out)
		TYPE(IntegralEquationOperator),TARGET,INTENT(INOUT)::ie_op
		COMPLEX(REALPARM),POINTER,INTENT(IN)::v_in(:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::v_out(:)
		COMPLEX(REALPARM),POINTER::field_in(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_out(:,:,:,:)
		COMPLEX(REALPARM),POINTER::fft_buff(:,:,:,:)

		COMPLEX(REALPARM)::d1,d2
		COMPLEX(REALPARM)::asiga,dsig,gsig
		REAL(8)::time1,time2
		INTEGER::IERROR
		INTEGER ::Ix,Iy,Iz,Ic,Ixy,N
		time1=GetTime()
		field_in(1:ie_op%Nx,1:ie_op%Ny_loc,1:ie_op%Nz,1:3)=>v_in
		field_out(1:ie_op%Nx,1:ie_op%Ny_loc,1:ie_op%Nz,1:3)=>v_out
                fft_buff(1:2*ie_op%Nx,1:ie_op%Ny_loc,1:ie_op%Nz,1:3)=>ie_op%DFD%field_in
                N=4*ie_op%Nx*ie_op%Ny
		!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(Ix,Ic,Iz,Iy,asiga,dsig,gsig)
#ifndef IBM_Bluegene
                !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
#else
                !$OMP DO SCHEDULE(GUIDED) 
#endif
                DO Iz=1,ie_op%Nz
                        DO Iy=1,ie_op%Ny_loc
                                DO Ix=1,ie_op%Nx
                                        asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
                                        dsig=ie_op%csiga(Ix,Iy,Iz)-ie_op%csigb(Iz)
                                        gsig=dsig*asiga

                                        fft_buff(Ix,Iy,Iz,EX)=&
                                        &field_in(Ix,Iy,Iz,EX)*gsig
                                ENDDO
                        ENDDO
                ENDDO
                !$OMP END DO
#ifndef IBM_Bluegene
                !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
#else
                !$OMP DO SCHEDULE(GUIDED) 
#endif
                DO Iz=1,ie_op%Nz
                        DO Iy=1,ie_op%Ny_loc
                                DO Ix=1,ie_op%Nx
                                        asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
                                        dsig=ie_op%csiga(Ix,Iy,Iz)-ie_op%csigb(Iz)
                                        gsig=dsig*asiga

                                        fft_buff(Ix,Iy,Iz,EY)=&
                                        &field_in(Ix,Iy,Iz,EY)*gsig
                                ENDDO
                        ENDDO
                ENDDO
                !$OMP END DO
#ifndef IBM_Bluegene
                !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
#else
                !$OMP DO SCHEDULE(GUIDED) 
#endif
                DO Iz=1,ie_op%Nz
                        DO Iy=1,ie_op%Ny_loc
                                DO Ix=1,ie_op%Nx
                                        asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
                                        dsig=ie_op%csiga(Ix,Iy,Iz)-ie_op%csigb(Iz)
                                        gsig=dsig*asiga

                                        fft_buff(Ix,Iy,Iz,EZ)=&
                                        &field_in(Ix,Iy,Iz,EZ)*gsig
                                ENDDO
                        ENDDO
                ENDDO
                !$OMP END DO
                !$OMP WORKSHARE
                fft_buff(ie_op%Nx+1:,:,:,:)=C_ZERO
                !$OMP END WORKSHARE
		!$OMP END PARALLEL

		CALL	MULT_IE_OP(ie_op)
		!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Iy,Ix,Iz,d1,d2,asiga)
#ifndef IBM_Bluegene
                !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
#else
                !$OMP DO SCHEDULE(GUIDED) 
#endif
		DO Iz=1,ie_op%Nz
			DO Iy=1,ie_op%Ny_loc
				DO Ix=1,ie_op%Nx
					asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
					d1=field_in(Ix,Iy,Iz,EX)*asiga
					d2=d1-fft_buff(Ix,Iy,Iz,EX)/(ie_op%dz(Iz))/N
					field_out(Ix,Iy,Iz,EX)=d2*ie_op%sqsigb(Iz)
				ENDDO
			ENDDO
		ENDDO
		!$OMP END DO
#ifndef IBM_Bluegene
                !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
#else
                !$OMP DO SCHEDULE(GUIDED) 
#endif
		DO Iz=1,ie_op%Nz
			DO Iy=1,ie_op%Ny_loc
				DO Ix=1,ie_op%Nx
					asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
					d1=field_in(Ix,Iy,Iz,EY)*asiga
					d2=d1-fft_buff(Ix,Iy,Iz,EY)/(ie_op%dz(Iz))/N
					field_out(Ix,Iy,Iz,EY)=d2*ie_op%sqsigb(Iz)
				ENDDO
			ENDDO
		ENDDO
		!$OMP END DO
#ifndef IBM_Bluegene
                !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
#else
                !$OMP DO SCHEDULE(GUIDED) 
#endif
		DO Iz=1,ie_op%Nz
			DO Iy=1,ie_op%Ny_loc
				DO Ix=1,ie_op%Nx
					asiga=C_TWO*ie_op%sqsigb(Iz)/(ie_op%csiga(Ix,Iy,Iz)+CONJG(ie_op%csigb(Iz)))
					d1=field_in(Ix,Iy,Iz,EZ)*asiga
					d2=d1-fft_buff(Ix,Iy,Iz,EZ)/(ie_op%dz(Iz))/N
					field_out(Ix,Iy,Iz,EZ)=d2*ie_op%sqsigb(Iz)
				ENDDO
			ENDDO
		ENDDO
		!$OMP END DO
		!$OMP END PARALLEL
		time2=GetTime()
		ie_op%counter%mult_num=ie_op%counter%mult_num+1
		ie_op%counter%apply=ie_op%counter%apply+time2-time1
	END SUBROUTINE
	SUBROUTINE APPLY_IE_ZEROS(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op

		!$OMP PARALLEL	 DEFAULT(SHARED)
		!$OMP WORKSHARE
                ie_op%DFD%field_in=C_ZERO
		!$OMP END WORKSHARE
		!$OMP ENDPARALLEL
		CALL MULT_IE_OP(ie_op)
	ENDSUBROUTINE
	SUBROUTINE MULT_IE_OP(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		INTEGER::Ixy,Nz
		REAL(8)::time1,time2,time3,time4
		COMPLEX(REALPARM),POINTER::p_send(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER::p_recv(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER::p_in(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER::p_out(:,:,:,:,:)
                INTEGER(MPI_CTL_KIND)::NxHalf,xfirst,xlast,data_length
                INTEGER(MPI_CTL_KIND)::proc,recv_start
		INTEGER(MPI_CTL_KIND)::requests(2),IERROR
		COMPLEX(REALPARM)::MirrorY
		time1=GetTime()
		CALL IE_OP_FFTW_FWD_PREPARED(ie_op)
		time2=GetTime()
		ie_op%counter%mult_fftw=ie_op%counter%mult_fftw+time2-time1

                NxHalf=ie_op%Nx/2
                data_length=3*ie_op%Nz*ie_op%Ny_loc*ie_op%Nx

                IF (ie_op%real_space) THEN
                        p_send=>ie_op%field_inT(:,:,:,:,NxHalf:)
                        p_recv=>ie_op%field_outT(:,:,:,:,NxHalf:)
                        MirrorY=C_ONE
                ELSE
                        p_send=>ie_op%field_inT(:,:,:,:,0:)
                        p_recv=>ie_op%field_outT(:,:,:,:,0:)
                        MirrorY=-C_ONE
                ENDIF

                proc=ie_op%partner

                CALL  MPI_IRECV(p_recv,data_length,MPI_DOUBLE_COMPLEX,proc,proc,&
                         &ie_op%ie_comm,requests(1),IERROR)
                CALL  MPI_ISEND(p_send,data_length,MPI_DOUBLE_COMPLEX,proc,ie_op%me,&
                         &ie_op%ie_comm,requests(2),IERROR)
                p_in=>ie_op%field_inT
                p_out=>ie_op%field_outT


                IF (ie_op%real_space) THEN
                        CALL VERTICAL_MULT_GENERAL(ie_op,1,NxHalf-1,MirrorY,p_in,p_out)
                        CALL VERTICAL_MULT_GENERAL_ZERO(ie_op,MirrorY,p_in,p_out)
                ELSE
                        CALL VERTICAL_MULT_GENERAL(ie_op,NxHalf,ie_op%Nx-1,MirrorY,p_in,p_out)
                ENDIF

                CALL MPI_WAITALL(2, requests,  MPI_STATUSES_IGNORE, IERROR)

                MirrorY=-MirrorY
                IF (ie_op%real_space) THEN
                       p_in(1:,1:,1:,1:,0:)=>p_recv 
                       p_out(1:,1:,1:,1:,0:)=>p_send 
                        CALL VERTICAL_MULT_GENERAL(ie_op,1,NxHalf-1,MirrorY,p_in,p_out)
                        CALL VERTICAL_MULT_GENERAL_ZERO(ie_op,MirrorY,p_in,p_out)
                ELSE
                       p_in(1:,1:,1:,1:,NxHalf:)=>p_recv 
                       p_out(1:,1:,1:,1:,NxHalf:)=>p_send 
                        CALL VERTICAL_MULT_GENERAL(ie_op,NxHalf,ie_op%Nx-1,MirrorY,p_in,p_out)
                ENDIF
                CALL  MPI_IRECV(p_recv,data_length,MPI_DOUBLE_COMPLEX,proc,proc,&
                         &ie_op%ie_comm,requests(1),IERROR)
                CALL  MPI_ISEND(p_send,data_length,MPI_DOUBLE_COMPLEX,proc,ie_op%me,&
                         &ie_op%ie_comm,requests(2),IERROR)

                CALL MPI_WAITALL(2, requests,  MPI_STATUSES_IGNORE, IERROR)

		time3=GetTime()
		CALL IE_OP_FFTW_BWD_PREPARED(ie_op)
		time4=GetTime()
		ie_op%counter%mult_fftw_b=ie_op%counter%mult_fftw_b+time4-time3
		ie_op%counter%mult_zgemv=ie_op%counter%mult_zgemv+time3-time2

	ENDSUBROUTINE
!----------------------------------------------------------------------------------------------------------------!
        SUBROUTINE VERTICAL_MULT_GENERAL(ie_op,xfirst,xlast,MirrorY,p_in,p_out)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		COMPLEX(REALPARM),POINTER,INTENT(IN)::p_in(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::p_out(:,:,:,:,:)
                COMPLEX(REALPARM),INTENT(IN)::MirrorY
                INTEGER,INTENT(IN)::xfirst,xlast

                INTEGER::Ix,Iy

		!$OMP PARALLEL	PRIVATE(Ix,Iy) DEFAULT(SHARED)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ix=xfirst,xlast
                        DO Iy=1,ie_op%Ny_loc
                                CALL IE_OP_ZGEMV_SYMM_DOUBLE(ie_op,S_EXX,EX,EX,Ix,Iy,C_ONE,C_ONE,C_ZERO,p_in,p_out)

                                CALL IE_OP_ZGEMV_SYMM_DOUBLE(ie_op,S_EXY,EY,EX,Ix,Iy,MirrorY,-MirrorY,C_ONE,p_in,p_out)
                                CALL IE_OP_ZGEMV_SYMM_DOUBLE(ie_op,S_EYX,EX,EY,Ix,Iy,MirrorY,-MirrorY,C_ZERO,p_in,p_out)

                                CALL IE_OP_ZGEMV_SYMM_DOUBLE(ie_op,S_EYY,EY,EY,Ix,Iy,C_ONE,C_ONE,C_ONE,p_in,p_out)

                                CALL IE_OP_ZGEMV_SYMM_DOUBLE(ie_op,S_EZZ,EZ,EZ,Ix,Iy,C_ONE,C_ONE,C_ZERO,p_in,p_out)

                                CALL IE_OP_ZGEMV_ASYM_DOUBLE(ie_op,'N',A_EXZ,EZ,EX,Ix,Iy,C_ONE,-C_ONE,C_ONE,p_in,p_out)
                                CALL IE_OP_ZGEMV_ASYM_DOUBLE(ie_op,'T',A_EZX,EX,EZ,Ix,Iy,-C_ONE,C_ONE,C_ONE,p_in,p_out)

                                CALL IE_OP_ZGEMV_ASYM_DOUBLE(ie_op,'N',A_EYZ,EZ,EY,Ix,Iy,MirrorY,MirrorY,C_ONE,p_in,p_out)
                                CALL IE_OP_ZGEMV_ASYM_DOUBLE(ie_op,'T',A_EZY,EY,EZ,Ix,Iy,-MirrorY,-MirrorY,C_ONE,p_in,p_out)
                        ENDDO
                ENDDO
               !$OMP END DO
       	!$OMP END  PARALLEL
        ENDSUBROUTINE

        SUBROUTINE VERTICAL_MULT_GENERAL_ZERO(ie_op,MirrorY,p_in,p_out)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		COMPLEX(REALPARM),POINTER,INTENT(IN)::p_in(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::p_out(:,:,:,:,:)
                COMPLEX(REALPARM),INTENT(IN)::MirrorY
                INTEGER::Iy

		!$OMP PARALLEL	PRIVATE(Iy) DEFAULT(SHARED)
                !$OMP DO SCHEDULE(GUIDED)
                DO Iy=1,ie_op%Ny_loc
                                CALL IE_OP_ZGEMV_SYMM_ZERO(ie_op,S_EXX,EX,EX,Iy,C_ONE,C_ZERO,p_in,p_out)

                                CALL IE_OP_ZGEMV_SYMM_ZERO(ie_op,S_EXY,EY,EX,Iy,MirrorY,C_ONE,p_in,p_out)
                                CALL IE_OP_ZGEMV_SYMM_ZERO(ie_op,S_EYX,EX,EY,Iy,MirrorY,C_ZERO,p_in,p_out)

                                CALL IE_OP_ZGEMV_SYMM_ZERO(ie_op,S_EYY,EY,EY,Iy,C_ONE,C_ONE,p_in,p_out)


                                CALL IE_OP_ZGEMV_SYMM_ZERO(ie_op,S_EZZ,EZ,EZ,Iy,C_ONE,C_ZERO,p_in,p_out)

                                CALL IE_OP_ZGEMV_ASYM_ZERO(ie_op,'N',A_EXZ,EZ,EX,Iy,C_ONE,C_ONE,p_in,p_out)
                                CALL IE_OP_ZGEMV_ASYM_ZERO(ie_op,'T',A_EZX,EX,EZ,Iy,-C_ONE,C_ONE,p_in,p_out)

                                CALL IE_OP_ZGEMV_ASYM_ZERO(ie_op,'N',A_EYZ,EZ,EY,Iy,MirrorY,C_ONE,p_in,p_out)
                                CALL IE_OP_ZGEMV_ASYM_ZERO(ie_op,'T',A_EZY,EY,EZ,Iy,-MirrorY,C_ONE,p_in,p_out)

                ENDDO
              !$OMP END DO
       	!$OMP END  PARALLEL
        ENDSUBROUTINE


	SUBROUTINE IE_OP_ZGEMV_SYMM_ZERO(ie_op,Tc,c_in,c_out,Iy,ALPHA,BETA,p_in,p_out)
			TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_in(:,:,:,:,:)
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_out(:,:,:,:,:)
			INTEGER,INTENT(IN)::Tc,c_in,c_out,Iy
			COMPLEX(REALPARM),INTENT(IN)::ALPHA,BETA
                        
			CALL ZSPMV('U',ie_op%Nz,ALPHA,ie_op%G_symmT(:,Tc,Iy,ZERO),&
			&p_in(:,1,c_in,Iy,ZERO),ONE,BETA,p_out(:,1,c_out,Iy,ZERO),ONE)
	END SUBROUTINE
	SUBROUTINE IE_OP_ZGEMV_ASYM_ZERO(ie_op,TRANS,Tc,c_in,c_out,Iy,ALPHA,BETA,p_in,p_out)
			TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_in(:,:,:,:,:)
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_out(:,:,:,:,:)
			CHARACTER,INTENT(IN)::TRANS(*)
			INTEGER,INTENT(IN)::Tc,c_in,c_out,Iy
			COMPLEX(REALPARM),INTENT(IN)::ALPHA,BETA
			CALL ZGEMV(TRANS,ie_op%Nz,ie_op%Nz,ALPHA,ie_op%G_asymT(:,:,Tc,Iy,ZERO),&
                        &ie_op%Nz,p_in(:,1,c_in,Iy,ZERO),ONE,BETA,p_out(:,1,c_out,Iy,ZERO),ONE)
	END SUBROUTINE
	SUBROUTINE IE_OP_ZGEMV_SYMM_DOUBLE(ie_op,Tc,c_in,c_out,Ix,Iy,ALPHA1,ALPHA2,BETA,p_in,p_out)
			TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_in(:,:,:,:,:)
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_out(:,:,:,:,:)
			INTEGER,INTENT(IN)::Tc,c_in,c_out,Ix,Iy
			COMPLEX(REALPARM),INTENT(IN)::ALPHA1,ALPHA2,BETA
                        
			CALL ZSPMV('U',ie_op%Nz,ALPHA1,ie_op%G_symmT(:,Tc,Iy,Ix),&
			&p_in(:,1,c_in,Iy,Ix),ONE,BETA,p_out(:,1,c_out,Iy,Ix),ONE)
			CALL ZSPMV('U',ie_op%Nz,ALPHA2,ie_op%G_symmT(:,Tc,Iy,Ix),&
			&p_in(:,2,c_in,Iy,Ix),ONE,BETA,p_out(:,2,c_out,Iy,Ix),ONE)
	END SUBROUTINE

	SUBROUTINE  IE_OP_ZGEMV_ASYM_DOUBLE(ie_op,TRANS,Tc,c_in,c_out,Ix,Iy,ALPHA1,ALPHA2,BETA,p_in,p_out)
			TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_in(:,:,:,:,:)
                        COMPLEX(REALPARM),POINTER,INTENT(IN)::p_out(:,:,:,:,:)
			CHARACTER,INTENT(IN)::TRANS(*)
			INTEGER,INTENT(IN)::Tc,c_in,c_out,Ix,Iy
			COMPLEX(REALPARM),INTENT(IN)::ALPHA1,ALPHA2,BETA
			CALL ZGEMV(TRANS,ie_op%Nz,ie_op%Nz,ALPHA1,ie_op%G_asymT(:,:,Tc,Iy,Ix),&
                        &ie_op%Nz,p_in(:,1,c_in,Iy,Ix),ONE,BETA,p_out(:,1,c_out,Iy,Ix),ONE)
			CALL ZGEMV(TRANS,ie_op%Nz,ie_op%Nz,ALPHA2,ie_op%G_asymT(:,:,Tc,Iy,Ix),&
                        &ie_op%Nz,p_in(:,2,c_in,Iy,Ix),ONE,BETA,p_out(:,2,c_out,Iy,Ix),ONE)
	END SUBROUTINE
!--------------------------------------------------------------------------------------------------------------------!       
#ifdef gggg
        SUBROUTINE VERTICAL_MULT_UNIFORM(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
                TYPE(LOCAL_OMP_FFT_DATA),POINTER::LFFT
                COMPLEX(REALPARM),POINTER::data_in(:,:),data_out(:,:)
                INTEGER::Ixy,TN,Ix,M,I,J,Nz,Nz2,K
                INTEGER,PARAMETER::REZ=EZ+1

                INTEGER,PARAMETER::TEX=2*EX-1
                INTEGER,PARAMETER::HEX=2*EX

                INTEGER,PARAMETER::TEY=2*EY-1
                INTEGER,PARAMETER::HEY=2*EY

                INTEGER,PARAMETER::TEZ=2*EZ-1
                INTEGER,PARAMETER::HEZ=2*EZ

                INTEGER,PARAMETER::T=1
                INTEGER,PARAMETER::H=2
                LFFT=>ie_op%LFFT
                Nz=ie_op%Nz
                Nz2=2*Nz
		!$OMP PARALLEL	PRIVATE(Ixy,Ix,M,I,J,TN,data_in,data_out) DEFAULT(SHARED)
                TN=OMP_GET_THREAD_NUM()
                data_in(1:2*Nz,1:6)=>LFFT%data_in(:,TN)
                data_out(1:2*Nz,1:6)=>LFFT%data_out(:,TN)
		!$OMP DO SCHEDULE(GUIDED) 
		DO I=1,2*ie_op%NxNy_loc
                        Ix=MODULO(I-1,2*ie_op%Nx)
                        M=(I-Ix-1)/2
                        Ix=MIN(Ix,2*ie_op%Nx-Ix)
                        IF (Ix==ie_op%Nx) THEN
                                ie_op%field_out3(:,:,I)=C_ZERO
                                CYCLE
                        ENDIF

                        Ixy=Ix+M+1

                       data_in=C_ZERO
                       data_in(1:Nz,EX)=ie_op%field_in3(:,EX,I)
                       data_in(1:Nz,EY)=ie_op%field_in3(:,EY,I)
                       data_in(1:Nz,EZ)=ie_op%field_in3(:,EZ,I)
                       DO K=1,Nz
                               data_in(K,REZ)=ie_op%field_in3(Nz-K+1,EZ,I)
                       ENDDO
                       CALL CALCULATE_FORWARD_AT_THREAD(LFFT,TN)

!---------------------------------------------     EX  ----------------------------------------                               
                       data_in(:,TEX)=ie_op%G_symm5(:,T,S_EXX,Ixy)*data_out(:,EX)
                       data_in(:,TEX)=ie_op%G_symm5(:,T,S_EXY,Ixy)*data_out(:,EY)+data_in(:,TEX)

                       data_in(:,HEX)=ie_op%G_symm5(:,H,S_EXX,Ixy)*data_out(:,EX)
                       data_in(:,HEX)=ie_op%G_symm5(:,H,S_EXY,Ixy)*data_out(:,EY)+data_in(:,HEX)

                       data_in(:,HEX)=ie_op%G_asym4(:,H,A_EXZ,Ixy)*data_out(:,EZ)+data_in(:,HEX)
                       data_in(:,HEX)=-ie_op%G_asym4(:,T,A_EXZ,Ixy)*data_out(:,REZ)+data_in(:,HEX)
!---------------------------------------------     EY  ----------------------------------------                               
                       data_in(:,TEY)=ie_op%G_symm5(:,T,S_EXY,Ixy)*data_out(:,EX)
                       data_in(:,TEY)=ie_op%G_symm5(:,T,S_EYY,Ixy)*data_out(:,EY)+data_in(:,TEY)

                       data_in(:,HEY)=ie_op%G_symm5(:,H,S_EXY,Ixy)*data_out(:,EX)
                       data_in(:,HEY)=ie_op%G_symm5(:,H,S_EYY,Ixy)*data_out(:,EY)+data_in(:,HEY)

                       data_in(:,HEY)=ie_op%G_asym4(:,H,A_EYZ,Ixy)*data_out(:,EZ)+data_in(:,HEY)
                       data_in(:,HEY)=-ie_op%G_asym4(:,T,A_EYZ,Ixy)*data_out(:,REZ)+data_in(:,HEY)

!---------------------------------------------     EZ  ----------------------------------------                               
                       data_in(:,TEZ)=ie_op%G_asym4(:,T,A_EZX,Ixy)*data_out(:,EX)
                       data_in(:,TEZ)=ie_op%G_asym4(:,T,A_EZY,Ixy)*data_out(:,EY)+data_in(:,TEZ)
                       data_in(:,TEZ)=ie_op%G_symm5(:,T,S_EZZ,Ixy)*data_out(:,EZ)+data_in(:,TEZ)

                       data_in(:,HEZ)=ie_op%G_asym4(:,H,A_EZX,Ixy)*data_out(:,EX)
                       data_in(:,HEZ)=ie_op%G_asym4(:,H,A_EZY,Ixy)*data_out(:,EY)+data_in(:,HEZ)
                       data_in(:,HEZ)=ie_op%G_symm5(:,H,S_EZZ,Ixy)*data_out(:,EZ)+data_in(:,HEZ)

!-----------------------------------------------------------------------------------------                               
                       CALL CALCULATE_BACKWARD_AT_THREAD(LFFT,TN)

                       DO K=1,Nz
                                ie_op%field_out3(K,EX,I)=data_out(K,TEX)+data_out(Nz-K+1,HEX)
                       ENDDO
                       DO K=1,Nz
                                ie_op%field_out3(K,EY,I)=data_out(K,TEY)+data_out(Nz-K+1,HEY)
                       ENDDO
                       DO K=1,Nz
                                ie_op%field_out3(K,EZ,I)=data_out(K,TEZ)+data_out(Nz-K+1,HEZ)
                       ENDDO
 
		ENDDO
		!$OMP END DO
		!$OMP END  PARALLEL
        ENDSUBROUTINE
#endif
END MODULE  
