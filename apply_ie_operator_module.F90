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


	SUBROUTINE APPLY_PRECOND(int_eq,v_in,v_out)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::int_eq
		COMPLEX(REALPARM),POINTER,INTENT(IN)::v_in(:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::v_out(:)
		COMPLEX(REALPARM),POINTER::field_in(:,:,:,:)
		COMPLEX(REALPARM),POINTER::field_out(:,:,:,:)
		COMPLEX(REALPARM)::d1,d2
		COMPLEX(REALPARM)::asiga,dsig,gsig
		REAL(8)::time1,time2
		INTEGER::IERROR
		INTEGER ::Ix,Iy,Iz,Ic,Ixy
		time1=GetTime()
		field_in(1:int_eq%Nz,1:3,1:int_eq%Nx,1:int_eq%Ny_loc)=>v_in
		field_out(1:int_eq%Nz,1:3,1:int_eq%Nx,1:int_eq%Ny_loc)=>v_out
		DO Iy=1,int_eq%Ny_loc
		!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(Ix,Ic,Iz,asiga,dsig,gsig)
			!$OMP DO SCHEDULE(GUIDED)
			DO Ix=1,int_eq%Nx
				DO Ic=1,3
					DO Iz=1,int_eq%Nz
						asiga=C_TWO*int_eq%sqsigb(Iz)/(int_eq%csiga(Iz,Ix,Iy)+CONJG(int_eq%csigb(Iz)))
						dsig=int_eq%csiga(Iz,Ix,Iy)-int_eq%csigb(Iz)
						gsig=dsig*asiga
						int_eq%field_in4(Iz,Ic,Ix,Iy)=&
						&field_in(Iz,Ic,Ix,Iy)*gsig
					ENDDO
				ENDDO
			ENDDO
			!$OMP END DO
			!$OMP DO SCHEDULE(GUIDED)
			DO Ix=int_eq%Nx+1,int_eq%Nx_loc
				DO Ic=1,3
					DO Iz=1,int_eq%Nz
						int_eq%field_in4(Iz,Ic,Ix,Iy)=C_ZERO
					ENDDO
				ENDDO
			ENDDO
			!$OMP END DO
		!$OMP END PARALLEL
		ENDDO
		CALL	MULT_IE_OP(int_eq)
		DO Iy=1,int_eq%Ny_loc
			!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(Ix,Ic,Iz,d1,d2,asiga)
			!$OMP DO SCHEDULE(GUIDED)
			DO Ix=1,int_eq%Nx
				DO Ic=1,3
					DO Iz=1,int_eq%Nz
						asiga=C_TWO*int_eq%sqsigb(Iz)/(int_eq%csiga(Iz,Ix,Iy)+CONJG(int_eq%csigb(Iz)))
						d1=field_in(Iz,Ic,Ix,Iy)*asiga
						d2=d1-int_eq%field_out4(Iz,Ic,Ix,Iy)/(int_eq%dz(Iz))
						field_out(Iz,Ic,Ix,Iy)=d2*int_eq%sqsigb(Iz)
					ENDDO
				ENDDO
			ENDDO
			!$OMP END DO
			!$OMP END PARALLEL
		ENDDO
		time2=GetTime()
		int_eq%counter%mult_num=int_eq%counter%mult_num+1
		int_eq%counter%apply=int_eq%counter%apply+time2-time1
	END SUBROUTINE
	SUBROUTINE APPLY_IE_ZEROS(int_eq)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::int_eq

		!$OMP PARALLEL	 DEFAULT(SHARED)
		!$OMP WORKSHARE
		int_eq%field_in4=C_ZERO
		!$OMP END WORKSHARE
		!$OMP ENDPARALLEL
		CALL MULT_IE_OP(int_eq)
	ENDSUBROUTINE
	SUBROUTINE MULT_IE_OP(ie_op)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		INTEGER::Ixy,Nz
		REAL(8)::time1,time2,time3,time4
		time1=GetTime()
		CALL IE_OP_FFTW_FWD(ie_op)
		time2=GetTime()
		ie_op%counter%mult_fftw=ie_op%counter%mult_fftw+time2-time1
		!$OMP PARALLEL	PRIVATE(Ixy) DEFAULT(SHARED)
		!$OMP DO SCHEDULE(GUIDED) 
		DO Ixy=1,ie_op%NxNy_loc
			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EXX,EX,EX,Ixy,C_ONE,C_ZERO)
			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EXY,EY,EX,Ixy,C_ONE,C_ONE)

			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EYX,EX,EY,Ixy,C_ONE,C_ZERO)
			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EYY,EY,EY,Ixy,C_ONE,C_ONE)

			CALL IE_OP_ZGEMV_SYMM(ie_op,S_EZZ,EZ,EZ,Ixy,C_ONE,C_ZERO)

			CALL IE_OP_ZGEMV_ASYM(ie_op,'N',A_EXZ,EZ,EX,Ixy,C_ONE,C_ONE)
			CALL IE_OP_ZGEMV_ASYM(ie_op,'T',A_EZX,EX,EZ,Ixy,-C_ONE,C_ONE)

			CALL IE_OP_ZGEMV_ASYM(ie_op,'N',A_EYZ,EZ,EY,Ixy,C_ONE,C_ONE)
			CALL IE_OP_ZGEMV_ASYM(ie_op,'T',A_EZY,EY,EZ,Ixy,-C_ONE,C_ONE)

		ENDDO
		!$OMP END DO
		!$OMP END  PARALLEL
		time3=GetTime()
		CALL IE_OP_FFTW_BWD(ie_op)
		time4=GetTime()
		ie_op%counter%mult_fftw_b=ie_op%counter%mult_fftw_b+time4-time3
		ie_op%counter%mult_zgemv=ie_op%counter%mult_zgemv+time3-time2
	ENDSUBROUTINE
	SUBROUTINE IE_OP_ZGEMV_SYMM(ie_op,Tc,c_in,c_out,I,ALPHA,BETA)
			TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
			INTEGER,INTENT(IN)::Tc,c_in,c_out,I
			COMPLEX(REALPARM),INTENT(IN)::ALPHA,BETA
			CALL ZSPMV('U',ie_op%Nz,ALPHA,ie_op%G_symm4(:,Tc,I),&
			&ie_op%field_in3(:,c_in,I),ONE,BETA,ie_op%field_out3(:,c_out,I),ONE)
	END SUBROUTINE
	SUBROUTINE IE_OP_ZGEMV_ASYM(ie_op,TRANS,Tc,c_in,c_out,I,ALPHA,BETA)
			TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
			CHARACTER,INTENT(IN)::TRANS(*)
			INTEGER,INTENT(IN)::Tc,c_in,c_out,I
			COMPLEX(REALPARM),INTENT(IN)::ALPHA,BETA
			CALL ZGEMV(TRANS,ie_op%Nz,ie_op%Nz,ALPHA,ie_op%G_asym4(:,:,Tc,I),ie_op%Nz,&
			&ie_op%field_in3(:,c_in,I),ONE,BETA,ie_op%field_out3(:,c_out,I),ONE)
	END SUBROUTINE
END MODULE  
