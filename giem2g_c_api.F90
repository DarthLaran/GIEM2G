MODULE GIEM2G_C_API
	USE CONST_MODULE
	USE, INTRINSIC :: iso_c_binding
	USE DATA_TYPES_MODULE

        USE INTEGRAL_EQUATION_MODULE
	USE Calc_IE_Tensor_Module
	USE IE_SOLVER_MODULE

	USE CONTINUATION_FUNCTION_MODULE
	USE Calc_RC_Tensor_Module
	USE MPI_MODULE
	IMPLICIT NONE


        CONTAINS
!-------------------------------- Create Structures ------------------!
                FUNCTION GIEM2G_CREATE_BKG(N,p_sig,p_thick)RESULT(RES)
                        INTEGER(C_INT),INTENT(IN)::N
                        TYPE(C_PTR),INTENT(IN)::p_sig,p_thick
                        TYPE(C_PTR)::RES
        		TYPE(BKG_DATA_TYPE),POINTER::bkg
                        REAL(REALPARM),POINTER::sig(:), thick(:)
			INTEGER::I
                        ALLOCATE(bkg)
                        RES=C_LOC(bkg)
			bkg%Nl=N
			ALLOCATE (bkg%sigma(N), bkg%thick(1:N-1), bkg%depth(1:N-1),bkg%csig(1:N-1))
                        CALL C_F_POINTER(p_sig,sig,(/N/))
                        bkg%sigma=sig
        		IF (N>1) THEN
                                CALL C_F_POINTER(p_thick,thick,(/N-1/))
                                bkg%thick=thick
        			bkg%depth(1)=bkg%thick(1)
        			DO I=2,N-1
	        			bkg%depth(I)=bkg%depth(I-1)+bkg%thick(I)
		        	END DO
        		ENDIF
                END FUNCTION
                FUNCTION GIEM2G_CREATE_ANOMALY(Nx,Ny,Nz,dx,dy,c_z,c_bkg)RESULT(RES)
        		INTEGER(C_INT),INTENT(IN) :: Nx,Ny,Nz
                        REAL(REALPARM),INTENT(IN)::dx,dy
                        TYPE(C_PTR),INTENT(IN)::c_z,c_bkg
                        TYPE(C_PTR)::RES
        		TYPE(ANOMALY_TYPE),POINTER::anomaly
        		TYPE(BKG_DATA_TYPE),POINTER::bkg
                        ALLOCATE(anomaly)
                        RES=C_LOC(anomaly)
                        CALL C_F_POINTER(c_bkg,bkg)
        		anomaly%Nx=Nx
	        	anomaly%Ny=Ny
	        	anomaly%Nz=Nz
	        	anomaly%dx=dx
	        	anomaly%dy=dy
			ALLOCATE(anomaly%z(0:Nz),anomaly%Lnumber(0:Nz))
			anomaly%Lnumber(1:Nz)=GetLayer((anomaly%z(0:Nz-1)+anomaly%z(1:Nz))/2d0,bkg)
			anomaly%Lnumber(0)=GetLayer(anomaly%z(0)*(1d0-1d-7),bkg)
                END FUNCTION
                FUNCTION GIEM2G_CREATE_RECIEVERS(Nr,c_x,c_y,c_z,c_anomaly,c_bkg)RESULT(RES)
                        INTEGER(C_INT),INTENT(IN)::Nr! Number of patterns!
                        TYPE(C_PTR),INTENT(IN)::c_z,c_x,c_y,c_anomaly,c_bkg
			TYPE(C_PTR)::RES
	        	REAL(REALPARM),POINTER::zr(:),xshift(:),yshift(:)
	        	TYPE (RECEIVER_TYPE),POINTER::recvs(:)
                	TYPE(ANOMALY_TYPE),POINTER::anomaly
        		TYPE(BKG_DATA_TYPE),POINTER::bkg
			INTEGER::I
                        CALL C_F_POINTER(c_bkg,bkg)
                        CALL C_F_POINTER(c_anomaly,anomaly)
                        ALLOCATE(recvs(Nr))
                        RES=C_LOC(recvs)
                        CALL C_F_POINTER(c_z,zr,(/Nr/))
                        CALL C_F_POINTER(c_x,xshift,(/Nr/))
                        CALL C_F_POINTER(c_y,yshift,(/Nr/))
        		DO I=1,Nr
	        		recvs(I)%zrecv=zr(I)
		        	recvs(I)%x_shift=xshift(I)
		        	recvs(I)%y_shift=yshift(I)
		        ENDDO
                	CALL PrepareRecvs(recvs,anomaly,bkg)
        	END FUNCTION 
                FUNCTION GIEM2G_CREATE_IE_OPERATOR(c_anomaly,wcomm,threads_ok) RESULT(RES)
                        TYPE(C_PTR),INTENT(IN)::c_anomaly
			INTEGER(MPI_CTL_KIND),INTENT(IN)::wcomm
                        LOGICAL,INTENT(IN)::threads_ok
                	TYPE(IntegralEquation),POINTER:: int_eq
                	TYPE(ANOMALY_TYPE),POINTER::anomaly
			TYPE(C_PTR)::RES
                        ALLOCATE(int_eq)
                        RES=C_LOC(int_eq)
                        CALL C_F_POINTER(c_anomaly,anomaly)
                	CALL PrepareIntegralEquation(int_eq,anomaly,wcomm,threads_ok)
                	IF (int_eq%real_space) THEN
				ALLOCATE(int_eq%siga(int_eq%Nz,int_eq%Nx,int_eq%Ny_loc))
                	ENDIF
                END FUNCTION
                FUNCTION GIEM2G_CREATE_RC_OPERATOR(c_anomaly,Nr,c_recvs,wcomm,threads_ok) RESULT(RES)
                        TYPE(C_PTR),INTENT(IN)::c_anomaly
                        TYPE(C_PTR),INTENT(IN)::c_recvs
                        INTEGER(C_INT),INTENT(IN)::Nr
			INTEGER(MPI_CTL_KIND),INTENT(IN)::wcomm
                        LOGICAL,INTENT(IN)::threads_ok
			TYPE(C_PTR)::RES
                	TYPE (RECEIVER_TYPE),POINTER::recvs(:)
                	TYPE(RC_OPERATOR),POINTER::rc_op
                	TYPE(ANOMALY_TYPE),POINTER::anomaly
                        ALLOCATE(rc_op)
                        RES=C_LOC(rc_op)
                        CALL C_F_POINTER(c_anomaly,anomaly)
                        CALL C_F_POINTER(c_recvs,recvs,(/Nr/))
                	CALL PrepareContinuationOperator(rc_op,anomaly,recvs,wcomm,threads_ok)
                	IF (rc_op%real_space) THEN
				ALLOCATE(rc_op%siga(rc_op%Nz,rc_op%Nx,rc_op%Ny_loc))
                	ENDIF
                END FUNCTION
!------------------------------------- Delete structures --------------------!
                SUBROUTINE GIEM2G_DeleteIE_OP(c_ie)
                        TYPE(C_PTR),INTENT(IN)::c_ie
                	TYPE(IntegralEquation),POINTER:: int_eq
                        CALL C_F_POINTER(c_ie,int_eq)
                	CALL DeleteIE_OP(int_eq)
                END SUBROUTINE
                SUBROUTINE GIME2G_DeleteRC_OP(c_rc)
                        TYPE(C_PTR),INTENT(IN)::c_rc
                	TYPE(RC_OPERATOR),POINTER:: rc_op
                        CALL C_F_POINTER(c_rc,rc_op)
                	CALL DeleteRC_OP(RC_OP)
                END SUBROUTINE
!-----------------------------------------------------------------------------------!
!---------------------------- Calculations--------------------------------!
                SUBROUTINE GIEM2G_CALC_IE_TENSOR(c_ie,c_bkg,c_anomaly,freq,IE_Threshold)
                        TYPE(C_PTR),INTENT(IN)::c_ie,c_bkg,c_anomaly
                        REAL(REALPARM),INTENT(IN)::freq,IE_Threshold

                	TYPE(IntegralEquation),POINTER:: int_eq
        		TYPE(ANOMALY_TYPE),POINTER::anomaly
        		TYPE(BKG_DATA_TYPE),POINTER::bkg
			INTEGER::Nz
                        CALL C_F_POINTER(c_bkg,bkg)
                        CALL C_F_POINTER(c_anomaly,anomaly)
                        CALL C_F_POINTER(c_ie,int_eq)
                	CALL SetSigb(int_eq,anomaly,bkg)
        		CALL Set_Freq(bkg,freq)
			Nz=anomaly%Nz
                        IF (anomaly%Nz == int_eq%Nz )THEN
                                int_eq%dz=anomaly%z(1:)-anomaly%z(0:Nz-1)
	                	CALL CalcIntegralGreenTensor(int_eq,bkg,anomaly,IE_Threshold)
                        ELSE
                                PRINT*, 'Different shapes of anomaly and IE  operator'
                        ENDIF
                END SUBROUTINE

                SUBROUTINE GIEM2G_CALC_RC_TENSOR(c_rc,c_bkg,c_anomaly,freq,RC_Threshold)
                        TYPE(C_PTR),INTENT(IN)::c_rc,c_bkg,c_anomaly
                        REAL(REALPARM),INTENT(IN)::freq,RC_Threshold
                	TYPE(RC_OPERATOR),POINTER:: rc_op
        		TYPE(ANOMALY_TYPE),POINTER::anomaly
        		TYPE(BKG_DATA_TYPE),POINTER::bkg
                        CALL C_F_POINTER(c_bkg,bkg)
                        CALL C_F_POINTER(c_anomaly,anomaly)
                        CALL C_F_POINTER(c_rc,rc_op)
                	CALL SetSigbRC(rc_op,anomaly,bkg)
        		CALL Set_Freq(bkg,freq)
                	CALL CalcRecalculationGreenTensor(rc_op,bkg,anomaly,RC_Threshold)
                END SUBROUTINE

                SUBROUTINE GIEM2G_CALC_FFT_OF_IE_TENSOR(c_ie)
                        TYPE(C_PTR),INTENT(IN)::c_ie
                	TYPE(IntegralEquation),POINTER:: int_eq
                        CALL C_F_POINTER(c_ie,int_eq)
        		CALL CalcFFTofIETensor(int_eq)
                END SUBROUTINE
                SUBROUTINE GIEM2G_CALC_FFT_OF_RC_TENSOR(c_rc)
                        TYPE(C_PTR),INTENT(IN)::c_rc
                	TYPE(RC_OPERATOR),POINTER:: rc_op
                        CALL C_F_POINTER(c_rc,rc_op)
        		CALL CalcFFTofRCTensor(rc_op)
                END SUBROUTINE
!----------------------------------------------------------------------------------------------------------!
		SUBROUTINE GIEM2G_SET_ANOMALY_FOR_IE(c_ie,c_siga)
                        TYPE(C_PTR),INTENT(IN)::c_ie,c_siga
                	TYPE(IntegralEquation),POINTER:: int_eq
			INTEGER::siga_shape(3)
                        CALL C_F_POINTER(c_ie,int_eq)
			siga_shape=(/int_eq%Nz,int_eq%Nx,int_eq%Ny_loc/)
			IF (int_eq%real_space) THEN
	                        CALL C_F_POINTER(c_siga,int_eq%siga,siga_shape)
			ELSE
				int_eq%siga=>NULL()
			ENDIF
		END SUBROUTINE
		SUBROUTINE GIEM2G_SET_ANOMALY_FOR_RC(c_rc,c_siga)
                        TYPE(C_PTR),INTENT(IN)::c_rc,c_siga
                	TYPE(RC_OPERATOR),POINTER:: rc_op
			INTEGER::siga_shape(3)
                        CALL C_F_POINTER(c_rc,rc_op)
			siga_shape=(/rc_op%Nz,rc_op%Nx,rc_op%Ny_loc/)
			IF (rc_op%real_space) THEN
	                        CALL C_F_POINTER(c_siga,rc_op%siga,siga_shape)
			ELSE
				rc_op%siga=>NULL()
			ENDIF
		END SUBROUTINE
!----------------------------------------------------------------------------------------------------------!
		SUBROUTINE  GIEM2G_Solve_Equation(c_ie,c_en,c_eav,maxit,misfit,fgmres_buf,gmres_buf)
                        !!! ATTENTION c_en is normal electric field averaged over the domains
                        TYPE(C_PTR),INTENT(IN)::c_ie,c_en,c_eav
        		REAL(REALPARM)::misfit
	        	INTEGER,INTENT(IN)::maxit
	        	INTEGER,INTENT(IN)::gmres_buf
	        	INTEGER,INTENT(IN)::fgmres_buf
                	TYPE (FGMRES_CTL_TYPE)::fgmres_ctl
                	TYPE(IntegralEquation),POINTER:: int_eq
                        COMPLEX(REALPARM),POINTER::En(:),Esol(:)
                        CALL C_F_POINTER(c_ie,int_eq)
                        IF (int_eq%real_space) THEN
                                CALL C_F_POINTER(c_en,En,(/int_eq%Nloc/))
                                int_eq%rhs=En
                        ENDIF
                        fgmres_ctl%fgmres_maxit=maxit
                        fgmres_ctl%misfit=misfit
                        fgmres_ctl%fgmres_buf=fgmres_buf
                        fgmres_ctl%gmres_buf=gmres_buf
			CALL SolveEquation(int_eq,fgmres_ctl)
                        IF (int_eq%real_space) THEN
                                CALL C_F_POINTER(c_eav,Esol,(/int_eq%Nloc/))
                                Esol=int_eq%solution
                        ENDIF
                END SUBROUTINE

		SUBROUTINE  GIEM2G_Solve_Equation_with_Guess(c_ie,c_en,c_guess,c_eav,maxit,misfit,fgmres_buf,gmres_buf)
                        !!! ATTENTION c_en is normal electric field averaged over the domains
                        TYPE(C_PTR),INTENT(IN)::c_ie,c_en
                        TYPE(C_PTR),INTENT(IN)::c_guess,c_eav
        		REAL(REALPARM)::misfit
	        	INTEGER,INTENT(IN)::maxit
	        	INTEGER,INTENT(IN)::gmres_buf
	        	INTEGER,INTENT(IN)::fgmres_buf
                	TYPE (FGMRES_CTL_TYPE)::fgmres_ctl
                	TYPE(IntegralEquation),POINTER:: int_eq
                        COMPLEX(REALPARM),POINTER::En(:),Esol(:),Eguess(:,:,:,:)
                        INTEGER::guess_shape(4)
                        CALL C_F_POINTER(c_ie,int_eq)
                        IF (int_eq%real_space) THEN
                                CALL C_F_POINTER(c_en,En,(/int_eq%Nloc/))
                                int_eq%rhs=En
                                guess_shape=(/int_eq%Nz,3,int_eq%Nx,int_eq%Ny_loc/)
                                CALL C_F_POINTER(c_guess,Eguess,guess_shape)
                        ELSE
                                Eguess=>NULL()
                        ENDIF
                        fgmres_ctl%fgmres_maxit=maxit
                        fgmres_ctl%misfit=misfit
                        fgmres_ctl%fgmres_buf=fgmres_buf
                        fgmres_ctl%gmres_buf=gmres_buf
			CALL SolveEquation(int_eq,fgmres_ctl,Eguess)
                        IF (int_eq%real_space) THEN
                                CALL C_F_POINTER(c_eav,Esol,(/int_eq%Nloc/))
                                Esol=int_eq%solution
                        ENDIF
                END SUBROUTINE

	        SUBROUTINE GIEM2G_ReCalculation (c_rc,c_eav,c_Ea,c_Ha)
                        TYPE(C_PTR),INTENT(IN)::c_rc,c_eav,c_Ea,c_Ha
                	TYPE(RC_OPERATOR),POINTER::rc_op
                        COMPLEX(REALPARM),POINTER::Eint(:,:,:,:)
                        COMPLEX(REALPARM),POINTER::Ea(:,:,:,:),Ha(:,:,:,:)
                        INTEGER::shape_in(4),shape_out(4)
                        CALL C_F_POINTER(c_rc,rc_op)
                        IF (rc_op%real_space) THEN
                                shape_in=(/rc_op%Nz,3,rc_op%Nx,rc_op%Ny_loc/)
                                shape_out=(/rc_op%Nr,3,rc_op%Nx,rc_op%Ny_loc/)
                                CALL C_F_POINTER(c_eav,Eint,shape_in)
                                CALL C_F_POINTER(c_Ea,Ea,shape_out)
                                CALL C_F_POINTER(c_Ha,Ha,shape_out)
                        ELSE
                                Ea=>NULL()
                                Ha=>NULL()
                                Eint=>NULL()
                        ENDIF
			CALL ReCalculation(rc_op,Eint,Ea,Ha)
                END SUBROUTINE
!-----------------------------------------------------------------------------------!
!------------------------------------- Guru Interface ))---------------------!!!!!
                
                SUBROUTINE GIEM2G_GET_IE_MATRICIES_POINTERS(c_ie,c_gsymm,c_gasym)
                        TYPE(C_PTR),INTENT(IN)::c_ie
                        TYPE(C_PTR),INTENT(OUT)::c_gsymm,c_gasym
                	TYPE(IntegralEquation),POINTER:: int_eq
                        CALL C_F_POINTER(c_ie,int_eq)
                        c_gsymm=int_eq%pG_symm
                        c_gasym=int_eq%pG_asym
                END SUBROUTINE
                SUBROUTINE GIEM2G_GET_RC_MATRICIES_POINTERS(c_rc,c_e,c_h)
                        TYPE(C_PTR),INTENT(IN)::c_rc
                        TYPE(C_PTR),INTENT(OUT)::c_e,c_h
                	TYPE(RC_OPERATOR),POINTER:: rc_op
                        CALL C_F_POINTER(c_rc,rc_op)
                        c_e=rc_op%pG_E
                        c_h=rc_op%pG_H
                END SUBROUTINE
END MODULE
