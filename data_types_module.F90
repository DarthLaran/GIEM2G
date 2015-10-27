MODULE Data_Types_Module
	USE CONST_MODULE
	IMPLICIT NONE
	PRIVATE
	TYPE BKG_DATA_TYPE
		INTEGER :: Nl
		REAL (REALPARM),POINTER:: sigma(:), thick(:), depth(:)
		COMPLEX(REALPARM),POINTER::csigma(:)
		REAL (REALPARM) :: freq,omega
		COMPLEX (RealParm),POINTER :: k(:),k2(:)
		COMPLEX(REALPARM)::iwm
	ENDTYPE

	TYPE ANOMALY_TYPE
		INTEGER :: Nx,Ny,Nz
		INTEGER::Ny_loc
		REAL(REALPARM)::dx,dy
		REAL(REALPARM),POINTER::z(:)
		INTEGER,POINTER::Lnumber(:)
		COMPLEX(REALPARM),POINTER::csiga(:,:,:)
		REAL(REALPARM),POINTER::siga(:,:,:)
	ENDTYPE

	TYPE RECEIVER_TYPE
		REAL(REALPARM)::zrecv,x_shift,y_shift
		INTEGER::recv_layer
		INTEGER:: anom_cell !anom_cell<1 recv is upper anomaly , anom_cell >Nz recv is lower anomaly
	END TYPE

	TYPE FGMRES_CTL_TYPE
		REAL(REALPARM)::misfit
		INTEGER::fgmres_maxit
		INTEGER::gmres_buf
		INTEGER::fgmres_buf
		INTEGER::ort_type
	END TYPE

	INTERFACE GetLayer
		MODULE PROCEDURE GetLayerSingle,GetLayerArray
	END INTERFACE

	PUBLIC:: BKG_DATA_TYPE,ANOMALY_TYPE, RECEIVER_TYPE,FGMRES_CTL_TYPE

	PUBLIC::GetLayer,SET_FREQ
	PUBLIC::AllocateSiga,SliceAnomaly
	PUBLIC::PrepareRecvs
	CONTAINS
	FUNCTION GetLayerArray(z,bkg) RESULT(lnum)
		REAL (KIND = RealParm), INTENT(IN) ::z(:)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		INTEGER ::lnum(SIZE(z))
		INTEGER ::I
		DO I=1,SIZE(z)
			lnum(I)=GetLayerSingle(z(I),bkg)
		ENDDO
	END FUNCTION
	FUNCTION GetLayerSingle(z,bkg) RESULT(lnum)
		REAL (KIND = RealParm), INTENT(IN) ::z
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		INTEGER ::lnum
		IF (ABS(z)<1d-15) THEN
			lnum=0
			RETURN
		ENDIF
		
		IF (z<0) THEN
			lnum=-1
			RETURN
		ENDIF
		
		IF (bkg%Nl==1) THEN
			lnum=1
			RETURN
		ENDIF
			
		lnum=1
		IF (z>bkg%depth(bkg%Nl-1)) THEN
			lnum=bkg%Nl
			RETURN
		ELSE
			DO WHILE (z>bkg%depth(lnum))
				lnum=lnum+1
			ENDDO
		ENDIF
	END FUNCTION

	SUBROUTINE SET_FREQ(bkg,f)
		TYPE (BKG_DATA_TYPE),INTENT(INOUT)::bkg
		REAL(KIND=RealParm),INTENT(IN)::f
		bkg%freq=f
		bkg%omega=f*2D0*PI
		bkg%csigma(0)=-C_IONE*bkg%omega*EPS0
		bkg%csigma(1:)=bkg%sigma-C_IONE*bkg%omega*EPS0
		bkg%k2=C_IONE*bkg%csigma*MU0*bkg%omega
		bkg%k=SQRT(bkg%k2)
		bkg%iwm=C_IONE*MU0*bkg%omega
	END SUBROUTINE
	SUBROUTINE AllocateSiga(anomaly)
			TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
				IF (ASSOCIATED(anomaly%siga)) DEALLOCATE(anomaly%siga)
				ALLOCATE(anomaly%siga(anomaly%Nz,anomaly%Nx,anomaly%Ny_loc))
	ENDSUBROUTINE 


	SUBROUTINE SliceAnomaly(anomaly,bkg,top,thick,N)
			TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
			TYPE (BKG_DATA_TYPE),TARGET,INTENT(IN)::bkg
			INTEGER,INTENT(IN)::N
			REAL(REALPARM),INTENT(IN)::top,thick
			INTEGER::I
			anomaly%Nz=N
			IF (ASSOCIATED(anomaly%z))	DEALLOCATE(anomaly%z)
			IF (ASSOCIATED(anomaly%Lnumber)) DEALLOCATE(anomaly%Lnumber)
			ALLOCATE(anomaly%z(0:N),anomaly%Lnumber(0:N))
			anomaly%z=(/(top+I*thick/N,I=0,N)/)
			anomaly%Lnumber(1:N)=GetLayer((anomaly%z(0:N-1)+anomaly%z(1:N))/2d0,bkg)
			anomaly%Lnumber(0)=GetLayer(anomaly%z(0)*(1d0-1d-4),bkg)
	ENDSUBROUTINE

	SUBROUTINE PrepareRecvs(recvs,anomaly,bkg)
			TYPE (RECEIVER_TYPE),POINTER,INTENT(INOUT)::recvs(:)
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			TYPE (BKG_DATA_TYPE),TARGET,INTENT(IN)::bkg
			INTEGER::I,Nr,l
			REAL(REALPARM)::zrecv
			Nr=SIZE(recvs)
			DO I=1,Nr
				zrecv=recvs(I)%zrecv
				recvs(I)%recv_layer=MAX(GetLayer(zrecv,bkg),1)
				IF (zrecv<=anomaly%z(0)) THEN
					recvs(I)%anom_cell=0
					
				ELSEIF (zrecv>=anomaly%z(anomaly%Nz)) THEN
					recvs(I)%anom_cell=anomaly%Nz+1
				ELSE
					l=1
					DO WHILE (zrecv>anomaly%z(l))
						l=l+1
					ENDDO
					recvs(I)%anom_cell=l
				ENDIF
			ENDDO
	END SUBROUTINE
ENDMODULE

