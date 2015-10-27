MODULE MPI_SAVELOAD_MODULE
	USE CONST_MODULE
	USE DATA_TYPES_MODULE
	USE MPI_MODULE
	IMPLICIT NONE
CONTAINS
	SUBROUTINE LoadBackground(bkg,comm,background)
		TYPE(BKG_DATA_TYPE),INTENT(INOUT)::bkg
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::background
		INTEGER(MPI_CTL_KIND)::me,IERROR
		INTEGER::N,I
		CALL MPI_COMM_RANK(comm, me, IERROR)
		IF (me==0) THEN
			OPEN(077,file=background)
			READ(077,*) N
			bkg%Nl=N
			ALLOCATE (bkg%sigma(N), bkg%thick(1:N-1), bkg%depth(1:N-1),bkg%csigma(0:N))
			READ(077,*) bkg%sigma
			READ(077,*) bkg%thick
			CLOSE(077)
		ENDIF
		CALL MPI_BCAST(N,1,MPI_INTEGER, 0, comm, IERROR)
		IF (me/=0) THEN 
			bkg%Nl=N
			ALLOCATE (bkg%sigma(N), bkg%thick(1:N-1), bkg%depth(1:N-1),bkg%csigma(0:N))
		ENDIF
		ALLOCATE (bkg%k(0:N),bkg%k2(0:N))
		CALL MPI_BCAST(bkg%sigma,N,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		IF (N>1) THEN
			CALL MPI_BCAST(bkg%thick,N-1,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
			bkg%depth(1)=bkg%thick(1)
			DO I=2,N-1
				bkg%depth(I)=bkg%depth(I-1)+bkg%thick(I)
			END DO
		ENDIF
		IF (me==0) THEN
			PRINT*,"Background has been loaded"
			PRINT*, "Number of layers:    ", N
			PRINT*, "Layers conductivity: ", bkg%sigma
			IF (N/=1) THEN
				PRINT*, "Layers thickness:    ", bkg%thick
			ENDIF
		ENDIF
	ENDSUBROUTINE
	SUBROUTINE LoadAnomalyShape(anomaly,bkg,comm,anomaly_shape,loadz)
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
	        TYPE (BKG_DATA_TYPE),TARGET,INTENT(IN)::bkg
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::anomaly_shape
		LOGICAL,OPTIONAL,INTENT(IN)::loadz
		INTEGER(MPI_CTL_KIND)::IERROR,me
		INTEGER::Nx,NY,Nz
		REAL(REALPARM)::dx,dy
		LOGICAL::lz
		lz=.TRUE.
		IF (PRESENT(loadz))  lz=loadz
		CALL MPI_COMM_RANK(comm, me, IERROR)
		IF (me==0) THEN
			OPEN(078,file=anomaly_shape)
			READ(078,*) Nx,Ny,Nz
			READ(078,*) dx,dy
		ENDIF
		CALL MPI_BCAST(Nx,1,MPI_INTEGER, 0, comm, IERROR)
		CALL MPI_BCAST(Ny,1,MPI_INTEGER, 0, comm, IERROR)
		CALL MPI_BCAST(Nz,1,MPI_INTEGER, 0, comm, IERROR)
		CALL MPI_BCAST(dx,1,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		CALL MPI_BCAST(dy,1,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		anomaly%Nx=Nx
		anomaly%Ny=Ny
		anomaly%Nz=Nz
		anomaly%dx=dx
		anomaly%dy=dy
		anomaly%z=>NULL()
		anomaly%Lnumber=>NULL()
		IF (lz) THEN
			ALLOCATE(anomaly%z(0:Nz),anomaly%Lnumber(0:Nz))
			IF (me==0) THEN
				READ(078,*) anomaly%z
			ENDIF
			CALL MPI_BCAST(anomaly%z,Nz+1,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
			anomaly%Lnumber(1:Nz)=GetLayer((anomaly%z(0:Nz-1)+anomaly%z(1:Nz))/2d0,bkg)
			anomaly%Lnumber(0)=GetLayer(anomaly%z(0)*(1d0-1d-7),bkg)
			IF (me==0) THEN
				IF  (GetLayer(anomaly%z(0)*(1d0-1d-7),bkg)/=GetLayer(anomaly%z(0)*(1d0+1d-7),bkg)) THEN
					PRINT*, 'z(0) is on layers border, is not it?'
				ENDIF
			ENDIF

		ENDIF
		IF (me==0) CLOSE(078)
	ENDSUBROUTINE
        SUBROUTINE LoadAnomalySigmaList(fname,comm,anom_list,Na)
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::fname
        CHARACTER(*),POINTER,INTENT(INOUT)::anom_list(:)
        CHARACTER(len=1024)::tmp
        INTEGER,INTENT(OUT)::Na
		INTEGER(MPI_CTL_KIND)::IERROR,me
		INTEGER::I,l,io,Nl
		CALL MPI_COMM_RANK(comm, me, IERROR)
		
        IF (me==0) THEN
				OPEN(078,file=fname)
                READ(078,*) Na
				PRINT*, 'Number of anomalies:', Na
               ALLOCATE(anom_list(Na))
               DO I=1,Na
                    READ(078,'(A)',SIZE=l,advance='NO',IOSTAT=io) tmp
                    anom_list(I)=tmp(1:l)
                   PRINT*,trim(anom_list(I))
		       ENDDO
					CLOSE(078)
	       ENDIF
			CALL MPI_BCAST(Na,1,MPI_INTEGER, 0, comm, IERROR)
		    IF (me/=0) 	ALLOCATE(anom_list(Na))
			Nl=Na*LEN(anom_list(1))
			CALL MPI_BCAST(anom_list,Nl,MPI_CHARACTER, 0, comm, IERROR)
        ENDSUBROUTINE

	SUBROUTINE LoadAnomalySigma(anomaly,comm,anomaly_sig)
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::anomaly_sig
		REAL(REALPARM),ALLOCATABLE::siga1(:,:,:)
		INTEGER(MPI_CTL_KIND)::IERROR,me,csize
		INTEGER::Nc,io,I
		INTEGER(MPI_CTL_KIND):: REC_STATUS(MPI_STATUS_SIZE)
		REAL(8)::time1,time2
		time1=MPI_WTIME()
		CALL MPI_COMM_RANK(comm, me, IERROR)
		CALL MPI_COMM_SIZE(comm, csize, IERROR)
		Nc=anomaly%Nx*anomaly%Nz*anomaly%Ny_loc
		IF (me==0) THEN
			PRINT*,'Loading sigma from ', anomaly_sig
			OPEN(078,file=anomaly_sig)
			ALLOCATE(siga1(anomaly%Nz,anomaly%Nx,anomaly%Ny_loc))
			READ(078,*) anomaly%siga
			DO I=1,csize-1
				READ(078,*) siga1
				CALL MPI_SEND(siga1, Nc, MPI_DOUBLE_PRECISION, I, I, comm, IERROR)
			ENDDO
			CLOSE(078)
		ELSE
			CALL MPI_RECV(anomaly%siga, Nc, MPI_DOUBLE_PRECISION,0, me,  comm,REC_STATUS, IERROR)
		ENDIF
		CALL MPI_BARRIER(comm,IERROR)
		time2=MPI_WTIME()
		IF (me==0) THEN
			DEALLOCATE(siga1)
			PRINT "(A, ES10.2E3, A)", 'Sigma loaded in' ,time2-time1, 's' 
		ENDIF
	ENDSUBROUTINE
	SUBROUTINE LoadFrequencies(freq,comm,frequencies)
		REAL(REALPARM),INTENT(INOUT),POINTER::freq(:)
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::frequencies
		INTEGER(MPI_CTL_KIND)::me,IERROR
		INTEGER:: Nfreq
		IF (ASSOCIATED(freq)) THEN
			DEALLOCATE(freq)
			freq=>NULL()
		ENDIF
		CALL MPI_COMM_RANK(comm, me, IERROR)
		IF (me==0) THEN
			OPEN(078,file=frequencies)
			READ(078,*) Nfreq
			ALLOCATE(freq(Nfreq))
			READ(078,*) freq
			CLOSE(078)
		ENDIF
		CALL MPI_BCAST(Nfreq,1,MPI_INTEGER, 0, comm, IERROR)
		IF (me/=0) THEN
			ALLOCATE(freq(Nfreq))
		ENDIF
		CALL MPI_BCAST(freq,Nfreq,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
	ENDSUBROUTINE

	SUBROUTINE LoadRecievers(recvs,comm,recievers)
		TYPE (RECEIVER_TYPE),POINTER,INTENT(INOUT)::recvs(:)
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::recievers
		INTEGER(MPI_CTL_KIND):: me,IERROR
		INTEGER:: Nr,I
		REAL(REALPARM),ALLOCATABLE::zr(:),xshift(:),yshift(:)
		IF (ASSOCIATED(recvs)) THEN
			DEALLOCATE(recvs)
			recvs=>NULL()
		ENDIF
		CALL MPI_COMM_RANK(comm, me, IERROR)
		IF (me==0) THEN
			OPEN(078,file=recievers,form='FORMATTED')
			READ(078,*) Nr
			ALLOCATE(zr(Nr),xshift(Nr),yshift(Nr))
			DO I=1,Nr
				READ(078,*) zr(I),xshift(I),yshift(I)
			ENDDO
			CLOSE(078)
		ENDIF
		CALL MPI_BCAST(Nr,1,MPI_INTEGER, 0, comm, IERROR)
		IF (me/=0) THEN
			ALLOCATE(zr(Nr),xshift(Nr),yshift(Nr))
		ENDIF
		CALL MPI_BCAST(zr,Nr,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		CALL MPI_BCAST(xshift,Nr,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		CALL MPI_BCAST(yshift,Nr,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		ALLOCATE(recvs(Nr))
		DO I=1,Nr
			recvs(I)%zrecv=zr(I)
			recvs(I)%x_shift=xshift(I)
			recvs(I)%y_shift=yshift(I)
		ENDDO
		DEALLOCATE(zr,xshift,yshift)
	ENDSUBROUTINE
	SUBROUTINE LoadFGMRES_Ctl(fgmres_ctl,comm,fname)
		TYPE (FGMRES_CTL_TYPE),INTENT(OUT)::fgmres_ctl
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::fname
		INTEGER :: buf_length,maxit,maxit_precond,ort_type
		REAL(REALPARM)::misfit
		INTEGER(MPI_CTL_KIND):: me,IERROR
		CALL MPI_COMM_RANK(comm, me, IERROR)
		IF (me==0) THEN
			OPEN(078,file=fname)
			READ(078,*) buf_length,maxit,maxit_precond,misfit, ort_type
			CLOSE(078)
		ENDIF
		CALL MPI_BCAST(buf_length,1,MPI_INTEGER, 0, comm, IERROR)
		CALL MPI_BCAST(ort_type,1,MPI_INTEGER, 0, comm, IERROR)
		CALL MPI_BCAST(maxit,1,MPI_INTEGER, 0, comm, IERROR)
		CALL MPI_BCAST(maxit_precond,1,MPI_INTEGER, 0, comm, IERROR)
		CALL MPI_BCAST(misfit,1,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		fgmres_ctl%misfit=misfit
		fgmres_ctl%fgmres_maxit=maxit
		fgmres_ctl%fgmres_buf=buf_length
		fgmres_ctl%gmres_buf=maxit_precond
		fgmres_ctl%ort_type=ort_type
	ENDSUBROUTINE
	SUBROUTINE LoadThreshold(IE_Threshold,RC_Threshold,comm,fname)
		REAL(REALPARM),INTENT(OUT)::IE_Threshold,RC_Threshold
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::fname
		REAL(REALPARM)::tmp(2)
		INTEGER(MPI_CTL_KIND):: me,IERROR
		CALL MPI_COMM_RANK(comm, me, IERROR)
		IF (me==0) THEN
			OPEN(078,file=fname)
			READ(078,*) tmp
			CLOSE(078)
		ENDIF
		CALL MPI_BCAST(tmp,2,MPI_DOUBLE_PRECISION, 0, comm, IERROR)
		IE_Threshold=tmp(1)
		RC_Threshold=tmp(2)
	ENDSUBROUTINE
	SUBROUTINE SaveField(E,fname,me,t)
		COMPLEX(REALPARM),INTENT(INOUT)::E(:,:,:,:)
		CHARACTER(len=*), INTENT(IN)::fname
		INTEGER,INTENT(IN)::me
		INTEGER, OPTIONAL,INTENT(IN)::t
		INTEGER::t1
		CHARACTER(len=6)::fnum
		WRITE (fnum,'(I5.5)') me
		IF (PRESENT(t)) THEN
			t1=t
		ELSE
			t1=-1
		ENDIF
			IF (t1>0) THEN
				OPEN(UNIT = 11, STATUS='replace',FILE=fname//trim(fnum)//'.dat')
				WRITE (UNIT =11,FMT='(E25.16E3, E25.16E3)') E 
				CLOSE(11) 

				OPEN(UNIT = 11, STATUS='replace',FILE=fname//trim(fnum)//'.bin', form='unformatted')  
				WRITE (UNIT =11) E 
				CLOSE(11) !FMT='(E25.16E3 E25.16E3)'
			ELSEIF (t1==0)THEN
				OPEN(UNIT = 11, STATUS='replace',FILE=fname//trim(fnum)//'.dat')  
				WRITE (UNIT =11,FMT='(E25.16E3, E25.16E3)') E 
				CLOSE(11) 
			ELSE
				OPEN(UNIT = 11, STATUS='replace',FILE=fname//trim(fnum)//'.bin', form='unformatted')  
				WRITE (UNIT =11) E 
				CLOSE(11) !FMT='(E25.16E3 E25.16E3)'
			ENDIF
	END SUBROUTINE
	SUBROUTINE LoadField(E,me,fname)
		COMPLEX(REALPARM),INTENT(INOUT)::E(:,:,:,:)
		CHARACTER(len=*), INTENT(IN)::fname
		INTEGER,INTENT(IN)::me
		CHARACTER(len=6)::fnum
		WRITE (fnum,'(I5.5)') me
			OPEN(UNIT = 11, STATUS='old',FILE=fname//trim(fnum)//'.bin', form='unformatted')  
			READ (UNIT =11) E 
			CLOSE(11) 
	END SUBROUTINE

	SUBROUTINE SaveFieldOneFile(E,fname,t,convert,comm)
		COMPLEX(REALPARM),INTENT(INOUT)::E(:,:,:,:)
		CHARACTER(len=*), INTENT(IN)::fname
		LOGICAL,INTENT(IN)::convert
		INTEGER,INTENT(IN)::t
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		COMPLEX(REALPARM),ALLOCATABLE::Etmp1(:,:,:,:),Etmp2(:,:,:,:)
		CHARACTER(len=6)::fnum
		CHARACTER(len=4)::ftype
		INTEGER::Np,I
		INTEGER::SH(4),Ny
		INTEGER(MPI_CTL_KIND)::me,csize,IERROR
		CALL MPI_COMM_RANK(comm, me, IERROR)
		CALL MPI_COMM_SIZE(comm, csize, IERROR)

		IF (me==0) THEN
		     IF (t==0) THEN
				ftype='.dat'
				OPEN(UNIT = 11,STATUS='replace',FILE=fname//ftype)
		     ELSE
				ftype='.bin'
				OPEN(UNIT = 11,STATUS='replace',FILE=fname//ftype,form='unformatted')
		     ENDIF
		     SH=SHAPE(E)
			 Ny=SH(4)*csize
		     ALLOCATE(Etmp1(SH(1),3,SH(3),Ny),Etmp2(SH(3),Ny,SH(1),3))
		ENDIF
		
		CALL  MPI_GATHER(E, SIZE(E), MPI_DOUBLE_COMPLEX, Etmp1,	SIZE(E), MPI_DOUBLE_COMPLEX,0, comm, IERROR)
		IF (me==0)	THEN
				IF (convert) THEN
					CALL  ConvertInNa(Etmp1,Etmp2)
					IF (t==0) THEN
					       WRITE (UNIT =11,FMT='(E25.16E3, E25.16E3)') Etmp2 
					ELSE
						WRITE(UNIT=11) Etmp2
					ENDIF
				ELSE
					IF (t==0) THEN
					       WRITE (UNIT =11,FMT='(E25.16E3, E25.16E3)') Etmp1 
					ELSE
						WRITE(UNIT=11) Etmp1
					ENDIF
				ENDIF
			CLOSE(11)
		ENDIF	
			
	END SUBROUTINE

#ifdef old_output
	SUBROUTINE SaveOutputOneFile_Old(Ea,Et,Ha,Ht,anomaly,recvs,freq,comm,fname)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Ea(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Et(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Ha(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Ht(:,:,:,:)
		TYPE (RECEIVER_TYPE),POINTER,INTENT(IN)::recvs(:)
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::freq
		INTEGER,INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::fname
		CHARACTER(LEN=*), PARAMETER  :: output_fmt = "(A1, 28ES20.10E3)"
		CHARACTER(LEN=*), PARAMETER  :: title_fmt = "(A1, 28A20)"
		INTEGER::Ir,Nr,me, IERROR
		INTEGER::fsize
		INTEGER::Ix,Iy
		COMPLEX(REALPARM),DIMENSION(1:SIZE(recvs),EX:EZ,1:anomaly%Nx,1:anomaly%Ny)::Ea1,Et1
		COMPLEX(REALPARM),DIMENSION(1:SIZE(recvs),HX:HZ,1:anomaly%Nx,1:anomaly%Ny)::Ha1,Ht1
		REAL(REALPARM)::x,y
		CALL MPI_COMM_RANK(comm, me, IERROR)
		Nr=SIZE(recvs)
		fsize=Nr*3*anomaly%Nx*anomaly%Ny_loc
		CALL  MPI_GATHER(Ea, fsize, MPI_DOUBLE_COMPLEX, Ea1,	fsize, MPI_DOUBLE_COMPLEX,0, comm, IERROR)
		CALL  MPI_GATHER(Et, fsize, MPI_DOUBLE_COMPLEX, Et1,	fsize, MPI_DOUBLE_COMPLEX,0, comm, IERROR)
		CALL  MPI_GATHER(Ha, fsize, MPI_DOUBLE_COMPLEX, Ha1,	fsize, MPI_DOUBLE_COMPLEX,0, comm, IERROR)
		CALL  MPI_GATHER(Ht, fsize, MPI_DOUBLE_COMPLEX, Ht1,	fsize, MPI_DOUBLE_COMPLEX,0, comm, IERROR)
		IF (me/=0) THEN
			CALL MPI_BARRIER(comm,IERROR)
			RETURN
		ENDIF
		OPEN(UNIT = 77,STATUS='replace',FILE=fname//'.dat')
		WRITE (UNIT =77,FMT=title_fmt) "%","x_recv","y_recv","z_recv","frequency",&
											&"Re Ex^a","Im Ex^a", "Re Ey^a", "Im Ey^a",&
											&"Re Ez^a","Im Ez^a", "Re Hx^a", "Im Hx^a",&
											&"Re Hy^a","Im Hy^a", "Re Hz^a", "Im Hz^a",&
											&"Re Ex^t","Im Ex^t", "Re Ey^t", "Im Ey^t",&
											&"Re Ez^t","Im Ez^t", "Re Hx^t", "Im Hx^t",&
											&"Re Hy^t","Im Hy^t", "Re Hz^t", "Im Hz^t"

		DO Ir=1,Nr
			DO Iy=1,anomaly%Ny
				y=Iy*anomaly%dy-anomaly%dy/2d0+recvs(Ir)%y_shift
				DO Ix=1,anomaly%Nx
					x=Ix*anomaly%dx-anomaly%dx/2d0+recvs(Ir)%x_shift
					WRITE (UNIT =77,FMT=output_fmt) " ",x,y,recvs(Ir)%zrecv,freq,&
														& Ea1(Ir,EX,Ix,Iy),Ea1(Ir,EY,Ix,Iy),&
														& Ea1(Ir,EZ,Ix,Iy),Ha1(Ir,HX,Ix,Iy),&
														& Ha1(Ir,HY,Ix,Iy),Ha1(Ir,HZ,Ix,Iy),&
														& Et1(Ir,EX,Ix,Iy),Et1(Ir,EY,Ix,Iy),&
														& Et1(Ir,EZ,Ix,Iy),Ht1(Ir,HX,Ix,Iy),&
														& Ht1(Ir,HY,Ix,Iy),Ht1(Ir,HZ,Ix,Iy)
				ENDDO
			ENDDO
		ENDDO
		CLOSE(77)
		CALL MPI_BARRIER(comm,IERROR)
	ENDSUBROUTINE
#endif
	SUBROUTINE SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freq,comm,fname)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Ea(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Et(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Ha(:,:,:,:)
		COMPLEX(REALPARM),POINTER,INTENT(IN)::Ht(:,:,:,:)
		TYPE (RECEIVER_TYPE),POINTER,INTENT(IN)::recvs(:)
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::freq
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::fname
		CHARACTER(LEN=*), PARAMETER  :: output_fmt = "(A1, 28ES20.10E3)"
		CHARACTER(LEN=*), PARAMETER  :: title_fmt = "(A1, 28A20)"
		INTEGER(MPI_CTL_KIND):: REC_STATUS(MPI_STATUS_SIZE)
		INTEGER(MPI_CTL_KIND)::me, IERROR,csize
		INTEGER::Ir,Nr
		INTEGER::fsize,fshape(4),fsize2
		INTEGER::Ix,Iy,I,Ic
		COMPLEX(REALPARM),POINTER::Ea1(:,:,:,:),Et1(:,:,:,:)
		COMPLEX(REALPARM),POINTER::Ha1(:,:,:,:),Ht1(:,:,:,:)
		REAL(REALPARM)::x,y
		
		CALL MPI_COMM_RANK(comm, me, IERROR)
		CALL MPI_COMM_SIZE(comm, csize, IERROR)
		Nr=SIZE(recvs)
		fsize=SIZE(Ea)
		IF (me==0) THEN
			Ic=0
			fshape=SHAPE(Ea)
			ALLOCATE(Ea1(fshape(1),fshape(2),fshape(3),fshape(4)))
			ALLOCATE(Et1(fshape(1),fshape(2),fshape(3),fshape(4)))
			ALLOCATE(Ha1(fshape(1),fshape(2),fshape(3),fshape(4)))
			ALLOCATE(Ht1(fshape(1),fshape(2),fshape(3),fshape(4)))
			OPEN(UNIT = 77,STATUS='replace',FILE=fname//'.dat')
			WRITE (UNIT =77,FMT=title_fmt) "%","x_recv","y_recv","z_recv","frequency",&
							&"Re Ex^a","Im Ex^a", "Re Ey^a", "Im Ey^a",&
							&"Re Ez^a","Im Ez^a", "Re Hx^a", "Im Hx^a",&
							&"Re Hy^a","Im Hy^a", "Re Hz^a", "Im Hz^a",&
							&"Re Ex^t","Im Ex^t", "Re Ey^t", "Im Ey^t",&
							&"Re Ez^t","Im Ez^t", "Re Hx^t", "Im Hx^t",&
							&"Re Hy^t","Im Hy^t", "Re Hz^t", "Im Hz^t"

			DO Ir=1,Nr
				DO Iy=1,anomaly%Ny_loc
					y=Iy*anomaly%dy-anomaly%dy/2d0+recvs(Ir)%y_shift
						DO Ix=1,anomaly%Nx
							x=Ix*anomaly%dx-anomaly%dx/2d0+recvs(Ir)%x_shift
							WRITE (UNIT =77,FMT=output_fmt) " ",x,y,recvs(Ir)%zrecv,freq,&
											& Ea(Ir,EX,Ix,Iy),Ea(Ir,EY,Ix,Iy),&
											& Ea(Ir,EZ,Ix,Iy),Ha(Ir,HX,Ix,Iy),&
											& Ha(Ir,HY,Ix,Iy),Ha(Ir,HZ,Ix,Iy),&
											& Et(Ir,EX,Ix,Iy),Et(Ir,EY,Ix,Iy),&
											& Et(Ir,EZ,Ix,Iy),Ht(Ir,HX,Ix,Iy),&
											& Ht(Ir,HY,Ix,Iy),Ht(Ir,HZ,Ix,Iy)
							Ic=Ic+1
						ENDDO
				ENDDO
			ENDDO
		ENDIF

		IF (me/=0) THEN
			CALL MPI_SEND(Ea, fsize, MPI_DOUBLE_COMPLEX,0,me,comm, IERROR)
			CALL MPI_SEND(Ha, fsize, MPI_DOUBLE_COMPLEX, 0, me+csize,comm, IERROR)
			CALL MPI_SEND(Et, fsize, MPI_DOUBLE_COMPLEX, 0, me+2*csize,comm, IERROR)
			CALL MPI_SEND(Ht, fsize, MPI_DOUBLE_COMPLEX, 0, me+3*csize,comm, IERROR)
			CALL MPI_BARRIER(comm,IERROR)
			RETURN
		ENDIF
		DO I=1,csize-1
			CALL MPI_RECV(Ea1, fsize, MPI_DOUBLE_COMPLEX, I, I,comm,REC_STATUS, IERROR)
			CALL MPI_RECV(Ha1, fsize, MPI_DOUBLE_COMPLEX, I, I+csize,comm, REC_STATUS,IERROR)
			CALL MPI_RECV(Et1, fsize, MPI_DOUBLE_COMPLEX, I, I+2*csize,comm, REC_STATUS,IERROR)
			CALL MPI_RECV(Ht1, fsize, MPI_DOUBLE_COMPLEX, I, I+3*csize,comm, REC_STATUS,IERROR)
			DO Ir=1,Nr
				DO Iy=1,anomaly%Ny_loc
					y=(Iy+I*anomaly%Ny_loc)*anomaly%dy-anomaly%dy/2d0+recvs(Ir)%y_shift
					DO Ix=1,anomaly%Nx
						x=Ix*anomaly%dx-anomaly%dx/2d0+recvs(Ir)%x_shift
						WRITE (UNIT =77,FMT=output_fmt) " ",x,y,recvs(Ir)%zrecv,freq,&
										& Ea1(Ir,EX,Ix,Iy),Ea1(Ir,EY,Ix,Iy),&
										& Ea1(Ir,EZ,Ix,Iy),Ha1(Ir,HX,Ix,Iy),&
										& Ha1(Ir,HY,Ix,Iy),Ha1(Ir,HZ,Ix,Iy),&
										& Et1(Ir,EX,Ix,Iy),Et1(Ir,EY,Ix,Iy),&
										& Et1(Ir,EZ,Ix,Iy),Ht1(Ir,HX,Ix,Iy),&
										& Ht1(Ir,HY,Ix,Iy),Ht1(Ir,HZ,Ix,Iy)
						Ic=Ic+1
					ENDDO
				ENDDO
			ENDDO
		ENDDO
		CLOSE(77)
		DEALLOCATE(Ea1,Ha1,Et1,Ht1)
		CALL MPI_BARRIER(comm,IERROR)
	ENDSUBROUTINE
	SUBROUTINE SaveOutputSeparate(Ea,Et,Ha,Ht,anomaly,recvs,freq,tag,offset,fname)
		COMPLEX(REALPARM),INTENT(IN)::Ea(1:,EX:,1:,1:)
		COMPLEX(REALPARM),INTENT(IN)::Et(1:,EX:,1:,1:)
		COMPLEX(REALPARM),INTENT(IN)::Ha(1:,HX:,1:,1:)
		COMPLEX(REALPARM),INTENT(IN)::Ht(1:,HX:,1:,1:)
		INTEGER,INTENT(IN)::offset
		TYPE (RECEIVER_TYPE),POINTER,INTENT(IN)::recvs(:)
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::freq
		INTEGER,INTENT(IN)::tag
		CHARACTER(*),INTENT(IN)::fname
		CHARACTER(LEN=*), PARAMETER  :: output_fmt = "(A1, 28ES20.10E3)"
		CHARACTER(LEN=*), PARAMETER  :: title_fmt = "(A1, 28A20)"
		CHARACTER(LEN=6):: ftag
		INTEGER::Ir,Nr,me, IERROR
		INTEGER::fsize
		INTEGER::Ix,Iy,Nyloc,s(4)
		REAL(REALPARM)::x,y
		Nr=SIZE(recvs)
		s=SHAPE(Ea)
		Nyloc=s(4)
		WRITE (ftag,'(I5.5)') tag
		OPEN(UNIT = 277,STATUS='replace',FILE=fname//trim(ftag)//'.dat')
		WRITE (UNIT =277,FMT=title_fmt) "%","x_recv","y_recv","z_recv","frequency",&
											&"Re Ex^a","Im Ex^a", "Re Ey^a", "Im Ey^a",&
											&"Re Ez^a","Im Ez^a", "Re Hx^a", "Im Hx^a",&
											&"Re Hy^a","Im Hy^a", "Re Hz^a", "Im Hz^a",&
											&"Re Ex^t","Im Ex^t", "Re Ey^t", "Im Ey^t",&
											&"Re Ez^t","Im Ez^t", "Re Hx^t", "Im Hx^t",&
											&"Re Hy^t","Im Hy^t", "Re Hz^t", "Im Hz^t"

		DO Ir=1,Nr
			DO Iy=1,Nyloc
				y=(Iy+offset)*anomaly%dy-anomaly%dy/2d0+recvs(Ir)%y_shift
				DO Ix=1,anomaly%Nx
					x=Ix*anomaly%dx-anomaly%dx/2d0+recvs(Ir)%x_shift
					WRITE (UNIT =277,FMT=output_fmt) " ",x,y,recvs(Ir)%zrecv,freq,&
														& Ea(Ir,EX,Ix,Iy),Ea(Ir,EY,Ix,Iy),&
														& Ea(Ir,EZ,Ix,Iy),Ha(Ir,HX,Ix,Iy),&
														& Ha(Ir,HY,Ix,Iy),Ha(Ir,HZ,Ix,Iy),&
														& Et(Ir,EX,Ix,Iy),Et(Ir,EY,Ix,Iy),&
														& Et(Ir,EZ,Ix,Iy),Ht(Ir,HX,Ix,Iy),&
														& Ht(Ir,HY,Ix,Iy),Ht(Ir,HZ,Ix,Iy)
				ENDDO
			ENDDO
		ENDDO
		CLOSE(277)
	ENDSUBROUTINE
	SUBROUTINE SaveIESolutionSeparate(Eint,anomaly,freq,tag,offset,fname)
		COMPLEX(REALPARM),INTENT(IN)::Eint(1:,EX:,1:,1:)
		INTEGER,INTENT(IN)::offset
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::freq
		INTEGER,INTENT(IN)::tag
		CHARACTER(*),INTENT(IN)::fname
		CHARACTER(LEN=*), PARAMETER  :: output_fmt = "(A1, 10ES20.10E3)"
		CHARACTER(LEN=*), PARAMETER  :: title_fmt = "(A1, 10A20)"
		CHARACTER(LEN=6):: ftag
		INTEGER::Ir,Nx,Nz,me, IERROR
		INTEGER::fsize
		INTEGER::Ix,Iy,Iz,Nyloc,s(4)
		REAL(REALPARM)::x,y,z
		s=SHAPE(Eint)
		Nyloc=s(4)
		WRITE (ftag,'(I5.5)') tag
		OPEN(UNIT = 277,STATUS='replace',FILE=fname//trim(ftag)//'.dat')
		WRITE (UNIT =277,FMT=title_fmt) "%","x_c","y_c","z_c","frequency",&
						&"Re Ex","Im Ex", "Re Ey", "Im Ey",&
						&"Re Ez","Im Ez"
		Nx=anomaly%Nx
		DO Iy=1,Nyloc
			y=(Iy+offset)*anomaly%dy-anomaly%dy/2d0
			DO Ix=1,Nx
				x=Ix*anomaly%dx-anomaly%dx/2d0
				DO Iz=1,anomaly%Nz
				z=(anomaly%z(Iz)+anomaly%z(Iz-1))/2;
				WRITE (UNIT =277,FMT=output_fmt) " ",x,y,z,freq,&
								& Eint(Iz,EX,Ix,Iy),Eint(Iz,EY,Ix,Iy),&
								& Eint(Iz,EZ,Ix,Iy)
				ENDDO
			ENDDO
		ENDDO
		CLOSE(277)
	ENDSUBROUTINE
	SUBROUTINE SaveIESolutionSeparateBinary(Eint,tag,fname)
		COMPLEX(REALPARM),INTENT(IN)::Eint(1:,EX:,1:,1:)
		INTEGER,INTENT(IN)::tag
		CHARACTER(*),INTENT(IN)::fname
		CHARACTER(LEN=6):: ftag
		INTEGER::Ir,Nx,Nz,me, IERROR
		INTEGER::fsize
		INTEGER::Ix,Iy,Iz,Nyloc,s(4)
		REAL(REALPARM)::x,y,z
		WRITE (ftag,'(I5.5)') tag
!		OPEN(UNIT = 277,STATUS='replace',FILE=fname//trim(ftag)//'.bin',form='binary',access="stream")
		OPEN(UNIT = 277,STATUS='replace',FILE=fname//trim(ftag)//'.bin',form='unformatted',access="stream")
		WRITE (UNIT =277) Eint
		CLOSE(277)
	ENDSUBROUTINE
	SUBROUTINE SaveIESolutionOneFileBinary(Eint,comm,fname)
		COMPLEX(REALPARM),INTENT(IN)::Eint(1:,EX:,1:,1:)
		INTEGER(MPI_CTL_KIND),INTENT(IN)::comm
		CHARACTER(*),INTENT(IN)::fname
		INTEGER(MPI_CTL_KIND)::me, IERROR,csize
		INTEGER::fsize,fshape(4)
		INTEGER(MPI_CTL_KIND):: REC_STATUS(MPI_STATUS_SIZE)
		INTEGER::Ix,Iy,I
		COMPLEX(REALPARM),POINTER::Eint1(:,:,:,:)
		REAL(REALPARM)::x,y
		
		CALL MPI_COMM_RANK(comm, me, IERROR)
		CALL MPI_COMM_SIZE(comm, csize, IERROR)
		fsize=SIZE(Eint)
		IF (me==0) THEN
			fshape=SHAPE(Eint)
			ALLOCATE(Eint1(fshape(1),fshape(2),fshape(3),fshape(4)))
!			OPEN(UNIT = 77,STATUS='replace',FILE=fname//'.bin',form='binary',access="stream")
			OPEN(UNIT = 77,STATUS='replace',FILE=fname//'.bin',form='unformatted',access="stream")
			WRITE (UNIT =77) Eint
		ENDIF

		IF (me/=0) THEN
			CALL MPI_SEND(Eint, fsize, MPI_DOUBLE_COMPLEX,0,me,comm, IERROR)
			CALL MPI_BARRIER(comm,IERROR)
			RETURN
		ENDIF
		DO I=1,csize-1
			CALL MPI_RECV(Eint1, fsize, MPI_DOUBLE_COMPLEX, I, I,comm, REC_STATUS,IERROR)
			WRITE (UNIT =77) Eint1
		ENDDO
		CLOSE(77)
		DEALLOCATE(Eint1)
		CALL MPI_BARRIER(comm,IERROR)
	ENDSUBROUTINE
	SUBROUTINE LoadIESolutionSeparateBinary(Eint,tag,fname)
		COMPLEX(REALPARM),INTENT(INOUT)::Eint(1:,EX:,1:,1:)
		INTEGER,INTENT(IN)::tag
		CHARACTER(*),INTENT(IN)::fname
		CHARACTER(LEN=6):: ftag
		INTEGER::Ir,Nx,Nz,me, IERROR
		INTEGER::fsize
		INTEGER::Ix,Iy,Iz,Nyloc,s(4)
		REAL(REALPARM)::x,y,z
		WRITE (ftag,'(I5.5)') tag
!		OPEN(UNIT = 277,STATUS='old',FILE=fname//trim(ftag)//'.bin',form='binary',access="stream")
		OPEN(UNIT = 277,STATUS='old',FILE=fname//trim(ftag)//'.bin',form='unformatted',access="stream")
		READ (UNIT =277) Eint
		CLOSE(277)
	ENDSUBROUTINE
	SUBROUTINE ConvertNaIn(E1,E2) !Convert EM field array  E1 from "naive" format E1(x,y,z,c) to internal format E2(z,c,x,y)
		COMPLEX(REALPARM),INTENT(IN)::E1(:,:,:,:)
		COMPLEX(REALPARM),INTENT(OUT)::E2(:,:,:,:)
		INTEGER ::S1(4),S2(4),Iz,Ic
		S1=SHAPE(E1)
		S2=SHAPE(E2)
		IF (S1(3)==S2(1)) THEN
			!$OMP PARALLEL PRIVATE(Iz,Ic) DEFAULT(SHARED)
        			!$OMP DO  SCHEDULE (GUIDED)!COLLAPSE(2) 
				DO Ic=EX,EZ
					DO Iz=1,S2(1)
						E2(Iz,Ic,:,:)=E1(:,:,Iz,Ic)
					ENDDO
				ENDDO
	       		!$OMP END DO
			!$OMP END PARALLEL
		ELSE
			PRINT*, 'CONVERT ERROR'
		ENDIF
	
	END SUBROUTINE
	SUBROUTINE ConvertInNa(E1,E2) !Convert EM field array E1 from internal format E1(z,c,x,y) to "naive" E2(x,y,z,c)
		COMPLEX(REALPARM),INTENT(IN)::E1(:,:,:,:)
		COMPLEX(REALPARM),INTENT(OUT)::E2(:,:,:,:)
		INTEGER ::S1(4),S2(4),Ic,Iz
		S1=SHAPE(E1)
		S2=SHAPE(E2)
		IF (S1(1)==S2(3)) THEN
			!$OMP PARALLEL PRIVATE(Iz,Ic) DEFAULT(SHARED)
		        	!$OMP DO   SCHEDULE (GUIDED) !COLLAPSE(2)
				DO Ic=EX,EZ
					DO Iz=1,S1(1)
						E2(:,:,Iz,Ic)=E1(Iz,Ic,:,:)
					ENDDO
				ENDDO
			        !$OMP END DO
			!$OMP END PARALLEL
		ELSE
			PRINT*, 'CONVERT ERROR'
		ENDIF
	
	END SUBROUTINE
END MODULE
