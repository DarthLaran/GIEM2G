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

MODULE CHECK_MEMORY
	USE MPI_MODULE
	USE LOGGER_MODULE
	 IMPLICIT NONE
         REAL,PARAMETER::PAGE_SIZE=4.0/1024/1024
CONTAINS
	SUBROUTINE CHECK_MEM(me,master,comm)
		INTEGER(MPI_CTL_KIND),INTENT(IN)::me,master,comm
#ifndef NO_MEMORY_CHECKING
		CHARACTER(LEN=20) :: FILE
		INTEGER :: SIZE=0, RESIDENT=0, SHARE=0, TEXT=0, LIB=0, DATA_LOC=0, DT=0
		INTEGER:: MAX_USAGE,TOTAL_USAGE
		INTEGER:: IERR=0
		CHARACTER(LEN=2048)::message 

		FILE='/proc/self/statm'
		OPEN(UNIT=1,FILE=FILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
		READ(UNIT=1,FMT=*,IOSTAT=IERR) SIZE, RESIDENT, SHARE, TEXT, LIB, DATA_LOC, DT
		CLOSE(UNIT=1)

		IF (IERR < 0) THEN
			WRITE(*,*)'PROBLEM READING /PROC/SELF/STATM'
		ENDIF
		CALL MPI_REDUCE(DATA_LOC,MAX_USAGE,1,MPI_INT,MPI_MAX,master,comm,IERR)
		CALL MPI_REDUCE(DATA_LOC,TOTAL_USAGE,1,MPI_INT,MPI_SUM,master,comm,IERR)
		WRITE (message,'(A, 1ES10.2E2, A)')  "Total memory usage:", TOTAL_USAGE*PAGE_SIZE, " GB"
	        CALL LOGGER(TRIM(message)) 
		WRITE (message,'(A, 1ES10.2E2, A)')  "Maximum memory usage at one node:", MAX_USAGE*PAGE_SIZE, " GB"
	        CALL LOGGER(TRIM(message)) 
#endif
	END SUBROUTINE 
END
