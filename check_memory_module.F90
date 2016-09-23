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
         REAL,PARAMETER::KB2GB=1.0/1024/1024
CONTAINS
	SUBROUTINE CHECK_MEM(me,master,comm)
		INTEGER(MPI_CTL_KIND),INTENT(IN)::me,master,comm
#ifndef NO_MEMORY_CHECKING
		CHARACTER(LEN=20) :: FILE
		INTEGER :: SIZE=0, RESIDENT=0, SHARE=0, TEXT=0, LIB=0, DATA_LOC=0, DT=0
		INTEGER:: MAX_USAGE,TOTAL_USAGE
		INTEGER:: IERR=0
		CHARACTER(LEN=2048)::message 
                
                DATA_LOC=OBTAIN_PHYSAL_MEMORY()
                
		CALL MPI_REDUCE(DATA_LOC,MAX_USAGE,1,MPI_INT,MPI_MAX,master,comm,IERR)
		CALL MPI_REDUCE(DATA_LOC,TOTAL_USAGE,1,MPI_INT,MPI_SUM,master,comm,IERR)
		WRITE (message,'(A, 1ES10.2E2, A)')  "Total physical memory usage:", TOTAL_USAGE*KB2GB, " GB"
	        CALL LOGGER(TRIM(message)) 
		WRITE (message,'(A, 1ES10.2E2, A)')  "Maximum physical memory usage at one node:", MAX_USAGE*KB2GB, " GB"
	        CALL LOGGER(TRIM(message)) 
#endif
	END SUBROUTINE


        FUNCTION OBTAIN_PHYSAL_MEMORY() RESULT(MEM)!RESULT IS IN KB!!
                INTEGER::MEM
                CHARACTER(LEN=256)::line
                CHARACTER(LEN=*),PARAMETER::stat_file='/proc/self/status'
                INTEGER::IERR,FD
                OPEN(UNIT=FD,FILE=stat_file,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
		IF (IERR < 0) THEN
			WRITE(*,*)'PROBLEM READING /proc/self/status !!!'
		ENDIF
                DO
                        READ(UNIT=FD,FMT='(A)') line
                        IF (line(1:5)=="VmRSS") THEN
                                READ(line(7:),*) MEM
                                EXIT
                        ENDIF
                ENDDO
                CLOSE(FD)
        ENDFUNCTION
END
