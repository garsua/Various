!--------------------------------------------------------------------

PROGRAM row2column

! -------------------------------------------------------------------
! This program converts data written in a row to data in a column
! -------------------------------------------------------------------
! Coded by Victor M. Garcia-Suarez
! Universidad de Oviedo & CINN
! e-mail: garciavictor@uniovi.es
! Date:	October 2015
!--------------------------------------------------------------------

implicit none

integer:: i,j,ierr
!double precision:: ist,ifi
double precision:: vcd
logical:: fileexist

!call cpu_time(ist)

inquire(file='row.csv',exist=fileexist)
if(.not.fileexist) then
  write(6,*)'ERROR: There is no row.csv'
  stop
endif
open(12,file='row.csv',status='unknown')
open(13,file='col.csv',status='unknown')
ierr=0
j=1
do while(ierr.eq.0)
  read(12,*,iostat=ierr)(vcd,i=1,j)
  j=j+1
  rewind(12)
  if(ierr.eq.0)write(13,*)vcd
enddo
close(13)
close(12)

!call cpu_time(ifi)
!write(6,'(f12.8)')ifi-ist

end program row2column

!--------------------------------------------------------------------
