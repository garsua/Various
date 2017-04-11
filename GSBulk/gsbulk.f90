!--------------------------------------------------------------------

PROGRAM GSBulk

! -------------------------------------------------------------------
! This program calculates the conductance, Seebeck coefficiente and
! ofther quantities of a bulk material usin the Boltzmann T. E.
! -------------------------------------------------------------------
! Coded by Victor M. Garcia-Suarez
! Universidad de Oviedo & CINN
! e-mail: garciavictor@uniovi.es
! Date:	September 2015
!--------------------------------------------------------------------

use ioabc_m

implicit none

!call cpu_time(ist)

call io_data()

if(abcode.eq.'qe')then
  call qe_data()
else if(abcode.eq.'vasp')then
  call vasp_data()
else
  write(6,'(a)')'  ERROR: Unknown ab-initio code'
  stop
endif

!call checks()

! Print alat and # of atms and sps -
if(rgeom)then
  call get_unit_number(ugm)
  open(ugm,file='geo.csv',status='unknown')
  write(ugm,'(i8,a,i8,a,f14.8)')nat,',',nspec,',',alat
  close(ugm)
endif

! Calculate density of states ------
if(doscalc)then
  call get_unit_number(udos)
  open(udos,file='dos.dat',status='unknown')
  nel=0.d0
  if(fermiref)refe=ef
  en=emin
  allocate(dos(nspin))
  do while(en.le.emax)
    dos=0.d0
    do i=1,nk
      do j=1,nek
        do ispin=1,nspin
          dos(ispin)=dos(ispin)+kcart(4,i)/((en-enk(j,i,ispin))**2.d0+&
            delta**2.d0)
        enddo
      enddo
    enddo
    dos=dos*sdeg*delta/pi
    write(udos,*)en-refe,(dos(ispin),ispin=1,nspin)
    !nel=nel+dos*estep*fermif(en-ef,temp)
    en=en+estep
  enddo
  deallocate(dos)
  close(udos)
endif

! Find highest ev and lowest ec ----
if(gapcalc.or.evecalc)then
  hv=-1.d6
  lc=1.d6
  do i=1,nk
    do ispin=1,nspin
      j=1
      do while(enk(j,i,ispin)-ef.lt.0.d0)
        j=j+1
      enddo
      if(j.gt.1)then
        if(enk(j-1,i,ispin).gt.hv)hv=enk(j-1,i,ispin)
        if(enk(j,i,ispin).lt.lc)lc=enk(j,i,ispin)
      endif
    enddo
  enddo
endif

! Calculate and print gap ----------
if(gapcalc)then
  gap=lc-hv
  if(gap.le.0.1d0)then
    aa='M'
  else if(gap.le.2.d0)then
    aa='S'
  else if(gap.le.5.d0)then
    aa='I'
  else
    aa='L'
  endif
  call get_unit_number(ugp)
  open(ugp,file='gap.csv',status='unknown')
  write(ugp,'(f14.8,a,3x,a)')gap,',',aa
  close(ugp)
endif

! Calculate DOS vector -------------
if(evecalc)then
  newef=(lc+hv)/2.d0
  refvmin=newef+evecmin
  refvmax=newef+evecmax
  estep=(refvmax-refvmin)/real(evecpts)
  allocate(vecene(evecpts))
  vecene=0.d0
  do i=1,nk
    do ispin=1,nspin
      eprev=refvmin
      en=refvmin+estep
      j=1
      l=1
      do while(en.le.refvmax)
        do while(enk(l,i,ispin).lt.en.and.l.le.nek)
          if(enk(l,i,ispin).ge.eprev)&
            vecene(j)=vecene(j)+kcart(4,i)
          l=l+1
        enddo
        eprev=en
        en=en+estep
        j=j+1
      enddo
    enddo
  enddo
  !nel=0
  !do i=1,evecpts
  !  nel=nel+vecene(i)
  !  print*,i,vecene(i)
  !enddo
  !print*,nel
  call get_unit_number(uve)
  open(uve,file='eve.csv',status='unknown')
  write(uve,'(f12.4)',advance='no')vecene(1)
  do i=2,evecpts
    write(uve,'(a,f12.4)',advance='no')',',vecene(i)
  enddo
  deallocate(vecene)
  close(uve)
endif

if(allocated(enk))deallocate(enk)
if(allocated(nev))deallocate(nev)
if(allocated(kcart))deallocate(kcart)
if(allocated(spec))deallocate(spec)
if(allocated(coord))deallocate(coord)
if(allocated(revect))deallocate(revect)
if(allocated(vect))deallocate(vect)

!call cpu_time(ifi)
!write(6,'(f12.8)')ifi-ist

end program GSBulk

!--------------------------------------------------------------------
