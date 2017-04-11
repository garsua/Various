!--------------------------------------------------------------------

module ioabc_m

implicit none
integer,parameter           :: dp=selected_real_kind(14,100)
integer                     :: countarg,ierr,where,uin,uwout,ugm,udos,ugp,uve,&
                               nline
integer                     :: i,j,l,ispin,m,kss,nat,nspec,nspin,nk,nrow,&
                               ncol,ninline,nek,evecpts
integer,allocatable         :: nev(:),iopt(:)
real(dp),parameter          :: pi=3.14159265d0
real(dp),parameter          :: a0=0.52917721d0
!real(dp)                    :: ist,ifi
real(dp)                    :: emin,emax,estep,delta,refe,alat,ucvol,eln,twk,&
                               ee,ef,temp,fermif,nel,te,en,sdeg,hv,lc,gap,&
                               newef,evecmin,evecmax,refvmin,refvmax,eprev
real(dp),allocatable        :: vect(:,:),revect(:,:),auxv(:),coord(:,:),&
                               kcart(:,:),kfrac(:,:),enk(:,:,:),dos(:),&
                               vecene(:)
character(len=1)            :: aa
character(len=2)            :: efref
character(len=10)           :: abcode
character(len=50)           :: dmy,fwout,logfile
character(len=85)           :: cvalue,intpar
character(len=256)          :: line
character(len=3),allocatable:: spec(:),sp1(:)
logical                     :: fileexist
logical                     :: rtestdat,fermiref,rgeom,quitl,doscalc,&
                               gapcalc,evecalc
external                    :: intpar,ninline,fermif
!external                    :: cart_to_frac,frac_to_cart

contains

!--------------------------------------------------------------------

subroutine io_data()

! Read from the input arguments -----
countarg = command_argument_count()
if (countarg<1) then
  write(6,'(a)')'  ERROR: Too few arguments'
  write(6,'(a)')'  USAGE: gsbulk scf-out_file'
  stop
endif
call getarg(1,fwout)

! Read data from GSBulk input file --
inquire(file='gsb.in',exist=fileexist)
if(.not.fileexist)then
  write(6,'(a)')'  WARNING: There is no gsb.in file'
  write(6,'(a)')'  Default ab-initio code: QE'
  abcode='qe'
  write(6,'(a)')'  Default value of Emin = -20.0 eV'
  emin=-20.d0
  write(6,'(a)')'  Default value of Emax = 20.0 eV'
  emax=20.d0
  write(6,'(a)')'  Default value of Estep = 0.1 eV'
  estep=0.1d0
  write(6,'(a)')'  Default value of delta = 0.1 eV'
  delta=0.1d0
  write(6,'(a)')'  Default value of energy reference = 0.0 eV'
  refe=0.d0
  write(6,'(a)')'  Default value of temperature = 0.1 K'
  temp=0.1
endif

if(fileexist)then

  call get_unit_number(uin)
  open(uin,file='gsb.in',status='old')

  call read_in(uin,'ab code',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  Default ab-initio code: QE'
    abcode='qe'
  else
    read(cvalue,*)abcode
  endif

  call read_in(uin,'read testdat',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  No test data'
    rtestdat=.false.
  else
    read(cvalue,*)rtestdat
  endif

  call read_in(uin,'dos calc',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  DOS not calculated by default'
    doscalc=.false.
  else
    read(cvalue,*)doscalc
  endif

  if(doscalc)then

    call read_in(uin,'min energy',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of Emin = -20.0 eV'
      emin=-20.d0
    else
      read(cvalue,*)emin
    endif
  
    call read_in(uin,'max energy',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of Emax = 20.0 eV'
      emax=20.d0
    else
      read(cvalue,*)emax
    endif
  
    call read_in(uin,'energy step',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of Estep = 0.1 eV'
      estep=0.1d0
    else
      read(cvalue,*)estep
    endif
  
    call read_in(uin,'delta',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of delta = 0.1 eV'
      delta=0.1d0
    else
      read(cvalue,*)delta
    endif
  
    fermiref=.false.
    call read_in(uin,'ref energy',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of energy reference = 0.0 eV'
      refe=0.d0
    else
      read(cvalue,*,iostat=ierr)refe
      if(ierr.ne.0)then
        read(cvalue,*)efref
        if(efref.eq.'ef')then
          fermiref=.true.
        else
          write(6,'(a)')'  ERROR:  Incorrect value of the Fermi energy'
          stop
        endif
      endif
    endif

  endif

  call read_in(uin,'evec calc',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  Energy vector not calculated by default'
    evecalc=.false.
  else
    read(cvalue,*)evecalc
  endif

  if(evecalc)then

    call read_in(uin,'min vecenergy',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of VecEmin = -10.0 eV'
      evecmin=-10.d0
    else
      read(cvalue,*)evecmin
    endif
  
    call read_in(uin,'max vecenergy',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of VecEmax = 10.0 eV'
      evecmax=10.d0
    else
      read(cvalue,*)evecmax
    endif
  
    call read_in(uin,'evec points',cvalue)
    if(cvalue.eq.'')then
      write(6,'(a)')'  Default value of Vec points = 32'
      evecpts=32
    else
      read(cvalue,*)evecpts
    endif
    
  endif

  call read_in(uin,'temperature',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  Default value of temperature = 0.1 K'
    temp=0.1
  else
    read(cvalue,*)temp
  endif
 
  call read_in(uin,'read geometry',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  Geometry not read'
    rgeom=.false.
  else
    read(cvalue,*)rgeom
  endif

  call read_in(uin,'calculate gap',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  Gap not calculated'
    gapcalc=.false.
  else
    read(cvalue,*)gapcalc
  endif

  close(uin)

endif

! Open out file ---------------------
logfile=trim(fwout)
inquire(file=logfile,exist=fileexist)
if(.not.fileexist)then
  write(6,'(2a)')'  ERROR: There is no ',logfile
  stop
endif
call get_unit_number(uwout)
open(uwout,file=logfile,status='old')

end subroutine io_data

!--------------------------------------------------------------------

subroutine qe_data()

! Read data from QE out file --------

if(rtestdat)then

  call read_out(uwout,'unit-cell volume          =',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  WARNING: No value of U. C. volume in the output file'
  else
    read(cvalue,*)ucvol
  endif

  call read_out(uwout,'number of electrons       =',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  WARNING: No value of number of e- in the output file'
  else
    read(cvalue,*)eln
  endif

endif

call read_out(uwout,'number of Kohn-Sham states=',cvalue)
if(cvalue.eq.'')then
  write(6,'(a)')'  WARNING: No value of KS states (bands) in the output file'
else
  read(cvalue,*)kss
endif

if(rgeom)then

  call read_out(uwout,'lattice parameter (alat)  =',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  ERROR: No value of alat in the output file'
    stop
  else
    read(cvalue,*)alat
  endif
  alat=alat*a0

  call read_out(uwout,'number of atoms/cell      =',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  ERROR: No value of # of atoms in the  output file'
    stop
  else
    read(cvalue,*)nat
  endif

  call read_out(uwout,'number of atomic types    =',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  WARNING: No value of # of elements in the output file'
  else
    read(cvalue,*)nspec
  endif

endif

! Read spin components -------------
rewind(uwout)
where=0
ierr=0
do while(where.eq.0.and.ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'SPIN')
enddo
if(ierr.ne.0)then
  nspin=1
else
  nspin=2
endif

! Read Fermi energy and Temp. ------
call read_out(uwout,'the Fermi energy is',cvalue)
if(cvalue.eq.'')then
  write(6,'(a)')'  ERROR: No value of Fermi energy in the output file'
  write(6,'(a)')'  USAGE: You should probably include occupations and degauss'
  stop
else
  read(cvalue,*)ef
endif

if(rgeom)then

! Read vectors ----------------------
  rewind(uwout)
  where=0
  ierr=0
  do while(where.eq.0.and.ierr.eq.0)
    read(uwout,'(a)',iostat=ierr)line
    where=index(line,'crystal axes')
  enddo
  if(ierr.ne.0)then
    write(6,'(a)')'  ERROR: No lattice vectors in the ouput file'
    stop
  endif
  allocate(vect(3,3))
  do i=1,3
    read(uwout,'(a)')line
    where=index(line,'=')
    line=line(where:len(line))
    cvalue=intpar(line)
    read(cvalue,*)(vect(j,i),j=1,3)
  enddo
  rewind(uwout)
  where=0
  ierr=0
  do while(where.eq.0.and.ierr.eq.0)
    read(uwout,'(a)',iostat=ierr)line
    where=index(line,'reciprocal axes')
  enddo
  if(ierr.ne.0)then
    write(6,'(a)')'  ERROR: No reciprocal lattice vectors in the ouput file'
    stop
  endif
  allocate(revect(3,3))
  do i=1,3
    read(uwout,'(a)')line
    where=index(line,'=')
    line=line(where:len(line))
    cvalue=intpar(line)
    read(cvalue,*)(revect(j,i),j=1,3)
  enddo

! Read coordinates ------------------
  rewind(uwout)
  where=0
  ierr=0
  do while(where.eq.0.and.ierr.eq.0)
    read(uwout,'(a)',iostat=ierr)line
    where=index(line,'positions (alat units)')
  enddo
  if(ierr.ne.0)then
    write(6,'(a)')'  ERROR: No coordinates in the ouput file'
    stop
  endif
  allocate(spec(nat),coord(3,nat))
  do i=1,nat
    read(uwout,*)l,spec(i)
    backspace(uwout)
    read(uwout,'(a)')line
    where=index(line,'=')
    line=line(where:len(line))
    cvalue=intpar(line)
    read(cvalue,*)(coord(j,i),j=1,3)
  enddo
  coord=coord*alat*a0

endif

! Read k points ---------------------
rewind(uwout)
where=0
ierr=0
do while(where.eq.0.and.ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'number of k points=')
enddo
if(ierr.ne.0)then
  write(6,'(a)')'  ERROR: No k points in the ouput file'
  stop
endif
where=index(line,'=')
read(line(where+1:len(line)),*)nk
allocate(kcart(4,nk))

where=0
ierr=0
do while(where.eq.0.and.ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'cart. coord.')
enddo
if(ierr.ne.0)then
  write(6,'(a)')'  ERROR: No cart. coord. k points in the ouput file'
  stop
endif
do i=1,nk
  read(uwout,'(a)')line
  where=index(line,'=')
  cvalue=intpar(line(where+1:len(line)))
  read(cvalue,*,iostat=ierr)(kcart(j,i),j=1,3)
  where=index(line,'wk =')
  read(line(where+4:len(line)),*)kcart(4,i)
enddo

! Read eigenvalues ------------------
rewind(uwout)
where=0
ierr=0
nline=1
do while(where.eq.0.and.ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'End of')
  nline=nline+1
enddo
if(ierr.ne.0)then
  write(6,'(a)')'  ERROR: The self-consistent calc. has not finished'
  stop
endif
where=0
do while(where.eq.0)
  read(uwout,'(a)')line
  where=index(line,'k =')
enddo
read(uwout,'(a)')dmy
ierr=0
nrow=-1
do while(ierr.eq.0)
  read(uwout,*,iostat=ierr)ee
  nrow=nrow+1
enddo
if(nrow.eq.0)then
  write(6,'(a)')'  ERROR: No eigenvalues printed in the output file'
  stop
endif
allocate(nev(nrow))

rewind(uwout)
do i=1,nline
  read(uwout,'(a)')dmy
enddo
where=0
do while(where.eq.0)
  read(uwout,'(a)')line
  where=index(line,'k =')
enddo
read(uwout,'(a)')dmy
do i=1,nrow
  read(uwout,'(a)')line
  nev(i)=ninline(line)
enddo
nek=sum(nev)

allocate(enk(nek,nk,nspin))

where=0
ierr=0
rewind(uwout)
do i=1,nline
  read(uwout,'(a)')dmy
enddo
l=1
ispin=1
do while(ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'k =')
  if(where.ne.0)then
    m=1
    do i=1,nrow
      read(uwout,*)(enk(m+j-1,l,ispin),j=1,nev(i))
      m=m+nev(i)
    enddo
    where=0
    if(l.eq.nk/2.and.nspin.eq.2)ispin=2
    l=l+1
  endif
enddo

! Spin degeneracy -/-----------------
sdeg=1.d0

end subroutine qe_data

!--------------------------------------------------------------------

subroutine vasp_data()

! Read data from VASP out file ------

if(rtestdat)then

  call read_out(uwout,'volume of cell',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  WARNING: No value of U. C. volume in the output file'
  else
    read(cvalue,*)ucvol
  endif

  call read_out(uwout,'NELECT',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  WARNING: No value of number of e- in the output file'
  else
    read(cvalue,*)eln
  endif

endif

call read_out(uwout,'NBANDS',cvalue)
if(cvalue.eq.'')then
  write(6,'(a)')'  WARNING: No value of KS states (bands) in the output file'
else
  read(cvalue,*)kss
endif

if(rgeom)then

  call read_out(uwout,'ALAT',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  ERROR: No value of alat in the output file'
    stop
  else
    read(cvalue,*)alat
  endif

  call read_out(uwout,'NIONS',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  ERROR: No value of # of atoms in the  output file'
    stop
  else
    read(cvalue,*)nat
  endif

  rewind(uwout)
  nspec=0
  where=0
  ierr=0
  do while(where.eq.0.and.ierr.eq.0)
    read(uwout,'(a)',iostat=ierr)line
    where=index(line,'VRHFIN')
    if(where.ne.0)then
      nspec=nspec+1
      where=0
    endif
  enddo

endif

! Read spin components -------------
call read_out(uwout,'ISPIN',cvalue)
if(cvalue.eq.'')then
  write(6,'(a)')'  ERROR: No value of spin components in the output file'
  stop
else
  read(cvalue,*)nspin
endif

! Read Fermi energy and Temp. ------
call read_out(uwout,'E-fermi',cvalue)
if(cvalue.eq.'')then
  write(6,'(a)')'  ERROR: No value of Fermi energy in the output file'
  stop
else
  read(cvalue,*)ef
endif

if(rgeom)then

! Read vectors ----------------------
  rewind(uwout)
  where=0
  ierr=0
  do while(where.eq.0.and.ierr.eq.0)
    read(uwout,'(a)',iostat=ierr)line
    where=index(line,'direct lattice vectors')
  enddo
  if(ierr.ne.0)then
    write(6,'(a)')'  ERROR: No lattice or reciprocal vectors in the out file'
    stop
  endif
  allocate(vect(3,3))
  allocate(revect(3,3))
  do i=1,3
    read(uwout,*)(vect(j,i),j=1,3),(revect(j,i),j=1,3)
  enddo
  vect=vect/alat
  revect=revect*alat

! Read coordinates ------------------
  rewind(uwout)
  where=0
  ierr=0
  do while(where.eq.0.and.ierr.eq.0)
    read(uwout,'(a)',iostat=ierr)line
    where=index(line,'ions in cartesian')
  enddo
  if(ierr.ne.0)then
    write(6,'(a)')'  ERROR: No coordinates in the ouput file'
    stop
  endif
  allocate(spec(nat),coord(3,nat))
  do i=1,nat
    read(uwout,*)(coord(j,i),j=1,3)
  enddo

  allocate(sp1(nspec))
  rewind(uwout)
  where=0
  ierr=0
  j=0
  do while(where.eq.0.and.ierr.eq.0)
    read(uwout,'(a)',iostat=ierr)line
    where=index(line,'VRHFIN')
    if(where.ne.0)then
      where=0
      l=0
      do while(line(l:l).ne.':')
        l=l+1
        if(line(l:l).eq.'=')m=l 
      enddo
      j=j+1
      sp1(j)=line(m+1:l-1)
    endif
  enddo
  allocate(iopt(nspec))
  call read_out(uwout,'ions per type',cvalue)
  if(cvalue.eq.'')then
    write(6,'(a)')'  ERROR: No value of # of ion types in the output file'
    stop
  else
    read(cvalue,*)(iopt(i),i=1,nspec)
  endif
  j=0
  do i=1,nspec
    spec(j+1:j+iopt(i))=sp1(i)
    j=j+iopt(i)
  enddo
  deallocate(iopt)
  deallocate(sp1)

endif

! Read k points ---------------------
rewind(uwout)
where=0
ierr=0
do while(where.eq.0.and.ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'NKPTS')
enddo
if(ierr.ne.0)then
  write(6,'(a)')'  ERROR: No k points in the ouput file'
else
  cvalue=line(where+7:len(line))
  read(cvalue,*)nk
endif
allocate(kcart(4,nk))

rewind(uwout)
where=0
ierr=0
do while(where.eq.0.and.ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'Following cartesian')
enddo
if(ierr.ne.0)then
  write(6,'(a)')'  ERROR: No cart. coord. k points in the ouput file'
  stop
endif
read(uwout,*)aa
twk=0
do i=1,nk
  read(uwout,*)(kcart(j,i),j=1,4)
  twk=twk+kcart(4,i)
enddo
kcart(4,:)=kcart(4,:)/twk

! Read eigenvalues ------------------
rewind(uwout)
where=0
ierr=0
nline=1
do while(where.eq.0.and.ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'E-fermi')
  nline=nline+1
enddo
if(ierr.ne.0)then
  write(6,'(a)')'  ERROR: The self-consistent calc. has not finished'
  stop
endif
nek=kss

allocate(enk(nek,nk,nspin))

where=0
ierr=0
l=1
m=0
ispin=1
do while(ierr.eq.0)
  read(uwout,'(a)',iostat=ierr)line
  where=index(line,'k-point')
  if(where.ne.0)then
    read(uwout,'(a)')dmy
    do i=1,nek
      read(uwout,*)m,enk(i,l,ispin)
    enddo
    where=0
    if(l.eq.nk)then
      l=0
      if(nspin.eq.2)ispin=2
    endif
    l=l+1
  endif
enddo
if(m.eq.0)then
  write(6,'(a)')'  ERROR: No eigenvalues printed in the output file'
  stop
endif

! Spin degeneracy -/-----------------
sdeg=1.d0
if(nspin.eq.1)sdeg=2.d0

end subroutine vasp_data

!--------------------------------------------------------------------

subroutine checks()

! Check: cell volume ---------------
!if(rgeom)then
!  allocate(auxv(3))
!  call vect_prod(vect(:,1),vect(:,2),auxv)
!  print*,alat**3.d0*dot_product(auxv,vect(:,3))
!  if(rtestdat)print*,ucvol
!  deallocate(auxv)
!endif

! Check: cart and frac k-points ----
!if(rgeom)then
!  allocate(kfrac(4,nk))
!  do i=1,nk
!    call cart_to_frac(vect,kcart(1:3,i),kfrac(1:3,i))
!    write(6,'(3f12.8)')kfrac(1:3,i)
!    kcart(1:3,i)=0.d0
!    call frac_to_cart(revect,kfrac(1:3,i),kcart(1:3,i))
!    write(6,'(3f12.8)')kcart(1:3,i)
!  enddo
!  deallocate(kfrac)
!endif

! Check: sum of kweights -----------
!twk=0.d0
!do i=1,nk
!  twk=twk+kcart(4,i)
!enddo
!print*,'w',twk

! Check: number of e- --------------
!nel=0.d0
!do i=1,nk
!  do j=1,nek
!    do ispin=1,nspin
!      nel=nel+kcart(4,i)*fermif(enk(j,i,ispin)-ef,temp)
!    enddo
!  enddo
!enddo
!nel=nel*sdeg
!print*,nel
!if(rtestdat)print*,eln

! Check: number of eigenvalues -----
!print*,kss,nek

! Check: total energy (DC) ---------
!te=0.d0
!do i=1,nk
!  do j=1,nek
!    te=te+kcart(4,i)*enk(j,i,1)*fermif(enk(j,i,1)-ef,temp)
!  enddo
!enddo
!print*,te

end subroutine checks

!--------------------------------------------------------------------

end module ioabc_m

!--------------------------------------------------------------------
