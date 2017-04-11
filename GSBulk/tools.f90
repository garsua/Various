!--------------------------------------------------------------------

subroutine get_unit_number(lun)

implicit none
integer,intent(out) :: lun
integer             :: iostat
logical             :: used

do lun=10,99
  inquire(unit=lun,opened=used,iostat=iostat)
  if (iostat.ne.0)used=.true.
  if (.not.used)return
enddo
stop '  ERROR: Cannot get unit number'

end subroutine get_unit_number

!--------------------------------------------------------------------

subroutine read_in(uwin,descript,cvalue)

implicit none
integer,intent(in)            :: uwin
integer                       :: ierr,where,ll
character(len=1)              :: aa
character(len=256)            :: line
character*(*),intent(in)      :: descript
character(len=85),intent(out) :: cvalue

rewind(uwin)
where=0
ierr=0
do while(where.eq.0.and.ierr.eq.0)
  read(uwin,'(a)',iostat=ierr)line
  where=index(line,descript)
  if(where.ne.0)then
    ll=len(descript)
    aa=line(where+ll:where+ll)
    if(aa.ne.''.and.aa.ne.'='.and.aa.ne.':')where=0
  endif
enddo
ll=0
do while(aa.eq.'')
  ll=ll+1
  aa=line(ll:ll)
enddo
if(ierr.ne.0)then
  cvalue=''
else if(line(ll:ll).ne.descript(1:1))then
  cvalue=''
else
  where=index(line,'=')+1
  if(where.eq.1)where=index(line,':')+1
  if(where.eq.1)where=len(descript)+1
  read(line(where:len(line)),'(a40)')cvalue
endif

end subroutine read_in

!--------------------------------------------------------------------

subroutine read_out(uout,descript,cvalue)

implicit none
integer,intent(in)            :: uout 
integer                       :: ierr,where,ll
character(len=1)              :: aa
character(len=256)            :: line
character*(*),intent(in)      :: descript
character(len=85),intent(out) :: cvalue

rewind(uout)
where=0
ierr=0
do while(where.eq.0.and.ierr.eq.0)
  read(uout,'(a)',iostat=ierr)line
  where=index(line,descript)
  if(where.ne.0)then
    ll=len(descript)
    aa=line(where+ll:where+ll)
    if(aa.ne.''.and.aa.ne.'='.and.aa.ne.':')where=0
  endif
enddo
if(ierr.ne.0)then
  cvalue=''
else
  ll=index(line,descript)+len(descript)+2
  do while(aa.eq.'')
    aa=line(ll:ll)
    if(aa.eq.'='.or.aa.eq.':')aa=''
    ll=ll+1
  enddo
  read(line(ll-1:len(line)),'(a40)')cvalue
endif

end subroutine read_out

!--------------------------------------------------------------------

pure function intpar(line) result (cvalue)

implicit none
integer                       :: ni,nf
character(len=1)              :: aa
character(len=256),intent(in) :: line
character(len=85)             :: cvalue

ni=1
nf=0
aa=line(1:1)
do while(aa.ne.')')
  nf=nf+1
  aa=line(nf:nf)
  if(aa.eq.'(')ni=nf
enddo
cvalue=line(ni+1:nf-1)

end function intpar

!--------------------------------------------------------------------

subroutine vect_prod(vect1,vect2,vectpr)

implicit none
integer,parameter    :: dp=selected_real_kind(14,100)
integer              :: i,i1,i2
real(dp),intent(in)  :: vect1(3),vect2(3)
real(dp),intent(out) :: vectpr(3)

do i=1,3
  i1=mod(i,3)+1
  i2=mod(i+1,3)+1
  vectpr(i)=vect1(i1)*vect2(i2)-vect1(i2)*vect2(i1)
enddo

end subroutine vect_prod

!--------------------------------------------------------------------

subroutine cart_to_frac(revect,cart,frac)

implicit none
integer,parameter    :: dp=selected_real_kind(14,100)
integer              :: i
real(dp),intent(in)  :: revect(3,3),cart(3)
real(dp),intent(out) :: frac(3)

do i=1,3
  frac(i)=revect(1,i)*cart(1)+revect(2,i)*cart(2)+revect(3,i)*cart(3)
enddo

end subroutine cart_to_frac

!--------------------------------------------------------------------

subroutine frac_to_cart(vect,frac,cart)

implicit none
integer,parameter    :: dp=selected_real_kind(14,100)
integer              :: i
real(dp),intent(in)  :: vect(3,3),frac(3)
real(dp),intent(out) :: cart(3)

do i=1,3
  cart(i)=frac(1)*vect(i,1)+frac(2)*vect(i,2)+frac(3)*vect(i,3)
enddo

end subroutine frac_to_cart

!--------------------------------------------------------------------

integer function ninline(line)

implicit none
integer                :: i
character(len=1)       :: aa
character(len=256)     :: line
logical                :: kread

kread=.true.
ninline=0
do i=1,256
  aa=line(i:i)
  if(aa.ne.''.and.kread)then
    ninline=ninline+1
    kread=.false.
  else if(aa.eq.'')then
    kread=.true.
  endif
enddo

end function ninline

!--------------------------------------------------------------------

real(selected_real_kind(14,100)) function fermif(ee,tt)

implicit none
integer,parameter    :: dp=selected_real_kind(14,100)
real(dp),intent(in)  :: ee,tt
real(dp),parameter   :: evkb=11604.51928d0
real(dp)             :: fexp

fexp=ee*evkb/tt
fermif=1.d0/(1.d0+exp(fexp))

end function fermif

!--------------------------------------------------------------------
