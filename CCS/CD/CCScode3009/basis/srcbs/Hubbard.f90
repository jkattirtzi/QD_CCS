Module Hubbard
use BSETMK, only:ndim,nbf, hbar, gam, m,w
use BSETMK, only:csbasisfn
use BSETMK, only:OVRLMAT

character(LEN=20)::popname='nopop.dat'
contains

!  Subroutine Nchange(bs,change,intdim)
!    Implicit none
!    type(csbasisfn),dimension(:),intent(in)::bs
!    complex(kind=8),intent(out)::change
!    complex(kind=8), dimension(:,:),allocatable::omega
!    integer::i,j
!    integer,intent(in)::intdim
!    allocate(omega(nbf,nbf))
!    call OVRLMAT(omega,bs)
!    do i=1,nbf
!       do j=1,nbf
!        change = change+bs(i)%d*exp((cmplx(0.0d0,1.0d0)/hbar)*(bs(i)%s-bs(j)%s))*(dconjg(bs(j)%d))&
!                  *omega(j,i)*dconjg(bs(j)%z(intdim))*bs(i)%z(intdim)
!       end do
!    end do
!    deallocate(omega)
! End Subroutine Nchange


  Subroutine ReadPopOp(LINE1)
    Implicit None
    Character(LEN=10)::LINE1
    integer::ierrM
    ierrM=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierrM)
    IF (ierrM .ne. 0) then
       Print *, 'error opening inham in ReadOutOp'
       STOP
    End IF
    read(30,*,iostat=ierrM)LINE1 ! unit 30 was chosen randomly
    do while (ierrM==0)
       if(LINE1 == "Pop") then
          backspace(30)
          read(30,*,iostat=ierrM)LINE1,popname
          !Print *, LINE1, popname
          IF (ierrM .ne. 0) then
             Print *, 'error reading Trajname in ReadOutOp'
             STOP
          END IF
       End if
       read(30,*,iostat=ierrM)LINE1
    end do
    IF (.not. (popname=="nopop.dat")) then
       call popouthead
    END IF
  End Subroutine ReadPopOp



  Subroutine popouthead
    Implicit none
    Open(UNIT=38, FILE=popname, STATUS='NEW')
    Write(38,*), '#********************Population Output File********************'
    write(38,*), '#**************************************************************'
    Write(38,*), '#time',' time value',' n1 pop', 'n2 pop', 'sum narray'
    close(38)
  End Subroutine Popouthead


  Subroutine popwrite(t,narray)
    Implicit none
    real(kind=8),intent(in)::t
    !real(kind=8),intent(in)::n1
    !real(kind=8),intent(in)::n2
    real(kind=8),dimension(:),intent(in)::narray
    Open(Unit=38, File=popname, STATUS='OLD', ACCESS='APPEND')
    write(38,*), t,narray, sum(narray)
    close(38)
  End Subroutine popwrite



  Subroutine Dim2Npop(t,bs)
    Implicit None
    type(csbasisfn),dimension(:),intent(in)::bs
    !integer::intdim1, intdim2
    integer,dimension(:),allocatable::intdimarray    
    !real(kind=8)::n1
    !real(kind=8)::n2
    real(kind=8),dimension(:),allocatable::narray    
    real(kind=8),intent(in)::t
    integer::intdim, i
    IF (.not. (popname=="nopop.dat")) then
       !intdim1=1
       !n1=NchangeFcn(bs,intdim1)
       !intdim2=2
       !n2=NchangeFcn(bs,intdim2)
       allocate(narray(ndim))
     !  allocate(intdimarray(ndim))
       do i=1,ndim
           intdim=i
           narray(i)=NchangeFcn(bs,intdim)
       End do   
    !   Print *, 'narrary2', narray(2)
    !   STOP 
       call popwrite(t,narray)
       deallocate(narray)
     !  deallocate(intdimarray)
    End IF
  End Subroutine Dim2Npop


  Function NchangeFcn(bs,intdim)
    Implicit none
    type(csbasisfn),dimension(:),intent(in)::bs
    complex(kind=8)::NchangeFcn
    complex(kind=8), dimension(:,:),allocatable::omega
    integer::i,j
    integer,intent(in)::intdim
    !Print *, 'bsiz2', bs(1)%z(2)
    !STOP
    NchangeFcn=0.0d0
    allocate(omega(nbf,nbf))
    call OVRLMAT(omega,bs)
    do i=1,nbf
       do j=1,nbf
        !Print *, 'i', i, 'j', j
        NchangeFcn = NchangeFcn+(bs(i)%d*exp((cmplx(0.0d0,1.0d0)/hbar)*(bs(i)%s-bs(j)%s))*(dconjg(bs(j)%d))&
                  *omega(j,i)*dconjg(bs(j)%z(intdim))*bs(i)%z(intdim))
       end do
    end do
    deallocate(omega)
 End Function NchangeFcn

 ! This function does that NchangeFcn but without the do loop

 Function NpopBf(bf,intdim)
   Implicit None
   type(csbasisfn),intent(in)::bf
   real(kind=8)::NpopBf
   integer,intent(in)::intdim
  ! Npopbf=bf%d*dconjg(bf%d)*exp((cmplx(0.0d0,1.0d0)/hbar)*(bf%s-bf%s))*bf%z(intdim)*bf%z(intdim)
   Npopbf=dconjg(bf%z(intdim))*bf%z(intdim)
  End Function
End Module Hubbard
