Module Hubbard
use BSET2, only:ndim,nbf, hbar, gam, m,w
use BSET2, only:csbasisfn
use BSET2, only:OVRLMAT

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


  Subroutine popwrite(t,narray, ncdarray)
    Implicit none
    real(kind=8),intent(in)::t
    real(kind=8),dimension(:),intent(in)::narray
    real(kind=8),dimension(:),intent(in)::ncdarray
    Open(Unit=38, File=popname, STATUS='OLD', ACCESS='APPEND')
    write(38,*), t,narray, sum(narray), ncdarray, sum(ncdarray)
    close(38)
  End Subroutine popwrite



  Subroutine Dim2Npop(t,bs,c)
    Implicit None
    type(csbasisfn),dimension(:),intent(in)::bs
    integer,dimension(:),allocatable::intdimarray    
    real(kind=8),dimension(:),allocatable::narray, ncdarray    
    real(kind=8),intent(in)::t
    integer::i
    complex(kind=8),dimension(:)::c
    IF (.not. (popname=="nopop.dat")) then
       allocate(narray(ndim))
       allocate(ncdarray(ndim))
       do i=1,ndim
           narray(i)=NchangeFcn(bs,i)
           ncdarray(i)=NpopCDin(bs,c,i)
       End do   
       call popwrite(t,narray,ncdarray)
       deallocate(narray)
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

  Function NpopCDin(bs,c,intdim)
    Implicit none
    type(csbasisfn),dimension(:),intent(in)::bs
    complex(kind=8),dimension(:)::c
    complex(kind=8)::NpopCDin
    integer::intdim, i
    complex(kind=8),dimension(:,:),allocatable::omega
    allocate(omega(nbf,nbf))
    call OVRLMAT(omega,bs)
    NpopCDin=cmplx(0.0d0,0.0d0)
    do i=1,nbf
        NpopCDin= NpopCDin +(dconjg(c(i))*bs(i)%d*exp(bs(i)%s)*bs(i)%z(intdim)*dconjg(bs(i)%z(intdim)))
    end do    
    deallocate(omega)
  END function

    
End Module Hubbard
