!%********************************************************************************************************************!!
!%      Author: John Kattirtzi, Date: 20/08/09
!%                                                      Module: BSET2
!%      This module was originally intended to be part of the main
!%      program. It has now been modified so as to be part of a 
!%      smaller program that configures a basis set. 
!%      This module consist of subroutines for:
!%      allocate the memory of the basis set
!%      initialise s
!%      define basis function type- overload operators
!%      The overlap matrix
!%      Going from C-TO-D initially, but then going from D TO C
!%      Calculating the norm
!%      Note that the following comments were wrriten as part of 
!%      BSET module and may now be outdated
!%********************************************************************************************************************!!

MODULE BSET2
!$********************************************************************************************************************!!
!$                                              Derived Type: csbasisfn
!$
!$      This is a derived type which will define each individual basis function
!$      An array of these types (basis functions) will therefore give the basis
!$      set. 
!$      The type consists of a complex array of zs which will depend on n*
!$      dimensions (ndim). 
!$      It will also consist of an amplitude(complex) and an action variable
!$      (real).
!$********************************************************************************************************************!!

  type csbasisfn
     complex(kind=8),dimension(:),allocatable::z
     complex(kind=8)::d
     real   (kind=8)::s
  end type csbasisfn
!$********************************************************************************************************************!!
!$                      Basis Set Parameters:
!$      gam defines the gamma in the expression for z:
!$      z=SQRT(gam/2)*q+i/hbarSQRT(1/2gam)*p
!$      m defines the mass 
!$      w defines the angular frequency (omega)
!$      hbar is hbar and is initialised to 1.0
!$      nbf is the number of basis functions
!$      ndim is the number of dimensions
!$      The interface operators are needed to overload the operators
!$
!$      Note- think about initialising the others to 1 too and not reading them
!$
!$********************************************************************************************************************!!
  real(kind=8),dimension(:),allocatable::gam
  real(kind=8),dimension(:),allocatable::m
  real(kind=8)::hbar=1.0d0
  real(kind=8),dimension(:),allocatable::w ! omega- angular frequency
!  real(kind=8)::ma
  ! type parameters for memory can go here:
  integer::nbf   ! number of basis functions
  integer::ndim ! number of dimensions/degrees of freedom-might have this with nbf, skr suggests there
  character(LEN=5)::symop! determines whether basis set should be symmetrical
  interface operator (+)
     Module procedure add
  end interface

  interface operator (*)
     Module procedure mult
  end interface

contains

!$********************************************************************************************************************!!
!$                                              Functions: OVERLOAD
!$      The functions add and mult are required to overload the operators
!$      The function add, adds two basis functions togther
!$      It adds the corresponding components of the basis functions
!$      The function mult, multiplies a real number by a basis function
!$      These functions are required for the integration routine
!$      
!$********************************************************************************************************************!!
  Function add(bs1,bs2)
    Implicit none
    type(csbasisfn),dimension(:),intent(in)::bs1,bs2
    type(csbasisfn),dimension(size(bs1))::add
    integer::i
    call IBASISSET(add)! allocates memory for the z array
    If( (.not. size(bs1) ==size(bs2)) .OR. (size(bs1)==0 .OR. size(bs2) ==0)) then
            Print *, 'basis sets added arent equal sizes'
    End If
    do i=1,nbf
       add(i)%s=bs1(i)%s+bs2(i)%s
       add(i)%d=bs1(i)%d+bs2(i)%d
       !Print *, size(add%z), 'size'
       add(i)%z(1:ndim)=bs1(i)%z(1:ndim)+bs2(i)%z(1:ndim)
    end do
  END FUNCTION add

  Function mult(r,bs1)
    Implicit none
    type(csbasisfn),intent(in),dimension(:)::bs1
    real(kind=8),intent(in)::r               ! this will help check kind is always the same
    type(csbasisfn),dimension(size(bs1))::mult
    integer::i
    call IBASISSET(mult)! allocates memory for the z array
    do i=1,nbf
       mult(i)%s=bs1(i)%s*r
       mult(i)%d=bs1(i)%d*r
       mult(i)%z(1:ndim)=bs1(i)%z(1:ndim)*r
    end do
  End Function mult

  SUBROUTINE INBASISFN(bf)
    ! subroutine that will allocate each bf%z
    IMPLICIT NONE
    integer::n
    TYPE(csbasisfn), intent(inout)::bf
    ! n=ndim -is useful if you are having it in type
    !Print *, 'ndim in function', ndim
    allocate(bf%z(ndim))
  END SUBROUTINE INBASISFN

  SUBROUTINE IBASISSET(bs)
    ! subroutine that will do a loop to allocate the basis function array
    IMPLICIT NONE
    type(csbasisfn),dimension(:),intent(inout)::bs
    integer:: i
    do i=1,size(bs)
       call inbasisfn(bs(i))! this allocates each bf%z
    end do
  END SUBROUTINE IBASISSET
  ! need to allocate bs array before calling this routine

  SUBROUTINE ALCBASIS(bs)
    IMPLICIT NONE
    type(csbasisfn),dimension(:),allocatable::bs
    !integer::n ! if try to delcare this as in/inout/out get an error- not sure
    allocate(bs(nbf))
    CALL IBASISSET(bs)! calls another routine that allocates memory to each basis fn
    !n=bs(1)%ndim
  END SUBROUTINE ALCBASIS
!$********************************************************************************************************************!!
!$                                              Subroutine: ALCPAR
!$      This allocates the parameters that are ndim arrays- gamma, mass,angular frequency
!$********************************************************************************************************************!!

  Subroutine ALCPAR
    Implicit none
    allocate(gam(ndim))
    allocate(m(ndim))
    allocate(w(ndim))
    gam(1:ndim)=1.0d0    ! I don't have an input file for them now
     m(1:ndim)=1.0d0
     w(1:ndim)=1.0d0
  End Subroutine ALCPAR

!$********************************************************************************************************************!!
!$                                              Subroutine: BSGAMCHK
!$      This subroutine checks that the values of gamma, mass and omega inputed
!$      are sensible, i.e gamma inputed = m*w/hbar  
!$********************************************************************************************************************!!

  Subroutine BSGAMCHK
    Implicit none
    real(kind=8),dimension(:),allocatable::gamout
    integer::i
    allocate(gamout(size(gam)))
    do i=1,ndim
       gamout(i)=m(i)*w(i)/hbar
       IF (gamout(i) .ne. gam(i)) then
          Print *, 'WARNING gamma isnt mw/hbar in dimension', i
          Print *, 'gam', gam(i), 'gamout', gamout(i)
          Print *, 'm',m(i), 'wfrq', w(i), 'hbar', hbar
       end if
    end do
  End Subroutine BSGAMCHK
!$********************************************************************************************************************!!
!$                                              Subroutine IntS
!$      This subroutine will initialise the action in the basis set to 0.0d0
!$********************************************************************************************************************!!
  Subroutine IntS(bs)
    ! subroutine to initialise the action
    IMPLICIT NONE
    type(csbasisfn),dimension(:)::bs
    if (size(bs)==0) then
       Print *, 'error in size bs in IntS', size(bs)
    end if
    bs(1:nbf)%s=0.0d0
  End Subroutine IntS


!$********************************************************************************************************************!!
!$                                              Function OVRLELMT
!$      This functions gives the overlap element for two basis functions
!$      For two basis functions with multidimensions the results is:
!$      OVERLAPELMT=PRODUCT_dim( <bf1_d1 | bf2_d1> <bf1_d2|bf2_d2>
!       ...<bf1_dn|bf2_dn>) where d is the dimension index
!$********************************************************************************************************************!!

  FUNCTION OVRLPELMT(bf1,bf2)
    IMPLICIT NONE
    ! a function to calculate  the overlap between two matrix elements
    type(csbasisfn),intent(in)::bf1,bf2
    !complex(kind=8)::bf1,bf2
    complex(kind=8)::OVRLPELMT
    complex(kind=8)::expnt
    ! new variables for checking
    complex(kind=8)::z1,z2,z1c,z2c, z1cz2,z22,z11
    integer::i,j
    OVRLPELMT=cmplx(1.0d0,0.0d0)
    do i=1,ndim
        !do j=1,ndim
       z1=bf1%z(i)
       z2=bf2%z(i)
       z1c=dconjg(z1)
       z2c=dconjg(z2)
       z1cz2=z1c*z2
       z11=z1c*z1
       z22=z2*z2c
    !Print *, 'z2', z2,'z1cz2', z1cz2, 'z11', z11, 'z22', z22
       OVRLPELMT=OVRLPELMT*cdexp(z1cz2 - 0.5d0*z11- 0.5d0*z22)
       !end do
    end do
    If (size(bf1%z) .ne. size(bf2%z)) then
      Print *, 'error in sizes of basis fns in OVERLELMT'
    End if
  END FUNCTION OVRLPELMT

!$********************************************************************************************************************!!
!$                                              Function OVRLMAT
!$      This functions will give the overlap matrix using the OVRLELMT function
!$********************************************************************************************************************!!

  SUBROUTINE OVRLMAT(omega,bs)
 !   subroutine calculates the overlap matrix from the basis set array
    IMPLICIT NONE
    type(csbasisfn),dimension(:),intent(in)::bs
    complex(kind=8),dimension(:,:),intent(out)::omega
    integer ::i,j,n, nbf2
    nbf2=nbf**2 ! nbf is global variable
    If (size(omega) .ne. nbf2 ) then
      Print *, 'error in size of omega, stop!'
      Print *, 'size of omega', size(omega)
      STOP
    End if
    n=size(bs)
    if (n .ne. nbf) then 
       Print *, "error in OVRLMAT routine"
    end if
    do i=1,n-1
       omega(i,i)=cmplx(1.0d0,0.0d0)! ensures diagonal elements are 1           
       do j=i+1,n
          omega(i,j)=OVRLPELMT(bs(i),bs(j))! uses the overlap function written
          omega(j,i)=dconjg(omega(i,j))
       end do
    end do
    omega(n,n)=cmplx(1.0d0,0.0d0)! ensures last diagonal element is 1
    !Print *, 'omega', omega
  END SUBROUTINE OVRLMAT
!$********************************************************************************************************************!!
!$                                              Subroutine CTODin
!$      This subroutine takes in the input c and gives a value of d. 
!$      However, this can only be called before the trajectory is run.
!$      This is because the phase factor (e^(i/hbar)*si-sj) is ignored.
!$      This function can therefore be used to calculate initial c and initial
!$      norm but should not be used when the trajectory is being calculated.
!$      The equation that is solved is:
!$      c=OMEGA*d
!$      where OMEGA is the overlap matrix and c and d are arrays. 
!$      To solve this LAPACK is needed (ZGESV). 
!$********************************************************************************************************************!!
  SUBROUTINE CTODin(c,bs)
    ! subroutine to go from c to d initially - phase =1
    ! check for a phase factor if using not initially
    implicit none
    integer i,j,k, LDA,LDB,N, NRHS, INFO
    integer, DIMENSION(nbf)::IPIV
    complex(kind=8), DIMENSION(:,:),allocatable :: omega
    COMPLEX(kind=8), DIMENSION(:,:),allocatable::zlapin
    COMPLEX(kind=8), DIMENSION(:), INTENT(IN)::c
    complex(kind=8), DIMENSION(size(c))::cin
    complex(kind=8), dimension(size(c))::cout
    !type(csbasisfn), DIMENSION(:), allocatable::bs
    type(csbasisfn), dimension(:)::bs
    !     Declare character needed to read input
    !INTERFACE
    include "./LAPACK/ZGESV.f90"
    !END INTERFACE
    allocate(zlapin(nbf,nbf))
    allocate(omega(nbf,nbf))
    if ((size(c) == 0) .OR. (size(bs) ==0) .OR. (size(omega)==0)) then
       Print *, 'WARNING (in CTODin basisset mod )arrays assumed not allocated- will prob get segmentation fault!'
       Print *, 'size c', size(c), 'size bs', size(bs), 'omega', size(omega)
    end if
    call OVRLMAT(omega,bs)
   ! Print *, 'initial ovrlmat', omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !     routine to calculate d from c
    !     really need to put this in a seperate subroutine according to Dmitry
    LDA=nbf
    LDB=nbf
    NRHS=1 ! This was one before needs checking!!! could be 1 or ndim?
    N=nbf
    INFO=0
    IPIV=0
    zlapin=omega
    cin=c
    !Print *, 'zlapin',zlapin
    !Print *, 'cin', cin
    CALL ZGESV(N, NRHS, zlapin, LDA, IPIV, cin, LDB, INFO) ! be weary of compiler issues
    IF (INFO .ne. 0) then
       Print *, 'WARNING: INFO in CTODin', INFO
    END IF
    bs%d=cin
   ! cout=MATMUL(omega,bs%d)
   ! do i=1,nbf
    !  IF (cout(i) .ne. c(i)) then 
    !    Print *, 'error with LAPACK- cout neq cin in CTODin'
    !  END IF
   ! end do
    deallocate(zlapin)
    deallocate(omega)
    !Print *, 'dint from c', d
    !     Check that zlap*d =c
    !      c=MATMUL(zlap, d)
    !      Print *, 'c', c
    !      Print *, 'if c is the same as cin before then ds are sorted'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE CTODin

!$********************************************************************************************************************!!
!$                                              Subroutine DTOC
!$    This subroutine solves the simpler equation:
!$    c=Omega*D*e(i/hbar(si-sj)) where Omega and D are both known
!$    It takes phase into account and so can be used at any point
!$********************************************************************************************************************!!
  
  SUBROUTINE DTOC(c,bs)
    ! subroutine to go from d to c
    Implicit none
    complex(kind=8), dimension(:,:),allocatable::omega
    complex(kind=8), dimension(:), intent(out)::c 
    type(csbasisfn), dimension(:),intent(in)::bs
    integer::  j,k, l,m
    allocate(omega(nbf,nbf))
    if ((size(c) == 0) .OR. (size(bs) ==0)) then
       Print *, 'WARNING (in DTOC basisset mod )arrays assumed not allocated- will prob get segmentation fault!'
       Print *, 'size c', size(c), 'size bs', size(bs) 
    end if
    call OVRLMAT(omega,bs)! modified this after written prog.f90 so might need to change this
    !do l=1,nbf
    !    do m=1,nbf
    !           Print *, 'omega', omega(l,m)
    !    end do
    !end do
    c(1:nbf)=0.0d0
    do j=1,nbf
       do k=1,nbf
           c(j)=c(j)+(omega(j,k)*bs(k)%d*exp((cmplx(0.0d0,1.0d0)/hbar)*(bs(k)%s-bs(j)%s)))
        !  Print *, 'omegajk', omega(j,k)
       end do   
    !   Print *, 'c', c(j)
    end do
    deallocate(omega)
  End SUBROUTINE DTOC

!$********************************************************************************************************************!!
!$                                                      Subroutine Norml
!$      This suboroutine takes in a c value and a d value (from bs) and
!$      calculates the normal from it. 
!$      It can be used at any time but c must be calculated first. 
!$      norm =Sum_i (c*_i d_i_)
!$********************************************************************************************************************!!

  Subroutine Norml(c,norm,bs)
    ! subroutine to calc the norm
    Implicit None
    type(csbasisfn),dimension(:),intent(in)::bs
    complex(kind=8),dimension(:),intent(in)::c
    real(kind=8),intent(out)::norm
    if ((size(c) == 0) .OR. (size(bs) ==0)) then
       Print *, 'WARNING (in Norm1 basisset mod )arrays assumed not allocated- will prob get segmentation fault!'
       Print *, 'size c', size(c), 'size bs', size(bs) 
    end if
 
    norm=sum(dconjg(c(1:nbf))*bs(1:nbf)%d)
  !Print *,'in norml,c', c(1:nbf),'d', bs(1:nbf)%d
  End Subroutine Norml

!$********************************************************************************************************************!!
!$                                              Subroutine Norm2
!$      This subroutine calculates the norm without the need to calculate c
!$      prior to this. 
!$      It solves norm=sum d_i e(i/hbar (si-sj)) d*_j omega_ji
!$      This is equivalent to the expression for norm calculated with c e
!$********************************************************************************************************************!!

  Subroutine Norm2(bs,norm)
    Implicit none
    type(csbasisfn),dimension(:),intent(in)::bs
    real(kind=8),intent(out)::norm
    real(kind=8),dimension(nbf)::normarray
    complex(kind=8), dimension(:,:),allocatable::omega
    integer::h,i,j
    allocate(omega(nbf,nbf))
    call OVRLMAT(omega,bs)
    norm=0.0d0
    do i=1,nbf
       do j=1,nbf
          norm=norm+bs(i)%d*exp((cmplx(0.0d0,1.0d0)/hbar)*(bs(i)%s-bs(j)%s))*(dconjg(bs(j)%d))*omega(j,i)
       end do
    end do
    deallocate(omega)
  End Subroutine Norm2

!$********************************************************************************************************************!!
!$   
!$ A subroutine that will read in a basis set, z and c!
!$********************************************************************************************************************!!
  
  Subroutine ReadBSIn(bs,c)
    Implicit None
    type(csbasisfn),dimension(:), allocatable,intent(out)::bs
    complex(kind=8),dimension(:),allocatable::c
    real(kind=8):: cr, ci
    real(kind=8)::zr, zi
    integer::ierr, i, j, k
    character(LEN=10)::Line1
    ierr=0
    OPEN(72, FILE='basisset.dat', STATUS='OLD',iostat=ierr)
    IF (ierr .ne. 0) then
        Print *, 'error reading basisset.dat'
    END IF
    read(72,*,iostat=ierr)Line1, ndim
    If (ierr .ne. 0) then
       Print *, 'error reading ndim'
    End If
    Print *, Line1, ndim
    read(72,*, iostat=ierr), Line1, nbf
    If (ierr .ne. 0) then
       Print *, 'error reading ndim'
    End If
    Print *, Line1, nbf
    allocate(c(nbf))
    call ALCBASIS(bs)
    do i=1,nbf
       read(72,*,iostat=ierr), Line1, k
       IF ((ierr .ne. 0 ) .OR. (k .ne. i)) then
          Print *, 'error in reading basis function number'
          print *, 'i', i, 'k', k, 'iostat', ierr
          stop
       End IF
       read(72,*,iostat=ierr), Line1,cr, ci
       IF (.not. (Line1 =='C')) then
           Print *, 'error associated with reading c', Line1
           STOP
       END IF
       c(i)=cmplx(cr,ci)
       If (ierr .ne. 0) then
          Print *, 'error reading ndim'
       End If
       do j=1,ndim
          read(72,*,iostat=ierr), zr, zi
          bs(i)%z(j)=cmplx(zr,zi)
          If (ierr .ne. 0) then
             Print *, 'error reading ndim'
          End If
       end do
    end do
    Print *, 'calling CTODin'
    Print *, 'c', c
    call CTODin(c,bs)
    Print *, 'finished bs initialisation'
    close(72)
  End Subroutine ReadBSIn

END MODULE BSET2
