!%********************************************************************************************************************!!
!%      Author: John Kattirtzi, Date: 31/05/09
!%                                                      Module: basisset
!%
!%      This module sets up the basis set before the trajectory is run. 
!%      It will read in the parameters and initialise the basis set. 
!%      The basis set will be in a derived type.
!%      The operators for the basis set will be overloaded here. 
!%      This module will specifically consist of
!%      functions to overload parameters                    done
!%      subroutines for:
!%      reading in the parameters                            done
!%      initialising z                                       done
!%      initialising s                                       done
!%      initialising d, including:                           done
!%      1. calc zp                                          done-check
!%      2. calc c initial                                   done 
!%      3. calc overlap initial                             done                           
!%      4. go from c to d initial                           done
!%      need to overload operators here                     done
!%      functions that will go from z to p and q            done
!%      subroutine for norm                                 done                         
!%        consider: 
!%      subroutine to initialise everything
!%
!%      Note: Above was written a while back..check it is still true
!%********************************************************************************************************************!!

MODULE BSETMK
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

!$********************************************************************************************************************!!
!$                                              Subroutine READBSPAR
!$      This subroutine will read the parameters for the basis set
!$      It will read:
!$      The number of dimensions/degrees of freedom (ndim)
!$      The number of basis functions (nbf)
!$      The value of hbar            - maybe it shouldn't read this and set it to 1.0d0
!$      Note- these parameters are not in arrays. The next subroutine deals with that
!$********************************************************************************************************************!!
  SUBROUTINE READBSPAR(LINE1)
    IMPLICIT NONE
    Character(LEN=100)::LINE1
    integer::ierr=0
    OPEN(UNIT=27, FILE='input.dat',STATUS='OLD', iostat=ierr)
    If (ierr .ne. 0) then
       Print *, 'error in opening input.dat file'
    End If
    read(27,*,iostat=ierr)LINE1
    do while (ierr==0)
       if(LINE1== "ndim") then
          backspace(27)
          read(27,*,iostat=ierr)LINE1,ndim ! could write a subroutine to put this in the type
          if(ierr.ne.0) then
             Print *,  "Error reading ndim"
          endif
       else if (LINE1=="hbar") then
          backspace(27)
          read(27,*,iostat=ierr)LINE1,hbar
          if(ierr.ne.0) then
             Print *,  "Error reading hbar"
          endif
       else if (LINE1=="nbf") then
          backspace(27)
          read(27,*,iostat=ierr)LINE1,nbf
          if(ierr.ne.0) then
             Print *,  "Error reading nbf"
          endif
       else if (LINE1=="sym") then
         backspace(27)
         read(27,*,iostat=ierr)LINE1,symop
         if(ierr.ne.0) then
           Print *, "Error symmetry option"
         endif
       end if
       read(27,*,iostat=ierr) LINE1
    end do
    close(27)
  END SUBROUTINE READBSPAR

!$********************************************************************************************************************!!
!$                                              Subroutine ReadArr
!$      This subroutine will read in the parameters that are in arrays of ndim,
!$      i.e:
!$      m-  mass
!$      gam-gamma
!$      w-  omega
!$********************************************************************************************************************!!

  SUBROUTINE READARR(LINE2)
    IMPLICIT NONE
    Character(LEN=100)::LINE2
    integer::ierr
    integer::i,k
    integer::j
    OPEN(27, FILE='input.dat', STATUS='OLD', iostat=ierr)
    read(27,*,iostat=ierr)LINE2
    do while (ierr==0)
       if(LINE2=="gamma")then
          backspace(27)
          do k=1,ndim
             read(27,*,iostat=ierr)LINE2,j,gam(k)
             if(.not. j==k) then
                print *,'error in i and j in gamma'
                Print *, 'check you have the correct number of gammas and their labels'
             end if
          end do
          if (ierr.ne.0) then
             Print *, "Error read gamma"
          end if
       else if (LINE2=="wfreq")then
          backspace(27)
          do k=1,ndim
             read(27,*,iostat=ierr)LINE2,j,w(k)
             if(.not. j==k) then
                print *,'error in i and j in wfreq'
                Print *, 'check you have the correct number of ws and their labels'
             end if
          end do
          if (ierr.ne.0) then
             Print *, "Error read wfreq"
          end if
       else if (LINE2=="wfreq")then
          backspace(27)
          do k=1,ndim
             read(27,*,iostat=ierr)LINE2,j,w(k)
             if(.not. j==k) then
                print *,'error in i and j in wfreq'
                Print *, 'check you have the correct number of ws and their labels'
             end if
          end do
          if (ierr.ne.0) then
             Print *, "Error read wfreq"
          end if
       else if (LINE2=="mass")then
          backspace(27)
          do k=1,ndim
             read(27,*,iostat=ierr)LINE2,j,m(k)
             if(.not. j==k) then
                print *,'error in i and j in mass'
                Print *, 'check you have the correct number of ms and their labels'
             end if
          end do
          if (ierr.ne.0) then
             Print *, "Error read mass"
          end if
       end if
       
       !write a subroutine to put this in the type
       read(27,*,iostat=ierr) LINE2
    end do
    close(27)
    
  END SUBROUTINE READARR
!$********************************************************************************************************************!!
!$                                              Subroutines:INBASISFN,IBASISSET,ALCBASIS
!$      The aim for these subroutines is to allocate the entire basis set
!$      The z array in the basis function is allocated in INBASISFN
!$      The IBASISSET routine calls the INBASISFN for each bf in the basis set
!$      The ALCBASIS routine allocates the bs as an array of n basis functions
!$      and calls IBASISSET to allocate each basis function within it. 
!$
!$      Note- there is a way to do it so that this is done in 1 routine but
!$      having them split up helps for afterwards 
!$********************************************************************************************************************!!
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
!$                                              Subroutine Readz
!$      This subroutine will read in the parameters that are needed to
!$      initialise z. These include:
!$      myseed- seed needed in random number generator- a value of 0 automatically generates one
!$      outp- gives a unit file number where the seed will be written to
!$      mup - gives the mean for the random gaussian distribution of p (momentum)
!$      muq - gives the mean for the random gaussian distribution of q (positions)
!$      sigp- gives the standard deviation for the random gaussian distb of p
!$      sigq- gives the standard deviation for the random gaussian distb of q
!$      LINE1 reads in the name of parameter
!$********************************************************************************************************************!!
  Subroutine Readz(myseed,outp,mup,muq,sigp,sigq, LINE1)
    Implicit none
    integer(kind=8),intent(out)::myseed
    integer, intent(inout)::outp   ! sort this out when writing output section
    real(kind=8),intent(inout),dimension(:),allocatable::mup
    real(kind=8),intent(inout),dimension(:),allocatable::muq
    real(kind=8),intent(inout)::sigp
    real(kind=8),intent(inout)::sigq
    character(LEN=100)::LINE1
    integer::ierr=0
    integer::i,j
    IF (.not. allocated(mup)) then
       Print *, 'WARNING not allocated mup in readz'
    END IF
    IF (.not. allocated(muq)) then
       Print *, 'WARNING not allocated muq in readz'
    END IF

    OPEN(27, FILE='input.dat', STATUS='OLD',iostat=ierr)

    If (ierr .ne. 0) then
       Print *, 'error in opening input.dat file'
    End If
    
    do while (ierr==0)
       if( LINE1=="outp") then
          backspace(27)
          read(27,*,iostat=ierr)LINE1, outp
          if (ierr.ne.0) then
             Print *, 'ierr error in outp', ierr
          end if
       else if (LINE1=="meanp") then
          backspace(27)
          do i=1,ndim
             read(27,*,iostat=ierr)LINE1,j,mup(i)
             if (j.ne.i) then
                Print *, 'error in mup i and j'
             end if
             if(ierr.ne.0) then
                Print *,  "Error(ierr) reading meanp"
             endif
          end do
       else if (LINE1=="meanq") then
          backspace(27)
          do i=1,ndim
             read(27,*,iostat=ierr)LINE1,j,muq(i)
             if (j.ne.i) then
                Print *, 'error in muq i and j'
             end if
             if(ierr.ne.0)then
                Print *, "Error(ierr) reading muq"
             end if
          end do
       else if (LINE1=="sigp")then
          backspace(27)
          read(27,*,iostat=ierr)LINE1,sigp
          if(ierr.ne.0) then
             Print *,  "Error reading sigp"
          endif
       else if (LINE1=="sigq")then
          backspace(27)
          read(27,*,iostat=ierr)LINE1,sigq
          if(ierr.ne.0) then
             Print *,  "Error (ierr) reading sigq", ierr
          endif
       else if (LINE1=="SEED")then
          backspace(27)
          read(27,*,iostat=ierr)LINE1,myseed
          if(ierr.ne.0) then
             Print *,  "Error (ierr) reading seed"
          endif
       end if
       read(27,*,iostat=ierr)LINE1! exit loop so ierr should not equal to 0 here
    end do
    close(27)
  END Subroutine Readz

!$********************************************************************************************************************!!
!$                                              Subroutine IntZ
!$      This subroutine will initialise the z in the basis function
!$      ZBQLINI is called first to SEED the random number routine
!$      ZBQLNOR is then called to generate a random number in a gaussian distb
!$      ZBQLNOR needs to call the mean and standard dev for the gaussian
!$      outp gives the unit file name that the SEED will be printed to
!$      z=SQRT(gam/2)*q+i/hbar*SQRT(1/2GAM)p
!$********************************************************************************************************************!!

  SUBROUTINE IntZ(myseed,outp,muq,sigq,mup,sigp,bs,ZBQLNOR)
    IMPLICIT NONE
    type(csbasisfn),dimension(:),intent(inout)::bs
    real(kind=8),external::ZBQLNOR
    type(csbasisfn)::bf
    type(csbasisfn)::bfconjz
    integer(kind=8),intent(in)::myseed
    integer, intent(in)::outp   ! sort this out when writing output section
    real(kind=8),intent(in),dimension(:),allocatable::mup
    real(kind=8),intent(in),dimension(:),allocatable::muq
    real(kind=8),intent(in)::sigp
    real(kind=8),intent(in)::sigq
    integer::i,j, nbfhalf 
    real(kind=8)::RLnbf, Modnbf
    RLnbf=REAL(nbf)
    IF (.not. allocated(mup)) then
       Print *, 'WARNING not allocated mup in intz'
    END IF
    IF (.not. allocated(muq)) then
       Print *, 'WARNING not allocated muq in intz'
    END IF
    CALL INBASISFN(bfconjz)
    CALL INBASISFN(bf)
    CALL ZBQLINI(myseed,outp)! this will automatically write the SEED to an outp
    nbfhalf=0.5d0*nbf
    IF ((symop == "NO") .OR. (nbf==1)) then
       IF (symop =="YES") then
          Print *, 'This isnt a warning, it is more of a reminder'
          Print *, 'You have asked for a symmetrical one basis fnc basis set'
          Print *, 'Either use more basis function or put sym NO'
       END IF   
       do i=1,nbf
          bf=bs(i)
          do j=1,ndim
             bf%z(j)=cmplx(SQRT(gam(j)/2.0d0)*(ZBQLNOR(muq(j),sigq))&               
                  ,(1.0d0/hbar*SQRT(1.0d0/(2.0d0*gam(j)))*ZBQLNOR(mup(j),sigp)))
             !   Print *, 'bf%z(j)', bf%z(j), 'bf', i, 'dim', j
          end do
          bs(i)=bf   
       end do
    ELSE IF (symop =="YES") then
       Print *, 'Symmetry in the basis set has been specified'
       Modnbf=mod(RLnbf,2.0d0)
        IF (Modnbf > 0.0d0) then
        Print *, 'WARNING How can I ensure a symmetrical odd numbered basis set?'
        STOP
       END IF      
       do i=1,nbfhalf
          bf=bs(i)
          do j=1,ndim
             bf%z(j)=cmplx(SQRT(gam(j)/2.0d0)*(ZBQLNOR(muq(j),sigq))&               
                  ,(1.0d0/hbar*SQRT(1.0d0/(2.0d0*gam(j)))*ZBQLNOR(mup(j),sigp)))
          end do
          bs(i)=bf
          bfconjz%s=bf%s
          bfconjz%d=bf%d
          bfconjz%z(1:ndim)=dconjg(bf%z(1:ndim))
          bs(i+nbfhalf)=bfconjz
       end do
    ELSE 
       Print *, 'Havent understood symop'
       Print *, 'Should be YES or NO in input.dat'
       Print *, 'Symop=', Symop
    END IF
    !do i=1,nbf
    !Print *,'bs(i)%z', bs(i)%z
    !end do
    close(27)
    deallocate (bf%z)
    deallocate (bfconjz%z)
  End Subroutine IntZ
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
!$                                              Subroutine ZPS
!$      This function will give the initial z (z0) that gets projected onto the
!$      basis set. 
!$      It consists of the mean of the gaussian distb for p and q
!$********************************************************************************************************************!!
 
  FUNCTION ZPS(muq,mup)
    ! A function that will give the exponent z in the zp
    IMPLICIT NONE
    real(kind=8),dimension(:),intent(in)::muq
    real(kind=8),dimension(:),intent(in)::mup
    complex(kind=8),dimension(ndim)::ZPS
    integer::i
    if ((size(muq) .ne. size(mup) ) .OR. (size(muq) ==0) .OR. (size(mup) == 0))  then
        Print *, 'error in size means ZPS function'
        Print *, 'size muq', size(muq), 'size mup', size(mup)
    end if
    ZPS(1:ndim)=cmplx(SQRT(gam(1:ndim)/2.0d0)*muq(1:ndim),1.0d0/hbar*SQRT(1.0d0/2.0d0*gam(1:ndim))* mup(1:ndim)) 
  End FUNCTION ZPS

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
!$                                                      Subroutine ReadChk
!$      ReadChk is a subroutine which reads through the input file and check to
!$      see if the number of lines that have been read are what is expected.
!$      There should be ndim values for:
!$      mass
!$      gamma
!$      angular frequency
!$      mean p
!$      mean q
!$      And the following is also expected to be read:
!$      ndim, hbar, nbf, SEED, outp, sigp, sigqa = 7
!$      This routine also checks that these terms are written correctly 
!$      and will tell you if it hasn't understood something. 
!$********************************************************************************************************************!!

  Subroutine ReadChk(Line1)
    ! subroutine that will read input again and check understood everything
    Implicit none
    character(LEN=100)::LINE1
    integer::ierr, er, nlines
    ierr=0
    OPEN(27, FILE='input.dat', STATUS='OLD',iostat=ierr)
    If (ierr .ne. 0) then
       Print *, 'error in opening input.dat file'
    End If
    read(27,*,iostat=ierr)LINE1
    er=0
    do while (ierr == 0)
       if (LINE1 =='ndim') then
          er=er+1    
       else if (LINE1=='hbar')then
          er=er+1
       else if (LINE1=='nbf') then
          er=er+1
       else if (LINE1=='sym') then
          er=er+1   
       else if (LINE1=='mass')then
          er=er+1
       else if (LINE1=='gamma')then
          er=er+1
       else if (LINE1=='wfreq')then
          er=er+1     
       else if (LINE1=='SEED') then
          er=er+1
       else if (LINE1=='outp') then
          er=er+1
       else if (LINE1=='meanp')then
          er=er+1
       else if (LINE1=='meanq') then
          er=er+1
       else if (LINE1=='sigp') then
          er=er+1
       else if (LINE1=='sigq') then
          er=er+1
       else  if ( .not. (LINE1(1:1) =="#")) then! Catches things not read
          Print *, "Option ",trim(LINE1) ," ignored"
          Print *, "I haven't understood this. Comment it out!"
       end if
    read(27,*,iostat=ierr) LINE1
    end do
    nlines=8+5*ndim
    if (er .ne. nlines) then
       Print *, 'I have read:', er,'lines, I expected', nlines,'from input.dat'
       Print *, 'ndim= ', ndim, 'nbf=', nbf, 'STOP!'
       STOP
    end if
    !   Print *, 'er=',er
  End Subroutine ReadChk

!$********************************************************************************************************************!!
!$                                              Subroutine IntD
!$      IntD initialises D in the basis set.
!$      It does this by first calculating c using the overlap of zp -propagated
!$      z with the z in the basis function. 
!$      After c is calculated the CTODin routine is called which gives the value
!$      of d. 
!$      It worked out easier to do calc c in this routine rather than using
!$      other routines to do this for me. 
!$********************************************************************************************************************!!


  Subroutine IntC(c, mup, muq, bs)
    ! subroutine to initialise d
    IMPLICIT NONE
    type(csbasisfn), dimension(:)::bs
    complex(kind=8), dimension(:), intent(out)::c
    complex(kind=8), dimension(:),allocatable::zp
    real(kind=8),dimension(:),allocatable,intent(in)::mup
    real(kind=8),dimension(:),allocatable,intent(in)::muq
    integer::i
    If (.not. allocated(mup)) then
    Print *, 'not allocated mup in IntD'
    End if
    IF (.not. allocated(muq)) then
    Print *, 'not allocated mup in IntD'
    End If

    allocate(zp(ndim))
    zp(:)=ZPS(muq,mup)
    do i=1,nbf
         c(i)=product(exp((dconjg(bs(i)%z(1:ndim))*zp(1:ndim))-(0.5d0*dconjg(bs(i)%z(1:ndim))*bs(i)%z(1:ndim))-&      ! Im not sure about the sum term here  
              &(0.5d0*dconjg(zp(1:ndim))*zp(1:ndim))))                                                                   ! important for ndim but not when eq 1
    end do
    deallocate(zp)
  END Subroutine IntC
     
  Subroutine ReadTime(Line1,t0,tmax,dt)
    real(kind=8)::t0
    real(kind=8)::tmax
    real(kind=8)::dt
    character(LEN=100)::Line1
    integer::ierr=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
    read(30,*,iostat=ierr)LINE1 ! unit 30 was chosen randomly
    If (ierr .ne. 0) then
       Print *, 'WARNING: Not reading Time in Int mod'
       Print *, 'ierr', ierr
    End If
    do while (ierr==0)
       if(Line1 == "t0") then
          backspace(30)
          read(30,*,iostat=ierr)Line1,t0
          If (ierr .ne. 0) then
             Print *, 'Error reading in t0 '
             Print *, 'ierr', ierr
          End If
       else if(Line1 == "tmax") then
          backspace(30)
          read(30,*,iostat=ierr)Line1,tmax
          If (ierr .ne. 0) then
             Print *, 'Error reading in tmax '
             Print *, 'ierr', ierr
          End If
       else if(Line1 == "dt") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,dt
          If (ierr .ne. 0) then
             Print *, 'Error reading in dt '
             Print *, 'ierr', ierr
          End If
       end if
       
       read(30,*, iostat=ierr)LINE1
    end do
    close(30)
  End Subroutine ReadTime
  


END MODULE BSETMK
