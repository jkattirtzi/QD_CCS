!%********************************************************************************************************************!!
!%  Author: John Kattirtzi, Date: 17th April                                                                          !!   
!%                                                         Module Hamiltonian Specific                                !!
!%                                                                                                                    !!
!%      The purpose of this module is to have a list of functions that define Hord and                                !!
!%      DhordDz1 for different systems. The systems that I am considering at the moment are:                          !! 
!%      System:                                                                                                       !!
!%      The free particle                                                        HordFree                             !!
!%      The Harmonic Oscillator                                                  HordHarm                             !!
!%      The Morse Oscillator                                                     HordMorse                            !!
!%      Hubbard Model- Needs to be considered and defined! - do the maths first                                       !!
!%                                                                                                                    !!
!%      The Hord function will provide the value for the matrix element Hord_ij                                       !!
!%      This module contains the core of the physics that is needed.                                                  !!
!%                                                                                                                    !!
!%      It would be a good idea to write a subroutine to read in the parameters needed                                !!
!%      in this module to make it more modular.                                                                       !!
!%                                                                                                                    !!
!%      Subroutines needed:                                                                                           !!
!%      Hord for each system                                                             done                         !!
!%      DhordDz1 for each system                                                         done                         !!
!%      Read in the parameters for the system                                            done                         !!
!%      Allocate memory for each parameter in an array                                   done                         !!
!%                                                                                                                    !!
!%      Links with other modules:                                                                                     !!
!%      The basisset module will be used in this module as ndim,hbar,gamma are required. To                           !!
!%      ensure the parameters that are read do not get mixed up between the two modules                               !!
!%      (e.g omega) read from a different input file (call it inham.dat).                                             !!
!%      The HamGen module will need to use this module to work too.                                                   !!   
!%                                                                                                                    !!   
!%                                                                                                                    !!
!%********************************************************************************************************************!!


MODULE HamiltonianSpec
use BSETMK, only:ndim, hbar, gam, m,w
use BSETMK, only:csbasisfn

!! need to use BSET as it has parameter ndim in it


! I will write this in terms of z but it could be written in terms of p and q with the function from basisset mod     !!
! I'm going to assume that the parameters required will be one per degree of freedom 
!$********************************************************************************************************************!!
!$                                                         Global Parameters
!$      It would be nice in terms of efficiency to have a small number of
!$      parameters that can be used globally in this module rather than having specific
!$      parameers for each system. I.e instead of having m_free, m_harm and m_morse
!$      (mass in each system) I will have one mass parameter called m and will use
!$      this in each system. Parameters explained:
!$      m: mass of the degree of freedom
!$      k: Force Constant (used in Harmonic Oscillator, could also use in Morse)
!$      D: Well Depth in Morse Oscillator
!$      a: exp parameter in Morse a=SQRT(k/2D)
!$      These parameters are arrays of (ndim=n dimensions/degrees of freedom)
!$      Each degree of freedom might have a different value
!$      Hubbard Model parameters (not sure what these mean exactly :---------------------- fill in when you do!
!$      eps: epsilon = -0.5
!$      nu:  nu=1
!$      g:   g=-1/10 or -3/10
!$********************************************************************************************************************!!
!real(kind=8),dimension(:),allocatable::mh
real(kind=8),dimension(:),allocatable::kh
real(kind=8),dimension(:),allocatable::Dh
real(kind=8),dimension(:),allocatable::ah
real(kind=8)::eps=0.0d0
real(kind=8)::nu=0.0d0
real(kind=8)::g=0.0d0
contains
!subroutine to read in parameters
!              free particle
!              harmonic oscillator
!              Morse Oscillator
!              Hubbard model

!$*******************************************************************************************************************!!
!$                                                      Subroutine ALCHamPar
!$      This subroutine will allocate the memory for the global parameters array
!$      defined above. The ndim is taken from the basisset module
!$*******************************************************************************************************************!!
  Subroutine ALCHamPar
    ! No arguments are required as gloabal parameters are allocated
    IMPLICIT NONE
 !  allocate(mh(ndim))! not sure if this needs allocating as mass already read in basisset
    allocate(kh(ndim))
    allocate(Dh(ndim))
    allocate(ah(ndim))
  End Subroutine ALCHamPar

!$*******************************************************************************************************************!!
!$                                                      Subroutine  ReadHamSpec
!$      This subroutine will read in the global parameters from the input file
!$      The parameters have to be grouped together in the input file
!$      Line1 is called as an argument like all the reading rubroutines
!$      This is so that i can write a seperate subroutine at the end and call
!$      the other readin subroutines with the same line1 argument.  
!$*******************************************************************************************************************!!
  Subroutine ReadHamSpec(LINE1, NameSys)
    Implicit None
    Character(LEN=100)::LINE1
    Character(LEN=20)::NameSys
    integer::ierr, i,j
    ierr=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
    read(30,*,iostat=ierr)LINE1 ! unit 30 was chosen randomly
    do while (ierr==0)

       if (NameSys=="Notread") then
          Print *, 'WARNING you are reading parameters before system name!!'
          Print *, 'exit the loop now'
          ierr=10
       end if
 !      if(LINE1== "mass") then
 !         backspace(30)
 !         do i=1,ndim
 !            read(30,*,iostat=ierr)LINE1,j,mh(i)
 !            if (.not. j==i) then
 !              ! statement below checks mass 2 wont go in mass(2)     
 !               print *, 'error in i and j mass reading'
 !               Print *, 'check you have the correct number of masses'
 !            end if
 !         end do
 !         if(ierr.ne.0) then
 !            Print *, "Error reading m"
 !         end if
       if (NameSys=="Morse") then

          if(LINE1=="alpha") then ! morse oscillator parameter
             backspace(30)
             do i=1,ndim
                read(30,*,iostat=ierr)LINE1,j,ah(i)
                if (.not. j==i) then
                   Print *, 'error in i and j in alpha'
                end if
             end do
             if (ierr.ne.0) then
                Print *, 'error in alpha ierr'
             end if
          else if(LINE1=="MorseD") then ! Morse well
             backspace(30)
             do i=1,ndim
                read(30,*,iostat=ierr)LINE1,j,Dh(i)
                if (.not. j==i) then
                   Print *, 'error in i and j in DH'
                   Print *, 'i=', i, 'j=', j, 'DH', Dh
                end if
             end do
             if (ierr.ne.0)then
                Print *, 'error in ierr in MorseD'
             end if
          end if
          
       else if(NameSys=="Harmonic") then
          if(LINE1=="ForceK")then! Harmonic parameter
             backspace(30)
             do i=1,ndim
                read(30,*,iostat=ierr)LINE1,j,kh(i)
          Print *, 'reading harmonic parameters', kh(i)
                if (.not. j==i) then
                   Print *, 'error in i and j forcek'
                end if
             if (ierr.ne.0) then
                Print *, 'error in ierr in forcek'
             end if
             end do
          end if
       
       else if (NameSys=='Hubbard') then
          if (ndim .ne. 2) then
             Print *, 'WARNING ndim ne 2'
             Print *, 'For Hubbard ndim=2!!!, STOPPING NOW'
             STOP
          End if
          if(LINE1=="HubEps") then
             backspace(30)
             read(30,*,iostat=ierr)LINE1, eps
             if (ierr .ne. 0) then
              Print *, 'error in reading HubEps', ierr
             end if
          else if (LINE1=="HubG") then
            backspace(30)
            read(30,*,iostat=ierr)LINE1, g    
            if (ierr .ne. 0) then
              Print *, 'error in reading HubG', ierr  
            end if
          else if (LINE1=='HubNu') then
            backspace(30)
            read(30,*,iostat=ierr)LINE1, nu
            if (ierr .ne. 0) then
               Print *,'error in reading HubNu',ierr
            end if
          end if
       end if
       read(30,*,iostat=ierr) LINE1
    end do
    close(30)
  END Subroutine ReadHamSpec

!$*******************************************************************************************************************!!
!$                                                      Subroutine HarmChk
!$      This subroutine will take in the arguments for the parameters of the Harmonic Oscillator and check
!$      that the force constant=w^2*m
!$      w angular frequency read in the basis set module
!$      m mass read in the basis set module
!$      kh - force constant read in the subroutine above (ReadHamSpec)
!$*******************************************************************************************************************!!

  Subroutine HarmChk
    Implicit none
    complex(kind=8),dimension(size(w))::wout
    integer::i
    wout(1:ndim)=SQRT(kh(1:ndim)/m(1:ndim))
    do i =1,ndim
       if (wout(i) .ne. w(i)) then
          Print *, 'WARNING, HarmChk in HamSpec, w neq SQRT(k/m)'
          Print *, 'Dimension', i, 'w',w(i), 'wout', w(i), 'k', kh(i),'mass', m(i)
       end if
    end do
  End Subroutine HarmChk

!$*******************************************************************************************************************!!
!$                                                      Function HordFree
!$      This function describes the free particle
!$      The arguments for this value are bf1cz and bf2
!$      bf1cz is the complex conjugate of the first bf(i.e z_i*)
!$      Need to calculate the complex conjg of first basis function before using this function
!$      bf2 is the second z value                     (i.e z_j)    
!$      Note:
!$      I'm not exactly sure about the multidimensionality- check this
!$ 
!$*******************************************************************************************************************!!

  Function HordFree(bf1cz,bf2)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    !complex(kind=8),intent(in)::z1c,z2
    complex(kind=8)::HordFree
    complex(kind=8)::HordFree2
    complex(kind=8)::bfbit
    !HordFree=sum(-1.0d0/(4.0d0*m(1:ndim)*hbar*gam(1:ndim))*((bf1cz%z(1:ndim)-bf2%z(1:ndim))**2))
  !  HordFree=sum(-1.0d0/(4.d0*(m(1:ndim)**2.0d0)*w(1:ndim))*(bf2%z(1:ndim)-bf1cz%z(1:ndim))**2.0d0)
    !HordFree=sum(-0.25d0*(bf2%z(1:ndim)-bf1cz%z(1:ndim))**2.0d0 + 0.25)
    HordFree=sum(-0.25d0*(bf2%z(1:ndim)**2.0d0+bf1cz%z(1:ndim)**2.0d0-2.0d0*bf1cz%z(1:ndim)*bf2%z(1:ndim))+0.25)
 !  +1.0d0/(4.0d0*(m(1:ndim)**2.0d0)*w(1:ndim)))
    !IF (HordFree .ne. HordFree2) then 
    !        Print *, 'Error in Hord Free'
    !        Print *, 'Hord Free', HordFree, '2', HordFree2, 'w', w
    !End if
    ! need to make sure you have :
    ! m
    ! hbar
    ! gamma
  end Function HordFree


 ! Function HordFree(bf:q

!$*******************************************************************************************************************!!
!$                                                      Function HordHarm
!$      This function describes the harmonic oscillator
!$      The arguments for this value are bf1cz and bf2
!$      bf1cz is the complex conjugate of the first bf (i.e z_i*)
!$      bf2 is the second z value                      (i.e z_j)    
!$      Note:
!$      I'm not exactly sure about the multidimensionality- check this
!$ 
!$*******************************************************************************************************************!!



  Function HordHarm(bf1cz,bf2)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
!    complex(kind=8),intent(in)::z1c,z2
    complex(kind=8)::HordHarm
!    integer::i
!    HordHarm=cmplx(0.0d0,0.0d0)
!    do i=1,ndim
!        HordHarm= HordHarm + (hbar*w(i)*(bf2%z(i)*bf1cz%z(i)+0.5d0))
!    end do
    !Print *, '', hbar
    !Print *, 'zstar', bf1cz%z(1:ndim), 'z', bf2%z(1:ndim)
    HordHarm=sum(hbar*w(1:ndim)*(bf2%z(1:ndim)*bf1cz%z(1:ndim)+0.5d0))! classical system no uncertainty
    !HordHarm=sum(hbar**2*(kh(1:ndim)/gam(1:ndim))*(bf1cz%z(1:ndim)*bf2%z(1:ndim)+0.5d0))
    !need hbar, gamma, m
  End Function HordHarm





!$*******************************************************************************************************************!!
!$                                                      Function HordMorse
!$      This function describes the Morse Oscillator
!$      The arguments for this value are bf1cz and bf2
!$      bf1cz is the complex conjugate of the first bf (i.e z_i*)
!$      bf2 is the second z value                    (i.e z_j) 
!$      ah is alpha- the exponent in the Morse potential
!$      Dh is the potential well
!$      gam is gamma in the basis set
!$      These parameters are global so do not need to be declared
!$      Note:
!$      I'm not exactly sure about the multidimensionality- check this
!$ 
!$*******************************************************************************************************************!!

  Function HordMorse2(bf1cz,bf2)
    Implicit None
    type(csbasisfn),intent(in)::bf1cz,bf2
    complex(kind=8)::HordMorse2
    complex(kind=8)::Dbit
    Dbit=sum(Dh(1:ndim)*(exp(ah(1:ndim)**2.0d0)*exp(-SQRT(2.0d0)*ah(1:ndim)*(bf1cz%z(1:ndim)+bf2%z(1:ndim)))&
    -2.0d0*exp((ah(1:ndim)**2.0d0)*0.25)*exp(-(ah(1:ndim)/SQRT(2.0d0))*(bf1cz%z(1:ndim)+bf2%z(1:ndim)))+1.0d0))
    HordMorse2=HordFree(bf1cz,bf2)+Dbit
  End Function HordMorse2


!$*******************************************************************************************************************!!
!$                                                      Function HordHubb
!$      This function describes the Hubbard Model
!$      It is a two dimensional hamiltonian and (for now) only two dimensions
!$      will  be programmed. 
!$      bf1cz is the complex conjugate of the first bf
!$      bf2 is the second basis function
!$*******************************************************************************************************************!!

  Function HordHubb(bf1cz,bf2)
    Implicit None
    type(csbasisfn),intent(in)::bf1cz,bf2
    complex(kind=8)::HordHubb
    complex(kind=8)::epbit, nubit,gbit
    !IF ((eps .eq. 0.0d0) .OR. (nu .eq. 0.0d0) .OR. (g .eq. 0.0d0)) then
    !   Print *, 'warning possible error in Hubbard parameters'
    !   Print *, 'eps', eps, 'nu', nu, 'g', g
    !   !STOP
    !END IF
    epbit=eps*(bf1cz%z(1)*bf2%z(1)-bf1cz%z(2)*bf2%z(2))
    nubit=nu*(bf1cz%z(1)*bf2%z(2)+bf1cz%z(2)*bf2%z(1))
    gbit=g*(bf1cz%z(1)*bf2%z(1)+bf1cz%z(1)*bf1cz%z(1)*bf2%z(1)*bf2%z(1)&
         +bf1cz%z(2)*bf2%z(2)+bf1cz%z(2)*bf1cz%z(2)*bf2%z(2)*bf2%z(2))
    HordHubb=epbit+nubit+gbit
  End Function HordHubb


! write DHord functions
!$********************************************************************************************************************!!
!$                                                      DHordDz1c Functions
!$      These functions are written for the three systems with the same
!$      argument.
!$      I think that DHordDz1c(bf1cz,bf2)=DhordDz2(bf2,bf1cz) -check this but should be ok
!$      so I haven't programmed the right hand side.
!$
!$      Modification on 23/04/09 - I'm changing this function so that i is an
!$      integer argument and I am removing the sum. This makes things easier
!$      when calculating the delta2 Hord in the GenHam module       
!$********************************************************************************************************************!!

  Function DHordDz1cHubb(bf1cz,bf2,i)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    integer::i
    complex(kind=8)::DHordDz1cHubb
    IF (i .eq. 1) then
       DHordDz1cHubb=eps*bf2%z(1)+nu*bf2%z(2)+&
            g*(bf2%z(1)+2*bf1cz%z(1)*bf2%z(1)*bf2%z(1))
    Else IF (i .eq. 2) then
       DHordDz1cHubb=eps*(-1.0d0*bf2%z(2))+nu*bf2%z(1)+&
            g*(bf2%z(2)+2*bf1cz%z(2)*bf2%z(2)*bf2%z(2))
    Else IF (i .ne. 1 .AND. i .ne. 2) then
       Print *, 'Error in i, i=', i
       Print *, 'For Hubbar i=2 max'
       STOP
    End IF
  End Function DHordDz1cHubb



  Function DHordDz1cFree(bf1cz,bf2,i)
  ! DHordDz1c(bf1cz,bf2)=DhordDz2(bf2,bf1cz)- might need to check this but should be ok
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    complex(kind=8)::DHordDz1cFree
    integer,intent(in)::i
    integer::j
    do j =1,ndim
       IF (bf1cz%z(j) .ne. dconjg(bf2%z(j))) then
          Print *, ' bf1cz ne bf2c, j=',j
       END IF
    end do
    !DHordDz1cFree=1.0d0/(2.0d0*m(i)*gam(i)*hbar)*(bf2%z(i)-bf1cz%z(i))
    DHordDz1cFree=-0.50d0*(bf1cz%z(i)-bf2%z(i))
  End Function DHordDz1cFree


  Function DHordDz1cHarm(bf1cz,bf2,i)
    Implicit none
    type(csbasisfn),intent(in)::bf2
    type(csbasisfn)::bf1cz ! think I'll get an error if this is an input, for this system seems independent of this
    integer,intent(in)::i
    complex(kind=8)::DHordDZ1cHarm
    integer::j
    do j =1,ndim
         IF (bf1cz%z(j) .ne. dconjg(bf2%z(j))) then
            Print *, ' bf1cz ne bf2c, j=',j
         END IF
    end do
     DHordDz1cHarm=hbar*w(i)*bf2%z(i)
  End Function DHordDz1cHarm

  Function DHordDz1cMorse2(bf1cz,bf2,i)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    integer,intent(in)::i
    complex(kind=8)::DHordDz1cMorse2
    complex(kind=8)::Dbit
    Dbit=Dh(i)*(exp(ah(i)**2.0d0)*((-1.0d0)*SQRT(2.0d0)*ah(i))*exp(-SQRT(2.0d0)*ah(i)*(bf1cz%z(i)+bf2%z(i)))&
    -2.0d0*exp((ah(i)**2.0d0)*0.25)*(-ah(i)/SQRT(2.0d0))*exp(-(ah(i)/SQRT(2.0d0))*(bf1cz%z(i)+bf2%z(i))))
    DHordDz1cMorse2=DHordDz1cFree(bf1cz,bf2,i) + DBit
  End Function DHordDz1cMorse2

!$********************************************************************************************************************!!
!$                                                      DHordDz1 Functions
!$      These functions do the same as the DHordDz1c functions but for the first argument.
!$      These are required for the delta2Hord function
!$********************************************************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! I don't need the functions below at the moment but its worth having for later
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Function DHordDz1Free(bf1cz,bf2,i)
    ! DHordDz1c(bf1cz,bf2)=DhordDz2(bf2,bf1cz)- might need to check this but should be ok
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    complex(kind=8)::DHordDz1Free
    integer,intent(in)::i
    integer::j
    do j =1,ndim
       IF (bf1cz%z(j) .ne. dconjg(bf2%z(j))) then
          Print *, ' bf1cz ne bf2c, j=',j
       END IF
    end do
    DHordDz1Free=-0.5d0*(bf2%z(i)-bf1cz%z(i))
  End Function DHordDz1Free
  
  Function DHordDz1Harm(bf1cz,bf2,i)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    complex(kind=8)::DHordDz1Harm
    integer,intent(in)::i
    integer::j
    do j=1,ndim
       IF (bf1cz%z(j) .ne. dconjg(bf2%z(j))) then
          Print *, 'bf1cz ne bf2 in DhordDz1,j=', j
       END IF
    end do
    DHordDz1Harm= hbar*w(i)*bf1cz%z(i)
  END FUNCTION
  
END MODULE HamiltonianSpec
