!%********************************************************************************************************************!!
!% Author: John Kattirtzi, Date: 21/04/09      
!%                                    Module General Hamiltonian
!%
!%      This module is written to read in the type of Hamiltonian and to
!%      assign the correct one to Hord (i.e HordFree, HordHarm, HordMorse)
!%      Subroutines/functions that are needed
!%      read in name of hamiltonian needed
!%      take in name and assign the correct Hord
!%      take in name and assign the correct DHorDz1c
!%      calculate Delta2
!%      calculate total energy that needs to be conserved (classical)
!%      calculate total quantum energy
!%      calculate expectation energy (quantum)
!%
!%      Links with other modules:
!%      Need to link with HamiltonianSpec as that has Hord for each system
!%      Need to link with the integration module as derivs routine calls Hord
!%      Might need to link with BSET to get parameters like ndig
!$
!%********************************************************************************************************************!!

Module GenHam
use BSETMK, only:csbasisfn!, ndim, nbf
use BSETMK, only:ndim,nbf, DTOC! does it get ndim from hamspec? it might be worth investigating as it seems to get it anyway
use BSETMK, only:INBASISFN
!use BSET, only: gam, hbar ! needs them for the Echeck mod
use HamiltonianSpec
!$********************************************************************************************************************!!
!$      Global Variables
!$ Name is a global variable that will read in the name of the system (Free_Particle, Harmonic, Morse)
!$********************************************************************************************************************!!

Character(LEN=20)::NameSys='Notread'
!real(kind=8)::Ecut=0.0d0
!integer ntries
contains


!$********************************************************************************************************************!!
!$                                                      Subroutine ReadSysName
!$      This subroutine will read in the name of the system
!$      It will work by doing a do loop where it searches for System: Name
!$      where System: is a character called Line1 and Name is the actual name of the system
!$********************************************************************************************************************!!


  Subroutine ReadSysName(LINE1)
    Implicit None
    Character(Len=10)::LINE1
    !System is the character that should read System: whilst Name is the name of the system
    integer::ierr
    ierr=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
    read(30,*,iostat=ierr)LINE1 ! unit 30 was chosen randomly
    do while (ierr==0)
       if(LINE1 == "System:") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,NameSys
       end if
       read(30,*, iostat=ierr)LINE1
    end do
    close(30)
  End Subroutine ReadSysName

!$********************************************************************************************************************!!
!$                                                      Function H_ord
!$      This is a function that will take in the name of the system read in the routine above and select the
!$      appropriate H_ord. The advantage of doing it this way is that you can use H_ord in an equation without having
!$      to call a subroutine. 
!$      Obviously have to call ReadSysName first otherwise you won't have a system
!$      I'm going to use IF statements here rather than CASE statements
!$      I've decided this becausse it will mean I don't have to write case(1) when I want to use it 
!$      and I can have an if statement if nothing is read to tell me that
!$      Arguments:
!$      bf1cz and bf2- these are the same as the arguments in Hamiltonian Specific module
!$********************************************************************************************************************!!

  Function H_ord(bf1cz,bf2)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    complex(kind=8)::H_ord
    if (NameSys == 'Notread')then ! Name was initialised to notread 
       Print *, 'havent read in System Name in, have you called ReadSysName? '
    else if (NameSys =='Free') then 
       H_ord=HordFree(bf1cz,bf2)
    else if (NameSys=='Harmonic') then
       H_ord=HordHarm(bf1cz,bf2)
    else if (NameSys=='Morse') then
       H_ord=HordMorse2(bf1cz,bf2)
    else if (NameSys=='Hubbard') then
       H_ord=HordHubb(bf1cz,bf2)     
    else if (.not. (NameSys=='Morse' .OR. NameSys=='Harmonic' .OR. NameSys=='Free')) then
       Print *, 'WARNING, I havent understood name of system!!!!!!!!!!!'
    end if
  End Function H_ord

!$********************************************************************************************************************!!
!$                                                      Function DH_ord_Dz1c
!$      This does the same as the H_ord function but for DHord/Dz1*
!$      Name is global so isn't an argument
!$      I like how each of the functions takes the same two arguments
!$********************************************************************************************************************!!

  Function DH_ord_Dz1c(bf1cz,bf2,i)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    integer,intent(in)::i
    complex(kind=8)::DH_ord_Dz1c
    if (NameSys == 'Notread')then ! Name was initialised to notread 
       Print *, 'havent read in System Name in, have you called ReadSysName? '
    else if (NameSys =='Free') then 
       DH_ord_Dz1c=DHordDz1cFree(bf1cz,bf2,i)
    else if (NameSys=='Harmonic') then
       DH_ord_Dz1c=DHordDz1cHarm(bf1cz,bf2,i)
    else if (NameSys=='Morse') then
       DH_ord_Dz1c= DHordDz1cMorse2(bf1cz,bf2,i) 
    else if (NameSys=='Hubbard') then
       DH_ord_Dz1c=DHordDz1cHubb(bf1cz,bf2,i)     
    else if (.not. (NameSys=='Morse' .OR. NameSys=='Harmonic' .OR. NameSys=='Free')) then
       Print *, 'WARNING, I havent understood name of system!!!!!!!!!!!'
    end if
  End Function DH_ord_Dz1c

!$********************************************************************************************************************!!
!$                                                      Function DHordDz1
!$      This function gives the derivative of z by dz whilst the function above gives it by z*
!$********************************************************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!I don't actually need the function below anymore
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Function DH_ord_Dz1(bf1cz,bf1,i)
!    Implicit none
!    type(csbasisfn),intent(in)::bf1cz,bf1
!    integer,intent(in)::i
!    complex(kind=8)::DH_ord_Dz1
!    if (NameSys == 'Notread')then ! Name was initialised to notread 
!       Print *, 'havent read in System Name in, have you called ReadSysName? '
!    else if (NameSys =='Free') then
!       DH_ord_Dz1=DHordDz1Free(bf1cz,bf1,i)
!    else if (NameSys=='Harmonic') then
!       !Print *, 'havent written this yet'
!       !STOP
!       DH_ord_Dz1=DHordDz1Harm(bf1cz,bf1,i)
!    else if (NameSys=='Morse') then
!       !Print *, 'havent written this yet'
!       !STOP
!    !   DH_ord_Dz1= DHordDz1Morse(bf1cz,bf1,i)
!    else if (.not. (NameSys=='Morse' .OR. NameSys=='Harmonic' .OR. NameSys=='Free')) then
!       Print *, 'WARNING, I havent understood name of system!!!!!!!!!!!'
!    end if
!  End Function DH_ord_Dz1

!$********************************************************************************************************************!!
!$                                                      Function Delta2Hord
!$      This value uses the Hord and DhordDz values calculated above
!$      It will be used when calculating the amplitudes
!$      This function looks a bit more complicated then others because it is
!$      split into two parts. 
!$      The DHordBit does the sum of the DH_ord function over all the degrees of
!$      freedom. This could not be done in the DH_ord function itself as the
!$      DH_ord(zic,zj) must be multiplied by zi-zj then summed over all the
!$      degrees of freedom. For the H_ord function the sum is written in the
!$      function as it is simpler and this doesn't need to be multiplied.  
!$********************************************************************************************************************!!


  Function Delta2Hord(bfi,bfj)
    Implicit none
    type(csbasisfn),intent(in)::bfi,bfj
    type(csbasisfn)::bfic, bfjc
    integer::i,j
    complex(kind=8)::Delta2Hord
    complex(kind=8)::DHordBit
    CALL INBASISFN(bfic)
    CALL INBASISFN(bfjc)
    bfic%z(1:ndim)=conjg(bfi%z(1:ndim))
    bfjc%z(1:ndim)=conjg(bfj%z(1:ndim))
    DHordBit=cmplx(0.0d0,0.0d0)
    !do i=1,ndim  !                      - might need to change this for 2d prob
       do j=1,ndim
          DHordBit=DHordBit+(DH_ord_Dz1c(bfjc,bfj,j)*(bfic%z(j)-bfjc%z(j)))
          !DHordBit=DHordBit+(DH_ord_Dz1c(bfjc,bfj,j)*(bfjc%z(j)-bfic%z(i)))
       end do   
    !end do
    !Delta2Hord=H_ord(bfic,bfi)-H_ord(bfic,bfj)+DHordBit
    Delta2Hord=(H_ord(bfic,bfj)-H_ord(bfjc,bfj)-DHordBit)
    !Print *, 'Delta2Hord', Delta2Hord
    deallocate(bfic%z)
    deallocate(bfjc%z)
    ! The DH_ord_Dz1c should actually be by dz2c so might have to check if this
    ! makes a difference
  End Function Delta2Hord


!       Write Energy functions:

!$********************************************************************************************************************!!
!$                                                      Function Energyclass
!$      This function is to calculate the classical energy for a basis set(bs)
!$      It will take the whole basis set as the argument
!$      and sum over all the degrees of freedom and basis functions
!$      I think the expression for the classical energy is:
!$      classical energy=Sum_i Hord(zic,zi)
!$      I think it doesn't involve zj and there won't be any overlap
!$********************************************************************************************************************!!

  Function Energyclass(bs)
    Implicit none
    type(csbasisfn),dimension(:), intent(in)::bs
    type(csbasisfn)::bf1cz,bf2
    integer::i
    complex(kind=8)::Energyclass ! makes sense for this to be a real number so check if imag is small
    !real(kind=8)::Energyclass
    Energyclass=cmplx(0.0d0,0.0d0)
    call INBASISFN(bf1cz)
    call INBASISFN(bf2)
    do i=1,nbf
       bf2=bs(i)
       bf1cz%z(1:ndim)=conjg(bf2%z(1:ndim))
       Energyclass=Energyclass + (H_ord(bf1cz,bf2))
    end do
    deallocate(bf1cz%z)
    deallocate(bf2%z)
  End Function Energyclass

!  Function EnergyclassBF(bf)
!    Implicit none
!    type(csbasisfn)::bf, bf1c
!    complex(kind=8)::EnergyclassBF
!    call INBASISFN(bf1c)
!    bf1c%s=bf%s
!    bf1c%d=bf%d
!    bf1c%z(1:ndim)=dconjg(bf%z(1:ndim))
!    EnergyclassBF=H_ord(bf1c,bf)
!    deallocate(bf1c%z) 
!  End Function EnergyclassBF

!$********************************************************************************************************************!!
!$                                                      Function Energyquantum
!$      This function calculates the total quantum energy of the basis set (bs)
!$      It will take the whole basis set as the argument, bs also has the
!$      amplitudes in it. 
!$      The expression I will programme is:
!$      Energy=Sum_i Sum_j( D_iconjg * D_j*Hord(zic,zj)
!$********************************************************************************************************************!!

  Function Energyquantum(bs)
    Implicit none
    type(csbasisfn),dimension(:),allocatable,intent(in)::bs
    type(csbasisfn)::bf1cz,bf2
    integer::i,j
    complex(kind=8)::Energyquantum
    call INBASISFN(bf1cz)
    call INBASISFN(bf2)
    Energyquantum=cmplx(0.0d0,0.0d0)
    do i=1,nbf
       bf1cz%z(1:ndim)=conjg(bs(i)%z(1:ndim))
       do j=1,nbf
          bf2=bs(j)
          Energyquantum=Energyquantum + conjg(bs(i)%d)*bs(j)%d*H_ord(bf1cz,bf2)
       end do
    end do
    deallocate(bf1cz%z)
    deallocate(bf2%z)
  End Function Energyquantum

!$********************************************************************************************************************!!
!$                                                      Function Energyexpct
!$      This function calculates the quantum expectation energy of the basis set
!$      It will take the whole bais set (bs) as argument, like the other energy
!$      functions. 
!$      The expression I will programme is:
!$      Total Energy Expectation: Sum_i C_ic * Di * Hord(zic,zi)
!$      This argument needs C so need to call DTOC routine in basis set
!$      I have changed the DTOC routine so that it calls overlap rather than
!$      needing to call overlap in this function too
!$      This keeps the overlap up to date. 
!$********************************************************************************************************************!!
  
  Function Energyexpct(bs)
    ! this function is like the classical one (one do loop) but with amplitudes
    Implicit none
    type(csbasisfn),dimension(:),allocatable,intent(in)::bs
    type(csbasisfn)::bf1cz,bf2
    integer::i
    complex(kind=8),dimension(:),allocatable::c
    complex(kind=8)::Energyexpct
    allocate(c(nbf))
    call INBASISFN(bf1cz)
    call INBASISFN(bf2)
    call DTOC(c,bs)
    Energyexpct=cmplx(0.0d0,0.0d0)
    do i=1,nbf
       bf2=bs(i)
       bf1cz%z(1:ndim)=conjg(bf2%z(1:ndim))
       Energyexpct= Energyexpct + c(i)*bs(i)%d*H_ord(bf1cz,bf2)
    end do
    deallocate(bf1cz%z)
    deallocate(bf2%z)
    deallocate(c)
  End Function Energyexpct

End Module GenHam
