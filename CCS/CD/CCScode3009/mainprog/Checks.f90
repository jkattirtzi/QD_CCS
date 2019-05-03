!%********************************************************************************************************************!!
!% Author: John Kattirtzi, Date: 11/06/09      
!%                                    Module Checks
!%
!%      This module is written to check any ares of the programme
!%      particular focus should be on the basis set. 
!%      There are only two subroutines at the moment:
!%      1. Read in a cut off energy value and number of tries for routine to do this
!%      2. Reinitialise z so it is lower than cutoff
!%
!%      Links with other modules:
!%      Need to link with GenHam for Hamiltonian
!%      Might need to link with BSET to get parameters like ndim and also derived type
!%
!%********************************************************************************************************************!!

Module Check
use BSET2!, only:INBASISFN
!use BSET, only: gam, hbar ! needs them for the Echeck mod
!use BSET, ty
use GenHam, only:H_ord
use GenHam, only:NameSys
use Hubbard, only:Npopbf
!$********************************************************************************************************************!!
!$      Global Variables
!$ Ecut is the cut off value for the energy of a basis function
!$ ntries is the number of tries the subroutine has to get this
!$********************************************************************************************************************!!

real(kind=8)::Ecut=0.0d0
real(kind=8)::Emin=0.0d0
real(kind=8)::Popmax=0.0d0
real(kind=8)::Popmin=0.0d0
character(LEN=10)::PCheck="NO"
integer ntries
contains

!$********************************************************************************************************************!!
!$                                              Subroutine ReadEcut
!!! This subroutine does not belong here and neither does the one below, didn't know where else to put them
!!! Worry about this later when you have more time
!$  This subroutine reads in the Ecut value and the ntries value.
!$  Ecut is a real number
!$  It is a parameter that determines the greatest energy for a basis function
!$  ntries is the number of tries per basis function the subroutine tries to get an energy for a basis function
!$  lower than Ecut
!$********************************************************************************************************************!!

  Subroutine ReadEcut(LINE1)
    Implicit None
    Character(LEN=10)::LINE1
    integer::ierr
    ierr=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
    read(30,*,iostat=ierr)LINE1 
    do while (ierr==0)
       if(LINE1 == "Ebfmax") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Ecut
          if (ierr .ne. 0) then
             Print *,' error reading Ebfmax'
          End if
       else if (LINE1=="Ebfmin") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Emin
          if (ierr .ne. 0) then
             Print *, 'error reading Ebfmin'
          End if
       else if (LINE1=="PopChk") then  
         backspace(30)
        read(30,*,iostat=ierr), LINE1, PCheck
        if (ierr .ne.0) then
           Print *, 'error reading PopChk'
        End IF
       else if(LINE1 == "Popmax") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Popmax
          Print *, 'Popmax=', Popmax
          if (ierr .ne. 0) then
             Print *,' error reading Popmax'
          End if
       else if (LINE1=="Popmin") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Popmin
          Print *, 'Popmin=', Popmin
          if (ierr .ne. 0) then
             Print *, 'error reading Popmin'
          End if

       else if(LINE1== "Ntries") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,ntries
          if (ierr .ne. 0) then
              Print *,' error reading ntries'
          End if
       end if
       read(30,*, iostat=ierr)LINE1
    end do
    If ((Ecut .lt. Emin) .OR. (Popmax .lt. Popmin)) then
       Print *, 'Emax is lt Emin OR Popmax is lt Popmin'
       Print *, 'Emax:', Ecut, 'Emin:', Emin, 'Popmax:', Popmax, 'Popmin:', Popmin
       STOP
    END IF   
    close(30)
  End Subroutine ReadEcut


!$********************************************************************************************************************!!
!$                                              Subroutine Echeck (now CheckInt)
!$  This subroutine calculates the Hord(bf1*,bf1) for a basis function.
!$  If it is higher than a cut off point it produces a different basis function to replace it. 
!$  ntries is the number in the input file
!$  nbfcount is the number of basis functions lower than Ecut
!$  nout is an integer to get out of the do loop
!$  ntimes counts the number of times the basis function has been intialised
!$  It's max number can be ntries. 
!$  ntries is per basis function 
!$********************************************************************************************************************!!


  Subroutine CheckInt(bs,mup,muq,sigp,sigq, ZBQLNOR)
    Implicit None
    type(csbasisfn),dimension(:)::bs
    type(csbasisfn)::bfcz,bf
    real(kind=8),intent(in),dimension(:),allocatable::mup
    real(kind=8),intent(in),dimension(:),allocatable::muq
    real(kind=8),intent(in)::sigp
    real(kind=8),intent(in)::sigq
    real(kind=8),external::ZBQLNOR
    integer::i
    integer::j, k
    real(kind=8)::Totpop, Pop1, Pop2
    complex(kind=8)::EHord
    integer::nbfcount, ntimes, nout
    integer::counter
    call INBASISFN(bfcz)
    call INBASISFN(bf)
    IF (Ecut .eq. 0 .OR. Emin .eq. 0) then
               Print *, 'Ebfmax and Ebfmin needs to be specified in inham'
 !               STOP
 !   END IFa
    ELSE IF (Ecut .ne. 0 .AND. Emin .ne. 0) then
    i=1
    nbfcount=1
    do while (nbfcount .le. nbf)
       nout=0
       bf=bs(i)
       bfcz%z(1:ndim)=dconjg(bf%z(1:ndim))
       ntimes=0
       do while (nout .eq. 0)
           ntimes=ntimes+1
           IF (ntimes .gt. ntries) then
                Print *, 'ntimes greater than ntries'
                Print *, 'Im giving up here'
                STOP
                nout=1
            END IF
            call Popcheck(bs,bfcz,bf,nbfcount, ntimes, i,Totpop, Pop1, Pop2, muq,mup, sigp, sigq, ZBQLNOR)
            call Echeck2(bs,bfcz,bf,nbfcount, ntimes, i,EHord, muq,mup, sigp, sigq, ZBQLNOR)
            call Getout(nout,nbfcount,Totpop, Pop1, Pop2, EHord )!- need to write a subroutine that will check if
            !    can accept the basis function and get out of the loop
            IF (.not.((nout .eq. 0 ).OR. (nout .eq. 1))) then
                Print *, 'Error with nout, nout=', nout 
            END IF
       end do
    nbfcount=nbfcount+1
    i=i+1
    end do
    END IF
    deallocate(bf%z)
    deallocate(bfcz%z)
  End Subroutine CheckInt

  Subroutine Getout(nout,nbfcount,Totpop, Pop1,Pop2, EHord)
    Implicit None
    integer,intent(out)::nout
    integer,intent(in)::nbfcount
    real(kind=8),intent(in)::Totpop, Pop1, Pop2
    complex(kind=8),intent(in)::EHord
    IF ((NameSys=="Hubbard") .AND. (.not.( PCheck=="NO" ))) then
       ! Print *, 'In Getout in the Hubbard pop one, Popmin', Popmin, 'Popmax',Popmax
       IF (Popmax .gt. Totpop .AND. Popmin .lt. Totpop .AND. Ecut .gt. Real(EHord) .AND. Emin .lt. Real(EHord)) then
          Print *, 'Accepted bf for Energy AND Pop'
          Print *, 'Using bf', nbfcount, 'Total Pop', Totpop, 'pop1', Pop1, 'pop2', Pop2,'Ehord', Real(Ehord)
          nout=1
       End IF
    ELSE 
       !Print *, 'Not in the pop get out'     
       IF (Ecut .gt. Real(EHord) .AND. Emin .lt. Real(EHord)) then
          Print *, 'Accepted bf for Energy'
          Print *, 'Using bf', nbfcount, 'Ehord', Real(Ehord)
          nout=1
       End IF
    End IF
  End Subroutine Getout
    
  Subroutine Popcheck(bs,bfcz,bf, nbfcount, ntimes,i,Totpop, Pop1, Pop2, muq,mup, sigp, sigq, ZBQLNOR)
    Implicit None
    type(csbasisfn),dimension(:)::bs
    type(csbasisfn)::bfcz,bf
    integer:: i, j, nbfcount, ntimes
    real(kind=8),intent(in),dimension(:),allocatable::mup
    real(kind=8),intent(in),dimension(:),allocatable::muq
    real(kind=8),intent(in)::sigp
    real(kind=8),intent(in)::sigq
    real(kind=8),external::ZBQLNOR
    real(kind=8)::Totpop, Pop1, Pop2
    integer::dim1=1
    integer::dim2=2
    complex(kind=8),dimension(:),allocatable::c
    allocate(c(nbf))
    IF (NameSys=="Hubbard") then
       Print *, 'Checking int pops'
       Pop1=Npopbf(bf,dim1)
       Pop2=Npopbf(bf,dim2)
       Totpop=Pop1+Pop2
       IF ((Totpop .gt. Popmax) .OR. (Totpop .lt. Popmin)) then
          Print *, 'Rejecting Bf because of pop, ntime=', ntimes, 'for nbf=', nbfcount, 'Totpop', Totpop
          Print *, 'Popmax=', Popmax, 'Popmin=', Popmin, 'ntries=', ntries
          do j=1,ndim
             bf%z(j)=cmplx(SQRT(gam(j)/2.0d0)*(ZBQLNOR(muq(j),sigq))&
                  ,(1.0d0/hbar*SQRT(1.0d0/(2.0d0*gam(j)))*ZBQLNOR(mup(j),sigp)))
          !call IntD(c, mup, muq, bs)
          end do
          bfcz%z(1:ndim)=dconjg(bf%z(1:ndim))
          bs(i)%z(1:ndim)=bf%z(1:ndim)
       End IF
    End IF
    deallocate(c)
 End Subroutine Popcheck

  Subroutine Echeck2(bs,bfcz,bf, nbfcount, ntimes,i,EHord, muq,mup, sigp, sigq, ZBQLNOR)
    Implicit None
    type(csbasisfn),dimension(:)::bs
    type(csbasisfn)::bfcz,bf
    integer:: i, j, nbfcount, ntimes
    real(kind=8),intent(in),dimension(:),allocatable::mup
    real(kind=8),intent(in),dimension(:),allocatable::muq
    real(kind=8),intent(in)::sigp
    real(kind=8),intent(in)::sigq
    real(kind=8),external::ZBQLNOR
    complex(kind=8)::EHord
    EHord=H_ord(bfcz,bf)
    !Print *, 'In Echecks2'
    IF ((Ecut .le. Real(EHord)) .OR. (Emin .gt. Real(EHord))) then
       Print *, 'Rejecting bf because of energy, ntime=', ntimes, 'for nbf=', nbfcount, 'EHord', EHord
       Print *, 'Emax=', Ecut, 'Emin=', Emin, 'ntries=', ntries
       do j=1,ndim
          bf%z(j)=cmplx(SQRT(gam(j)/2.0d0)*(ZBQLNOR(muq(j),sigq))&
               ,(1.0d0/hbar*SQRT(1.0d0/(2.0d0*gam(j)))*ZBQLNOR(mup(j),sigp)))
       end do
       bfcz%z(1:ndim)=dconjg(bf%z(1:ndim))
       bs(i)%z(1:ndim)=bf%z(1:ndim)
    ELSE
    !  Print *, 'bf fine energy wise'   
    END IF
  End Subroutine Echeck2


End Module Check
