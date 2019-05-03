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
use BSETMK!, only:INBASISFN
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
          !Print *, 'Popmax=', Popmax
          if (ierr .ne. 0) then
             Print *,' error reading Popmax'
          End if
       else if (LINE1=="Popmin") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Popmin
          !Print *, 'Popmin=', Popmin
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

End Module Check
