!%********************************************************************************************************************!!
!%                                                      Program TestEverything
!%      This is a program this is written to test the modules written
!%      It will test:
!%      Initialisation of the basis set
!%      Running the trajectory (Integration)
!%      If the classical energy is being conserved
!%      If the norm is being conserved
!%********************************************************************************************************************!!             


Program TestEverything
use BSET2
use HamiltonianSpec                                          
use GenHam                                                   ! namesys is declared here   
use Integration
use OutMod
use ACFMOD
!use Check
use Hubbard
Implicit none
type(csbasisfn),dimension(:),allocatable::bs                 ! basis set
character(LEN=100)::LINE1                                    ! character variable that is read in
character(LEN=100)::LINE2
! parameters for reading in z (except for LINE1 above):
integer(kind=8):: myseed
integer::outp
!real(kind=8):: sigp,sigq
real(kind=8),dimension(:), allocatable::mup, muq
!real(kind=8),external::ZBQLNOR
! extra variables for initialising d
complex(kind=8), dimension(:),allocatable::c
complex(kind=8), dimension(:),allocatable::zp
complex(kind=8), dimension(:,:),allocatable::omega
complex(kind=8)::psi
! variable for norm
real(kind=8)::norm
! variables for hamiltonian tests
complex(kind=8):: classical
complex(kind=8):: quantum
complex(kind=8):: total
! variables for integration test
real(kind=8)::t0
real(kind=8)::dt
real(kind=8)::tmax
real(kind=8)::t
type(csbasisfn),dimension(:),allocatable::bsout
type(csbasisfn),dimension(:),allocatable::dbsdt
integer::i
! variable for output
real(kind=8)::t2
integer::n
integer::j,ierr
integer::dn
! variable for ACF
complex(kind=8)::ACFval
integer::jcount
real(kind=8) :: starttime,stoptime
Print *, 'Hello, programme is starting'
call CPU_TIME(starttime)
!$********************************************************************************************************************!!
!$                                                      Output file
!$      Read in the values of dn
!$      dn refers to the increment you want it printed
!$********************************************************************************************************************!!
!ierr=0
Call ReadDn(dn)
call ReadOutOp(LINE1)
!call TrajHead
!call ConsvHead
!$********************************************************************************************************************!!
!$                                                      Initialise basis set
!$      Read in the input file
!$      Allocate the memory
!$      Initialise z,s, d -call initialisation routines
!$********************************************************************************************************************!!

!call READBSPAR(LINE1)! reads in n* dimensions(ndim), n* basis functions(nbf) etc

! Allocate the memory for the array for the mean values for p and q- need to do this before reading them into the array
!allocate(mup(ndim))
!allocate(muq(ndim))
!allocate(c(nbf))
!allocate(gam(ndim))
!call alcbasis(bs)! allocates the memory for the basis set array but also for the z array in the basis set and for gamma!
! Read in mean values for p and q -used for propagated wavefunction
!call READARR(LINE1)  ! reads in nbf parameters: gamma - could read in mass and w too 
!call Readz(myseed,outp,mup,muq,sigp,sigq,LINE1)
!call ReadChk(Line1)

!Allocate the memory for the basis set before initialising 
! allocate memory for overlap omega and calculate it- need it to initialise D
  
! initialise S(action) and D(amplitude),-arrays allocated in subroutines
Print *, 'Read In Basis Set'
call ReadBSIN(bs,c)
call ALCPAR
call BSGAMCHK
!call IntZ(myseed,outp,muq,sigq,mup,sigp,bs,ZBQLNOR)
!call ReadEcut(Line1)  ! reads cut off point for energy
!$********************************************************************************************************************!!
!$                                                      Test Hamiltonian module
!$    Read Hamin.dat file (need to change hamiltonian module to use gamma from basis set)
!$    Test by calculating energy of system
!$********************************************************************************************************************!!
CALL ALCHamPar
CALL ReadSysName(LINE1)
Print *, 'Namesys ', NameSys
CALL  ReadHamSpec(LINE1, NameSys)
If (NameSys=='Harmonic') then
    Print *, 'checking harmonic parameters'
    call HarmChk
End If
!********************************************************************************************************************!!
!  Date: 11.06.08
!  I need to redistribute things for Echeck- it is calling H_ord without knowing which Hord to call!
!********************************************************************************************************************!!
call IntS(bs)
!Print *, 'Starting IntD'
!call IntD(c, mup, muq,bs)
!Print *, 'IntD', bs%d
!call CheckInt(bs,mup,muq,sigp,sigq, ZBQLNOR) 
!call IntD(c, mup, muq,bs) ! need to initialise d again IFF d has changed in checks- probably could find a better way to do this, think about it
call CTODin(c,bs)  
call Norml(c,norm,bs)
Print *,'Initial norm',  norm
call Norm2(bs,norm)
Print *, 'initial norm2', norm
!*********************************************************************************************************************!!
!                                            Energy Testing
!*********************************************************************************************************************!!
classical=Energyclass(bs)
quantum= Energyexpct(bs)
total= Energyquantum(bs)
Print *,'*************************************************************************'
Print *, 'classical energy ', classical
Print *, 'quantum expectation', quantum
Print *, 'quantum total', total
Print *,'*************************************************************************'
!$********************************************************************************************************************!!
!$                                                      Integration part testing
!$      Read in the time parameters
!$      write a do loop to calculate the trajectory
!$      In the do loop print off:
!$      the classical energy
!$      the norm
!$      maybe also look at the quantum and expectation energy
!$********************************************************************************************************************!!
call ReadTime(Line1,t0,tmax,dt)
Print *, 'number of basis functions', nbf, 'dimensions', ndim
Print *, 'initial time', t0, 'tmax', tmax, 'dt', dt
Print *,'*************************************************************************'
!do i=1,nbf
!   Print *, 'bf', i, 'initial z',bs(i)%z(1:ndim)
!end do
call alcbasis(bsout)
call alcbasis(dbsdt)
t=t0
t2=t
bsout(:)=bs(:)
n=0
!STOP
!$********************************************************************************************************************!!
!$                                              ACF testing
!$      Need to read in if ACF should be stored
!$      If it should be stored then call the routine to store it and print in an
!$      output file. 
!$      After have done this write a seperate program that will read in the ACF
!$      output and do the FT.  
!$********************************************************************************************************************!!
!Print *, 'in prog acf: ', LINE1
!call ACFHead
call ReadPopOp(Line1)
call ReadACF(LINE1,mup,muq)
call ACFHead
!Print *, 'in prog, LINE1::', Line1
! testing the general write out file
!call WriteOutp(bs,mup, muq, sigp,sigq, t0, tmax, dt,norm, classical, quantum, total)
call ParName(Line1,bs, t0, tmax, dt,norm,classical,quantum, total )
jcount=1
Print *, 'starting main do loop'
     do while (ABS(t)<=ABS(tmax))
         !  open(unit=20,FILE='Tproga.dat',STATUS='OLD', ACCESS='APPEND')
         !  do i=1,nbf
         !  write(20,*)t, Real(bs(i)%z), AIMAG(bs(i)%z)
         !  end do
        IF (dn ==0) then
             call Genoutp(bs,t,c)
        ELSE IF (NINT(t)==n) then
           !Writes the output
           call GenOutp(bs,t,c)
           n=n+dn
        End IF
        !        Print *, 'in main prog calling derivs' 
        ! ACFval and j need to be sorted/declared and j has to be given a value
        call Dim2Npop(t,bs, c)
        IF (.not. (ACFname=="NoACFoutin.dat")) then
            call ACFrtn(bs,mup,muq,ACFval,jcount,t)
        END IF    
        jcount=jcount+1
        call derivs(t,bs,dbsdt)
        !Print *, 'after derivs'
        !classical=Energyclass(dbsdt)
        !  Print *, 'energy of derivs', classical
        call RK4ccs(bs,dbsdt,t,dt,bsout)
       ! STOP
        ! do i=1,nbf
        !     Print *, ' traj', bsout(i)%z(1:ndim)
        ! end do
        
        bs=bsout
        !classical=Energyclass(bsout)
        !Print *, 'classical bsout', classical 
        t=t+dt
     end do
     jcount=jcount-1
     Print *, 'j', jcount
     call ACFTail(jcount)
     !Print *, 'w', w
     deallocate(bsout)
     deallocate(dbsdt)
          !IF (size(muq) > 0 .AND. size(mup) > 0) then
     IF (allocated(muq)) then
       deallocate(muq)
     end if
     If (allocated(mup)) then
      deallocate(mup)
     end if       
    !      Print *, 'deallocating mup, muq', size(mup), size(muq)        
    ! deallocate(mup)
    ! deallocate(muq)
    ! end if
     deallocate(c)
!$********************************************************************************************************************!!
!$                                                      End of test
!$********************************************************************************************************************!!
    call CPU_TIME(stoptime)
    Print *,  '************************************************************************'
    Print *,  '************************************************************************'
    Print *, 'Time taken:', (stoptime-starttime)/60.0d0, 'mins'
    Print *, 'Outputs:',Trajname, Consvname, Outpname , popname, bfoutname
    Print *, 'program terminated ok'
    Print *,  '************************************************************************'
    Print *,  '************************************************************************'
   End Program TestEverything
