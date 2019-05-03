!%********************************************************************************************************************!!
!%Author: John Kattirtzi, Date: 03/06/09                                         
!%                                                   Module: Outmod
!%      This is a module that will have subroutines written which will
!%      determine the outputs that are produced. 
!%      The outputs will be written in the current directory for now but
!%      this will probably be changed later 
!%      It would probably be a good idea to have an inout file to keep
!%      input and output files      
!%      Subroutines required:
!%      Subroutine to read dn- determines the number of lines produced in output file
!%      Subroutine to read option of trajectories or conservation(energy+norm)      
!%      Subroutine to print conservation
!%      Subroutine to print norm
!%      You will need an output for autocorrelation function but worry
!%      about this when you do that!      
!%********************************************************************************************************************!!

 Module Outmod
   USE BSET2
   USE GenHam
   USE ACFMOD
   integer::ierrM=0
!   Character(LEN=10)::Op1
!   Character(LEN=10)::Op2
    Character(LEN=20)::Trajname='No.Tout.dat'
    Character(LEN=20)::Consvname='No.Cout.dat'
    character(LEN=20)::Outpname='No.outp.dat'
    character(LEN=20)::bfoutname='No.bfout.dat'
   contains

!$********************************************************************************************************************!!
!$                                              Subroutine ReadDn
!$  This subroutine will read Dn
!$  Dn is the incrememnet in the do loop for the output of the Traj and Consv file
!$  It must be an integer
!$  A larger value of Dn will mean that FEWER lines are written as the increment is larger
!$********************************************************************************************************************!!

     Subroutine ReadDn(Dn)
       Implicit None
       integer,intent(out)::Dn
       character(LEN=10)::LINE1
       OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierrM)
       IF (ierrM .ne. 0) then
          Print *, 'error opening inham in ReadDn'
          STOP
       End IF
       read(30,*,iostat=ierrM)LINE1 ! unit 30 was chosen randomly
       do while (ierrM==0)
          if(LINE1 == "dn") then
             backspace(30)
             read(30,*,iostat=ierrM)LINE1,dn
             IF (ierrM .ne. 0) then
                Print *, 'error reading Dn in ReadDN'
                STOP
             END IF
          end if
          read(30,*,iostat=ierrM)LINE1
       end do
       close(30)
     End Subroutine ReadDn

!$********************************************************************************************************************!!
!$                                              Subroutine ReadOutOp
!$  
!$  This subroutine will read the names of the Traj and Consv file
!$  If these are not specified in the inham.dat file then the ones set at the beginning of the module will be used.
!$  This will mean that no Traj and Consv file will be produced.
!$  The second part of this subroutine calls the Head routine which will write the header or give a warning if no output
!$  The advantage of having this part here is that these routines need to be called once only as does the first part. 
!$********************************************************************************************************************!!

     Subroutine ReadOutOp(LINE1)
       Implicit None
       Character(LEN=10)::LINE1
       OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierrM)
       IF (ierrM .ne. 0) then
          Print *, 'error opening inham in ReadOutOp'
          STOP
       End IF
       read(30,*,iostat=ierrM)LINE1 ! unit 30 was chosen randomly
       do while (ierrM==0)
          if(LINE1 == "Traj") then
             backspace(30)
             read(30,*,iostat=ierrM)LINE1,Trajname
             !Print *, LINE1, Trajname
             IF (ierrM .ne. 0) then
                Print *, 'error reading Trajname in ReadOutOp'
                STOP
             END IF
          else if (LINE1 == "Consv") then
             backspace(30)
             read(30,*,iostat=ierrM)LINE1,Consvname
          else if (LINE1 == "Bfout") then
             backspace(30)
             read(30,*,iostat=ierrM)LINE1,bfoutname
             !Print *, LINE1, Consvname
          End if
          read(30,*,iostat=ierrM)LINE1
       end do

      ! part 2: calling the header and warning if no output
       IF (.not. (Trajname=="No.Tout.dat")) then
               call TrajHead
       END IF
       IF (.not. (Consvname=="No.Cout.dat")) then
               call ConsvHead
       END IF
       IF (.not. (bfoutname=='No.bfout.dat')) then
               call bfoutHead
       END IF        

     !  If (Consvname=="No.Cout.dat" .AND. Trajname=="No.Tout.dat")then
     !     !Print *, 'Im giving you a warning but change if I shouldnt!!!'
     !     Print *, "Traj =", Trajname, "Consv", Consvname
     !     Print *, 'WARNING, no output'
     !  End IF
                                                                
       close(30)
     End Subroutine ReadOutOp

!$********************************************************************************************************************!!
!$                                              Subroutine GenOup
!$  This subroutine is called from the main program
!$  It calls the other subroutines which write the output lines
!$********************************************************************************************************************!!

     Subroutine GenOutp(bs,t,c)
       Implicit None
       type(csbasisfn),dimension(:),intent(in)::bs
       real(kind=8),intent(in)::t
       complex(kind=8),dimension(:)::c
       IF (.not. (Trajname=="No.Tout.dat")) then
          Call WriteTraj(bs,t,c)
       END IF
       IF (.not. (Consvname=="No.Cout.dat")) then
          call WriteConsv(bs,t)
       END IF 
       IF (.not. (bfoutname=="No.bfout.dat")) then
          call Writebfout(bs,t)
       END IF
     End Subroutine GenOutp

!$********************************************************************************************************************!!
!$                                              Subroutine Writebfout
!$      This subroutine bfout will write out all components of the z basis function
!$      i.e: time, z, a =dexp(iS)
!$********************************************************************************************************************!!
     Subroutine Writebfout(bs,t)
       Implicit None
       type(csbasisfn),dimension(:),intent(in)::bs
       real(kind=8),intent(in)::t
       integer::i,nbfhalf,j
       nbfhalf=nbf/2
       j=1+nbfhalf
       Open(UNIT=33, FILE=bfoutname, STATUS='OLD', ACCESS='APPEND')
        Write(33,*),t,'1',Real(bs(1)%z), AIMAG(bs(1)%z),Real(bs(1)%d*exp(bs(1)%s*cmplx(0.0d0,1.0d0))) ,&
             AIMAG(bs(1)%d*exp(bs(1)%s*cmplx(0.0d0,1.0d0)))
        Write(33,*),t,j,Real(bs(j)%z),AIMAG(bs(j)%z),Real(bs(j)%d*exp(bs(j)%s*cmplx(0.0d0,1.0d0))),&
              AIMAG(bs(j)%d*exp(bs(j)%s*cmplx(0.0d0,1.0d0)))
       close(33)
     End Subroutine Writebfout
!$********************************************************************************************************************!!
!$                                              Subroutine WriteTraj
!$  This subroutine writes out 1 line for the trajectory at a specific point in time
!$********************************************************************************************************************!!

     Subroutine WriteTraj(bs,t,c)
       Implicit None
       type(csbasisfn),dimension(:),intent(in)::bs
       real(kind=8),intent(in)::t
       integer::i
       complex(kind=8),dimension(:)::c
       Open(UNIT=31, FILE=Trajname, STATUS='OLD', ACCESS='APPEND')
      ! IF (symop == 'YES') then
      !    do i=1,nbf
      !        IF (i .eq.1 .OR. i .eq. 2) then
      !            Write(31,*),t,i,Real(bs(i)%z), AIMAG(bs(i)%z)
      !        END IF
      !    end do
      ! ELSE
        do i=1,nbf
           Write(31,*),i,t,Real(bs(i)%z), AIMAG(bs(i)%z), Real(bs(i)%d), AIMAG(bs(i)%d), bs(i)%s!, Real(c(i)), AIMAG(c(i))
        end do
       !END IF
        Close(31)
      End Subroutine WriteTraj
    
!$********************************************************************************************************************!!
!$                                              Subroutine TrajHead
!$  This subroutines writes a header for the traj output file
!$********************************************************************************************************************!!

      Subroutine TrajHead
        Implicit None
        Open(UNIT=31, FILE=Trajname, STATUS='NEW')
        Write(31,*), '#********************Trajectory Output File********************'
        write(31,*), '#**************************************************************'
        Write(31,*), '#The trajectory is written as:'
        Write(31,*), '# bf index,time value, real z(ndim), imag z(ndim), real d, imag d, action s'
        Write(31,*),'#'
        close(31)
      End Subroutine TrajHead

!$********************************************************************************************************************!!
!$                                              Subroutine ConsvHead
!$  This subroutines writes a header for the traj output file
!$********************************************************************************************************************!!

     Subroutine ConsvHead
        Implicit None
        Open(unit=32, File=Consvname, STATUS='New')
        Write(32,*), '********************Conservation Output File********************'
        write(32,*), '****************************************************************'
        Write(32,*), 'The classical energy of the basis set and the norm are written out against time:'
        Write(32,*),''
        close(32)
      End Subroutine ConsvHead
      
!$********************************************************************************************************************!!
!$                                              Subroutine bfoutHead
!$      This subroutine will write the header for the bfout output file
!$********************************************************************************************************************!!
    
      Subroutine bfoutHead
        Implicit None
        Open(unit=33, File=bfoutname, STATUS='New')
        Write(33,*), '#******************Basis Function Output File*******************'
        write(33,*), '#***************************************************************'
        Write(33,*), '# The columns printed are:'
        write(33,*), '# Time, basis function index(i),REAL(z(i)), IMAG(z(i)), REAL(d(i)), IMAG(d(i)),s(i), phase'  
        Write(33,*),''
        close(33)
      End Subroutine bfoutHead


!$********************************************************************************************************************!!
!$                                              Subroutine WriteConsv
!$  This subroutines writes out the norm and the classical energy which should be conserved
!$  The classical energy is for the whole basis set 
!$********************************************************************************************************************!!

      Subroutine WriteConsv(bs,t)
        Implicit None
        type(csbasisfn),dimension(:),intent(in)::bs
        real(kind=8),intent(in)::t
        real(kind=8)::norm
        complex(kind=8)::classical
        call Norm2(bs,norm)
        classical=Energyclass(bs)
        Open(unit=32,File=Consvname, STATUS='OLD', ACCESS='Append')
        Write(32,*),t,classical, norm
      End Subroutine WriteConsv

      Subroutine ParName(LINE1,bs,t0, tmax, dt,norm0,class,expt, qttotal )
        Implicit none
        character(len=10)::LINE1
        !real(kind=8),dimension(:)::meanp, meanq
        !real(kind=8)::sigp,sigq
        real(kind=8):: t0, tmax, dt, norm0
        complex(kind=8)::class, expt, qttotal
        integer::i
        type(csbasisfn),dimension(:)::bs
        integer::ierrP, ierr
        !integer(kind=8)::myseed
        ierr=0
        ierrP=0
        !ierrM=0
        !Print *, 'in parname'
        OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierrP)
        rewind(30)
        !Print *, 'opened inham'
        IF (ierrP .ne. 0) then
           Print *, 'error opening inham in ReadOutOp'
           STOP
        End IF
        !LINE1='Nothing!'
        !read(30,*,iostat=ierrP)LINE1 
        !     Print *, 'LINE1 in OUTMOD::', LINE1
        !     IF (ierr .ne. 0) then
        !        Print *, 'error pre loop in par name', ierr
        !        STOP
        !     End IF
        do while (ierrP==0)
           if(LINE1 == "OutPname") then
             backspace(30)
             read(30,*,iostat=ierrP)LINE1,Outpname
             !Print *, LINE1, Outpname
             IF (ierrP .ne. 0) then
                Print *, 'error reading Outpname'
                STOP
             END IF
          end if
          read(30,*,iostat=ierrP)LINE1
       end do
       close(30)
       IF (Outpname .ne.'No.outp.dat') then
           call writeOutp(bs, t0, tmax, dt,norm0, class,expt, qttotal)
           !!Print *, 'calling write out'
        END IF
     end Subroutine ParName

     Subroutine WriteOutp(bs, t0, tmax, dt,norm0, class, expt, qttotal)
       Implicit none
!       real(kind=8),dimension(:)::meanp, meanq
!       real(kind=8)::sigp,sigq
       real(kind=8):: t0, tmax, dt, norm0
       complex(kind=8)::class, expt, qttotal
       integer::i
!       integer(kind=8)::myseed
       type(csbasisfn),dimension(:)::bs 
       OPEN(UNIT=50, FILE=Outpname, STATUS='NEW')
       write(50,*),'# This is a general output file'
        write(50,*)'it wil contain the initial conditions of the basis set '
        write(50,*), 'Parameters from input.dat'
        write(50,*), 'Number of degrees of freedom/dimensions', ndim
        write(50,*), 'Number of basis functions', nbf
        write(50,*), 'units'
        write(50,*), 'hbar', hbar, 'm', m, 'angular frequency', w, 'gamma', gam
        write(50,*), 'initial wavefunction'
        !write(50,*), 'mean p', meanp, 'mean q', meanq, 'sigp', sigp, 'sigq',sigq 
        !write(50,*), 'seed', myseed
        write(50,*), 'Parameters from inham.dat'
        write(50,*), 'System:',NameSys
        write(50,*), 'Hamiltonian Parameters'
        IF (NameSys =="Morse")then
           write(50,*), 'Morse D', Dh, 'Morse alpha', ah
        ELSE IF (NameSys =="Harmonic")then
           write(50,*), 'Force Constant', kh 
        ELSE IF (NameSys =="Hubbard") then
           write(50,*) 'Hubbard: g', g, 'nu', nu, 'eps', eps
        END IF 
        write(50,*) 'Initial Time', t0
        write(50,*), 'Max Time', tmax
        write(50,*), 'Time step', dt
        write(50,*), 'Initial basis set'
        write(50,*), 'Initial z'
        do i=1,nbf
           write(50,*), 'basis function', i, 'z',  bs(i)%z(1:ndim)
        end do
        !do i=1,nbf
        !    write(50,*), bs(i)
        !end do
        write(50, *) 'Initial d'
        do i=1,nbf
           write(50,*), 'basis function', i, 'd', bs(i)%d
        end do   
        write(50,*), 'Initial Calculated Properties:'
        write(50,*), 'Initial norm', norm0
        write(50,*), 'Initial classical energy', class
        write(50,*), 'Initial quantum total energy', qttotal
        write(50,*), 'Initial quantum expectation energy',expt 
       close(50)
     End Subroutine WriteOutp
    End Module Outmod
