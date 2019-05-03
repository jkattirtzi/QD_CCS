!%********************************************************************************************************************!!
!%                                              Module FT1
!%     Author: John Kattirtzi, Date: 31/07/09
!%     Version:2
!!     Note: This is a new version because the other one was hard to understand.
!%  This module deals with the Fourier Transformation of the Articulation Function.
!%  It will be used for a main program, that isn't part of the CCS code
!%  It will require:
!%  A subroutine to read the inham input file
!%  A subroutine to read the ACF input file
!%  A function to the do FT
!%  A function for the window function before the FT
!%  Outputs for the FT and the decay output
!%********************************************************************************************************************!!


Module FT1
character(LEN=5)::Wdwf='NO' !Winop
real(kind=8)::Emax
real(kind=8)::E0
character(LEN=20)::decayname='No.decay.out'
character(LEN=20)::FTname='No.ft.out'
character(LEN=20)::ACFname
real(kind=8)::dt
real(kind=8)::Tcut
integer(kind=8)::jmax
complex(kind=8)::Img=cmplx(0.0d0,1.0d0)
contains

!$********************************************************************************************************************!!
!$                                              Subroutine ReadFT1
!$  This is the first subroutine that reads in variables for the FT.
!$  It will read in from the inham.dat - might change this name.
!$  It will read in:
!$  Winop (Wdwf)- the option of a window function (Wdwf=No or not no)
!$  FTEmax(Emax)- the max value of E 
!$  FTE0 (E0)- the min value of E
!$  DecOut(decayname)- the name of the decay output file- ACF*window
!$  FTOut(FTname)- a name for the FTout file
!$  ACF (ACFname)-a name for the ACF input file (out from CCS code)
!$  dt (dt) -time step
!$  FTcut (Tcut) - cut off time
!$  MAXJ - (Jmax)- max number of lines in ACF
!$********************************************************************************************************************!!
  Subroutine ReadFT1(Line1)
    Implicit None
    Character(LEN=10)::LINE1
    integer::ierr
    ierr=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
    read(30,*,iostat=ierr)LINE1
    do while (ierr==0)
       if(LINE1=="Winop")then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Wdwf
          if (ierr .ne. 0) then
             Print *, 'error reading Wdwf'
          end if
      else if(LINE1=="FTEmax") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Emax
          if(ierr .ne. 0) then
             Print *, 'error reading ETmax'
          end if     
       else if (LINE1=="FTE0") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,E0
          if (ierr .ne. 0) then
             Print *, 'error reading E0 in ReadFT'
             STOP
          end if
      else if (LINE1=="DecOut") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,decayname
          if (ierr .ne. 0) then
              Print *, 'error reading ACFdecay'
          end if 
       else if(LINE1 == "FTOut") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,FTname
          if (ierr .ne. 0) then
             Print *, 'error reading FTOP in ReadFT'
             STOP
          end if
       else if(LINE1 == "ACF") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,ACFname
          if (ierr .ne. 0) then
             Print *, 'error reading FTOP in ReadFT'
             STOP
          end if
       else if (LINE1=="dt") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,dt
          if (ierr .ne. 0) then
                Print *, 'error reading dt in ReadFT'
                STOP
          end if
       else if (LINE1=="FTcut") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,Tcut
          if (ierr .ne. 0) then
             Print *, 'error reading FTcut in ReadFT'
             STOP
          end if
       else if (LINE1=="MAXJ") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,jmax
          if (ierr .ne. 0) then
               Print *, 'error reading jmax in ReadFT'
               STOP
          end if
       End if
       read(30,*, iostat=ierr)LINE1
    end do
    close(30)
  End Subroutine ReadFT1

!$********************************************************************************************************************!!
!$                                              Subroutine ReadFT2
!$  This is the second subroutine for the FT module, it will read from the ACF file. 
!$  It will first read in MaxJ- an integer from the ACF that counts the number of lines in the ACF printed
!$  MaxJ is then used to allocate an array for time and a complex array for ACF
!$  To read in the complex array the real and imaginary parts are read seperately at first then the cmplx function 
!$  is used. This might not be the most efficient method but it does work. 
!$********************************************************************************************************************!!
  Subroutine ReadFT2(Line1,t,ACFarray)! is it worth having Line1 as an argument to be consistent?
    Implicit None
    Character(LEN=10)::LINE1, LINE2, LINE3, LINE4
    integer::ierr
    real(kind=8),dimension(:),allocatable::t
    complex(kind=8),dimension(:),allocatable::ACFarray
    real(kind=8)::RACF,IMGACF
    integer::j
    integer::k
    LOGICAL::Header
    ! ACFHead is a subroutine that needs to be called somewhere and as this is
    ! the point where the ACF decay output will also be written it is worth
    ! having it here. 
    call ACFHead
    ierr=0
    OPEN(UNIT=40, FILE=ACFname,STATUS='OLD', iostat=ierr)
    IF (ierr .ne. 0) then
        Print *,'iostat non zero in READFT2'
        Print *, 'error in opening ACFname file'
    END IF
    allocate(t(jmax))
    allocate(ACFarray(jmax))
    IF (.not. Wdwf=="NO") then
       Print *, 'Window function will be used, Wdwf= ', Wdwf
    END IF
    Header=.TRUE. 
    do while (ierr==0 .AND. Header)
         read(40,*,iostat=ierr)LINE1
         IF (LINE1(1:1)=="#") then
            print *, 'A comment'
         ELSE
         Header=.FALSE. 
         END IF
    END DO
    backspace(40)
    j=1
    read(40,*,iostat=ierr)k,t(j),RACF,IMGACF
    do while (ierr==0 .AND. j .le. jmax)
             IF (ierr .ne. 0) then 
                     Print *, 'error reading ACF'
                     STOP
             END IF
              IF (Wdwf=="NO") then
                   ACFarray(j)=cmplx(RACF,IMGACF)  
              ELSE IF (.not.(Wdwf=="NO")) then 
                  ACFarray(j)=cmplx(RACF,IMGACF)*COS2Fcn(t(j),Tcut)! change this function!!!!!!!!!!!!!!!!
                  call ACFDecayOut(t(j),ACFarray(j))  
              ELSE 
              Print *, 'hasnt gone in the correct if statement'
               STOP    
              END IF
             IF (j .ne. k) then
                Print *,'error in ReadFT2'
                Print *, 'j .ne. k', 'j=',j, 'k=', k
                STOP
             End IF
             j=j+1
             IF (j .le. jmax) then
                 read(40,*,iostat=ierr)k,t(j),RACF,IMGACF
             ENDIF
    End do
    j=j-1
    IF (j .ne. jmax) then
      Print *, 'WARNING end of j isnt eq jmax'
      STOP
    END IF
  End Subroutine ReadFT2

!$********************************************************************************************************************!!
!$                                              Subroutine ACFdecayHead
!$      This subroutine will write the header for the ACFdecayi ouputput file. 
!$      
!$********************************************************************************************************************!!
  Subroutine ACFHead
    Implicit None
    IF (.not. (decayname=="No.decay.out")) then
       OPEN(Unit=43,File=decayname,STATUS='NEW')
       Write(43,*), '#********************ACF Decay OUTPUT********************'
       Write(43,*), '#This output gives the ACF after it has been multiplied by a function (COSSQ)'
       write(43,*), '#**************************************************'
       write(43,*), '# time          Real (ACF*Fcn)          Imag (ACF*Fcn)'
       close(43)
    End IF
  End Subroutine ACFHead

!$********************************************************************************************************************!!
!$                                              Subroutine ACFDecayOut
!$      This subroutine will write one line in the output for the ACFDecayOut 
!$********************************************************************************************************************!!

  Subroutine ACFDecayOut(t, ACFvalue)
  Implicit None
  real(kind=8)::t
  complex(kind=8)::ACFvalue
  IF (.not. (decayname=="No.decay.out")) then
  OPEN(Unit=43,File=decayname,STATUS='OLD', ACCESS='APPEND')
  write(43,*), t, Real(ACFvalue), AIMAG(ACFvalue)
  close(43)
  ELSE
  Return
  End IF
  End Subroutine ACFDecayOut


!$********************************************************************************************************************!!
!$                                              Subroutine FTRtn
!$  This subroutine writes the header for the FToutput file and it writes the FT too.
!$  The do loop calls the FTFCn function and then this is written directly. 
!$********************************************************************************************************************!!

  Subroutine FTRtn(tarray,acfarray)
    Implicit None
    complex(kind=8)::FT
    real(kind=8)::E
    real(kind=8),dimension(:)::tarray
    complex(kind=8),dimension(:)::acfarray
    integer::ierr
    ierr=0
    E=E0
    IF (.not. (FTname =="No.ft.out")) then
       OPEN(UNIT=42,File=FTname,iostat=ierr)
       ! call ierr function
       Write(42,*), '#*******************FT OUTPUT********************'
       Write(42,*), '#', 'E', E, 'Emax', Emax, 'dt', dt
       write(42,*), '#', 'E  ', '   REAL FT(E)   ', '    IMAG FT(E)    '
       do while(E<Emax)
          FT=FTFCN(E,acfarray,tarray)
          ! write out or call subroutine to write out
          !Write(42,*),'E',E,'Real FT(E)',Real(FT),'IMAG FT(E)', AIMAG(FT)
          Write(42,*),E,Real(FT), AIMAG(FT)
          ! the FT value
          E=E+dt 
       End do
       close(42)
    ELSE 
       Return
    EndIF
  End Subroutine FTRtn

!$********************************************************************************************************************!!
!$                                              Function FTFcn
!$  This is a function that will return the Fourier Transformation for a specified E and t value.
!$  The equation is:
!$  FT(E)=Sum_i exp(IEt)ACF(t)DT
!$  This means that the function ACF(t) becomes a function F(E). 
!$********************************************************************************************************************!!
  Function FTFcn(E,acfarray,tarray)
    Implicit None
    complex(kind=8)::FTFcn
    real(kind=8),dimension(:),intent(in)::tarray! needs to be dim (k) ortherwise get seg fault
    complex(kind=8),dimension(:),intent(in)::acfarray
    integer::j
    real(kind=8)::E
    FTFcn=cmplx(0.0d0,0.0d0)
    do j=1,jmax
       FTFcn=FTFcn+exp(Img*E*tarray(j))*acfarray(j)*dt
    end do

  End Function FTFcn

!$********************************************************************************************************************!!
!$                                              Function COS2Fcn
!$  This function will return a value that can be used in the Fourier Transformation Routine.
!$  This will be used to help converge the ACF 
!$********************************************************************************************************************!!
  Function COS2Fcn(t,Tmax)
  Implicit none
  complex(kind=8)::COS2Fcn
  real(kind=8)::Pi=3.141592653589793238462643383279502884197
  real(kind=8)::t
  real(kind=8)::Tmax
  IF (Tmax .ge. t) then 
     COS2Fcn=COS(2.0d0*Pi*t/(4.0d0*Tmax))**2.0d0
  Else IF (Tmax .lt. t) then
     COS2Fcn=0.0d0
  End IF
  End Function

  !Function EFcn(t,Tmax)
  !Implicit none
  !complex(kind=8)::EFcn
  !real(kind=8)::Pi=3.141592653589793238462643383279502884197
  !real(kind=8)::t
  !real(kind=8)::Tmax
  !IF (Tmax .ge. t) then 
  !   EFcn=exp(-t/Tmax)
  !Else If (Tmax .lt. t) then
  !   EFcn=0.0d0
  !End IF
  !End Function


End Module FT1
