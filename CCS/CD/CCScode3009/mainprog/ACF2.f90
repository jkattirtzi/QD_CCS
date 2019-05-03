!%********************************************************************************************************************!!
!%                                              Module: ACF
!%      Author: John Kattirtzi, Date: 03/06/09
!%      
!%      This module will consist of the subroutines related to the Autocorrelation
!%      Function.
!%      A subroutine to read if ACF will be calculated
!%      A subroutine to read if FT will be calculated
!%      A subroutine to read in E0, Emax
!%      Need to write output for ACF and FT(E) either here or in output
!%      Output needs to show Frank-Condon Spectrum
!%********************************************************************************************************************!!
Module ACFMOD

USE BSET2
complex(kind=8)::Img=cmplx(0.0d0,1.0d0)
!Character(LEN=10)::ACFOP,FTOP
Character(LEN=20)::ACFname="NoACFoutin.dat"
!Character(LEN=20)::FTname="NoFTout.dat"
contains

!$********************************************************************************************************************!!
!$                                              Subroutine ZPS
!$      This function would be better suited to the basisset module.
!$      However, it is only needed to calculate the ACF and so I have put it here for now
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

  Subroutine Readmpmq(Line1,mup,muq)
    Implicit None
    real(kind=8),intent(inout),dimension(:),allocatable::mup,muq
    integer::i,j
    Character(Len=10)::LINE1
    integer::ierr
    ierr=0
    allocate(mup(ndim))
    allocate(muq(ndim))
    OPEN(UNIT=27, FILE='input.dat',STATUS='OLD', iostat=ierr)
    read(27,*,iostat=ierr)LINE1 
    do while (ierr==0)
       if(LINE1 == "meanp") then
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
       end if
       read(27,*,iostat=ierr)LINE1
    end do
  END Subroutine Readmpmq

!$********************************************************************************************************************!!
!$                                              Subroutine ReadACF
!$    This subroutine will read the name of the ACFoutput file
!$********************************************************************************************************************!!

  Subroutine ReadACF(LINE1,mup,muq)
    Character(Len=10)::LINE1
    real(kind=8),intent(inout),dimension(:),allocatable::mup,muq
    integer::ierr
    ierr=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
    print *, 'in read acf before read', ierr
    rewind(30)
    read(30,*,iostat=ierr)LINE1 
    IF (ierr .ne. 0) then
        Print *, 'error reading line1 in readacf', LINE1
    END IF
    do while (ierr==0)
       if(LINE1 == "ACF") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,ACFname
          call Readmpmq(LINE1,mup,muq)
       End if
       read(30,*, iostat=ierr)LINE1
    end do
    Print *, 'ACF option', ACFname
    close(30)
  End Subroutine ReadACF

!$********************************************************************************************************************!!
!$                                              Function ACF
!$  This function will return the value of the autocorrelation function
!$  The equation is:
!$  ACF=Sum_i d_i*exp(I/hbar*s_i)*(PROD_ndim<z_0|z_i>)
!$********************************************************************************************************************!!

  
  Function ACF(bs,muq, mup)
    Implicit none
    type(csbasisfn),dimension(:)::bs
    complex(kind=8),dimension(:),allocatable::zp
    real(kind=8),dimension(:),intent(in)::mup
    real(kind=8),dimension(:),intent(in)::muq
    complex(kind=8)::ACF
    complex(kind=8)::ovrlpbit
    integer::i,j
    allocate(zp(ndim))
    zp(:)=ZPS(muq,mup)
  !  Print *, 'zp', zp, 'd', bs%d, 's', bs%s
    ACF=cmplx(0.0d0,0.0d0)
    !ovrlpbit=cmplx(1.0d0,0.0d0)
    do i=1,nbf
          ovrlpbit=Product(exp(dconjg(zp(1:ndim))*bs(i)%z(1:ndim)-0.5d0*dconjg(zp(1:ndim))*zp(1:ndim)& ! when doing multidim might need to reset ovrlpbit before this step but i think this should be fine
               & -0.5d0*dconjg(bs(i)%z(1:ndim))*bs(i)%z(1:ndim)))
       ACF=ACF+bs(i)%d*exp(Img/hbar*bs(i)%s)*ovrlpbit
    end do
    deallocate(zp)
  End  Function ACF
  

  Function ACF2H(bs,muq,mup)
    Implicit None
    type(csbasisfn),dimension(:)::bs
    type(csbasisfn)::bfkstar
    real(kind=8),dimension(:),intent(in)::mup
    real(kind=8),dimension(:),intent(in)::muq
    complex(kind=8)::ACF2H
    integer::j,k
    allocate(bfkstar%z(ndim))
    do j=1,nbf
       do k=1,nbf
          bfkstar%d=dconjg(bs(k)%d)
          bfkstar%z=dconjg(bs(k)%z)
          bfkstar%s=(bs(k)%s)
          ACF2H=ACF2H+bs(k)%d*bs(j)%d*exp(Img*(bs(k)%s+bs(j)%s))*OVRLPELMT(bfkstar,bs(j))
       end do
    end do
    deallocate(bfkstar%z)
  End Function ACF2H

!$********************************************************************************************************************!!
!$                                              Subroutine ACFHead
!$  This subroutine will write the header to the ACF out file
!$********************************************************************************************************************!!


  Subroutine ACFHead
  Implicit None
  IF (.not. (ACFname=="NoACFoutin.dat")) then
     OPEN(Unit=35,File=ACFname,STATUS='NEW')
     Write(35,*), '#********************ACF OUTPUT********************'
     Write(35,*), '#This output gives the ACF as a function of t'
     write(35,*), '#**************************************************'
     write(35,*),'#'
     write(35,*), '# Assuming the wave function to be Gaussian (ground state CS)'
     write(35,*), '#','j', 't', 'Real ACF value', 'Imag ACF value'
     close(35)
  End IF
  End Subroutine ACFHead

!$********************************************************************************************************************!!
!$                                              Subroutine ACFtail
!$  This subroutine will write the end of the ACF out file
!$********************************************************************************************************************!!

  Subroutine ACFTail(j)
  Implicit None
  integer::j
  IF  (.not. (ACFname=="NoACFoutin.dat")) then
     !Open(unit=35,File=ACFname, STATUS='OLD', ACCESS='APPEND')
     Open(unit=35,File=ACFname, STATUS='OLD', POSITION='APPEND')
     Write(35,*),'#','MAXJ', j
     !Write(35,'(MAXJ)',REC=5) j
     !close(35)
     Write(35,*)'#END of ACF OUTPUT'
     Write(35,*)'#***************************************************'
     close(35)
  End IF
  End Subroutine ACFTail
!$********************************************************************************************************************!!
!$                                              Subroutine ACFrtn
!$  This subroutine calls the ACF function and writes out the ACF value for a specific time
!$  providing that the ACFfile is required. 
!$********************************************************************************************************************!!
  Subroutine ACFrtn(bs,mup,muq,ACFval,j,t)
    Implicit None
    type(csbasisfn),dimension(:),intent(in)::bs
    real(kind=8),dimension(:),intent(in)::mup
    real(kind=8),dimension(:),intent(in)::muq
    real(kind=8),intent(in)::t
    complex(kind=8),intent(out)::ACFval
    integer::j
    IF  (.not. (ACFname=="NoACFoutin.dat")) then
        OPEN(UNIT=35, FILE=ACFname,STATUS='OLD',ACCESS='APPEND')
        !ACFval=ACF(bs,muq,mup)
        ACFval=ACF2H(bs,muq,mup)
        Write(35,*), j,t,REAL(ACFval), AIMAG(ACFval)
        close(35)
    END IF
  End Subroutine ACFrtn

END MODULE ACFMOD

!  Subroutine FTRtn(acfarray,dt,k,FT)
!    complex(kind=8),intent(out)::FT







!       else if(LINE1=="FT") then
!          backspace(30)
!          read(30,*,iostat=ierr)LINE1,FTOP
!       else if(LINE1=="E0")then
!          backspace(30)
!          read(30,*,iostat=ierr)LINE1,E0
!       else if(LINE1=="Emax")then
!          backspace(30)
!          read(30,*,iostat=ierr)LINE1,Emax
!    IF (ACFOP =="ON" .AND. FTOP=="OFF") then
!       Print *, 'You want me to calculated FT without ACF?'
!       Print *, 'I think you are insane, change this !!!'
!       STOP
!    End IF
