!%********************************************************************************************************************!!
!%      Author: John Kattirtzi, Date 23/04/09
!%                                                      Module Integration
!%      
!%      This module deals with the Integration part of the program. 
!%      It will consist of two subroutines:
!%      1. Derivs routine- a routine that will give the derivatives of
!%              the basis set      
!%      2. RK4ccs- a 4th order Runge-kutta subroutine that will use the
!%      type made in the bset module
!%      The rk4 routine in numerical recipes uses assert_eq function to
!%      check that the arrays are the same size- can't copy this easily
!%      so might use an IF statement instead
!%      
!%      Links with other modules:
!%      Need to link with BSET to understand type and overload the
!%      operators
!%      Need to link with program/module that calls the integration
!%      Need to link with GenHam to get the Hord and DHord that is needed 
!%        
!%      Global parameters for this module:
!%      Img=0.0d0,1.0d0
!%********************************************************************************************************************!!
  Module Integration
    !use BSET, only: csbasisfn,ndim,nbf
    use BSET2!, only: operator(+) , operator(*)  
    use GenHam
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Global parameters:
    complex(kind=8),PARAMETER::Img=cmplx(0.0d0,1.0d0)
    contains

!$********************************************************************************************************************!!
!$                                                      Subroutine ReadTime
!$      This subroutine will read in the time parameters:
!$      t0- initial time
!$      tmax- maximum time you want the trajectory to run for
!$      dt- time increment
!$********************************************************************************************************************!!
  
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


!$********************************************************************************************************************!!
!$                                                      Subroutine Derivs
!$      This subroutine will take in the basis set and give out the derivative
!$      of the z,s and ds. 
!$      The arguments are:
!$      t- time
!$,     bs- basis set(input)
!$      dbsdt- derivative of basis set(output)
!$       
!$      Equation for zdot:
!$      zdot=-i/hbar*dhord(zic,zi)/dzic
!$      Equation for Sdot=L
!$      Sdot=L=ihbar/2.0sum_ndim((zc*zdot-zdotc*z)-Hord(zc,z))
!$
!$********************************************************************************************************************!!

    Subroutine Derivs(t,bs,dbsdt)
      Implicit none
      type(csbasisfn),dimension(:),intent(in)::bs
      type(csbasisfn)::bf1cz,bf2, bf1
      type(csbasisfn),dimension(:),intent(inout)::dbsdt! need to have this as in for s
      type(csbasisfn),dimension(:),allocatable::dbsdtc
      complex(kind=8),dimension(:),allocatable::ddot ! by this i mean nbf, not sure if will work
      real(kind=8)::t
      integer::i,j
      ! The IF statement below will check that bs and dbsdt are same size
      !If (.not. size(bs) ==size(dbsdt)) then
      !   Print *, 'input and output arrays are different sizes in derivs'
      !   STOP
      !End If
      !!!!!!!!!!!!!!!!!!!!!Z derivs below:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call INBASISFN(bf1cz)
      call INBASISFN(bf2)
      call ALCBASIS(dbsdtc)
      !call INBASISFN(dbsdt)! will want to move this before derivs is called 
      !call INBASISFN(dbsdtc)
      do i=1,nbf
         bf1=bs(i)                                                ! assign the zi and zic bits
         bf1cz=bf1 ! might not need this check later
         bf1cz%z(1:ndim)=dconjg(bf1%z(1:ndim))
         do j=1,ndim
            dbsdt(i)%z(j)=-(Img/hbar)*DH_ord_Dz1c(bf1cz,bf1,j)     ! DHord takes in a j for a degree of freedom
!            Print *, 'called deriv z'
         end do
      end do
      !!!!!!!!!!!!!!!!!!!!!S derivs below:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Print *, 'ACT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     ! Print *,' in act', Img, hbar, 'time', t
      do i=1,nbf
         bf1=bs(i)
         bf1cz%z(1:ndim)=dconjg(bf1%z(1:ndim))
         dbsdt(i)%s=Img*0.5d0*sum(bf1cz%z(1:ndim)*dbsdt(i)%z(1:ndim)-dconjg(dbsdt(i)%z(1:ndim))*bf1%z(1:ndim))
         dbsdt(i)%s=dbsdt(i)%s - H_ord(bf1cz,bf1)
      end do
      !Print *, 'deriv action', dbsdt(1:nbf)%s
!      !!!!!!!!!!!!!!!!!!!!!D derivs below:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(ddot(nbf))
      call Amp2d(bs,ddot)
      dbsdt(1:nbf)%d=ddot(1:nbf)
      !dbsdt(1:nbf)%s=0.0d0
      !do i=1,nbf
      !  do j=1,ndim
      !    dbsdt(i)%z(j)=0.0d0
      !  end do
      !end do
      !Print *, 'dbsdt', dbsdt(1:nbf)%d
      deallocate(ddot)
    End Subroutine Derivs

!$********************************************************************************************************************!!
!$                                                      Subroutine Amp2d
!$      This subroutine will do what Ampsd was supposed to do but failed
!$      I have decided to keep this routine out of the derivs routine to keep the derivs routine more simple
!$      Arguments:
!$      bs- input-  bs(i)%d will be used
!$      ddot-output- which will go into dbsdt(i)%d
!$      external:
!$      need to call overlap matrix and  Hord (Delta2Hord) functions
!$      need to call LAPACK
!$      Equation solved in element form:
!$      Sum_i D_i^dot *e^i/hbar*S_i*Omega_ji=i/hbar*Sum_i d_i e^(i/hbar*S_i)*Omega_ji*D2Hordji
!$      In LAPACK I will solve
!$      ==> Omega*x=b
!$      where Omega is the overlap matrix
!$      x is the vector ddot*exp(i/hbar*Si)                     - ddot *exp is a scalar multiplication
!$      b is the vector:
!$      b=OmegaD2*diexp(i/hbar*Si)                              -Matrix vector multiplication
!$      where OmegaD2 is a matrix where the elements are constructed from a scalar multiplication of the ij elements 
!$      diexp(i/hbar*si) is a vector where the element di is multiplied by exp(i/hbar*si) - also called vector a!
!$      After solved with LAPACK need to divide x by exp(i/hbar*Si) to get ddot
!$********************************************************************************************************************!!

    Subroutine Amp2d(bs,ddot)
      Implicit none
      type(csbasisfn),dimension(:),intent(in)::bs
      type(csbasisfn)::bf1,bf2
      complex(kind=8),dimension(:),intent(out)::ddot
      complex(kind=8),dimension(:),allocatable::a
      complex(kind=8),dimension(:),allocatable::b
      complex(kind=8),dimension(:,:),allocatable::OMEGA
      complex(kind=8),dimension(:,:),allocatable::D2
      complex(kind=8),dimension(:,:),allocatable::OMD2
      complex(kind=8),dimension(nbf,nbf)::E ! new phase matrix
      integer::i,j
      ! variables for LAPACK
      integer, dimension(nbf) ::IPIV
      integer:: LDA,LDB, N, NRHS, INFO
      complex(kind=8),dimension(nbf,nbf)::Omegain
      complex(kind=8),dimension(nbf)::bout
      !INTERFACE
       include "./LAPACK/ZGESV.f90"
      !END INTERFACE
      allocate(OMEGA(nbf,nbf))
      allocate(OMD2(nbf,nbf))
      allocate(D2(nbf,nbf))
      allocate(a(nbf))
      allocate(b(nbf))
      call OVRLMAT(OMEGA,bs)
      CALL INBASISFN(bf1)
      CALL INBASISFN(bf2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Right Hand Side!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! sort out right hand side matrix
      do i=1,nbf
         bf2=bs(i)
         do j=1,nbf
            bf1=bs(j)
            IF ( i .eq. j) then
               D2(j,i)=cmplx(0.0d0, 0.0d0)
            ELSE IF (i .ne. j) then
                D2(j,i)=Delta2Hord(bf1,bf2)
            END IF
         end do
      end do
      !Print *, 'D2', D2(1:nbf,1:nbf)
      do i=1,nbf 
         do j=1,nbf
                IF (j .ne. i) then
                        E(j,i)=0.0d0
                Else if (j .eq. i) then
                        E(i,j)=exp(Img/hbar*(bs(i)%s))!-bs(j)%s))  - find it confusing as it matches what i have in notes but differs to paper
                        !Print *, 'e', E(i,j)
                End IF
         end do
      end do
      do i=1,nbf
         do j=1,nbf
            OMD2(j,i)=OMEGA(j,i)*D2(j,i)
         end do
      end do
      ! sort out right hand side vector a
      do j=1,nbf
        a(j)=bs(j)%d*E(j,j)                              
      end do
      ! sort out right hand side vector b that goes into lapack
     ! do i =1,nbf
     ! do j=1,nbf
     ! Print *, 'Omega', OMEGA(j,i) !Real(i)*2.0d0+Real(j)
     ! end do
     ! end do
      b(:)=MATMUL(OMD2,a)
      b(:)=-1.0d0*(Img/hbar)*b(:)
      !Print *, 'b before LAPACK', b
      ! matrix equation now is : Omega * x=b , where x is ddot*exp(i/hbar*si) and b is vector above- solve with LAPACK
      ! variables for LAPACK****************************************************************
      omegain(:,:)=OMEGA(:,:)
      LDA=nbf
      LDB=nbf
      NRHS=1
      N=nbf
      INFO=0
      IPIV=0
      !Print *, 'LDA', LDA, 'LDB', LDB,'NRHS', NRHS, 'N', N, 'INFO', INFO, 'IPIV', IPIV
      !***********************************************************************************
      CALL ZGESV(N, NRHS, omegain , LDA, IPIV, b , LDB, INFO )
      bout(:)=b(:)
      If (INFO .ne. 0) then
         Print *, 'WARNING INFO =', INFO
      End If
      ! ddot=b/exp(i/hbar*si)
      do j=1,nbf
        ddot(j)=bout(j)/E(j,j)
      end do
      b(:)=MATMUL(Omega,bout)! omegain gets changed by LAPACK so need to use OMEGA
      !Print *, 'b after LAPACK', b
      !STOP
      deallocate(D2)
      deallocate(OMD2)
      deallocate(OMEGA)
      deallocate(b)
      deallocate(a)
      deallocate(bf1%z)
      deallocate(bf2%z)
    End Subroutine Amp2d


!$********************************************************************************************************************!!
!$                                                      Subroutine RK4ccs
!$      This subroutine will be heavily based on the one in the numberical
!$      recipes book, p1297 (second volume) but it will work for the type
!$      defined in the basis set module.
!$      Double precision will be used instead of SP and DP to be consistent
!$      No interface will be needed as the derivs routine will also be in this
!$      module. 
!$      The arguments are:
!$      y - this is bs
!$      dydx- this is bsdot (the deriv of bs)
!$      x- time t- I will actually change this to t so not to be confusing
!$      h- time step, I will call this dt 
!$      yout- the integrated output of the bs
!$      The RK4 routine has derivs as an argument but I dont think it is needed
!$      if it is in the same module.
!$********************************************************************************************************************!!

    Subroutine RK4ccs(y,dydx,x,h,yout)
      Implicit none
      type(csbasisfn),dimension(:),intent(in)::y,dydx
      type(csbasisfn),dimension(:),intent(out)::yout
      real(kind=8),intent(in)::x,h
      real(kind=8)::h6,hh,xh
      type(csbasisfn),dimension(:),allocatable::dym,dyt,yt
      !if(.not. (size(y) == size (dydx) .OR. size(y)== size(yout) .or. size(yout)==(size(dydx)) ))then
      !   Print *, 'error in size of bs arrays in RK4ccs'
      !end if
      !Print *, ' starting derivs'
      CALL ALCBASIS(dym)
      CALL ALCBASIS(dyt)
      CALL ALCBASIS(yt)
      !Print *, 'allocated memory in derivs'
      hh=h*0.5d0
      h6=h/6.0d0
      xh=x+hh
      yt(:)=hh*dydx(:)
      yt(:)=yt(:)+y(:)
      !yt=y+tt*dydx ! first step
      !Print *, 'in rk4 calling derivs'
      call derivs(xh,yt,dyt)!second
      ! equation was: yt=y+tt*dyt but had to split it up
      yt(:)=hh*dyt(:)!+y
      yt(:)=yt(:)+y(:)
      !Print *, 'in rk4 calling derivs'
      call derivs(xh,yt,dym)! third step
      !equation was: yt= y+tt*dym
      yt(:)=h*dym(:) ! changed this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      yt(:)=yt(:)+y(:)
      dym(:)=dyt(:)+dym(:)
      !Print *, 'in rk4 calling derivs'
      call derivs(x+h,yt,dyt)! fourth step
      ! equation was: yout=y+t6*(dydx+dyt+2.0d0*dym)! accumulate the proper weights
      dym(:)=2.0d0*dym(:)
      yout(:)=dydx(:)+dyt(:)
      yout=yout(:)+dym(:)
      yout=h6*yout(:)
      yout=yout(:)+y(:)
      DEALLOCATE(dym)
      DEALLOCATE(dyt)
      DEALLOCATE(yt)
      !Print *, 'finished derivs'
    End Subroutine RK4ccs
      

  End Module Integration

