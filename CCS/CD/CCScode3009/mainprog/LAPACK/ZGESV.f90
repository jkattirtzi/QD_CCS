Interface 
        SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        INTEGER :: N
        INTEGER :: NRHS
        INTEGER :: IPIV(*)
        INTEGER :: LDA, LDB
        COMPLEX*16:: A(LDA,*), B(LDB,*)
        INTEGER :: INFO
        END SUBROUTINE
END INTERFACE 
        
!*
!*  -- LAPACK driver routine (version 3.1) --
!*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!*     November 2006
!*
!*     .. Scalar Arguments ..
!      INTEGER            INFO, LDA, LDB, N, NRHS
!*     ..
!*     .. Array Arguments ..
!      INTEGER            IPIV( * )
!      COMPLEX*16         A( LDA, * ), B( LDB, * )

!END SUBROUTINE ZGESV

