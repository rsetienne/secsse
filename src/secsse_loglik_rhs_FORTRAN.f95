! Helper function: 
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE secsse_fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
     
        DO I = 1, DIMP
          II = II + 1
          vec(I) = parms(II)
        ENDDO
        
      END SUBROUTINE secsse_fill1d

!==========================================================================
! module with declarations
!==========================================================================

      MODULE secsse_dimmod

      ! length of the vector -  decided in R-code
      INTEGER  :: N
      
      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)                                
      
      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE secsse_dimmod

!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE secsse_initmod (steadyparms)
      USE secsse_dimmod 

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 1  ! constant-length parameters
      
      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block 

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector       
      N = INT(parms(1) + 1e-6)  

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)  
      ALLOCATE(P(N * N))

      initialised = .FALSE.
       
      END SUBROUTINE secsse_initmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE secsse_runmod (neq, t, Conc, dConc, yout, ip)
      USE secsse_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii

      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*), FF1, FF2


! parameters - named here
      DOUBLE PRECISION rn
      COMMON /XCBPar/rn


! local variables
      CHARACTER(len=80) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL secsse_fill1d(P, N * N / 4 + N, yout, ii)   ! ii is updated in Fill1D
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

!  R code
!  dE<- mus-(lambdas + mus + Q %*% (rep(1,d))) * Es + lambdas * Es * Es + ( Q %*% Es )

      ! Es

      DO I = 1, N/2
        FF1 = P(N/2 + I) - (P(I) + P(N/2 + I)) * Conc(I)
        FF2 = P(I) * Conc(I) * Conc(I)
        dConc(I) = FF1 + FF2
        DO II = 1, N/2
           FF1 = Conc(II) - Conc(I)
           FF2 = P(N + (I - 1) * N/2 + II) * FF1
           dConc(I) = dConc(I) + FF2
        ENDDO
      ENDDO

!  R code
!  dD<- -(lambdas + mus + Q %*% (rep(1,d))) * Ds + 2 * lambdas * Es * Ds + ( Q %*% Ds )

      ! Ds
      DO I = 1, N/2
        FF1 = (-P(I) - P(N/2 + I) + 2 * P(I) * Conc(I)) * Conc(N/2 + I)
        dConc(N/2 + I) = FF1
        DO II = 1, N/2
           FF1 = Conc(N/2 + II) - Conc(N/2 + I)
           FF2 = P(N + (I - 1) * N/2 + II) * FF1
           dConc(N/2 + I) = dConc(N/2 + I) + FF2
        ENDDO
      ENDDO

      END SUBROUTINE secsse_runmod
      
