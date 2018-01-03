C  ---------------------------------------------------------------------
      SUBROUTINE PNDLHF ( F, X, N, IORD, HES, LD )
C  ---------------------------------------------------------------------
C
C  Description:                              PNDL user interface routine.
C                                            ---------------------------
C    This is a simple interface to the main differentiation
C    routine PNDLHFA.
C    Given a multidimensional function (F), this routine returns the
C    Hessian matrix (HES) by applying a numerical differentiation
C    formula according to the desired order of accuracy (IORD).
C    This routine does not support bounds on the variables. Fatal error
C    messages are printed on the standard output device.
C
C  Input arguments:
C    F          The function to be differentiated. Must be declared as:
C                 FUNCTION F ( X, N )
C                 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                 DIMENSION X(N)
C    X          Array containing the function variables.
C    N          The number of variables for this function.
C    IORD       Order of the numerical derivative.
C               Possible values: 1, 2.
C  Output arguments:
C    HES        Array containing the resulting Hessian matrix.
C               Note that only the lower triangular part (plus the
C               diagonal elements) is returned.
C    LD         Leading dimension of matrix HES.
C
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      DIMENSION X(N), HES(LD,N)
      EXTERNAL PNDLHFA
      COMMON /DUMMY/ NOCx, IERRx
      include 'torcf.h'
C
      PARAMETER ( BIG = 1.0D300 )
C
      DIMENSION XL(N), XU(N), UH(N)
C
C     Use the machine accuracy.
      FEPS = 0.0D0
C
C     Print fatal errors.
      IPRINT = 1
C
C     Set lower/upper bounds to large numbers (=>no bounds).
C     Set stepsizes to zero (=>default steps).
      DO 10,I=1,N
         XL(I) = -BIG
         XU(I) =  BIG
         UH(I) = 0.0D0
10    CONTINUE
C
C     Call the main routine.
      CALL PNDLHFA(F,X,N,XL,XU,UH,FEPS,IORD,IPRINT,HES,LD,NOC,IERR)
      END SUBROUTINE PNDLHF
