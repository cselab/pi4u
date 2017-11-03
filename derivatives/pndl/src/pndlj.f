C  ---------------------------------------------------------------------
	SUBROUTINE PNDLJ ( RSD, X, N, M, IORD, FJ, LD )
C  ---------------------------------------------------------------------
C
C  Description:                              PNDL user interface routine.
C                                            ---------------------------
C    This is a simple interface to the main differentiation
C    routine PNDLJA.
C    Given a multidimensional function that is written a sum of squared
C    terms (residuals):
C                  M
C          F(X) = Sum R_i(X)**2
C                 i=1
C    this routine returns the Jacobian matrix (FJ) by applying a
C    numerical differentiation formula according to the desired order
C    of accuracy (IORD).
C    This routine does not support bounds on the variables. Fatal error
C    messages are printed on the standard output device.
C
C  Input arguments:
C    RSD        A subroutine that returns the residuals (F).
C               Must be declared as:
C                 SUBROUTINE RSD ( X, N, M, F )
C                 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                 DIMENSION X(N), F(M)
C    X          Array containing the function variables.
C    N          The number of variables.
C    M          The number of squared terms.
C    IORD       Order of the numerical derivative.
C               Possible values: 1, 2.
C
C  Output arguments:
C    FJ         The Jacobian matrix. FJ(i,j) = dRi/dXj
C    LD         Leading dimension of matrix FJ.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	EXTERNAL RSD
	DIMENSION X(N), FJ(LD,N)
	EXTERNAL PNDLJA
	COMMON /DUMMY/ NOCx, IERRx
	include 'torcf.h'
C
	PARAMETER ( BIG = 1.0D300 )
C
	DIMENSION XL(N), XU(N), UH(N)
C
C  Use the machine accuracy.
	FEPS = 0.0D0
C
C  Print fatal errors.
	IPRINT = 1
C
C  Set lower/upper bounds to large numbers (=>no bounds).
C  Set stepsizes to zero (=>default steps).
	DO 10,I=1,N
		XL(I) = -BIG
		XU(I) =  BIG
		UH(I) = 0.0D0
10	CONTINUE
C
C  Call the main routine.
	CALL PNDLJA(RSD,X,N,M,XL,XU,UH,FEPS,IORD,
     &              IPRINT,FJ,LD,NOCx,IERRx)
c	call torc_task(PNDLJA, 1, 14,
c     &      1, MPI_INTEGER, CALL_BY_VAD,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      LD*N, MPI_DOUBLE_PRECISION, CALL_BY_RES,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_RES,
c     &      1, MPI_INTEGER, CALL_BY_RES,
c     &      RSD,X,N,M,XL,XU,UH,FEPS,IORD,IPRINT,FJ,LD,NOCx,IERRx)
C
	END
