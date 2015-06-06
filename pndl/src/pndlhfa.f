C  ---------------------------------------------------------------------
	SUBROUTINE PNDLHFA ( F, X, N, XL, XU, UH, FEPS, IORD, IPRINT, HES, 
     &                    LD, NOC, IERR ) 
C  ---------------------------------------------------------------------
C
C  Description:                              PNDL user interface routine.
C                                            ---------------------------
C    Given a multidimensional function (F), this routine returns the 
C    Hessian matrix (HES) by applying a numerical differentiation
C    formula according to the desired order of accuracy (IORD).
C
C  Input arguments:
C    F          The function to be differentiated. Must be declared as:
C                 FUNCTION F ( X, N )
C                 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                 DIMENSION X(N)
C    X          Array containing the function variables.
C    N          The number of variables for this function.
C    XL         Array containing lower bounds on the variables.
C    XU         Array containing upper bounds on the variables.
C    UH         Array containing user-speficied steps.
C    FEPS       An estimation of the relative accuracy with which 
C               the function is evaluated.
C               If FEPS=0 the machine accurracy is computed and used
C               instead.
C    IORD       Order of the numerical derivative.
C               Possible values: 1, 2.
C    IPRINT     This option controls the amount of printout from the
C               routine. Note that all output appears on the standard
C               output device. Possible values are:
C                 IPRINT=0 -> No printout at all.
C                 IPRINT=1 -> Fatal error messages are printed.
C                 IPRINT=2 -> Warning messages are printed.
C                 IPRINT=3 -> Detailed information is printed (the
C                             formula that was used, differentiation
C                             steps etc).
C    LD         Leading dimension of matrix HES.
C
C  Output arguments:
C    HES        Array containing the resulting Hessian matrix.
C               Note that only the lower triangular part (plus the
C               diagonal elements) is returned.
C    NOC        Number of calls to the function F.
C    IERR       Error indicator. Possible values are:
C                 IERR=0 -> No errors at all.
C                 IERR=1 -> The supplied IORD is incorrect.
C                 IERR=2 -> The supplied N is less than 1.
C                 IERR=3 -> Some of the supplied upper bounds (XU) 
C                           are less than the lower bounds (XL).
C                 IERR=4 -> Some of the supplied values in X do
C                           not lie inside the bounds.
C                 IERR=5 -> The supplied value of FEPS is incorrect
C                           (less than 0 or greater than 1).
C                 IERR=6 -> The supplied IPRINT is incorrect.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	EXTERNAL F
	DIMENSION X(N), XL(N), XU(N), UH(N), HES(LD,N)
	EXTERNAL PNDL2FF, PNDL2FQ
	include 'torcf.h'
C
	DIMENSION TMP3(N), TMP4(N), ITMP(N)
C
C  Informative message (routine starts executing).
	CALL PNDLMSG(IPRINT,'PNDLHFA',0)
	IF (IPRINT.EQ.3) WRITE (*,10) N, FEPS, IORD
C
C  Check validity of arguments.
	CALL PNDLARG(X,N,XL,XU,FEPS,IPRINT,IERR)
	IF (IERR.NE.0) GOTO 100
C
C  Print bounds.
	CALL PNDLBND(X,XL,XU,N,IPRINT)
C
C  Call a routine to do the actual computation according to the
C  desired order of accuracy.
	IF (IORD.EQ.1) THEN
		CALL PNDL2FF(F,X,N,XL,XU,UH,FEPS,IPRINT,HES,LD,NOC)
c	call torc_task(PNDL2FF, 1, 11,
c     &      1, MPI_INTEGER, CALL_BY_VAD,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      LD*N, MPI_DOUBLE_PRECISION, CALL_BY_RES,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_RES,
c     &      F,X,N,XL,XU,UH,FEPS,IPRINT,HES,LD,NOC)
	ELSE IF (IORD.EQ.2) THEN
		CALL PNDL2FQ(F,X,N,XL,XU,UH,FEPS,IPRINT,HES,LD,NOC)
c	call torc_task(PNDL2FQ, 1, 11,
c     &      1, MPI_INTEGER, CALL_BY_VAD,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      LD*N, MPI_DOUBLE_PRECISION, CALL_BY_RES,
c     &      1, MPI_INTEGER, CALL_BY_VAL,
c     &      1, MPI_INTEGER, CALL_BY_RES,
c     &      F,X,N,XL,XU,UH,FEPS,IPRINT,HES,LD,NOC)
	ELSE
		WRITE (*,20) IORD
		IERR = 1
	END IF
100	CALL PNDLMSG(IPRINT,'PNDLHFA',1)
C
C  Informative message (routine starts executing).
10	FORMAT (' PNDL: ',3X,'N = ',I6,3X,'FEPS = ',1PG21.14,3X,'IORD = ',
     &        I2)
20	FORMAT (/' PNDL: Error: Incorrect order (IORD) ',I6/)
	END
