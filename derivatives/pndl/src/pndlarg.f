C  ---------------------------------------------------------------------
	SUBROUTINE PNDLARG ( X, N, XL, XU, FEPS, IPRINT, IERR )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Checks the input arguments, making sure that all of them
C    have reasonable values. The output variable IERR is set 
C    accordingly.
C
C  Input arguments:
C    X          Array containing the function variables.
C    N          The number of variables.
C    XL         Array containing lower bounds on the variables.
C    XU         Array containing upper bounds on the variables.
C    IPRINT     Printout control option. Possible values are:
C                 IPRINT=0    -> No printout at all.
C                 IPRINT.NE.0 -> Print error messages.
C
C  Output arguments:
C    IERR       Error indicator. Possible values are:
C                 IERR=0 -> No errors at all.
C                 IERR=2 -> The supplied N is less than 1.
C                 IERR=3 -> Some of the supplied upper bounds (XU)
C                           are less than the lower bounds (XL).
C                 IERR=4 -> Some of the supplied values in X do
C                           not lie inside the bounds.
C                 IERR=5 -> The supplied value of FEPS is incorrect
C                           (less than 0 or greater than 1).
C                 IERR=6 -> The supplied IPRINT is incorrect.
C
C  Input / Output arguments:
C    FEPS       On input an estimation of the relative accuracy with 
C               which the function is evaluated.
C               If FEPS=0 the machine accurracy is computed and returned
C               in FEPS.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), XL(N), XU(N)
C
	LOGICAL PR
C
	IERR = 0
	PR = IPRINT.GT.0
C
C  Check the printout control option.
	IF (IPRINT.LT.0 .OR. IPRINT.GT.3) THEN
		WRITE (*,100) IPRINT
		IERR = 6
	END IF
C
C  Make sure that the number of variables is at least 1.
	IF (N.LT.1) THEN
		IF (PR) WRITE (*,200) N
		IERR = 2
		RETURN
	END IF
C
C  Make sure the lower bound is less than the upper bound.
	DO 10,I=1,N
		IF (XU(I).LT.XL(I)) THEN
			IF (PR) WRITE (*,300) I
			IERR = 3
		END IF
10	CONTINUE
	IF (IERR.NE.0) RETURN
C
C  Check whether the variables lie inside the supplied bounds.
	DO 20,I=1,N
		XI = X(I)
		IF (XI.GE.XU(I) .OR. XI.LE.XL(I)) THEN
			IF (PR) WRITE (*,400) I
			IERR = 4
		END IF
20	CONTINUE
	IF (IERR.NE.0) RETURN
C
C  Make sure FEPS is valid.
	IF (FEPS.LT.0.0D0 .OR. FEPS.GT.1.0D0) THEN
		IF (PR) WRITE (*,500) FEPS
		IERR = 5
	END IF
C
C  When FEPS is set to zero we must use the machine accuracy instead.
	IF (FEPS.EQ.0.0D0) THEN
C
C  Compute the machine accuracy.
		CALL PNDLACC(COMACC)
		FEPS = COMACC
		IF (IPRINT.EQ.3) WRITE (*,600) COMACC
	END IF
C
100	FORMAT (/' PNDL: Error: Incorrect print level (IPRINT) ',I6/)
200	FORMAT (/' PNDL: Error: Incorrect number of variables (N) ',I6/)
300	FORMAT (/' PNDL: Error: The upper bound (XU) is less than the ',
     &        'lower bound (XL) for variable ',I6/)
400	FORMAT (/' PNDL: Error: Variable ',I6,' is out of the specified ',
     &        'bounds'/)
500	FORMAT (/' PNDL: Error: Incorrect function accuracy (FEPS) ',
     &        1PG21.14/)
600	FORMAT (/' PNDL: Setting FEPS to computer accuracy: ',1PG21.14)
	END
