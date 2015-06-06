	PROGRAM TEST1
C  ---------------------------------------------------------------------
C  PNDL test program.
C  -----------------
C  Computes the first derivative of several one-dimensional
C  functions, with all possible orders of accuracy.
C  The computation takes place at several different points, taking
C  into account all possible situations of point-bound proximity.
C  Relative errors against the analyticaly known derivatives are
C  evaluated and printed. Since the amount of output is quite large,
C  an overall failure indication is printed, based on the comparison
C  of each relative error with the FAIL threshold (defined below).
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INCLUDE 'torcf.h'
	EXTERNAL F
C
C  Number of test points.
	PARAMETER ( NPT = 4 )
C
C  Number of test functions.
	PARAMETER ( NFUN = 1 ) 
C
C  Derivative computations with relative errors above 
C  this threshold are considered "FAILed".
	PARAMETER ( FAIL = 1.0D-6 )
C
C  BND is the distance of a bound to the evaluation point
C  (in order to enforce point-bound proximity).
	PARAMETER ( BND = 1.0D-9 )
C
C  Printout control option.
	PARAMETER ( IPRINT = 0 )
C
C  A large number, to be used as lower/upper bound.
	PARAMETER ( BIG = 1.0D300 )
C
	DIMENSION X(NPT)
	CHARACTER*4 STATA, STATB
	CHARACTER*40 FNAME
	COMMON / SELECT / ISEL
C
C  The points at which the derivatives will be tested.
	DATA X / -4.0D0, -2.0D0, 0.0D0, 2.0D0 /
C
c	CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, IPROV, IERR)
	CALL PNDL_INIT()
	CALL TORC_INIT()

C	CALL MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
	IRANK = 0

C  The relative error in the evaluation of our test functions is
C  in the order of the machine accuracy.
	FEPS = 0.0D0
C  No stepsize specified 
	UH = 0.0D0
C
C  NFAIL is incremented each time the relative error is above 
C  the "FAIL" threshold,
	NFAIL = 0
C
C  Loop over all point-bound proximity situations.
	DO 40,M=-1,1
		IF (IRANK.EQ.0) THEN
		WRITE (*,*)
		CALL WLINE('=')
		IF (M.EQ.-1) THEN
			CALL CPRINT('X IS NEAR THE LOWER BOUND')
		ELSE IF (M.EQ.0) THEN
			CALL CPRINT('NO BOUNDS')
		ELSE IF (M.EQ.1) THEN
			CALL CPRINT('X IS NEAR THE UPPER BOUND')
		END IF
		CALL WLINE('=')
		ENDIF

C
C  Loop over all test functions.
		DO 30,IFUN=1,NFUN
			ISEL = IFUN
			IF (IRANK.EQ.0) THEN
			CALL NAMEIT(FNAME)
			ENDIF
C
C  Loop over all test points.
			DO 20,K=1,NPT,2
				IF (IRANK.EQ.0) THEN
				X1 = X(K)
				X2 = X(K+1)
				CALL WLINE('-')
				CALL CPRINT('F(X) = '//FNAME)
				WRITE (*,210) X1, X2
				CALL WLINE('-')
				WRITE (*,200)
				CALL WLINE('-')
				CALL GRAD(X1,1,GA)
				CALL GRAD(X2,1,GB)
				ENDIF
C				CALL PNDL_BARRIER()
C
C  Loop over all available orders of accuracy.
				DO 10,I=0,2
					IF (IRANK.EQ.0) THEN
					IORD = 2**I
					IF (M.EQ.-1) THEN
						XL = X1 - BND
						XU = BIG


						CALL PNDLGA(F,X1,1,XL,XU,UH,FEPS,IORD,IPRINT,
     &					           A1,NOCA,IERR)



						XL = X2-BND
						XU = BIG


						CALL PNDLGA(F,X2,1,XL,XU,UH,FEPS,IORD,IPRINT,
     &					           B1,NOCA,IERR)


					ELSE IF (M.EQ.0) THEN


						CALL PNDLG(F,X1,1,IORD,A1)





						CALL PNDLG(F,X2,1,IORD,B1)


					ELSE
						XL = -BIG
						XU = X1+BND


						CALL PNDLGA(F,X1,1,XL,XU,UH,FEPS,IORD,IPRINT,
     &					           A1,NOCA,IERR)



						XL = -BIG
						XU = X2+BND


						CALL PNDLGA(F,X2,1,XL,XU,UH,FEPS,IORD,IPRINT,
     &					           B1,NOCA,IERR)


					END IF
					ENDIF
C					CALL PNDL_BARRIER()
					
					IF (IRANK.EQ.0) THEN
					CALL RELERR(A1,GA,EA1)
					CALL RELERR(B1,GB,EB1)
					IF (EA1.GT.FAIL) THEN
						STATA = 'FAIL'
						NFAIL = NFAIL+1
					ELSE
						STATA = 'OK'
					END IF
					IF (EB1.GT.FAIL) THEN
						STATB = 'FAIL'
						NFAIL = NFAIL+1
					ELSE
						STATB = 'OK'
					END IF
					WRITE (*,100) IORD, A1, EA1, STATA, B1, EB1, STATB
					ENDIF
10				CONTINUE
				IF (IRANK.EQ.0) THEN
				WRITE (*,110) 'Exact:', GA, GB
				ENDIF
20			CONTINUE
30		CONTINUE
40	CONTINUE
C
C  Print the overall failure counter.
	IF (IRANK.EQ.0) THEN
	WRITE (*,*)
	WRITE (*,300) NFAIL
	WRITE (*,310) FAIL
	WRITE (*,*)
	ENDIF
C
	CALL TORC_FINALIZE()
	CALL PNDL_FINALIZE()
C	CALL MPI_FINALIZE()
C
100	FORMAT ('  ',I1,'  :',1X,1PG21.14,1X,1PE7.1,1X,A4,3X,
     &                            1PG21.14,1X,1PE7.1,1X,A4)
110	FORMAT (A6,1X,1PG21.14,16X,1PG21.14)
200	FORMAT ('Order',7X,'Derivative',7X,'Rel. error',
     &         10X,'Derivative',7X,'Rel. error')
210	FORMAT (14X,'X = ',F5.1,27X,'X = ',F5.1)
300	FORMAT (1X,'NUMBER OF FAILURES: ',I4)
310	FORMAT (1X,'FAILURE THRESHOLD:  ',1PE8.1)
C
	END
C  ---------------------------------------------------------------------
	SUBROUTINE RELERR ( FN, FA, RE )
C  ---------------------------------------------------------------------
C
C  Description:
C    Computes the relative error between a numericaly and an analyticaly
C    evaluated derivative.
C
C  Input arguments:
C    FN         Numerical value.
C    FA         Analytical value.
C
C  Output arguments:
C    RE         The relative error.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	RE = ABS(FN-FA)/MAX(1.D0,ABS(FA))
	END
C  ---------------------------------------------------------------------
	FUNCTION F ( X, N )
C  ---------------------------------------------------------------------
C
C  Description:
C    The test function. According to a selector (ISEL), it may assume
C    the value of several different functions.
C
C  Input arguments:
C    X          Array containing the values of the variables.
C    N          Number of variables.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N)
	COMMON / SELECT / ISEL
C
C	CALL TORC_SLEEP(100)
	IF (ISEL.EQ.1) THEN
		F = SIN(X(1))
	ELSE IF (ISEL.EQ.2) THEN
		F = EXP(X(1))
	ELSE IF (ISEL.EQ.3) THEN
		F = X(1)**2*SIN(X(1))
	ELSE IF (ISEL.EQ.4) THEN
		F = 7.D0*EXP(-X(1)**2/8.0D0)
	ELSE
		WRITE (*,*) 'ERROR: Incorrect ISEL'
		STOP
	END IF
C
	END
C  ---------------------------------------------------------------------
	SUBROUTINE GRAD ( X, N, G )
C  ---------------------------------------------------------------------
C
C  Description:
C    Computes the analytic gradient vector (G). According to a
C    selector (ISEL), it calculates the gradient of several different
C    functions.
C
C  Input arguments:
C    X          Array containing the values of the variables.
C    N          Number of variables.
C
C  Output arguments:
C    G          Array contating the first partial derivatives of the
C               test function.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), G(N)
	COMMON / SELECT / ISEL
C
	IF (ISEL.EQ.1) THEN
		G(1) = COS(X(1))
	ELSE IF (ISEL.EQ.2) THEN
		G(1) = EXP(X(1))
	ELSE IF (ISEL.EQ.3) THEN
		G(1) = 2.0D0*X(1)*SIN(X(1))+X(1)**2*COS(X(1))
	ELSE IF (ISEL.EQ.4) THEN
		G(1) = -7.0D0*X(1)*EXP(-X(1)**2/8.0D0)/4.0d0
	ELSE
		WRITE (*,*) 'ERROR: Incorrect ISEL'
		STOP
	END IF
C
	END
C  ---------------------------------------------------------------------
	SUBROUTINE NAMEIT ( FNAME )
C  ---------------------------------------------------------------------
C
C  Description:
C    Assigns a name (or description/formula) to each of the test 
C    functions, according to a selector (ISEL).
C
C  Output arguments:
C    FNAME      The name assigned to the function.
C
C  ---------------------------------------------------------------------
	CHARACTER*(*) FNAME
	COMMON / SELECT / ISEL
C
	IF (ISEL.EQ.1) THEN
		FNAME = 'SIN(X)'
	ELSE IF (ISEL.EQ.2) THEN
		FNAME = 'EXP(X)'
	ELSE IF (ISEL.EQ.3) THEN
		FNAME = 'X**2 * SIN(X)'
	ELSE IF (ISEL.EQ.4) THEN
		FNAME = '7 * EXP(-X**2/8)'
	ELSE
		WRITE (*,*) 'ERROR: Incorrect ISEL'
		STOP
	END IF
C
	END
C  ---------------------------------------------------------------------
	SUBROUTINE WLINE ( CH )
C  ---------------------------------------------------------------------
C
C  Description:
C    Prints a line on the screen using the character specified in
C    variable CH.
C
C  Input arguments:
C    CH         The character to be used for printing the line.
C
C  ---------------------------------------------------------------------
	CHARACTER CH
C
	WRITE (*,10) (CH,I=1,79)
10	FORMAT (79A1)
	END
C  ---------------------------------------------------------------------
	SUBROUTINE CPRINT ( MSG )
C  ---------------------------------------------------------------------
C
C  Description:
C    Prints a message (MSG) centered, on an 80 character line.
C
C  Input arguments:
C    MSG        The message to be printed.
C
C  ---------------------------------------------------------------------
C  Arguments:
	CHARACTER MSG*(*)
C  ---------------------------------------------------------------------
C  Local variables:
	CHARACTER LINE*80
C  ---------------------------------------------------------------------
C
	LE = LENGTH(MSG)
	IB = (74-LE)/2
	IBLE = IB+LE
	LINE = ' '
	LINE(IB+1:IBLE) = MSG(1:LE)
	WRITE (*,'(A)') LINE(1:IBLE+5)
C
	END
C  ---------------------------------------------------------------------
	FUNCTION LENGTH ( STRING )
C  ---------------------------------------------------------------------
C
C  Description:
C    Returns the position of the last non-blank character in
C    the input string, ie. its effective length.
C
C  Input arguments:
C    STRING    A character variable whose length we seek.
C
C  ---------------------------------------------------------------------
C  Arguments:
	CHARACTER STRING*(*)
C  ---------------------------------------------------------------------
C  Local variables:
	CHARACTER C
C  ---------------------------------------------------------------------
C
	DO 10,I=LEN(STRING),1,-1
		C = STRING(I:I)
		IF (C.NE.' ') THEN
			LENGTH = I
			RETURN
		END IF
10	CONTINUE
	LENGTH = 0
C
	END
