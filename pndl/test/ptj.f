	PROGRAM TESTJ
C  ---------------------------------------------------------------------
C  PNDL test program.
C  -----------------
C  Computes the Jacobian matrix of several one-dimensional
C  functions, with all possible orders of accuracy.
C  The computation takes place at several different points, taking
C  into account all possible situations of point-bound proximity.
C  Relative errors against the analyticaly known derivatives are
C  evaluated and printed. Since the amount of output is quite large,
C  an overall failure indication is printed, based on the comparison
C  of each relative error with the FAIL threshold (defined below).
C  The test functions are of the form:
C                  M
C          F(X) = Sum R_i(X)**2
C                 i=1
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INCLUDE 'torcf.h'
	EXTERNAL RSD
C
C  Number of test points.
	PARAMETER ( NPT = 5 )
C
C  Number of test functions.
	PARAMETER ( NFUN = 1 )
C
C  The number of squared terms in each test function.
	PARAMETER ( MR = 3 )
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
	DIMENSION FJ(MR,1), FJA(MR,1)
	CHARACTER*4 STAT
	CHARACTER*40 FNAME
	COMMON / SELECT / ISEL
C
C  The points at which the derivatives will be tested.
	DATA X / -4.0D0, -2.0D0, 0.0D0, 2.0D0, 4.0D0 /
C
C	CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, IPROV, IERR)
	CALL PNDL_INIT()
	CALL TORC_REGISTER_TASK(RSD)
	CALL TORC_INIT()

c	CALL MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
	IRANK = 0

C
C  The relative error in the evaluation of our test functions is
C  in the order of the machine accuracy.
	FEPS = 0.0D0
C No stepsize specified
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
			DO 20,K=1,NPT
				IF (IRANK.EQ.0) THEN
				X1 = X(K)
				CALL WLINE('-')
				CALL CPRINT(FNAME)
				WRITE (*,210) X1
				CALL WLINE('-')
				WRITE (*,200)
				CALL WLINE('-')
				CALL JACAN(X1,1,MR,FJA,MR)
				ENDIF
C
C  Loop over all available orders of accuracy.
				DO 10,I=0,1
					IF (IRANK.EQ.0) THEN
					IORD = 2**I
					IF (M.EQ.-1) THEN
						XL = X1 - BND
						XU = BIG
					ELSE IF (M.EQ.0) THEN
						XL = -BIG
						XU = BIG
					ELSE
						XL = -BIG
						XU = X1+BND
					END IF

					CALL PNDLJA(RSD,X1,1,MR,XL,XU,UH,FEPS,IORD,
     &				           IPRINT,FJ,MR,NOC,IERR)
					ENDIF
C					CALL PNDL_BARRIER()

					IF (IRANK.EQ.0) THEN
					CALL RELERR(FJ(1,1),FJA(1,1),E1)
					CALL RELERR(FJ(2,1),FJA(2,1),E2)
					CALL RELERR(FJ(3,1),FJA(3,1),E3)
					RE = MAX(E1,E2,E3)
					IF (RE.GT.FAIL) THEN
						STAT = 'FAIL'
						NFAIL = NFAIL+1
					ELSE
						STAT = 'OK'
					END IF
					WRITE (*,100) IORD, FJ(1,1), FJ(2,1), FJ(3,1),
     &				              RE, STAT
					ENDIF
10				CONTINUE
				IF (IRANK.EQ.0) THEN
				WRITE (*,110) 'Exact:', FJA(1,1), FJA(2,1), FJA(3,1)
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
C	CALL MPI_FINALIZE()
C
100	FORMAT ('  ',I1,'  :',1X,3(1PG19.12,1X),1PE7.1,1X,A4)
110	FORMAT (A6,1X,3(1PG19.12,1X))
200	FORMAT ('Order',8X,'J(1,1)',14X,'J(2,1)',
     &         13X,'J(3,1)',6X,'Max. rel. error')
210	FORMAT (34X,'X = ',F5.1)
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
	SUBROUTINE RSD ( X, N, M, R )
C  ---------------------------------------------------------------------
C
C  Description:
C    Implements the test function which has the form:
C                  M
C          F(X) = Sum R_i(X)**2
C                 i=1
C    According to a selector (ISEL), it may assume the value of several
C    different functions.
C    Note that instead of a single value, the values of the residuls
C    R_i are returned.
C
C  Input arguments:
C    X          Array containing the values of the variables.
C    N          Number of variables.
C    M          Number of squared terms.
C
C  Output arguments:
C    R          Array containing the values of the residuals.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), R(M)
	COMMON / SELECT / ISEL
C
C	CALL TORC_SLEEP(100)
	IF (ISEL.EQ.1) THEN
		R(1) = SIN(X(1))
		R(2) = EXP(X(1))
		R(3) = X(1)*COS(X(1))
	ELSE
		WRITE (*,*) 'ERROR: Incorrect ISEL'
		STOP
	END IF
C
	END
C  ---------------------------------------------------------------------
	SUBROUTINE JACAN ( X, N, M, FJ, LD )
C  ---------------------------------------------------------------------
C  Description:
C    Returns the analytic value of the Jacobian matrix. It may return 
C    the Jacobian of several different functions according to a
C    selector (ISEL).
C
C  Input arguments:
C    X          Array containing the values of the variables.
C    N          Number of variables.
C    M          Number of squared terms.
C    LD         LEading dimension of the JAcobian (FJ).
C
C  Output arguments:
C    FJ         The Jacobian matrix. FJ(i,j) = dRi/dXj
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), FJ(LD,N)
	COMMON / SELECT / ISEL
C
	IF (ISEL.EQ.1) THEN
		FJ(1,1) = COS(X(1))
		FJ(2,1) = EXP(X(1))
		FJ(3,1) = COS(X(1))-X(1)*SIN(X(1))
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
		FNAME = 'F1=SIN(X)   F2=EXP(X)   F3=X*COS(X)'
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
