	PROGRAM TEST2
C  ---------------------------------------------------------------------
C  PNDL test program.
C  -----------------
C  Computes the second derivative matrix (Hessian) of several
C  two-dimensional functions, with all possible orders of accuracy.
C  The computation takes place at several different points, taking
C  into account all possible situations of point-bound proximity.
C  Relative errors against the analyticaly known derivatives are
C  evaluated and printed. Since the amount of output is quite large,
C  an overall failure indication is printed, based on the comparison
C  of each relative error with the FAIL threshold (defined below).
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INCLUDE 'torcf.h'
	EXTERNAL F, GRAD
C
C  Number of test points.
	PARAMETER ( NPT = 1 )
C
C  Number of test functions.
	PARAMETER ( NFUN = 1 )
C
C  Derivative computations with relative errors above
C  this threshold are considered "FAILed".
	PARAMETER ( FAIL = 1.0D-4 )
C
C  BND is the distance of a bound to the evaluation point
C  (in order to enforce point-bound proximity).
	PARAMETER ( BND = 1.0D-9 )
C
C  Printout control option.
	PARAMETER ( IPRINT = 0 )
C
C  Dimensionality of the test functions.
	PARAMETER ( N = 2 )
C
C  A large number, to be used as lower/upper bound.
	PARAMETER ( BIG = 1.0D300 )
C
	DIMENSION X(N,NPT), XL(N), XU(N), UH(N)
	DIMENSION H(N,N), HA(N,N)
	CHARACTER CFG
	CHARACTER*4 STAT
	CHARACTER*40 FNAME
	COMMON / SELECT / ISEL
C
C  The points at which the derivatives will be tested.
	DATA X / -1.0D0, 1.0D0 /
C
C	CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, IPROV, IERR)
	CALL PNDL_INIT()
	CALL TORC_REGISTER_TASK(GRAD)
	CALL TORC_INITF()

c	CALL MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
	IRANK = 0

C  The relative error in the evaluation of our test functions is
C  in the order of the machine accuracy.
	FEPS = 0.0D0
C  No stepsizes specified.
	UH(1) = 0.0D0
	UH(2) = 0.0D0
C
C  NFAIL is incremented each time the relative error is above
C  the "FAIL" threshold,
	NFAIL = 0
C
C  Loop over all point-bound proximity situations for the 1st variable.
	DO 40,M1=-1,1
C
C  Loop over all point-bound proximity situations for the 2nd variable.
		DO 50,M2=-1,1
			IF (IRANK.EQ.0) THEN
			WRITE (*,*)
			CALL WLINE('=')
			IF (M1.EQ.-1) THEN
				CALL CPRINT('VARIABLE X1 IS NEAR THE LOWER BOUND')
			ELSE IF (M1.EQ.0) THEN
				CALL CPRINT('NO BOUNDS ON VARIABLE X1')
			ELSE IF (M1.EQ.1) THEN
				CALL CPRINT('VARIABLE X1 IS NEAR THE UPPER BOUND')
			END IF
			IF (M2.EQ.-1) THEN
				CALL CPRINT('VARIABLE X2 IS NEAR THE LOWER BOUND')
			ELSE IF (M2.EQ.0) THEN
				CALL CPRINT('NO BOUNDS ON VARIABLE X2')
			ELSE IF (M2.EQ.1) THEN
				CALL CPRINT('VARIABLE X2 IS NEAR THE UPPER BOUND')
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
					X1 = X(1,K)
					X2 = X(2,K)
					CALL WLINE('-')
					CALL CPRINT('F(X) = '//FNAME)
					WRITE (*,210) X1, X2
					CALL WLINE('-')
					WRITE (*,200)
					CALL WLINE('-')
					CALL HESS(X(1,K),N,HA,N)
					ENDIF
C
C  Test using function/gradient values.
					DO 15,IFG=1,2
C
C  Loop over all available orders of accuracy.
						DO 10,I=0,1
							IF (IRANK.EQ.0) THEN
							IORD = 2**I
							IF (M1.EQ.-1) THEN
								XL(1) = X1 - BND
								XU(1) = BIG
							ELSE IF (M1.EQ.0) THEN
								XL(1) = -BIG
								XU(1) = BIG
							ELSE
								XL(1) = -BIG
								XU(1) = X1+BND
							END IF
							IF (M2.EQ.-1) THEN
								XL(2) = X2 - BND
								XU(2) = BIG
							ELSE IF (M2.EQ.0) THEN
								XL(2) = -BIG
								XU(2) = BIG
							ELSE
								XL(2) = -BIG
								XU(2) = X2+BND
							END IF
							IF (IFG.EQ.1) THEN


							CALL PNDLHFA(F,X(1,K),N,XL,XU,UH,FEPS,
     &							       IORD,IPRINT,H,N,NOC,IERR)


								CFG = 'F'
							ELSE


							CALL PNDLHGA(GRAD,X(1,K),N,XL,XU,UH,FEPS,
     &							       IORD,IPRINT,H,N,NOC,IERR)


								CFG = 'G'
							END IF
							ENDIF
C							CALL PNDL_BARRIER()
							IF (IRANK.EQ.0) THEN
							CALL RELERR(H(1,1),HA(1,1),E11)
							CALL RELERR(H(2,1),HA(2,1),E21)
							CALL RELERR(H(2,2),HA(2,2),E22)
							RE = MAX(E11,E21,E22)
							IF (RE.GT.FAIL) THEN
								STAT = 'FAIL'
								NFAIL = NFAIL+1
							ELSE
								STAT = 'OK'
							END IF
							WRITE (*,100) CFG, IORD, H(1,1), H(2,1),
     &						              H(2,2), RE, STAT
							ENDIF
10						CONTINUE
15					CONTINUE
					IF (IRANK.EQ.0) THEN
					WRITE (*,110) 'Exact:', HA(1,1), HA(2,1), HA(2,2)
					ENDIF
20				CONTINUE
30			CONTINUE
50		CONTINUE
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
100	FORMAT (' ',A,I1,'  :',1X,3(1PG19.12,1X),1PE7.1,1X,A4)
110	FORMAT (A6,1X,3(1PG19.12,1X))
200	FORMAT ('Order',7X,'H(1,1)',14X,'H(2,1)',
     &         14X,'H(2,2)',6X,'Max. rel. error')
210	FORMAT (14X,'X1 = ',F5.1,27X,'X2 = ',F5.1)
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
		F = X(1)**2*SIN(X(2))+EXP(X(1))*COS(X(2))**2
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
C	CALL TORC_SLEEP(100)
	IF (ISEL.EQ.1) THEN
		G(1) = 2.0D0*X(1)*SIN(X(2))+EXP(X(1))*COS(X(2))**2
		G(2) = X(1)**2*COS(X(2))-2.0D0*EXP(X(1))*COS(X(2))*SIN(X(2))
	ELSE
		WRITE (*,*) 'ERROR: Incorrect ISEL'
		STOP
	END IF
C
	END
C  ---------------------------------------------------------------------
	SUBROUTINE HESS ( X, N, HES, LD )
C  ---------------------------------------------------------------------
C
C  Description:
C    Computes the analytical  Hessian matrix (HES). According to a
C    selector (ISEL), it calculates the Hessian of several different
C    functions.
C
C  Input arguments:
C    X          Array containing the values of the variables.
C    N          Number of variables.
C    LD         Leading dimension of matrix HES.
C
C  Output arguments:
C    HES        Array containing the Hessian matrix. Note that only
C               the lower triangular part (plus the diagonal elements)
C               is returned.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), HES(LD,N)
	COMMON / SELECT / ISEL
C
	IF (ISEL.EQ.1) THEN
		HES(1,1) = 2.0D0*SIN(X(2))+EXP(X(1))*COS(X(2))**2
		HES(2,2) = -X(1)**2*SIN(X(2))
     &	           -2.0D0*EXP(X(1))*(COS(X(2))**2-SIN(X(2))**2)
		HES(2,1) = 2.0D0*X(1)*COS(X(2))
     &	          -2.0D0*EXP(X(1))*COS(X(2))*SIN(X(2))
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
		FNAME = 'X1**2*SIN(X2) + EXP(X1)*COS(X2)**2'
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
