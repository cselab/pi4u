C  ---------------------------------------------------------------------
	SUBROUTINE PNDL2GQ ( GRD, X, N, XL, XU, UH, FEPS, IPRINT, HES,
     &                    LD, IGC )
C     &                    LD, GW, G, IWFIX, IGC )
C  ---------------------------------------------------------------------
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the Hessian matrix using O(h_i*h_j)
C    formulae. Selection of the appropriate formula is based on the
C    proximity of the variable with its bounds.
C    Instead of calling the function directly, this routine uses the
C    first partial derivatives of the function.
C    Care is taken to avoid function evaluations at the same point.
C
C  Input arguments:
C    GRD        A subroutine that returns the gradient vector (G),
C               given the values of the variables (X).
C               It must be declared as:
C                 SUBROUTINE GRD ( X, N, G )
C                 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                 DIMENSION X(N), G(N)
C    X          Array containing the function variables.
C    N          The number of variables.
C    XL         Array containing lower bounds on the variables.
C    XU         Array containing upper bounds on the variables.
C    UH         Array containing user-speficied steps.
C    FEPS       An estimation of the relative accuracy with which
C               the function is evaluated. FEPS must be non-zero.
C    IPRINT     This option controls the amount of printout from the
C               routine. Note that all output appears on the standard
C               output device. Possible values are:
C                 IPRINT=0 -> No printout at all.
C                 IPRINT=1 -> Fatal error messages are printed.
C                 IPRINT=2 -> Warning messages are printed.
C                 IPRINT=3 -> Detailed information is printed (the
C                             formula that was used, differentiation
C                             steps and the resulting gradient
C                             vector).
C    LD         Leading dimension of matrix HES.
C
C  Output arguments:
C    HES        Array containing the resulting Hessian matrix.
C               Note that only the lower triangular part (plus the
C               diagonal elements) is returned.
C    IGC        Number of calls to the function F.
C
C  Work spaces:
C    GW         Real work space of length NxN.
C    G          Real work space of length N.
C    IWFIX      Integer work space of length N.
C
C  Notes:
C    The main differentiation formula uses central differences, however
C    when one or more variables are near the bounds, forward or backward
C    difference formulae are used. The later require evaluation of the
C    gradient at the given point X. To avoid uneccessary evaluations
C    when more than one variable are near the bounds, a simple caching
C    mechanism is implemented. Variable IGAT indicates whether the
C    gradient has been evaluated at the given point X.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), XL(N), XU(N), UH(N), HES(LD,N)
C
	DIMENSION GW(N,N), G(N), IWFIX(N)
C
	DIMENSION XX(N)
	DIMENSION IFB(N)
	DIMENSION HH(N)
	INTEGER K
C                                                  
        EXTERNAL GRD, SPAWNTASKG
C
	LOGICAL ILB, IUB, PR2, PR3
	CHARACTER*2 FORM
C
	DATA ONE   / 1.D0 /
C
	DO K=1,N
	XX(K) = X(K)
	END DO
C
	CON = FEPS**(1.0D0/3.0D0)
C

C	!! FIRST PASS !!

C  Analytic gradient call counter.
	IGC = 0
C
C  Flag to indicate whether the gradient has been evaluated at the
C  current point. The gradient at the current point is neccessary 
C  only if forward/backward formulae are to be used.
	IGAT = 0
C
	PR2 = .FALSE.
	PR3 = .FALSE.

C  Loop over all variables.
	IF (PR3) WRITE (*,200)

	DO 100,I=1,N

		IWFIX(I) = 1
		TI = XX(I)
		IF (UH(I).EQ.0.0D0) THEN
			HI = CON*MAX(ABS(TI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		HH(I) = HI
		TIP = TI+HI
		TIM = TI-HI
		TWOHI = 2.0D0*HI
		TIP2 = TI + TWOHI
		TIM2 = TI - TWOHI
		ILB = TIM.GT.XL(I)
		IUB = TIP.LT.XU(I)
C
C  Basic case: the variable plus/minus the differentiation step is 
C  inside the margins. Use central differences.
		IF (ILB .AND. IUB) THEN
			XX(I) = TIP
			CALL SPAWNTASKG(GRD, XX, N, HES(1,I))
			IGC = IGC + 1

			XX(I) = TIM
			CALL SPAWNTASKG(GRD, XX, N, GW(1,I))
			IGC = IGC + 1

c			IFB(I) = 0
c			FORM = 'CD'
C
C  Cannot use central differences since the variable plus the 
C  differentiation step exceeds the upper bound. However we can 
C  use the three point backward difference formula.
		ELSE IF (.NOT.IUB .AND. TIM2.GT.XL(I)) THEN
C
C  We need the gradient at the current point. Check whether it has been
C  evaluated already.

			IF (IGAT.EQ.0) THEN
				CALL SPAWNTASKG(GRD, XX, N, G)
				IGC = IGC + 1
				IGAT = 1
			END IF

C
			XX(I) = TIM
			CALL SPAWNTASKG(GRD, XX, N, HES(1,I))
			IGC = IGC + 1


C
			XX(I) = TIM2
			CALL SPAWNTASKG(GRD, XX, N, GW(1,I))
			IGC = IGC + 1

c			IFB(I) = -1
c			FORM = 'BD'
C
C  Cannot use central or backward differences. However we can use the 
C  three point forward difference formula.
		ELSE IF (.NOT.ILB .AND. TIP2.LT.XU(I)) THEN
C  We need the gradient at the current point. Check whether it has been
C  evaluated already.

			IF (IGAT.EQ.0) THEN
				CALL SPAWNTASKG(GRD, XX, N, G)
				IGC = IGC + 1
				IGAT = 1
			END IF

C
			XX(I) = TIP
			CALL SPAWNTASKG(GRD, XX, N, HES(1,I))
			IGC = IGC + 1

C
			XX(I) = TIP2
			CALL SPAWNTASKG(GRD, XX, N, GW(1,I))
			IGC = IGC + 1

c			IFB(I) = 1
c			FORM = 'FD'
		ELSE
C  The variable is confined in a very small interval. (An "almost fixed" 
C  variable). Numerical differentiation is not possible.
			IWFIX(I) = 0
			FORM = '  '
			IFB(I) = 2
			IF (PR2) WRITE (*,220) I
		END IF
		XX(I) = TI
		IF (PR3) WRITE (*,210) I, FORM, TI, HI
100	CONTINUE

c	CALL torc_enable_stealing()
	CALL torc_waitall()
c	CALL torc_disable_stealing()

C	!! SECOND PASS !!
	PR2 = IPRINT.GE.2
	PR3 = IPRINT.EQ.3

C  Loop over all variables.
	IF (PR3) WRITE (*,200)

	DO 101,I=1,N

		IWFIX(I) = 1
		TI = XX(I)
		IF (UH(I).EQ.0.0D0) THEN
			HI = CON*MAX(ABS(TI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		HH(I) = HI
		TIP = TI+HI
		TIM = TI-HI
		TWOHI = 2.0D0*HI
		TIP2 = TI + TWOHI
		TIM2 = TI - TWOHI
		ILB = TIM.GT.XL(I)
		IUB = TIP.LT.XU(I)
C
C  Basic case: the variable plus/minus the differentiation step is 
C  inside the margins. Use central differences.
		IF (ILB .AND. IUB) THEN

			DO 20,J=1,N
20				HES(J,I) = (HES(J,I)-GW(J,I))/TWOHI

			IFB(I) = 0
			FORM = 'CD'
C
C  Cannot use central differences since the variable plus the 
C  differentiation step exceeds the upper bound. However we can 
C  use the three point backward difference formula.
		ELSE IF (.NOT.IUB .AND. TIM2.GT.XL(I)) THEN
C
C  We need the gradient at the current point. Check whether it has been
C  evaluated already.

			DO 50,J=1,N
50				HES(J,I) = (-HES(J,I) + 0.25D0*GW(J,I))



			IFB(I) = -1
			FORM = 'BD'
C
C  Cannot use central or backward differences. However we can use the 
C  three point forward difference formula.
		ELSE IF (.NOT.ILB .AND. TIP2.LT.XU(I)) THEN
C  We need the gradient at the current point. Check whether it has been
C  evaluated already.

			DO 90,J=1,N
90				HES(J,I) = (HES(J,I) - 0.25D0*GW(J,I))

			IFB(I) = 1
			FORM = 'FD'
		ELSE
C  The variable is confined in a very small interval. (An "almost fixed" 
C  variable). Numerical differentiation is not possible.
			IWFIX(I) = 0
			FORM = '  '
			IFB(I) = 2
			IF (PR2) WRITE (*,220) I
		END IF
		XX(I) = TI
		IF (PR3) WRITE (*,210) I, FORM, TI, HI
101	CONTINUE

C	!! THIRD PASS !!
C  Perform division using H()
	DO 2,I=1,N
		IF (IFB(I).EQ.-1) THEN
			HALFHI = 0.5D0*HH(I)
			DO 21,J=1,N
			HES(J,I) = ( 0.75D0*G(J)+HES(J,I))/HALFHI
21			CONTINUE
		ELSE IF (IFB(I).EQ.1) THEN
			HALFHI = 0.5D0*HH(I)
			DO 22,J=1,N
			HES(J,I) = (-0.75D0*G(J)+HES(J,I))/HALFHI
22			CONTINUE
		END IF
2	CONTINUE


C
C  Make sure the Hessian is symmetric. (The result is stored in
C  the lower triangular part).
	DO 110,J=1,N
		IWJ = IWFIX(J)
		DO 120,I=J,N
		  HES(I,J) = IWJ*IWFIX(I)*(HES(I,J)+HES(J,I))/2.0D0
120		CONTINUE
110	CONTINUE
C

200	FORMAT (/' PNDL:',' Index',1X,'Formula',7X,'X_i',19X,'Step_i')
210	FORMAT (' PNDL:',I6,3X,A,3X,1PG21.14,1X,1PG21.14)
220	FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
	END 
