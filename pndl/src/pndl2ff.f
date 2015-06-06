C  ---------------------------------------------------------------------
	SUBROUTINE PNDL2FF ( F, X, N, XL, XU, UH, FEPS, IPRINT, HES, LD,
     &                    NOC)
C, ISTIP, WTIP )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the Hessian matrix using O(h_i)+O(h_j)
C    formulae.
C    Selection of the appropriate formula is based on the proximity
C    of the variable with its bounds. Care is taken to avoid function
C    evaluations at the same point.
C
C  Input arguments:
C    F          The function to be differentiated. Must be declared as:
C                 FUNCTION F ( X, N )
C                 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                 DIMENSION X(N)
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
C                             steps etc).
C    LD         Leading dimension of matrix HES.
C
C  Output arguments:
C    HES        Array containing the resulting Hessian matrix.
C               Note that only the lower triangular part (plus the
C               diagonal elements) is returned.
C    NOC        Number of calls to the function F.
C
C  Work spaces:
C    ISTIP      Integer work space of length N.
C    WTIP       Real work space of length N.
C
C    Notes:
C      A simple caching mechanism is implemented, to avoid uneccessary
C      function evaluations at the same points. This might happen
C      when one ore more variables are near the bounds. The caching
C      mechanism uses the folowing arrays:
C        ISTIP(I) =  0 -> No value in the cache.
C                 =  1 -> F(X_I+) is available in the cache.
C                 = -1 -> F(X_I-) is available in the cache.
C        WTIP(I)  = F(X_I+) when ISTIP(I)=1
C                 = F(X_I-) when ISTIP(I)=-1
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), XL(N), XU(N), UH(N), HES(LD,N)
C
	DIMENSION ISTIP(N), WTIP(N)
C
	DIMENSION XX(N)
	DIMENSION FV(8*N*N)
	INTEGER FC
	INTEGER K
C                                                  
        EXTERNAL F, SPAWNTASK
C
	LOGICAL ILB, IUB, JLB, JUB, PR2, PR3
	CHARACTER*5 FORM
C
	DATA ONE   / 1.0D0 /
	DATA THREE / 3.0D0 /
C
	DO K=1,N
	XX(K) = X(K)
	END DO
C
	CON = FEPS**(ONE/THREE)
	
C	!! FIRST PASS !!
	FC = 1
C  Function call counter.
	NOC = 1
C
	PR2 = .FALSE.
	PR3 = .FALSE.
C
	CALL SPAWNTASK(F, XX, N, FV(FC))
	FC = FC + 1

C
C  Calculate the diagonal elements.
	IF (PR3) WRITE (*,30)
	DO 100,I=1,N
		ISTIP(I) = 0
		TI = XX(I)
		IF (UH(I).EQ.0.0D0) THEN
			HI = CON*MAX(ABS(TI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		TIP = TI+HI
		TIM = TI-HI
		TIP2 = TI+2.0D0*HI
		TIM2 = TI-2.0D0*HI
C  ---------------------------------------------------------------------
C  X_I++ is below the upper bound.
		IF (TIP2.LT.XU(I)) THEN
			XX(I) = TIP
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1
			ISTIP(I) = 1
			WTIP(I) = F1
C
			XX(I) = TIP2
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1
C
			NOC = NOC + 2
			HES(I,I) = (F00 - 2.0D0*F1 + F2)/(HI**2)
			FORM = 'FD-FD'
C  ---------------------------------------------------------------------
C  X_I-- is above the lower bound.
		ELSE IF (TIM2.GT.XL(I)) THEN
			XX(I) = TIM
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1
			ISTIP(I) = -1
			WTIP(I) = F1
C
			XX(I) = TIM2
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1
C
			NOC = NOC + 2
			HES(I,I) = (F00 - 2.0D0*F1 + F2)/(HI**2)
			FORM = 'BD-BD'
		ELSE
C  ---------------------------------------------------------------------
C  X_I is confined in a very small interval.
			HES(I,I) = 0.0D0
			FORM = '     '
			IF (PR2) WRITE (*,50) I
		END IF
		XX(I) = TI
		IF (PR3) WRITE (*,40) I, I, FORM, HI, HI, HES(I,I)
100	CONTINUE
C
C  Calculate off-diagonal elements.
	DO 10,J=1,N-1
		TJ = XX(J)
		IF (UH(J).EQ.0.0D0) THEN
			HJ = CON*MAX(ABS(TJ),ONE)
		ELSE
			HJ = UH(J)
		ENDIF
		TJP = TJ+HJ
		TJM = TJ-HJ
		JLB = TJM.GT.XL(J)
		JUB = TJP.LT.XU(J)
		DO 20,I=J+1,N
			TI = XX(I)
			IF (UH(I).EQ.0.0D0) THEN
				HI = CON*MAX(ABS(TI),ONE)
			ELSE
				HI = UH(I)
			ENDIF
			TIP = TI+HI
			TIM = TI-HI
			ILB = TIM.GT.XL(I)
			IUB = TIP.LT.XU(I)
C  ---------------------------------------------------------------------
C  X_I+ below the upper bound, X_J+ below the upper bound.
			IF (IUB .AND. JUB) THEN
				IF (ISTIP(I).EQ.1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				XX(I) = TIP
				XX(J) = TJP
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				NOC = NOC + 1
C
			HES(I,J) = (F00 - F10 - F01 + F11)/(HI*HJ)
				FORM = 'FD-FD'
C  ---------------------------------------------------------------------
C  X_I+ below the upper bound, X_J- above the lower bound.
			ELSE IF (IUB .AND. JLB) THEN
				IF (ISTIP(I).EQ.1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.-1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				XX(I) = TIP
				XX(J) = TJM
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				NOC = NOC + 1
C
			HES(I,J) = (-F00 + F10 + F01 - F11)/(HI*HJ)
				FORM = 'FD-BD'
C  ---------------------------------------------------------------------
C  X_I- above the lower bound, X_J+ below the upper bound.
			ELSE IF (ILB .AND. JUB) THEN
				IF (ISTIP(I).EQ.-1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				XX(I) = TIM
				XX(J) = TJP
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				NOC = NOC + 1
C
			HES(I,J) = (-F00 + F10 + F01 - F11)/(HI*HJ)
				FORM = 'BD-FD'
C  ---------------------------------------------------------------------
C  X_I- above the lower bound, X_J- above the lower bound.
			ELSE IF (ILB .AND. JLB) THEN
				IF (ISTIP(I).EQ.-1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.-1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC + 1
				END IF
C
				XX(I) = TIM
				XX(J) = TJM
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				NOC = NOC + 1
C
			HES(I,J) = (F00 - F10 - F01 + F11)/(HI*HJ)
				FORM = 'BD-BD'
			ELSE
C  ---------------------------------------------------------------------
C  X_I or X_J is confined in a very small interval.
				HES(I,J) = 0.0D0
				FORM = '     '
				IF (PR3) WRITE (*,60) I, J
			END IF
			XX(I) = TI
			XX(J) = TJ
			IF (PR3) WRITE (*,40) I, J, FORM, HI, HJ, HES(I,J)
20		CONTINUE
10	CONTINUE
C

c	CALL torc_enable_stealing()
	CALL torc_waitall()
c	CALL torc_disable_stealing()

C	!! SECOND PASS !!
	FC = 1
C  Function call counter.
c	NOC = 1
C
	PR2 = IPRINT.GE.2
	PR3 = IPRINT.EQ.3
C
	F00 = FV(FC)
	FC = FC + 1
C
C  Calculate the diagonal elements.
	IF (PR3) WRITE (*,30)
	DO 101,I=1,N
		ISTIP(I) = 0
		TI = XX(I)
		IF (UH(I).EQ.0.0D0) THEN
			HI = CON*MAX(ABS(TI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		TIP = TI+HI
		TIM = TI-HI
		TIP2 = TI+2.0D0*HI
		TIM2 = TI-2.0D0*HI
C  ---------------------------------------------------------------------
C  X_I++ is below the upper bound.
		IF (TIP2.LT.XU(I)) THEN
			XX(I) = TIP
			F1 = FV(FC)
			FC = FC + 1
			ISTIP(I) = 1
			WTIP(I) = F1
C
			XX(I) = TIP2
			F2 = FV(FC)
			FC = FC + 1
C
c			NOC = NOC + 2
			HES(I,I) = (F00 - 2.0D0*F1 + F2)/(HI**2)
			FORM = 'FD-FD'
C  ---------------------------------------------------------------------
C  X_I-- is above the lower bound.
		ELSE IF (TIM2.GT.XL(I)) THEN
			XX(I) = TIM
			F1 = FV(FC)
			FC = FC + 1
			ISTIP(I) = -1
			WTIP(I) = F1
C
			XX(I) = TIM2
			F2 = FV(FC)
			FC = FC + 1
C
c			NOC = NOC + 2
			HES(I,I) = (F00 - 2.0D0*F1 + F2)/(HI**2)
			FORM = 'BD-BD'
		ELSE
C  ---------------------------------------------------------------------
C  X_I is confined in a very small interval.
			HES(I,I) = 0.0D0
			FORM = '     '
			IF (PR2) WRITE (*,50) I
		END IF
		XX(I) = TI
		IF (PR3) WRITE (*,40) I, I, FORM, HI, HI, HES(I,I)
101	CONTINUE
C
C  Calculate off-diagonal elements.
	DO 11,J=1,N-1
		TJ = XX(J)
		IF (UH(J).EQ.0.0D0) THEN
			HJ = CON*MAX(ABS(TJ),ONE)
		ELSE
			HJ = UH(J)
		ENDIF
		TJP = TJ+HJ
		TJM = TJ-HJ
		JLB = TJM.GT.XL(J)
		JUB = TJP.LT.XU(J)
		DO 21,I=J+1,N
			TI = XX(I)
			IF (UH(I).EQ.0.0D0) THEN
				HI = CON*MAX(ABS(TI),ONE)
			ELSE
				HI = UH(I)
			ENDIF
			TIP = TI+HI
			TIM = TI-HI
			ILB = TIM.GT.XL(I)
			IUB = TIP.LT.XU(I)
C  ---------------------------------------------------------------------
C  X_I+ below the upper bound, X_J+ below the upper bound.
			IF (IUB .AND. JUB) THEN
				IF (ISTIP(I).EQ.1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIP
					F10 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJP
					F01 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				XX(I) = TIP
				XX(J) = TJP
				F11 = FV(FC)
				FC = FC + 1
c				NOC = NOC + 1
C
			HES(I,J) = (F00 - F10 - F01 + F11)/(HI*HJ)
				FORM = 'FD-FD'
C  ---------------------------------------------------------------------
C  X_I+ below the upper bound, X_J- above the lower bound.
			ELSE IF (IUB .AND. JLB) THEN
				IF (ISTIP(I).EQ.1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIP
					F10 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.-1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJM
					F01 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				XX(I) = TIP
				XX(J) = TJM
				F11 = FV(FC)
				FC = FC + 1
c				NOC = NOC + 1
C
			HES(I,J) = (-F00 + F10 + F01 - F11)/(HI*HJ)
				FORM = 'FD-BD'
C  ---------------------------------------------------------------------
C  X_I- above the lower bound, X_J+ below the upper bound.
			ELSE IF (ILB .AND. JUB) THEN
				IF (ISTIP(I).EQ.-1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIM
					F10 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJP
					F01 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				XX(I) = TIM
				XX(J) = TJP
				F11 = FV(FC)
				FC = FC + 1
c				NOC = NOC + 1
C
			HES(I,J) = (-F00 + F10 + F01 - F11)/(HI*HJ)
				FORM = 'BD-FD'
C  ---------------------------------------------------------------------
C  X_I- above the lower bound, X_J- above the lower bound.
			ELSE IF (ILB .AND. JLB) THEN
				IF (ISTIP(I).EQ.-1) THEN
					F10 = WTIP(I)
				ELSE
					XX(I) = TIM
					F10 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				IF (ISTIP(J).EQ.-1) THEN
					F01 = WTIP(J)
				ELSE
					XX(I) = TI
					XX(J) = TJM
					F01 = FV(FC)
					FC = FC + 1
c					NOC = NOC + 1
				END IF
C
				XX(I) = TIM
				XX(J) = TJM
				F11 = FV(FC)
				FC = FC + 1
c				NOC = NOC + 1
C
			HES(I,J) = (F00 - F10 - F01 + F11)/(HI*HJ)
				FORM = 'BD-BD'
			ELSE
C  ---------------------------------------------------------------------
C  X_I or X_J is confined in a very small interval.
				HES(I,J) = 0.0D0
				FORM = '     '
				IF (PR3) WRITE (*,60) I, J
			END IF
			XX(I) = TI
			XX(J) = TJ
			IF (PR3) WRITE (*,40) I, J, FORM, HI, HJ, HES(I,J)
21		CONTINUE
11	CONTINUE

C
30	FORMAT (/' PNDL:',3X,'I',6X,'J',1X,'Formula',7X,'Step_i',13X,
     &        'Step_j',12X,'HES_ij')
40	FORMAT (' PNDL:',I4,3X,I4,2X,A,2X,1PG18.11,1X,1PG18.11,1X,
     &        1PG18.11)
50	FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
60	FORMAT (' PNDL: Warning: Variable ',I4,' or ',I4,' is confined ',
     &        'in a very small interval.'
     &        /' PNDL:          Setting derivative to zero.')
	END
