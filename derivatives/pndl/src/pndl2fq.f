C  ---------------------------------------------------------------------
	SUBROUTINE PNDL2FQ ( F, X, N, XL, XU, UH, FEPS, IPRINT, HES, LD,
     &                    NOC)
C    ISTIP, WTIP, WTIP2 )
C  ---------------------------------------------------------------------
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the Hessian matrix using O(h_i*h_j)
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
C    WTIP2      Real work space of length N.
C
C    Notes:
C      A simple caching mechanism is implemented, to avoid uneccessary
C      function evaluations at the same points. This might happen
C      when one ore more variables are near the bounds. The caching
C      mechanism uses the folowing arrays:
C        ISTIP(I) =  0 -> No value in the cache.
C                 =  1 -> F(X_I+), F(X_I++) are available in the cache.
C                 = -1 -> F(X_I-), F(X_I--) are available in the cache.
C                 =  2 -> F(X_I+), F(X_I-) are available in the cache.
C        WTIP(I)  = F(X_I+)  when ISTIP(I)= 1  or ISTIP(I)=2
C                 = F(X_I-)  when ISTIP(I)=-1
C        WTIP2(I) = F(X_I++) when ISTIP(I)= 1
C                 = F(X_I--) when ISTIP(I)=-1
C                 = F(X_I-)  when ISTIP(I)= 2
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  Arguments:
	DIMENSION X(N), XL(N), XU(N), UH(N), HES(LD,N)
C
	DIMENSION ISTIP(N), WTIP(N), WTIP2(N)
C
	DIMENSION FV(4*N*N)
	INTEGER FC
	DIMENSION XX(N)
C                                                  
        EXTERNAL F, SPAWNTASK
C
	LOGICAL ILB, IUB, JLB, JUB, PR2, PR3
	CHARACTER*5 FORM
C
	DATA ONE  / 1.0D0 /
	DATA FOUR / 4.0D0 /
C
	DO K=1,N
	XX(K) = X(K)
	END DO
C
	CON = FEPS**(ONE/FOUR)
C
	PR2 = .FALSE.
	PR3 = .FALSE.

C	!! FIRST PASS !!
	FC = 1
	CALL SPAWNTASK(F, XX, N, FV(FC))
	FC = FC+1
C
C  Function call counter.
	NOC = 1
C
C  Calculate the diagonal elements.


	DO 200,I=1,N
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
		TIP3 = TI+3.0D0*HI
		TIM3 = TI-3.0D0*HI
C  ---------------------------------------------------------------------
C  X_I+, X_I- are inside the bounds.
		IF (TIM.GT.XL(I) .AND. TIP.LT.XU(I)) THEN



			XX(I) = TIP
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC+1

C
			XX(I) = TIM
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

C
			ISTIP(I) = 2
			WTIP(I) = FP
			WTIP2(I) = FM
C
			NOC = NOC + 2
			FORM = 'CD-CD'
C  ---------------------------------------------------------------------
C  X_I+++ below the upper bound.
		ELSE IF (TIP3.LT.XU(I)) THEN



			XX(I) = TIP
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

C
			XX(I)= TIP2
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

C
			XX(I)= TIP3
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1


C
			ISTIP(I) = 1
			WTIP(I) = F1
			WTIP2(I) = F2
C
			NOC = NOC + 3
			HES(I,I) = (2.0D0*F00-5.0D0*F1+4.0D0*F2-F3)/HI**2
			FORM = 'FD-FD'
C  ---------------------------------------------------------------------
C  X_I--- above the lower bound.
		ELSE IF (TIM3.GT.XL(I)) THEN



			XX(I) = TIM
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

C
			XX(I)= TIM2
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

C
			XX(I)= TIM3
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1


C
			ISTIP(I) = -1
			WTIP(I) = F1
			WTIP2(I) = F2
C
			NOC = NOC + 3
			HES(I,I) = (2.0D0*F00-5.0D0*F1+4.0D0*F2-F3)/HI**2
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
200	CONTINUE
C
C  Calculate off-diagonal elements.


	DO 10,J=1,N-1
		TJ = XX(J)
		IF (UH(J).EQ.0D0) THEN
			HJ = CON*MAX(ABS(TJ),ONE)
		ELSE
			HJ = UH(J)
		ENDIF
		TJP = TJ+HJ
		TJM = TJ-HJ
		TJP2 = TJ+2.0D0*HJ
		TJM2 = TJ-2.0D0*HJ
		JLB = TJM.GT.XL(J)
		JUB = TJP.LT.XU(J)
		DO 20,I=J+1,N
			TI = XX(I)
			IF (UH(I).EQ.0D0) THEN
				HI = CON*MAX(ABS(TI),ONE)
			ELSE
				HI = UH(I)
			ENDIF
			TIP = TI+HI
			TIM = TI-HI
			TIP2 = TI+2.0D0*HI
			TIM2 = TI-2.0D0*HI
			ILB = TIM.GT.XL(I)
			IUB = TIP.LT.XU(I)
C  ---------------------------------------------------------------------
C  X_I+ inside the bounds, X_J+ inside the bounds.
			IF (ILB .AND. JLB .AND. IUB .AND. JUB) THEN



				XX(I) = TIP
				XX(J) = TJP
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
C


				XX(J) = TJM
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
C

				XX(I) = TIM
				XX(J) = TJP
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
C

				XX(J) = TJM
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
C


				NOC = NOC + 4
				HES(I,J) = (FPP-FPM+FMM-FMP)/(4.0D0*HI*HJ)
				FORM = 'CD-CD'
C  ---------------------------------------------------------------------
C  X_I- below the lower bound, X_J- below the lower bound.
			ELSE IF (.NOT.ILB .AND. .NOT.JLB) THEN
C  X_I++ below the upper bound, X_J++ below the upper bound.
				IF (TIP2.LT.XU(I) .AND. TJP2.LT.XU(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01 = WTIP(J)
					ELSE
						XX(J) = TJP
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1

					END IF
C

					IF (IST.EQ.1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJP2
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1

					END IF
C

					XX(I) = TIP
					IST = ISTIP(I)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F10 = WTIP(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1

					END IF
C

					XX(J)= TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J)= TJP2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIP2
					IF (IST.EQ.1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1

					END IF
C

					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJP2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (9.0D0*F00+3.0D0*F02-12.0D0*F01
     &				            +16.0D0*F11-4.0D0*F12-12.0D0*F10
     &				            +3.0D0*F20+F22-4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'FD-FD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I,J
				END IF
C  ---------------------------------------------------------------------
C  X_I- below the lower bound,X_J+ above the upper bound.
			ELSE IF (.NOT.ILB .AND. .NOT.JUB) THEN
C  X_I++ below the upper bound,X_J-- above the lower bound.
				IF (TIP2.LT.XU(I) .AND. TJM2.GT.XL(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01 = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01 = WTIP2(J)
					ELSE
						XX(J) = TJM
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.-1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJM2
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(I) = TIP
					IST = ISTIP(I)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F10 = WTIP(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(J)= TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J)= TJM2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIP2
					IF (IST.EQ.1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJM2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (-9.0D0*F00-3.0D0*F02+12.0D0*F01
     &				            -16.0D0*F11+4.0D0*F12+12.0D0*F10
     &				            -3.0D0*F20-F22+4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'FD-BD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I,J
				END IF
C  ---------------------------------------------------------------------
C  X_I+ above the upper bound, X_J- below the lower bound.
			ELSE IF (.NOT.IUB .AND. .NOT.JLB) THEN
C  X_I-- above the lower bound, X_I++ below the upper bound.
				IF (TIM2.GT.XL(I) .AND. TJP2.LT.XU(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01 = WTIP(J)
					ELSE
						XX(J) = TJP
						XX(J) = TJP
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJP2
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F10 = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F10 = WTIP2(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(J)= TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J)= TJP2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIM2
					IF (IST.EQ.-1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJP2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (-9.0D0*F00-3.0D0*F02+12.0D0*F01
     &				            -16.0D0*F11+4.0D0*F12+12.0D0*F10
     &				            -3.0D0*F20-F22+4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'BD-FD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+ above the upper bound, X_J+ above the upper bound.
			ELSE IF (.NOT.IUB .AND. .NOT.JUB) THEN
C  X_I-- above the lower bound, X_J-- above the lower bound.
				IF (TIM2.GT.XL(I) .AND. TJM2.GT.XL(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01 = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01 = WTIP2(J)
					ELSE
						XX(J) = TJM
						FV(FC) = F(X,N)
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.-1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJM2
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F10 = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F10 = WTIP2(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(J)= TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J)= TJM2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIM2
					IF (IST.EQ.-1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1

					XX(J) = TJM2
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (9.0D0*F00+3.0D0*F02-12.0D0*F01
     &				            +16.0D0*F11-4.0D0*F12-12.0D0*F10
     &				            +3.0D0*F20+F22-4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'BD-BD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I- below the lower bound, X_J+, X_J- inside the bounds.
			ELSE IF (.NOT.ILB .AND. JUB .AND. JLB) THEN
C  X_I++ below the upper bound.
				IF (TIP2.LT.XU(I)) THEN



					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(J)
					ELSE
						XX(J) = TJM
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(J)
					ELSE
						XX(J) = TJP
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(I) = TIP
					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIP2
					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (4.0D0*(F11P-F11M) - 3.0D0*(F01P-F01M)
     &				           -(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'FD-CD'
				ELSE
C  X_I is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+ above the upper bound, X_J+, X_J- inside the bounds.
			ELSE IF (.NOT.IUB .AND. JUB .AND. JLB) THEN
C  X_I-- above the lower bound.
				IF (TIM2.GT.XL(I)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(J)
					ELSE
						XX(J) = TJM
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(J)
					ELSE
						XX(J) = TJP
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIM2
					XX(J) = TJM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (-4.0D0*(F11P-F11M) + 3.0D0*(F01P-F01M)
     &				           +(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'BD-CD'
				ELSE
C  X_I is confined in a very small intrval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+, X_I- inside the bounds, X_J- below the lower bound.
			ELSE IF (.NOT.JLB .AND. IUB .AND. ILB) THEN
C  X_J++ below the upper bound.
				IF (TJP2.LT.XU(J)) THEN




					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(I)
					ELSE
						XX(I) = TIM
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(I)
					ELSE
						XX(I) = TIP
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					XX(J) = TJP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJP2
					XX(I) = TIM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (4.0D0*(F11P-F11M) - 3.0D0*(F01P-F01M)
     &				           -(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'CD-FD'
				ELSE
C  X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+, X_I- inside the bounds, X_J+ above the upper bound.
			ELSE IF (.NOT.JUB .AND. IUB .AND. ILB) THEN
C  X_J-- above the lower bound.
				IF (TJM2.GT.XL(J)) THEN




					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(I)
					ELSE
						XX(I) = TIM
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(I)
					ELSE
						XX(I) = TIP
						CALL SPAWNTASK(F, XX, N, FV(FC))
						FC = FC + 1

						NOC = NOC + 1


					END IF
C

					XX(J) = TJM
					XX(I) = TIM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(J) = TJM2
					XX(I) = TIM
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C

					XX(I) = TIP
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
C



					NOC = NOC + 4
					HES(I,J) = (-4.0D0*(F11P-F11M) + 3.0D0*(F01P-F01M)
     &				           +(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'CD-BD'
				ELSE
C  X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
			ELSE
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
	PR2 = IPRINT.GE.2
	PR3 = IPRINT.EQ.3

	FC = 1

	F00 = FV(FC)
	FC = FC + 1

C
C
C  Calculate the diagonal elements.
	IF (PR3) WRITE (*,30)

	DO 201,I=1,N
		ISTIP(I) = 0
		TI = XX(I)
		IF (UH(I).EQ.0D0) THEN
			HI = CON*MAX(ABS(TI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		TIP = TI+HI
		TIM = TI-HI
		TIP2 = TI+2.0D0*HI
		TIM2 = TI-2.0D0*HI
		TIP3 = TI+3.0D0*HI
		TIM3 = TI-3.0D0*HI
C  ---------------------------------------------------------------------
C  X_I+, X_I- are inside the bounds.
		IF (TIM.GT.XL(I) .AND. TIP.LT.XU(I)) THEN



			XX(I) = TIP
			FP = FV(FC)
			FC = FC + 1

C
			XX(I) = TIM
			FM = FV(FC)
			FC = FC + 1

C
			ISTIP(I) = 2
			WTIP(I) = FP
			WTIP2(I) = FM
C
c			NOC = NOC + 2
			HES(I,I) = (FP+FM-2.0D0*F00)/HI**2
			FORM = 'CD-CD'
C  ---------------------------------------------------------------------
C  X_I+++ below the upper bound.
		ELSE IF (TIP3.LT.XU(I)) THEN



			XX(I) = TIP
			F1 = FV(FC)
			FC = FC + 1
C
			XX(I)= TIP2
			F2 = FV(FC)
			FC = FC + 1
C
			XX(I)= TIP3
			F3 = FV(FC)
			FC = FC + 1

C
			ISTIP(I) = 1
			WTIP(I) = F1
			WTIP2(I) = F2
C
c			NOC = NOC + 3
			HES(I,I) = (2.0D0*F00-5.0D0*F1+4.0D0*F2-F3)/HI**2
			FORM = 'FD-FD'
C  ---------------------------------------------------------------------
C  X_I--- above the lower bound.
		ELSE IF (TIM3.GT.XL(I)) THEN



			XX(I) = TIM
			F1 = FV(FC)
			FC = FC + 1
C
			XX(I)= TIM2
			F2 = FV(FC)
			FC = FC + 1
C
			XX(I)= TIM3
			F3 = FV(FC)
			FC = FC + 1

C
			ISTIP(I) = -1
			WTIP(I) = F1
			WTIP2(I) = F2
C
c			NOC = NOC + 3
			HES(I,I) = (2.0D0*F00-5.0D0*F1+4.0D0*F2-F3)/HI**2
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
201	CONTINUE
C
C  Calculate off-diagonal elements.



	DO 11 J=1,N-1
		TJ = XX(J)
		IF (UH(J).EQ.0D0) THEN
			HJ = CON*MAX(ABS(TJ),ONE)
		ELSE
			HJ = UH(J)
		ENDIF
		TJP = TJ+HJ
		TJM = TJ-HJ
		TJP2 = TJ+2.0D0*HJ
		TJM2 = TJ-2.0D0*HJ
		JLB = TJM.GT.XL(J)
		JUB = TJP.LT.XU(J)
		DO 21,I=J+1,N
			TI = XX(I)
			IF (UH(I).EQ.0D0) THEN
				HI = CON*MAX(ABS(TI),ONE)
			ELSE
				HI = UH(I)
			ENDIF
			TIP = TI+HI
			TIM = TI-HI
			TIP2 = TI+2.0D0*HI
			TIM2 = TI-2.0D0*HI
			ILB = TIM.GT.XL(I)
			IUB = TIP.LT.XU(I)
C  ---------------------------------------------------------------------
C  X_I+ inside the bounds, X_J+ inside the bounds.
			IF (ILB .AND. JLB .AND. IUB .AND. JUB) THEN



				XX(I) = TIP
				XX(J) = TJP
				FPP = FV(FC)
				FC = FC + 1


				XX(J) = TJM
				FPM = FV(FC)
				FC = FC + 1

				XX(I) = TIM
				XX(J) = TJP
				FMP = FV(FC)
				FC = FC + 1


				XX(J) = TJM
				FMM = FV(FC)
				FC = FC + 1


c				NOC = NOC + 4
				HES(I,J) = (FPP-FPM+FMM-FMP)/(4.0D0*HI*HJ)
				FORM = 'CD-CD'
C  ---------------------------------------------------------------------
C  X_I- below the lower bound, X_J- below the lower bound.
			ELSE IF (.NOT.ILB .AND. .NOT.JLB) THEN
C  X_I++ below the upper bound, X_J++ below the upper bound.
				IF (TIP2.LT.XU(I) .AND. TJP2.LT.XU(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01 = WTIP(J)
					ELSE
						XX(J) = TJP
						F01 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1

					END IF
C

					IF (IST.EQ.1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJP2
						F02 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1

					END IF
C

					XX(I) = TIP
					IST = ISTIP(I)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F10 = WTIP(I)
					ELSE
						XX(J) = TJ
						F10 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1

					END IF
C

					XX(J)= TJP
					F11 = FV(FC)
					FC = FC + 1

					XX(J)= TJP2
					F12 = FV(FC)
					FC = FC + 1

					XX(I) = TIP2
					IF (IST.EQ.1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						F20 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1

					END IF
C

					XX(J) = TJP
					F21 = FV(FC)
					FC = FC + 1

					XX(J) = TJP2
					F22  = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (9.0D0*F00+3.0D0*F02-12.0D0*F01
     &				            +16.0D0*F11-4.0D0*F12-12.0D0*F10
     &				            +3.0D0*F20+F22-4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'FD-FD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I- below the lower bound, X_J+ above the upper bound.
			ELSE IF (.NOT.ILB .AND. .NOT.JUB) THEN
C  X_I++ below the upper bound, X_J-- above the lower bound.
				IF (TIP2.LT.XU(I) .AND. TJM2.GT.XL(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01 = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01 = WTIP2(J)
					ELSE
						XX(J) = TJM
						F01  = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.-1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJM2
						F02  = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(I) = TIP
					IST = ISTIP(I)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F10 = WTIP(I)
					ELSE
						XX(J) = TJ
						F10 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(J)= TJM
					F11 = FV(FC)
					FC = FC + 1

					XX(J)= TJM2
					F12 = FV(FC)
					FC = FC + 1

					XX(I) = TIP2
					IF (IST.EQ.1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						F20 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(J) = TJM
					F21 = FV(FC)
					FC = FC + 1

					XX(J) = TJM2
					F22  = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (-9.0D0*F00-3.0D0*F02+12.0D0*F01
     &				            -16.0D0*F11+4.0D0*F12+12.0D0*F10
     &				            -3.0D0*F20-F22+4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'FD-BD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+ above the upper bound, X_J- below the lower bound.
			ELSE IF (.NOT.IUB .AND. .NOT.JLB) THEN
C  X_I-- above the lower bound, X_I++ below the upper bound.
				IF (TIM2.GT.XL(I) .AND. TJP2.LT.XU(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01 = WTIP(J)
					ELSE
						XX(J) = TJP
						F01  = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJP2
						F02  = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F10 = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F10 = WTIP2(I)
					ELSE
						XX(J) = TJ
						F10 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(J)= TJP
					F11 = FV(FC)
					FC = FC + 1

					XX(J)= TJP2
					F12 = FV(FC)
					FC = FC + 1

					XX(I) = TIM2
					IF (IST.EQ.-1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						F20 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(J) = TJP
					F21  = FV(FC)
					FC = FC + 1

					XX(J) = TJP2
					F22  = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (-9.0D0*F00-3.0D0*F02+12.0D0*F01
     &				            -16.0D0*F11+4.0D0*F12+12.0D0*F10
     &				            -3.0D0*F20-F22+4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'BD-FD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+ above the upper bound, X_J+ above the upper bound.
			ELSE IF (.NOT.IUB .AND. .NOT.JUB) THEN
C  X_I-- above the lower bound, X_J-- above the lower bound.
				IF (TIM2.GT.XL(I) .AND. TJM2.GT.XL(J)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01 = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01 = WTIP2(J)
					ELSE
						XX(J) = TJM
						F01  = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.-1) THEN
						F02 = WTIP2(J)
					ELSE
						XX(J) = TJM2
						F02  = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F10 = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F10 = WTIP2(I)
					ELSE
						XX(J) = TJ
						F10 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(J)= TJM
					F11 = FV(FC)
					FC = FC + 1

					XX(J)= TJM2
					F12 = FV(FC)
					FC = FC + 1

					XX(I) = TIM2
					IF (IST.EQ.-1) THEN
						F20 = WTIP2(I)
					ELSE
						XX(J) = TJ
						F20 = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(J) = TJM
					F21  = FV(FC)
					FC = FC + 1
					XX(J) = TJM2
					F22  = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (9.0D0*F00+3.0D0*F02-12.0D0*F01
     &				            +16.0D0*F11-4.0D0*F12-12.0D0*F10
     &				            +3.0D0*F20+F22-4.0D0*F21)/
     &				           (4.0D0*HI*HJ)
					FORM = 'BD-BD'
				ELSE
C  X_I or X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I- below the lower bound, X_J+, X_J- inside the bounds.
			ELSE IF (.NOT.ILB .AND. JUB .AND. JLB) THEN
C  X_I++ below the upper bound.
				IF (TIP2.LT.XU(I)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(J)
					ELSE
						XX(J) = TJM
						F01M = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(J)
					ELSE
						XX(J) = TJP
						F01P = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(I) = TIP
					XX(J) = TJM
					F11M = FV(FC)
					FC = FC + 1

					XX(J) = TJP
					F11P = FV(FC)
					FC = FC + 1

					XX(I) = TIP2
					XX(J) = TJM
					F21M = FV(FC)
					FC = FC + 1

					XX(J) = TJP
					F21P = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (4.0D0*(F11P-F11M) - 3.0D0*(F01P-F01M)
     &				           -(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'FD-CD'
				ELSE
C  X_I is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+ above the upper bound, X_J+, X_J- inside the bounds.
			ELSE IF (.NOT.IUB .AND. JUB .AND. JLB) THEN
C  X_I-- above the lower bound.
				IF (TIM2.GT.XL(I)) THEN




					IST = ISTIP(J)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(J)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(J)
					ELSE
						XX(J) = TJM
						F01M = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(J)
					ELSE
						XX(J) = TJP
						F01P = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					XX(J) = TJM
					F11M =  FV(FC)
					FC = FC + 1

					XX(J) = TJP
					F11P = FV(FC)
					FC = FC + 1

					XX(I) = TIM2
					XX(J) = TJM
					F21M = FV(FC)
					FC = FC + 1

					XX(J) = TJP
					F21P = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (-4.0D0*(F11P-F11M) + 3.0D0*(F01P-F01M)
     &				           +(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'BD-CD'
				ELSE
C  X_I is confined in a very small intrval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+, X_I- inside the bounds, X_J- below the lower bound.
			ELSE IF (.NOT.JLB .AND. IUB .AND. ILB) THEN
C  X_J++ below the upper bound.
				IF (TJP2.LT.XU(J)) THEN




					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(I)
					ELSE
						XX(I) = TIM
						F01M = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(I)
					ELSE
						XX(I) = TIP
						F01P = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					XX(J) = TJP
					F11M = FV(FC)
					FC = FC + 1

					XX(I) = TIP
					F11P = FV(FC)
					FC = FC + 1

					XX(J) = TJP2
					XX(I) = TIM
					F21M = FV(FC)
					FC = FC + 1

					XX(I) = TIP
					F21P = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (4.0D0*(F11P-F11M) - 3.0D0*(F01P-F01M)
     &				           -(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'CD-FD'
				ELSE
C  X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
C  ---------------------------------------------------------------------
C  X_I+, X_I- inside the bounds, X_J+ above the upper bound.
			ELSE IF (.NOT.JUB .AND. IUB .AND. ILB) THEN
C  X_J-- above the lower bound.
				IF (TJM2.GT.XL(J)) THEN




					IST = ISTIP(I)
					IF (IST.EQ.-1) THEN
						F01M = WTIP(I)
					ELSE IF (IST.EQ.2) THEN
						F01M = WTIP2(I)
					ELSE
						XX(I) = TIM
						F01M = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					IF (IST.EQ.1 .OR. IST.EQ.2) THEN
						F01P = WTIP(I)
					ELSE
						XX(I) = TIP
						F01P = FV(FC)
						FC = FC + 1
c						NOC = NOC + 1


					END IF
C

					XX(I) = TIM
					XX(J) = TJM
					F11M = FV(FC)
					FC = FC + 1

					XX(I) = TIP
					F11P = FV(FC)
					FC = FC + 1

					XX(J) = TJM2
					XX(I) = TIM
					F21M = FV(FC)
					FC = FC + 1

					XX(I) = TIP
					F21P = FV(FC)
					FC = FC + 1



c					NOC = NOC + 4
					HES(I,J) = (-4.0D0*(F11P-F11M) + 3.0D0*(F01P-F01M)
     &				           +(F21P-F21M)) / (4.0D0*HI*HJ)
					FORM = 'CD-BD'
				ELSE
C  X_J is confined in a very small interval.
					HES(I,J) = 0.0D0
					FORM = '     '
					IF (PR3) WRITE (*,60) I, J
				END IF
			ELSE
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
CCCCCCCC

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
