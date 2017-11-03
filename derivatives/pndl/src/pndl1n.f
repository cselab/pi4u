C  ---------------------------------------------------------------------
	SUBROUTINE PNDL1N ( F, X, N, XL, XU, UH, FEPS, IPRINT, G, NOC )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the gradient vector using O(h**4)
C    formulae.
C    Selection of the appropriate formula is based on the proximity
C    of the variable with its bounds.
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
C                             steps and the resulting gradient
C                             vector).
C
C  Output arguments:
C    G          Array containing the resulting gradient vector.
C    NOC        Number of calls to the function F.
C
C  Notes:
C    The main differentiation formula uses central differences, however
C    when one or more variables are near the bounds, forward or backward
C    difference formulae are used. The later require evaluation of the
C    function at the given point X. To avoid uneccessary evaluations
C    when more than one variable are near the bounds, a simple caching
C    mechanism is implemented. Variable NOF0 indicates whether the
C    function value F0=F(X) has been evaluated.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N), XL(N), XU(N), UH(N), G(N)
C
	DIMENSION XX(N)
	DIMENSION FV(4*N)
	INTEGER FC
	INTEGER K
C                                                  
        EXTERNAL F, SPAWNTASK
C
	LOGICAL NOF0, PR2, PR3
	CHARACTER*2 FORM
C
	DATA ONE   / 1.0D0 /
	DATA TWO   / 2.0D0 /
	DATA FOUR  / 4.0D0 /
	DATA EIGHT / 8.0D0 /
C
	DO K=1,N
	XX(K) = X(K)
	END DO
C
	CON = FEPS**0.2D0
	
C	!! FIRST PASS !!
	FC = 1
	NOC = 0
	NOF0 = .TRUE.
	PR2 = .FALSE.
	PR3 = .FALSE.
C
	IF (PR3) WRITE (*,40)
	DO 10,I=1,N
		XI = X(I)
		IF (UH(I).EQ.0D0) THEN
			HI = CON*MAX(ABS(XI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		XP = XI+HI
		XM = XI-HI
		X2P = XI+TWO*HI
		X2M = XI-TWO*HI
		GI = 0.0D0
		IF (X2M.GT.XL(I) .AND. X2P.LT.XU(I)) THEN
			XX(I) = XP
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1
			
			XX(I) = XM
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

			GI = GI + 8.0D0*(FP-FM)

			XX(I) = X2P
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

			XX(I) = X2M
			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

			GI = GI - (F2P-F2M)
			GI = GI / (12.0D0*HI)
			NOC = NOC+4
			FORM = 'CD'
		ELSE IF (X2M.LT.XL(I)) THEN
			X4P = XI+FOUR*HI
			X8P = XI+EIGHT*HI
			IF (X8P.LT.XU(I)) THEN
				IF (NOF0) THEN
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
			
					NOC = NOC+1
					NOF0 = .FALSE.
				END IF
C
				XX(I) = XP
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI + 512.0D0*(FP-F0)
C
				XX(I) = X2P
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI-224.0D0*(F2P-F0)
C
				XX(I) = X4P
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI+28.0D0*(F4P-F0)
C
				XX(I) = X8P
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI-(F8P-F0)
C
				GI = GI / (168.0D0*HI)
				NOC = NOC+4
				FORM = 'FD'
			ELSE
				GI = 0.0D0
				IF (PR2) WRITE (*,20) I
			END IF
		ELSE IF (X2P.GT.XU(I)) THEN
			X4M = XI-FOUR*HI
			X8M = XI-EIGHT*HI
			IF (X8M.GT.XL(I)) THEN
				IF (NOF0) THEN
					CALL SPAWNTASK(F, XX, N, FV(FC))
					FC = FC + 1
					NOC = NOC+1
					NOF0 = .FALSE.
				END IF
C
				XX(I) = XM
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI-512.0D0*(FM-F0)
C
				XX(I) = X2M
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI+224.0D0*(F2M-F0)
C
				XX(I) = X4M
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI-28.0D0*(F4M-F0)
C
				XX(I) = X8M
				CALL SPAWNTASK(F, XX, N, FV(FC))
				FC = FC + 1
				GI = GI+(F8M-F0)
C
				GI = GI / (168.0D0*HI)
				NOC = NOC+4
				FORM = 'BD'
			ELSE
				GI = 0.0D0
				IF (PR2) WRITE (*,20) I
			END IF
		ELSE
			GI = 0.0D0
			IF (PR2) WRITE (*,20) I
		END IF
		G(I) = GI
		XX(I) = XI
		IF (PR3) WRITE (*,30) I, FORM, XI, HI, GI
10	CONTINUE

c	CALL torc_enable_stealing()
	CALL torc_waitall()
c	CALL torc_disable_stealing()

C	!! SECOND PASS !!
	FC = 1
c	NOC = 0
	NOF0 = .TRUE.
	PR2 = IPRINT.GE.2
	PR3 = IPRINT.EQ.3
C
	IF (PR3) WRITE (*,40)
	DO 11,I=1,N
		XI = XX(I)
		IF (UH(I).EQ.0D0) THEN
			HI = CON*MAX(ABS(XI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		XP = XI+HI
		XM = XI-HI
		X2P = XI+TWO*HI
		X2M = XI-TWO*HI
		GI = 0.0D0
		IF (X2M.GT.XL(I) .AND. X2P.LT.XU(I)) THEN
			XX(I) = XP
			FP = FV(FC)
			FC = FC + 1
			XX(I) = XM
			FM = FV(FC)
			FC = FC + 1
			GI = GI + 8.0D0*(FP-FM)
			XX(I) = X2P
			F2P = FV(FC)
			FC = FC + 1
			XX(I) = X2M
			F2M = FV(FC)
			FC = FC + 1
			GI = GI - (F2P-F2M)
			GI = GI / (12.0D0*HI)
c			NOC = NOC+4
			FORM = 'CD'
		ELSE IF (X2M.LT.XL(I)) THEN
			X4P = XI+FOUR*HI
			X8P = XI+EIGHT*HI
			IF (X8P.LT.XU(I)) THEN
				IF (NOF0) THEN
					F0 = FV(FC)
					FC = FC + 1
c					NOC = NOC+1
					NOF0 = .FALSE.
				END IF
C
				XX(I) = XP
				FP = FV(FC)
				FC = FC + 1
				GI = GI + 512.0D0*(FP-F0)
C
				XX(I) = X2P
				F2P = FV(FC)
				FC = FC + 1
				GI = GI-224.0D0*(F2P-F0)
C
				XX(I) = X4P
				F4P = FV(FC)
				FC = FC + 1
				GI = GI+28.0D0*(F4P-F0)
C
				XX(I) = X8P
				F8P = FV(FC)
				FC = FC + 1
				GI = GI-(F8P-F0)
C
				GI = GI / (168.0D0*HI)
c				NOC = NOC+4
				FORM = 'FD'
			ELSE
				GI = 0.0D0
				IF (PR2) WRITE (*,20) I
			END IF
		ELSE IF (X2P.GT.XU(I)) THEN
			X4M = XI-FOUR*HI
			X8M = XI-EIGHT*HI
			IF (X8M.GT.XL(I)) THEN
				IF (NOF0) THEN
					F0 = FV(FC)
					FC = FC + 1
c					NOC = NOC+1
					NOF0 = .FALSE.
				END IF
C
				XX(I) = XM
				FM = FV(FC)
				FC = FC + 1
				GI = GI-512.0D0*(FM-F0)
C
				XX(I) = X2M
				F2M = FV(FC)
				FC = FC + 1
				GI = GI+224.0D0*(F2M-F0)
C
				XX(I) = X4M
				F4M = FV(FC)
				FC = FC + 1
				GI = GI-28.0D0*(F4M-F0)
C
				XX(I) = X8M
				F8M = FV(FC)
				FC = FC + 1
				GI = GI+(F8M-F0)
C
				GI = GI / (168.0D0*HI)
c				NOC = NOC+4
				FORM = 'BD'
			ELSE
				GI = 0.0D0
				IF (PR2) WRITE (*,20) I
			END IF
		ELSE
			GI = 0.0D0
			IF (PR2) WRITE (*,20) I
		END IF
		G(I) = GI
		XX(I) = XI
		IF (PR3) WRITE (*,30) I, FORM, XI, HI, GI
11	CONTINUE
C
20	FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
30	FORMAT (' PNDL:',I6,3X,A,3X,1PG21.14,1X,1PG21.14,1X,1PG21.14)
40	FORMAT (/' PNDL:',' Index',1X,'Formula',7X,'X_i',19X,'Step_i',17X,
     &        'G_i')
	END 
