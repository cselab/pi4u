C  ---------------------------------------------------------------------
	SUBROUTINE PNDL1F ( F, X, N, XL, XU, UH, FEPS, IPRINT, G, NOC )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the gradient vector using O(h)
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
	LOGICAL PR2, PR3
	CHARACTER*2 FORM
C
	DATA ONE  / 1.0D0 /
C
	INCLUDE 'torcf.h'
C
	CON = SQRT(FEPS)

C
	DO K=1,N
	XX(K) = X(K)
	END DO
	
C	!! FIRST PASS !!
	PR2 = .FALSE.
	PR3 = .FALSE.

	FC = 1
	CALL SPAWNTASK(F, XX, N, FV(FC))
    
	FC = FC + 1
	
	NOC = 1
C
	IF (PR3) WRITE (*,40)
	DO 10,I=1,N
		XI = XX(I)
		IF (UH(I).EQ.0D0) THEN
			HI = CON*MAX(ABS(XI),ONE)
		ELSE
			HI = UH(I)
		ENDIF
		XP = XI+HI
		XM = XI-HI
		IF (XP.LT.XU(I)) THEN
			XX(I) = XP

			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

			G(I) = (FP-F0)/HI
			NOC = NOC+1
			FORM = 'FD'
		ELSE IF (XM.GT.XL(I)) THEN
			XX(I) = XM

			CALL SPAWNTASK(F, XX, N, FV(FC))
			FC = FC + 1

			G(I) = (F0-FM)/HI
			NOC = NOC+1
			FORM = 'BD'
		ELSE
			G(I) = 0.0D0
			FORM = '  '
			IF (PR2) WRITE (*,20) I
		END IF
		XX(I) = XI
		IF (PR3) WRITE (*,30) I, FORM, XI, HI, G(I)
10	CONTINUE

c	CALL torc_enable_stealing()
	CALL torc_waitall()
c	CALL torc_disable_stealing()

C
C	!! SECOND PASS !!
	PR2 = IPRINT.GE.2
	PR3 = IPRINT.EQ.3

	FC = 1
	F0 = FV(FC)
	FC = FC + 1
c	NOC = 1
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
		IF (XP.LT.XU(I)) THEN
			XX(I) = XP
			FP = FV(FC)
			FC = FC + 1

			G(I) = (FP-F0)/HI
c			NOC = NOC+1
			FORM = 'FD'
		ELSE IF (XM.GT.XL(I)) THEN
			XX(I) = XM
			FM =  FV(FC)
			FC = FC + 1
			G(I) = (F0-FM)/HI
c			NOC = NOC+1
			FORM = 'BD'
		ELSE
			G(I) = 0.0D0
			FORM = '  '
			IF (PR2) WRITE (*,20) I
		END IF
		XX(I) = XI
		IF (PR3) WRITE (*,30) I, FORM, XI, HI, G(I)
11	CONTINUE
C
20	FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
30	FORMAT (' PNDL:',I6,3X,A,3X,1PG21.14,1X,1PG21.14,1X,1PG21.14)
40	FORMAT (/' PNDL:',' Index',1X,'Formula',7X,'X_i',19X,'Step_i',17X,
     &        'G_i')
	END



C  ---------------------------------------------------------------------
	SUBROUTINE SPAWNTASK (F, XX, N, RES)
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	EXTERNAL F, TASKFUN
	DIMENSION XX(N)
	include 'torcf.h' 
    
	call torc_task(TASKFUN, 0, 4,
     &	    1, MPI_INTEGER, CALL_BY_VAD,
     &	    N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &	    1, MPI_INTEGER, CALL_BY_VAL,
     &	    1, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &	    F, XX, N, RES)

c	CALL TASKFUN(F, XX, N, RES)

	RETURN
	END 

    
C  ---------------------------------------------------------------------
	SUBROUTINE TASKFUN (F1, XX, N, RES)
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	EXTERNAL F1
	DIMENSION XX(N)

	RES = F1(XX, N)
	RETURN
	END 

