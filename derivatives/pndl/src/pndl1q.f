C  ---------------------------------------------------------------------
      SUBROUTINE PNDL1Q ( F, X, N, XL, XU, UH, FEPS, IPRINT, G, NOC )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the gradient vector using O(h**2)
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
      DIMENSION FV(4*N)
      INTEGER FC
      DIMENSION XX(N)
C
          EXTERNAL F, SPAWNTASK
C
      LOGICAL NOF0, PR2, PR3
      CHARACTER*2 FORM
C
      DATA ONE   / 1.0D0 /
C
      DO K=1,N
      XX(K) = X(K)
      END DO
C
      CON = FEPS**(1.0D0/3.0D0)

C   !! FIRST PASS !!

      NOC = 0
      NOF0 = .TRUE.
      PR2 = .FALSE.
      PR3 = .FALSE.
      FC = 1
C
      IF (PR3) WRITE (*,40)

c    PRINT *, 'PNDL1Q'

      DO 10,I=1,N
          XI = XX(I)
          IF (UH(I).EQ.0D0) THEN
              HI = CON*MAX(ABS(XI),ONE)
          ELSE
              HI = UH(I)
          ENDIF
          XP = XI+HI
          XM = XI-HI
          IF (XM.GT.XL(I) .AND. XP.LT.XU(I)) THEN



              XX(I) = XP
                  CALL SPAWNTASK(F, XX, N, FV(FC))
              FC = FC + 1

              XX(I) = XM
                  CALL SPAWNTASK(F, XX, N, FV(FC))
              FC = FC + 1

              G(I) = (FP-FM)/(2.0D0*HI)
              NOC = NOC+2
              FORM = 'CD'
          ELSE IF (XM.LT.XL(I)) THEN
              XPP = XI+2.0D0*HI
              IF (XPP.LT.XU(I)) THEN

                  IF (NOF0) THEN
                          CALL SPAWNTASK(F, XX, N, FV(FC))
                      FC = FC + 1
                      NOC = NOC+1
                      NOF0 = .FALSE.
                  END IF




                  XX(I) = XP
                      CALL SPAWNTASK(F, XX, N, FV(FC))
                  FC = FC + 1

                  XX(I) = XPP
                      CALL SPAWNTASK(F, XX, N, FV(FC))
                  FC = FC + 1

                  G(I) = (4.0D0*FP-3.0D0*F0-FPP)/(2.0D0*HI)
                  NOC = NOC+2
                  FORM = 'FD'
              ELSE
                  G(I) = 0.0D0
                  FORM = '  '
                  IF (PR2) WRITE (*,20) I
              END IF
          ELSE IF (XP.GT.XU(I)) THEN
              XMM = XI-2.0D0*HI
              IF (XMM.GT.XL(I)) THEN

                  IF (NOF0) THEN
                          CALL SPAWNTASK(F, XX, N, FV(FC))
                      FC = FC + 1

                      NOC = NOC+1
                      NOF0 = .FALSE.
                  END IF




                  XX(I) = XM
                      CALL SPAWNTASK(F, XX, N, FV(FC))
                  FC = FC + 1

                  XX(I) = XMM
                      CALL SPAWNTASK(F, XX, N, FV(FC))
                  FC = FC + 1

                  G(I) = -(4.0D0*FM-3.0D0*F0-FMM)/(2.0D0*HI)
                  NOC = NOC+2
                  FORM = 'BD'
              ELSE
                  G(I) = 0.0D0
                  FORM = '  '
                  IF (PR2) WRITE (*,20) I
              END IF
          ELSE
              G(I) = 0.0D0
              FORM = '  '
              IF (PR2) WRITE (*,20) I
          END IF
          XX(I) = XI
          IF (PR3) WRITE (*,30) I, FORM, XI, HI, G(I)
10    CONTINUE

c    CALL torc_enable_stealing()
      CALL torc_waitall()
c    CALL torc_disable_stealing()

C   !! SECOND PASS !!

      FC = 1
c    NOC = 0
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
          IF (XM.GT.XL(I) .AND. XP.LT.XU(I)) THEN



              XX(I) = XP
              FP = FV(FC)
              FC = FC + 1

              XX(I) = XM
              FM = FV(FC)
              FC = FC + 1

              G(I) = (FP-FM)/(2.0D0*HI)
c            NOC = NOC+2
              FORM = 'CD'
          ELSE IF (XM.LT.XL(I)) THEN
              XPP = XI+2.0D0*HI
              IF (XPP.LT.XU(I)) THEN

                  IF (NOF0) THEN
                      F0 = FV(FC)
                      FC = FC + 1
c                    NOC = NOC+1
                      NOF0 = .FALSE.
                  END IF




                  XX(I) = XP
                  FP = FV(FC)
                  FC = FC + 1

                  XX(I) = XPP
                  FPP =  FV(FC)
                  FC = FC + 1

                  G(I) = (4.0D0*FP-3.0D0*F0-FPP)/(2.0D0*HI)
c                NOC = NOC+2
                  FORM = 'FD'
              ELSE
                  G(I) = 0.0D0
                  FORM = '  '
                  IF (PR2) WRITE (*,20) I
              END IF
          ELSE IF (XP.GT.XU(I)) THEN
              XMM = XI-2.0D0*HI
              IF (XMM.GT.XL(I)) THEN

                  IF (NOF0) THEN
                      F0 =  FV(FC)
                      FC = FC + 1
c                    NOC = NOC+1
                      NOF0 = .FALSE.
                  END IF




                  XX(I) = XM
                  FM = FV(FC)
                  FC = FC + 1

                  XX(I) = XMM
                  FMM =  FV(FC)
                  FC = FC + 1

                  G(I) = -(4.0D0*FM-3.0D0*F0-FMM)/(2.0D0*HI)
c                NOC = NOC+2
                  FORM = 'BD'
              ELSE
                  G(I) = 0.0D0
                  FORM = '  '
                  IF (PR2) WRITE (*,20) I
              END IF
          ELSE
              G(I) = 0.0D0
              FORM = '  '
              IF (PR2) WRITE (*,20) I
          END IF
          XX(I) = XI
          IF (PR3) WRITE (*,30) I, FORM, XI, HI, G(I)
11    CONTINUE

C
20    FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
30    FORMAT (' PNDL:',I6,3X,A,3X,1PG21.14,1X,1PG21.14,1X,1PG21.14)
40    FORMAT (/' PNDL:',' Index',1X,'Formula',7X,'X_i',19X,'Step_i',17X,
     &        'G_i')
      END
