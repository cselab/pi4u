C  ---------------------------------------------------------------------
      SUBROUTINE PNDL2GF ( GRD, X, N, XL, XU, UH, FEPS, IPRINT, HES,
     &                    LD, IGC )
C     &                    LD, GW, G, IWFIX, IGC )
C  ---------------------------------------------------------------------
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the Hessian matrix using O(h_i)+O(h_j)
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
C    GW         Real work space of length N.
C    G          Real work space of length N.
C    IWFIX      Integer work space of length N.
C
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N), XL(N), XU(N), UH(N), HES(LD,N)
C
      DIMENSION GW(N), G(N), IWFIX(N)
C
      DIMENSION IFB(N)
      DIMENSION HH(N)
      DIMENSION XX(N)
      INTEGER K
C
      EXTERNAL GRD, SPAWNTASKG
C
      LOGICAL PR2, PR3
      CHARACTER*2 FORM
C
      DATA ONE / 1.0D0 /
C
      DO K=1,N
      XX(K) = X(K)
      END DO
C
      CON = SQRT(FEPS)
C
C  Analytic gradient call counter.
      IGC = 0
C
      PR2 = IPRINT.GE.2
      PR3 = IPRINT.EQ.3
C
C  Get the gradient at the given point.
      CALL SPAWNTASKG(GRD, XX, N, G)
      IGC = IGC + 1
C
C  Loop over all variables (H-columns)
      IF (PR3) WRITE (*,30)
      DO 100,I=1,N
C
          IWFIX(I) = 1
          TI = XX(I)
          IF (UH(I).EQ.0.0D0) THEN
              HI = CON*MAX(ABS(TI),ONE)
          ELSE
              HI = UH(I)
          ENDIF
          HH(I) = HI
          TIP = TI + HI
          TIM = TI - HI
C
C  Basic case: the variable plus the differentiation step does not
C  exceed the upper bound. Use forward differences.
          IF (TIP.LT.XU(I)) THEN
              XX(I) = TIP
              CALL SPAWNTASKG(GRD, XX, N, HES(1,I))
              IGC = IGC + 1
              IFB(I) = 1
              FORM = 'FD'

C
C  The variable plus the differentiation step exceeds the upper bound.
C  Check if we can use backward differences.
          ELSE IF (TIM.GT.XL(I)) THEN
              XX(I) = TIM
              CALL SPAWNTASKG(GRD, XX, N, HES(1,I))
              IGC = IGC + 1
              IFB(I) = -1
              FORM = 'BD'
          ELSE
C
C  The variable is confined in a very small interval. (An "almost fixed"
C  variable). Numerical differentiation is not possible.
              IWFIX(I) = 0
              IFB(I) = 0
              FORM = '  '
              IF (PR2) WRITE (*,70) I
          END IF
          XX(I) = TI
          IF (PR3) WRITE (*,60) I, FORM, TI, HI
100    CONTINUE
C

c    CALL torc_enable_stealing()
      CALL torc_waitall()
c    CALL torc_disable_stealing()

C    !! ~SECOND PASS !!
C  Perform division using H()
          DO 2,J=1,N
              IF (IFB(J).EQ.1) THEN
                  DO 21,I=1,N
                      HES(I,J) = (HES(I,J)-G(I))/HH(J)
21                CONTINUE
              ELSE IF (IFB(J).EQ.-1) THEN
                  DO 22,I=1,N
                      HES(I,J) = (G(I)-HES(I,J))/HH(J)
22                CONTINUE
              END IF
2        CONTINUE



C  Make sure the Hessian is symmetric. (The result is stored in
C  the lower triangular part).
      DO 40,J=1,N
          IWJ = IWFIX(J)
          DO 50,I=J,N
              HES(I,J) = IWJ*IWFIX(I)*(HES(I,J)+HES(J,I))/2.0D0
50        CONTINUE
40    CONTINUE
C
30    FORMAT (/' PNDL:',' Index',1X,'Formula',7X,'X_i',19X,'Step_i')
60    FORMAT (' PNDL:',I6,3X,A,3X,1PG21.14,1X,1PG21.14)
70    FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
      END

C  ---------------------------------------------------------------------
      SUBROUTINE SPAWNTASKG (GRD, XX, N, RES)
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL GRD
      DIMENSION XX(N), RES(N)
      INCLUDE 'torcf.h'

      CALL torc_taskf(GRD, 0, 3,
     &                N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &                1, MPI_INTEGER,          CALL_BY_VAL,
     &                N, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &                XX, N, RES)

      RETURN
      END

