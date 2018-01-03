C  ---------------------------------------------------------------------
      SUBROUTINE PNDLJF ( RSD, X, N, M, XL, XU, UH, FEPS, IPRINT, FJ,
     &                   LD, NOC)
c     &                  LD, NOC, F0 )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the Jacobian matrix (FJ) using O(h)
C    formulae.
C    Selection of the appropriate formula is based on the proximity
C    of the variable with its bounds.
C
C  Input arguments:
C    RSD        A subroutine that returns the residuals (F).
C               Must be declared as:
C                 SUBROUTINE RSD ( X, N, M, F )
C                 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                 DIMENSION X(N), F(M)
C    X          Array containing the function variables.
C    N          The number of variables.
C    M          The number of squared terms.
C    XL         Array containing lower bounds on the variables.
C    XU         Array containing upper bounds on the variables.
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
C    LD         Leading dimension of matrix FJ.
C
C  Output arguments:
C    FJ         The Jacobian matrix. FJ(i,j) = dRi/dXj
C    NOC        Number of calls to subroutine RSD.
C
C  Work spaces:
C    F0         Real work space of length M.
C
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N), XL(N), XU(N), UH(N), FJ(LD,N)
C
      DIMENSION F0(M)
C
      DIMENSION XX(N)
      DIMENSION IFB(N)
      DIMENSION HH(N)
      INTEGER K
C
      EXTERNAL RSD, SPAWNTASKJ
C
      LOGICAL PR2, PR3
      CHARACTER*2 FORM
C
      DATA ONE  / 1.0D0 /
C
      DO K=1,N
      XX(K) = X(K)
      END DO
C
      CON = SQRT(FEPS)

      CALL SPAWNTASKJ(RSD, XX, N, M, F0)

      NOC = 1
      PR2 = IPRINT.GE.2
      PR3 = IPRINT.EQ.3
C
      IF (PR3) WRITE (*,50)


      DO 10,I=1,N
          XI = XX(I)
          IF (UH(I).EQ.0.0D0) THEN
              HI = CON*MAX(ABS(XI),ONE)
          ELSE
              HI = UH(I)
          ENDIF
          HH(I) = HI
          XP = XI+HI
          XM = XI-HI
C  ---------------------------------------------------------------------
C  X_I+ below the upper bound.
          IF (XP.LT.XU(I)) THEN
              XX(I) = XP
              CALL SPAWNTASKJ(RSD, XX, N, M, FJ(1,I))
              NOC = NOC+1
c            DO 20,K=1,M
c20                FJ(K,I) = (FJ(K,I)-F0(K))/HI
              IFB(I) = 1
              FORM = 'FD'
C  ---------------------------------------------------------------------
C  X_I- above the lower bound.
          ELSE IF (XM.GT.XL(I)) THEN
              XX(I) = XM
              CALL SPAWNTASKJ(RSD, XX, N, M, FJ(1,I))
              NOC = NOC+1
c            DO 30,K=1,M
c30                FJ(K,I) = (F0(K)-FJ(K,I))/HI
              IFB(I) = -1
              FORM = 'BD'
C  ---------------------------------------------------------------------
C  X_I is confined in a very small interval.
          ELSE
              DO 40,K=1,M
40                FJ(K,I) = 0.0D0
              FORM = ' '
              IFB(I) = 0
              IF (PR2) WRITE (*,70)
          END IF
          XX(I) = XI
          IF (PR3) WRITE (*,60) I, FORM, XI, HI
10    CONTINUE
C
c    CALL torc_enable_stealing()
      CALL torc_waitall()
c    CALL torc_disable_stealing()

C  Perform division using H()
      DO 2,I=1,N
          IF (IFB(I).EQ.1) THEN
              DO 3,K=1,M
                   FJ(K,I) = (FJ(K,I)-F0(K))/HH(I)
3            CONTINUE
          ELSE IF (IFB(I).EQ.-1) THEN
              DO 4,K=1,M
                   FJ(K,I) = (F0(K)-FJ(K,I))/HH(I)
4            CONTINUE
          END IF
2     CONTINUE


C
50    FORMAT (/' PNDL:',' Index',1X,'Formula',7X,'X_i',19X,'Step_i')
60    FORMAT (' PNDL:',I6,3X,A,3X,1PG21.14,1X,1PG21.14)
70    FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
      END

C  ---------------------------------------------------------------------
      SUBROUTINE SPAWNTASKJ (RSD, XX, N, M, RES)
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RSD
      DIMENSION XX(N), RES(N)
      INCLUDE 'torcf.h'

      CALL torc_taskf(RSD, 0, 4,
     &                N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &                1, MPI_INTEGER,          CALL_BY_VAL,
     &                1, MPI_INTEGER,          CALL_BY_VAL,
     &                M, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &                XX, N, M, RES)

      RETURN
      END SUBROUTINE SPAWNTASKJ

