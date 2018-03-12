C  ---------------------------------------------------------------------
      SUBROUTINE PNDLJQ ( RSD, X, N, M, XL, XU, UH, FEPS, IPRINT, FJ,
     &                   LD, NOC)
c     &                  LD, NOC, F0, F1 )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Computes approximations to the Jacobian matrix (FJ) using O(h**2)
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
C    F1         Real work space of length MxN.
C
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N), XL(N), XU(N), UH(N), FJ(LD,N)
C
      DIMENSION F0(M), F1(M,N)
C
      DIMENSION XX(N)
      DIMENSION IFB(N)
      DIMENSION HH(N)
      INTEGER K
C
      EXTERNAL RSD, SPAWNTASKJ
C
      LOGICAL NOF0, PR2, PR3
      CHARACTER*2 FORM
C
      DATA ONE   / 1.0D0 /
      DATA THREE / 3.0D0 /
      DATA FOUR  / 4.0D0 /
C
      DO K=1,N
      XX(K) = X(K)
      END DO
C
      CON = FEPS**(1.0D0/3.0D0)
C
      NOC = 0
      NOF0 = .TRUE.
      PR2 = IPRINT.GE.2
      PR3 = IPRINT.EQ.3
C
      IF (PR3) WRITE (*,90)

C    !! FIRST PASS !!
      DO 10,I=1,N
          XI = XX(I)
          IF (UH(I).EQ.0.0D0) THEN
              HI = CON*MAX(ABS(XI),ONE)
          ELSE
              HI = UH(I)
          ENDIF
          HH(I) = HI
          TWOHI = 2.0D0*HI
          XP = XI+HI
          XM = XI-HI
C  ---------------------------------------------------------------------
C  X_I-, X_I+ inside the bounds.
          IF (XM.GT.XL(I) .AND. XP.LT.XU(I)) THEN

              XX(I) = XP
              CALL SPAWNTASKJ(RSD,XX,N,M,FJ(1,I))

              XX(I) = XM
              CALL SPAWNTASKJ(RSD,XX,N,M,F1(1,I))

c            DO 20,K=1,M
c20                FJ(K,I) = (FJ(K,I)-F1(K,I))/TWOHI
              NOC = NOC+2
              IFB(I) = 0
              FORM = 'CD'
C  ---------------------------------------------------------------------
C  X_I- below the lower bound.
          ELSE IF (XM.LT.XL(I)) THEN
              XPP = XI+2.0D0*HI
C  X_I++ below the upper bound.
              IF (XPP.LT.XU(I)) THEN

                  IF (NOF0) THEN
                      CALL SPAWNTASKJ(RSD,XX,N,M,F0)
c                    do K=1, M
c                    F0(K) = -THREE * F0(K)
c                    end do
                      NOC = NOC+1
                      NOF0 = .FALSE.
                  END IF


                  XX(I) = XP
                  CALL SPAWNTASKJ(RSD,XX,N,M,FJ(1,I))
c                do K=1, M
c                    FJ(K, I) = FOUR * FJ(K, I)
c                end do

                  XX(I) = XPP
                  CALL SPAWNTASKJ(RSD,XX,N,M,F1(1,I))
c                do K=1, M
c                    F1(K,I) = -F1(K,I)
c                end do

c                DO 30,K=1,M
c30                    FJ(K,I) = (FJ(K,I)+F0(K)+F1(K,I))/TWOHI
                  NOC = NOC+2
                  IFB(I) = 1
                  FORM = 'FD'
C  X_I is confined in a very small interval.
              ELSE
                  DO 40,K=1,M
40                    FJ(K,I) = 0.0D0
                  IFB(I) = 2
                  FORM = '  '
                  IF (PR2) WRITE (*,100) I
              END IF
C  ---------------------------------------------------------------------
C  X_I+ above the upper bound.
          ELSE IF (XP.GT.XU(I)) THEN
              XMM = XI-2.0D0*HI
C  X_I-- above the lower bound.
              IF (XMM.GT.XL(I)) THEN

                  IF (NOF0) THEN
                      CALL SPAWNTASKJ(RSD,XX,N,M,F0)
c                    do K=1, M
c                    F0(K) = THREE * F0(K)
c                    end do

                      NOC = NOC+1
                      NOF0 = .FALSE.
                  END IF

                  XX(I) = XM
                  CALL SPAWNTASKJ(RSD,XX,N,M,FJ(1,I))
c                do K=1, M
c                    FJ(K, I) = -FOUR * FJ(K, I)
c                end do

                  XX(I) = XMM
                  CALL SPAWNTASKJ(RSD,XX,N,M,F1(1,I))
c                do K=1, M
c                    F1(K,I) = F1(K,I)
c                end do

c                DO 50,K=1,M
c50                    FJ(K,I) = (FJ(K,I)+F0(K)+F1(K,I))/TWOHI

                  NOC = NOC+2
                  FORM = 'BD'
                  IFB(I) = -1
C  X_I is confined in a very small interval.
              ELSE
                  DO 60,K=1,M
60                    FJ(K,I) = 0.0D0
                  FORM = '  '
                  IFB(I) = 2
                  IF (PR2) WRITE (*,100) I
              END IF
C  ---------------------------------------------------------------------
C  X_I is confined in a very small interval.
          ELSE
              DO 70,K=1,M
70                FJ(K,I) = 0.0D0
              FORM = '  '
              IF (PR2) WRITE (*,100) I
          END IF
          XX(I) = XI
          IF (PR3) WRITE (*,110) I, FORM, XI, HI
10    CONTINUE

c    CALL torc_enable_stealing()
      CALL torc_waitall()
c    CALL torc_disable_stealing()

C    !! ~SECOND PASS !!
C  Perform division using H()
      DO 2,I=1,N
         TWOHI = 2.0D0*HH(I)
         IF (IFB(I).EQ.0) THEN
            DO K=1,M
               FJ(K,I) = (FJ(K,I)-F1(K,I))/TWOHI
            END DO
         ELSE IF (IFB(I).EQ.1) THEN
            DO K=1,M
               FJ(K,I) = (FOUR * FJ(K,I) -THREE * F0(K) - F1(K,I))/TWOHI
            END DO
         ELSE IF (IFB(I).EQ.-1) THEN
            DO K=1,M
               FJ(K,I) = (-FOUR * FJ(K,I) + THREE * F0(K) +
     &                   F1(K,I))/TWOHI
            END DO
        END IF
2     CONTINUE

C
90    FORMAT (/' PNDL:',' Index',1X,'Formula',7X,'X_i',19X,'Step_i')
100   FORMAT (' PNDL: Warning: Variable ',I6,' is confined in a very ',
     &        'small interval.'
     &        /' PNDL:          Setting derivative to zero.')
110   FORMAT (' PNDL:',I6,3X,A,3X,1PG21.14,1X,1PG21.14)
      END
