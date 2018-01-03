C  ---------------------------------------------------------------------
      SUBROUTINE PNDLJA ( RSD, X, N, M, XL, XU, UH, FEPS, IORD, IPRINT,
     &                   FJ, LD, NOC, IERR )
C  ---------------------------------------------------------------------
C
C  Description:                              PNDL user interface routine.
C                                            ---------------------------
C    Given a multidimensional function that is written a sum of squared
C    terms (residuals):
C                  M
C          F(X) = Sum R_i(X)**2
C                 i=1
C    this routine returns the Jacobian matrix (FJ) by applying a
C    numerical differentiation formula according to the desired order
C    of accuracy (IORD).
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
C    UH         Array containing user-speficied steps.
C    FEPS       An estimation of the relative accuracy with which
C               the function is evaluated.
C               If FEPS=0 the machine accurracy is computed and used
C               instead.
C    IORD       Order of the numerical derivative.
C               Possible values: 1, 2.
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
C    IERR       Error indicator. Possible values are:
C                 IERR=0 -> No errors at all.
C                 IERR=1 -> The supplied IORD is incorrect.
C                 IERR=2 -> The supplied N is less than 1.
C                 IERR=3 -> Some of the supplied upper bounds (XU)
C                           are less than the lower bounds (XL).
C                 IERR=4 -> Some of the supplied values in X do
C                           not lie inside the bounds.
C                 IERR=5 -> The supplied value of FEPS is incorrect
C                           (less than 0 or greater than 1).
C                 IERR=6 -> The supplied IPRINT is incorrect.
C                 IERR=7 -> Not used
C                 IERR=8 -> Not used
C                 IERR=9 -> The supplied number of squared terms (M)
C                           is less than 1.
C
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RSD
      DIMENSION X(N), XL(N), XU(N), UH(N), FJ(LD,N)
      EXTERNAL PNDLJF, PNDLJQ
      INCLUDE 'torcf.h'

C     Informative message (routine starts executing).
      CALL PNDLMSG(IPRINT,'PNDLJA',0)
      IF (IPRINT.EQ.3) WRITE (*,10) N, M, FEPS, IORD

C     Check validity of arguments.
      CALL PNDLARG(X,N,XL,XU,FEPS,IPRINT,IERR)
      IF (IERR.NE.0) GOTO 100

C     Print bounds.
      CALL PNDLBND(X,XL,XU,N,IPRINT)

C     Make sure M is at least 1.
      IF (M.LT.1) THEN
          IF (IPRINT.GT.0) WRITE (*,20) M
          IERR = 9
          GOTO 100
      END IF

C     Call a routine to do the actual computation according to the
C     desired order of accuracy.
      IF (IORD.EQ.1) THEN
         CALL PNDLJF(RSD,X,N,M,XL,XU,UH,FEPS,IPRINT,FJ,LD,NOC)
      ELSE IF (IORD.EQ.2) THEN
         CALL PNDLJQ(RSD,X,N,M,XL,XU,UH,FEPS,IPRINT,FJ,LD,NOC)
      ELSE
          IF (IPRINT.GT.0) WRITE (*,30) IORD
          IERR = 1
      END IF
C
100   CALL PNDLMSG(IPRINT,'PNDLJA',1)
C
10    FORMAT (' PNDL: ',3X,'N = ',I6,3X,'M = ',I6,3X,'FEPS = ',1PG21.14,
     &        3X,'IORD = ',I2)
20    FORMAT (/' PNDL: Error: Number of functions (M) is less than 1'/)
30    FORMAT (/' PNDL: Error: Incorrect order (IORD) ',I6/)
      END SUBROUTINE PNDLJA
