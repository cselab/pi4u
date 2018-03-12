C  ---------------------------------------------------------------------
      SUBROUTINE PNDLBND ( X, XL, XU, N, IPRINT )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Prints the values of the function variables (X), along with
C    their lower (XL) and upper (XU) bounds. The routine takes into
C    account the printout control option IPRINT.
C
C  Input arguments:
C    X          Array containing the function variables.
C    XL         Array containing lower bounds on the variables.
C    XU         Array containing upper bounds on the variables.
C    N          The number of variables.
C    IPRINT     Printout control option. Possible values are:
C                 IPRINT.NE.3 -> No printout at all.
C                 IPRINT=3    -> Print variables and bounds.
C
C  ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N), XL(N), XU(N)
C
      IF (IPRINT.NE.3) RETURN
C
      WRITE (*,100)
      DO 10,I=1,N
          WRITE (*,200) I, XL(I), X(I), XU(I)
10    CONTINUE
C
100    FORMAT (/' PNDL: ',' Index',9X,'XL_i',20X,'X_i',19X,'XU_i')
200    FORMAT (' PNDL: ',I6,2X,1PG21.14,2X,1PG21.14,2X,1PG21.14)
      END

