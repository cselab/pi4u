C  ---------------------------------------------------------------------
	SUBROUTINE PNDLACC ( COMACC )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    This routine returns an estimation of the relative machine 
C    accuracy (COMACC).
C
C  Output arguments:
C    COMACC     The relative machine accuracy.
C
C  ---------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
	TEN = 10.0D0
	ONE = 1.0D0
	DO 10,I=1,1000
		COMACC = TEN**(-I)
		X1 = ONE+COMACC
		X2 = ONE-COMACC
		IF (X1.EQ.ONE .OR. X2.EQ.ONE) THEN
			COMACC = TEN**(-I+1)
			RETURN
		END IF
10	CONTINUE
	WRITE (*,*) 'ERROR: Cannot determine machine accuracy'
	STOP 'PNDL ERROR'
C
	END
