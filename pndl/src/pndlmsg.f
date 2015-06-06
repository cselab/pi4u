C  ---------------------------------------------------------------------
	SUBROUTINE PNDLMSG ( IPRINT, NAME, IEND )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Prints a message indicating the beggining or termination
C    of a library routine, taking into account the printout
C    control option IPRINT.
C
C  Input arguments:
C    IPRINT     Printout control option.
C    NAME       Name of the routine that starts/ends.
C    IEND       IEND=0 -> Indicates that the routine starts.
C               IEND=1 -> Indicates that the routine ends.
C
C  ---------------------------------------------------------------------
	CHARACTER*(*) NAME
C
	IF (IPRINT.NE.3) RETURN
C
	IF (IEND.EQ.0) THEN
		WRITE (*,10) NAME
	ELSE
		WRITE (*,20) NAME
	END IF
C
10	FORMAT (/' PNDL: ',75('-')
     &        /' PNDL: ',22X,'Start of SUBROUTINE ',A
     &        /' PNDL: ',75('-'))
20	FORMAT (/' PNDL: ',75('-')
     &        /' PNDL: ',23X,'End of SUBROUTINE ',A
     &        /' PNDL: ',75('-'))
	END
