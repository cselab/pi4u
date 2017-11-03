C  ---------------------------------------------------------------------
	SUBROUTINE PNDL_INIT ( )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Initializes the underlying tasking library.
C
C  ---------------------------------------------------------------------
C
	INCLUDE 'torcf.h'
C
C	CALL torc_init()
	EXTERNAL TASKFUN
C       GRD, RSD

	call torc_register_task(TASKFUN)
c	call torc_register_task(GRD)
c	call torc_register_task(RSD)


	END

C  ---------------------------------------------------------------------
	SUBROUTINE PNDL_FINALIZE ( )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Shutdowns the underlying tasking library.
C
C  ---------------------------------------------------------------------
C
	INCLUDE 'torcf.h'
C

	END

