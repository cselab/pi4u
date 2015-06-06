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
	CALL torc_waitall()
	CALL torc_finalize()

	END

C  ---------------------------------------------------------------------
	SUBROUTINE PNDL_BARRIER ( )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Task-aware SPMD barrier.
C
C  ---------------------------------------------------------------------
C
	INCLUDE 'torcf.h'
C
c	CALL torc_enable()
c	CALL torc_waitall()
	CALL torc_spmd_barrier()
c	CALL torc_disable()

	END

C  ---------------------------------------------------------------------
	SUBROUTINE PNDL_TASKWAIT ( )
C  ---------------------------------------------------------------------
C
C  Description:                                    PNDL internal routine.
C                                                  ---------------------
C    Task-aware SPMD barrier.
C
C  ---------------------------------------------------------------------
C
	INCLUDE 'torcf.h'
C
	CALL torc_waitall()

	END
