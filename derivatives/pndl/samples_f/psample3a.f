c
c     A single PNDL CALL issued by rank 0
c
      PROGRAM ptest
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      PARAMETER (nn=2)
      DIMENSION x1(nn), g1(nn)
      DIMENSION x2(nn), g2(nn)
      DIMENSION x3(nn), g3(nn)
      INTEGER n
      EXTERNAL f
      EXTERNAL pndlg
      INCLUDE 'torcf.h'

      CALL pndl_init()
      CALL torc_register_task(pndlg)
      CALL torc_initf()
c   ...
      x1(1) = 1.0d0
      x1(2) = 1.1d0
      iord = 2
      n = nn


      CALL torc_taskf(pndlg, 1, 5,
     &               1, MPI_INTEGER,          CALL_BY_VAD,
     &               N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &               1, MPI_INTEGER,          CALL_BY_VAL,
     &               1, MPI_INTEGER,          CALL_BY_VAL,
     &               n, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &               f, x1, n, iord, g1)

      x2(1) = 2.0d0
      x2(2) = 2.1d0
      iord = 2

      CALL torc_taskf(pndlg, 1, 5,
     &               1, MPI_INTEGER,          CALL_BY_VAD,
     &               N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &               1, MPI_INTEGER,          CALL_BY_VAL,
     &               1, MPI_INTEGER,          CALL_BY_VAL,
     &               n, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &               f, x2, n, iord, g2)

      x3(1) = 3.0d0
      x3(2) = 3.1d0
      iord = 2

      CALL torc_taskf(pndlg, 1, 5,
     &               1, MPI_INTEGER,          CALL_BY_VAD,
     &               N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &               1, MPI_INTEGER,          CALL_BY_VAL,
     &               1, MPI_INTEGER,          CALL_BY_VAL,
     &               n, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &               f, x3, n, iord, g3)

      CALL torc_enable_stealing()
      CALL torc_waitall()
c   ...
      PRINT *, 'point ', (x1(i), i=1, n)
      PRINT *, 'gradient ', (g1(i), i=1, n)

      PRINT *, 'point ', (x2(i), i=1, n)
      PRINT *, 'gradient ', (g2(i), i=1, n)

      PRINT *, 'point ', (x3(i), i=1, n)
      PRINT *, 'gradient ', (g3(i), i=1, n)

c   ...

      CALL torc_finalize()
      CALL pndl_finalize()
      END
c----------------------------------------------------------

      FUNCTION f(x,n)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      DIMENSION x(n)
      CALL torc_sleep(100)
      f = x(1)*COS(x(2))+x(2)*COS(x(1))
      END
