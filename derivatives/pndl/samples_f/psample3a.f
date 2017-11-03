c
c	A single PNDL call issued by rank 0
c
	program ptest
	implicit double precision (a-h, o-z)
	parameter (nn=2)
	dimension x1(nn), g1(nn)
	dimension x2(nn), g2(nn)
	dimension x3(nn), g3(nn)
	integer n
	external f
	external pndlg
	include 'torcf.h'

	call pndl_init()
	call torc_register_task(pndlg)
	call torc_init()
c   ...
	x1(1) = 1.0d0
	x1(2) = 1.1d0
	iord = 2
	n = nn


        call torc_task(pndlg, 1, 5,
     &      1, MPI_INTEGER, CALL_BY_VAD,
     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &      1, MPI_INTEGER, CALL_BY_VAL,
     &      1, MPI_INTEGER, CALL_BY_VAL,
     &      n, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &      f,x1,n,iord,g1)

	x2(1) = 2.0d0
	x2(2) = 2.1d0
	iord = 2

        call torc_task(pndlg, 1, 5,
     &      1, MPI_INTEGER, CALL_BY_VAD,
     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &      1, MPI_INTEGER, CALL_BY_VAL,
     &      1, MPI_INTEGER, CALL_BY_VAL,
     &      n, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &      f,x2,n,iord,g2)

	x3(1) = 3.0d0
	x3(2) = 3.1d0
	iord = 2

        call torc_task(pndlg, 1, 5,
     &      1, MPI_INTEGER, CALL_BY_VAD,
     &      N, MPI_DOUBLE_PRECISION, CALL_BY_VAL,
     &      1, MPI_INTEGER, CALL_BY_VAL,
     &      1, MPI_INTEGER, CALL_BY_VAL,
     &      n, MPI_DOUBLE_PRECISION, CALL_BY_RES,
     &      f,x3,n,iord,g3)

	call torc_enable_stealing()
	call torc_waitall()
c   ...
	print *, 'point ', (x1(i), i=1, n)
	print *, 'gradient ', (g1(i), i=1, n)

	print *, 'point ', (x2(i), i=1, n)
	print *, 'gradient ', (g2(i), i=1, n)

	print *, 'point ', (x3(i), i=1, n)
	print *, 'gradient ', (g3(i), i=1, n)

c   ...

	call torc_finalize()
	call pndl_finalize()
	end
c----------------------------------------------------------

	function f(x,n)
	implicit double precision (a-h, o-z)
	dimension x(n)
	call torc_sleep(100)
	f = x(1)*cos(x(2))+x(2)*cos(x(1))
	end
