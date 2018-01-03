c
c	A single PNDL call issued by rank 0
c
	program ptest
	implicit double precision (a-h, o-z)
	parameter (n=2)
	dimension x1(n), g1(n)
	dimension x2(n), g2(n)
	dimension x3(n), g3(n)
	external f
	include 'torcf.h'

	call pndl_init()
	call torc_initf()
c   ...
	x1(1) = 1.0d0
	x1(2) = 1.1d0
	iord = 2
	call pndlg ( f, x1, n, iord, g1 )

	x2(1) = 2.0d0
	x2(2) = 2.1d0
	iord = 2
	call pndlg ( f, x2, n, iord, g2 )


	x3(1) = 3.0d0
	x3(2) = 3.1d0
	iord = 2
	call pndlg ( f, x3, n, iord, g3 )

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
