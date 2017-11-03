AC_INIT(pndl, 1.0.0,,pndl)
AC_CONFIG_AUX_DIR(auxdir)
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE


noFFLAGS=""
if test "x${FFLAGS}" = x; then
        noFFLAGS=yes
fi

noCFLAGS=""
if test "x${CFLAGS}" = x; then
        noCFLAGS=yes
fi

AC_CONFIG_HEADER(src/pndl_config.h)
AC_PROG_CC(mpicc)
AC_PROG_F77(mpif77)
AC_PROG_SED
AC_LANG_C
AC_PROG_CPP

dnl AC_LANG(Fortran 77)
AC_F77_WRAPPERS
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AM_PROG_AS

G77OPT=""
if test x$G77 = xyes; then
	if test "x${noFFLAGS}" = xyes; then
		FFLAGS="-g -O0"
	fi
#	G77OPT="-O0"
fi

if test x$GCC = xyes; then
	if test "x${noCFLAGS}" = xyes; then
		CFLAGS="-g -O0"
	fi
#	CFLAGS="${CFLAGS} -O0"
fi


AC_ARG_WITH(mpi, 
	[  --with-mpi=<dir>        Location of the MPI installation],
	if test x$withval = xyes; then
		AC_MSG_ERROR(--with-mpi must be given a pathname)
	else
		CFLAGS="${CFLAGS} -I${withval}/include"
		FFLAGS="${FFLAGS} -I${withval}/include"
		LDFLAGS="${LDFLAGS} -L${withval}/lib"
		MPILIB="-lmpich"
	fi)

AC_ARG_WITH(mpilib, 
	[  --with-mpilib=<mpilib>        Name of MPI library],
	if test x$withval = xyes; then
		AC_MSG_WARN(--with-mpilib must be given a library. Using default: mpich)
		MPILIB="-lmpich"
	else
		MPILIB="-l${withval}"
	fi)


AC_ARG_WITH(mpiincdir, 
	[  --with-mpiincdir=<mpiincdir>        Directory of MPI header files],
	if test x$withval = xyes; then
		AC_MSG_ERROR(--with-mpiincdir must be given a directory)
	else
		CFLAGS="${CFLAGS} -I${withval}"
		FFLAGS="${FFLAGS} -I${withval}"
	fi)


AC_ARG_WITH(mpilibdir, 
	[  --with-mpilibdir=<mpilibdir>        Directory of MPI library],
	if test x$withval = xyes; then
		AC_MSG_ERROR(--with-mpilibdir must be given a directory)
	else
		LDFLAGS="${LDFLAGS} -L${withval}"
	fi)

AC_ARG_ENABLE(debug,
	AC_HELP_STRING([--enable-debug], [enable debugging (default:disabled)]),
	if eval "test x$enable_debug = xyes"; then
		DEBUGFLAG="-g"
	else
		DEBUGFLAG=""
	fi,
	DEBUGFLAG="")


LDFLAGS="${LDFLAGS} ${MPILIB}"
LDFLAGS="${LDFLAGS} -lpthread"

echo "Making with LDFLAGS=$LDFLAGS"

AC_SUBST(G77OPT)
AC_SUBST(DEBUG_FLAG)
AC_SUBST(MPICC)
AC_SUBST(MPIF77)
AC_SUBST(MPILIB)

# Check for MPI
AC_MSG_CHECKING([for MPI/Fortran77 support])
mpif_enabled=""
AC_LANG_PUSH(Fortran 77)
AC_LINK_IFELSE(
[
        PROGRAM test
	integer ERROR
	include 'mpif.h'
        print *, 'hello'
	call MPI_INIT(ERROR)
	call MPI_FINALIZE(ERROR)
        END
],
[mpif_enabled=yes], [mpif_enabled=no])
AC_LANG_POP(Fortran 77)

#echo $mpif_enabled

if test ".$mpif_enabled" = ".yes"; then
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
	AC_MSG_FAILURE([MPI/Fortran77 not enabled])
fi


AC_OUTPUT(Makefile src/Makefile test/Makefile)
