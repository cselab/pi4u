bin_PROGRAMS=pt1 pt2 ptj

AM_CFLAGS = @DEBUG_FLAG@ 

LLIBS = -L../src -lpndl -L../../../torc_lite/src -ltorc

pt1_SOURCES = pt1.f
pt1_LDADD =$(LLIBS) 
pt1_FFLAGS = -I../../../torc_lite/include

pt2_SOURCES = pt2.f
pt2_LDADD = $(LLIBS)
pt2_FFLAGS = -I../../../torc_lite/include

ptj_SOURCES = ptj.f
ptj_LDADD = $(LLIBS)
ptj_FFLAGS = -I../../../torc_lite/include
