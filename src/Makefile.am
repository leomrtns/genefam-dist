AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libgenefamdist.la $(AM_LDFLAGS)

bin_PROGRAMS = distance
distance_SOURCES = main_distance.c 
