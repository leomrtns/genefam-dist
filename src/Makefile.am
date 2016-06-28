AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libgenefamdist.la $(AM_LDFLAGS)

bin_PROGRAMS = distmatrix_genetree_sptree generate_spr_trees

distmatrix_genetree_sptree_SOURCES = distmatrix_genetree_sptree.c 
generate_spr_trees_SOURCES = generate_spr_trees.c
