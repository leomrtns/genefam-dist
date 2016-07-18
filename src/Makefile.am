AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libgenefamdist.la $(AM_LDFLAGS)

bin_PROGRAMS = gf_distmatrix_genetree_sptree gf_generate_spr_trees gf_spr_distance

gf_distmatrix_genetree_sptree_SOURCES = distmatrix_genetree_sptree.c 
gf_generate_spr_trees_SOURCES = generate_spr_trees.c
gf_spr_distance_SOURCES = spr_distance.c
