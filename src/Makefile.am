AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libgenefamdist.la $(AM_LDFLAGS)

bin_PROGRAMS = test_python_functions gf_concatenate_trees gf_find_best_trees \
							 gf_distmatrix_genetree_sptree gf_generate_spr_trees gf_spr_distance gf_distsignal_genetree_sptree 


test_python_functions_SOURCES = test_python_functions.c
gf_concatenate_trees_SOURCES = concatenate_trees.c
gf_find_best_trees_SOURCES = find_best_trees.c
gf_distmatrix_genetree_sptree_SOURCES = distmatrix_genetree_sptree.c 
gf_distsignal_genetree_sptree_SOURCES = distsignal_genetree_sptree.c 
gf_generate_spr_trees_SOURCES = generate_spr_trees.c
gf_spr_distance_SOURCES = spr_distance.c
