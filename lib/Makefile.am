AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS)

common_headers = hashtable.h lowlevel.h random_number_gen.h random_number.h nexus_common.h \
                 topology_common.h topology_space.h topology_mrca.h bipartition.h topology_build.h \
                 topology_splitset.h prob_distribution.h empirical_frequency.h genefam_dist.h

common_src     = hashtable.c lowlevel.c random_number_gen.c random_number.c nexus_common.c \
                 topology_common.c topology_space.c topology_mrca.c bipartition.c topology_build.c \
                 topology_splitset.c prob_distribution.c empirical_frequency.c genefam_dist.c

lib_LTLIBRARIES = libgenefamdist.la
otherincludedir = $(includedir)/genefam
otherinclude_HEADERS = config.h $(common_headers)
libgenefamdist_la_SOURCES = $(common_src) $(common_headers)

