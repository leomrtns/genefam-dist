#include <biomcmc.h> 

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, j, distance;  
  splitset split;

  time0 = clock ();
  if ((argc != 2) && (argc != 3)) 
    biomcmc_error ( "USAGE\n SPR distance between neighbor trees from file: %s <tree file>\nOR\n random SPR: %s <n_leaves> <n_simulations>", 
                    basename (argv[0]), basename(argv[0]));

  if (argc == 2) {
    topology_space gt = NULL;
    // read and order nexus_tree 
    gt = read_topology_space_from_file (argv[1], NULL);
    split = new_splitset_dSPR (gt->tree[0]->nleaves);
    time1 = clock (); fprintf (stderr, "read timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    time0 = time1;
    fprintf (stderr, "gene    trees (unique) = %d (%d)\n", gt->ntrees, gt->ndistinct);

    for (i = 0; i < gt->ndistinct - 1; i++) {
      distance = dSPR_topology (gt->distinct[i], gt->distinct[i+1], split);
      printf ("%9d %9d %9d %9d %9d\n", i, distance, split->hdist, split->rf, split->spr_extra);
    }
    del_topology_space (gt);
  }
  else {
    topology origtree, randtree; /* these trees DON'T have taxa labels */
    int n_leaves, n_iter, truedist, sprdist;

    biomcmc_random_number_init(0ULL);
    sscanf (argv[1], " %d ", &n_leaves);
    sscanf (argv[2], " %d ", &n_iter);
    origtree = new_topology (n_leaves);
    randtree = new_topology (n_leaves);
    split = new_splitset_dSPR (n_leaves);

    for (i=0; i < n_iter; i++) {
      randomize_topology (origtree);
      copy_topology_from_topology (randtree, origtree);
      truedist = biomcmc_rng_unif_int (n_leaves - 3) + 1;
      printf ("%9d ", truedist); fflush(stdout);
      for (j = 0; j < truedist; j++) topology_apply_spr_unrooted (randtree, false);
      sprdist = dSPR_topology (origtree, randtree, split);
      printf ("%9d %9d %9d %9d\n", sprdist, split->hdist, split->rf, split->spr_extra);
    }
    biomcmc_random_number_finalize(); /* free the global variable */
    del_topology (randtree);
    del_topology (origtree);
    del_splitset (split);
  }

  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}
