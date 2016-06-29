#include <genefam_dist.h> // timestamp 2016.06.28 

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, j, n_leaves, n_iter, n_spr;
  char *s;
  splitset split;

  time0 = clock ();
  if (argc != 4) {
    fprintf (stderr, "Generate random trees with a given number of SPR operations apart\n");
    fprintf (stderr, "USAGE: %s <n_leaves> <n_simulations> <n_SPR>\n",  basename (argv[0])); return EXIT_FAILURE;
  }

  topology origtree, randtree; /* these trees don't have taxa labels */

  biomcmc_random_number_init(0ULL);
  sscanf (argv[1], " %d ", &n_leaves);
  sscanf (argv[2], " %d ", &n_iter);
  sscanf (argv[3], " %d ", &n_spr);
  split = new_splitset_dSPR (n_leaves);

  origtree = new_topology (n_leaves);
  randtree = new_topology (n_leaves);

  randomize_topology (randtree);
  s = topology_to_string_create_name (randtree, NULL); /* second parameter is vector with branch lengths */
  printf ("#NEXUS\nBegin trees;\ntree PAUP_1 = ");
  printf ("%s\n",s); fflush(stdout); free (s);

  for (i=0; i < n_iter; i++) {
    copy_topology_from_topology (origtree, randtree);
    for (j = 0; j < n_spr; j++) topology_apply_spr_unrooted (randtree, false);

    s = topology_to_string_create_name (randtree, NULL); /* second parameter is vector with branch lengths */
    printf ("tree PAUP_%d = %s\n", i+2, s); fflush(stdout); free (s);
  }
  printf ("End;\n");
  biomcmc_random_number_finalize(); /* free the global variable */
  del_topology (randtree);
  del_topology (origtree);
  del_splitset (split);

  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}
