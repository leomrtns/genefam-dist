#include <genefam_dist.h>  // timestamp 2016.06.25 

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i;  
  topology_space gt = NULL;
  splitset split;

  biomcmc_random_number_init(0ULL);
  time0 = clock ();

  if (argc != 2) { 
    fprintf (stderr, "Calculates SPR distance between consecutive trees\n");
    fprintf (stderr, " USAGE: %s <nexus tree file>\n", basename (argv[0])); return EXIT_FAILURE; 
  }

  // read and order nexus_tree 
  gt = read_topology_space_from_file (argv[1], NULL, false); // false -> neglects root location
  split = new_splitset_dSPR (gt->tree[0]->nleaves);

  time1 = clock (); fprintf (stderr, "read timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;
  fprintf (stderr, "gene    trees (unique) = %d (%d)\n", gt->ntrees, gt->ndistinct);

  printf ("treeID   dSPR     dSPR_extra    RF   Hdist\n"); 
  for (i=0; i < gt->ntrees - 1; i++) {
    dSPR_topology (gt->tree[i], gt->tree[i+1], split);
    printf ("%9d %9d %9d %9d %9d\n", i, split->spr, split->spr_extra, split->rf, split->hdist);
  }

  del_splitset (split);
  del_topology_space (gt);
  biomcmc_random_number_finalize(); /* free the global variable */

  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
}

