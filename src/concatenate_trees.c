#include <genefam_dist.h>  // timestamp 2017.02.07 

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, j;  
  topology_space gt1 = NULL, gt2 = NULL;
  splitset split;

  biomcmc_random_number_init(0ULL);
  time0 = clock ();

  if (argc != 3) { 
    fprintf (stderr, "Concatenate two nexus tree files into one\n"); // TODO: include burnin, thinning, and remove outliers  
    fprintf (stderr, " USAGE: %s <first tree file> <second tree file>", basename (argv[0])); return EXIT_FAILURE; 
  }

  // read and order nexus_tree 
  gt1 = read_topology_space_from_file (argv[1], NULL);
  time1 = clock (); fprintf (stderr, "read timing 1: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;

  gt2 = read_topology_space_from_file (argv[2], gt1->taxlabel_hash);
  time1 = clock (); fprintf (stderr, "read timing 1: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;

  fprintf (stderr, "tree file 1 (unique) = %d (%d)\n", gt1->ntrees, gt1->ndistinct);
  fprintf (stderr, "tree file 1 (unique) = %d (%d)\n", gt2->ntrees, gt2->ndistinct);

  // We should test if external hashtable messes with translate numbers 
  merge_topology_spaces (gt1, gt2, (double)(gt1->ntrees)/(double)(gt2->ntrees));
  save_topology_space_to_trprobs_file (gt1, "concat.tre");


  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  del_topology_space (gt1);
  del_topology_space (gt2);
  biomcmc_random_number_finalize(); /* free the global variable */

  return EXIT_SUCCESS;
}

