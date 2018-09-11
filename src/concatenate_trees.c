#include <genefam_dist.h>  // timestamp 2017.02.07 

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int burnin = 0, thin = 1;  
  double credible = 1.;
  topology_space gt1 = NULL, gt2 = NULL;

  biomcmc_random_number_init(0ULL);
  time0 = clock ();

  if (argc != 6) { 
    fprintf (stderr, "Concatenate two nexus tree files into one\n"); 
    fprintf (stderr, "   USAGE: %s <first tree file> <second tree file> <burnin> <thin> <credible interval>\n", basename (argv[0])); 
    fprintf (stderr, "(where 'thin' is the frequency at which we sample, and 'credible'<=1 is the cutoff for the cummulative posterior)\n"); 
    return EXIT_FAILURE; 
  }
  sscanf (argv[3], " %d ", &burnin);
  sscanf (argv[4], " %d ", &thin);
  sscanf (argv[5], " %lf ", &credible);
  if (credible < 0.00001) credible = 0.00001;

  // read and order nexus_tree 
  gt1 = read_topology_space_from_file_with_burnin_thin (argv[1], NULL, burnin, thin, true); // true-> assume species trees, respect root location
  time1 = clock (); fprintf (stderr, "read timing 1: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;

  gt2 = read_topology_space_from_file_with_burnin_thin (argv[2], gt1->taxlabel_hash, burnin, thin, true);
  time1 = clock (); fprintf (stderr, "read timing 2: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;

  fprintf (stderr, "tree file 1 (unique) = %d (%d)\n", gt1->ntrees, gt1->ndistinct);
  fprintf (stderr, "tree file 2 (unique) = %d (%d)\n", gt2->ntrees, gt2->ndistinct);

  // We should test if external hashtable messes with translate numbers 
  merge_topology_spaces (gt1, gt2, (double)(gt1->ntrees)/(double)(gt2->ntrees), true); // true --> preserves rooting
  save_topology_space_to_trprobs_file (gt1, "concat.tre", credible);


  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  del_topology_space (gt1);
  del_topology_space (gt2);
  biomcmc_random_number_finalize(); /* free the global variable */

  return EXIT_SUCCESS;
}

