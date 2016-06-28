#include <genefam_dist.h>  // timestamp 2016.06.25 

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, j;  
  topology_space gt = NULL, st = NULL;
  splitset split;

  biomcmc_random_number_init(0ULL);
  time0 = clock ();

  if (argc != 3) { 
    fprintf (stderr, "Calculates set of distances between all tree pairs\n");
    fprintf (stderr, " USAGE: %s <nexus gene tree file> <nexus species tree file>", basename (argv[0])); return EXIT_FAILURE; 
  }

  // read and order nexus_tree 
  gt = read_topology_space_from_file (argv[1], NULL);
  st = read_topology_space_from_file (argv[2], NULL);

  time1 = clock (); fprintf (stderr, "read timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;
  fprintf (stderr, "gene    trees (unique) = %d (%d)\n", gt->ntrees, gt->ndistinct);
  fprintf (stderr, "species trees (unique) = %d (%d)\n", st->ntrees, st->ndistinct);

  split = create_splitset_dSPR_genespecies (gt->distinct[0], st->distinct[0]);
  printf ("genetree sptree  ndups  nloss ndeepcoals dSPR_lower dSPR_extra  RF   Hdist\n"); 
  for (i=0; i < st->ntrees; i++) {
    for (j = 0; j < gt->ntrees; j++) { 
      init_tree_recon_from_species_topology (gt->tree[j], st->tree[i]);
      dSPR_gene_species (gt->tree[j], st->tree[i], split);
      printf ("%6d %6d   %5d %5d %5d   %5d %5d %5d %6d\n", j, i, gt->tree[j]->rec->ndups, gt->tree[j]->rec->nloss, gt->tree[j]->rec->ndcos, 
              split->spr, split->spr_extra, split->rf, split->hdist); 
    }
  }

  del_splitset (split);

  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  del_topology_space (gt);
  del_topology_space (st);
  biomcmc_random_number_finalize(); /* free the global variable */

  return EXIT_SUCCESS;
}

