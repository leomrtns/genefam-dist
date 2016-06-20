#include <genefam_dist.h>

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, j;  
  topology_space gt = NULL, st = NULL;
  splitset split;
  int dSPR;
  FILE *stream;

  biomcmc_random_number_init(0ULL);
  time0 = clock ();

  if (argc != 3) 
    biomcmc_error ( " USAGE: %s <gene tree file> <species tree file>", basename (argv[0]));

  // read and order nexus_tree 
  gt = read_topology_space_from_file (argv[1], NULL);
  st = read_topology_space_from_file (argv[2], NULL);

  time1 = clock (); fprintf (stderr, "read timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;
  fprintf (stderr, "gene    trees (unique) = %d (%d)\n", gt->ntrees, gt->ndistinct);
  fprintf (stderr, "species trees (unique) = %d (%d)\n", st->ntrees, st->ndistinct);
  
//  for (i=0; i < gt->ntrees; i++) printf ("genetree:: %8d  %8d  %12.8lf\n", i, gt->tree[i]->id, gt->freq[ gt->tree[i]->id ]);
//  for (i=0; i < st->ntrees; i++) printf ("  sptree:: %8d  %8d  %12.8lf\n", i, st->tree[i]->id, st->freq[ st->tree[i]->id ]);

  stream = biomcmc_fopen ("pairwise.txt", "w");
  split = create_splitset_dSPR_genespecies (gt->distinct[0], st->distinct[0]);
  fprintf (stream, "genetree sptree  ndups  nloss ndeepcoals dSPR  RF   Hdist\n"); 
  for (i=0; i < st->ntrees; i++) {
    for (j = 0; j < gt->ntrees; j++) { 
      init_tree_recon_from_species_topology (gt->tree[j], st->tree[i]);
      dSPR_gene_species (gt->tree[j], st->tree[i], split);
      dSPR = split->spr + split->spr_extra;
      fprintf (stream, "%6d %6d   %5d %5d %5d   %5d %5d %6d\n", j, i, gt->tree[j]->rec->ndups, gt->tree[j]->rec->nloss, gt->tree[j]->rec->ndcos, dSPR, split->rf, split->hdist); 
    }
  }

  del_splitset (split);
  fclose (stream);

  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  del_topology_space (gt);
  del_topology_space (st);
  biomcmc_random_number_finalize(); /* free the global variable */

  return EXIT_SUCCESS;
}

