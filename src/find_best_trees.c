#include <genefam_dist.h> 
/* based on guenomu/maxtree but also implementing elements from guenomu in ML mode */

#define NDISTS 6

typedef struct pooled_matrix_struct* pooled_matrix;
typedef struct genetree_struct* genetree;
typedef struct gene_sptrees_struct* gene_sptrees;

/* idea: store bad  topologies and reject those close to them (RF between sptrees cheaper) <- too many?
 * 'bad' in the sense that appeared in proposal (i.e. not random) but rejected (worse than any in ratchet?); we can even output these at end */

// unlike genefam_dist library, here we do not have ref_counter 
struct pooled_matrix_struct
{
  int n_sptrees, n_species, n_sets_per_gene, next; 
  spdist_matrix *d, d_total, this_gene_spdist;
  distance_matrix square_matrix;
};

struct genetree_struct
{
  topology tree;
  double minmax[2 * NDISTS], dist[NDISTS]; // values are integers; might have a variable "scale" to allow for quick or no scaling 
  splitset split;
};

struct gene_sptrees_struct
{ // we may need to have something like genetree_struct since trees must have pointer to distances/score (to sort in locu) or we use empfreq
  int n_genes, n_ratchet, n_proposal;
  int next_proposal, next_free; // idx to ratchet elements
  genetree *gene;
  topology *ratchet, *proposal; // may need vector of distances
  double best_score; // I'm forgetting other vars like minmax etc...
};

void print_usage (char *progname);
pooled_matrix new_pooled_matrix (int n_sets, int n_species);
void del_pooled_matrix (pooled_matrix pool);
void update_pooled_matrix_from_gene_tree (pooled_matrix pool, topology gene_topol, int *idx_gene_to_species);
void finalise_pooled_matrix (pooled_matrix pool, char_vector sptaxa);
void find_maxtree_and_add_to_topol_space (spdist_matrix dist, pooled_matrix pool, topology_space tsp, char_vector sptaxa, int tree_method, bool use_within_gf_means);
topology_space maxtrees_from_subsamples (char_vector sptaxa, char **arg_filename, int arg_nfiles, char_vector gfilename);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  char_vector species, genefiles;
  topology_space tsp;

  time0 = clock ();
  if (argc < 3) print_usage (argv[0]);
  biomcmc_random_number_init(0);

  species = new_char_vector_from_file (argv[1]);
  genefiles = new_char_vector (1); // updated by maxtrees_from_subsamples()
  char_vector_remove_duplicate_strings (species); /* duplicate names */
  char_vector_longer_first_order (species);
  tsp = maxtrees_from_subsamples (species, argv + 2, argc - 2, genefiles); 
  save_topology_space_to_trprobs_file (tsp, "patristic.tre", 1.);
 

  del_char_vector (species);
  del_char_vector (genefiles);
  del_topology_space (tsp);
  biomcmc_random_number_finalize();

	time1 = clock (); fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

	return EXIT_SUCCESS;
}

void
print_usage (char *progname)
{
  fprintf (stderr, "\tUSAGE: %s <species file> <gene tree 1> ... <gene tree N>\n\n", basename (progname));
  fprintf (stderr, "where  < species file >      is the name of a text file with the names of the species\n");
  fprintf (stderr, "       < gene tree ?? >      are file names for the gene family trees\n");
  fprintf (stderr, "\n  - remember that the species names should be found within the gene names\n");
  exit (EXIT_FAILURE);
}

topology_space
maxtrees_from_subsamples (char_vector sptaxa, char **arg_filename, int arg_nfiles, char_vector gfilename)
{
  int i, *idx_gene_to_sp = NULL, valid_species_size = 0;
  topology_space tsp, genetre;
  pooled_matrix pool;
  /* order species names from longer names to shorter (so that EColi is searched only after EColiII, e.g.) */
  empfreq ef = new_empfreq_sort_decreasing (sptaxa->nchars, sptaxa->nstrings, 1); /* 1=size_t (0=char, 2=int)*/

  pool = new_pooled_matrix ((int)(arg_nfiles/100), sptaxa->nstrings);

  for (i = 0; i < arg_nfiles; i++) { /* scan through gene files */
    genetre  = read_topology_space_from_file (arg_filename[i], NULL);/* read i-th gene family */
    
    idx_gene_to_sp = (int *) biomcmc_realloc ((int*) idx_gene_to_sp, genetre->distinct[0]->nleaves * sizeof (int));/* for each gene leaf, index of species */
    index_sptaxa_to_genetaxa (sptaxa, genetre->taxlabel, idx_gene_to_sp, ef);/* map species names to gene names and store into idx[] */
    valid_species_size = prepare_spdistmatrix_from_gene_species_map (pool->this_gene_spdist, idx_gene_to_sp, genetre->distinct[0]->nleaves);
    if (valid_species_size > 4) {
      char_vector_add_string (gfilename, arg_filename[1]);
      update_pooled_matrix_from_gene_tree (pool, genetre->distinct[0], idx_gene_to_sp);
    }
    del_topology_space (genetre); 
  }

  finalise_pooled_matrix (pool, sptaxa);

  tsp = new_topology_space ();
  tsp->taxlabel = sptaxa; sptaxa->ref_counter++; /* sptaxa is shared here as well */
  for (i =0; i < 3; i++) {
    find_maxtree_and_add_to_topol_space (pool->d_total, pool, tsp, sptaxa, i, true);  // mean length within gene family (locus)
    find_maxtree_and_add_to_topol_space (pool->d_total, pool, tsp, sptaxa, i, false); // min lengths within gene family
  }
  for (i = 0; i < pool->n_sptrees; i++) {
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, sptaxa, 0, true); // bionj/upgma/singlelinkage, means/mins
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, sptaxa, 1, false); // upgma on mins 
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, sptaxa, 2, false); // single linkage on mins 
  }
  del_empfreq (ef);
  del_pooled_matrix (pool);
  if (idx_gene_to_sp) free (idx_gene_to_sp);
  return tsp;
}

/* pooled_matrix_struct functions */

pooled_matrix
new_pooled_matrix (int n_sets, int n_species)
{
  int i;
  pooled_matrix pool;

  pool = (pooled_matrix) biomcmc_malloc (sizeof (struct pooled_matrix_struct));

  pool->n_sptrees = n_sets;
  pool->n_species = n_species;
  pool->n_sets_per_gene = n_sets/5;
  pool->next = 0; // vector is shuffled whenever next arrives at last element

  if (pool->n_sptrees < 1) pool->n_sptrees = 1; 
  if (pool->n_sets_per_gene < 1) pool->n_sets_per_gene = 1;

  pool->this_gene_spdist = new_spdist_matrix (n_species);
  pool->d_total = new_spdist_matrix (n_species);
  pool->square_matrix = new_distance_matrix (n_species);
  pool->d = (spdist_matrix*) biomcmc_malloc (pool->n_sptrees * sizeof (spdist_matrix));
  for (i = 0; i < pool->n_sptrees; i++) {
    pool->d[i] = new_spdist_matrix(n_species);
    zero_all_spdist_matrix (pool->d[i]);
  }
  printf ("Creating %d distance matrices with each gene family participating in %d random ones\n", pool->n_sptrees, pool->n_sets_per_gene);
  return pool;
}

void
del_pooled_matrix (pooled_matrix pool)
{
  int i;
  if (! pool) return;
  if (pool->d) {
    for (i = pool->n_sptrees-1; i >=0; i--) del_spdist_matrix (pool->d[i]);
    free (pool->d);
  }
  del_spdist_matrix (pool->this_gene_spdist);
  del_spdist_matrix (pool->d_total);
  del_distance_matrix (pool->square_matrix);
  free (pool);
  return;
}

int
next_element_from_pooled_matrix (pooled_matrix pool)
{
  if (++pool->next < pool->n_sptrees) return pool->next;

  int i, j;
  spdist_matrix pivot;
  for (i = pool->n_sptrees - 1; i > 0; i--) {  // Knuth shuffle
    j = biomcmc_rng_unif_int (i+1); // including i
    pivot = pool->d[j];
    pool->d[j] = pool->d[i];
    pool->d[i] = pivot;
  }
  pool->next = 0;
  return pool->next;
}

void
update_pooled_matrix_from_gene_tree (pooled_matrix pool, topology gene_topol, int *idx_gene_to_species)
{
  int i, j;
  distance_matrix genedist = new_distance_matrix_for_topology (gene_topol->nleaves);
  // here we can play with branch lengths (now it's NULL)
  fill_distance_matrix_from_topology (genedist, gene_topol, NULL, true); /* true=store in upper diagonal */
  fill_spdistmatrix_from_gene_dists (pool->this_gene_spdist, genedist, idx_gene_to_species, true);
  for (i = 0; i < pool->n_sets_per_gene; i++) {
    j = next_element_from_pooled_matrix (pool); // shuffles pool->d[] when necessary 
    update_spdistmatrix_from_spdistmatrix (pool->d[j], pool->this_gene_spdist);
  }
  update_spdistmatrix_from_spdistmatrix (pool->d_total, pool->this_gene_spdist);
  del_distance_matrix (genedist);
}

void
finalise_pooled_matrix (pooled_matrix pool, char_vector sptaxa)
{
  int i;
  finalise_spdist_matrix (pool->d_total);
  if (pool->d_total->n_missing) fprintf (stderr, "OBS: %d species pair combinations never appear on same gene familiy\n", pool->d_total->n_missing);
  for (i = 0; i < sptaxa->nstrings; i++) if (!pool->d_total->species_present[i])
    fprintf (stderr, "OBS: species \'%s\' never appears in data set; consider removing it from name list\n", sptaxa->string[i]);
  for (i = 0; i < pool->n_sptrees; i++) {
    finalise_spdist_matrix (pool->d[i]);
    complete_missing_spdist_from_global_spdist (pool->d[i], pool->d_total);
  }
}

void
find_maxtree_and_add_to_topol_space (spdist_matrix dist, pooled_matrix pool, topology_space tsp, char_vector sptaxa, int tree_method, bool use_within_gf_means)
{
  topology maxtree;
  copy_spdist_matrix_to_distance_matrix_upper (dist, pool->square_matrix, use_within_gf_means);
  maxtree = new_topology (sptaxa->nstrings);
  maxtree->taxlabel = sptaxa; sptaxa->ref_counter++; /* sptaxa is pointed here and at the calling function */
  if (tree_method == 0) bionj_from_distance_matrix (maxtree, pool->square_matrix); 
  else if (tree_method == 1) upgma_from_distance_matrix (maxtree, pool->square_matrix, false); // false -> upgma, true -> single linkage 
  else upgma_from_distance_matrix (maxtree, pool->square_matrix, true); // false -> upgma, true -> single linkage 
  add_topology_to_topology_space_if_distinct (maxtree, tsp, true); // true -> consider root location
  return;
}

/* genetree_struct functions */

genetree 
new_genetree (topology g_tree, topology s_tree)
{
  int i;
  genetree gt;
  gt = (genetree) biomcmc_malloc (sizeof (struct genetree_struct));
  gt->tree = g_tree; gt->tree->ref_counter++; // point to element of topology_space (usually)

  gt->split = create_splitset_dSPR_genespecies (gt->tree, s_tree);
  init_tree_recon_from_species_topology (gt->tree, s_tree);
  for (i = 0; i < NDISTS; i++) {
    gs->minmax[i] = 1.e35;
    gs->minmax[i+NDISTS] = -1.e35;
    gs->dist[i] = -1.;
  }
  del_char_vector (g_tree->taxlabel); // should just decrease ref_counter and then delete together with topology_space, later
  return gt;
}

void
del_genetree (genetree gt)
{
  if (!gt) return;
  del_splitset (split);
  del_topology (gt->tree);
  free (gt);
  return;
}

gene_sptrees
new_gene_sptrees (int n_genes, int n_ratchet, int n_proposal)
{
  gene_sptrees gs;
  gs = (gene_sptrees) biomcmc_malloc (sizeof (struct gene_sptrees_struct));
  gs->n_genes = n_genes;
  gs->n_ratchet = n_ratchet;
  gs->n_proposal = n_proposal;
  gs->next_proposal = 0; // idx of topol to fill proposal (and therefore to be modified)
  gs->next_free = 0; // idx of currently worse tree (which is just after best tree in a ratchet)
  gs->best_score = 1.e35;
  gs->gene = (genetree*) biomcmc_malloc (gs->n_genes * sizeof (genetree*)); 
  gs->ratchet  = (topology*) biomcmc_malloc (gs->n_ratchet * sizeof (topology*)); 
  gs->proposal = (topology*) biomcmc_malloc (gs->n_proposal * sizeof (topology*));
  // for i in ... malloc
  return gs;
}

void
del_gene_sptrees (gene_sptrees gs)
{
  int i;
  if (!gs) return;
  for (i = gs->n_genes - 1; i >= 0 ; i--) del_genetree (gs->gene);
  for (i = gs->n_ratchet- 1; i >= 0 ; i--) del_topology (gs->ratchet);
  for (i = gs->n_proposal - 1; i >= 0 ; i--) del_topology (gs->proposal);
  free (gs);
  return;
}



