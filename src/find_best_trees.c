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
  double minmax[NDISTS], dist[NDISTS]; // values are integers; might have a variable "scale" to allow for quick or no scaling 
  splitset split;
};

struct gene_sptrees_struct
{
  int n_genes, n_ratchet, n_proposal;
  int next_avail; // idx to ratchet elements 
  genetree *gene;
  topology *ratchet, *proposal; // may need vector of distances
  double *ratchet_score, best_score; // I'm forgetting other vars like minmax etc. but for now we will use simple score (gene dists are rescaled)
};

void print_usage (char *progname);
topology_space maxtrees_from_subsamples (char_vector sptaxa, char **arg_filename, int arg_nfiles, char_vector gfilename);

pooled_matrix new_pooled_matrix (int n_sets, int n_species);
void del_pooled_matrix (pooled_matrix pool);
void update_pooled_matrix_from_gene_tree (pooled_matrix pool, topology gene_topol, int *idx_gene_to_species);
void finalise_pooled_matrix (pooled_matrix pool, char_vector sptaxa);
void find_maxtree_and_add_to_topol_space (spdist_matrix dist, pooled_matrix pool, topology_space tsp, int tree_method, bool use_within_gf_means);

genetree new_genetree (topology g_tree, topology s_tree);
double calculate_genetree_distance (genetree gt, topology s_tree, bool update_minmax);
void del_genetree (genetree gt);

gene_sptrees new_gene_sptrees (int n_genes, int n_ratchet, int n_proposal, topology_space sptree);
topology init_topol_for_new_gene_sptrees (topology original);
void initialise_ratchet_from_topol_space (gene_sptrees gs, topology_space sptree);
void initialise_gene_sptrees_with_genetree_filenames (gene_sptrees gs, char_vector gene_files);
void del_gene_sptrees (gene_sptrees gs);
void initial_sorting_of_gene_sptrees_ratchet (gene_sptrees gs, topology_space tsp);
double calculate_species_tree_score();
void improve_gene_sptrees_ratchet (gene_sptrees gs, int iteration);
void generate_new_gene_sptrees_proposal_trees (gene_sptrees gs, int idx, int iteration);
void add_topol_to_gene_sptrees_ratchet (gene_sptrees gs, int idx_proposal, double score);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i,j;
  char_vector species, genefiles;
  topology_space tsp;
  gene_sptrees gs;
  FILE *stream;
  char *s;

  time0 = clock ();
  if (argc < 3) print_usage (argv[0]);
  biomcmc_random_number_init(0);

  species = new_char_vector_from_file (argv[1]);
  genefiles = new_char_vector (1); // updated by maxtrees_from_subsamples()
  char_vector_remove_duplicate_strings (species); /* duplicate names */
  char_vector_longer_first_order (species); // order species names from longest to shortest (speed up gene/spnames comparison) 
  tsp = maxtrees_from_subsamples (species, argv + 2, argc - 2, genefiles); 
  save_topology_space_to_trprobs_file (tsp, "patristic.tre", 1.);

  i = (int)(genefiles->nstrings/100); 
  if (i < 10) i = 10;
  gs = new_gene_sptrees (i, 100, 50, tsp); // i = n_genes, n_ratchet, n_proposal
  initialise_gene_sptrees_with_genetree_filenames (gs, genefiles);
  initial_sorting_of_gene_sptrees_ratchet (gs, tsp);
  for (i=0; i < 10; i++) {
    improve_gene_sptrees_ratchet (gs, i);
    j = gs->next_avail+1; if (j == gs->n_ratchet) j = 0;
    s = topology_to_string_by_name (gs->ratchet[j], NULL); /* second parameter is vector with branch lengths */
    printf ("[score %lf] %s;\n",gs->ratchet_score[j],s); fflush(stdout); free (s);
  }

  stream = biomcmc_fopen ("best_trees.tre", "w");
  fprintf (stream, "#NEXUS\nBegin trees;\n");
  for (i = 0; i < gs->n_ratchet; i++) {
    s = topology_to_string_by_name (gs->ratchet[i], NULL); 
    fprintf (stream, "tree PAUP_%d = %s;\n", i, s); free (s);
  }
  fprintf (stream, "End;\n");
  fclose (stream);

  del_char_vector (species);
  del_char_vector (genefiles);
  del_topology_space (tsp);
  del_gene_sptrees (gs);
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

  pool = new_pooled_matrix ((int)(arg_nfiles/10), sptaxa->nstrings);

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
    find_maxtree_and_add_to_topol_space (pool->d_total, pool, tsp, i, true);  // mean length within gene family (locus)
    find_maxtree_and_add_to_topol_space (pool->d_total, pool, tsp, i, false); // min lengths within gene family
  }
  for (i = 0; i < pool->n_sptrees; i++) {
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, 0, true); // bionj/upgma/singlelinkage, means/mins
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, 1, false); // upgma on mins 
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, 2, false); // single linkage on mins 
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
  pool->n_sets_per_gene = n_sets/20;
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
find_maxtree_and_add_to_topol_space (spdist_matrix dist, pooled_matrix pool, topology_space tsp, int tree_method, bool use_within_gf_means)
{
  topology maxtree;
  copy_spdist_matrix_to_distance_matrix_upper (dist, pool->square_matrix, use_within_gf_means);
  maxtree = new_topology (tsp->taxlabel->nstrings);
  maxtree->taxlabel = tsp->taxlabel; tsp->taxlabel->ref_counter++; /* sptaxa is pointed here and at the calling function */
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
  // init_tree_recon_from_species_topology() would order sptree names everytime (and already calculates distances...)
  if (!gt->tree->rec) gt->tree->rec = new_reconciliation (g_tree->nleaves, s_tree->nleaves); // if gene has rec, it better be from same sptree size!
  index_sptaxa_to_reconciliation (s_tree->taxlabel, g_tree->taxlabel, gt->tree->rec);
  for (i = 0; i < NDISTS; i++) {
    gt->minmax[i] = -1.e35;  // this is NOT minmax, but only MAX (currently we dont use min, only max)
    gt->minmax[i+NDISTS] = -1.e35;
    gt->dist[i] = -1.;
  }
  del_char_vector (g_tree->taxlabel); // should just decrease ref_counter and then delete later, with topology_space
  return gt;
}

double
calculate_genetree_distance (genetree gt, topology s_tree, bool update_minmax)
{
  int k;
  double g_score = 0.;
  gene_tree_reconcile (gt->tree, s_tree);
//  dSPR_gene_species_rf (gt->tree, s_tree, gt->split); // does not compute hdist and dSPR, only RF
  dSPR_gene_species (gt->tree, s_tree, gt->split); 
  gt->dist[0] = (double) gt->tree->rec->ndups;
  gt->dist[1] = (double) gt->tree->rec->nloss;
  gt->dist[2] = (double) gt->tree->rec->ndcos;
  gt->dist[3] = (double) gt->split->rf;
  gt->dist[4] = (double) gt->split->hdist;
  gt->dist[5] = (double) gt->split->spr + gt->split->spr_extra;  // order same as guenomu, NOT same as genefam_dist or treesignal.py
  if (update_minmax) for (k=0; k < NDISTS; k++) { // we update ONLY MAX currently
    if (gt->minmax[k] <  gt->dist[k]) gt->minmax[k] = gt->dist[k];
  }
  for (k=0; k < NDISTS; k++) g_score += biomcmc_log1p( (double)(gt->dist[k])/(double)(gt->minmax[k])); 
//  for (k=0; k < NDISTS; k++) g_score += biomcmc_log1p( (double)(gt->dist[k])); 
  return g_score;
}

void
del_genetree (genetree gt)
{
  if (!gt) return;
  del_splitset (gt->split);
  del_topology (gt->tree);
  free (gt);
  return;
}

gene_sptrees
new_gene_sptrees (int n_genes, int n_ratchet, int n_proposal, topology_space sptree)
{
  int i;
  gene_sptrees gs;
  gs = (gene_sptrees) biomcmc_malloc (sizeof (struct gene_sptrees_struct));
  gs->n_genes = n_genes;
  gs->n_ratchet = n_ratchet;
  gs->n_proposal = n_proposal;
  gs->best_score = 1.e35;
  if (gs->n_ratchet < 10) gs->n_ratchet = 10; // ratchet[0] has best score
  if (gs->n_proposal < 10) gs->n_proposal = 10;
  gs->next_avail = gs->n_ratchet - 1; // idx of currently worse tree (which is just before best tree in a ratchet)

  gs->gene = (genetree*) biomcmc_malloc (gs->n_genes * sizeof (genetree)); 
  gs->ratchet  = (topology*) biomcmc_malloc (gs->n_ratchet * sizeof (topology)); 
  gs->proposal = (topology*) biomcmc_malloc (gs->n_proposal * sizeof (topology));
  gs->ratchet_score  = (double*) biomcmc_malloc (gs->n_ratchet * sizeof (double)); 

  for (i=0; i < gs->n_ratchet; i++) {
    gs->ratchet[i] = init_topol_for_new_gene_sptrees (sptree->distinct[0]);
    randomize_topology (gs->ratchet[i]); // should be filled outside this function (to recycle optimal sptrees from previous iteration)
  }
  for (i=0; i < gs->n_proposal; i++) gs->proposal[i] = init_topol_for_new_gene_sptrees (sptree->distinct[0]);
  for (i=0; i < gs->n_genes; i++) gs->gene[i] = NULL; // must be created later, when we receive the gene topol_spaces 

  initialise_ratchet_from_topol_space (gs, sptree);
  return gs;
}

topology
init_topol_for_new_gene_sptrees (topology original)
{
  topology tree = new_topology(original->nleaves);
  tree->taxlabel = original->taxlabel;
  tree->taxlabel->ref_counter++;
  new_mrca_for_topology (tree);
  return tree;
}

void
initialise_ratchet_from_topol_space (gene_sptrees gs, topology_space sptree)
{
  int i, n_trees_to_copy = sptree->ndistinct;
  if (n_trees_to_copy > gs->n_ratchet) n_trees_to_copy = gs->n_ratchet;
  for (i=0; i < n_trees_to_copy; i++) copy_topology_from_topology (gs->ratchet[gs->n_ratchet - i -1], sptree->distinct[i]); 
  // the remaining trees (from n_trees_to_copy to n_ratchet) are random (if it's first time this function is called) 
  // or are best from previous iteration (first on ratchet) 
  return;
}

void
initialise_gene_sptrees_with_genetree_filenames (gene_sptrees gs, char_vector gene_files)
{ // sample gene families, and calculate max values 
  int i, j, *idx;
  topology_space genetre;

  idx = (int*) biomcmc_malloc (gs->n_genes * sizeof (int)); 
  for (i=0; i < gs->n_genes; i++) idx[i] = i;

  for (i=0; i < gs->n_genes; i++) {
    j = biomcmc_rng_unif_int (gs->n_genes - i);
    genetre = read_topology_space_from_file (gene_files->string[ idx[j] ], NULL);/* read i-th gene family */
    idx[j] = idx[gs->n_genes - i -1]; // avoids replacement
    gs->gene[i] = new_genetree (genetre->distinct[0], gs->proposal[0]);
    del_topology_space (genetre);
    for (j = 0; j < 32; j++) {
      randomize_topology (gs->proposal[0]);
      calculate_genetree_distance (gs->gene[i], gs->proposal[0], true);
    }
  }
  if (idx) free (idx);
  return;
}

void
del_gene_sptrees (gene_sptrees gs)
{
  int i;
  if (!gs) return;
  if (gs->gene) {
    for (i = gs->n_genes - 1; i >= 0 ; i--) del_genetree (gs->gene[i]);
    free (gs->gene);
  }
  if (gs->ratchet) {
    for (i = gs->n_ratchet- 1; i >= 0 ; i--) del_topology (gs->ratchet[i]);
    free (gs->ratchet);
  }
  if (gs->proposal) {
    for (i = gs->n_proposal - 1; i >= 0 ; i--) del_topology (gs->proposal[i]);
    free (gs->proposal);
  }
  if (gs->ratchet_score) free (gs->ratchet_score);
  free (gs);
  return;
}

void
initial_sorting_of_gene_sptrees_ratchet (gene_sptrees gs, topology_space tsp)
{ // first step before start optimisation: sort ratchet 
  topology *pivot; // temporary vector with original location of topologies
  empfreq_double ef;
  double *score_vec = NULL;
  int i;

  if ((!tsp) || (tsp->ndistinct <= gs->n_ratchet)) {
    for (i=0; i < gs->n_ratchet; i++) gs->ratchet_score[i] = calculate_species_tree_score (gs, gs->ratchet[i]);
    /* Find the score of each tree and sort their scores from lowest (best to highest (worst)  */ 
    ef = new_empfreq_double_sort_increasing (gs->ratchet_score, gs->n_ratchet);
    pivot = (topology*) biomcmc_malloc (gs->n_ratchet * sizeof (topology)); 
    /* reorder ratchet, by creating temporary pointers with previous location */
    for (i=0; i < gs->n_ratchet; i++) pivot[i] = gs->ratchet[i]; //vector version of tmp=a;a=b;b=tmp
    for (i=0; i < gs->n_ratchet; i++) gs->ratchet[i] = pivot[ ef->d[i].idx ]; 
    if (pivot) free (pivot);
  }
  else {
    score_vec  = (double*) biomcmc_malloc (tsp->ndistinct * sizeof (double)); 
    for (i=0; i < tsp->ndistinct; i++) score_vec[i] = calculate_species_tree_score (gs, tsp->distinct[i]);
    /* Find the score of each tree and sort their scores from lowest (best to highest (worst)  */ 
    ef = new_empfreq_double_sort_increasing (score_vec, tsp->ndistinct);
    for (i=0; i < gs->n_ratchet; i++) copy_topology_from_topology (gs->ratchet[i], tsp->distinct[ ef->d[i].idx ]);
    if (score_vec) free (score_vec);
  }

  for (i=0; i < gs->n_ratchet; i++) gs->ratchet_score[i] = ef->d[i].freq; // doesnt need pivot since ef->freqs came from ratchet_score[]
  gs->best_score = ef->d[0].freq; // lowest (best) score is first element (may need to be recalc every time new proposals are created, using new minima)
  for (i=0; i < gs->n_ratchet; i++) printf ("DEBUG::score:: %lf\n", gs->ratchet_score[i]);

  del_empfreq_double (ef);
  return;
}

double
calculate_species_tree_score (gene_sptrees gs, topology s_tree)
{
  int i;
  double score = 0.;
  for (i = 0; i < gs->n_genes; i++) score += calculate_genetree_distance (gs->gene[i], s_tree, false);
  return score;
} 

void
improve_gene_sptrees_ratchet (gene_sptrees gs, int iteration)
{ 
  int i,j;
  double score = 0;
  for (j = 0; j < gs->n_ratchet; j++) { 
    generate_new_gene_sptrees_proposal_trees (gs, j, iteration);
    for (i = 0; i < gs->n_proposal; i++) {
      score = calculate_species_tree_score (gs, gs->proposal[i]);
      if (score < gs->best_score) add_topol_to_gene_sptrees_ratchet (gs, i, score);
    }
  }
  return;
}

void
generate_new_gene_sptrees_proposal_trees (gene_sptrees gs, int idx, int iteration)
{
  int i, j, n_2 = gs->n_proposal/2;
  for (i = 0; i < gs->n_proposal; i++) copy_topology_from_topology (gs->proposal[i], gs->ratchet[idx]);
  if (!(iteration%4)) {
    for (i = 0; i < n_2; i++) topology_apply_rerooting (gs->proposal[i], false);
    for (i = n_2; i < gs->n_proposal; i++) topology_apply_shortspr (gs->proposal[i], false);
  }
  if (!((iteration+1)%4)) {
    for (i = 0; i < n_2; i++) topology_apply_spr_unrooted (gs->proposal[i], false);
    for (i = n_2; i < gs->n_proposal; i++) topology_apply_nni (gs->proposal[i], false);
  }
  if (!((iteration+2)%4)) {
    for (i = 0; i < gs->n_proposal; i++) topology_apply_spr_unrooted (gs->proposal[i], false);
    for (i = 0; i < n_2; i++) topology_apply_nni (gs->proposal[i], false);
  }
  if (!((iteration+3)%4)) {
    for (i = 0; i < gs->n_proposal; i++) { 
      topology_apply_shortspr (gs->proposal[i], false); 
      topology_apply_spr_unrooted (gs->proposal[i], false); 
    }
  }
  // if rooted trees are the same, randomise one of them (assuming random tree won't be the same anymore...)
  for (i = 0; i < gs->n_proposal; i++) for (j = 0; j < i; j++) if (topology_is_equal (gs->proposal[i], gs->proposal[j])) randomize_topology (gs->proposal[j]); 
  return;
}

void
add_topol_to_gene_sptrees_ratchet (gene_sptrees gs, int idx_proposal, double score)
{
  topology pivot = gs->proposal[idx_proposal];
  gs->proposal[idx_proposal] = gs->ratchet[gs->next_avail];
  gs->ratchet[gs->next_avail] = pivot;
  gs->ratchet_score[gs->next_avail] = score;

  printf ("DEBUG::bestfound::%d  \t%lf\n",gs->next_avail, score);
  gs->next_avail--;
  if (gs->next_avail < 0) gs->next_avail = gs->n_ratchet - 1;
  gs->best_score = score;
  return;
}
