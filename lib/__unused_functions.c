
/* from topology_space.c - unused since we must reorder all leaf numbers; easier to use external_hash when *reading* tree file */
void 
reorder_topol_space_leaves_from_hash (topology_space tsp, hashtable external_hash)
{
  int *order, i;
  order = (int*) biomcmc_malloc (2 * tsp->taxlabel->nstrings * sizeof (int)); /* two vecs, second is temporary */
  for (i=0; i < tsp->taxlabel->nstrings; i++) {
    /* map order in which taxlabels appear originally - where hashtable came from, e.g. the alignment file */
    order[i] = lookup_hashtable (external_hash, tsp->taxlabel->string[i]);
    if (order[i] < 0) {
      del_topology_space (tsp);
      biomcmc_error ( "tree label %s not found in global hashtable\n", tsp->taxlabel->string[i]); 
    }
  }
  // STOPHERE lines below (tree is nexus) must be replaced by changing all topol pointers 
  //original_order = (*order) + tsp->taxlabel->nstrings;
  //for (i=0; i < tree->nleaves; i++) original_order[i] = tree->leaflist[i]->id;
  //for (i=0; i < tree->nleaves; i++) tree->leaflist[i]->id = (*order)[ original_order[i] ];

  del_hashtable (tsp->taxlabel_hash);
  tsp->taxlabel_hash = external_hash;
  external_taxhash->ref_counter++; /* since we are sharing the hashfunction */
  // reorder taxlabels to conform to hashtable 
  char_vector_reorder_strings (tsp->taxlabel, order);
}
