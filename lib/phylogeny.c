/* 
 * This file is part of guenomu, a hierarchical Bayesian procedure to estimate the distribution of species trees based
 * on multi-gene families data.
 * Copyright (C) 2009  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * Guenomu is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "phylogeny.h"

node_likelihood new_node_likelihood (int n_cat, int n_pat, int n_state, int n_cycle);
void            del_node_likelihood (node_likelihood l, int n_cat, int n_pat);
lk_vector new_lk_vector (int n_cat, int n_pat, int n_state);
void      del_lk_vector (lk_vector u, int n_cat, int n_pat);

void init_evolution_model_parameters (evolution_model m, double kappa, double alpha, double beta, double *pi);
void init_eigenvectors_from_eq_frequencies (double **z1, double **z2, double *pi);

phylogeny
new_phylogeny_from_alignment (alignment align, int n_cat, int n_state, int n_cycle, distance_matrix external_dist)
{
  int i, j, k, l;
  phylogeny phy;
  distance_matrix dist;
  double alpha, beta;

  if (!align->is_aligned) biomcmc_error ("can't build a phylogeny, sequences not aligned");

  if (external_dist == NULL) dist = new_distance_matrix_from_alignment (align);
  else dist = external_dist;

  if (dist->mean_K2P_dist < 1e-16) biomcmc_error ("average pairwise distance is too small");
  if (dist->var_K2P_dist  < 1e-16) biomcmc_error ("variance in pairwise distances is too small");
  beta  = dist->mean_K2P_dist / dist->var_K2P_dist;
  alpha = beta * dist->mean_K2P_dist;

  phy = new_phylogeny (align->ntax, n_cat, align->npat, n_state, n_cycle);
  init_evolution_model_parameters (phy->model, dist->mean_R, alpha, beta, dist->freq);
  phy->nsites = align->nchar; /* original number of sites (may be used to calculate some constant like gamma rates) */

  phy->align_filename = align->filename; /* inherit original file name information */
  align->filename = NULL;

  for (i = 0; i < phy->ntax; i++) { /* store "trivial likelihood" at first rate category lk[0] */
    store_likelihood_info_at_leaf (phy->l[i]->d[0]->lk[0], align->character->string[i], align->npat, n_state);
    for (j=1; j<n_cat; j++) for (k=0; k < phy->npat; k++) { 
      for (l=0; l < 4; l++) phy->l[i]->d[0]->lk[j][k][l] = phy->l[i]->d[0]->lk[0][k][l]; /* copy to other categories */
      phy->l[i]->d[0]->lnmax[j][k] = 0.; /* log (scale factor) is zero since tips are already scaled */
    }
  }

  for (i = 0; i < phy->npat; i++) phy->weight[i] = (double) align->pattern_freq[i];

  if (external_dist == NULL) del_distance_matrix (dist);

  return phy;
}

phylogeny
new_phylogeny (int n_tax, int n_cat, int n_pat, int n_state, int n_cycle)
{
  int i;
  phylogeny phy;

  /* n_tax, n_pat   => stored here
   * n_cat, n_state => stored in evolution_model (n_cat = nrates)
   * n_cycle        => stored in node_likelihood */
  phy = (phylogeny) biomcmc_malloc (sizeof (struct phylogeny_struct));
  phy->npat = n_pat;
  phy->ntax = n_tax;
  phy->nsites = phy->npat; /* temporary assumption that all patterns are distinct (fixed by calling function) */
  phy->nnodes = 2 * n_tax - 1; 
  phy->lk_current = phy->lk_proposal = phy->lk_accepted = 0.;
  phy->align_filename = NULL;

  phy->l = (node_likelihood*) biomcmc_malloc ((phy->nnodes) * sizeof (node_likelihood));
  phy->weight   = (double*) biomcmc_malloc (n_pat * sizeof (double)); /* frequency of pattern */
  phy->pat_lnLk = (double*) biomcmc_malloc (n_pat * sizeof (double)); /* log likelihood of pattern */

  phy->model = new_evolution_model (n_cat, n_state);

  /* internal nodes must have at least one extra partial likelihood vectors (for proposal state) */
  for (i = 0; i < n_tax; i++)  phy->l[i] = new_node_likelihood (n_cat, n_pat, n_state, 1); /* leaf */
  for (; i < phy->nnodes; i++) phy->l[i] = new_node_likelihood (n_cat, n_pat, n_state, n_cycle + 2); /* internal node */

  return phy;
}

void
del_phylogeny (phylogeny phy)
{
  int i;
  if (!phy) return;
  if (phy->weight)         free (phy->weight);
  if (phy->pat_lnLk)       free (phy->pat_lnLk);
  if (phy->align_filename) free (phy->align_filename);
  if (!phy->model) biomcmc_error ("I cannot deallocate phylogenetic memory since I lost the model");
  if (phy->l) {
    for (i = phy->nnodes - 1; i >= 0; i--) del_node_likelihood (phy->l[i], phy->model->nrates, phy->npat);
    free (phy->l);
  }
  del_evolution_model (phy->model);
  free (phy);
}

void
phylogeny_order_accepted_lk_vector (phylogeny phy)
{
  int i, j;
  lk_vector lk;

  for (i = phy->ntax; i < phy->nnodes; i++) { /* internal nodes only */
    lk = phy->l[i]->d_accepted; /* zeroes the ratchet (so that accepted->next->next = d[2] etc.) */
    for (j = 0; j < phy->l[i]->n_cycle; j++, lk = lk->next) phy->l[i]->d[j] = lk;
  }
}

void 
phylogeny_link_accepted_to_current (phylogeny phy)
{
  int i;
  for (i = phy->ntax; i < phy->nnodes; i++) phy->l[i]->d_accepted = phy->l[i]->d_current; /* update lk ring */
  phy->lk_accepted = phy->lk_current;
}	

void 
phylogeny_link_current_to_accepted (phylogeny phy)
{
  int i;
  for (i = phy->ntax; i < phy->nnodes; i++) phy->l[i]->d_current = phy->l[i]->d_accepted; /* update lk ring */
  phy->lk_current = phy->lk_accepted;
}	

node_likelihood
new_node_likelihood (int n_cat, int n_pat, int n_state, int n_cycle)
{
  int i;
  node_likelihood l;
  l =	(node_likelihood) biomcmc_malloc (sizeof (struct node_likelihood_struct));
  /* one vector for current state, one for proposal update and one for accepted */
  l->n_cycle = n_cycle;
  l->u = (lk_vector*) biomcmc_malloc ((l->n_cycle) * sizeof (lk_vector));
  l->d = (lk_vector*) biomcmc_malloc ((l->n_cycle) * sizeof (lk_vector));

  for (i=0; i < l->n_cycle; i++) { /* doubly-linked circular list */
    l->u[i] = new_lk_vector (n_cat, n_pat, n_state);
    l->d[i] = new_lk_vector (n_cat, n_pat, n_state);
  }

  l->u[0]->prev = l->u[l->n_cycle-1];
  l->d[0]->prev = l->d[l->n_cycle-1];
  l->u[l->n_cycle-1]->next = l->u[0]; 
  l->d[l->n_cycle-1]->next = l->d[0]; 
  for (i=0; i < l->n_cycle; i++) {
    if (!l->u[i]->prev) l->u[i]->prev = l->u[i-1];
    if (!l->d[i]->prev) l->d[i]->prev = l->d[i-1];
    if (!l->u[i]->next) l->u[i]->next = l->u[i+1];
    if (!l->d[i]->next) l->d[i]->next = l->d[i+1];
  }

  l->u_current = l->u_accepted = l->u[0];
  l->d_current = l->d_accepted = l->d[0];
  l->d_proposal = l->d[0];

  return l;
}

void
del_node_likelihood (node_likelihood l, int n_cat, int n_pat)
{
  int i;
  if (!l) return;
  if (l->u) {
    for (i = l->n_cycle - 1; i >= 0; i--) del_lk_vector (l->u[i], n_cat, n_pat);
    free (l->u);
  }
  if (l->d) {
    for (i = l->n_cycle - 1; i >= 0; i--) del_lk_vector (l->d[i], n_cat, n_pat);
    free (l->d);
  }
  free (l);
}

lk_vector
new_lk_vector (int n_cat, int n_pat, int n_state)
{
  int j, k;
  lk_vector u;

  u = (lk_vector) biomcmc_malloc (sizeof (struct lk_vector_struct));
  u->prev = u->next = NULL;

  u->lk    = (double***) biomcmc_malloc ((n_cat) * sizeof (double**));
  u->lnmax = (double**)  biomcmc_malloc ((n_cat) * sizeof (double*));
  for (j = 0; j < n_cat; j++) { 
    u->lk[j]    = (double**) biomcmc_malloc (n_pat * sizeof (double*));
    u->lnmax[j] = (double*)  biomcmc_malloc (n_pat * sizeof (double));
    for (k = 0; k < n_pat; k++) u->lk[j][k] = (double*) biomcmc_malloc (n_state * sizeof (double));
  }

  return u;
}

void
del_lk_vector (lk_vector u, int n_cat, int n_pat)
{
  if (!u) return;
  if (u->lk) {
    int j, k;
    for (j = n_cat - 1; j >= 0; j--) { 
      if (u->lk[j]) {
        for (k = n_pat - 1; k >= 0; k--) if (u->lk[j][k]) free (u->lk[j][k]);
        free (u->lk[j]);
      }
      if (u->lnmax[j]) free (u->lnmax[j]);
    }
    if (u->lk)    free (u->lk);
    if (u->lnmax) free (u->lnmax);
  }
  free (u);
}

evolution_model
new_evolution_model (int n_cat, int n_state) /*n_state MUST BE 4 */
{
  int i, j;
  evolution_model m;
  m = (evolution_model) biomcmc_malloc (sizeof (struct evolution_model_struct));
  m->nrates  = n_cat;
  m->n_state = n_state; 
  m->kappa = m->alpha = m->beta = 1.; /* arbitrary values */

  m->rate = (double*)   biomcmc_malloc (n_cat * sizeof (double));
  m->Q    = (double***) biomcmc_malloc (n_cat * sizeof (double**));
  for (i = 0; i < n_cat; i++) {
    m->Q[i]  = (double**) biomcmc_malloc (n_state * sizeof (double*));
    for (j = 0; j < n_state; j++)
      m->Q[i][j]  = (double*) biomcmc_malloc (n_state * sizeof (double));
  }

  m->pi  = (double*) biomcmc_malloc ((n_state + 2) * sizeof (double)); /* pi[4] = pi_Y; pi[5] = pi_R */
  for (i = 0; i < n_state; i++) m->pi[i] = 0.25;
  for (; i < n_state + 2; i++) m->pi[i] = 0.5; /* arbitrary values */

  m->psi = (double*) biomcmc_malloc (n_state * sizeof (double));

  m->z1 = (double**) biomcmc_malloc (n_state * sizeof (double*));
  m->z2 = (double**) biomcmc_malloc (n_state * sizeof (double*));
  for (i = 0; i < n_state; i++) {
    m->z1[i] = (double*) biomcmc_malloc (n_state * sizeof (double));
    m->z2[i] = (double*) biomcmc_malloc (n_state * sizeof (double));
  }

  return m;
}

void
del_evolution_model (evolution_model m)
{
  int i, j;

  if (!m) return;
  if (m->rate) free (m->rate);
  if (m->pi)   free (m->pi);
  if (m->psi)  free (m->psi);
  if (m->z1) { for (i = m->n_state - 1; i >= 0; i--) if (m->z1[i]) free (m->z1[i]); free (m->z1); }
  if (m->z2) { for (i = m->n_state - 1; i >= 0; i--) if (m->z2[i]) free (m->z2[i]); free (m->z2); }
  if (m->Q) {
    for (i = m->nrates - 1; i >= 0; i--) if (m->Q[i]) { 
      for (j = m->n_state - 1; j >= 0; j--) if (m->Q[i][j]) free (m->Q[i][j]);
      free (m->Q[i]);
    }
    free (m->Q);
  }
  free (m);
}

void
init_evolution_model_parameters (evolution_model m, double kappa, double alpha, double beta, double *pi)
{
  int i;

  for (i = 0; i < m->n_state; i++) m->pi[i] = pi[i];
  m->pi[4] = pi[1] + pi[3]; /* pi_Y = pi_C + pi_T */ /* ASSUMES n_state = 4 */
  m->pi[5] = pi[0] + pi[2]; /* pi_R = pi_A + pi_G */

  /* initialize left and right eigenvectors (just need to be done once since eq. freqs. don't change) */
  init_eigenvectors_from_eq_frequencies (m->z1, m->z2, m->pi);

  m->kappa = kappa;
  m->alpha = alpha;
  m->beta  = beta;
  /* update psi */
  update_model_eigenvalues_from_kappa (m, m->kappa);
  /* calculate rates for each category *//*FIXME: what if rates too low (maybe reescale?) */
  biomcmc_discrete_gamma (m->alpha, m->beta, m->rate, m->nrates);
  /* update Q matrix for each rate */
  update_Q_matrix_from_average_rate (m, m->rate);
}

void
copy_evolution_model (evolution_model to, evolution_model from, bool copy_Qmatrix)
{
  int i, j, k;
  to->kappa = from->kappa;
  to->alpha = from->alpha;
  to->beta  = from->beta;

  for (i = 0; i < from->n_state + 2; i++) to->pi[i] = from->pi[i];
  for (i = 0; i < from->nrates; i++) to->rate[i] = from->rate[i];

  if (copy_Qmatrix) for (i = 0; i < from->nrates; i++) 
    for (j = 0; j < from->n_state; j++) for (k = 0; k < from->n_state; k++) to->Q[i][j][k] = from->Q[i][j][k];

  for (j = 0; j < from->n_state; j++) {
    to->psi[j] = from->psi[j];
    for (k = 0; k < from->n_state; k++) {
      to->z1[j][k] = from->z1[j][k];
      to->z2[j][k] = from->z2[j][k];
    }
  }
}


void
init_eigenvectors_from_eq_frequencies (double **z1, double **z2, double *pi)
{
  /* pi[4] = pi_Y = pi_C + pi_T (0123 -> ACGT)
   * pi[5] = pi_R = pi_A + pi_G  */
  z1[0][0] =  pi[0];
  z1[0][1] =  pi[1];
  z1[0][2] =  pi[2];
  z1[0][3] =  pi[3];
  z1[1][0] = -pi[0] * pi[4];
  z1[1][1] =  pi[1] * pi[5];
  z1[1][2] = -pi[2] * pi[4];
  z1[1][3] =  pi[3] * pi[5];
  z1[3][0] =  z1[3][2] = z1[2][1] = z1[2][3] = 0.;
  z1[3][3] =  z1[2][0] = 1.;
  z1[3][1] =  z1[2][2] = -1.;

  z2[0][0] =  z2[0][1] = z2[0][2] = z2[0][3] = 1.;
  z2[1][0] = -1./pi[5];
  z2[1][1] =  1./pi[4];
  z2[1][2] = -1./pi[5];
  z2[1][3] =  1./pi[4];
  z2[3][0] =  z2[3][2] = z2[2][1] = z2[2][3] = 0.;
  z2[3][1] = -pi[3]/pi[4];
  z2[3][3] =  pi[1]/pi[4];
  z2[2][0] =  pi[2]/pi[5];
  z2[2][2] = -pi[0]/pi[5];
}

void
update_model_eigenvalues_from_kappa (evolution_model m, double kappa)
{  /* (double *psi, double *pi, double *kappa) */
  double k = kappa * ((m->pi[0]*m->pi[2]) + (m->pi[1]*m->pi[3])) + (m->pi[4]*m->pi[5]);
  //  k = (2. * k)/(2. + kappa);
  k = 0.5/k;
  m->psi[0] = 0.;
  m->psi[1] = k;
  m->psi[2] = ((kappa * m->pi[5]) + m->pi[4]) * k;
  m->psi[3] = ((kappa * m->pi[4]) + m->pi[5]) * k;
}

void
update_Q_matrix_from_average_rate (evolution_model m, double *lambda)
{  /* (double **Q, double **z1, double **z2, double *psi, double lambda) */
  int i, j, k, cat;
  for (cat = 0; cat < m->nrates; cat++) for (i=0; i < m->n_state; i++) for (j=0; j < m->n_state;j++) {
    m->Q[cat][j][i] = 0.;
    for (k=0; k < m->n_state; k++) m->Q[cat][j][i] += (m->z1[k][i] * m->z2[k][j])/(1. + (m->psi[k] * lambda[cat]));
  }
}


