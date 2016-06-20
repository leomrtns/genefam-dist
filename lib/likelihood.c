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

#include "likelihood.h"

const int LikScaleFrequency = 20;

/*! \brief main function that calculates log(likelihood) for changed nodes (called by high-level functions */
void calculate_ln_likelihood_proposal (phylogeny phy, topology tre);

/* real calculation (posterior distribution, using data) */
/*! \brief ln(likelihood) of topology, updating all internal nodes */ 
void ln_likelihood_real (phylogeny phy, topology tre);
/*! \brief ln(likelihood) of topology, based on changed nodes by dynamically updating lk_vector */
void ln_likelihood_moved_branches_real (phylogeny phy, topology tre);
/*! \brief ln(likelihood) of topology, based on changed nodes by statically updating lk_vector (under calling 
 * function control) */
void ln_likelihood_moved_branches_at_lk_vector_real (phylogeny phy, topology tre, int idx);

/* pointer to real or dummy likelihood calculation (defined in likelihood.h as external) */
void (*ln_likelihood) (phylogeny phy, topology tre) = &ln_likelihood_real;
void (*ln_likelihood_moved_branches) (phylogeny phy, topology tre) = &ln_likelihood_moved_branches_real;
void (*ln_likelihood_moved_branches_at_lk_vector) (phylogeny phy, topology tre, int idx) = &ln_likelihood_moved_branches_at_lk_vector_real;

/* dummy likelihood calculations (prior distribution, always zero) */
void
ln_likelihood_dummy (phylogeny phy, topology tre) 
{ phy->lk_proposal = 0.; (void) tre; }
void
ln_likelihood_moved_branches_dummy (phylogeny phy, topology tre) 
{ phy->lk_proposal = 0.; (void) tre; }
void
ln_likelihood_moved_branches_at_lk_vector_dummy (phylogeny phy, topology tre, int idx) 
{ phy->lk_proposal = 0.; (void) tre; (void) idx; }

void
set_likelihood_to_prior (void)
{
  ln_likelihood = &ln_likelihood_dummy;
  ln_likelihood_moved_branches = &ln_likelihood_moved_branches_dummy;
  ln_likelihood_moved_branches_at_lk_vector = &ln_likelihood_moved_branches_at_lk_vector_dummy;
}

void
set_likelihood_to_posterior (void)
{
  ln_likelihood = &ln_likelihood_real;
  ln_likelihood_moved_branches = &ln_likelihood_moved_branches_real;
  ln_likelihood_moved_branches_at_lk_vector = &ln_likelihood_moved_branches_at_lk_vector_real;
}

void
ln_likelihood_real (phylogeny phy, topology tre)
{ /* current --> proposal (=current->next) --> current */
  int i, cat, pat, s1, s2;
  double LikSite, lkMax, lkl, lkr, *left, *right, sum_of_lnLk = 0.;
  bool scale;

  if (!tre->traversal_updated) update_topology_traversal (tre);

  /*  "collapse(2)" only avail in openMP3.0 -- it collapses two for loops into one */
#ifdef _OPENMP
#pragma omp parallel for shared(phy,tre) \
  private(pat,cat,i,s1,s2,lkl,lkr,left,right,LikSite,lkMax,scale) reduction (+:sum_of_lnLk)
#endif
  for (pat = 0; pat < phy->npat; pat++) { 
    for (cat = 0; cat < phy->model->nrates; cat++) {
      for (i = 0; i < tre->nleaves - 2; i++) { /* skip postorder[nleaves-2] which is root node */
        scale = !(tre->postorder[i]->level % LikScaleFrequency); 
        lkMax = 0.; /* maximum partial likelihood for this node/category/pattern */

        left  = phy->l[tre->postorder[i]->left->id ]->d_current->next->lk[cat][pat];
        right = phy->l[tre->postorder[i]->right->id]->d_current->next->lk[cat][pat];
          
        phy->l[tre->postorder[i]->id]->d_current->next->lnmax[cat][pat] = 
        (phy->l[tre->postorder[i]->left->id ]->d_current->next->lnmax[cat][pat] +
         phy->l[tre->postorder[i]->right->id]->d_current->next->lnmax[cat][pat]);

        for (s1 = 0; s1 < 4; s1++) {
          lkl = lkr = 0.0;
          for (s2 = 0; s2 < 4; s2++) {
            lkl += phy->model->Q[cat][s1][s2] * left[s2];
            lkr += phy->model->Q[cat][s1][s2] * right[s2];
          }
          LikSite = phy->l[tre->postorder[i]->id]->d_current->next->lk[cat][pat][s1] = lkr * lkl;
          if (scale && (LikSite > lkMax)) lkMax = LikSite; /* may slow down since compared many times, most of them "false" */
        }
          
        if (scale) {
          /* scale the partial likelihoods to avoid underflow: unlike Yang's suggestion (JMolEvol.2000.423) we scale
           * each pattern, while he suggested over all patterns/sites. Each rate category is treated independently. */
          if (lkMax <= 0.) biomcmc_error ("underflow: all partial likelihoods are <= 0.");

          /* reescale only a few times since it is computationally expensive. Note that 
           * left->split->n_ones >= right->split->n_ones always (by design of update_topology_traversal() ) */
          phy->l[tre->postorder[i]->id]->d_current->next->lnmax[cat][pat] += log (lkMax);
          for (s1 = 0; s1 < 4; s1++) phy->l[tre->postorder[i]->id]->d_current->next->lk[cat][pat][s1] /= lkMax;
        }

      } // for (i < tre->nleaves) 

      /* root node is superfluous: the site likelihood is calculated between root->left and root->right */
      LikSite = 0.;
      left  = phy->l[tre->root->left->id ]->d_current->next->lk[cat][pat];
      right = phy->l[tre->root->right->id]->d_current->next->lk[cat][pat];
      /* lkMax will contain the sum of all scaling factors in log scale */
      lkMax = (phy->l[tre->root->left->id ]->d_current->next->lnmax[cat][pat] + 
               phy->l[tre->root->right->id]->d_current->next->lnmax[cat][pat]);

      for (s1 = 0; s1 < 4; s1++) for (s2 = 0; s2 < 4; s2++) /* likelihood at root for pattern */
        LikSite += phy->model->pi[s1] * left[s1] * phy->model->Q[cat][s1][s2] * right[s2];

      /* log likelihood of pattern, averaged over discretized rates */
      if (!cat) phy->pat_lnLk[pat] = log (LikSite) + lkMax; /* logspace_add(A,B) = log(exp(A)+exp(B)) below */
      else      phy->pat_lnLk[pat] = biomcmc_logspace_add (phy->pat_lnLk[pat], log (LikSite) + lkMax);
    } // for (category)

    /* phylogenetic log likelihood over sites (weighted patterns), summed through parallel reduction */
    sum_of_lnLk += phy->pat_lnLk[pat] * phy->weight[pat];
  } // for (pattern) 

  /* log (phy->model->nrates) is irreleveant in MCMC since it is a constant. It's here for completeness */
  phy->lk_proposal = sum_of_lnLk - ((double) (phy->nsites) * log ((double) phy->model->nrates));
}

void 
accept_likelihood (phylogeny phy, topology tre)
{ /* doesn't need topology actually; it's here just for consistency with ohter equivalent functions */
  int i;
  /* when working with u_done preorder update should come here */

  phy->lk_current = phy->lk_proposal;

  for (i=tre->nleaves; i < tre->nnodes; i++) { /* topology and phylogeny share same ids */
    phy->l[i]->d_current = phy->l[i]->d_current->next; /* update lk ring */
  }
}	

void
ln_likelihood_moved_branches_real (phylogeny phy, topology tre)
{ /* current --> proposal (=current->next) --> current ; final acceptance from calling function */
  int i;

  if (!tre->traversal_updated) update_topology_traversal (tre);
  if ((!tre->n_undone) && tre->root->left->d_done && tre->root->right->d_done) return; 

  for (i = 0; i < tre->n_undone; i++) { /* scan all nodes with d_done = false, but updating only children */ 
    if (tre->undone[i]->left->d_done)
      phy->l[ tre->undone[i]->left->id  ]->d_proposal = phy->l[ tre->undone[i]->left->id  ]->d_current;
    else
      phy->l[ tre->undone[i]->left->id  ]->d_proposal = phy->l[ tre->undone[i]->left->id  ]->d_current->next;

    if (tre->undone[i]->right->d_done)
      phy->l[ tre->undone[i]->right->id  ]->d_proposal = phy->l[ tre->undone[i]->right->id  ]->d_current;
    else
      phy->l[ tre->undone[i]->right->id  ]->d_proposal = phy->l[ tre->undone[i]->right->id  ]->d_current->next;
  }

  calculate_ln_likelihood_proposal (phy, tre);
}

void 
accept_likelihood_moved_branches (phylogeny phy, topology tre)
{
  int i;
  /* when working with u_done preorder update should come here */

  phy->lk_current = phy->lk_proposal;

  for (i = 0; i < tre->n_undone; i++) { /* topology and phylogeny share same ids */
    phy->l[tre->undone[i]->id]->d_current = phy->l[tre->undone[i]->id]->d_proposal; /* update lk ring */
    tre->undone[i]->d_done = true; /* update tree d_done (since we still don't use u_done) */
  }
}


void
ln_likelihood_moved_branches_at_lk_vector_real (phylogeny phy, topology tre, int idx)
{ /* accepted --> proposal (by array index, not linked list) --> accepted */
  int i;

  if (!tre->traversal_updated) update_topology_traversal (tre);
  if (!idx) biomcmc_error ("proposal likelihood will overwrite accepted (not your fault, it's a bug)");

  for (i = 0; i < tre->n_undone; i++) { /* scan all nodes with d_done = false, but updating only children */ 
    if (tre->undone[i]->left->d_done)
      phy->l[ tre->undone[i]->left->id  ]->d_proposal = phy->l[ tre->undone[i]->left->id  ]->d[0];
    else
      phy->l[ tre->undone[i]->left->id  ]->d_proposal = phy->l[ tre->undone[i]->left->id  ]->d[idx];

    if (tre->undone[i]->right->d_done)
      phy->l[ tre->undone[i]->right->id  ]->d_proposal = phy->l[ tre->undone[i]->right->id  ]->d[0];
    else
      phy->l[ tre->undone[i]->right->id  ]->d_proposal = phy->l[ tre->undone[i]->right->id  ]->d[idx];
  }

  calculate_ln_likelihood_proposal (phy, tre);
}

void 
accept_likelihood_moved_branches_at_lk_vector (phylogeny phy, topology tre, int idx, double likelihood)
{
  int i;
  /* when working with u_done preorder update should come here */

  phy->lk_accepted = likelihood;

  for (i = 0; i < tre->n_undone - 1; i++) { /* topology and phylogeny share same ids */
    phy->l[tre->undone[i]->id]->d_accepted = phy->l[tre->undone[i]->id]->d_current = 
    phy->l[tre->undone[i]->id]->d[idx]; /* update lk ring */
    tre->undone[i]->d_done = true; /* update tree d_done (since we still don't use u_done) */
  }
}


void
calculate_ln_likelihood_proposal (phylogeny phy, topology tre)
{ 
  int i, cat, pat, s1, s2;
  double LikSite, lkMax, lkl, lkr, *left, *right, sum_of_lnLk = 0.;
  bool scale;

#ifdef _OPENMP
#pragma omp parallel for  shared(phy,tre) \
  private(pat,cat,i,s1,s2,lkl,lkr,left,right,LikSite,lkMax,scale) reduction (+:sum_of_lnLk)
#endif
  for (pat = 0; pat < phy->npat; pat++) {
    for (cat = 0; cat < phy->model->nrates; cat++) {
      for (i = 0; i < tre->n_undone - 1; i++) { /* only nodes nodes that changed minus the root (n_undone -1)  */
        scale = !(tre->undone[i]->level % LikScaleFrequency); /* crude choice (the best would be distance from leaves) */ 

        left  = phy->l[ tre->undone[i]->left->id  ]->d_proposal->lk[cat][pat];
        right = phy->l[ tre->undone[i]->right->id ]->d_proposal->lk[cat][pat];

        phy->l[ tre->undone[i]->id        ]->d_proposal->lnmax[cat][pat] = 
        phy->l[ tre->undone[i]->left->id  ]->d_proposal->lnmax[cat][pat] +
        phy->l[ tre->undone[i]->right->id ]->d_proposal->lnmax[cat][pat];

        /* reescale only a few times since it is computationally expensive. Note that 
         * left->split->n_ones >= right->split->n_ones always (by design of update_topology_traversal() ) */
        if (scale) { /* reescaling is necessary */ 
          lkMax = 0.; /* maximum partial likelihood for this node/category/pattern */
          for (s1 = 0; s1 < 4; s1++) {
            lkl = lkr = 0.0;
            for (s2 = 0; s2 < 4; s2++) {
              lkl += phy->model->Q[cat][s1][s2] * left[s2];
              lkr += phy->model->Q[cat][s1][s2] * right[s2];
            }
            LikSite = phy->l[tre->undone[i]->id]->d_proposal->lk[cat][pat][s1] = lkr * lkl;
            if (LikSite > lkMax) lkMax = LikSite;
          }

          if (lkMax <= 0.) biomcmc_error ("underflow: all partial likelihoods are <= 0.");
          phy->l[tre->undone[i]->id]->d_proposal->lnmax[cat][pat] += log (lkMax);
          for (s1 = 0; s1 < 4; s1++)
            phy->l[tre->undone[i]->id]->d_proposal->lk[cat][pat][s1] /= lkMax;
        } 

        else { /* do not reescale: economize four comparisons (search for maximum) */
          for (s1 = 0; s1 < 4; s1++) {
            lkl = lkr = 0.0;
            for (s2 = 0; s2 < 4; s2++) {
              lkl += phy->model->Q[cat][s1][s2] * left[s2];
              lkr += phy->model->Q[cat][s1][s2] * right[s2];
            }
            phy->l[tre->undone[i]->id]->d_proposal->lk[cat][pat][s1] = lkr * lkl;
          }
        }
      } // for (i < tre->n_undone) 

      /* root node is superfluous: the site likelihood is calculated between root->left and root->right.
       * By design the heavier node (more nodes) is on the left */
      LikSite = 0.;
      
      left  = phy->l[ tre->root->left->id  ]->d_proposal->lk[cat][pat];
      right = phy->l[ tre->root->right->id ]->d_proposal->lk[cat][pat];

      lkMax = (phy->l[ tre->root->left->id  ]->d_proposal->lnmax[cat][pat] + 
               phy->l[ tre->root->right->id ]->d_proposal->lnmax[cat][pat]);

      for (s1 = 0; s1 < 4; s1++) for (s2 = 0; s2 < 4; s2++)
        LikSite += phy->model->pi[s1] * left[s1] * phy->model->Q[cat][s1][s2] * right[s2];

      if (!cat) phy->pat_lnLk[pat] = log (LikSite) + lkMax; /* logspace_add(A,B) = log(exp(A)+exp(B)) below */
      else      phy->pat_lnLk[pat] = biomcmc_logspace_add (phy->pat_lnLk[pat], log (LikSite) + lkMax);
    } // for (category)

    sum_of_lnLk += phy->pat_lnLk[pat] * phy->weight[pat];
  } // for (pattern) 

  /* log (phy->model->nrates) is irreleveant in MCMC since it is a constant. It's here for completeness */
  phy->lk_proposal = sum_of_lnLk - ((double) (phy->nsites) * log ((double) phy->model->nrates));
}
