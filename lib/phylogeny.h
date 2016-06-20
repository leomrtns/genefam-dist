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

/*! \file phylogeny.h
 *  \brief Phylogeny, defined here as topology plus likelihood vectors.
 */

#ifndef _biomcmc_phylogeny_h_
#define _biomcmc_phylogeny_h_

#include "alignment.h"
#include "prob_distribution.h"

typedef struct phylogeny_struct* phylogeny;
typedef struct evolution_model_struct* evolution_model;
typedef struct node_likelihood_struct* node_likelihood;
typedef struct lk_vector_struct* lk_vector;

/*! \brief Model parameters and likelihood vectors for one segment. */
struct phylogeny_struct
{
  node_likelihood *l; /*! \brief Likelihood structure (circular list of vectors) for each node. */
  int npat;   /*! \brief Sequence size (in number of distinct patterns, not sites). */
  int nsites; /*! \brief Sequence original size, in sites (but in likelihood we use basically npat) */
  int ntax;   /*! \brief Number of taxa. */
  int nnodes;	/*! \brief Number of nodes (internal and leaves). */
  double *weight; /*! \brief Frequency of each site pattern (same AGCT pattern). */
  evolution_model model; /*! \brief Evolutionary model parameters */
  double lk_current;	/*! \brief Current \f$ ln(L) \f$. Equivalent to likelihood_accepted within minisample */
  double lk_proposal;	/*! \brief Proposal \f$ ln(L) \f$. Ultimately subject to acceptance/rejection by MCMC.*/
  double lk_accepted;	/*! \brief Accepted \f$ ln(L) \f$. */
  double *pat_lnLk;   /*! \brief sitewise (pattern-wise, in fact) log of likelihood, marginalized over rates */
  char *align_filename;  /*! \brief name of original alignment file, without extension */ 
};

/* model parameters assuming integrated rate and discrete gamma */
struct evolution_model_struct
{
  double *rate,	 /*! \brief expected substitution rate (one for each gamma category) */
         ***Q,   /*! \brief Transition probability matrix (one 4x4 vector for each category) */
         kappa,  /*! \brief transition/transversion ratio \f$\kappa_i\f$ for HKY model */
         *pi,    /*! \brief Equilibrium base distribution */
         **z1,   /*! \brief Left eigenvector for HKY model (depends on pi[]) */
         **z2,   /*! \brief Right eigenvector for HKY model (depends on pi[] */ 
         *psi;   /*! \brief Eigenvalues for HKY model (function of kappa) */
  double alpha,  /*! \brief alpha from the discrete gamma (sitewise heterogeneity) E[x]=alpha/beta */
         beta;   /*! \brief beta from the discrete gamma (sitewise heterogeneity) */
  int nrates,    /*! \brief number of discrete rate categories */
      n_state;   /*! \brief number of states (4 for DNA, 64 for codon...) MUST BE 4 currently */
};

/*! \brief Partial Likelihood information for each node such that no calculation is necessary 
 * if new state is rejected. */
struct node_likelihood_struct
{
  lk_vector *u; /*! \brief Upstream partial likelihood linked list. */
  lk_vector *d; /*! \brief Dounstream partial likelihood linked list. */
  int n_cycle;  /*! \brief Likelihood linked list size (1 for leaves and minisampler size for internal nodes). */
  /*! \brief Upstream partial likelihood of current topology (inside Al-Awadi update, for instance). Proposal topologies are 
   * accessed through node_likelihood_struct::u_current->next. */
  lk_vector u_current;
  /*! \brief downstream partial likelihood of current topology (inside Al-Awadi update, for instance). Proposal topologies are 
   * accessed through node_likelihood_struct::d_current->next. */
  lk_vector d_current;
  lk_vector d_proposal; /*!< \brief precalculated proposal vector element */
  lk_vector u_accepted, /*!< \brief Upstream partial likelihood of last accepted topology (before update). */
            d_accepted; /*!< \brief Downstream partial likelihood of last accepted topology (before update). */
};

/*! \brief Circular linked list with partial likelihood information for a node. Its size is the largest between 
 * chain_data_struct::n_cycles and chain_data_struct::n_mini. */
struct lk_vector_struct
{
  double ***lk;   /*! \brief Partial likelihood values for each pattern,  gamma category and state (A,G,C,T). */
  double **lnmax; /*! \brief scaling factors following Yang's JMolEvol.2000.423 to avoid underflow */
  lk_vector next, prev; /*! \brief Double-linked circular list information */
};


phylogeny new_phylogeny_from_alignment (alignment align, int n_cat, int n_state, int n_cycle, distance_matrix external_dist);

phylogeny new_phylogeny (int n_tax, int n_cat, int n_pat, int n_state, int n_cycle);

void del_phylogeny (phylogeny phy);

/* "rewinds" the circular liknked list of partial likelihoods such that we access the proposal likelihoods through d[]
 * instead of d->next->...->next. IOW makes the equivalence between a linked list and an array */
void phylogeny_order_accepted_lk_vector (phylogeny phy);

/* update lk_vector elements by making d_accepted point to d_current (as when new state is definitely accepted) */
void phylogeny_link_accepted_to_current (phylogeny phy);

/* update lk_vector elements by making d_current point to d_accepted (used before complex proposals since ln_likelihood
 * functions generally work on d_current */
void phylogeny_link_current_to_accepted (phylogeny phy);

/*! \brief Phylogenetic evolutionary model parameters (for likelihood calculation) */
evolution_model new_evolution_model (int n_cat, int n_state);
void del_evolution_model (evolution_model m);

/*! \brief copy values from one evolution_model to another, possibly skipping the transition matrix */
void copy_evolution_model (evolution_model to, evolution_model from, bool copy_Qmatrix);

void update_model_eigenvalues_from_kappa (evolution_model m, double kappa);

/*! \brief  HKY model integrated over branch length FIXME: change to E[x]=1/lambda (redo calcs) 
 *
 *  Calculates \f$ Prob(j/i,\lambda) \f$ where \f$ i,j\in {A,C,G,T}\f$ and branch length is assumed to have 
 *  exponential distribution with mean \f${\lambda}\f$, leading to
 *  \f[ Prob(j/i,\lambda)=\int Prob(j/i,t)p(t/\lambda)\, \mathrm{d}t \f]
 *  Assuming
 *  \f[ p(t/\lambda)\mathrm{d}t=\frac1{\lambda} e^{\frac{t}{\lambda}}\mathrm{d}t \f]
 *  and the HKY model \f$ Prob(j/i)\f$ reduces to
 *  \f[ Prob(j/i,\lambda)=\sum_{i=1}^4\frac{Z_i Z_i^{-1}}{1 - \psi_i \lambda} \f]
 *  where \f$Z_i \f$ is the matrix of eigenvectors and \f$\psi_i \f$ are the eigenvalues for the HKY model. */
void update_Q_matrix_from_average_rate (evolution_model m, double *lambda);

#endif
