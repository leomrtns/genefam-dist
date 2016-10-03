#include <genefam_dist.h>  // timestamp 2016.09.27 

/* very few checks, since it's just for debugging module functions */
int
main (int argc, char **argv)
{
  clock_t time0, time1;
  char *s1 = "(((((T1:0.05992300866519096,(T2:0.0006180303152917179,T3:0.0006180303152917179):0.05930497834989924):0.043272898719766625,T4:0.10319590738495758):1.1239760857488426,(T5:0.4980243170498337,T6:0.4980243170498337):0.7291476760839662):1.281009272423232,T7:2.508181265557032):0.8827018107597968,T8:3.390883076316829):0.4625856443738513;";
  char *s2 = "(((T1:0.07863280033561784,T2:0.07863280033561784):0.13962732945264802,T3:0.21826012978826587):1.489139947849041,((T4:0.6897841786566293,(T5:0.6791207380820032,T6:0.6791207380820032):0.010663440574626088):0.24517159995111354,(T7:0.054232722534953585,T8:0.054232722534953585):0.8807230560727893):0.772444299029564):3.1394598688780517; ((T1:0.2512197737294074,T2:0.2512197737294074):0.7860571952762823,(((T3:0.004744266649772175,T4:0.004744266649772175):0.911730146173057,(T5:0.01730670374674797,T6:0.01730670374674797):0.8991677090760812):0.0813022597405953,(T7:0.568812647177929,T8:0.5688126471779292):0.4289640253854954):0.03950029644226517):4.045752587089589;";
  char *s3, *s4;
  double *dists;
  int i, n;

  (void) argc; (void) argv;
  time0 = clock ();

  printf ("T1: genefam_module_treesignal_fromtrees() from pre-defined strings\n DISTS = ");
  n = genefam_module_treesignal_fromtrees (s1, s2, &dists);
  for (i=0; i < n; i++) printf ("%lf ", dists[i]);
  if (dists) free (dists);

  printf ("\n\nT2: genefam_module_generate_sprtrees(100 leaves, 40 iterations, 1 SPR)  [ = sptree]\n");
  s3 = genefam_module_generate_spr_trees (100, 40, 1);
  printf (" STRING t2 = %s \n", s3);

  printf ("\nT3: genefam_module_generate_sprtrees(60 leaves, 1 iteration, 1 SPR) [ = genetree]\n");
  s4 = genefam_module_generate_spr_trees (60, 1, 1);
  printf (" STRING t3 = %s \n", s4);

  printf ("\nT4: genefam_module_treesignal_fromtrees(genetree, sptree) [only first gene_tree is used]\n DISTS = ");
  n = genefam_module_treesignal_fromtrees (s4, s3, &dists);
  for (i=0; i < n; i++) printf ("%lf ", dists[i]);
  if (dists) free (dists);
  printf ("\n(%d distances)\n", n);

  printf ("\nT5: genefam_module_treesignal_fromtrees_pvalue(genetree, sptree, 100 replicates) [only first gene_tree is used]\n DISTS = ");
  n = genefam_module_treesignal_fromtrees_pvalue (s4, s3, 100, &dists);
  for (i=0; i < n; i++) printf ("%lf ", dists[i]);
  if (dists) free (dists);
  printf ("\n(%d distances)\nend of tests\n", n);

  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
}

