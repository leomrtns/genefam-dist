#include <genefam_dist.h>  // timestamp 2016.09.27 

/* very few checks, since it's just for debugging module functions */
int
main (int argc, char **argv)
{
  clock_t time0, time1;
  char *s1 = "(((((T1:0.05992300866519096,(T2:0.0006180303152917179,T3:0.0006180303152917179):0.05930497834989924):0.043272898719766625,T4:0.10319590738495758):1.1239760857488426,(T5:0.4980243170498337,T6:0.4980243170498337):0.7291476760839662):1.281009272423232,T7:2.508181265557032):0.8827018107597968,T8:3.390883076316829):0.4625856443738513;";
  char *s2 = "(((T1:0.07863280033561784,T2:0.07863280033561784):0.13962732945264802,T3:0.21826012978826587):1.489139947849041,((T4:0.6897841786566293,(T5:0.6791207380820032,T6:0.6791207380820032):0.010663440574626088):0.24517159995111354,(T7:0.054232722534953585,T8:0.054232722534953585):0.8807230560727893):0.772444299029564):3.1394598688780517; ((T1:0.2512197737294074,T2:0.2512197737294074):0.7860571952762823,(((T3:0.004744266649772175,T4:0.004744266649772175):0.911730146173057,(T5:0.01730670374674797,T6:0.01730670374674797):0.8991677090760812):0.0813022597405953,(T7:0.568812647177929,T8:0.5688126471779292):0.4289640253854954):0.03950029644226517):4.045752587089589;";
  char *g_bug = "(Am_my_Ppal_EFA75307.1,(EE_ap_Ttra_AMSG_07526,((Am_is_Fflu_FL6WO5102SIE8Q,(Ex_ja_Rame_EC795785.1,Ex_ja_Jbah_EC686471.1)),(((Sr_rh_Cten_Contig4696,Am_di_Mspa_9162_1),(EE_is_Pbil_8550_1,(Sr_st_pram_fgenesh1pg.Cscaffold82000007_OG5_128421,Sr_st_Sdic_530743653))),(Sr_st_Tspa_5514_1,(((Pl_gr_ppat_egw1.95.21.1_OG5_128421,(Pl_gr_Atri_548840672,(Pl_gr_atha_NP193403_OG5_128421,Pl_gr_rcom_30101.m000378_OG5_128421))),(((((Ex_eu_tbru_Tb427tmp.211.1610_OG5_128421,(Ex_eu_Adea_528231608,Ex_eu_lmaj_LmjF.35.4590_OG5_128421)),Ex_eu_Ndes_7845_1),Ex_he_Pcos_5195_1),(Pl_rh_Ccho_546301620,Sr_rh_Crep_16086_1)),((Ex_he_Ngru_XP_002676866.1,((Ex_fo_glae_GLP155211_OG5_128421,Ex_fo_Ssal_558600419),(Sr_rh_Aspa_15156_1,Sr_rh_Rfil_ETO26683))),((Sr_st_Bhom_CBK20999.2,(Sr_is_Vbra_9211_1,((Sr_ch_Cvel_343870877,((Sr_ap_pfal_PFI1370c_OG5_128421,Sr_ap_cmur_CMU001440_OG5_128421),Sr_ap_tgon_TGME49025550_OG5_128421)),(((Sr_di_Gfol_48521_1,Sr_di_Shan_38824_1),((Sr_di_Pbah_48298_1,Sr_di_Lpol_346250803),((Sr_di_Pbei_59132_1,(Sr_st_Selo_11536_1,Sr_rh_Sspa_3228_1)),((Sr_di_Pbah_10608_1,Sr_di_Lpol_346274605),Sr_di_Gcat_43423_1)))),Sr_pe_Perk_XP_002787562.1)))),(Sr_ci_Slem_Contig5197.g5570,(Sr_ci_tthe_152.m00116_OG5_128421,Sr_ci_tthe_152.m00114_OG5_128421)))))),(((((((Op_me_Mlei_FC476215.1,Op_me_Ppil_295244175),(Op_me_nvec_estExtGenewiseH1.C1500035_OG5_128421,Op_me_nvec_egw.3160.3.1_OG5_128421)),((Op_fu_Npar_EIJ93167,Op_me_Aque_340373393),(((Op_fu_lbic_egwh1.11.122.1_OG5_128421,(((Op_fu_ncra_NCU03695T0_OG5_128421,Op_fu_afum_Afu1g15760_OG5_128421),Op_fu_scer_scers288cYNL169C_OG5_128421),(Op_fu_Bden_EGF83384.1,(Op_fu_spom_NP594463_OG5_128421,Op_fu_spom_NP595799_OG5_128421)))),(Op_fu_Mmel_328863906,Op_fu_Mglo_164662659)),Op_fu_Rirr_552918109))),((Op_me_Ctel_443694529,((Op_me_hsap_ENSP00000391739_OG5_128421,Op_me_cint_ENSCINP00000000439_OG5_128421),(Op_me_Skow_291229248,Op_me_Bflo_260833724))),((Pl_rh_Paer_23069_1,Op_me_cele_WBGene00015159_OG5_128421),(Op_me_sman_Smp021830_OG5_128421,Op_me_dmel_FBpp0083921_OG5_128421)))),((Op_ic_Cowc_EFW41470.1,(Op_ch_Mova_DC471793.1,Op_ch_Sros_514651458)),(Op_me_Hvul_221124716,Op_me_tadh_egw1.3.1254.1_OG5_128421))),Pl_rh_Gsul_545712519),((((Am_th_Tqua_Contig4779_281_reads_1162_bases,Am_di_Vexi_6493_1),(Am_di_Vrob_14231_1,Am_is_Sste_FL6WO5102QZJU6)),((Am_is_Sram_Contig47237,Am_di_Pess_4213),Am_is_Fnol_11646)),(Am_my_ddis_DDBG0276503_OG5_128421,(Am_is_Vant_Contig_1582_3353_reads_3468_bases,(Am_di_Naes_40059_1,Am_di_Patl_4453_1)))))))))),Am_my_ddis_DDBG0292748_OG5_128421);";
  char *s_bug = "(((((((((((((((Sr_st_Cspb,Sr_st_Aana),Sr_st_Ngad),Sr_st_Rmar),((Sr_st_Esil,Sr_st_Csub),Sr_st_Croe)),Ex_he_Smar),Ex_he_Pcos),Ex_he_Ngru),((((Ex_eu_Adea,Ex_eu_lmaj),Ex_eu_tbru),Ex_eu_Ndes),(Ex_eu_Egra,Ex_eu_Egym))),((((Ex_ja_Secu,Ex_ja_Rame),Ex_ja_Jbah),Ex_ja_Jlib),((Sr_ci_Slem,Sr_ci_tthe),Sr_ci_Cvir))),(((((((Sr_di_Pbah,Sr_di_Lpol),Sr_di_Gcat),(Sr_di_Shan,Sr_di_Gfol)),(Sr_rh_Sspa,Sr_di_Pbei)),Sr_pe_Perk),(((Sr_ap_pfal,Sr_ap_tgon),Sr_ap_cmur),(Sr_ch_Cvel,Sr_is_Vbra))),(((Sr_rh_Aspa,Sr_rh_Rfil),Sr_rh_Cten),(Sr_rh_Gsph,Sr_rh_Crep)))),((((((Ex_pa_Tfoe,Ex_pa_Hmel),Ex_pa_tvag),(EE_br_Stet,Sr_st_Tspa)),((Op_fu_ecun,Op_fu_Ncer),Op_fu_Npar)),(((Am_ar_Enut,Am_ar_ehis),(Ex_fo_glae,Ex_fo_Ssal)),((Sr_st_pram,Sr_st_Sdic),Sr_st_Bhom))),((((Sr_st_Tfra,Sr_st_Ptri),(Sr_st_Espi,Sr_st_tpse)),Sr_st_Bpac),((Sr_st_Spus,Sr_st_Selo),Sr_st_Ppar)))),((((((((Pl_gr_rcom,Pl_gr_atha),Pl_gr_Atri),Pl_gr_ppat),Pl_gr_Corb),Pl_gr_Mver),(((EE_is_Tsub,Pl_gr_Cvar),Pl_gr_vcar),(Pl_gr_Npyr,Pl_gr_micr))),(((Op_nu_Falb,EE_is_Pbil),EE_ka_Rtru),((EE_cr_Ccur,EE_cr_Rlen),EE_cr_Hand))),(((((Pl_rh_Pyez,Pl_rh_Ccho),Pl_rh_Rmar),((Pl_rh_Paer,Pl_rh_Ccoe),Pl_rh_Rmac)),Pl_rh_Gsul),((((EE_ha_Pcar,EE_ha_Ehux),EE_ha_Pant),EE_ha_Cpol),(EE_ha_Crho,EE_ha_Pspa))))),(((((((Op_me_Acal,Op_me_Ctel),Op_me_sman),(Op_me_dmel,Op_me_cele)),((Op_me_Skow,Op_me_Bflo),(Op_me_cint,Op_me_hsap))),(((Op_me_Ocar,Op_me_Aque),(Op_me_Mlei,Op_me_Ppil)),((Op_me_Hvul,Op_me_nvec),Op_me_tadh))),((Op_ch_Mova,Op_ch_Sros),(Op_ic_Cowc,Op_ic_Sarc))),(((((Op_fu_Mglo,Op_fu_lbic),Op_fu_Mmel),((Op_fu_ncra,Op_fu_scer),Op_fu_spom)),((Op_fu_Pspa,Op_fu_Bden),(Op_fu_Rory,Op_fu_Rirr))),((((EE_is_Ctri,Pl_gl_Cglo),(Pl_gl_Gnos,Pl_gl_Cpad)),Am_ar_Mbal),((Ex_ma_Mjak,Ex_ma_Mcal),EE_ap_Ttra))))),((((((Op_fu_afum,Am_th_Tqua),Am_is_Sste),Am_is_Vant),((Am_va_Usch,Am_di_Vrob),Am_va_Clsp)),((Am_di_Vexi,Am_di_Mspa),(Am_di_Patl,Am_di_Naes))),((((Am_tu_Hver,Am_di_Odes),Am_is_Pmon),((Am_hi_Gfon,Am_is_Sram),Am_di_Acas)),((Am_is_Fflu,Am_my_Ppol),(Am_is_Fnol,Am_di_Pess))))),(Am_my_Ppal,Am_my_ddis))";
  char *s3, *s4;
  double *dists;
  int i, n;

  (void) argc; (void) argv;
  time0 = clock ();

  printf ("T1: genefam_module_treesignal_fromtrees() from pre-defined strings\n DISTS = ");
  n = genefam_module_treesignal_fromtrees (s1, s2, &dists);
  for (i=0; i < n; i++) printf ("%lf ", dists[i]);
  if (dists) free (dists);

  printf ("\n\nT1.2: genefam_module_treesignal_fromtrees() from pre-defined strings (real example, with trifurcate unrooted gene tree)\n DISTS = ");
  n = genefam_module_treesignal_fromtrees (g_bug, s_bug, &dists);
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

