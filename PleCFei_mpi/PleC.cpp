//**********************************************************
//Caulobacter crescentus cell cycle model
//
//For model details please refer to 
//
//A Stochastic Spatiotemporal Model of a Response-Regulator Network 
//in the Caulobacter crescentus Cell Cycle
//Fei Li , Kartik Subramanian, Minghan Chen, John J. Tyson, Yang Cao
//
//Physical Biology, special issue of 2015 q-bio 
//
//coded by: Fei Li
//felix@cs.vt.edu
//Dec 4, 2015
//***********************************************************

#include <cstring>
using std::strcpy; 

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <cstdlib>
using std::rand;

#include <cmath>
using std::log;

#include <fstream>
using std::ifstream;
using std::ofstream;
using std::ios;

#include "PleC.h"


PleC::PleC():
  none(),

  gene_ori(BINNUM), gene_dup(BINNUM),
  mRNA_DivJ(BINNUM), mRNA_DivK(BINNUM), mRNA_PleC(BINNUM),
  mRNA_DivL(BINNUM), mRNA_CckA(BINNUM), mRNA_CtrA(BINNUM),

  DivL(BINNUM), DivL_f(BINNUM), DivL_st(BINNUM),
  DivLKp(BINNUM),

  CckA(BINNUM), CckA_st(BINNUM),
  CckA_phos(BINNUM), CckA_kin(BINNUM),
  CtrA(BINNUM), CtrAp(BINNUM),

  DivJ(BINNUM), DivJ_f(BINNUM), DivJ_st(BINNUM),
  DivJK(BINNUM), DivJKp(BINNUM),

  DivK(BINNUM), DivKp(BINNUM),

  PleC_f(BINNUM), PleC_st(BINNUM),
  PleCphos(BINNUM), PleCph1(BINNUM), PleCph2(BINNUM),
  PleCkin11(BINNUM), PleCkin12(BINNUM), PleCkin22(BINNUM),
  PleCkin0(BINNUM), PleCkin2(BINNUM), PleCkin4(BINNUM),
  PleCkin1(BINNUM), PleCkin3(BINNUM),
  PleC2p(BINNUM), PleC1p(BINNUM),
  PleCkin02p(BINNUM), PleCkin10p(BINNUM), PleCkin01p(BINNUM),
  PleCpt4(BINNUM), PleCpt2(BINNUM), 

  /*DivJ mRNA reactions*/
  syn_OrimrDivJ(gene_ori, mRNA_DivJ, 0.625, syn),
  syn_DupmrDivJ(gene_dup, mRNA_DivJ, 0.625, syn),
  deg_mrDivJ(mRNA_DivJ, none, 0.25, deg),
  lftdiff_mrDivJ(mRNA_DivJ, none, 1.0, diffusion_lft, -1),
  rgtdiff_mrDivJ(mRNA_DivJ, none, 1.0, diffusion_rgt, -1),
  /*DivK mRNA reactions*/
  syn_OrimrDivK(gene_ori, mRNA_DivK, 0.625, syn),
  syn_DupmrDivK(gene_dup, mRNA_DivK, 0.625, syn),
  deg_mrDivK(mRNA_DivK, none, 0.25, deg),
  lftdiff_mrDivK(mRNA_DivK, none, 1.0, diffusion_lft, -1),
  rgtdiff_mrDivK(mRNA_DivK, none, 1.0, diffusion_rgt, -1),
  /*PleC mRNA reactions*/
  syn_OrimrPleC(gene_ori, mRNA_PleC, 0.625, syn),
  syn_DupmrPleC(gene_dup, mRNA_PleC, 0.625, syn),
  deg_mrPleC(mRNA_PleC, none, 0.25, deg),
  lftdiff_mrPleC(mRNA_PleC, none, 1.0, diffusion_lft, -1),
  rgtdiff_mrPleC(mRNA_PleC, none, 1.0, diffusion_rgt, -1),
  /*DivL mRNA reactions*/
  syn_OrimrDivL(gene_ori, mRNA_DivL, 0.625, syn),
  syn_DupmrDivL(gene_dup, mRNA_DivL, 0.625, syn),
  deg_mrDivL(mRNA_DivL, none, 0.25, deg),
  lftdiff_mrDivL(mRNA_DivL, none, 1.0, diffusion_lft, -1),
  rgtdiff_mrDivL(mRNA_DivL, none, 1.0, diffusion_rgt, -1),
  /*CckA mRNA reactions*/
  syn_OrimrCckA(gene_ori, mRNA_CckA, 0.625, syn),
  syn_DupmrCckA(gene_dup, mRNA_CckA, 0.625, syn),
  deg_mrCckA(mRNA_CckA, none, 0.25, deg),
  lftdiff_mrCckA(mRNA_CckA, none, 1.0, diffusion_lft, -1),
  rgtdiff_mrCckA(mRNA_CckA, none, 1.0, diffusion_rgt, -1),
  /*CtrA mRNA reactions*/
  syn_OrimrCtrA(gene_ori, mRNA_CtrA, 6.25, syn),
  syn_DupmrCtrA(gene_dup, mRNA_CtrA, 6.25, syn),
  deg_mrCtrA(mRNA_CtrA, none, 0.25, deg),
  lftdiff_mrCtrA(mRNA_CtrA, none, 1.0, diffusion_lft, -1),
  rgtdiff_mrCtrA(mRNA_CtrA, none, 1.0, diffusion_rgt, -1),

  /*DivK, DivKp*/
  syn_DivK(mRNA_DivK, DivK, 0.025*scale, syn),
  deg_DivK(DivK,   none, 0.005, deg),
  deg_DivKp(DivKp, none, 0.005, deg),
  lftdiff_DivK(DivK, none, 10.0, diffusion_lft, -1),
  rgtdiff_DivK(DivK, none, 10.0, diffusion_rgt, -1),
  lftdiff_DivKp(DivKp, none, 10.0, diffusion_lft, -1),
  rgtdiff_DivKp(DivKp, none, 10.0, diffusion_rgt, -1),
  /*PleC*/
  syn_PleCf(mRNA_PleC, PleC_f, 0.05*scale, syn),
  deg_PleCf(PleC_f,  none, 0.05, deg),
  lftdiff_PleC(PleC_f, none, 1.0, diffusion_lft, -1),
  rgtdiff_PleC(PleC_f, none, 1.0, diffusion_rgt, -1),
  bnd_PleC(PleC_f, PleCphos, PleC_st, 1.0, sticky),
  ubd_PleC(PleCphos, PleC_f, 0.5, first_order),
  /*PleCph1, PleCph2*/
  bndKp_PleCphos(PleCphos, DivKp, PleCph1, none, 5.0/scale, combination),
  bndK_PleCphos(PleCphos, DivK, PleCph2, none, 0.05/scale, combination),

  ubdKp_PleCph1(PleCph1, PleCphos, DivKp, 5.0, decomposition),
  ubdK_PleCph2(PleCph2,  PleCphos, DivK, 5.0, decomposition),

  phos_PleCph2(PleCph2, PleCph1, 0.005, first_order),
  deph_PleCph1(PleCph1, PleCph2, 10.0, first_order),

  disp_PleCph1(PleCph1, DivK,  PleCph2, DivKp, 0.016/scale, displacement),
  disp_PleCph2(PleCph2, DivKp, PleCph1, DivK,  1.6/scale, displacement),
  /*PleCkin11, PleCkin12, PleCkin22*/
  bndKp_PleCph1(PleCph1, DivKp, PleCkin11, none, 5.0/scale, combination),
  bndK_PleCph1(PleCph1, DivK, PleCkin12,  none, 0.016/scale, combination),
  bndKp_PleCph2(PleCph2, DivKp, PleCkin12, none, 1.6/scale, combination),
  bndK_PleCph2(PleCph2, DivK, PleCkin22, none, 0.016/scale, combination),

  ubdKp_PleCkin11(PleCkin11, PleCph1, DivKp, 2.5, decomposition),
  ubdKp_PleCkin12(PleCkin12, PleCph2, DivKp, 1.6e-4,decomposition),
  ubdK_PleCkin12(PleCkin12,  PleCph1, DivK,  1.6e-4, decomposition),
  ubdK_PleCkin22(PleCkin22,  PleCph2, DivK,  1.6e-8, decomposition),
  /*PleCkin0, PleCkin2, PleCkin4*/
  phos_PleCkin11(PleCkin11, PleCkin0, 2.5, first_order),
  phos_PleCkin12(PleCkin12, PleCkin2, 5.0, first_order),
  phos_PleCkin22(PleCkin22, PleCkin4, 5.0, first_order),

  deph_PleCkin0(PleCkin0, PleCkin11, 5.0, first_order),
  deph_PleCkin2(PleCkin2, PleCkin12, 5.0, first_order),
  deph_PleCkin4(PleCkin4, PleCkin22, 5.0, first_order),
  /*PleCkin1, PleCkin3*/
  ubdKp_PleCkin0(PleCkin0, PleCkin1, DivKp, 0.16, decomposition),
  ubdKp_PleCkin2(PleCkin2, PleCkin3, DivKp, 0.16, decomposition),
  ubdK_PleCkin2(PleCkin2,  PleCkin1, DivK,  1.6e-3, decomposition),
  ubdK_PleCkin4(PleCkin4,  PleCkin3, DivK,  1.6e-3, decomposition),

  bndKp_PleCkin1(PleCkin1, DivKp, PleCkin0, none, 5.0/scale, combination),
  bndK_PleCkin1(PleCkin1, DivK, PleCkin2,   none, 5.0/scale, combination),
  bndKp_PleCkin3(PleCkin3, DivKp, PleCkin2, none, 5.0/scale, combination),
  bndK_PleCkin3(PleCkin3, DivK, PleCkin4,   none, 5.0/scale, combination),
  /*PleC2p*/
  ubdKp_PleCkin1(PleCkin1, PleC2p, DivKp, 0.16, decomposition),
  ubdK_PleCkin3(PleCkin3,  PleC2p, DivK,  1.6e-3, decomposition),

  bndKp_PleC2p(PleC2p, DivKp, PleCkin1, none, 5.0/scale, combination),
  bndK_PleC2p(PleC2p, DivK, PleCkin3,   none, 5.0/scale, combination),

  deph_PleC2p(PleC2p, PleCphos, 5.0, first_order),

  /*PleC1p*/
  ubdKp_PleCkin10p(PleCkin10p, PleC1p, DivKp, 0.16, decomposition),
  ubdKp_PleCkin01p(PleCkin01p, PleC1p, DivKp, 0.16, decomposition),
  ubdK_PleCkin02p(PleCkin02p,  PleC1p, DivK,  0.0016, decomposition),

  bndKp_PleC1pt(PleC1p, DivKp, PleCkin10p, none, 5.0/scale, combination),
  bndKp_PleC1ps(PleC1p, DivKp, PleCkin01p, none, 5.0/scale, combination),
  bndK_PleC1ps(PleC1p, DivK, PleCkin02p,   none, 5.0/scale, combination),
 
  deph_PleC1p(PleC1p, PleCphos, 5.0, first_order),

  /*auto-phosphation*/
  autoph_PleCkin3(PleCkin3, PleCkin10p, 5.0, first_order),
  autoph_PleCkin2(PleCkin2, PleCpt2,    5.0, first_order),
  autoph_PleCkin4(PleCkin4, PleCpt4,    5.0, first_order),
  autoph_PleCpt4(PleCpt4,   PleCkin11,  5.0, first_order),

  deauto_PleCkin10p(PleCkin10p, PleCkin3, 0.16, first_order),
  deauto_PleCpt2(PleCpt2,  PleCkin2,  0.16, first_order),
  deauto_PleCpt4(PleCpt4,  PleCkin4,  0.16, first_order),
  deauto_PleCkin11(PleCkin11, PleCpt4,0.0755, first_order),

  /*PleCKin10p, PleCKin01p, PleCkin20p, PleCkin02p*/
  deph_PleCkin1t(PleCkin1, PleCkin01p, 5.0, first_order),  
  deph_PleCkin3t(PleCkin3, PleCkin02p, 5.0, first_order),  

  phos_PleCkin01p(PleCkin01p, PleCkin1, 5.0, first_order),  
  phos_PleCkin02p(PleCkin02p, PleCkin3, 5.0, first_order),  

  /*PleCpt2, PleCpt4*/
  bndKp_PleCkin01p(PleCkin01p, DivKp, PleCpt2, none, 5.0/scale, combination),
  bndKp_PleCkin02p(PleCkin02p, DivKp, PleCpt4, none, 5.0/scale, combination),

  ubdKp_PleCpt2s(PleCpt2, PleCkin01p, DivKp, 0.16, decomposition),
  ubdKp_PleCpt4(PleCpt4,  PleCkin02p, DivKp, 0.16, decomposition),
  /*PleC complex degradation*/
  deg_PleCphos(PleCphos, none, 0.05, deg),
  deg_PleCph1(PleCph1, DivKp,  0.05, first_order),
  deg_PleCph2(PleCph2,   none, 0.05, deg),
  deg_PleCkin11(PleCkin11, DivKp, DivKp, 0.05, decomposition),
  deg_PleCkin12(PleCkin12, DivKp, 0.05, first_order),
  deg_PleCkin22(PleCkin22,  none, 0.05, deg),
  deg_PleCkin0(PleCkin0, DivKp, DivKp, 0.05, decomposition),
  deg_PleCkin2(PleCkin2, DivKp, 0.05, first_order),
  deg_PleCkin4(PleCkin4, none, 0.05, deg),
  deg_PleCkin1(PleCkin1, DivKp, 0.05, first_order),
  deg_PleCkin3(PleCkin3,  none, 0.05, deg),
  deg_PleC2p(PleC2p, none, 0.05, deg),
  deg_PleC1p(PleC1p, none, 0.05, deg),
  deg_PleCkin02p(PleCkin02p,  none, 0.05, deg),
  deg_PleCkin10p(PleCkin10p, DivKp, 0.05, first_order),
  deg_PleCkin01p(PleCkin01p, DivKp, 0.05, first_order),
  deg_PleCpt4(PleCpt4,  DivKp,  0.05, first_order),
  deg_PleCpt2(PleCpt2,  DivKp, DivKp, 0.05, decomposition),

  /*DivJ*/
  syn_DivJf(mRNA_DivJ, DivJ_f, 0.0125*scale, syn), //New rate, only half of the original
  deg_DivJf(DivJ_f, none, 0.05, deg),
  lftdiff_DivJ(DivJ_f, none, 10.0, diffusion_lft, -1),
  rgtdiff_DivJ(DivJ_f, none, 10.0, diffusion_rgt, -1),

  bnd_DivJ(DivJ_f, DivJ, DivJ_st, 1.0, sticky),
  ubd_DivJ(DivJ, DivJ_f, 0.0, first_order),
  deg_DivJ(DivJ, none, 0.05, deg),

  bndK_DivJ(DivJ, DivK, DivJK, none, 5.0/scale, combination),
  ubdK_DivJK(DivJK, DivJ, DivK, 0.0016, decomposition),
  phos_DivJK(DivJK,  DivJKp, 5.0, first_order), //New rate, half of original
  deph_DivJKp(DivJKp,DivJK,  0.16, first_order),
  ubdKp_DivJKp(DivJKp, DivJ, DivKp, 1.0, decomposition),
  bndKp_DivJ(DivJ, DivKp, DivJKp, none, 5.0/scale, combination),

  deg_DivJK(DivJK, none, 0.05, deg),
  deg_DivJKp(DivJKp, none, 0.05, deg),

  phos_DivK(DivK,  DivKp, 0.0, first_order),//New reaction, base phosphorylation
  deph_DivKp(DivKp,DivK,  0.0, first_order),//new reaction
  /*CckA*/
  syn_CckA(mRNA_CckA, CckA, 0.0125*scale, syn),
  deg_CckA(CckA, none, 0.05, deg),
  lftdiff_CckA(CckA, none, 10.0, diffusion_lft, -1),
  rgtdiff_CckA(CckA, none, 10.0, diffusion_rgt, -1),

  bnd_CckA(CckA, CckA_phos, CckA_st, 1.0, sticky),
  ubd_CckA(CckA_phos, CckA, 0.1, first_order),
  ph2k_CckA(CckA_phos, CckA_kin, DivL, 10.0, hill, 0.5*scale),
  k2ph_CckA(CckA_kin,  CckA_phos, 1.0, first_order),
  deg_CckAkin(CckA_kin,   none, 0.05, deg),
  deg_CckAphos(CckA_phos, none, 0.05, deg),
  /*CtrA*/
  syn_CtrA(mRNA_CtrA, CtrA, 0.025*scale, syn),
  deg_CtrA(CtrA,   none, 0.05, deg),
  deg_CtrAp(CtrAp, none, 0.05, deg),
  lftdiff_CtrA(CtrA,  none, 10.0, diffusion_lft, -1),
  rgtdiff_CtrA(CtrA,  none, 10.0, diffusion_rgt, -1),
  lftdiff_CtrAp(CtrAp,none, 10.0, diffusion_lft, -1),
  rgtdiff_CtrAp(CtrAp,none, 10.0, diffusion_rgt, -1),

  phos_CtrA(CtrA,  CtrAp, CckA_kin, 600.0/scale, catalytic),
  deph_CtrAp(CtrAp, CtrA, CckA_phos, 600.0/scale, catalytic),
  /*DivL */
  syn_DivLf(mRNA_DivL, DivL_f, 0.0125*scale, syn),
  deg_DivLf(DivL_f, none, 0.05, deg),
  lftdiff_DivL(DivL_f, none, 10.0, diffusion_lft, -1),
  rgtdiff_DivL(DivL_f, none, 10.0, diffusion_rgt, -1),
  bnd_DivL(DivL_f, DivL, DivL_st, 1.0, sticky),
  ubd_DivL(DivL, DivL_f, 0.1, first_order),
  deg_DivL(DivL, none, 0.05, deg),
  
  bndKp_DivL(DivL, DivKp, DivLKp, none, 2.5/scale, combination),
  ubdKp_DivLKp(DivLKp, DivL, DivKp, 0.5, decomposition),
  deg_DivLKp(DivLKp, none, 0.05, deg){

  setLength(1.30);

  //*********************************************************
  //45 Species total
  //*********************************************************
  Species_List.push_back(& DivJ);
  Species_List.push_back(& DivJ_f);
  Species_List.push_back(& DivJ_st);
  Species_List.push_back(& DivJK);
  Species_List.push_back(& DivJKp);

  Species_List.push_back(& DivL);
  Species_List.push_back(& DivL_f);
  Species_List.push_back(& DivL_st);
  Species_List.push_back(& DivLKp);

  Species_List.push_back(& CckA);
  Species_List.push_back(& CckA_st);
  Species_List.push_back(& CckA_phos);
  Species_List.push_back(& CckA_kin);
  Species_List.push_back(& CtrA);
  Species_List.push_back(& CtrAp);

  Species_List.push_back(& DivK);
  Species_List.push_back(& DivKp);

  Species_List.push_back(& PleC_f);
  Species_List.push_back(& PleC_st);
  Species_List.push_back(& PleCphos);
  Species_List.push_back(& PleCph1);
  Species_List.push_back(& PleCph2);
  Species_List.push_back(& PleCkin11);
  Species_List.push_back(& PleCkin12);
  Species_List.push_back(& PleCkin22);
  Species_List.push_back(& PleCkin0);
  Species_List.push_back(& PleCkin2);
  Species_List.push_back(& PleCkin4);
  Species_List.push_back(& PleCkin1);
  Species_List.push_back(& PleCkin3);
  Species_List.push_back(& PleC2p);
  Species_List.push_back(& PleC1p);
  Species_List.push_back(& PleCkin02p);
  Species_List.push_back(& PleCkin10p);
  Species_List.push_back(& PleCkin01p);
  Species_List.push_back(& PleCpt4);
  Species_List.push_back(& PleCpt2);

  Species_List.push_back(& mRNA_DivJ); 
  Species_List.push_back(& mRNA_DivK); 
  Species_List.push_back(& mRNA_PleC); 
  Species_List.push_back(& mRNA_DivL); 
  Species_List.push_back(& mRNA_CckA); 
  Species_List.push_back(& mRNA_CtrA); 

  Species_List.push_back(& gene_ori); 
  Species_List.push_back(& gene_dup); 

  //***********************************************************
  //125+24+6+8 reactions
  //***********************************************************
  Reaction_List.push_back(& lftdiff_DivK);
  Reaction_List.push_back(& lftdiff_DivKp);
  Reaction_List.push_back(& lftdiff_PleC);
  Reaction_List.push_back(& lftdiff_DivJ);
  Reaction_List.push_back(& lftdiff_DivL);
  Reaction_List.push_back(& lftdiff_CckA);
  Reaction_List.push_back(& lftdiff_CtrA);
  Reaction_List.push_back(& lftdiff_CtrAp);

  Reaction_List.push_back(& rgtdiff_DivK);
  Reaction_List.push_back(& rgtdiff_DivKp);
  Reaction_List.push_back(& rgtdiff_PleC);
  Reaction_List.push_back(& rgtdiff_DivJ);
  Reaction_List.push_back(& rgtdiff_DivL);
  Reaction_List.push_back(& rgtdiff_CckA);
  Reaction_List.push_back(& rgtdiff_CtrA);
  Reaction_List.push_back(& rgtdiff_CtrAp);
  /*syn-dig-diff DivJ mRNA*/
  Reaction_List.push_back(& syn_OrimrDivJ);
  Reaction_List.push_back(& syn_DupmrDivJ);
  Reaction_List.push_back(& deg_mrDivJ);
  Reaction_List.push_back(& lftdiff_mrDivJ);  
  Reaction_List.push_back(& rgtdiff_mrDivJ);  
  /*syn-dig-diff DivK mRNA*/
  Reaction_List.push_back(& syn_OrimrDivK);
  Reaction_List.push_back(& syn_DupmrDivK);
  Reaction_List.push_back(& deg_mrDivK);
  Reaction_List.push_back(& lftdiff_mrDivK);  
  Reaction_List.push_back(& rgtdiff_mrDivK);  
  /*syn-dig-diff PleC mRNA*/
  Reaction_List.push_back(& syn_OrimrPleC);
  Reaction_List.push_back(& syn_DupmrPleC);
  Reaction_List.push_back(& deg_mrPleC);
  Reaction_List.push_back(& lftdiff_mrPleC);  
  Reaction_List.push_back(& rgtdiff_mrPleC);  
  /*syn-dig-diff DivL mRNA*/
  Reaction_List.push_back(& syn_OrimrDivL);
  Reaction_List.push_back(& syn_DupmrDivL);
  Reaction_List.push_back(& deg_mrDivL);
  Reaction_List.push_back(& lftdiff_mrDivL);  
  Reaction_List.push_back(& rgtdiff_mrDivL);  
  /*syn-dig-diff CckA mRNA*/
  Reaction_List.push_back(& syn_OrimrCckA);
  Reaction_List.push_back(& syn_DupmrCckA);
  Reaction_List.push_back(& deg_mrCckA);
  Reaction_List.push_back(& lftdiff_mrCckA);  
  Reaction_List.push_back(& rgtdiff_mrCckA);  
  /*syn-dig-diff CtrA mRNA*/
  Reaction_List.push_back(& syn_OrimrCtrA);
  Reaction_List.push_back(& syn_DupmrCtrA);
  Reaction_List.push_back(& deg_mrCtrA);
  Reaction_List.push_back(& lftdiff_mrCtrA);  
  Reaction_List.push_back(& rgtdiff_mrCtrA);  
  /*DivJ*/
  Reaction_List.push_back(& syn_DivJf);
  Reaction_List.push_back(& deg_DivJf);
  Reaction_List.push_back(& bnd_DivJ);
  Reaction_List.push_back(& ubd_DivJ);
  Reaction_List.push_back(& deg_DivJ);
  /*DivJ + DivK <-> DivJK <-> DivJKp <-> DivJ + DivKp*/
  Reaction_List.push_back(& bndK_DivJ);
  Reaction_List.push_back(& ubdK_DivJK);
  Reaction_List.push_back(& phos_DivJK);
  Reaction_List.push_back(& deph_DivJKp);
  Reaction_List.push_back(& ubdKp_DivJKp);
  Reaction_List.push_back(& bndKp_DivJ);
  
  Reaction_List.push_back(& deg_DivJK);
  Reaction_List.push_back(& deg_DivJKp);

  Reaction_List.push_back(& phos_DivK);
  Reaction_List.push_back(& deph_DivKp);
  /*CckA*/
  Reaction_List.push_back(& syn_CckA);
  Reaction_List.push_back(& deg_CckA);
  Reaction_List.push_back(& bnd_CckA);
  Reaction_List.push_back(& ubd_CckA);
  Reaction_List.push_back(& ph2k_CckA);
  Reaction_List.push_back(& k2ph_CckA);
  Reaction_List.push_back(& deg_CckAkin);
  Reaction_List.push_back(& deg_CckAphos);
  /*CtrA*/
  Reaction_List.push_back(& syn_CtrA);
  Reaction_List.push_back(& deg_CtrA);
  Reaction_List.push_back(& deg_CtrAp);
  Reaction_List.push_back(& phos_CtrA);
  Reaction_List.push_back(& deph_CtrAp);
  /*DivL */
  Reaction_List.push_back(& syn_DivLf);
  Reaction_List.push_back(& deg_DivLf);
  Reaction_List.push_back(& bnd_DivL);
  Reaction_List.push_back(& ubd_DivL);
  Reaction_List.push_back(& deg_DivL);
  
  Reaction_List.push_back(& bndKp_DivL);
  Reaction_List.push_back(& ubdKp_DivLKp);
  Reaction_List.push_back(& deg_DivLKp);
  /*DivK, DivKp*/
  Reaction_List.push_back(& syn_DivK);
  Reaction_List.push_back(& deg_DivK);
  Reaction_List.push_back(& deg_DivKp);
  /*PleC*/
  Reaction_List.push_back(& syn_PleCf);
  Reaction_List.push_back(& deg_PleCf);
  Reaction_List.push_back(& bnd_PleC);
  Reaction_List.push_back(& ubd_PleC);
  /*PleCph1, PleCph2*/
  Reaction_List.push_back(& bndKp_PleCphos);
  Reaction_List.push_back(& bndK_PleCphos);
  Reaction_List.push_back(& ubdKp_PleCph1);
  Reaction_List.push_back(& ubdK_PleCph2);

  Reaction_List.push_back(& phos_PleCph2);
  Reaction_List.push_back(& deph_PleCph1);
  Reaction_List.push_back(& disp_PleCph1);
  Reaction_List.push_back(& disp_PleCph2);
  /*PleCkin11, PleCkin12, PleCkin22*/
  Reaction_List.push_back(& bndKp_PleCph1);
  Reaction_List.push_back(& bndK_PleCph1);
  Reaction_List.push_back(& bndKp_PleCph2);
  Reaction_List.push_back(& bndK_PleCph2);
  Reaction_List.push_back(& ubdKp_PleCkin11);
  Reaction_List.push_back(& ubdKp_PleCkin12);
  Reaction_List.push_back(& ubdK_PleCkin12);
  Reaction_List.push_back(& ubdK_PleCkin22);
  /*PleCkin11p, PleC12p, PleC21p, PleC22p*/
  Reaction_List.push_back(& phos_PleCkin11);
  Reaction_List.push_back(& phos_PleCkin12);
  Reaction_List.push_back(& phos_PleCkin22);

  Reaction_List.push_back(& deph_PleCkin0);
  Reaction_List.push_back(& deph_PleCkin2);
  Reaction_List.push_back(& deph_PleCkin4);
  /*PleCkin1, PleCkin3*/
  Reaction_List.push_back(& ubdKp_PleCkin0);
  Reaction_List.push_back(& ubdKp_PleCkin2);
  Reaction_List.push_back(& ubdK_PleCkin2);
  Reaction_List.push_back(& ubdK_PleCkin4);

  Reaction_List.push_back(& bndKp_PleCkin1);
  Reaction_List.push_back(& bndK_PleCkin1);
  Reaction_List.push_back(& bndKp_PleCkin3);
  Reaction_List.push_back(& bndK_PleCkin3);
  /*PleC2p*/
  Reaction_List.push_back(& ubdKp_PleCkin1);
  Reaction_List.push_back(& ubdK_PleCkin3);

  Reaction_List.push_back(& bndKp_PleC2p);
  Reaction_List.push_back(& bndK_PleC2p);

  Reaction_List.push_back(& deph_PleC2p);
  /*PleC1p*/
  Reaction_List.push_back(& ubdKp_PleCkin10p);
  Reaction_List.push_back(& ubdKp_PleCkin01p);
  Reaction_List.push_back(& ubdK_PleCkin02p);

  Reaction_List.push_back(& bndKp_PleC1pt);
  Reaction_List.push_back(& bndKp_PleC1ps);
  Reaction_List.push_back(& bndK_PleC1ps);
 
  Reaction_List.push_back(& deph_PleC1p);
  /*auto-phosphation*/
  Reaction_List.push_back(& autoph_PleCkin3);
  Reaction_List.push_back(& autoph_PleCkin2);
  Reaction_List.push_back(& autoph_PleCkin4);
  Reaction_List.push_back(& autoph_PleCpt4);

  Reaction_List.push_back(& deauto_PleCkin10p);
  Reaction_List.push_back(& deauto_PleCpt2);
  Reaction_List.push_back(& deauto_PleCpt4);
  Reaction_List.push_back(& deauto_PleCkin11);
  /*PleCKin10p, PleCKin01p, PleCkin20p, PleCkin02p*/
  Reaction_List.push_back(& deph_PleCkin1t);
  Reaction_List.push_back(& deph_PleCkin3t);

  Reaction_List.push_back(& phos_PleCkin01p);
  Reaction_List.push_back(& phos_PleCkin02p);
  /*PleCkin11p, PleCkin12p, PleCkin21p, PleCkin22p*/
  Reaction_List.push_back(& bndKp_PleCkin01p);
  Reaction_List.push_back(& bndKp_PleCkin02p);

  Reaction_List.push_back(& ubdKp_PleCpt2s);
  Reaction_List.push_back(& ubdKp_PleCpt4);
  /*Degradation*/
  Reaction_List.push_back(& deg_PleCphos);
  Reaction_List.push_back(& deg_PleCph1);
  Reaction_List.push_back(& deg_PleCph2); 
  Reaction_List.push_back(& deg_PleCkin11); 
  Reaction_List.push_back(& deg_PleCkin12);
  Reaction_List.push_back(& deg_PleCkin22);
  Reaction_List.push_back(& deg_PleCkin0);
  Reaction_List.push_back(& deg_PleCkin2);
  Reaction_List.push_back(& deg_PleCkin4);
  Reaction_List.push_back(& deg_PleCkin1);
  Reaction_List.push_back(& deg_PleCkin3);
  Reaction_List.push_back(& deg_PleC2p);
  Reaction_List.push_back(& deg_PleC1p);
  Reaction_List.push_back(& deg_PleCkin02p);
  Reaction_List.push_back(& deg_PleCkin10p);
  Reaction_List.push_back(& deg_PleCkin01p);
  Reaction_List.push_back(& deg_PleCpt4);
  Reaction_List.push_back(& deg_PleCpt2);

}


//****************************************
//Cell size growth
//***************************************
void PleC::binLength_update(double tau){ len=len+len*mu*tau; }
void PleC::setLength(double leng){ len=leng/BINNUM; }
double PleC::getLength(){ return len;}


void PleC::DNAorigination(){
  gene_ori.setInit(); gene_dup.setInit(); 
  gene_ori.population[BINNUM*4/5] = 1; gene_ori.total_population = 1; 
}

void PleC::DNAreplication(){
  gene_dup.population[BINNUM/5-1] = 1; gene_dup.total_population = 1; 
}


void PleC::InitStage(){
  setLength(1.30);
  for(int i=0; i<Species_List.size(); i++) Species_List[i]->setInit();
}

//////////////////////////////////////////////////////////////
////void PleC::setSwarmer(char fname[]){
////  ifstream datafile(fname, ios::in);
////
////  for(int i=0; i<Species_List.size(); ++i)
////    for(int k = 0; k<Species_List[i]->size(); ++k)
////      datafile>>Species_List[i]->population[k];
////
////  for(int i=0; i<Species_List.size(); ++i)
////      datafile>>Species_List[i]->total_population;
////
////  datafile.close();
////
////  setLength(1.30);
////}
////////////////////////////////////////////////////////////////

void PleC::initSpatial(){
  PleC_st.setInit();
  for(int i=BINNUM*9/10; i<BINNUM; i++)PleC_st.population[i]=1;
  PleC_st.total_population = BINNUM/10;

  DivJ_st.setInit();
  DivL_st.uniFill();
  CckA_st.uniFill();
}

void PleC::introDivJ(){
  DivJ_st.setInit();
  for(int i=BINNUM*9/10; i<BINNUM; i++)DivJ_st.population[i]=1;
  DivJ_st.total_population = BINNUM/10;

  phos_DivK.setReactionRate(0.05);//New reaction, base phosphorylation
  deph_DivKp.setReactionRate(0.01);//new reaction
}

void PleC::clearPleC(){
  PleC_st.setInit();

  CckA_st.setInit();
  for(int i=BINNUM*9/10; i<BINNUM; i++)CckA_st.population[i]=1;
  CckA_st.total_population = BINNUM/10;

}

void PleC::checkPoint3(){
  PleC_st.setInit();
  for(int i=0; i<BINNUM/10; i++)PleC_st.population[i] = 1;
  PleC_st.total_population = BINNUM/10;

  DivL_st.setInit();
  for(int i=0; i<BINNUM/10; i++)DivL_st.population[i] = 1;
  DivL_st.total_population = BINNUM/10;

  CckA_st.setInit();
  for(int i=0; i<BINNUM/10; i++)CckA_st.population[i] = 1;
  for(int i=BINNUM*9/10; i<BINNUM; i++)CckA_st.population[i]=1;
  CckA_st.total_population = 2*BINNUM/10;  
}



void PleC::DivJtranslation(double foldup){

  double rt_fu = syn_DivJf.getReactionRate()*foldup;
  syn_DivJf.setReactionRate(rt_fu);
}

void PleC::DivKOverExpression(int foldup){

  double rt_fu = syn_DivK.getReactionRate()*foldup;
  syn_DivK.setReactionRate(rt_fu);
}


void PleC::compartmentization(){
 //*********************************************
 //Barriers for mRNA
 //*********************************************
 lftdiff_mrDivJ.setBarrier(BINNUM/2);
 rgtdiff_mrDivJ.setBarrier(BINNUM/2);
 
 lftdiff_mrDivK.setBarrier(BINNUM/2);
 rgtdiff_mrDivK.setBarrier(BINNUM/2);

 lftdiff_mrPleC.setBarrier(BINNUM/2);
 rgtdiff_mrPleC.setBarrier(BINNUM/2);
 
 lftdiff_mrDivL.setBarrier(BINNUM/2);
 rgtdiff_mrDivL.setBarrier(BINNUM/2);
 
 lftdiff_mrCckA.setBarrier(BINNUM/2);
 rgtdiff_mrCckA.setBarrier(BINNUM/2);
 
 lftdiff_mrCtrA.setBarrier(BINNUM/2);
 rgtdiff_mrCtrA.setBarrier(BINNUM/2);
 
 //***********************************************
 //Barrier for proteins
 //***********************************************
 lftdiff_DivK.setBarrier(BINNUM/2);
 rgtdiff_DivK.setBarrier(BINNUM/2);
 lftdiff_DivKp.setBarrier(BINNUM/2);
 rgtdiff_DivKp.setBarrier(BINNUM/2);

 lftdiff_PleC.setBarrier(BINNUM/2);
 rgtdiff_PleC.setBarrier(BINNUM/2);
 
 lftdiff_DivJ.setBarrier(BINNUM/2);
 rgtdiff_DivJ.setBarrier(BINNUM/2);
 
 lftdiff_CtrA.setBarrier(BINNUM/2);
 rgtdiff_CtrA.setBarrier(BINNUM/2);
 lftdiff_CtrAp.setBarrier(BINNUM/2);
 rgtdiff_CtrAp.setBarrier(BINNUM/2);

 lftdiff_CckA.setBarrier(BINNUM/2);
 rgtdiff_CckA.setBarrier(BINNUM/2);
 
 lftdiff_DivL.setBarrier(BINNUM/2);
 rgtdiff_DivL.setBarrier(BINNUM/2);

}

void PleC::cellFixed(){  
  mu=0.0; 
  syn_OrimrDivJ.setReactionRate(0.0);
  syn_DupmrDivJ.setReactionRate(0.0);
}

void PleC::cellGrowth(){ 
  mu=0.005577; 
  syn_OrimrDivJ.setReactionRate(0.625);
  syn_DupmrDivJ.setReactionRate(0.625);
}

////
////void PleC::InitFirings(){
////  firings.clear();
////  firings.resize(Reaction_List.size(), 0.0);
////}
////
void PleC::printFirings(string filename){

  ofstream ofile(filename.c_str(), ios::out|ios::app);

  for(int i=0; i<12; i++) ofile<<firings[i]<<"\t";
  ofile<<endl;

  ofile.close();
}


void PleC::printPopulation(){

  for(int i=0; i<Species_List.size(); ++i)
    for(int k = 0; k<Species_List[i]->size(); ++k)
      cout<<Species_List[i]->population[k]<<'\t';

  for(int i=0; i<Species_List.size(); ++i)
      cout<<Species_List[i]->total_population<<'\t';

//  for(int i=0; i<firings.size(); ++i)
//      cout<<firings[i]<<"      ";
  cout<<len<<endl;
}


void PleC::printPopulation(string filename){

  ofstream ofile(filename.c_str(), ios::out|ios::app);

  for(int i=0; i<Species_List.size(); ++i)
    for(int k = 0; k<Species_List[i]->size(); ++k)
      ofile<<Species_List[i]->population[k]<<'\t';

  for(int i=0; i<Species_List.size(); ++i)
      ofile<<Species_List[i]->total_population<<'\t';

  ofile<<len<<endl;

  ofile.close();

}


void PleC::printPropensity(string filename){

  ofstream ofile(filename.c_str(), ios::out|ios::app);

  for(int i=0; i<Reaction_List.size(); ++i)
    ofile<<Reaction_List[i]->get_propensity()<<endl;
  ofile<<a0<<endl;  

  ofile.close();

}

void PleC::printPleC(){
  cout<<"Model PleC has "<<Species_List.size()<<" species and "
      <<Reaction_List.size()<<" reactions"<<endl;

}

double PleC::uniRandNum(){

  double r1; 
  do{r1=1.0*rand()/RAND_MAX;}while(r1<=0 || r1>=1);

  return r1;
}
/////////////////////////////////////////////////////////////////
//void popZ::population_decrease(){
//  for(int i=0; i<Reaction_List.size(); ++i)Reaction_List[i]->population_update(length);
//   
//}
///////////////////////////////////////////////////////////////
void PleC::propensity_update(int r_ita, int b_ita){

  switch(Reaction_List[r_ita]->getReactionType()){
   case diffusion_lft:
    for(int i=0; i<Reaction_List.size(); ++i)Reaction_List[i]->propensity_update(b_ita, b_ita-1, len);
    break;

   case diffusion_rgt:
    for(int i=0; i<Reaction_List.size(); ++i)Reaction_List[i]->propensity_update(b_ita, b_ita+1, len);
    break;

   default: 
    for(int i=0; i<Reaction_List.size(); ++i)Reaction_List[i]->propensity_update(b_ita, len);
   break; 
  }

  a0=0.0;
  for(int i=0; i<Reaction_List.size(); ++i)a0+=Reaction_List[i]->get_propensity();
}


void PleC::propensity_calculation(){

  for(int i=0; i<Reaction_List.size(); ++i)Reaction_List[i]->cal_propensity(len);

  a0=0.0;
  for(int i=0; i<Reaction_List.size(); ++i)a0+=Reaction_List[i]->get_propensity();
}


void PleC::SSA(double t,string fname, double timeStep){

  double timeTracker=0.0;
  double Timetag=timeStep;

  double tau;
  int rule_ita, bin_ita;
  double ra0;

  //**********************************************
  //Get ready to the SSA 
  //**********************************************
  propensity_calculation();
  
  while(timeTracker < t){
    tau=-1.0/a0*log(uniRandNum());
    timeTracker+=tau;
    //***************************************
    //Keep the population at every min
    //***************************************
    while(timeTracker > Timetag){
      printPopulation(fname);   Timetag += timeStep;
    }

    ra0 = a0 * uniRandNum();

    rule_ita=0;
    while(ra0 > Reaction_List[rule_ita]->get_propensity()){
      ra0 = ra0 - Reaction_List[rule_ita]->get_propensity();
      rule_ita++; 
    }
////    cout<<rule_ita<<'\t'<<ra0<<'\t'<<len<<'\t';
    
    if(rule_ita>=16) firings[Reaction_List[rule_ita]->getReactionType()]++;
    bin_ita = Reaction_List[rule_ita]->population_update(ra0, len);
////    cout<<bin_ita<<endl;
   
    if(mu != 0)binLength_update(tau);
    propensity_update(rule_ita, bin_ita);
 }

 cout<<timeTracker<<endl;
}


void PleC::clearFire(){
    for(int i=0;i<12;i++){
      firings[i] = 0;
    } 
}

