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


#ifndef POPZ_H
#define POPZ_H

#include <vector>
using std::vector;

#include <string>
using std::string; 

#include "Reaction.h"
#include "CatReaction.h"
#include "DecmpstReaction.h"
#include "CmbReaction.h"
#include "HillReaction.h"
#include "LftDiffusion.h"
#include "RgtDiffusion.h"

class PleC{
 public:
  
  PleC();

  void DivJtranslation(double);
  void DivKOverExpression(int);

  double getLength();
  void setLength(double);
  void binLength_update(double);
  void clearFire();

////  void InitFirings();
  void InitStage();
  void setSwarmer(char[]);

  void initSpatial();
  void introDivJ();
  void clearPleC();
  void checkPoint3();
  void compartmentization();

  void cellFixed();
  void cellGrowth(); 
  void DNAorigination();
  void DNAreplication();

  void propensity_calculation();
  void propensity_update(int, int);

  void SSA(double, string, double);

  void printPopulation();
  void printPopulation(string);
  void printPropensity(string);
  void printPleC();
  void printFirings(string);


 private:
  static const int scale = 1000;
  static const int BINNUM = 50;
  Species none;

  Species gene_ori, gene_dup;
  Species mRNA_DivJ, mRNA_DivK, mRNA_PleC;
  Species mRNA_DivL, mRNA_CckA, mRNA_CtrA;

  Species DivJ, DivJ_f, DivJ_st;
  Species DivJK, DivJKp;

  Species DivL, DivL_f, DivL_st;
  Species DivLKp;

  Species CckA, CckA_st;
  Species CckA_phos, CckA_kin;
  Species CtrA, CtrAp;

  Species DivK, DivKp;

  Species PleC_f, PleC_st;
  Species PleCphos, PleCph1, PleCph2;
  Species PleCkin11, PleCkin12, PleCkin22;
  Species PleCkin0, PleCkin2, PleCkin4;
  Species PleCkin1, PleCkin3;
  Species PleC2p, PleC1p;
  Species PleCkin02p, PleCkin10p, PleCkin01p;
  Species PleCpt4, PleCpt2;

  /*mRNA reactions*/
  Reaction syn_OrimrDivJ, syn_DupmrDivJ, deg_mrDivJ;
  Reaction syn_OrimrDivK, syn_DupmrDivK, deg_mrDivK;
  Reaction syn_OrimrPleC, syn_DupmrPleC, deg_mrPleC;
  Reaction syn_OrimrDivL, syn_DupmrDivL, deg_mrDivL;
  Reaction syn_OrimrCckA, syn_DupmrCckA, deg_mrCckA;
  Reaction syn_OrimrCtrA, syn_DupmrCtrA, deg_mrCtrA;
  LftDiffusion lftdiff_mrDivJ, lftdiff_mrDivK, lftdiff_mrPleC;
  RgtDiffusion rgtdiff_mrDivJ, rgtdiff_mrDivK, rgtdiff_mrPleC;
  LftDiffusion lftdiff_mrDivL, lftdiff_mrCckA, lftdiff_mrCtrA;
  RgtDiffusion rgtdiff_mrDivL, rgtdiff_mrCckA, rgtdiff_mrCtrA;

  /*Protein Diffusion*/
  LftDiffusion lftdiff_DivJ, lftdiff_DivK, lftdiff_DivKp, lftdiff_PleC;
  RgtDiffusion rgtdiff_DivJ, rgtdiff_DivK, rgtdiff_DivKp, rgtdiff_PleC;
  LftDiffusion lftdiff_DivL, lftdiff_CckA, lftdiff_CtrAp, lftdiff_CtrA;
  RgtDiffusion rgtdiff_DivL, rgtdiff_CckA, rgtdiff_CtrAp, rgtdiff_CtrA;

  /*DivK, DivKp*/
  Reaction syn_DivK, deg_DivK, deg_DivKp;
  /*PleC*/
  Reaction syn_PleCf, deg_PleCf, ubd_PleC;
  CatReaction bnd_PleC;
  /*degradation PleC cmplx*/
  DecmpstReaction deg_PleCkin11, deg_PleCkin0, deg_PleCpt2;
  Reaction deg_PleCphos, deg_PleCph1, deg_PleCph2;
  Reaction deg_PleCkin12, deg_PleCkin22;
  Reaction deg_PleCkin2, deg_PleCkin4, deg_PleCkin1, deg_PleCkin3;
  Reaction deg_PleC2p, deg_PleC1p;
  Reaction deg_PleCkin02p, deg_PleCkin10p, deg_PleCkin01p, deg_PleCpt4;

  /*PleCph1, PleCph2*/
  CmbReaction bndKp_PleCphos, bndK_PleCphos;
  DecmpstReaction ubdKp_PleCph1, ubdK_PleCph2;
  Reaction phos_PleCph2, deph_PleCph1;
  CmbReaction disp_PleCph1, disp_PleCph2;
  /*PleCkin11, PleCkin12, PleCkin22*/
  CmbReaction bndKp_PleCph1, bndK_PleCph1, bndKp_PleCph2, bndK_PleCph2;
  DecmpstReaction ubdKp_PleCkin11, ubdKp_PleCkin12, ubdK_PleCkin12, ubdK_PleCkin22;
  /*PleCkin0, PleC2, PleC4*/
  Reaction phos_PleCkin11, phos_PleCkin12, phos_PleCkin22;
  Reaction deph_PleCkin0, deph_PleCkin2, deph_PleCkin4;
  /*PleCkin1, PleCkin3*/
  DecmpstReaction ubdKp_PleCkin0, ubdKp_PleCkin2, ubdK_PleCkin2, ubdK_PleCkin4;
  CmbReaction bndKp_PleCkin1, bndK_PleCkin1, bndKp_PleCkin3, bndK_PleCkin3;
  /*PleC2p*/
  DecmpstReaction ubdKp_PleCkin1, ubdK_PleCkin3;
  CmbReaction bndKp_PleC2p, bndK_PleC2p;
  Reaction deph_PleC2p;
  /*PleC1p*/
  DecmpstReaction ubdKp_PleCkin10p, ubdKp_PleCkin01p, ubdK_PleCkin02p;
  CmbReaction bndKp_PleC1pt, bndKp_PleC1ps, bndK_PleC1ps;
 
  Reaction deph_PleC1p;
  /*auto-phosphation*/
  Reaction autoph_PleCkin3, autoph_PleCkin2, autoph_PleCkin4, autoph_PleCpt4;
  Reaction deauto_PleCkin10p, deauto_PleCpt2, deauto_PleCpt4, deauto_PleCkin11;
  /*PleCKin10p, PleCKin01p, PleCkin20p, PleCkin02p*/
  Reaction deph_PleCkin1t, deph_PleCkin3t;
  Reaction phos_PleCkin01p, phos_PleCkin02p;
  /*PleCkin11p, PleCkin12p, PleCkin21p, PleCkin22p*/
  CmbReaction bndKp_PleCkin01p, bndKp_PleCkin02p;
  DecmpstReaction ubdKp_PleCpt2s, ubdKp_PleCpt4;

  /*DivJ*/
  Reaction syn_DivJf, deg_DivJf;
  CatReaction bnd_DivJ; 
  Reaction ubd_DivJ, deg_DivJ;

  CmbReaction bndK_DivJ, bndKp_DivJ;
  DecmpstReaction ubdK_DivJK, ubdKp_DivJKp;
  Reaction phos_DivJK, deph_DivJKp;
  Reaction deg_DivJK, deg_DivJKp;

  Reaction phos_DivK, deph_DivKp;
  /*CckA*/
  Reaction syn_CckA, deg_CckA;
  CatReaction bnd_CckA;
  Reaction ubd_CckA;
  HillReaction ph2k_CckA;
  Reaction k2ph_CckA, deg_CckAkin, deg_CckAphos;
  /*CtrA*/
  Reaction syn_CtrA, deg_CtrA, deg_CtrAp;
  CatReaction phos_CtrA, deph_CtrAp;
  /*DivL */
  Reaction syn_DivLf, deg_DivLf;
  CatReaction bnd_DivL;
  Reaction ubd_DivL, deg_DivL;
  
  CmbReaction bndKp_DivL; DecmpstReaction ubdKp_DivLKp;
  Reaction deg_DivLKp;

  //******************************************************
  //Reaction and Species List
  //******************************************************
  vector<Species *> Species_List;
  vector<Reaction *> Reaction_List;

  //vector<double> firings;
  double a0;
  double len, mu;
  long long int firings[12];
  double uniRandNum();

};

#endif
