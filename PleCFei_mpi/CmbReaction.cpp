#include "Reaction.h"
#include "CmbReaction.h"

CmbReaction :: CmbReaction(Species& r1, Species& r2, Species& p1, Species& p2, double rt, Reaction_type tp) : Reaction(r1, p1, rt, tp), reactant2(r2), product2(p2){

  bin_prop = new int[reactant.size()];
  total_bin_prop = 0;

}


void CmbReaction :: cal_propensity(double h){
  for(int i=0; i<reactant.size(); i++)
    bin_prop[i] = reactant.population[i] * reactant2.population[i];

  total_bin_prop = 0;
  for(int i=0; i<reactant.size(); i++) total_bin_prop += bin_prop[i];

  propensity = total_bin_prop*rate/h;
}


void CmbReaction::propensity_update(int bin_num, double h){
  total_bin_prop -= bin_prop[bin_num];
  bin_prop[bin_num] = reactant.population[bin_num]*reactant2.population[bin_num];
  total_bin_prop += bin_prop[bin_num];

  propensity = total_bin_prop*rate/h;

}


void CmbReaction::propensity_update(int bin_num1, int bin_num2, double h){
  total_bin_prop -= (bin_prop[bin_num1]+bin_prop[bin_num2]);
  bin_prop[bin_num1] = reactant.population[bin_num1]*reactant2.population[bin_num1];
  bin_prop[bin_num2] = reactant.population[bin_num2]*reactant2.population[bin_num2];
  total_bin_prop += (bin_prop[bin_num1]+bin_prop[bin_num2]);

  propensity = total_bin_prop*rate/h;

}


int CmbReaction::population_update(double a_res, double h){
  double bin_ra0 = a_res*h/rate;

  int b_ita=0;
  while(bin_ra0>bin_prop[b_ita]){
    bin_ra0 -= bin_prop[b_ita];  b_ita++;
  }

  reactant.population[b_ita] --;  reactant.total_population --;
  reactant2.population[b_ita] --; reactant2.total_population --;
  product.population[b_ita] ++;  product.total_population ++;
  
  if(type == displacement){
    product2.population[b_ita] ++;  product2.total_population ++;
  }
  return b_ita;
}

