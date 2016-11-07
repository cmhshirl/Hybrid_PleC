#include "Reaction.h"
#include "DecmpstReaction.h"


DecmpstReaction :: DecmpstReaction(Species& r1, Species& p1, Species& p2, double r, Reaction_type tp):Reaction(r1, p1, r, tp), product2(p2){}


void DecmpstReaction::cal_propensity(double h){
  propensity = rate * reactant.total_population;
}

void DecmpstReaction::propensity_update(int bin_num, double h){
  propensity = rate * reactant.total_population;
}

void DecmpstReaction::propensity_update(int bin_num1, int bin_num2, double h){
  propensity = rate * reactant.total_population;
}

int DecmpstReaction::population_update(double a_res, double h){

  int bin_ita = reactant.randBinSelect(a_res/rate);

  reactant.population[bin_ita] --;  reactant.total_population --;
  product.population[bin_ita] ++;   product.total_population ++;
  product2.population[bin_ita] ++;  product2.total_population ++;
    
  return bin_ita;
}

