#include "Reaction.h"
#include "RgtDiffusion.h"


RgtDiffusion :: RgtDiffusion(Species& r1, Species& p1, double r, Reaction_type tp, int bar):Reaction(r1, p1, r, tp){barrier = bar;}

void RgtDiffusion::setBarrier(int bar){barrier = bar; }

void RgtDiffusion::cal_propensity(double h){

  dif_population = reactant.total_population-reactant.population[reactant.size()-1];

  if(barrier > 0 && barrier < reactant.size()) dif_population -= reactant.population[barrier-1];

  propensity = rate/(h*h) * dif_population;

}

void RgtDiffusion::propensity_update(int bin_num, double h){
  cal_propensity(h);
}


void RgtDiffusion::propensity_update(int bin_num1, int bin_num2, double h){
  cal_propensity(h);
}


int RgtDiffusion::population_update(double a_res, double h){

  double bin_ra0 = a_res*h*h/rate;

  int b_ita=0;
  while(bin_ra0 > reactant.population[b_ita]){
    bin_ra0 -= reactant.population[b_ita]; 
    b_ita++;
    if(b_ita == barrier-1) b_ita++;
  }

  reactant.population[b_ita] --;  reactant.population[b_ita+1] ++;  
  return b_ita;    
}

