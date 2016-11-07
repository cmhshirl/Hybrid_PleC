#include "Reaction.h"
#include "LftDiffusion.h"


LftDiffusion :: LftDiffusion(Species& r1, Species& p1, double r, Reaction_type tp, int bar):Reaction(r1, p1, r, tp){  barrier = bar; }


void LftDiffusion::setBarrier(int bar){ barrier = bar; }


void LftDiffusion::cal_propensity(double h){

  dif_population = reactant.total_population-reactant.population[0];
  if(barrier > 0 && barrier < reactant.size()) dif_population -= reactant.population[barrier];

  propensity = rate/(h*h) * dif_population;

}

void LftDiffusion::propensity_update(int bin_num, double h){
  cal_propensity(h);
}

void LftDiffusion::propensity_update(int bin_num1, int bin_num2, double h){
  cal_propensity(h);
}


int LftDiffusion::population_update(double a_res, double h){

  double bin_ra0 = a_res*h*h/rate;

  int b_ita=1;
  while(bin_ra0 > reactant.population[b_ita]){
    if(b_ita != barrier)bin_ra0 -= reactant.population[b_ita]; 
    b_ita++;
    if(b_ita == barrier) b_ita++;
  }

  reactant.population[b_ita] --;
  reactant.population[b_ita-1] ++;  

  return b_ita;  
}

