#include <cstdlib>
using std::rand;


#include "Reaction.h"
#include "HillReaction.h"

HillReaction :: HillReaction(Species& r1, Species& p1, Species& au, double rt, Reaction_type tp, double km)
  :Reaction(r1, p1, rt, tp), aux(au){
  
  Kmdl = km*km*km*km;

  bin_prop = new double[reactant.size()];
  total_bin_prop = 0.0;

  epsilon = 0.01;
}


double HillReaction :: KmdlInPopulation(double h){

  return Kmdl*h*h*h*h;
}

//********************************************************
//Smooth DivL for the 5 neighboring bins
//
//********************************************************
void HillReaction :: cal_propensity(double h){
  double Kmdl_ppl = Kmdl*h*h*h*h;

  double meanPPL = 0.0;
  
  for(int i=0; i<10; i++){
    meanPPL = 0.0;
    for(int j=0; j<aux.size()/10; j++) meanPPL += aux.population[i*aux.size()/10+j];
    meanPPL = meanPPL/(aux.size()/10);
    meanPPL = meanPPL*meanPPL*meanPPL*meanPPL;
    for(int j=0; j<aux.size()/10; j++) bin_prop[i*aux.size()/10+j] = meanPPL;
  }


  for(int i=0; i<reactant.size(); i++){ 
    bin_prop[i] = reactant.population[i]*bin_prop[i]/(Kmdl_ppl+bin_prop[i]);
  }

  total_bin_prop = 0.0;
  for(int i=0; i<reactant.size(); i++) total_bin_prop += bin_prop[i];

  propensity = rate*total_bin_prop;

}
//*********************************************************
//*********************************************************
//New Hill function
// gamma = (epsilon * Km^4 + DivL^4)/(Km^4+DivL^4)
//
//*********************************************************
//void HillReaction :: cal_propensity(double h){
//  double Kmdl_ppl = Kmdl*h*h*h*h;
//
//  double aux_ppl = 0;
//  
//  for(int i=0; i<reactant.size(); i++){ 
//    aux_ppl = aux.population[i]*aux.population[i]*aux.population[i]*aux.population[i];
//    bin_prop[i] = reactant.population[i]*(Kmdl_ppl * epsilon + aux_ppl)/(Kmdl_ppl + aux_ppl);
//  }
//
//  total_bin_prop = 0.0;
//  for(int i=0; i<reactant.size(); i++) total_bin_prop += bin_prop[i];
//
//  propensity = rate*total_bin_prop;
//
//}
//*************************************************************

void HillReaction::propensity_update(int bin_num, double h){
  cal_propensity(h);
}

void HillReaction::propensity_update(int bin_num1, int bin_num2, double h){
  cal_propensity(h);
}


int HillReaction :: population_update(double a_res, double h){
  
  double resid = a_res/rate;
  int bin_ita=0;
  while(resid > bin_prop[bin_ita]){
    resid -= bin_prop[bin_ita]; bin_ita++;
  }

  reactant.population[bin_ita] --;  reactant.total_population --;
  product.population[bin_ita] ++;   product.total_population ++;
    
  return bin_ita;
}


