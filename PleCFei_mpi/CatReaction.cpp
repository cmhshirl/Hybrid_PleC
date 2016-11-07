#include <cstdlib>

#include <iostream>
using std::cout;
using std::endl;

#include "Reaction.h"
#include "CatReaction.h"

CatReaction :: CatReaction(Species& r1, Species& p1, Species& au,double rt, Reaction_type tp) : Reaction(r1, p1, rt, tp), aux(au){

  bin_prop = new int[reactant.size()];
  total_bin_prop = 0;
}


void CatReaction :: cal_propensity(double h){
  for(int i=0; i<reactant.size(); i++)
    bin_prop[i] = reactant.population[i] * aux.population[i];

  total_bin_prop=0;
  for(int i=0; i<reactant.size(); i++) total_bin_prop += bin_prop[i];

  switch(type){
   case catalytic:
    propensity = total_bin_prop*rate/h;
    break;

   case sticky:
    propensity = total_bin_prop*rate;
    break;

   default:  
    cout<<"Reaction_type : "<< type <<" not found "<<endl;
    exit(1);
  }

}

void CatReaction::propensity_update(int bin_num, double h){
  total_bin_prop -= bin_prop[bin_num];
  bin_prop[bin_num] = reactant.population[bin_num]*aux.population[bin_num];
  total_bin_prop += bin_prop[bin_num];

  switch(type){
   case catalytic:
    propensity = total_bin_prop*rate/h;
    break;

   case sticky:
    propensity = total_bin_prop*rate;
    break;

   default:  
    cout<<"Reaction_type : "<< type <<" not found "<<endl;
    exit(1);
  }

}


void CatReaction::propensity_update(int bin_num1, int bin_num2, double h){

  total_bin_prop -= (bin_prop[bin_num1]+bin_prop[bin_num2]);

  bin_prop[bin_num1] = reactant.population[bin_num1]*aux.population[bin_num1];
  bin_prop[bin_num2] = reactant.population[bin_num2]*aux.population[bin_num2];

  total_bin_prop += (bin_prop[bin_num1]+bin_prop[bin_num2]);

  switch(type){
   case catalytic:
    propensity = total_bin_prop*rate/h;
    break;

   case sticky:
    propensity = total_bin_prop*rate;
    break;

   default:  
    cout<<"Reaction_type : "<< type <<" not found "<<endl;
    exit(1);
  }
}

int CatReaction::population_update(double a_res, double h){
  double bin_ra0;

  if(type==catalytic) bin_ra0 = a_res*h/rate;
  else bin_ra0 = a_res/rate;

  int b_ita=0;
  while(bin_ra0>bin_prop[b_ita]){
    bin_ra0 -= bin_prop[b_ita];  b_ita++;
  }

  reactant.population[b_ita] --; reactant.total_population --;
  product.population[b_ita] ++;  product.total_population ++;

  return b_ita;
}

