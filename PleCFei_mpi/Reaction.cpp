#include <iostream>
using std::cout;
using std::endl;

#include <cstdlib>
using std::rand;

#include "Reaction.h"

Reaction::Reaction(Species& s1, Species& p1, double r, Reaction_type tp): reactant(s1), product(p1), rate(r), type(tp){}
  

double Reaction::get_propensity(){
  return propensity;
}


void Reaction::setReactionRate(double rt){
  rate = rt;
}


double Reaction::getReactionRate(){
  return rate;
}


Reaction_type Reaction::getReactionType(){
  return type;
}


void Reaction::propensity_update(int bin_num, double h){
  cal_propensity(h);
}

void Reaction::propensity_update(int bin_num1, int bin_num2, double h){
  cal_propensity(h);
}

void Reaction::cal_propensity(double h){

  switch(type){

   case zero:
    propensity = rate * h * product.size();
    break;

   case syn:
   case deg:
   case first_order:
    propensity = rate * reactant.total_population;
    break;

   default:  
    cout<<"Reaction_type : "<< type <<" not found "<<endl;
    exit(1);
  }
}

int Reaction::population_update(double a_res, double h){

  int bin_ita;

  switch(type){

   case first_order: //A -> B
    bin_ita = reactant.randBinSelect(a_res/rate);
    reactant.population[bin_ita] --;  reactant.total_population --;
    product.population[bin_ita] ++;   product.total_population ++;
    break;
    
   case deg: //A -> null
    bin_ita = reactant.randBinSelect(a_res/rate);
    reactant.population[bin_ita] --;  reactant.total_population --;
    break;
    
   case syn: // null -> A; aux
    bin_ita = reactant.randBinSelect(a_res/rate);
    product.population[bin_ita] ++;  product.total_population ++;
    break;  
  
   case zero: // null -> A
    bin_ita = (int)(a_res/(rate*h));
    product.population[bin_ita] ++;  product.total_population ++;
    break;  
  
   default:
    cout<<"Reaction type : "<<type<<" not found"<<endl;
    exit(1);
    break;
  }

  return bin_ita;
}

