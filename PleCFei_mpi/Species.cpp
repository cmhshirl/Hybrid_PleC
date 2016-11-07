#include <new>
#include <cstdlib>
using std::rand;

#include "Species.h"

Species::Species():bin_num(1), total_population(0) 
  {population = new int[1]; population[0]=0;}


Species::Species(int bin): bin_num(bin){
    total_population=0;  population = new int[bin];
    for(int i=0; i<bin; i++) population[i] = 0;
}
 
void Species::setInit(){
    for(int i=0; i<bin_num; i++) population[i] = 0;
    total_population = 0; 
}

void Species::uniFill(){
    for(int i=0; i<bin_num; i++) population[i] = 1;
    total_population = bin_num; 
}

int Species::size(){return bin_num;}

int Species::randBinSelect(double bin_ra0){

  int b_ita=0;
  while(bin_ra0 > population[b_ita]){
    bin_ra0 -= population[b_ita];
    b_ita++; 
   }

  return b_ita;

}


