#ifndef CMBREACTION_H
#define CMBREACTION_H

#include "Species.h"
#include "Reaction.h"

class CmbReaction : public Reaction{
 public:
  CmbReaction(Species&, Species&, Species&, Species&, double, Reaction_type);
  virtual void cal_propensity(double);

  virtual void propensity_update(int, double);  
  virtual void propensity_update(int, int, double);  
  virtual int population_update(double, double);  

 private: 
   Species &reactant2, &product2;
   int *bin_prop;
   int total_bin_prop;;
};

#endif
