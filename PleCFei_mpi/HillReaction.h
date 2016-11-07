#ifndef HILLREACTION_H
#define HILLREACTION_H

#include "Species.h"
#include "Reaction.h"

class HillReaction : public Reaction {

 public:
  HillReaction(Species&, Species&, Species&, double, Reaction_type, double);

  virtual int population_update(double, double);  

  virtual void propensity_update(int, double);  
  virtual void propensity_update(int, int, double);  
  virtual void cal_propensity(double);

 private:
  Species &aux;
  double KmdlInPopulation(double);
  double epsilon, Kmdl;
  double *bin_prop;
  double total_bin_prop;
};

#endif


