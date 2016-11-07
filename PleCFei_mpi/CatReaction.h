#ifndef CATREACTION_H
#define CATREACTION_H

#include "Species.h"
#include "Reaction.h"

class CatReaction : public Reaction{
 public:
  CatReaction(Species&, Species&, Species&, double, Reaction_type);
  virtual void cal_propensity(double);
  virtual void propensity_update(int, double);  
  virtual void propensity_update(int, int, double);  
  virtual int population_update(double, double);  

 private: 
   Species& aux;
   int *bin_prop;
   int total_bin_prop;

};

#endif
