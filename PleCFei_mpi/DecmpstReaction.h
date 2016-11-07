#ifndef DECMPSTREACTION_H
#define DECMPSTREACTION_H

#include "Species.h"
#include "Reaction.h"

class DecmpstReaction : public Reaction{
 public:
  DecmpstReaction(Species&, Species&, Species&, double, Reaction_type);
  virtual void cal_propensity(double);

  virtual void propensity_update(int, double);  
  virtual void propensity_update(int, int, double);  
  virtual int population_update(double, double);  

 private: 
   Species &product2;

};

#endif

