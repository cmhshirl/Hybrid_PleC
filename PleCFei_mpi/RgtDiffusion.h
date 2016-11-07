#ifndef RGTDIFFUSION_H
#define RGTDIFFUSION_H

#include "Species.h"
#include "Reaction.h"

class RgtDiffusion : public Reaction{
 public:
  RgtDiffusion(Species&, Species&, double, Reaction_type, int);
  virtual void cal_propensity(double);

  virtual void propensity_update(int, double);  
  virtual void propensity_update(int, int, double);  
  virtual int population_update(double, double);  

  void setBarrier(int);
 private:
  int barrier;
  int dif_population;
};

#endif

