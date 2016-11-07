#ifndef REACTION_H
#define REACTION_H

#include "Species.h"

enum Reaction_type{
  zero, 
  syn, deg, first_order, decomposition, 
  sticky, catalytic, 
  displacement, combination, 
  hill,
  diffusion_lft, diffusion_rgt
};  

class Reaction{

 protected:
  Species &reactant, &product;
  double rate;
  Reaction_type type;
  double propensity;

 public:
  Reaction(Species&, Species&, double, Reaction_type);
  double get_propensity();
  Reaction_type getReactionType();

  double getReactionRate();
  void setReactionRate(double);

  virtual void cal_propensity(double);
  virtual void propensity_update(int, double);
  virtual void propensity_update(int, int, double);

  virtual int population_update(double, double);  
};

#endif

