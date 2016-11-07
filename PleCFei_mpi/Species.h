#ifndef SPECIES_H
#define SPECIES_H

class Species{

 private:
  int bin_num;
 public:
  int* population;
  int total_population;

  Species();
  Species(int);

  int size(); 
  void setInit();
  void uniFill();
  int randBinSelect(double);

};

#endif

