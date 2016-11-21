#include "compound.h"

int main()
{
  
  string name("nca40");
  // initial and readin information for 40Ca
  compound CN(&name);


  //specify level in 41Ca
  level lev(1.5,-1,5.);  //decaying level j,parity,excitation energy

  //find all levels in 40Ca for which decy is possible
  CN.decayElastic(lev);
}
