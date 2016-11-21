#include "twoFermi.h"
#include <iostream>

using namespace std;

int main()
{
  double Efermi = -13.21;
  double gap = 7.92;
  twoFermi FH(44.19,12.56,25.04,11.72,gap,Efermi);

  
  double e = 0.5;
  cout << e << " " << (FH.funct(e+.05)-FH.funct(e-.05))/.1 << " " 
       << FH.derFunct(e) << endl;
  e = 20.;
  cout << e << " " << (FH.funct(e+.05)-FH.funct(e-.05))/.1 << " " << 
      FH.derFunct(e) << endl;
  e = -25.;
  cout << e << " " << (FH.funct(e+.05)-FH.funct(e-.05))/.1 << " " 
       << FH.derFunct(e) << endl;
  e = -35.;
  cout << e << " " << (FH.funct(e+.05)-FH.funct(e-.05))/.1 << " " << 
       FH.derFunct(e) << endl;
  
  e = -17.;
  cout << e << " " << (FH.deltaV(e+.05)-FH.deltaV(e-.05))/.1 << " " << 
       FH.derDeltaV(e) << endl;


  for (int i=0;i<200;i++)
    {
      double E = (double)(i-100);
      double one = FH.funct(E);
      double two = FH.deltaV(E);
      double three = FH.derDeltaV(E);
      double four;
      if (E < Efermi) four = FH.derDeltaVHole(E);
      else four = FH.derDeltaVParticle(E);
      double five = FH.derDeltaV();
      double six = FH.derDeltaVParticleHole();

      FH.deltaVAboveBelow(E);
      double seven = FH.deltaVAbove();
      double eight = FH.deltaVBelow();
      double nine  = FH.derDeltaVAbove();
      double ten  = FH.derDeltaVBelow();
      

      cout << E << " " << one << " " << two << 
	" " << three << " " << seven << " " <<
        eight << " " << seven+eight << " " << nine << " " << ten << " " 
        << nine+ten <<endl;
    }

}
