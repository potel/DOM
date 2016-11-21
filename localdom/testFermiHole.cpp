#include "fermiHole.h"
#include <iostream>

using namespace std;

int main()
{
  double Efermi = -15.;
  fermiHole FH(1.,10.,20.,5.,4,Efermi);

  /*
  double e = 0.;
  cout << (FH.funct(e+.05)-FH.funct(e-.05))/.01 << FH.derFunct(e) << endl;
  e = 20.;
  cout << (FH.funct(e+.05)-FH.funct(e-.05))/.01 << FH.derFunct(e) << endl;
  e = -25.;
  cout << (FH.funct(e+.05)-FH.funct(e-.05))/.01 << FH.derFunct(e) << endl;
  e = -35.;
  cout << (FH.funct(e+.05)-FH.funct(e-.05))/.01 << FH.derFunct(e) << endl;
  */

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
      cout << E << " " << one << " " << two << 
	" " << three << " " << four << " " << five << " " << six << endl;
    }

}
