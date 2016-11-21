#include "twoGauss.h"
#include <iostream>

using namespace std;

int main()
{
  double Efermi = -15.;
  twoGauss TG(1.,10.,10.,15.,4,Efermi);

  for (int i=0;i<200;i++)
    {
      double E = (double)(i-100);
      double one = TG.funct(E);
      double two = TG.deltaV(E);
      double three = TG.derDeltaV(E);
      double four;
      if (E < Efermi) four = TG.derDeltaVHole(E);
      else four = TG.derDeltaVParticle(E);
      cout << E << " " << one << " " << two << 
	" " << three << " " << four << endl;
    }

}
