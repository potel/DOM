#include "standard.h"
#include "imaginaryForm.h"
#include <iostream>

using namespace std;

int main()
{
  double Efermi = -5.65;
  standard TG(20.649,11.0264,.05,4.,Efermi);
  imaginaryForm form(20.649,11.0264,.05,0.,Efermi,4,0,0.,0.);
  for (int i=0;i<400;i++)
    {
      double E = (double)(i-200);
      double one = TG.funct(E);
      double two = TG.deltaV(E);
      double three = TG.derDeltaV(E);
      double four;
      if (E < Efermi) four = TG.derDeltaVHole(E);
      else four = TG.derDeltaVParticle(E);


      double One = form.ImaginaryPotential(E);
      double Two = form.DeltaRealPotential(E);
      double Three = form.DerDeltaHole(E)+form.DerDeltaParticle(E);
      double Four;
      if (E < Efermi) Four = form.DerDeltaHole(E);
      else Four = Four = form.DerDeltaParticle(E);

      cout << E << " " << four << " " << Four << endl;

    }

}
