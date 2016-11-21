#include "reaction.h"
#include <iostream>

using namespace std;


int main()
{
  
  double Elab = 80.;
  
  reaction Reaction;

  Reaction.Zp = 0;
  Reaction.Z = 4;
  Reaction.A = 9;
  Reaction.DOM = 0;
  string *title = new string (" ");
  bool bool1 = 0;
  Reaction.scatter = new scat(Reaction.Zp,Reaction.Z,Reaction.A,bool1,title);
  Reaction.scatter->initIntegration();
  Reaction.Rc = 1.3*pow(Reaction.A,1./3.);
  Reaction.VHF = 38.5 - 0.145*Elab;
  Reaction.RHF = (1.447- 0.005*(Elab-20.))*pow(Reaction.A,1./3.);
  Reaction.aHF = 0.387;
  Reaction.Avolume = 7.5-0.02*(Elab-40.);
  Reaction.Rvolume = 1.368*pow(Reaction.A,1./3.);
  Reaction.avolume = 0.3;
  Reaction.Asurface = 16.226-0.1*(Elab-40.);
  Reaction.Rsurface = Reaction.Rvolume;
  Reaction.asurface = Reaction.avolume;
  Reaction.Vso = 0.;
  Reaction.Rso = Reaction.RHF;
  Reaction.aso = Reaction.aHF;
  Reaction.AWso = 0.;

  double Ecm = Reaction.energyLab2Cm(Elab);
  Reaction.loadOM();
  Reaction.InitializeForEcm(Ecm,Elab);
  Reaction.scatter->integrateWave();

  for (int i=0;i<10;i++)
    {
      double b = ((double)i+0.5)/Reaction.scatter->Kwave;
      Reaction.scatter->getSmatrixEikonal(b);
      cout << i << " " << Reaction.scatter->TransCoef(i,(double)i+0.5)
	   <<  " " << Reaction.scatter->getSmatrixEikonal(b) << endl;
    }
  

}
