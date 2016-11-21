#include "reaction.h"
#include <iostream>

using namespace std;


int main()
{
  
  double Elab = 24.;
  
  reaction Reaction;

  Reaction.Zp = 1;
  Reaction.Z = 70;
  Reaction.A = 160;
  Reaction.DOM = 0;
  string *title = new string (" ");
  bool bool1 = 0;
  Reaction.scatter = new scat(Reaction.Zp,Reaction.Z,Reaction.A,bool1,title);
  Reaction.scatter->initIntegration();
  Reaction.Rc = 1.25*pow(Reaction.A,1./3.);
  Reaction.VHF = 53.3 + -.55*Elab+ 0.4*Reaction.Z/pow(Reaction.A,(1./3.))+
       27.*(Reaction.A-2.*Reaction.Z)/Reaction.A;
  Reaction.RHF = 1.25*pow(Reaction.A,1./3.);
  Reaction.aHF = 0.65;
  Reaction.Avolume = 3.*pow(Reaction.A,1./3.);
  Reaction.Rvolume = 1.25*pow(Reaction.A,1./3.);
  Reaction.avolume = 0.47;
  Reaction.Asurface = 0.;
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
