#include <iostream>
#include <fstream>
#include "scatter.h"

using namespace std;

int main()
{
  int proton = 0;
  
  double Z = 29.;
  double A = 64.;
  double Zp;
  if (proton == 1) Zp = 1.;
  else Zp = 0.;

  scat scatter(Zp,Z,A,1);

  //determine sign for proton or neutron
  double sign;
  if (proton) sign = 1.;
  else sign = -1.;

  harteeFock HartreeFock;
  volume Volume;
  surface isoscalerSurface;
  surface isovectorSurface;
  spinOrbit SpinOrbit;

  PotPara VHF;
  PotPara Vvolume;
  PotPara Vsurface;
  PotPara Wvolume;
  PotPara Wsurface;
  PotPara Vspino;

  //zero the real volume
  Vvolume.init(0.,1.,1.);
  //zero the real surface
  Vsurface.init(0.,1.,1.);


  double const Efermi;
  double const r0 = 1.250;
  double const r00 = -.225;

  double const a0 = .690;
  double const rc = 1.24;
  double const rc0 = 0.12;

  double const Vso = 5.9;
  double const rso = 1.34;
  double const rso0 = -1.2;
  double const aso = 0.63;

  double const V0 = 52.9;
  double const Vt = 13.1;
  double const Ve = -.299;
  double const Wv0 = 7.8;
  double const Wve0 = 35.;
  double const Wvew = 16.;
  double const Ws0 = 10.0;
  double const Wst = 18.0;
  double const Wse0 = 36.;
  double const Wsew = 37.;

  double rw = 1.33;
  double rw0 = -.42;
  double aw = .69;

  double R0 = r0*pow(A,1./3.) + r00;
  double Rc = rc*pow(A,1./3.) + rc0;
  double Rso = rso*pow(A,1./3.) + rso0;
  double Rw = rw*pow(A,1./3.) + rw0;


  double Ec;
  if (proton) Ec = 1.73*Z/Rc;
  else Ec = 0.;

  Vspino.init(Vso,Rso,aso);


  int NstepsE = 200;
  int NstepsL = 11;
  ofstream FileUp("up.tl");
  ofstream FileDown("down.tl");

  FileUp << " " << NstepsL +1 << " " << NstepsE << endl;
  FileDown << " " << NstepsL + 1 << " " << NstepsE << endl;

  for (int i=0;i<NstepsE;i++)
    {
      double Ecm = ((double)i+0.5)/10.;
      double Elab = Ecm*(A+1.)/A;

      double V = V0 + sign*(A-2.*Z)/A + (Elab-Ec)*Ve;
      HartreeFock.load(R0,a0,V,0.,0.,0.,0.,0.,0.);



      double Wv = Wv0/(1.+exp((Wve0-(Elab-Ec))/Wvew));
      Wvolume.init(Wv,Rw,aw);

      double Ws = (Ws0+sign*Wst*(A-2.*Z)/A)/(1.+exp(((Elab-Ec)-Wse0)/Wsew));
      Wsurface.init(Ws,Rw,aw);

      scatter.loadPotential(&VHF,&Vvolume,&Vsurface,&Wvolume,
                             &Wsurface,&Vspino,Rc,Ecm,Elab);

      scatter.integrateWave();

       FileUp << Ecm << " ";
       FileDown << Ecm << " " ;
       for (int jj=0;jj<=NstepsL;jj++)
           {
	     FileUp << scatter.TransCoef(jj,(double)jj+.5)<< " " ;
	     FileDown << scatter.TransCoef(jj,(double)jj-.5)<< " ";  
           }
       FileUp << endl;
       FileDown << endl;      
      
    }
 
  return 1;
}
