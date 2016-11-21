#include <iostream>
#include <fstream>
#include "scat.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2S.h"
#include <fstream>

using namespace std;

int main()
{
  int proton = 1;
  
  double Z = 19.;
  double A = 40.;
  double ZZ;
  if (proton == 1) ZZ = Z;
  else ZZ = 0.;

  int const NstepsE = 200;
  int const NstepsL = 11;
  double Ecmarray[NstepsE];
  double Elabarray[NstepsE];

  double potR[NstepsE];
  double potS[NstepsE];
  double potV[NstepsE];
  double potR2[10];
  double potS2[10];
  double potV2[10];
  double Elabarray2[10];
  struct tt
  {
    double up[NstepsE];
    double down[NstepsE];
  };

  tt transm[NstepsL+1];
  double xsec[NstepsE];
  double xsec2[10];


  scat scatter(ZZ,A);

  //determine sign for proton or neutron
  double sign;
  if (proton) sign = 1.;
  else sign = -1.;

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



  ofstream FileUp("up.tl");
  ofstream FileDown("down.tl");

  FileUp << " " << NstepsL +1 << " " << NstepsE << endl;
  FileDown << " " << NstepsL + 1 << " " << NstepsE << endl;

  for (int i=0;i<NstepsE;i++)
    {
      double Ecm = ((double)i+0.5)/10.;
      double Elab = Ecm*(A+1.)/A;
      Ecmarray[i] = Ecm;
      Elabarray[i] = Elab;

      double V = V0 + sign*(A-2.*Z)/A + (Elab-Ec)*Ve;
      VHF.init(V,R0,a0);

      double Wv = Wv0/(1.+exp((Wve0-(Elab-Ec))/Wvew));
      Wvolume.init(Wv,Rw,aw);

      double Ws = (Ws0+sign*Wst*(A-2.*Z)/A)/(1.+exp(((Elab-Ec)-Wse0)/Wsew));
      Wsurface.init(Ws,Rw,aw);

      potR[i] = V;
      potV[i] = Wv;
      potS[i] = Ws;
      cout << Ecm << " " << V << " " << Ws << " " << Wv << endl;

      scatter.loadPotential(&VHF,&Vvolume,&Vsurface,&Wvolume,
                             &Wsurface,&Vspino,Rc,Ecm,Elab);

      scatter.integrateWave();
     if (i > 1) xsec[i] = scatter.AbsorptionXsection();
     else xsec[i] = 0.;
     //cout << Elabarray[i] << " " << xsec[i] << endl;

       FileUp << Ecm << " ";
       FileDown << Ecm << " " ;
       for (int jj=0;jj<=NstepsL;jj++)
           {
	     double ttup =  scatter.TransCoef(jj,(double)jj+.5);
	     double ttdown =  scatter.TransCoef(jj,(double)jj-.5);
	     FileUp <<ttup << " " ;
	     FileDown <<ttdown << " ";  
             transm[jj].up[i] = ttup;
	     transm[jj].down[i] = ttdown;
           }
       FileUp << endl;
       FileDown << endl;      
      
    }
 


  for (int i=0;i<10;i++)
    {
      double Ecm = 20. + (double)i*30./10.;
      double Elab = Ecm*(A+1.)/A;

      Elabarray2[i] = Elab;

      double V = V0 + sign*(A-2.*Z)/A + (Elab-Ec)*Ve;
      VHF.init(V,R0,a0);

      double Wv = Wv0/(1.+exp((Wve0-(Elab-Ec))/Wvew));
      Wvolume.init(Wv,Rw,aw);

      double Ws = (Ws0+sign*Wst*(A-2.*Z)/A)/(1.+exp(((Elab-Ec)-Wse0)/Wsew));
      Wsurface.init(Ws,Rw,aw);

      potR2[i] = V;
      potV2[i] = Wv;
      potS2[i] = Ws;

      scatter.loadPotential(&VHF,&Vvolume,&Vsurface,&Wvolume,
                             &Wsurface,&Vspino,Rc,Ecm,Elab);

      scatter.integrateWave();
       xsec2[i] = scatter.AbsorptionXsection();
    }
 


  TFile f("CH89.root","RECREATE");

  TCanvas absx("absx");
  TH2S hist1("hist1","CH89",10,0,50,10,0,1300);
  hist1.SetStats(kFALSE);
  hist1.GetXaxis()->SetTitle("E_{lab} [MeV]");
  hist1.GetYaxis()->SetTitle("#sigma_{abs} [mb]");
  hist1.Draw();


  TGraph *gxsec = new TGraph(NstepsE,Elabarray,xsec);
  gxsec->Draw("C");

  TGraph *gxsec2 = new TGraph(10,Elabarray2,xsec2);
  gxsec2->Draw("C");

  if (A == 40 && proton)
    {
      double energy[23];
      double sig[23];
      double error[23];
      double zero[23];
      ifstream File("ca40_reaction.xsec");
      string line;
      getline(File,line);
      getline(File,line);
      int N;
      File >> N;
      for (int i=0;i<N;i++)
	{
	  File >> energy[i] >> sig[i] >> error[i];
	  zero[i] = 0.;
	}
      File.close();
      TGraphErrors *exp = new TGraphErrors(N,energy,sig,zero,error);
      exp->SetMarkerStyle(21);
      exp->Draw("P");
    }


  absx.Write();
  //-------------------------------------------------
 TCanvas Ctl("trans");

  TH2S histTr("histTr","CH89",10,0,20,10,0,1);
  histTr.SetStats(kFALSE);
  histTr.GetXaxis()->SetTitle("Ek_{cm} [MeV]");
  histTr.GetYaxis()->SetTitle("Transmission Coefficient");
  histTr.Draw();
  
  TGraph *graphUP[NstepsL];
  TGraph *graphDOWN[NstepsL];
  for (int jj=0;jj<NstepsL;jj++)
    {
      
      graphUP[jj] = new TGraph(NstepsE,Ecmarray,transm[jj].up);
      graphUP[jj]->SetLineStyle(1);
      graphUP[jj]->SetLineColor(jj+1);
      graphUP[jj]->SetLineWidth(2);
      graphUP[jj]->Draw("C");
      
      if (jj >= 1)
	{

      graphDOWN[jj] = new TGraph(NstepsE,Ecmarray,transm[jj].down);
      graphDOWN[jj]->SetLineStyle(2);
      graphDOWN[jj]->SetLineColor(jj+1);
      graphDOWN[jj]->SetLineWidth(2);
      graphDOWN[jj]->Draw("C");

	}
      
    }
 
  Ctl.Write();

  TCanvas pot("pot");
  TH2S histp("histp","CH89",10,0,50,10,0,60);
  histp.SetStats(kFALSE);
  histp.Draw();

  TGraph *gpot = new TGraph(NstepsE,Elabarray,potR);
  gpot->SetLineWidth(2);
  gpot->Draw("C");
  gpot->DrawGraph(10,Elabarray2,potR2,"C");

  TGraph *gpots = new TGraph(NstepsE,Elabarray,potS);
  gpots->SetLineColor(2);
  gpots->SetLineWidth(2);
  gpots->Draw("C");
  gpots->DrawGraph(10,Elabarray2,potS2,"C");

  TGraph *gpotv = new TGraph(NstepsE,Elabarray,potV);
  gpotv->SetLineColor(3);
  gpotv->SetLineWidth(2);
  gpotv->Draw("C");
  gpotv->DrawGraph(10,Elabarray2,potV2,"C");
  


  pot.Write();



  if (proton != 1) return 1;

  TCanvas elas("elast");

  TH2S hist3("hist3","CH89",10,0,180,10,0.,5.6);
  hist3.SetStats(kFALSE);
  hist3.GetXaxis()->SetTitle("#theta_{cm} [deg]");
  hist3.GetYaxis()->SetTitle("#sigma/#sigma_{R}");

  hist3.Draw();
  double Ecm = 40.;
  double Elab = Ecm*(A+1.)/A;

  double V = V0 + sign*(A-2.*Z)/A + (Elab-Ec)*Ve;
  VHF.init(V,R0,a0);

  double Wv = Wv0/(1.+exp((Wve0-(Elab-Ec))/Wvew));
  Wvolume.init(Wv,Rw,aw);

  double Ws = (Ws0+sign*Wst*(A-2.*Z)/A)/(1.+exp(((Elab-Ec)-Wse0)/Wsew));
  Wsurface.init(Ws,Rw,aw);

   scatter.loadPotential(&VHF,&Vvolume,&Vsurface,&Wvolume,
                             &Wsurface,&Vspino,Rc,Ecm,Elab);

   scatter.integrateWave();
   
   double angle[180];
   double ratio[180];
   for (int i=0;i<178;i++)
     {
       double theta = 1. + (double)i;
       angle[i] = theta;
       theta *= pi/180.;
       ratio[i] =  scatter.DifferentialXsection(theta);
       //cout << ratio[i] << endl;
       ratio[i] /= scatter.Rutherford(theta);
     }
 
   TGraph *gelast = new TGraph(178,angle,ratio);
   gelast->Draw("C");
   

  elas.Write();

 
  f.Write();
  f.Close();

  return 1;
}
