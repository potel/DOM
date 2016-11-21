#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "TCanvas.h"
#include "TH2S.h"
//#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
//#include "gRoot.h"
using namespace std;


class line
{
public:
  float x[20];
  float y[20];
  int N;

  line();
  void add(float,float);
  void plotline(int,int);
  void plotpoints(int,int);
};

line::line()
{
  N = 0;
}

void line::add(float xx, float yy)
{
  x[N] = xx;
  y[N] = yy;
  N++;
}

void line::plotline(int icol, int istyle)
{
  if (N == 0) return;
  TGraph *graph = new TGraph(N,x,y);
  graph->SetLineColor(icol);
  graph->SetLineStyle(istyle);
  graph->Draw("L");
}

void line::plotpoints(int icol, int istyle)
{
  if (N == 0) return;
  TGraph *graph = new TGraph(N,x,y);
  graph->SetMarkerColor(icol);
  graph->SetMarkerStyle(istyle);
  graph->Draw("P");
}

int main()
{

gROOT->Reset();

TStyle * Sty = new TStyle("MyStyle","MyStyle");
Sty->SetOptTitle(0);
Sty->SetOptStat(0);
Sty->SetPalette(8,0);
Sty->SetCanvasColor(10);
Sty->SetCanvasBorderMode(0);
Sty->SetFrameLineWidth(3);
Sty->SetFrameFillColor(10);
Sty->SetPadColor(10);
Sty->SetPadTickX(1);
Sty->SetPadTickY(1);
Sty->SetPadBottomMargin(.2);
Sty->SetPadLeftMargin(.2);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.06,"xyz");
Sty->SetTitleOffset(1.4,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.025,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();



  line FermiN;
  line N0s12;
  line N0p32;
  line N0p12;
  line N0d52;
  line N1s12;
  line N0d32;
  line N0f72;
  line N0f52;
  line N1p12;
  line N1p32;
  line N0g92;
  line N1d52;

  line NN0d52;
  line NN0d32;
  line NN1s12;
  line NN0f72;
  line NN1p12;
  line NN1p32;
  line NN0f52;


  line FermiP;
  line P0s12;
  line P0p32;
  line P0p12;
  line P0d52;
  line P1s12;
  line P0d32;
  line P0f72;
  line P0f52;
  line P1p12;
  line P1p32;
  line P0g92;
  line P1d52;

  line PP0d52;
  line PP0d32;
  line PP1s12;
  line PP0f72;
  line PP1p12;
  line PP1p32;
  line PP0f52;

  float z,a,three,Efermi;
  //TFile *f = new TFile("alev.root","RECREATE");
  char stuff[100];
  ostringstream outstring;
  string name;
  ifstream file;
  int Nneut = 0;
  float j,energy,RMS,occup,spect,Delta;
  float error,RMSerror,spectError,DeltaError;
  int l,N,colour,Efit,Rfit,Dfit,sfit;

  
  // neutrons
  for (int A=36;A<74;A+=2)
    {
      outstring.str("");
      outstring << "nca" << A << ".inp";
      name = outstring.str();
      file.open(name.c_str());
      if (file.fail()) 
	{
	  cout << name << " not opened " << endl;
          file.clear();
          continue;
	}
      Nneut++;
      file >> three >> z >> a >> Efermi;
      file.close();
      file.clear();

      if (Efermi == 0.) 
	{
         outstring.str("");
         outstring << "nca" << A << ".fermi";
         name = outstring.str();
         file.open(name.c_str());
         if (file.fail()) 
	   {
	     cout << name << " not opened " << endl;
             file.clear();
             continue;
	   }
	 file >> Efermi;
	 file.close();
	 file.clear();
	}
      cout << Efermi << endl;
      FermiN.add((float)A,Efermi);

      outstring.str("");
      outstring << "nca" << A << ".level";
      name = outstring.str();
      file.open(name.c_str());
      if (file.fail()) 
	{
	  cout << name << " not opened " << endl;
          file.clear();
          continue;
	}
       file.getline(stuff,100);

       for (;;)
         {
           file >> N >> j >> l >> energy >> RMS >> occup 
                >> spect >> Delta>> colour;
           if (file.eof())break;

           if (l == 0 && j == 0.5 && N ==0 ) N0s12.add((float)A,energy);
           if (l == 0 && j == 0.5 && N ==1 ) N1s12.add((float)A,energy);
           if (l == 1 && j == 0.5 && N ==0 ) N0p12.add((float)A,energy);
           if (l == 1 && j == 1.5 && N ==0 ) N0p32.add((float)A,energy);
           if (l == 2 && j == 2.5 && N ==0 ) N0d52.add((float)A,energy);
           if (l == 2 && j == 1.5 && N ==0 ) N0d32.add((float)A,energy);
           if (l == 3 && j == 3.5 && N ==0 ) N0f72.add((float)A,energy);
           if (l == 3 && j == 2.5 && N ==0 ) N0f52.add((float)A,energy);
           if (l == 1 && j == 0.5 && N ==1 ) N1p12.add((float)A,energy);
           if (l == 1 && j == 1.5 && N ==1 ) N1p32.add((float)A,energy);
           if (l == 4 && j == 4.5 && N ==0 ) N0g92.add((float)A,energy);
           if (l == 2 && j == 2.5 && N ==1 ) N1d52.add((float)A,energy);
	 }
       file.close();
       file.clear();

     outstring.str("");
      outstring << "nca" << A << ".lev";
      name = outstring.str();
      file.open(name.c_str());
      if (file.fail()) 
	{
	  cout << name << " not opened " << endl;
	  file.clear();
          continue;
	}
       file.getline(stuff,100);
       file.getline(stuff,100);

       for (;;)
         {
            file >> energy >> error >> N >> j >> l >> colour >>Efit>> 
            RMS >> RMSerror >> Rfit >> Delta>>DeltaError >> Dfit >> 
            spect>>spectError>>sfit;

           if (file.eof())break;

           //if (l == 0 && j == 0.5 && N ==0 ) NN0s12.add((float)A,energy);
           if (l == 0 && j == 0.5 && N ==1 ) NN1s12.add((float)A,energy);
           //if (l == 1 && j == 0.5 && N ==0 ) NN0p12.add((float)A,energy);
           //if (l == 1 && j == 1.5 && N ==0 ) NN0p32.add((float)A,energy);
           if (l == 2 && j == 2.5 && N ==0 ) NN0d52.add((float)A,energy);
           if (l == 2 && j == 1.5 && N ==0 ) NN0d32.add((float)A,energy);
           if (l == 3 && j == 3.5 && N ==0 ) NN0f72.add((float)A,energy);
           if (l == 3 && j == 2.5 && N ==0 ) NN0f52.add((float)A,energy);
           //if (l == 3 && j == 2.5 && N ==0 ) NN0f52.add((float)A,energy);
           if (l == 1 && j == 0.5 && N ==1 ) NN1p12.add((float)A,energy);
           if (l == 1 && j == 1.5 && N ==1 ) NN1p32.add((float)A,energy);
           //if (l == 4 && j == 4.5 && N ==0 ) NN0g92.add((float)A,energy);
	 }
       file.close();
       file.clear();

    }
  cout << "proton" << endl;

  // protons
  for (int A=36;A<73;A+=2)
    {
      outstring.str("");
      outstring << "pca" << A << ".inp";
      name = outstring.str();
      file.open(name.c_str());
      if (file.fail()) 
	{
	  cout << name << " not opened " << endl;
          file.clear();
          continue;
	}
      Nneut++;
      file >> three >> z >> a >> Efermi;
      file.close();
      file.clear();
      if (Efermi == 0.) 
	{
         outstring.str("");
         outstring << "pca" << A << ".fermi";
         name = outstring.str();
         file.open(name.c_str());
         if (file.fail()) 
	   {
	     cout << name << " not opened " << endl;
             file.clear();
             continue;
	   }
	 file >> Efermi;
	 file.close();
	 file.clear();
	}
      cout << Efermi << endl;
      FermiP.add((float)A,Efermi);

      outstring.str("");
      outstring << "pca" << A << ".level";
      name = outstring.str();
      file.open(name.c_str());
      if (file.fail()) 
	{
	  cout << name << " not opened " << endl;
          file.clear();
          continue;
	}
       file.getline(stuff,100);

       for (;;)
         {
           file >> N >> j >> l >> energy >> RMS >> occup 
                >> spect >> Delta>> colour;
           if (file.eof())break;

           if (l == 0 && j == 0.5 && N ==0 ) P0s12.add((float)A,energy);
           if (l == 0 && j == 0.5 && N ==1 ) P1s12.add((float)A,energy);
           if (l == 1 && j == 0.5 && N ==0 ) P0p12.add((float)A,energy);
           if (l == 1 && j == 1.5 && N ==0 ) P0p32.add((float)A,energy);
           if (l == 2 && j == 2.5 && N ==0 ) P0d52.add((float)A,energy);
           if (l == 2 && j == 1.5 && N ==0 ) P0d32.add((float)A,energy);
           if (l == 3 && j == 3.5 && N ==0 ) P0f72.add((float)A,energy);
           if (l == 3 && j == 2.5 && N ==0 ) P0f52.add((float)A,energy);
           if (l == 1 && j == 0.5 && N ==1 ) P1p12.add((float)A,energy);
           if (l == 1 && j == 1.5 && N ==1 ) P1p32.add((float)A,energy);
           if (l == 4 && j == 4.5 && N ==0 ) P0g92.add((float)A,energy);
           if (l == 2 && j == 2.5 && N ==1 ) P1d52.add((float)A,energy);
	 }
       file.close();
       file.clear();
    outstring.str("");
      outstring << "pca" << A << ".lev";
      name = outstring.str();
      file.open(name.c_str());
      if (file.fail()) 
	{
	  cout << name << " not opened " << endl;
          file.clear();
          continue;
	}
       file.getline(stuff,100);
       file.getline(stuff,100);

       for (;;)
         {
            file >> energy >> error >> N >> j >> l >> colour >>Efit>> 
            RMS >> RMSerror >> Rfit >> Delta>>DeltaError >> Dfit >> 
            spect>>spectError>>sfit;

           if (file.eof())break;

           //if (l == 0 && j == 0.5 && N ==0 ) PP0s12.add((float)A,energy);
           if (l == 0 && j == 0.5 && N ==1 ) PP1s12.add((float)A,energy);
           //if (l == 1 && j == 0.5 && N ==0 ) PP0p12.add((float)A,energy);
           //if (l == 1 && j == 1.5 && N ==0 ) PP0p32.add((float)A,energy);
           if (l == 2 && j == 2.5 && N ==0 ) PP0d52.add((float)A,energy);
           if (l == 2 && j == 1.5 && N ==0 ) PP0d32.add((float)A,energy);
           if (l == 3 && j == 3.5 && N ==0 ) PP0f72.add((float)A,energy);
           if (l == 3 && j == 2.5 && N ==0 ) PP0f52.add((float)A,energy);
           //if (l == 3 && j == 2.5 && N ==0 ) PP0f52.add((float)A,energy);
           if (l == 1 && j == 0.5 && N ==1 ) PP1p12.add((float)A,energy);
           if (l == 1 && j == 1.5 && N ==1 ) PP1p32.add((float)A,energy);
           //if (l == 4 && j == 4.5 && N ==0 ) PP0g92.add((float)A,energy);
	 }
       file.close();
       file.clear();
    }

  TCanvas *canvas = new TCanvas("alev","",600,900);

  //TPad *pad1 = new TPad("pad1","",0,0,.59,1.);
  //TPad *pad2 = new TPad("pad2","",.41,0,1.,1.);
  TPad *pad1 = new TPad("pad1","",0,0.42,1,1.);
  TPad *pad2 = new TPad("pad2","",0,0.,1,.58);
  pad1->SetFillStyle(4000);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  TH2S *frame1 = new TH2S ("frame1","",10,36,72,10,-70,5);
  frame1->SetStats(kFALSE);
  //frame1->GetXaxis()->SetTitle("A");
  frame1->SetLabelSize(0.);
  frame1->GetXaxis()->CenterTitle();
  frame1->GetYaxis()->SetTitle("E_{n,l,j}[MeV]");
  frame1->GetYaxis()->CenterTitle();
  
  frame1->Draw();


  FermiN.plotline(1,2);
  N0s12.plotline(1,1);
  N0p12.plotline(2,1);
  N0p32.plotline(2,1);
  N0d52.plotline(3,1);
  N0d32.plotline(3,1);
  N1s12.plotline(1,1);
  N0f72.plotline(4,1);
  N0f52.plotline(4,1);
  N1p12.plotline(2,1);
  N1p32.plotline(2,1);
  N0g92.plotline(6,1);
  N1d52.plotline(3,1);

  NN0d52.plotpoints(3,20);
  NN0d32.plotpoints(3,21);
  NN1s12.plotpoints(1,21);
  NN0f72.plotpoints(4,21);
  NN0f52.plotpoints(4,21);
  NN1p12.plotpoints(2,20);
  NN1p32.plotpoints(2,21);
  TLatex text;
  text.DrawLatex(40,-54.,"(a) Neutron");
  text.SetTextSize(.04);
  text.DrawLatex(40,-65,"0s_{1/2}");
  text.DrawLatex(50,-21,"1s_{1/2}");

  text.SetTextColor(2);
  text.DrawLatex(38,-46,"0p_{3/2}");
  text.DrawLatex(40,-36,"0p_{1/2}");
  text.DrawLatex(57,-15.5,"1p");
  text.SetTextColor(3);
  text.DrawLatex(66,-30,"0d_{5/2}");
  text.DrawLatex(50,-14.,"0d_{3/2}");
  text.DrawLatex(64,2.,"1d_{5/2}");

  text.SetTextColor(4);
  text.DrawLatex(66,-14,"0f_{7/2}");
  text.DrawLatex(36.7,1.,"0f_{5/2}");
  text.SetTextColor(6);
  text.DrawLatex(45,2.,"0g_{9/2}");

  TArrow arrow;
  arrow.SetAngle(30.);
  arrow.SetFillColor(2);
  arrow.SetLineColor(2);
  arrow.DrawArrow(57.5,-12,57,-5,.03,"|>");
  arrow.DrawArrow(58,-12,58.5,-7,.03,"|>");

  pad2->cd();
  TH2S *frame2 = new TH2S ("frame2","",10,36,72,10,-70,5);
  frame2->SetStats(kFALSE);
  frame2->GetXaxis()->SetTitle("A");
  frame2->GetXaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitle("E_{n,l,j}[MeV]");
  frame2->GetYaxis()->CenterTitle();



  frame2->Draw();
  frame2->Draw();
  text.SetTextColor(1);
  text.SetTextSize(.05);
  text.DrawLatex(45,-54.,"(b) Proton");
  text.SetTextSize(.04);
  text.DrawLatex(40,-65,"0s_{1/2}");
  text.DrawLatex(36,-19,"1s_{1/2}");

  text.SetTextColor(2);
  text.DrawLatex(66,-60,"0p_{3/2}");
  text.DrawLatex(40,-28,"0p_{1/2}");
  text.DrawLatex(63,-8.7,"1p");
  text.SetTextColor(3);
  text.DrawLatex(66,-39,"0d_{3/2}");
  text.DrawLatex(50,-27.5,"0d_{5/2}");
  text.DrawLatex(66,-3,"1d_{5/2}");

  text.SetTextColor(4);
  text.DrawLatex(59.2,-7.7,"0f");

  text.SetTextColor(6);
  text.DrawLatex(66.4,-10,"0g_{9/2}");
 

  arrow.SetFillColor(1);
  arrow.SetLineColor(1);
  arrow.DrawArrow(38,-17.,38.,-8,.03,"|>");
  arrow.SetFillColor(2);
  arrow.SetLineColor(2);
  arrow.DrawArrow(64,-9.6,63.,-15.7,.02,"|>");
  arrow.DrawArrow(64.5,-9.6,65.5.,-14,.02,"|>");
  arrow.SetFillColor(4);
  arrow.SetLineColor(4);
  arrow.DrawArrow(59.5,-9,59,-18,.02,"|>");
  arrow.DrawArrow(60.5,-9,61.5,-14.,.02,"|>");
  arrow.SetFillColor(3);
  arrow.SetLineColor(3);
  arrow.DrawArrow(69,-38,69,-32,.02,"|>");

  FermiP.plotline(1,2);
  P0s12.plotline(1,1);
  P0p12.plotline(2,1);
  P0p32.plotline(2,1);
  P0d52.plotline(3,1);
  P0d32.plotline(3,1);
  P1s12.plotline(1,1);
  P0f72.plotline(4,1);
  P0f52.plotline(4,1);
  P1p12.plotline(2,1);
  P1p32.plotline(2,1);
  P0g92.plotline(6,1);
  P1d52.plotline(3,1);

  PP0d52.plotpoints(3,20);
  PP0d32.plotpoints(3,21);
  PP1s12.plotpoints(1,21);
  PP0f72.plotpoints(4,21);
  PP0f52.plotpoints(4,21);
  PP1p12.plotpoints(2,20);
  PP1p32.plotpoints(2,21);

  //canvas->Write();
  //f->Write();
  //f->Close();


  double xs[1] = {40};
  double ys[1] = {-50};
  double xe[1] = {0.};
  double ye[1] = {4};

  TGraphErrors *s = new TGraphErrors(1,xs,ys,xe,ye);
  s->SetMarkerStyle(21);
  s->Draw("P"); 


  double xp[1] = {40};
  double yp[1] = {-35.46};
  double xpe[1] = {0.};
  double ype[1] = {4};

  TGraphErrors *p = new TGraphErrors(1,xp,yp,xpe,ype);
  p->SetMarkerStyle(20);
  p->SetMarkerColor(2);
  p->SetLineColor(2);
  p->Draw("P"); 


}
