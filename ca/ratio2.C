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
Sty->SetPadBottomMargin(.18);
Sty->SetPadLeftMargin(.18);
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

 TCanvas canvas("ratio");

 canvas->SetLogx();
  TH2S * frame = new TH2S("frame","",10,5,300,10,-.08,.12);
  frame->GetXaxis()->SetTitle("E_{lab}");
  frame->GetYaxis()->
         SetTitle("(#sigma_{48}-#sigma_{40})/(#sigma_{48}+#sigma_{40})");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->CenterTitle();
  frame->Draw();


  ifstream file("nca40_48.data");
  string out;
  getline(file,out);
  int Nexp;
  file >> Nexp;

  double xexp[Nexp];
  double yexp[Nexp];
  double sexp[Nexp];
  double sxexp[Nexp];

  for (int i=0;i<Nexp;i++)
    {
      file >> xexp[i] >> yexp[i] >> sexp[i];
      sxexp[i] = 0.;
      //yexp[i] /= 200.;
      //sexp[i] /= 200.;
    }

  TGraphErrors gexp(Nexp,xexp,yexp,sxexp,sexp);
  gexp.SetMarkerStyle(20);
  gexp.Draw("P");

  file.close();
  file.clear();


  double xfit[200];
  double yfit[200];
  double yalpha[200];
  double yvso[200];
  double yaHF[200];
  double yalphaoverA[200];
  double yV[200];
  double yalphaNZ[200];
  double yvsoNZ[200];
  int Nfit = 0;
  file.open("ratioMin.dat");
  ifstream filealpha("ratioalpha.dat");
  ifstream filevso("ratiovso.dat");
  ifstream fileaHF("ratioaHF.dat");
  ifstream filealphaoverA("ratioalphaoverA.dat");
  ifstream fileV("ratioV.dat");
  ifstream filealphaNZ("ratioalphaNZ.dat");
  ifstream filevsoNZ("ratiovsoNZ.dat");
  float crap;
  for (;;)
    {
      file >> xfit[Nfit] >> yfit[Nfit];
      filealpha >> crap >> yalpha[Nfit];
      filevso >> crap >> yvso[Nfit];
      fileaHF >> crap >> yaHF[Nfit];
      filealphaoverA >> crap >> yalphaoverA[Nfit];
      fileV >> crap >> yV[Nfit];
      filealphaNZ >> crap >> yalphaNZ[Nfit];
      filevsoNZ >> crap >> yvsoNZ[Nfit];
      if (file.eof())break;
      if (file.bad())break;
      Nfit++;
    }
 
  TGraph gfit(Nfit,xfit,yfit);
  gfit.Draw("C");

  //20% increase in alpha for both
  TGraph galpha(Nfit,xfit,yalpha);
  galpha.SetLineColor(2);
  galpha.Draw("C");


  //set vso to zero from 6.
  TGraph gvso(Nfit,xfit,yvso);
  gvso.SetLineColor(3);
  gvso.Draw("C");

  //20% increase in aHF
  TGraph gaHF(Nfit,xfit,yaHF);
  gaHF.SetLineColor(4);
  gaHF.Draw("C");

  //alphaoverA set to zero from .08
  TGraph galphaoverA(Nfit,xfit,yalphaoverA);
  galphaoverA.SetLineColor(5);
  galphaoverA.Draw("C");


  //set volume1 to zero from 4.22, no asy dependence of volume
  TGraph gV(Nfit,xfit,yV);
  gV.SetLineColor(6);
  //gV.Draw("C");


  //set alphaNZ to zero from .17, no asy dependence of eff mass
  TGraph galphaNZ(Nfit,xfit,yalphaNZ);
  galphaNZ.SetLineColor(7);
  //galphaNZ.Draw("C");


  //set vsoNZ to zero from -1., no asy dependence of sin orbit
  TGraph gvsoNZ(Nfit,xfit,yvsoNZ);
  gvsoNZ.SetLineColor(8);
  //gvsoNZ.Draw("C");



  file.close();
  file.clear();

  TLine line;
  line.SetLineStyle(2);
  line.DrawLine(5.,.06069,300,.06069);
  return;
 
  double xfit[200];
  double yfit[200];
  int Nfit = 0;
  file.open("ratiofi.dat");
  for (;;)
    {
      file >> xfit[Nfit] >> yfit[Nfit];
      if (file.eof())break;
      if (file.bad())break;
      Nfit++;
    }
 
  TGraph gfit(Nfit,xfit,yfit);
  gfit.Draw("C");


  file.close();
  file.clear();
 
  double x0[200];
  double y0[200];
  int N0 = 0;
  file.open("ratio0.dat");
  for (;;)
    {
      file >> x0[N0] >> y0[N0];
      if (file.eof())break;
      if (file.bad())break;
      N0++;
    }
 
  TGraph g0(N0,x0,y0);
  g0.SetLineColor(2);
  g0.SetLineStyle(2);
  g0.Draw("C");
  


  file.close();
  file.clear();
 
  double xL[200];
  double yL[200];
  int NL = 0;
  file.open("ratioL.dat");
  for (;;)
    {
      file >> xL[NL] >> yL[NL];
      if (file.eof())break;
      if (file.bad())break;
      NL++;
    }
 
  TGraph gL(NL,xL,yL);
  gL.SetLineColor(4);
  gL.SetLineStyle(8);
  gL.Draw("C");

  TLatex text;
  text.SetNDC();
  text.DrawLatex(.53,.4,"fit to elastic d#sigma/d#Omega data");
  text.SetTextColor(2);
  text.DrawLatex(.53,.35,"zero");
  text.SetTextColor(4);
  text.DrawLatex(.53,.3,"opposite");  



}

