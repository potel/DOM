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
Sty->SetPadBottomMargin(.25);
Sty->SetPadLeftMargin(.2);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.07,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
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

 TCanvas can("asyi","",300,600);

//divide up canvas into three pads
 double overlap = .08;
 TPad *pad1  = new TPad("pad1","",0.,0.,1.,.5+overlap);
 TPad *pad2  = new TPad("pad2","",0.,.5-overlap,1.,1.);

pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

gPad->GetFrame()->SetLineWidth(1);
pad1->Draw();
pad2->Draw();


 pad1->cd();

  ifstream file("fitI5.inp");
  string out;

  for (int i=0;i<19;i++)
    {
     getline(file,out);
     cout << out << endl;
    }

  int i1,i2;
  double d1;
  double A40,Ap48,Ap42,Ap44,An48;
  double B40,Bp48,Bp42,Bp44,Bn48;
  file >> A40  >> i1 >> i2 >> d1 >> out;

  file >> B40  >> i1 >> i2 >> d1 >> out;

  file >> Ap48  >> i1 >> i2 >> d1 >> out;
  file >> Bp48  >> i1 >> i2 >> d1 >> out;
  file >> An48  >> i1 >> i2 >> d1 >> out;
  file >> Bn48  >> i1 >> i2 >> d1 >> out;
  file >> Ap42  >> i1 >> i2 >> d1 >> out;
  file >> Bp42  >> i1 >> i2 >> d1 >> out;  
  file >> Ap44  >> i1 >> i2 >> d1 >> out;
  file >> Bp44  >> i1 >> i2 >> d1 >> out;


  double x[5]={-.16666,0.,.0476,.090,.1667};

  TH2S frameA("frameA","",10,-.2,.2,10,0,40);
  frameA.Draw();


  double yA[5];
  yA[0] = An48;
  yA[1] = A40;
  yA[2] = Ap42;
  yA[3] = Ap44;
  yA[4] = Ap48;

  TGraph gA(5,x,yA);
  gA.SetMarkerStyle(21);
  gA.SetMarkerColor(2);
  gA.SetLineColor(2);
  gA.Draw("PL");


  pad2->cd();

  TH2S frameB("frameB","",10,-.2,.2,10,0,20);
  frameA.Draw();


  double yB[5];
  yB[0] = Bn48;
  yB[1] = B40;
  yB[2] = Bp42;
  yB[3] = Bp44;
  yB[4] = Bp48;

  TGraph gB(5,x,yB);
  gB.SetMarkerStyle(21);
  gB.SetMarkerColor(4);
  gB.SetLineColor(4);
  gB.Draw("PL");


}
