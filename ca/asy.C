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
Sty->SetPadBottomMargin(.15);
Sty->SetPadLeftMargin(.15);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.05,"xyz");
Sty->SetLabelOffset(0.01,"y");
Sty->SetLabelOffset(0.01,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.05,"xyz");
Sty->SetTitleOffset(1.,"y");
Sty->SetTitleOffset(1.,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.07,"xz");
Sty->SetTickLength(.05,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("asyCa","",800,800);

  double x[4];
  x[0] = (40.- 40.)/40.;
  x[1] = (42.- 40.)/42.;
  x[2] = (44.- 40.)/44.;
  x[3] = (48.- 40.)/48.;




  TH2S frame("frame","",10,-0.01,.22,10,0,15);
  frame.GetXaxis()->SetTitle("(N-Z)/A");
  frame.GetYaxis()->SetTitle("peak mag of Surface [MeV]");
  frame.Draw();

  double ywithXnp[4] = {6.72,8.96,10.59,9.24};
  TGraph gwithXnp(4,x,ywithXnp);
  gwithXnp.SetMarkerStyle(21);
  gwithXnp.SetMarkerSize(2);
  gwithXnp.Draw("P");



  double ywithXp[4]={6.64,8.72,10.52,8.85};

  TGraph gwithXp(4,x,ywithXp);
  gwithXp.SetMarkerStyle(25);
  gwithXp.SetMarkerSize(2);
  gwithXp.SetMarkerColor(2);
  gwithXp.Draw("P");


  double ywithoutXp[4]={6.7,8.76,10.36,8.94};

  TGraph gwithoutXp(4,x,ywithoutXp);
  gwithoutXp.SetMarkerStyle(20);
  gwithoutXp.SetMarkerSize(2);
  gwithoutXp.SetMarkerColor(3);
  gwithoutXp.Draw("P");
}
