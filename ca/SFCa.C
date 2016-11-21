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
Sty->SetPadBottomMargin(.16);
Sty->SetPadLeftMargin(.14);
Sty->SetPadTopMargin(.05);
Sty->SetPadRightMargin(.05);
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
Sty->SetTitleOffset(1.2,"y");
Sty->SetTitleOffset(1.2,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.025,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

TCanvas can("SFCa");
  TH2S frame("frame","",10,38,50,10,.2,1.);
  frame.GetXaxis()->SetTitle("A");
  frame.GetYaxis()->SetTitle("S/S_{IPM}");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.Draw();

  double xn[2]={40,48};
  double yn[2]={.719,.702};
  TGraph gn(2,xn,yn);
  gn.SetMarkerStyle(21);
  gn.SetMarkerColor(4);
  gn.SetLineColor(4);
  //gn.Draw("L");


  double xp[4]={40,42,44,48};
  double yp[4]={.722,.594,.586,.625};
  TGraph gp(4,xp,yp);
  gp.SetMarkerStyle(20);
  gp.SetMarkerColor(2);
  gp.SetMarkerSize(3);
  gp.SetLineColor(2);
  gp.Draw("L");


  double xpE[2]={40,48};
  double ypE[2]={.645,.565};
  double sypE[2]={.048,.040};
  double sxpE[2]={0.};

  TGraphErrors gpE(2,xpE,ypE,sxpE,sypE);
  gpE.SetMarkerStyle(20);
  gpE.SetMarkerColor(2);
  gpE.SetMarkerSize(2);
  gpE.SetLineColor(2);
  

  gpE.Draw("P");
  

  TLatex text;
  text.SetTextColor(2);
  text.SetNDC();
  //text.DrawLatex(.5,.46,"p (0d#frac{3}{2})");
  //text.DrawLatex(.78,.38,"(e,e'p)");
  text.SetTextColor(4);
  //text.DrawLatex(.5,.68,"n");
  //text.DrawLatex(.3,.68,"(0d#frac{3}{2})");
  //text.DrawLatex(.7,.68,"(0f#frac{7}{2})");
  text.SetTextColor(1);
  text.DrawLatex(.2,.25,"Ca isotopes");
}
