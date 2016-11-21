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
Sty->SetLabelSize(0.09,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.09,"xyz");
Sty->SetTitleOffset(1.,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("ImagAsyCa");


//divide up canvas into three pads
 double overlap = .1;
 double dist = (1.+3*overlap)/4.;
 TPad *pad1  = new TPad("pad1","",0.,0.,1.,.5+overlap);
 TPad *pad2  = new TPad("pad2","",0.,.5-overlap,1.,1.);
pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

gPad->GetFrame()->SetLineWidth(1);
pad1->Draw();
pad2->Draw();


 pad1->cd();
  TH2S frame("frame","",10,-200,200,10,0,19.9);
  frame.GetXaxis()->SetTitle("E-E_{F} [MeV]");
  frame.GetYaxis()->SetTitle("W [MeV]");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.SetStats(kFALSE);
  frame.Draw();


 pad2->cd();
  TH2S framev("framev","",10,-200,200,10,0,19.9);
  //framev.GetXaxis()->SetTitle("E [MeV]");
  framev.GetXaxis()->SetLabelSize(0.);
  framev.GetYaxis()->SetTitle("W [MeV]");
  framev.SetStats(kFALSE);
  framev.Draw();

  ifstream filepca40("pca40.poten");
  ifstream filepca48("pca48.poten");
  ifstream filenca40("nca40.poten");
  ifstream filenca48("nca48.poten");


  double xpca40[150];
  double surpca40scaler[150];
  double surpca40vector[150];
  double surpca40[150];
  double volpca40[150];
  double HFpca40[150];

  double xpca48[150];
  double surpca48scaler[150];
  double surpca48vector[150];
  double surpca48[150];
  double volpca48[150];
  double HFpca48[150];


  double xnca40[150];
  double surnca40scaler[150];
  double surnca40vector[150];
  double surnca40[150];
  double volnca40[150];
  double HFnca40[150];

  double xnca48[150];
  double surnca48scaler[150];
  double surnca48vector[150];
  double surnca48[150];
  double volnca48[150];
  double HFnca48[150];






  double zero,one,two,three,four, five,six,seven;

for (int i=0;i<150;i++) 
{

  filepca40 >> xpca40[i] >> surpca40scaler[i] >> surpca40vector[i] >> 
    volpca40[i] >> HFpca40[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surpca40[i] = surpca40scaler[i] + surpca40vector[i];

  filepca48 >> xpca48[i] >> surpca48scaler[i] >> surpca48vector[i] >> 
    volpca48[i] >> HFpca48[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surpca48[i] = surpca48scaler[i] + surpca48vector[i];


  filenca40 >> xnca40[i] >> surnca40scaler[i] >> surnca40vector[i] >> 
    volnca40[i] >> HFnca40[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surnca40[i] = surnca40scaler[i] + surnca40vector[i];

  filenca48 >> xnca48[i] >> surnca48scaler[i] >> surnca48vector[i] >> 
    volnca48[i] >> HFnca48[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surnca48[i] = surnca48scaler[i] + surnca48vector[i];



}
 pad1->cd();

 TGraph gpca40(150,xpca40,surpca40scaler);
 gpca40.SetLineColor(2);
 gpca40.Draw("L");


 pad2->cd();

 TGraph gpca48(150,xpca48,surpca48scaler);
 gpca48.SetLineColor(2);
 gpca48.Draw("L");


 pad1->cd();
 TGraph gnca40(150,xnca40,surnca40scaler);
 gnca40.SetLineColor(4);
 gnca40.Draw("L");

 pad2->cd();

 TGraph gnca48(150,xnca48,surnca48scaler);
 gnca48.SetLineColor(4);
 gnca48.Draw("L");




 pad1->cd();
 TLatex text;
 text.SetTextSize(.1);
 text.SetNDC();
 text.SetTextColor(1);
 text.DrawLatex(.25,.75,"^{40}Ca");
 pad2->cd();
 text.DrawLatex(.25,.75,"^{48}Ca");
 pad1->cd();
 text.SetTextSize(.08);
  text.DrawLatex(.75,.5,"volume");
  text.DrawLatex(.5,.6,"surface");
  pad2->cd();
  text.SetTextColor(2);
  text.DrawLatex(.59,.6,"p");
  text.SetTextColor(4);
  text.DrawLatex(.59,.39,"n");
  

 pad1->cd();

 TGraph gpca40vol(150,xpca40,volpca40);
 gpca40vol.SetLineColor(2);
 gpca40vol.Draw("L"); 

 pad2->cd();

 TGraph gpca48vol(150,xpca48,volpca48);
 gpca48vol.SetLineColor(2);
 gpca48vol.Draw("L");

 pad1->cd();
 TGraph gnca40vol(150,xnca40,volnca40);
 gnca40vol.SetLineColor(4);
 gnca40vol.Draw("L"); 

 pad2->cd();

 TGraph gnca48vol(150,xnca48,volnca48);
 gnca48vol.SetLineColor(4);
 gnca48vol.Draw("L");



}
