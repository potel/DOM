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

 TCanvas canvas("ComImag2Ca");


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
  TH2S frame("frame","",10,-200,200,10,0,14);
  frame.GetXaxis()->SetTitle("E-E_{F} [MeV]");
  frame.GetYaxis()->SetTitle("W_{surf} [MeV]");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.SetStats(kFALSE);
  frame.Draw();

  ifstream filepca48("pca48.poten");
  ifstream filepca44("pca44.poten");
  ifstream filepca42("pca42.poten");
  ifstream filenca48("nca48.poten");
  ifstream fileCa40("../ca/pca40.poten");

  double xpca48[150];
  double surpca48scaler[150];
  double surpca48vector[150];
  double surpca48[150];
  double volpca48[150];
  double HFpca48[150];
  double x[150];


  double xpca44[150];
  double surpca44scaler[150];
  double surpca44vector[150];
  double surpca44[150];
  double volpca44[150];
  double HFpca44[150];


  double xpca42[150];
  double surpca42scaler[150];
  double surpca42vector[150];
  double surpca42[150];
  double volpca42[150];
  double HFpca42[150];




  double xnca48[150];
  double surnca48scaler[150];
  double surnca48vector[150];
  double surnca48[150];
  double volnca48[150];
  double HFnca48[150];


  double xCa40[150];
  double surCa40scaler[150];
  double surCa40vector[150];
  double surCa40[150];
  double volCa40[150];
  double HFCa40[150];


  double one,two,three,four, five,six;

for (int i=0;i<150;i++) 
{
  filepca48 >> xpca48[i] >> surpca48scaler[i] >> surpca48vector[i] >> 
  volpca48[i] >> HFpca48[i] >> 
  one >> two >> three >> four >> five >> six;
 surpca48[i] = surpca48scaler[i] + surpca48vector[i];


  filepca44 >> xpca44[i] >> surpca44scaler[i] >> surpca44vector[i] >> 
  volpca44[i] >> HFpca44[i] >> 
  one >> two >> three >> four >> five >> six;
 surpca44[i] = surpca44scaler[i] + surpca44vector[i];


  filepca42 >> xpca42[i] >> surpca42scaler[i] >> surpca42vector[i] >> 
  volpca42[i] >> HFpca42[i] >> 
  one >> two >> three >> four >> five >> six;
 surpca42[i] = surpca42scaler[i] + surpca42vector[i];

  filenca48 >> xnca48[i] >> surnca48scaler[i] >> surnca48vector[i] >> 
  volnca48[i] >> HFnca48[i] >> 
  one >> two >> three >> four >> five >> six;
 surnca48[i] = surnca48scaler[i] + surnca48vector[i];


  fileCa40 >> xCa40[i] >> surCa40scaler[i] >> surCa40vector[i] >> 
  volCa40[i] >> HFCa40[i] >> 
  one >> two >> three >> four >> five >> six;
 surCa40[i] = surCa40scaler[i] + surCa40vector[i];
 // cout << x[i] << " " << surCa40scaler[i] << endl;
}

 TGraph gpca48(150,xpca48,surpca48);
 gpca48.SetLineColor(4);
 gpca48.Draw("L");


 TGraph gpca44(150,xpca44,surpca44);
 gpca44.SetLineColor(6);
 gpca44.Draw("L");

 TGraph gpca42(150,xpca42,surpca42);
 gpca42.SetLineColor(7);
 gpca42.Draw("L");

 TGraph gnca48(150,xnca48,surnca48);
 gnca48.SetLineColor(3);
 gnca48.Draw("L");
 
 TGraph gCa40(150,xCa40,surCa40);
 gCa40.SetLineColor(2);
 gCa40.Draw("L");

 TLatex text;
 text.SetTextSize(.08);
 text.SetNDC();
 text.SetTextColor(4);
 text.DrawLatex(.25,.7,"p+^{48}Ca");
 text.SetTextColor(6);
 text.DrawLatex(.25,.6,"p+^{44}Ca");
 text.SetTextColor(7);
 text.DrawLatex(.25,.5,"p+^{42}Ca");
 text.SetTextColor(2);
 text.DrawLatex(.25,.4,"n,p+^{40}Ca");
 text.SetTextColor(3);
 text.DrawLatex(.25,.3,"n+^{48}Ca");


 pad2->cd();
  TH2S framev("framev","",10,-200,200,10,0,12);
  //framev.GetXaxis()->SetTitle("E [MeV]");
  framev.GetXaxis()->SetLabelSize(0.);
  framev.GetYaxis()->SetTitle("W_{vol} [MeV]");
  framev.SetStats(kFALSE);
  framev.Draw();


 TGraph gpca48vol(150,xpca48,volpca48);
 gpca48vol.SetLineColor(4);
 gpca48vol.Draw("L");

 TGraph gpca44vol(150,xpca44,volpca44);
 gpca44vol.SetLineColor(6);
 gpca44vol.Draw("L");

 TGraph gpca42vol(150,xpca42,volpca42);
 gpca42vol.SetLineColor(7);
 gpca42vol.Draw("L");

 TGraph gnca48vol(150,xnca48,volnca48);
 gnca48vol.SetLineColor(3);
 gnca48vol.Draw("L");
 
 TGraph gCa40vol(150,xCa40,volCa40);
 gCa40vol.SetLineColor(2);
 gCa40vol.Draw("L");



}
