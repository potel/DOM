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

 TCanvas canvas("ComImagCa3");


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

  ifstream filepPb("pca48.poten");
  ifstream filenPb("nca48.poten");
  ifstream fileCa40("../ca/pca40.poten");
  ifstream filenCa40("../ca/nca40.poten");
  ifstream filepCa42("../ca/pca42.poten");
  ifstream filepCa44("../ca/pca44.poten");

  double xpPb[150];
  double surpPbscaler[150];
  double surpPbvector[150];
  double surpPb[150];
  double volpPb[150];
  double HFpPb[150];
  double x[150];


  double xnPb[150];
  double surnPbscaler[150];
  double surnPbvector[150];
  double surnPb[150];
  double volnPb[150];
  double HFnPb[150];


  double xCa40[150];
  double surCa40scaler[150];
  double surCa40vector[150];
  double surCa40[150];
  double volCa40[150];
  double HFCa40[150];

  double xnCa40[150];
  double surnCa40scaler[150];
  double surnCa40vector[150];
  double surnCa40[150];
  double volnCa40[150];
  double HFnCa40[150];


  double xpCa42[150];
  double surpCa42scaler[150];
  double surpCa42vector[150];
  double surpCa42[150];
  double volpCa42[150];
  double HFpCa42[150];

  double xpCa44[150];
  double surpCa44scaler[150];
  double surpCa44vector[150];
  double surpCa44[150];
  double volpCa44[150];
  double HFpCa44[150];


  double one,two,three,four, five,six;

for (int i=0;i<150;i++) 
{
  filepPb >> xpPb[i] >> surpPbscaler[i] >> surpPbvector[i] >> 
  volpPb[i] >> HFpPb[i] >> 
  one >> two >> three >> four >> five >> six;
 surpPb[i] = surpPbscaler[i] + surpPbvector[i];
 //cout << x[i] << " " << surpPb[i] << endl;


  filenPb >> xnPb[i] >> surnPbscaler[i] >> surnPbvector[i] >> 
  volnPb[i] >> HFnPb[i] >> 
  one >> two >> three >> four >> five >> six;
 surnPb[i] = surnPbscaler[i] + surnPbvector[i];


  fileCa40 >> xCa40[i] >> surCa40scaler[i] >> surCa40vector[i] >> 
  volCa40[i] >> HFCa40[i] >> 
  one >> two >> three >> four >> five >> six;
 surCa40[i] = surCa40scaler[i] + surCa40vector[i];


  filenCa40 >> xnCa40[i] >> surnCa40scaler[i] >> surnCa40vector[i] >> 
  volnCa40[i] >> HFnCa40[i] >> 
  one >> two >> three >> four >> five >> six;
 surnCa40[i] = surnCa40scaler[i] + surnCa40vector[i];

  filepCa42 >> xpCa42[i] >> surpCa42scaler[i] >> surpCa42vector[i] >> 
  volpCa42[i] >> HFpCa42[i] >> 
  one >> two >> three >> four >> five >> six;
 surpCa42[i] = surpCa42scaler[i] + surpCa42vector[i];


  filepCa44 >> xpCa44[i] >> surpCa44scaler[i] >> surpCa44vector[i] >> 
  volpCa44[i] >> HFpCa44[i] >> 
  one >> two >> three >> four >> five >> six;
 surpCa44[i] = surpCa44scaler[i] + surpCa44vector[i];

}

 TGraph gpPb(150,xpPb,surpPb);
 gpPb.SetLineColor(4);
 gpPb.Draw("L");


 TGraph gnPb(150,xnPb,surnPb);
 gnPb.SetLineColor(3);
 gnPb.Draw("L");
 
 TGraph gCa40(150,xCa40,surCa40);
 gCa40.SetLineColor(2);
 gCa40.Draw("L");

 TGraph gnCa40(150,xnCa40,surnCa40);
 gnCa40.SetLineColor(6);
 gnCa40.Draw("L");

 TGraph gpCa42(150,xpCa42,surpCa42);
 gpCa42.SetLineColor(7);
 gpCa42.Draw("L");


 TGraph gpCa44(150,xpCa44,surpCa44);
 gpCa44.SetLineColor(8);
 gpCa44.Draw("L");

 TLatex text;
 text.SetTextSize(.08);
 text.SetNDC();
 text.SetTextColor(4);
 text.DrawLatex(.25,.7,"p+^{48}Ca");
 text.SetTextColor(3);
 text.DrawLatex(.75,.7,"n+^{48}Ca");
 text.SetTextColor(2);
 text.DrawLatex(.75,.5,"p+^{40}Ca");
 text.SetTextColor(6);
 text.DrawLatex(.25,.5,"n+^{40}Ca");
 text.SetTextColor(7);
 text.DrawLatex(.75,.3,"p+^{42}Ca");
 text.SetTextColor(8);
 text.DrawLatex(.25,.3,"p+^{44}Ca");
 pad2->cd();
  TH2S framev("framev","",10,-200,200,10,0,12);
  //framev.GetXaxis()->SetTitle("E [MeV]");
  framev.GetXaxis()->SetLabelSize(0.);
  framev.GetYaxis()->SetTitle("W_{vol} [MeV]");
  framev.SetStats(kFALSE);
  framev.Draw();


 TGraph gpPbvol(150,xpPb,volpPb);
 gpPbvol.SetLineColor(4);
 gpPbvol.Draw("L");

 TGraph gnPbvol(150,xnPb,volnPb);
 gnPbvol.SetLineColor(3);
 gnPbvol.Draw("L");
 
 TGraph gCa40vol(150,xCa40,volCa40);
 gCa40vol.SetLineColor(2);
 gCa40vol.Draw("L");


 TGraph gnCa40vol(150,xnCa40,volnCa40);
 gnCa40vol.SetLineColor(6);
 gnCa40vol.Draw("L");



}
