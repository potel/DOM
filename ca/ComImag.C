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

 TCanvas canvas("ComImagCa");


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
  ifstream filepCa40("../ca/pca40.poten");
  ifstream filenCa40("../ca/nca40.poten");

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


  double xpCa40[150];
  double surpCa40scaler[150];
  double surpCa40vector[150];
  double surpCa40[150];
  double volpCa40[150];
  double HFpCa40[150];

  double xnCa40[150];
  double surnCa40scaler[150];
  double surnCa40vector[150];
  double surnCa40[150];
  double volnCa40[150];
  double HFnCa40[150];


  double zero,one,two,three,four, five,six,seven;

for (int i=0;i<150;i++) 
{
  filepPb >> xpPb[i] >> surpPbscaler[i] >> surpPbvector[i] >> 
    volpPb[i] >> HFpPb[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surpPb[i] = surpPbscaler[i] + surpPbvector[i];
 //cout << x[i] << " " << surpPb[i] << endl;


  filenPb >> xnPb[i] >> surnPbscaler[i] >> surnPbvector[i] >> 
    volnPb[i] >> HFnPb[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surnPb[i] = surnPbscaler[i] + surnPbvector[i];


  filepCa40 >> xpCa40[i] >> surpCa40scaler[i] >> surpCa40vector[i] >> 
    volpCa40[i] >> HFpCa40[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surpCa40[i] = surpCa40scaler[i] + surpCa40vector[i];
 // cout << x[i] << " " << surpCa40scaler[i] << endl;

  filenCa40 >> xnCa40[i] >> surnCa40scaler[i] >> surnCa40vector[i] >> 
    volnCa40[i] >> HFnCa40[i] >> zero >>
    one >> two >> three >> four >> five >> six >> seven;
 surnCa40[i] = surnCa40scaler[i] + surnCa40vector[i];
 // cout << x[i] << " " << surnCa40scaler[i] << endl;
}

 TGraph gpPb(150,xpPb,surpPbscaler);
 gpPb.SetLineColor(4);
 gpPb.Draw("L");


 TGraph gnPb(150,xnPb,surnPbscaler);
 gnPb.SetLineColor(3);
 gnPb.Draw("L");
 
 TGraph gpCa40(150,xpCa40,surpCa40scaler);
 gpCa40.SetLineColor(2);
 gpCa40.Draw("L");

 TGraph gnCa40(150,xnCa40,surnCa40scaler);
 gnCa40.SetLineColor(1);
 gnCa40.Draw("L");

 TLatex text;
 text.SetTextSize(.08);
 text.SetNDC();
 text.SetTextColor(4);
 text.DrawLatex(.25,.7,"p+^{48}Ca");
 text.SetTextColor(3);
 text.DrawLatex(.75,.7,"n+^{48}Ca");
 text.SetTextColor(2);
 text.DrawLatex(.25,.5,"p+^{40}Ca");
 text.SetTextColor(1);
 text.DrawLatex(.75,.5,"n+^{40}Ca");

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
 
 TGraph gpCa40vol(150,xpCa40,volpCa40);
 gpCa40vol.SetLineColor(2);
 gpCa40vol.Draw("L");


 TGraph gnCa40vol(150,xnCa40,volnCa40);
 gnCa40vol.SetLineColor(1);
 gnCa40vol.Draw("L");


 TGraph gpPbv(150,xpPb,surpPbvector);
 gpPbv.SetLineColor(4);
 gpPbv.Draw("L");


 TGraph gnPbv(150,xnPb,surnPbvector);
 gnPbv.SetLineColor(3);
 gnPbv.Draw("L");
 
 TGraph gpCa40v(150,xpCa40,surpCa40vector);
 gpCa40v.SetLineColor(2);
 gpCa40v.Draw("L");

 TGraph gnCa40v(150,xnCa40,surnCa40vector);
 gnCa40v.SetLineColor(1);
 gnCa40v.Draw("L");


}
