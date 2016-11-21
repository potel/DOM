{
gROOT->Reset();
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
Sty->SetLabelSize(0.05,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
Sty->SetTitleOffset(1.,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.025,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

TCanvas can("vca");

double A[4]= {40,42,44,48};
for (int i=0;i<4;i++)A[i] = (A[i] - 40.)/A[i];
double JW21[4] = {104.96,124.63,132.96,134.08};
double JW25[4] = {88.86,108.0,122.4,122.28};
double JW30[4] = {96.86,112.86,108.98,112.43};
double JW35[4] = {106.87,106.30,109.13,112.55};
double JW40[4] = {98.12,82.07,104.82,110.30};
double JW45[4] = {96.90,104.07,100.71,102.05};
double JW48[4] = {105.26,105.79,105.45,109.10};

for (int i=0;i<4;i++)
{
  JW21[i] += 60;
  JW25[i] += 50;
  JW30[i] += 40.;
    JW35[i] += 30.;
    JW40[i] += 20.;
    JW45[i] += 25.;
}

TF1 *lin21 = new TF1("lin21","pol1");
lin21->SetLineColor(1);
TF1 *lin25 = new TF1("lin25","pol1");
lin25->SetLineColor(2);
TF1 *lin40 = new TF1("lin40","pol1");
lin40->SetLineColor(6);
TF1 *lin48 = new TF1("lin48","pol1");
lin48->SetLineColor(4);


TGraph *graph21 = new TGraph(4,A,JW21);
graph21->SetMarkerStyle(21);
graph21->SetMarkerSize(2);
TGraph *graph25 = new TGraph(4,A,JW25);
graph25->SetLineColor(2);
graph25->SetMarkerStyle(22);
graph25->SetMarkerSize(2);
graph25->SetMarkerColor(2);
TGraph *graph30 = new TGraph(4,A,JW30);
graph30->SetLineColor(6);
TGraph *graph35 = new TGraph(4,A,JW35);
graph35->SetMarkerStyle(23);
graph35->SetMarkerColor(3);
graph35->SetMarkerSize(2);
graph35->SetLineColor(3);
TGraph *graph40 = new TGraph(4,A,JW40);
graph40->SetLineColor(5);
TGraph *graph45 = new TGraph(4,A,JW45);
graph45->SetLineColor(6);
graph45->SetMarkerColor(6);
graph45->SetMarkerSize(2);
TGraph *graph48 = new TGraph(4,A,JW48);
graph48->SetLineColor(7);
graph48->SetMarkerStyle(20);
graph48->SetMarkerColor(4);
graph48->SetMarkerSize(2);

TH2S frame ("frame","",10,-.01,.2,10,80,200);
frame.SetStats(kFALSE);
frame.Draw();
frame.GetYaxis()->SetTitle("J_{W} [MeV fm^{3}]");
frame.GetYaxis()->CenterTitle();
frame.GetXaxis()->SetTitle("(N-Z)/A");
frame.GetXaxis()->CenterTitle();

graph21->Draw("PL");
//graph21->Fit(lin21);
graph25->Draw("PL");
//graph25->Fit(lin25);
//graph30->Draw("PL");
//graph30->Fit("pol1");
graph35->Draw("PL");
//graph35->Fit("pol1");
//graph40->Draw("*");
//graph40->Fit("pol1");
graph45->Draw("PL");
//graph45->Fit(lin40);
graph48->Draw("PL");
//graph48->Fit(lin48);

TLatex text;
text.SetNDC();
text.SetTextColor(4);
text.DrawLatex(.65,.33,"E=48.4 MeV");
text.SetTextColor(6);
text.DrawLatex(.69,.45,"45 MeV (+25)");
text.SetTextColor(3);
text.DrawLatex(.69,.54,"35 MeV (+30)");
text.SetTextColor(2);
text.DrawLatex(.69,.67,"25 MeV (+50)");
text.SetTextColor(1);
text.DrawLatex(.69,.80,"21 MeV (+60)");

}
