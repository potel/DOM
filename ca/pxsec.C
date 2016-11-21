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
Sty->SetPadLeftMargin(.27);
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
Sty->SetTitleOffset(2.,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.07,"xz");
Sty->SetTickLength(.05,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("pxsecCa42");


//divide up canvas into three pads
 double overlap = .1;
 double dist = (1.+3*overlap)/4.;
 TPad *pad1  = new TPad("pad1","",0.,0.,.5+overlap,1.);
 TPad *pad2  = new TPad("pad2","",.5-overlap,0.,1.,1.);
pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

gPad->GetFrame()->SetLineWidth(1);
pad1->Draw();
pad2->Draw();



double All,isoScaler,isoVector,Volume,standard;
float lineWidth = 2.;
double ymax = 1300.;

 pad1->cd();
   TH2S frame ("frame","",10,0,99,10,0,ymax);
 frame->GetXaxis()->SetTitle("E_{lab} [MeV]");
 frame->GetYaxis()->SetTitle("#sigma_{react} [mb]");
 frame->GetYaxis()->CenterTitle();
 frame->Draw();

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.08);
 text.DrawLatex(.5,.35,"p+^{42}Ca");

 ifstream file("pca42r.dat");
 int Npoints;
 float x[300];
 float y[300];
 float error[300];
 float xerror[300];

 file >> Npoints;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> error[i];
  xerror[i] = 0.;
} 

TGraphErrors *pni58data = new TGraphErrors(Npoints,x,y,xerror,error);
pni58data->SetMarkerStyle(21);
pni58data->Draw("P");


//fit
file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> All >> isoScaler >> isoVector >> Volume;
} 

TGraph *pni58fit = new TGraph(Npoints,x,y);
pni58fit->SetLineWidth(lineWidth);
 pni58fit->SetLineColor(2);
pni58fit->Draw("C");
file.close();
file.clear();
//***************************************************************

 pad2->cd();
   TH2S frame2 ("frame2","",10,0,99,10,0,ymax);
 frame2->GetXaxis()->SetTitle("E_{lab} [MeV]");
 frame2->GetYaxis()->SetLabelSize(0.);
 frame2->Draw();
 text.DrawLatex(.5,.35,"p+^{44}Ca");
file.open("pca44r.dat");


file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> error[i];
  xerror[i] = 0.;
} 
TGraphErrors *pni60data = new TGraphErrors(Npoints,x,y,xerror,error);
pni60data->SetMarkerStyle(21);
pni60data->Draw("P");
//fit
file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> All >> isoScaler >> isoVector >> Volume;
} 

TGraph *pni60fit = new TGraph(Npoints,x,y);
pni60fit->SetLineWidth(lineWidth);
 pni60fit->SetLineColor(2);
pni60fit->Draw("C");
file.close();
file.clear();



}
