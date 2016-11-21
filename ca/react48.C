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

 TCanvas canvas("react48Ca");


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



double All,isoScaler,isoVector,Volume,standard;
float lineWidth = 2.;
double ymax = 1200.;

 pad1->cd();
 TH2S frame ("frame","",10,0,200,10,-100,ymax);
 frame->GetXaxis()->SetTitle("E_{lab} [MeV]");
 frame->GetXaxis()->CenterTitle();
 frame->GetYaxis()->SetTitle("#sigma_{react} [mb]");
 frame->GetYaxis()->CenterTitle();
 frame->Draw();

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.1);
 text.DrawLatex(.5,.35,"(b) p+^{48}Ca");

 ifstream file("pca48r.dat");
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

 cout << " 1 " << Npoints << endl;

TGraphErrors *p208pbdata = new TGraphErrors(Npoints,x,y,xerror,error);
p208pbdata->SetMarkerStyle(21);
p208pbdata->Draw("P");


//fit
file >> Npoints;
cout << "Npoints= " << Npoints << endl;
 float yS[300];
 float yA[300];
 float yVV[300];
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> All >> yS[i] >> isoVector >> Volume;
  yVV[i] = Volume;
  yA[i] = isoVector;
} 

TGraph *p208pbfit = new TGraph(Npoints,x,y);
p208pbfit->SetLineWidth(lineWidth);
p208pbfit->SetLineColor(2);

p208pbfit->Draw("L");


TGraph *p208pbfitS = new TGraph(Npoints,x,yS);
p208pbfitS->SetLineWidth(lineWidth);
p208pbfitS->SetLineColor(2);
p208pbfitS->SetLineStyle(8);
p208pbfitS->Draw("L");

TGraph *p208pbfitV = new TGraph(Npoints,x,yVV);
p208pbfitV->SetLineWidth(lineWidth);
p208pbfitV->SetLineColor(2);
p208pbfitV->SetLineStyle(2);
p208pbfitV->Draw("L");


TGraph *p208pbfitA = new TGraph(Npoints,x,yA);
p208pbfitA->SetLineWidth(lineWidth);
p208pbfitA->SetLineColor(4);
p208pbfitA->SetLineStyle(2);
p208pbfitA->Draw("L");




file.close();
file.clear();
double xOM[30];
double yOM[30];
double x_OM,y_OM;
int n_OM=0;
int ii;
file.open("pca48so.x");
string stuff;
for (;;)
  {
    file >> ii >> x_OM >> y_OM;
    if (file.eof()) break;
    if (file.bad()) break;
    if (x_OM < 0.) break;
    xOM[n_OM] = x_OM;
    yOM[n_OM] = y_OM;
    n_OM++;
    getline(file,stuff);
  }

TGraph gOM(n_OM,xOM,yOM);
gOM.SetMarkerStyle(20);
gOM.SetMarkerColor(3);
gOM.SetMarkerSize(.7);
gOM.Draw("P");

file.close();
file.clear();



//***************************************************************

 pad2->cd();
   TH2S frame2 ("frame2","",10,0,200,10,0,4000);
   //frame2->GetXaxis()->SetTitle("E_{lab} [MeV]");
   frame2->GetXaxis()->SetLabelSize(0.);
   frame2->GetYaxis()->SetTitle("#sigma [mb]");
   frame2->GetYaxis()->CenterTitle();
 frame2->Draw();
 text.DrawLatex(.5,.75,"(a) n+^{48}Ca");
text.SetTextColor(4);
text.DrawLatex(.8,.6,"#sigma_{tot}");
text.SetTextColor(2);
text.DrawLatex(.3,.38,"#sigma_{react}");
 file.open("nca48r.dat");
 int NNpoints;

 float xn[1];
 float yn[1];
 float errorn[1];
 float xerrorn[1];


 float xnf[300];
 float ynf[300];

 float xt[400];
 float yt[400];
 float errort[400];
 float xerrort[400];

 float xtf[400];
 float ytf[400];

/*
 file >> NNpoints;
for (int i=0;i<NNpoints;i++)
{
  file >> xn[i] >> yn[i] >> errorn[i];
  xerrorn[i] = 0.;
} 

 TGraphErrors n208Pba(NNpoints,xn,yn,xerrorn,errorn);
 n208Pba.SetMarkerStyle(20);
 n208Pba.SetMarkerColor(1);
 n208Pba.Draw("P");
*/
//fit
 int NfNpoints;
 file >>  NfNpoints;

 float ynS[300];
 float ynV[300];
for (int i=0;i<NfNpoints;i++)
{
  file >> xnf[i] >> ynf[i] >> All >> isoScaler >> isoVector >> Volume;
  ynS[i] = isoScaler;
  ynV[i] = Volume;
} 

 TGraph n208Pbaf(NfNpoints,xnf,ynf);
 n208Pbaf.SetLineColor(2);
n208Pbaf.SetLineWidth(lineWidth);
 n208Pbaf.Draw("C");


 TGraph n208PbafS(NfNpoints,xnf,ynS);
 n208PbafS.SetLineColor(2);
 n208PbafS.SetLineStyle(8);
n208PbafS.SetLineWidth(lineWidth);
 n208PbafS.Draw("C");


 TGraph n208PbafV(NfNpoints,xnf,ynV);
 n208PbafV.SetLineColor(2);
 n208PbafV.SetLineStyle(2);
n208PbafV.SetLineWidth(lineWidth);
 n208PbafV.Draw("C");



 int NTpoints;

file >> NTpoints;
for (int i=0;i<NTpoints;i++)
{
  file >> xt[i] >> yt[i] >> errort[i];
  xerrort[i] = 0.;
} 
TGraphErrors *n208pbtdata = new TGraphErrors(NTpoints,xt,yt,xerrort,errort);
n208pbtdata->SetMarkerStyle(20);
n208pbtdata->SetMarkerSize(.4);
 n208pbtdata->SetMarkerColor(1);
 n208pbtdata->SetLineColor(1);
n208pbtdata->Draw("P");

 int NTFpoints;
file >> NTFpoints;
for (int i=0;i<NTFpoints;i++)
{
  file >> xtf[i] >> ytf[i];
} 

TGraph *n208pbtfit = new TGraph(NTFpoints,xtf,ytf);
n208pbtfit->SetLineWidth(2);
 n208pbtfit->SetLineColor(4);
n208pbtfit->Draw("C");
file.close();
file.clear();

double xxOM[30];
double yyOM[30];
double ttOM[30];
int nnOM = 0;
double one,two,three,four,five;
file.open("nca48sat.x");
for(;;)
  {
    file >> ii >> one >> two >> three >> four >> five;
    if (file.eof()) break;
    if (file.bad()) break;
    if (one < 0.) break;
    xxOM[nnOM] = one;
    yyOM[nnOM] = two;
    ttOM[nnOM] = five;
    nnOM++;
    getline(file,stuff);
  }

TGraph ggOM(nnOM,xxOM,yyOM);
ggOM.SetMarkerStyle(20);
ggOM.SetMarkerColor(3);
ggOM.SetMarkerSize(.7);
ggOM.Draw("P");


TGraph gtOM(nnOM,xxOM,ttOM);
gtOM.SetMarkerStyle(21);
gtOM.SetMarkerColor(3);
gtOM.SetMarkerSize(.7);
gtOM.Draw("P");
}
