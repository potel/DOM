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
Sty->SetPadRightMargin(.05);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.09,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.1,"xyz");
Sty->SetTitleOffset(1.,"y");
Sty->SetTitleOffset(1.,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("pReactionCa","",400,700);


//divide up canvas into seven pads
 double lap = 0.08;
 TPad *pad4  = new TPad("pad4","",0.,             0.,1., 1./4.+lap);
 TPad *pad3  = new TPad("pad3","",0.,1./4.-1./3.*lap,1., 2./4.+2.*lap/3.);
 TPad *pad2  = new TPad("pad2","",0.,2./4.-2./3.*lap,1., 3./4.+1.*lap/3.);
 TPad *pad1  = new TPad("pad1","",0.,3./4.-lap,1., 1.);



pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);
pad3->SetFillStyle(4000);
pad4->SetFillStyle(4000);


gPad->GetFrame()->SetLineWidth(1);
pad1->Draw();
pad2->Draw();
pad3->Draw();
pad4->Draw();





 ifstream file;
 TLatex text;
 text.SetNDC();
 text.SetTextColor(1);
 text.SetTextSize(.1);
 double All, isoScaler,isoVector,Volume;
 int lineWidth = 2;
double one,two,three,four,five;
 int ii;


 pad1->cd();
   TH2S frame1 ("frame1","",10,0,200,10,0,1200);
   //frame1->GetXaxis()->SetTitle("E_{lab} [MeV]");
   frame1->GetXaxis()->SetLabelSize(0.);
   frame1->GetYaxis()->SetTitle("#sigma_{react} [mb]");
   frame1->GetYaxis()->CenterTitle();
 frame1->Draw();
 text.DrawLatex(.67,.75,"(a) p+^{40}Ca");


 file.open("/home/charity/DOM48/ca/pca40r.dat");
 int Npoints40;
 float x40[300];
 float y40[300];
 float error40[300];
 float xerror40[300];


 file >> Npoints40;
for (int i=0;i<Npoints40;i++)
{
  file >> x40[i] >> y40[i] >> error40[i];
  xerror40[i] = 0.;
} 



TGraphErrors *xx40 = new TGraphErrors(Npoints40,x40,y40,xerror40,error40);
xx40->SetMarkerStyle(21);
xx40->Draw("P");


//fit
file >> Npoints40;


for (int i=0;i<Npoints40;i++)
{
  file >> x40[i] >> y40[i] >> one >> two >> three >> four;
} 

TGraph *f40 = new TGraph(Npoints40,x40,y40);
f40->SetLineWidth(lineWidth);
f40->SetLineColor(2);

f40->Draw("L");
file.close();
file.clear();


//*****************************************************
 pad2->cd();
   TH2S frame2 ("frame2","",10,0,200,10,0,1200);
   //frame2->GetXaxis()->SetTitle("E_{lab} [MeV]");
   frame2->GetXaxis()->SetLabelSize(0.);
   frame2->GetYaxis()->SetTitle("#sigma_{react} [mb]");
   frame2->GetYaxis()->CenterTitle();
 frame2->Draw();
 text.DrawLatex(.67,.75,"(b) p+^{42}Ca");

 file.open("/home/charity/DOM48/ca/pca42r.dat");
 int Npoints42;
 float x42[300];
 float y42[300];
 float error42[300];
 float xerror42[300];


 file >> Npoints42;
for (int i=0;i<Npoints42;i++)
{
  file >> x42[i] >> y42[i] >> error42[i];
  xerror42[i] = 0.;
} 



TGraphErrors *xx42 = new TGraphErrors(Npoints42,x42,y42,xerror42,error42);
xx42->SetMarkerStyle(21);
xx42->Draw("P");


//fit
file >> Npoints42;


for (int i=0;i<Npoints42;i++)
{
  file >> x42[i] >> y42[i] >> one >> two >> three >> four;
} 

TGraph *f42 = new TGraph(Npoints42,x42,y42);
f42->SetLineWidth(lineWidth);
f42->SetLineColor(2);

f42->Draw("L");
file.close();
file.clear();



//*****************************************************
 pad3->cd();
   TH2S frame3 ("frame3","",10,0,200,10,0,1200);
   //frame3->GetXaxis()->SetTitle("E_{lab} [MeV]");
   frame3->GetXaxis()->SetLabelSize(0.);
   frame3->GetYaxis()->SetTitle("#sigma_{react} [mb]");
   frame3->GetYaxis()->CenterTitle();
 frame3->Draw();
 text.DrawLatex(.67,.75,"(c) p+^{44}Ca");

 file.open("/home/charity/DOM48/ca/pca44r.dat");
 int Npoints44;
 float x44[300];
 float y44[300];
 float error44[300];
 float xerror44[300];


 file >> Npoints44;
for (int i=0;i<Npoints44;i++)
{
  file >> x44[i] >> y44[i] >> error44[i];
  xerror44[i] = 0.;
} 



TGraphErrors *xx44 = new TGraphErrors(Npoints44,x44,y44,xerror44,error44);
xx44->SetMarkerStyle(21);
xx44->Draw("P");


//fit
file >> Npoints44;


for (int i=0;i<Npoints44;i++)
{
  file >> x44[i] >> y44[i] >> one >> two >> three >> four;
} 

TGraph *f44 = new TGraph(Npoints44,x44,y44);
f44->SetLineWidth(lineWidth);
f44->SetLineColor(2);

f44->Draw("L");
file.close();
file.clear();

//*****************************************************
 pad4->cd();
   TH2S frame4 ("frame4","",10,0,200,10,0,1200);
   frame4->GetXaxis()->SetTitle("E_{lab} [MeV]");
frame4->GetXaxis()->CenterTitle();
//frame4->GetXaxis()->SetLabelSize(0.);
   frame4->GetYaxis()->SetTitle("#sigma_{react} [mb]");
   frame4->GetYaxis()->CenterTitle();
 frame4->Draw();
 text.DrawLatex(.67,.75,"(d) p+^{48}Ca");
 file.open("/home/charity/DOM48/ca/pca48r.dat");
 int Npoints48;
 float x48[300];
 float y48[300];
 float error48[300];
 float xerror48[300];


 file >> Npoints48;
for (int i=0;i<Npoints48;i++)
{
  file >> x48[i] >> y48[i] >> error48[i];
  xerror48[i] = 0.;
} 



TGraphErrors *xx48 = new TGraphErrors(Npoints48,x48,y48,xerror48,error48);
xx48->SetMarkerStyle(21);
xx48->Draw("P");


//fit
file >> Npoints48;


for (int i=0;i<Npoints48;i++)
{
  file >> x48[i] >> y48[i] >> one >> two >> three >> four;
} 

TGraph *f48 = new TGraph(Npoints48,x48,y48);
f48->SetLineWidth(lineWidth);
f48->SetLineColor(2);

f48->Draw("L");
file.close();
file.clear();




}



