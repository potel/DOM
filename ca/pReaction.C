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
Sty->SetPadBottomMargin(.13);
Sty->SetPadTopMargin(.04);
Sty->SetPadLeftMargin(.2);
Sty->SetPadRightMargin(.04);
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
Sty->SetTitleOffset(1.8,"y");
Sty->SetTitleOffset(1.1,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.02,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();


 double offset;
 double delta = 300.;

 float All;
 int lineWidth = 2;


 TCanvas canvas("pReactionCa","",800,800);
   TH2S frame ("frame","",10,0,199,10,0,2100);
   //frame->GetXaxis()->SetTitle("E_{lab} [MeV]");
 frame->GetYaxis()->SetTitle("#sigma_{react} [mb]");
 frame->GetXaxis()->SetTitle("E_{lab} [MeV]");
 frame->GetYaxis()->CenterTitle();
 frame->GetXaxis()->CenterTitle();
 frame->Draw();




 TLatex text;
 text.SetTextAngle(-3);
 //text.SetNDC();
 text.SetTextSize(.05);
 text.SetTextColor(2);
 text.DrawLatex(150,610,"p+^{40}Ca");

 ifstream file("/home/charity/DOM48/ca/pca40r.dat");
 int Npoints;
 float x[300];
 float y[300];
 float Volume[300];
 float isoVector[300];
 float isoScaler[300];
 float error[300];
 float xerror[300];

 file >> Npoints;
for (int i=0;i<Npoints;i++)
{
  file  >> x[i] >> y[i] >> error[i];
  xerror[i] = 0.;
} 

TGraphErrors *pca40data = new TGraphErrors(Npoints,x,y,xerror,error);
pca40data->SetMarkerStyle(21);
pca40data->SetMarkerColor(2);
pca40data->SetLineColor(2);
pca40data->Draw("P");


//fit
file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> All >> isoScaler[i] >> isoVector[i] >> Volume[i];

} 




TGraph *pca40fit = new TGraph(Npoints,x,y);
pca40fit->SetLineWidth(lineWidth);
 pca40fit->SetLineColor(2);
pca40fit->Draw("C");

file.close();
file.clear();

 offset += delta;
 text.SetTextColor(4);
 text.DrawLatex(150,950,"p+^{42}Ca");
file.open("/home/charity/DOM48/ca/pca42r.dat");


file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> error[i];
  xerror[i] = 0.;
  y[i] += offset;
} 
TGraphErrors *pca42data = new TGraphErrors(Npoints,x,y,xerror,error);
pca42data->SetMarkerStyle(20);
pca42data->SetMarkerColor(4);
pca42data->SetLineColor(4);
pca42data->Draw("P");
//fit
file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> All >> isoScaler[i] >> isoVector[i] >> Volume[i];
  y[i]+= offset;
} 

TGraph *pca42fit = new TGraph(Npoints,x,y);
pca42fit->SetLineWidth(lineWidth);
 pca42fit->SetLineColor(4);
pca42fit->Draw("C");


file.close();
file.clear();


 offset+= delta;
 text->SetTextColor(2);
 text.DrawLatex(150,1280,"p+^{44}Ca");
file.open("/home/charity/DOM48/ca/pca44r.dat");


file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> error[i];
  xerror[i] = 0.;
y[i] += offset;
} 
TGraphErrors *pca44data = new TGraphErrors(Npoints,x,y,xerror,error);
pca44data->SetMarkerStyle(21);
pca44data->SetMarkerColor(2);
pca44data->SetLineColor(2);
pca44data->Draw("P");

//fit
file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> All >> isoScaler[i] >> isoVector[i] >> Volume[i];
  y[i] += offset;
} 

TGraph *pca44fit = new TGraph(Npoints,x,y);
pca44fit->SetLineWidth(lineWidth);
 pca44fit->SetLineColor(2);
pca44fit->Draw("C");



file.close();
file.clear();

 offset+= delta;
 text->SetTextColor(4);
 text.DrawLatex(150,1630,"p+^{48}Ca");
file.open("/home/charity/DOM48/ca/pca48r.dat");


file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> error[i];
  xerror[i] = 0.;
y[i] += offset;
} 
TGraphErrors *pca48data = new TGraphErrors(Npoints,x,y,xerror,error);
pca48data->SetMarkerStyle(21);
pca48data->SetMarkerColor(4);
pca48data->SetLineColor(4);
pca48data->Draw("P");

//fit
file >> Npoints;
cout << "Npoints= " << Npoints << endl;
for (int i=0;i<Npoints;i++)
{
  file >> x[i] >> y[i] >> All >> isoScaler[i] >> isoVector[i] >> Volume[i];
  y[i] += offset;
} 

TGraph *pca48fit = new TGraph(Npoints,x,y);
pca48fit->SetLineWidth(lineWidth);
 pca48fit->SetLineColor(4);
pca48fit->Draw("C");



file.close();
file.clear();



}

