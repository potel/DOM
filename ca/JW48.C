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
Sty->SetPadBottomMargin(.17);
Sty->SetPadLeftMargin(.23);
Sty->SetHistLineWidth(1);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(1);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.08,"xyz");
Sty->SetTitleOffset(1.2,"y");
Sty->SetTitleOffset(1.,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.05,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

  TCanvas canvas("JW48");
  double delta = .02;
  TPad * pad1 = new TPad("pad1","",0.,0.,.5+delta,1);
  TPad * pad2 = new TPad("pad2","",.5-delta,0,1,1);
  pad1->SetFillStyle(4000);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();

  int N,ii;
  double x[40];
  double xx[40];
  double y[40];
  double yy[40];
  ifstream file("nca48sat.Vint");
  string name;
  getline(file,name);
  cout << name << endl;
  file >> N >> ii; 
  cout << N << " " << ii << endl;

  string out;
  double E,jr,jw,vrms,wrms,jso,sorms;
  double avJw = 0.;
  double Efermi = -7.546;
  for (int i=0;i<N;i++)
    {
     
      file >> E >> jr >> jw >> vrms >> wrms >> jso >> sorms;
      getline(file,name);
      cout << E << " " << jw << endl;
      x[i] = E*48./49. - Efermi;
      xx[i] = E*48./49.;
      y[i] = jw;
      yy[i] = jr;
      avJw += jw;
    }
  file.close();
  file.clear();

  avJw /= (double)N;
  cout << "average JW for Neutrons = " << avJw << endl;


  TH2S frame("frame","",10,0,212,10,0,200);
  frame.GetXaxis()->SetTitle("E-E_{F} [MeV]");
  frame.GetYaxis()->SetTitle("J_{W}/A [MeV]");
  frame.SetStats(kFALSE);
  frame.Draw();

  TGraph gnca48(N,x,y);
  gnca48.SetMarkerStyle(20);
  gnca48.SetMarkerColor(3);
  gnca48.SetLineColor(3);

  gnca48.Draw("P");


  int NP;
  double xp[40];
  double xxp[40];
  double yp[40];
  double yyp[40];
  file.open("pca48sa.Vint");

  getline(file,name);
  cout << name << endl;
  file >> NP >> ii; 
  cout << NP << " " << ii << endl;
  avJw = 0.;
  double EfermiP = -5.79;
  double Cshift = .759*1.7*20./1.22/pow(48.,1./3.);
  cout << "Cshift= " << Cshift << endl;
  for (int i=0;i<NP;i++)
    {
     
      file >> E >> jr >> jw >> vrms >> wrms >> jso >> sorms;
      getline(file,name);
      xp[i] = E*48./49. - EfermiP;
      xxp[i] = E*48./49. - Cshift;
      yp[i] = jw;
      yyp[i] = jr;
      avJw += jw;
    }
  file.close();
  file.clear();
  avJw /= (double)NP;
  cout << "average JW for protons= " << avJw << endl;


  TGraph gpca48(NP,xp,yp);
  gpca48.SetMarkerStyle(21);
  gpca48.SetMarkerColor(4);
  gpca48.SetLineColor(4);

  gpca48.Draw("P");

  TLatex text;
  text.SetNDC();
  text.SetTextColor(3);
  text.DrawLatex(.7,.8,"n+^{48}Ca");
  text.SetTextColor(4);
  text.DrawLatex(.7,.7,"p+^{48}Ca");
  
  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(1);
  line.SetLineColor(3);
  line.DrawLine(7.2,0,7.2,200);
  line.SetLineColor(4);
  line.DrawLine(7.92,0,7.92,200);

  double xdom[200];
  double vdom[200];
  double wdom[200];
  
  int NN;
  ifstream domFile("nca48.Vinteg");
  domFile >> NN;
  getline(domFile,name);
  cout << NN << endl;
  for (int jj=0;jj<NN;jj++) 
    {
    getline(domFile,name);
    }
  domFile >> NN;
  for (int jj=0;jj<NN;jj++)
    {
      domFile >> xdom[jj] >> vdom[jj] >> wdom[jj];
      getline(domFile,name);
   }
  domFile.close();
  domFile.clear();

  TGraph ggw(NN,xdom,wdom);
  ggw.SetLineColor(3);
  ggw.SetLineWidth(1.);
  ggw.Draw("L");


  double xdomP[200];
  double vdomP[200];
  double wdomP[200];


  int PP;
  ifstream domFile("pca48.Vinteg");
  domFile >> PP;
  getline(domFile,name);
  cout << PP << endl;
  for (int jj=0;jj<PP;jj++) 
    {
    getline(domFile,name);
    }
  domFile >> PP;
  for (int jj=0;jj<PP;jj++)
    {
      domFile >> xdomP[jj] >> vdomP[jj] >> wdomP[jj];
      getline(domFile,name);
   }
  domFile.close();
  domFile.clear();

  TGraph ggwP(PP,xdomP,wdomP);
  ggwP.SetLineColor(4);
  ggwP.SetLineWidth(1.);
  ggwP.Draw("L");

  pad2->cd();
  TH2S frame2("frame2","",10,-20,200,10,0,550);
  frame2.GetXaxis()->SetTitle("E-#DeltaE_{C} [MeV]");
  frame2.GetYaxis()->SetTitle("J_{V}/A [MeV]");
  frame2.SetStats(kFALSE);
  frame2.Draw();    


  TGraph gpni58R(NP,xxp,yyp);
  gpni58R.SetMarkerStyle(21);
  gpni58R.SetMarkerColor(4);
  gpni58R.SetLineColor(4);

  gpni58R.Draw("P");


  TGraph gnni58R(N,xx,yy);
  gnni58R.SetMarkerStyle(20);
  gnni58R.SetMarkerColor(3);
  gnni58R.SetLineColor(3);

  gnni58R.Draw("P");

  for (int jj=0;jj<PP;jj++) xdomP[jj] += EfermiP - Cshift; 
  TGraph ggvP(PP,xdomP,vdomP);
  ggvP.SetLineColor(4);
  ggvP.SetLineWidth(1.);
  ggvP.Draw("L");

  for (int jj=0;jj<NN;jj++) xdom[jj] += Efermi; 
  TGraph ggv(NN,xdom,vdom);
  ggv.SetLineColor(3);
  ggv.SetLineWidth(1.);
  ggv.Draw("L");
}
