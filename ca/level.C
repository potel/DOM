


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

 TCanvas canvas("level","",400,500);

//divide up canvas into three pads
 double overlap = .05;
 double dist = (1.+3*overlap)/4.;
 TPad *pad2  = new TPad("pad2","",0.,0.,1.,.5+overlap/2.);
 TPad *pad3  = new TPad("pad3","",0.,.5-overlap,1.,1.);
pad2->SetFillStyle(4000);
pad3->SetFillStyle(4000);


pad2->Draw();
pad3->Draw();





  double ymax = 5.;
  double ymin = -25.;
  TH2S frame("frame","",10,0,1,10,ymin,ymax);
  frame.GetXaxis()->SetLabelSize(0.);
  frame.GetXaxis()->SetTickLength(0.);
  frame.GetXaxis()->SetAxisColor(0.);
  frame.GetYaxis()->SetTitle("E_{njl} [MeV]");
  frame.GetYaxis()->CenterTitle();




  string out;
  ifstream fileE;
  ifstream fileT;


  TLine line;
line.SetLineColor(1);
line.SetLineStyle(2);

TLatex text;
text.SetNDC();

//************************
 pad2->cd();

  TH2S frame2("frame2","",10,0,1,10,ymin,ymax);
  frame2.GetXaxis()->SetLabelSize(0.);
  frame2.GetXaxis()->SetTickLength(0.);
  frame2.GetXaxis()->SetAxisColor(0.);
  frame2.GetYaxis()->SetTitle("E_{njl} [MeV]");
  frame2.GetYaxis()->CenterTitle();
  frame2.Draw();



  fileE.close();
  fileE.clear();
  fileT.close();
  fileT.clear();

  cout << "pca40" << endl;
  fileE.open("../ca/pca40.lev");

  getline(fileE,out);
  getline(fileE,out);

  double E, error1,j,rms,error2,Delta,error3,Spect,error4,Occ;
  int n,l,icolor,Efit,Rfit,Dfit,Sfit;

  double energyE[20];
  int lE[20];
  int nE[20];
  double jE[20];  


  int NE = 0;
  for (;;)
    {
      fileE >> E >> error1 >> n >> j >> l >> icolor >> Efit >> rms >> error2 >> Rfit >> Delta >> error3 >> Dfit >> Spect >> error4 >> Sfit;
      if (fileE.eof())break;
      if (fileE.bad())break;

      energyE[NE] = E;
      lE[NE] = l;
      nE[NE] = n;
      jE[NE] = j;
      NE++;
    }

  fileT.open("../ca/pca40.level");
  getline(fileT,out);

  double energyT[40];
  int lT[40];
  int nT[40];
  double jT[40];  

  double startE = .05;
  double stopE = .22;
  double startT = .28;
  double stopT = .45;
  TLine line;

  int NT = 0;
  for (;;)
    {
      fileT >> n >> j >> l >> E >> rms >> Occ >> Spect >> Delta >> 
      error1 >> icolor;
      if (fileT.eof())break;
      if (fileT.bad())break;


      energyT[NT] = E;
      nT[NT] = n;
      lT[NT] = l;
      jT[NT] = j;
      NT++;
    }
  line.SetLineStyle(1);
  for (int i=0;i<NE;i++)
    {
      if (energyE[i] < ymin || energyE[i] > ymax) continue;
      line.SetLineColor(lE[i]+1);
      line.DrawLine(startE,energyE[i],stopE,energyE[i]);
      for (int k=0;k<NT;k++)
	{
          if (nE[i] == nT[k] && jE[i] == jT[k] && lE[i] == lT[k])
	    {
	    line.DrawLine(stopE,energyE[i],startT,energyT[k]);
	    line.DrawLine(startT,energyT[k],stopT,energyT[k]);
            }
	}
    }

  fileE.close();
  fileE.clear();
  fileT.close();
  fileT.clear();

  cout << "nca40" << endl;
  fileE.open("../ca/nca40.lev");

  getline(fileE,out);
  getline(fileE,out);

  double energyEn[20];
  int lEn[20];
  int nEn[20];
  double jEn[20];  

  int NEn = 0;
  for (;;)
    {
      fileE >> E >> error1 >> n >> j >> l >> icolor >> Efit >> rms >> error2 >> Rfit >> Delta >> error3 >> Dfit >> Spect >> error4 >> Sfit;
      if (fileE.eof())break;
      if (fileE.bad())break;

      energyEn[NEn] = E;
      lEn[NEn] = l;
      nEn[NEn] = n;
      jEn[NEn] = j;
      cout << lEn[NEn] << " " << nEn[NEn] << " " << jEn[NEn] << endl;
      NEn++;
    }


  cout << "theory" << endl;

  fileT.open("../ca/nca40.level");
  getline(fileT,out);

  double energyTn[40];
  int lTn[40];
  int nTn[40];
  double jTn[40];  


  int NTn = 0;
  for (;;)
    {
      fileT >> n >> j >> l >> E >> rms >> Occ >> Spect >> Delta >> 
      error1 >> icolor;
      if (fileT.eof())break;
      if (fileT.bad())break;


      energyTn[NTn] = E;
      nTn[NTn] = n;
      lTn[NTn] = l;
      jTn[NTn] = j;
      NTn++;
    }


  double startEn = .05 + .5;
  double stopEn = .22 + .5;
  double startTn = .28 + .5;
  double stopTn = .45 + .5;

  for (int i=0;i<NEn;i++)
    {
      if (energyEn[i] < ymin || energyEn[i] > ymax) continue;
      line.SetLineColor(lEn[i]+1);
      line.DrawLine(startEn,energyEn[i],stopEn,energyEn[i]);
      for (int k=0;k<NTn;k++)
	{
          if (nEn[i] == nTn[k] && jEn[i] == jTn[k] && lEn[i] == lTn[k])
	    {
	    line.DrawLine(stopEn,energyEn[i],startTn,energyTn[k]);
	    line.DrawLine(startTn,energyTn[k],stopTn,energyTn[k]);
            }
	}
    }





  ifstream file;
  double efermin,efermip;
ifstream file("../ca/pca40.inp");
file >> error1 >> error2 >> error3 >> efermip;
file.close();
file.clear();
file.open("../ca/nca40.inp");
file >> error1 >> error2 >> error3 >> efermin;

line.SetLineColor(1);
line.SetLineStyle(2);
line.DrawLine(startE,efermip,stopT,efermip);
line.DrawLine(startEn,efermin,stopTn,efermin);

text.DrawLatex(.32,.95,"p+^{40}Ca");
text.DrawLatex(.67,.2,"n+^{40}Ca"); 


 pad3->cd();

  TH2S frame3("frame3","",10,0,1,10,ymin,ymax);
  frame3.GetXaxis()->SetLabelSize(0.);
  frame3.GetXaxis()->SetTickLength(0.);
  frame3.GetXaxis()->SetAxisColor(0.);
  frame3.GetYaxis()->SetTitle("E_{njl} [MeV]");
  frame3.GetYaxis()->CenterTitle();
  frame3.Draw();



  fileE.close();
  fileE.clear();
  fileT.close();
  fileT.clear();

  cout << "pca48" << endl;

  fileE.open("../ca/pca48.lev");

  getline(fileE,out);
  getline(fileE,out);

  NE = 0;
  for (;;)
    {
      fileE >> E >> error1 >> n >> j >> l >> icolor >> Efit >> rms >> error2 >> Rfit >> Delta >> error3 >> Dfit >> Spect >> error4 >> Sfit;
      if (fileE.eof())break;
      if (fileE.bad())break;

      energyE[NE] = E;
      lE[NE] = l;
      nE[NE] = n;
      jE[NE] = j;
      NE++;
    }

  fileT.open("../ca/pca48.level");
  getline(fileT,out);


  NT = 0;
  for (;;)
    {
      fileT >> n >> j >> l >> E >> rms >> Occ >> Spect >> Delta >> 
      error1 >> icolor;
      if (fileT.eof())break;
      if (fileT.bad())break;


      energyT[NT] = E;
      nT[NT] = n;
      lT[NT] = l;
      jT[NT] = j;
      NT++;
    }
  line.SetLineStyle(1);
  for (int i=0;i<NE;i++)
    {
      line.SetLineColor(lE[i]+1);
      line.DrawLine(startE,energyE[i],stopE,energyE[i]);
      for (int k=0;k<NT;k++)
	{
          if (energyE[i] < ymin || energyE[i] > ymax) continue;
          if (nE[i] == nT[k] && jE[i] == jT[k] && lE[i] == lT[k])
	    {
	    line.DrawLine(stopE,energyE[i],startT,energyT[k]);
	    line.DrawLine(startT,energyT[k],stopT,energyT[k]);
            }
	}
    }

  fileE.close();
  fileE.clear();
  fileT.close();
  fileT.clear();

  cout << "Nca48" << endl;
  fileE.open("../ca/nca48.lev");

  getline(fileE,out);
  getline(fileE,out);


  NEn = 0;
  for (;;)
    {
      fileE >> E >> error1 >> n >> j >> l >> icolor >> Efit >> rms >> error2 >> Rfit >> Delta >> error3 >> Dfit >> Spect >> error4 >> Sfit;
      if (fileE.eof())break;
      if (fileE.bad())break;

      energyEn[NEn] = E;
      lEn[NEn] = l;
      nEn[NEn] = n;
      jEn[NEn] = j;
      NEn++;
    }

  fileT.open("../ca/nca48.level");
  getline(fileT,out);


  NTn = 0;
  for (;;)
    {
      fileT >> n >> j >> l >> E >> rms >> Occ >> Spect >> Delta >> 
      error1 >> icolor;
      if (fileT.eof())break;
      if (fileT.bad())break;


      energyTn[NTn] = E;
      nTn[NTn] = n;
      lTn[NTn] = l;
      jTn[NTn] = j;
      NTn++;
    }

  for (int i=0;i<NEn;i++)
    {
      line.SetLineColor(lEn[i]+1);
      line.DrawLine(startEn,energyEn[i],stopEn,energyEn[i]);
      for (int k=0;k<NTn;k++)
	{
          if (energyEn[i] < ymin || energyEn[i] > ymax) continue;
          if (nEn[i] == nTn[k] && jEn[i] == jTn[k] && lEn[i] == lTn[k])
	    {
	    line.DrawLine(stopEn,energyEn[i],startTn,energyTn[k]);
	    line.DrawLine(startTn,energyTn[k],stopTn,energyTn[k]);
            }
	}
    }





  file.close();
  file.clear();
ifstream file("../ca/pca48.inp");
file >> error1 >> error2 >> error3 >> efermip;
file.close();
file.clear();
file.open("../ca/nca48.inp");
file >> error1 >> error2 >> error3 >> efermin;

line.SetLineColor(1);
line.SetLineStyle(2);
line.DrawLine(startE,efermip,stopT,efermip);
line.DrawLine(startEn,efermin,stopTn,efermin);

text.DrawLatex(.32,.8,"p+^{48}Ca");
text.DrawLatex(.67,.3,"n+^{48}Ca"); 

}
