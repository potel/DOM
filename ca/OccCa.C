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
Sty->SetPadBottomMargin(.2);
Sty->SetPadLeftMargin(.2);
Sty->SetPadTopMargin(.05);
Sty->SetPadRightMargin(.05);
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
Sty->SetTitleOffset(1.4,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.025,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas can("OccCa","",400,600);

//divide up canvas into three pads
 double overlap = .06;
 double dist = (1.+3*overlap)/4.;
 TPad *pad1  = new TPad("pad1","",0.,0.,1.,.5+overlap);
 TPad *pad2  = new TPad("pad2","",0.,.5-overlap,1.,1.);
pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

gPad->GetFrame()->SetLineWidth(1);
pad1->Draw();
pad2->Draw();

 double Z,Zp,Ap;

 pad2->cd();

  TH2S frame("frame","",10,-60,20,10,0,1);
  //frame.GetXaxis()->SetTitle("E_{n,l,j}-E_{Fermi} [MeV]");
  frame.GetXaxis()->SetLabelSize(0.);
  frame.GetYaxis()->SetTitle("Occupancy");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.SetStats(kFALSE);
  frame.Draw();

  TLine line;
  line.SetLineStyle(2);
  line.DrawLine(0.,0.,0.,1.);


  ifstream file;
  file.open("pca48.inp");
  double Efermi48;
  file >> Z >> Zp >> Ap >> Efermi48;


  file.close();
  file.clear();
  file.open("pca48.level");
  string out;
  getline(file,out);

  int i1,i2,i3;
  double d1,d2,d3,d4,d5,d6,d7;

  double x48hole[60];
  double y48hole[60];
  double x48part[60];
  double y48part[60];
  int n48hole = 0;
  int n48part = 0;
  for(;;)
    {
      file >> i1 >> d1 >> i2 >> d2 >> d3 >> d4 >> d5 >> d6 >> d7 >> i3;
      if (file.eof())break;
      if (file.bad())break;
      if (d2 < Efermi48)
	{
         x48hole[n48hole] = d2 - Efermi48;
         y48hole[n48hole] = d4;
         n48hole++;
	}
      else
	{
         x48part[n48part] = d2 - Efermi48;
         y48part[n48part] = d4;
         n48part++;
	}

    }

  

  TGraph g48hole(n48hole,x48hole,y48hole);
  g48hole.SetMarkerStyle(20);
  g48hole.SetMarkerColor(4);
  g48hole.SetLineColor(4);
    g48hole.SetMarkerSize(1.2);
  g48hole.Draw("PL");

  TGraph g48part(n48part,x48part,y48part);
  g48part.SetMarkerStyle(20);
  g48part.SetMarkerColor(4);
  g48part.SetLineColor(4);
    g48part.SetMarkerSize(1.2);
  g48part.Draw("PL");


  file.close();
  file.clear();
  file.open("pca40.inp");
  double Efermi40;
  file >> Z >> Zp >> Ap >> Efermi40;

  file.close();
  file.clear();
  file.open("pca40.level");
  getline(file,out);

  double x40hole[60];
  double y40hole[60];
  double x40part[60];
  double y40part[60];
  int n40hole = 0;
  int n40part = 0;
  for(;;)
    {
      file >> i1 >> d1 >> i2 >> d2 >> d3 >> d4 >> d5 >> d6 >> d7 >> i3;
      if (file.eof())break;
      if (file.bad())break;
      if (d2 < Efermi40)
	{
         x40hole[n40hole] = d2 - Efermi40;
         y40hole[n40hole] = d4;
         n40hole++;
	}
      else
	{
         x40part[n40part] = d2 - Efermi40;
         y40part[n40part] = d4;
         n40part++;
	}
    }


  

  TGraph g40hole(n40hole,x40hole,y40hole);
  g40hole.SetMarkerStyle(21);
  g40hole.SetMarkerColor(2);
  g40hole.SetLineColor(2);
    g40hole.SetLineStyle(2);
    g40hole.SetMarkerSize(1.2);
  g40hole.Draw("PL");

  TGraph g40part(n40part,x40part,y40part);
  g40part.SetMarkerStyle(21);
  g40part.SetMarkerColor(2);
  g40part.SetLineColor(2);
  g40part.SetLineStyle(2);
    g40part.SetMarkerSize(1.2);
  g40part.Draw("PL");


  TLatex text;
  text.SetNDC();
  text.SetTextSize(.07);
  text.DrawLatex(.25,.31,"(a) protons");

  TLatex text2;
  text2.SetTextSize(.07);
  text2.SetTextColor(2);
  text2.DrawLatex(1.,.8,"^{40}Ca");
  text2.SetTextColor(4);
  text2.DrawLatex(-15,.67,"^{48}Ca");

  //**********************************************************
 pad1->cd();

  TH2S frame2("frame2","",10,-60,20,10,0,.99);
  frame2.GetXaxis()->SetTitle("E_{n,l,j}-E_{Fermi} [MeV]");
  frame2.GetXaxis()->CenterTitle();
  //frame2.GetYaxis()->SetLabelSize(0.);
  frame2.GetYaxis()->SetTitle("Occupancy");
  frame2.GetYaxis()->CenterTitle();
  frame2.SetStats(kFALSE);
  frame2.Draw();
  line.DrawLine(0.,0.,0.,1.);

  file.close();
  file.clear();
  file.open("nca48.inp");
  double EEfermi48;
  file >> Z >> Zp >> Ap >> EEfermi48;


  file.close();
  file.clear();
  file.open("nca48.level");
  getline(file,out);


  double xx48hole[60];
  double yy48hole[60];
  double xx48part[60];
  double yy48part[60];
  int nn48hole = 0;
  int nn48part = 0;
  for(;;)
    {
      file >> i1 >> d1 >> i2 >> d2 >> d3 >> d4 >> d5 >> d6 >> d7 >> i3;
      if (file.eof())break;
      if (file.bad())break;
      if (d2 < EEfermi48)
	{
         xx48hole[nn48hole] = d2 - EEfermi48;
         yy48hole[nn48hole] = d4;
         nn48hole++;
	}
      else
	{
         xx48part[nn48part] = d2 - EEfermi48;
         yy48part[nn48part] = d4;
         nn48part++;
	}

    }

  

  TGraph gg48hole(nn48hole,xx48hole,yy48hole);
  gg48hole.SetMarkerStyle(20);
  gg48hole.SetMarkerColor(4);
  gg48hole.SetLineColor(4);
    gg48hole.SetMarkerSize(1.2);
  gg48hole.Draw("PL");

  TGraph gg48part(nn48part,xx48part,yy48part);
  gg48part.SetMarkerStyle(20);
  gg48part.SetMarkerColor(4);
  gg48part.SetLineColor(4);
    gg48part.SetMarkerSize(1.2);
  gg48part.Draw("PL");


  file.close();
  file.clear();
  file.open("nca40.inp");
  double EEfermi40;
  file >> Z >> Zp >> Ap >> EEfermi40;

  file.close();
  file.clear();
  file.open("nca40.level");
  getline(file,out);

  double xx40hole[60];
  double yy40hole[60];
  double xx40part[60];
  double yy40part[60];
  int nn40hole = 0;
  int nn40part = 0;
  for(;;)
    {
      file >> i1 >> d1 >> i2 >> d2 >> d3 >> d4 >> d5 >> d6 >> d7 >> i3;
      if (file.eof())break;
      if (file.bad())break;
      if (d2 < EEfermi40)
	{
         xx40hole[nn40hole] = d2 - EEfermi40;
         yy40hole[nn40hole] = d4;
         nn40hole++;
	}
      else
	{
         xx40part[nn40part] = d2 - EEfermi40;
         yy40part[nn40part] = d4;
         nn40part++;
	}
    }


  

  TGraph gg40hole(nn40hole,xx40hole,yy40hole);
  gg40hole.SetMarkerStyle(21);
  gg40hole.SetMarkerColor(2);
  gg40hole.SetLineColor(2);
  gg40hole.SetLineStyle(2);
    gg40hole.SetMarkerSize(1.2);
  gg40hole.Draw("PL");

  TGraph gg40part(nn40part,xx40part,yy40part);
  gg40part.SetMarkerStyle(21);
  gg40part.SetMarkerColor(2);
  gg40part.SetLineColor(2);
  gg40part.SetLineStyle(2);
    gg40part.SetMarkerSize(1.2);
  gg40part.Draw("PL");

  text.DrawLatex(.25,.31,"(b) neutrons");


  TArrow arrow;
  arrow.SetAngle(30);
  arrow.SetLineColor(2);
  arrow.SetFillColor(2);
  arrow.SetLineStyle(2);
  pad1->cd();
  arrow.DrawArrow(15,.105,15,.820,.02,"<|>");
  arrow.SetLineStyle(1);
  arrow.DrawArrow(15,.105,15,.11,.02,"<|");
  arrow.DrawArrow(15,.81,15,.820,.02,"|>");
  arrow.SetLineColor(4);
  arrow.SetFillColor(4);

  arrow.DrawArrow(10,.162,10,.757,.02,"<|>");


}
