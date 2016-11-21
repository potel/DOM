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
Sty->SetPadBottomMargin(.15);
Sty->SetPadLeftMargin(.15);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.06,"xyz");
Sty->SetTitleOffset(1.2,"y");
Sty->SetTitleOffset(1.1,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.02,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("SpectF");

  TH2S frame("frame","",10,-80,30,10,0,1);
  frame.GetXaxis()->SetTitle("E [MeV]");
  frame.GetYaxis()->SetTitle("Spectroscpic Function");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.Draw();



  
  ifstream file("d32.dat");
  double x[500];
  double y[500];
  double xx,yy;
  int N=0;
  for (;;)
    {
      file >> xx >> yy;
      if (file.eof())break;
      if (file.bad())break;
      x[N] = xx;
      y[N] = yy*30;
      N++;
    }

  TGraph graph(N,x,y);
  
  graph.Draw("L same");

  float dd = .25;
  TLine line;
  line.DrawLine(-8.5-dd,0,-8.5-dd,.7);
  line.DrawLine(-8.5+dd,0,-8.5+dd,.7);
  line.DrawLine(-8.5-dd,.7,-8.5+dd,.7);

  line.SetLineStyle(2);
  line.DrawLine(-4.71,0.,-4.71,1.);

  TLatex text;
  text.SetNDC();
  text.DrawLatex(.18,.75,"0d#frac{3}{2} protons in ^{40}Ca");
  text.DrawLatex(.67,.8,"E_{Fermi}");

  return;
  TArrow arrow;
  arrow.SetAngle(30);
  arrow.SetLineWidth(1.);
  arrow.SetFillColor(2);
  arrow.DrawArrow(1,.3,-15,.2,.02,"|>");




}
