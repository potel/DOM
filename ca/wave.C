{gROOT->Reset();
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
Sty->SetPadBottomMargin(.18);
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
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("wave");



  TH2S frame("frame","",10,0,8,10,0,.45);
  frame.GetXaxis()->SetTitle("r [fm]");
  frame.GetYaxis()->SetTitle("|u_{n,l,j}|^{2} [fm^{-1}]");

  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.Draw();
  double xs[300];
  double ys[300];
  int Ns=0;

  ifstream file("s12.dat");
  double r,phi;
  for (;;)
    {
      file >> r >> phi;
      if (file.eof())break;
      if(file.bad())break;
      xs[Ns] = r;
      ys[Ns] = phi*phi;
      Ns++;
    }

  TGraph gs(Ns,xs,ys);
  gs.SetLineColor(2);
  gs.SetLineWidth(2);
  gs.Draw("C");

  file.close();
  file.clear();
  file.open("f72.dat");
  double xf[300];
  double yf[300];
  int Nf=0;

  ifstream file("f72.dat");

  for (;;)
    {
      file >> r >> phi;
      if (file.eof())break;
      if(file.bad())break;
      xf[Nf] = r;
      yf[Nf] = phi*phi;
      Nf++;
    }

  TGraph gf(Nf,xf,yf);
  gf.SetLineColor(4);
  gf.SetLineWidth(2);
  gf.Draw("C");

  TF1 *funct = new TF1("funct","[2]*exp((x-[0])/[1])/pow(1.+exp((x-[0])/[1]),2)",0,15);
  funct->SetParameter(0,4.248);
  funct->SetParameter(1,.6);
  funct->SetParameter(2,1.2);
  funct->SetLineStyle(9);
  funct->Draw("same");

  TLatex text;
  text.SetTextColor(4);
  text.DrawLatex(4.3,.4,"#nu f7/2");
  text.SetTextColor(2);
   text.DrawLatex(2.7,.04,"#pi s1/2");

   text.SetTextColor(3);
   text.DrawLatex(5.2,.2,"W_{surface}");

   text.SetTextColor(1);
   text.SetTextSize(.08);
   text.SetNDC();
   text.DrawLatex(.75,.8,"^{48}Ca");
}
