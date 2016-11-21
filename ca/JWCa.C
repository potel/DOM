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
Sty->SetPadBottomMargin(.16);
Sty->SetPadLeftMargin(.19);
Sty->SetPadTopMargin(.05);
Sty->SetPadRightMargin(.05);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.07,"xyz");
Sty->SetLabelOffset(0.01,"y");
Sty->SetLabelOffset(0.01,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
Sty->SetTitleOffset(1.3,"y");
Sty->SetTitleOffset(1.,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.03,"xz");
Sty->SetTickLength(.03,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("JwCa","",400,600);

//divide up canvas into three pads
 double overlap = .05;
 double dist = (1.+3*overlap)/4.;
 TPad *pad1  = new TPad("pad1","",0.,0.,1.,.5+overlap);
 TPad *pad2  = new TPad("pad2","",0.,.5-overlap,1.,1.);
pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

gPad->GetFrame()->SetLineWidth(1);
pad1->Draw();
pad2->Draw();

 pad2->cd();

 TH2S frame("frame","",10,110,126,10,100,140);
 //frame.GetXaxis()->SetTitle("A");
 frame.GetXaxis()->SetLabelSize(0.);
 frame.GetXaxis()->CenterTitle();
 frame.GetYaxis()->CenterTitle();
 frame.GetYaxis()->SetTitle("|J_{W}/A| [MeV fm^{3}]");
 frame.Draw();


 double x1[6]={112,116,118,120,122,124};

 //20 MeV data
 double y20[6]={114.1,116.1,123.7,130.4,134.2};
 double y16[6] = {114.142,116.519,127.682,123.592,125.319,129.578};
  TGraph g20(6,x1,y20);
 g20.SetMarkerStyle(21);
 g20.SetMarkerColor(2);
 g20.SetMarkerSize(2);
 g20.Draw("P");

  TGraph g16(6,x1,y16);
 g16.SetMarkerStyle(20);
 g16.SetMarkerColor(4);
 g16.SetMarkerSize(2);
 //g16.Draw("P");


 TLatex text;
 text.SetNDC();
 text.SetTextSize(.08);
 text.DrawLatex(.22,.25,"(a) Makofske et al. (1968)");

 text.SetTextColor(4);
 //text.DrawLatex(.6,.4,"E_{lab}=16 MeV");
 text.SetTextColor(1);
 text.DrawLatex(.22,.75,"E_{lab}=20 MeV");

 pad1->cd();
 TH2S frame("frame","",10,38,50,10,60,95);
 frame1.GetXaxis()->SetTitle("A");
 frame1.GetXaxis()->CenterTitle();
 frame1.GetYaxis()->CenterTitle();
 frame1.GetYaxis()->SetTitle("#sigma_{react}/A^{2/3} [mb]");
 frame1.Draw();

 //ELab = 25.3,24.8,24.9 Carlson
 double x1[4]={40,42,44,48};
 double y1[4]={876.,969.,1028,1058};
 double sx1[4] = {0.};
 double sy1[4]={33.,32.,28.,47.};
 //ELab = 30.3,29.7,30
 double y3[4] = {879.,934.,1002,999.};
 double sy3[4] = {26.,26.,25.,55};
 //Elab = 35
 double y4[4] = {853.,916.,951.,971.};
 double sy4[4] = {25.,25.,21.,32.};


 for (int i=0;i<4;i++)
   {
     y1[i] /= pow(x1[i],(2./3.));
     sy1[i] /= pow(x1[i],(2./3.));
     y3[i] /= pow(x1[i],(2./3.));
     sy3[i] /= pow(x1[i],(2./3.));
     y4[i] /= pow(x1[i],(2./3.));
     sy4[i] /= pow(x1[i],(2./3.));

   }


 TGraphErrors g1(4,x1,y1,sx1,sy1);
 g1.SetMarkerStyle(21);
 g1.SetMarkerSize(2);
 g1.Draw("PL");

 TGraphErrors g3(4,x1,y3,sx1,sy3);
 g3.SetMarkerStyle(20);
 g3.SetMarkerSize(2);
 g3.SetMarkerColor(4);
 g3.SetLineColor(4);
 //g3.Draw("PL");

 TGraphErrors g4(4,x1,y4,sx1,sy4);
 g4.SetMarkerStyle(20);
 g4.SetMarkerSize(2);
 g4.SetMarkerColor(3);
 g4.SetLineColor(3);
 //g4.Draw("PL");



 TLatex text;

 text.SetTextColor(1);
 text.DrawLatex(40,90,"E_{lab}#approx25 MeV");
 text.SetTextColor(4);
 //text.DrawLatex(40,1150,"E_{lab}=30 MeV");
 text.SetTextColor(3);
 //text.DrawLatex(40,1100,"E_{lab}=35 MeV");

 text.SetTextColor(1);
 text.DrawLatex(45,1150,"data from Carlson");



 text.DrawLatex(.22,.22,"(b) Carlson et al. (1995)");
}
