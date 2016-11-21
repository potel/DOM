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
Sty->SetPadLeftMargin(.15);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.01,"y");
Sty->SetLabelOffset(0.01,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.06,"xyz");
Sty->SetTitleOffset(1.3,"y");
Sty->SetTitleOffset(.8,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.03,"xz");
Sty->SetTickLength(.03,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("reactACa");

 TH2S frame("frame","",10,38,50,10,60,95);
 frame.GetXaxis()->SetTitle("A");
 frame.GetYaxis()->SetTitle("#sigma_{react}/A^{2/3} [mb]");
 frame.GetXaxis()->CenterTitle();
 frame.GetYaxis()->CenterTitle();

 frame.Draw();

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
 text.DrawLatex(40,1150,"E_{lab}=30 MeV");
 text.SetTextColor(3);
 text.DrawLatex(40,1100,"E_{lab}=35 MeV");

 text.SetTextColor(1);
 text.DrawLatex(45,1150,"data from Carlson");

}
