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
Sty->SetPadLeftMargin(.2);
Sty->SetHistLineWidth(1);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(1);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.05,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
Sty->SetTitleOffset(1.,"y");
Sty->SetTitleOffset(1.,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.05,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("doe","",400,600);


//divide up canvas into three pads
 double overlap = .11;

double All,isoScaler,isoVector,Volume,standard;
float lineWidth = 2.;
double ymax = 1300.;


 float x[200];
 float y[200];
 float sig[200];
 float xsig[200]={0.};
 float xx[200];
 float yy[200];



   gPad->SetLogy();
   TH2S frame ("frame","",10,0,180,10,2.,1.e6);
   frame->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 frame->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");
 frame->GetYaxis()->CenterTitle();
 frame->GetXaxis()->CenterTitle();
 frame->GetYaxis()->SetTitleOffset(1.3);
 frame->GetXaxis()->SetTickLength(.03);
 frame->GetYaxis()->SetTickLength(.03);
 frame->GetXaxis()->SetLabelOffset(.01);
 frame->Draw();

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.1);
 text.DrawLatex(.5,.8,"n+^{48}Ca");

 ifstream fxdata;
 ifstream fcdata;

 fxdata.open("../ca/nca48x.dat");
 fcdata.open("../ca/nca48c.dat");

 TGraphErrors graphE3[20];
 TGraph graph;
 graph.SetLineWidth(1.5);



 int Nth,Ndata;
 double Elab;
 bool first = 1;
 double scale = 1.;
 double scaleInc =30.;
 int ii = 0;
 for (;;)
   {
    
     fcdata >> Elab >> Nth;
     if (fcdata.eof()) break;
     for (int i=0;i<Nth;i++)
       {
	 fcdata >> xx[i] >> yy[i];
         yy[i] *= scale;
       }

    
     first = 0;
     fxdata >> Elab >> Ndata;
     for (int i=0;i<Ndata;i++)
       {
	 fxdata >> x[i] >> y[i] >> sig[i];
	 y[i] *= scale;
	 sig[i] *= scale;
       }
     double xmax = x[Ndata-1];
     int NthTry = Nth;


     for (int i=0;i<Nth;i++)
       {
         if (xx[i] > xmax)
	   {
             NthTry = i + 1;
	     break;
	   }
       }

     graphE3[ii] = new TGraphErrors(Ndata,x,y,xsig,sig);
     graphE3[ii]->SetMarkerStyle(20);
     graphE3[ii]->SetMarkerSize(.7);

     if (Elab < 10)
       {
	 graphE3[ii]->SetMarkerColor(1.);
	 graphE3[ii]->SetLineColor(1.);
	 graph.SetLineColor(1.);
       }
     else if (Elab < 20)
       {
	 graphE3[ii]->SetMarkerColor(2.);
	 graphE3[ii]->SetLineColor(2.);
	 graph.SetLineColor(2.);
       }
     else if (Elab < 40)
       {
	 graphE3[ii]->SetMarkerColor(3.);
	 graphE3[ii]->SetLineColor(3.);
	 graph.SetLineColor(3.);
       }
     else if (Elab < 100)
       {
	 graphE3[ii]->SetMarkerColor(4.);
	 graphE3[ii]->SetLineColor(4.);
	 graph.SetLineColor(4.);
       }
     else 
       {
	 graphE3[ii]->SetMarkerColor(6.);
	 graphE3[ii]->SetLineColor(6.);
	 graph.SetLineColor(6.);
       }




     graphE3[ii]->Draw("P");
     graph.DrawGraph(NthTry,xx,yy,"L");

     scale *= scaleInc;
     ii++;
   }


fxdata.close();
fxdata.clear();
fcdata.close();
fcdata.clear();

 text.SetTextSize(.05);
 text.DrawLatex(.6,.35,"8.0 MeV");
 text.SetTextColor(2);
 text.DrawLatex(.55,.54,"11.9 MeV(x30)");
 text.DrawLatex(.55,.73,"16.8 MeV(x900)");

}
