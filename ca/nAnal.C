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
 Sty->SetPadLeftMargin(.17);
Sty->SetHistLineWidth(1);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(1);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.07,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.09,"xyz");
Sty->SetTitleOffset(2.,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.05,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("nAnal","",400,600);


double All,isoScaler,isoVector,Volume,standard;
float lineWidth = 2.;
double ymax = 1300.;

 float Q;
 float x[200];
 float y[200];
 float sig[200];
 float xsig[200]={0.};
 float xx[200];
 float yy[200];



   TH2S frame ("frame","",10,0,180,10,0.,15.);
   frame->GetXaxis()->SetTitle("#theta_{cm} [deg]");
   frame->GetYaxis()->SetTitle("A");
   frame->GetYaxis()->CenterTitle();
   frame->GetYaxis()->SetTitleOffset(1. );
   frame->GetXaxis()->SetTickLength(.03);
   frame->GetYaxis()->SetTickLength(.03);
     frame->GetXaxis()->SetLabelOffset(-.0 );
   frame->GetXaxis()->SetTitleOffset( .7);
   frame->GetXaxis()->CenterTitle();     ;
 frame->Draw();



 TLine line;
 line.SetLineStyle(2);
 line.SetLineWidth(1);

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.08);

 text.DrawLatex(.6,.8,"n+^{40}Ca");

 ifstream fxdata("../ca/nca40a.dat");

 int Nth,Ndata;
 double Elab;
 TGraphErrors* graphE1[20];
 TGraph graph;
 graph.SetLineWidth(2.);
 bool first = 1;
 double scale = 1.;
 double scaleInc =2.;



 int ii = 0;
 for (;;)
   {
    
     fxdata >> Elab >> Ndata;
     if (fxdata.eof()) break;
     for (int i=0;i<Ndata;i++)
       {
	 fxdata >> x[i] >> y[i] >> sig[i];
         y[i] += scale;
       }
    
     first = 0;
     fxdata >> Nth;
     for (int i=0;i<Nth;i++)
       {
	 fxdata >> xx[i] >> yy[i] >> Q;
	 yy[i] += scale;
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


     graphE1[ii] = new TGraphErrors(Ndata,x,y,xsig,sig);
     graphE1[ii]->SetMarkerStyle(20);
     graphE1[ii]->SetMarkerSize(.6);

     if (Elab < 10)
       {
	 graphE1[ii]->SetMarkerColor(1.);
	 graphE1[ii]->SetLineColor(1.);
	 graph.SetLineColor(1.);
         line.SetLineColor(1.);
       }
     else if (Elab < 20)
       {
	 graphE1[ii]->SetMarkerColor(2.);
	 graphE1[ii]->SetLineColor(2.);
	 graph.SetLineColor(2.);
         line.SetLineColor(2.);
       }
     else if (Elab < 40)
       {
	 graphE1[ii]->SetMarkerColor(3.);
	 graphE1[ii]->SetLineColor(3.);
	 graph.SetLineColor(3.);
         line.SetLineColor(3.);
       }
     else if (Elab < 100)
       {
	 graphE1[ii]->SetMarkerColor(4.);
	 graphE1[ii]->SetLineColor(4.);
	 graph.SetLineColor(4.);
         line.SetLineColor(4.);
       }
     else 
       {
	 graphE1[ii]->SetMarkerColor(6.);
	 graphE1[ii]->SetLineColor(6.);
	 graph.SetLineColor(6.);
         line.SetLineColor(6.);
       }

     if (Ndata == 0) continue;
     cout << Ndata << " " << NthTry << endl;

     graphE1[ii]->Draw("P");
     graph.DrawGraph(NthTry,xx,yy,"L");
     line.DrawLine(0.,scale,180.,scale);
     scale += scaleInc;
     ii++;
   }


fxdata.close();
fxdata.clear();


 
}
