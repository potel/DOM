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
 Sty->SetPadLeftMargin(.17);
Sty->SetPadTopMargin(.05);
 Sty->SetPadRightMargin(.05);
Sty->SetHistLineWidth(1);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(1);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.07,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(-0.01,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.09,"xyz");
Sty->SetTitleOffset(1.,"y");
Sty->SetTitleOffset(.6,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.05,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("nAnalCaFe","",600,600);


//divide up canvas into three pads
 double overlap = .1;
 double dist = (1.+3*overlap)/4.;
 TPad *pad2  = new TPad("pad2","",0.,0.,.5+overlap/2.,1.);
 TPad *pad3  = new TPad("pad3","",.5-overlap/2.,0.,1.,1.);

pad2->SetFillStyle(4000);
pad3->SetFillStyle(4000);

gPad->GetFrame()->SetLineWidth(1);
pad2->Draw();
pad3->Draw();

 pad2->cd();
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



   TH2S frame ("frame","",10,0,180,10,0.,11.);
   frame->GetXaxis()->SetTitle("#theta_{cm} [deg]");
   frame->GetYaxis()->SetTitle("A");
   frame->GetYaxis()->CenterTitle();
   frame->GetXaxis()->SetTickLength(.03);
   frame->GetYaxis()->SetTickLength(.03);
   frame->GetXaxis()->CenterTitle();     ;
 frame->Draw();



 TLine line;
 line.SetLineStyle(2);
 line.SetLineWidth(1);

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.08);

 text.DrawLatex(.3,.88,"(a) n+^{40}Ca");

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

 pad3->cd();
   TH2S frame2 ("frame2","",10,0,180,10,0.,11.);
   frame2->GetXaxis()->SetTitle("#theta_{cm} [deg]");
   //frame2->GetYaxis()->SetTitle("A");
   frame2->GetYaxis()->CenterTitle();
   frame2->GetXaxis()->SetTickLength(.03);
   frame2->GetYaxis()->SetTickLength(.03);
   frame2->GetXaxis()->CenterTitle();     ;
   frame2->GetYaxis()->SetLabelSize(0.);
 frame2->Draw();


 text.DrawLatex(.3,.88,"(b) n+^{54}Fe");

 ifstream fxdata("../fe/nfe54a.dat");

 graph.SetLineWidth(2.);
 first = 1;
 scale = 1.;
 scaleInc =2.;



 ii = 0;
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
