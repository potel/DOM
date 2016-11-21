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
Sty->SetPadLeftMargin(.3);
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

 TCanvas canvas("pAnalCa42","",600,600);


//divide up canvas into three pads
 double overlap = .1;
 double dist = (1.+3*overlap)/4.;
 TPad *pad1  = new TPad("pad1","",0.,0.,.5+overlap,1.);
 TPad *pad2  = new TPad("pad2","",.5-overlap,0.,1.,1.);


pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);


gPad->GetFrame()->SetLineWidth(1);
pad1->Draw();
pad2->Draw();




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


 pad1->cd();
   TH2S frame ("frame","",10,0,180,10,0.,15.);
 frame->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 frame->GetYaxis()->SetTitle("A");
 frame->GetYaxis()->CenterTitle();
 frame->GetYaxis()->SetTitleOffset(1.5);
 frame->GetXaxis()->SetTickLength(.03);
 frame->GetYaxis()->SetTickLength(.03);
 frame->GetXaxis()->SetLabelOffset(-.01);
 frame->Draw();



 TLine line;
 line.SetLineStyle(2);
 line.SetLineWidth(1);

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.08);
 text.DrawLatex(.6,.8,"p+^{42}Ca");

 ifstream fxdata("pca42a.dat");

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


//**********************************************************
    cout << "pad2" << endl;
 pad2->cd();
   TH2S frame2 ("frame2","",10,0,180,10,0.,15.);
 frame2->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 frame2->GetYaxis()->SetLabelSize(0.);
 frame2->GetYaxis()->CenterTitle();
 frame2->GetYaxis()->SetTitleOffset(1.5);
 frame2->GetXaxis()->SetTickLength(.03);
 frame2->GetYaxis()->SetTickLength(.03);
 frame2->GetXaxis()->SetLabelOffset(-.01);
 frame2->Draw();



 text.DrawLatex(.6,.8,"p+^{44}Ca");

 ifstream fxdata("../ca/pca44a.dat");

 TGraphErrors* graphE2[20];
 first = 1;
 scale = 1.;


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

     graphE2[ii] = new TGraphErrors(Ndata,x,y,xsig,sig);
     graphE2[ii]->SetMarkerStyle(20);
     graphE2[ii]->SetMarkerSize(.6);

     if (Elab < 10)
       {
	 graphE2[ii]->SetMarkerColor(1.);
	 graphE2[ii]->SetLineColor(1.);
	 graph.SetLineColor(1.);
	 line.SetLineColor(1);
       }
     else if (Elab < 20)
       {
	 graphE2[ii]->SetMarkerColor(2.);
	 graphE2[ii]->SetLineColor(2.);
	 graph.SetLineColor(2.);
	 line.SetLineColor(2);
       }
     else if (Elab < 40)
       {
	 graphE2[ii]->SetMarkerColor(3.);
	 graphE2[ii]->SetLineColor(3.);
	 graph.SetLineColor(3.);
	 line.SetLineColor(3);
       }
     else if (Elab < 100)
       {
	 graphE2[ii]->SetMarkerColor(4.);
	 graphE2[ii]->SetLineColor(4.);
	 graph.SetLineColor(4.);
	 line.SetLineColor(4);
       }
     else 
       {
	 graphE2[ii]->SetMarkerColor(6.);
	 graphE2[ii]->SetLineColor(6.);
	 graph.SetLineColor(6.);
	 line.SetLineColor(6);
       }

     if (Ndata == 0) continue;
     cout << Ndata << " " << NthTry << endl;

     graph.DrawGraph(NthTry,xx,yy,"L");
     graphE2[ii]->Draw("P");
     line.DrawLine(0.,scale,180.,scale);

     scale += scaleInc;
     ii++;
   }


fxdata.close();
fxdata.clear();



 canvas.cd();
 text.SetTextSize(.05);
 text.DrawLatex(.45,.04,"#theta_{cm} [deg]");
 
}
