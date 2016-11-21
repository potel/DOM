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
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(.01,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.08,"xyz");
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

 TCanvas canvas("pRuthCa42","",800,600);


//divide up canvas into three pads
 double overlap = .15;
 double dist = (1.+3*overlap)/4.;
 TPad *pad2  = new TPad("pad2","",0.,0.,.5+overlap/2.,1.);
 TPad *pad3  = new TPad("pad3","",.5-overlap/2.,0.,1.,1.);
pad2->SetFillStyle(4000);
pad3->SetFillStyle(4000);
gPad->GetFrame()->SetLineWidth(1);
pad2->Draw();
pad3->Draw();


double All,isoScaler,isoVector,Volume,standard;
float lineWidth = 2.;
double ymax = 1300.;


 float x[200];
 float y[200];
 float sig[200];
 float xsig[200]={0.};
 float xx[200];
 float yy[200];



 pad2->cd();
   pad2->SetLogy();
   TH2S frame ("frame","",10,0,180,10,.1,1.e27);
   //frame->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 frame->GetYaxis()->SetTitle("#sigma/#sigma_{Ruth} ");
 //frame->GetYaxis()->SetLabelSize(0.);
 frame->GetYaxis()->CenterTitle();
 frame->GetXaxis()->SetTickLength(.03);
 frame->GetYaxis()->SetTickLength(.03);
 frame->GetXaxis()->SetLabelOffset(.01);
 frame->Draw();

 double scale = 2.5;
 double scaleInc =scale;

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.08);
 text.DrawLatex(.68,.8,"^{42}Ca");

 ifstream fxdata;
 fxdata.open("../ca/pca42xx.dat");

 int Nth,Ndata;
 double Elab;
 TGraph graphE;
 TGraph graph;

 graphE.SetMarkerStyle(20);
 graphE.SetMarkerSize(.5);
 bool first = 1;
 for (;;)
   {
     fxdata >> Nth;
     if (fxdata.eof()) break;
     for (int i=0;i<Nth;i++)
       {
	 fxdata >> xx[i] >> yy[i];
         xx[i] *= 180./3.14159;

       }

     first = 0;
     fxdata >> Ndata >> Elab;
     for (int i=0;i<Ndata;i++)
       {
	 fxdata >> x[i] >> y[i] >> sig[i];
         x[i] *= 180./3.14159;
       }
     xmax = x[Ndata-1];
     int NthTry = Nth;
     for (int i=0;i<Nth;i++)
       {
         if (xx[i] > xmax)
	   {
             NthTry = i + 1;
	     break;
	   }
       }

    if (Elab > 90.)
       {
	 cout << "scale" << endl;
	 for (int i=0;i<Ndata;i++) y[i] *=scale*scale*4.;
	 for (int i=0;i<NthTry;i++) yy[i] *=scale*scale*4.;
	 scale *= scaleInc;
        
       }

    if (Elab < 10)
       {
	 graphE.SetMarkerColor(1.);
	 graph.SetLineColor(1.);
       }
     else if (Elab < 20)
       {
	 graphE.SetMarkerColor(2.);
	 graph.SetLineColor(2.);
       }
     else if (Elab < 40)
       {
	 graphE.SetMarkerColor(3.);
	 graph.SetLineColor(3.);
       }
     else if (Elab < 100)
       {
	 graphE.SetMarkerColor(4.);
	 graph.SetLineColor(4.);
       }
     else 
       {
	 graphE.SetMarkerColor(6.);
	 graph.SetLineColor(6.);
       }

     graphE.DrawGraph(Ndata,x,y,"P");
     graph.DrawGraph(NthTry,xx,yy,"L");
   }


fxdata.close();
fxdata.clear();


//**********************************************************
 pad3->cd();
   pad3->SetLogy();
   TH2S frame2 ("frame2","",10,0,180,10,.1,1.e27);
   //frame2->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 //frame2->GetYaxis()->SetTitle("#sigma/#sigma_{Ruth} ");
 frame2->GetYaxis()->SetLabelSize(0.);
 frame2->GetYaxis()->CenterTitle();
 frame2->GetXaxis()->SetTickLength(.03);
 frame2->GetYaxis()->SetTickLength(.03);
 frame2->GetXaxis()->SetLabelOffset(.01);
 frame2->Draw();


 text.DrawLatex(.68,.8,"^{44}Ca");

 fxdata.open("../ca/pca44xx.dat");

 scale = scaleInc;
 graphE.SetMarkerStyle(20);
 graphE.SetMarkerSize(.5);
 bool first = 1;
 for (;;)
   {
     fxdata >> Nth;
     if (fxdata.eof()) break;
     for (int i=0;i<Nth;i++)
       {
	 fxdata >> xx[i] >> yy[i];
         xx[i] *= 180./3.14159;

       }

     first = 0;
     fxdata >> Ndata >> Elab;
     for (int i=0;i<Ndata;i++)
       {
	 fxdata >> x[i] >> y[i] >> sig[i];
         x[i] *= 180./3.14159;
       }
     xmax = x[Ndata-1];
     int NthTry = Nth;
     for (int i=0;i<Nth;i++)
       {
         if (xx[i] > xmax)
	   {
             NthTry = i + 1;
	     break;
	   }
       }


    if (Elab >= 10.)
       {
	 for (int i=0;i<Ndata;i++) y[i] *=100.;
	 for (int i=0;i<NthTry;i++) yy[i] *=100.;
       }


    if (Elab >= 20.)
       {
	 for (int i=0;i<Ndata;i++) y[i] *=1000.;
	 for (int i=0;i<NthTry;i++) yy[i] *=1000.;
       }

    if (Elab >= 40.)
       {
	 for (int i=0;i<Ndata;i++) y[i] *=1000.;
	 for (int i=0;i<NthTry;i++) yy[i] *=1000.;
       }



    if (Elab >= 100.)
       {
	 for (int i=0;i<Ndata;i++) y[i] *=scale*100.;
	 for (int i=0;i<NthTry;i++) yy[i] *=scale*100.;
	 scale *= scaleInc;
        
       }
     else 
       {
	 graphE.SetMarkerColor(6.);
	 graph.SetLineColor(6.);
       }

    if (Elab < 10)
       {
	 graphE.SetMarkerColor(1.);
	 graph.SetLineColor(1.);
       }
     else if (Elab < 20)
       {
	 graphE.SetMarkerColor(2.);
	 graph.SetLineColor(2.);
       }
     else if (Elab < 40)
       {
	 graphE.SetMarkerColor(3.);
	 graph.SetLineColor(3.);
       }
     else if (Elab < 100)
       {
	 graphE.SetMarkerColor(4.);
	 graph.SetLineColor(4.);
       }
     else if (Elab > 100)
       {
	 graphE.SetMarkerColor(6.);
	 graph.SetLineColor(6.);
       }

     graphE.DrawGraph(Ndata,x,y,"P");
     graph.DrawGraph(NthTry,xx,yy,"L");
   }


fxdata.close();
fxdata.clear();





 canvas.cd();
 text.SetTextSize(.06);
   text.DrawLatex(.47,.02,"#theta_{cm} [deg]");

   pad3.cd();
   text.SetTextSize(.04);
   text.DrawLatex(.4,.2,"E_{lab}<10 MeV");  
   text.SetTextColor(2);
   text.DrawLatex(.4,.36,"10<E_{lab}<20 MeV");  

   text.SetTextColor(3);
   text.DrawLatex(.4,.53,"20<E_{lab}<40 MeV");  

   text.SetTextColor(4);
   text.DrawLatex(.57,.74,"40<E_{lab}<100 MeV"); 

   pad2->cd();
   text.SetTextColor(6);
   //text.DrawLatex(.6,.7,"E_{lab}>100 MeV"); 
}
