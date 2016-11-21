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
Sty->SetPadBottomMargin(.12);
Sty->SetPadTopMargin(.05);
Sty->SetPadLeftMargin(.25);
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
Sty->SetTitleOffset(1.4,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.05,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("nAngCa","",600,600);

//divide up canvas into three pads
 double overlap = .15;
 double dist = (1.+3*overlap)/4.;
 TPad *pad1  = new TPad("pad1","",0.,0.,.5+overlap/2.,1.);
 TPad *pad2  = new TPad("pad2","",.5-overlap/2.,0.,1.,1.);

pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

pad1->Draw();
pad2->Draw();



double All,isoScaler,isoVector,Volume,standard;
float lineWidth = 2.;
double ymax = 1300.;


 float x[200];
 float y[200];
 float sig[200];
 float xsig[200]={0.};
 float xx[200];
 float yy[200];


 pad1->cd();
   pad1->SetLogy();
   TH2S frame ("frame","",10,0,180,10,1.,1.e22);
   //frame->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 frame->GetYaxis()->SetTitle("d#sigma/d#Omega");
 frame->GetYaxis()->CenterTitle();
 frame->GetXaxis()->SetTickLength(.03);
 frame->GetYaxis()->SetTickLength(.03);
 frame->Draw();

pad2->cd();
   pad2->SetLogy();
   TH2S frame2 ("frame2","",10,0,180,10,1.,1.e22);
   //frame2->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 //frame2->GetYaxis()->SetTitle("d#sigma/d#Omega");
 frame2->GetYaxis()->SetLabelSize(0.);
 frame2->GetXaxis()->SetTickLength(.03);
 frame2->GetYaxis()->SetTickLength(.03);
 frame2->Draw();

 pad1->cd();

 TLatex text;
 text.SetNDC();
 text.SetTextSize(.08);
 text.DrawLatex(.68,.8,"n+^{40}Ca");
 pad2->cd();
 text.DrawLatex(.4,.8,"n+^{40}Ca");
 pad1->cd();
 ifstream fxdata("nca40x.dat");
 ifstream fcdata("nca40c.dat");

 int Nth,Ndata;
 double Elab;
 TGraphErrors graphE1[30];
 TGraph graph;
 graph.SetLineWidth(1.5);
 bool first = 1;
 double scale = 1.;
 double scaleInc =30.;

 bool flag =1;

 int ii = 0;
 for (;;)
   {
    
     fcdata >> Elab >> Nth;
     if (flag && Elab >= 70.) 
       {
	 scale = 100000000.;
	flag = 0;
       }
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


     //remove larger error

     int nn = Ndata;
     int shift = 0;
     for (int i=0;i<nn;i++)
       {
         if (sig[i] >= y[i]) 
	   {
             cout << "caught" << endl;
	     Ndata--;
             shift++;
	   }
         else if (shift > 0)
	   {
             sig[i-shift] = sig[i];
             y[i-shift] = y[i];
             x[i-shift] = x[i];
	   }
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
       }
     else if (Elab < 20)
       {
	 graphE1[ii]->SetMarkerColor(2.);
	 graphE1[ii]->SetLineColor(2.);
	 graph.SetLineColor(2.);
       }
     else if (Elab < 40)
       {
	 graphE1[ii]->SetMarkerColor(3.);
	 graph.SetLineColor(3.);
       }
     else if (Elab < 100)
       {
	 graphE1[ii]->SetMarkerColor(4.);
	 graphE1[ii]->SetLineColor(4.);
	 graph.SetLineColor(4.);
       }
     else 
       {
	 graphE1[ii]->SetMarkerColor(6.);
	 graphE1[ii]->SetLineColor(6.);
	 graph.SetLineColor(6.);
       }

     if (Ndata == 0) continue;

     if (Elab < 70) pad1->cd();
     else pad2->cd();

     graphE1[ii].SetMarkerSize(.7);
     graphE1[ii].Draw("P");
     graph.DrawGraph(NthTry,xx,yy,"L");

     scale *= scaleInc;
     ii++;
   }


fxdata.close();
fxdata.clear();
fcdata.close();
fcdata.clear();





//******************************************************
    //pad3->cd();
    //pad3->SetLogy();
    /*
   TH2S frame3 ("frame3","",10,0,180,10,1.,1.e20);
   //frame3->GetXaxis()->SetTitle("#theta_{cm} [deg]");
 //frame3->GetYaxis()->SetTitle("d#sigma/d#Omega");
 frame3->GetYaxis()->SetLabelSize(0.);
 frame3->GetYaxis()->SetTitleOffset(1.5);
 frame3->GetXaxis()->SetTickLength(.03);
 frame3->GetYaxis()->SetTickLength(.03);
 frame3->GetXaxis()->SetLabelOffset(.01);
 frame3->Draw();
    */

 text.DrawLatex(.68,.3,"n+^{48}Ca");

 fxdata.open("../ca/nca48x.dat");
 fcdata.open("../ca/nca48c.dat");

 TGraphErrors graphE3[20];
 graph.SetLineWidth(1.5);



 first = 1;
 scale = 1.;
 ii = 0;
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


 canvas.cd();
 text.SetTextSize(.05);
 text.DrawLatex(.45,.03,"#theta_{cm} [deg]");

}
