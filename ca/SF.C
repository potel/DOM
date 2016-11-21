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
Sty->SetPadBottomMargin(.25);
Sty->SetPadLeftMargin(.2);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
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

 TCanvas canvas("SP");



 TH2S frame("frame","",10,-19,0,10,0,1);
 frame.GetXaxis()->SetTitle("E_{njl} [MeV]");
 frame.GetYaxis()->SetTitle("S_{njl} ");

 frame.GetXaxis()->CenterTitle();
 frame.GetYaxis()->CenterTitle();
 frame.Draw();


 TLatex text;
 text.SetTextColor(4);
 text.SetNDC();



 TH2S frame2("frame2","",10,-19,0,10,0,1);
 frame2.GetXaxis()->SetTitle("E_{njl} [MeV]");
 frame2.GetYaxis()->SetTitle("Spect. Factor");


 //frame2.GetYaxis()->SetLabelSize(0.);

 frame2.GetXaxis()->CenterTitle();
 frame2.Draw();

 ifstream file;
 string out;
 file.open("/home/charity/DOM48/ca/pca48.lev");
 getline(file,out);
 cout << out << endl;
 getline(file,out);

 double sCa48[10];
 double eCa48[10];
 double esCa48[10];
 double eeCa48[10];

 double one,two,three,four,five,six;
 double energy,spect,error;
 int i1,i2,i3,i4,i5,i6,i7;

 int nCa48 = 0;
 for (;;)
   {
     file >> energy >> one >> i1 >> two >> i2 >> i3 >> i4 >> three >> four >>
       i5 >> five >> six >> i6 >> spect >> error >> i7;
     if (file.bad())break;
     if (file.eof())break;
     if (spect > 0.)
       {
        eCa48[nCa48] = energy;
        sCa48[nCa48] = spect;
        esCa48[nCa48] = error;
        eeCa48[nCa48] = 0.;
        nCa48++;
       }
   }

 

 TGraphErrors gCa48(nCa48,eCa48,sCa48,eeCa48,esCa48);
 gCa48.SetMarkerStyle(20);
 gCa48.SetMarkerColor(2);
 gCa48.SetLineColor(2);
 gCa48.Draw("P");






 file.close();
 file.clear();
 file.open("/home/charity/DOM48/ca/pca48.level");
 getline(file,out);

 double eeeCa48[30];
 double ssCa48[30]; 
 int nnCa48 = 0;
 for(;;)
   {
     file >> i1 >> five >> i1 >> energy >> one >> two >> spect >> three >> four >> i3;
     if (file.eof())break;
     if (file.bad())break;

    
     if (energy < 0)
       {
        eeeCa48[nnCa48] = energy;
        ssCa48[nnCa48] = spect;
        nnCa48++;
       }
   }

 TGraph ggCa48(nnCa48,eeeCa48,ssCa48);
 ggCa48.SetLineColor(2);
 ggCa48.Draw("L");

 text.SetTextColor(2);
 text.DrawLatex(.3,.65,"^{48}Ca");





 file.close();
 file.clear();
 file.open("/home/charity/DOM48/ca/pca40.lev");
 getline(file,out);
 cout << out << endl;
 getline(file,out);

 double sCa40[10];
 double eCa40[10];
 double esCa40[10];
 double eeCa40[10];

 int nCa40 = 0;
 for (;;)
   {
     file >> energy >> one >> i1 >> two >> i2 >> i3 >> i4 >> three >> four >>
       i5 >> five >> six >> i6 >> spect >> error >> i7;
     if (file.bad())break;
     if (file.eof())break;
     if (spect > 0.)
       {
        eCa40[nCa40] = energy;
        sCa40[nCa40] = spect;
        esCa40[nCa40] = error;
        eeCa40[nCa40] = 0.;
        nCa40++;
       }
   }

 

 TGraphErrors gCa40(nCa40,eCa40,sCa40,eeCa40,esCa40);
 gCa40.SetMarkerStyle(20);
 gCa40.SetMarkerColor(3);
 gCa40.SetLineColor(3);
 gCa40.Draw("P");

 file.close();
 file.clear();
 file.open("/home/charity/DOM48/ca/pca40.level");
 getline(file,out);

 double eeeCa40[30];
 double ssCa40[30]; 
 int nnCa40 = 0;
 for(;;)
   {
     file >> i1 >> five >> i1 >> energy >> one >> two >> spect >> three >> four >> i3;
     if (file.eof())break;
     if (file.bad())break;

    
     if (energy < 100.)
       {
        eeeCa40[nnCa40] = energy;
        ssCa40[nnCa40] = spect;
        nnCa40++;
       }
   }

 TGraph ggCa40(nnCa40,eeeCa40,ssCa40);
 ggCa40.SetLineColor(3);
 ggCa40.Draw("L");

 text.SetTextColor(3);
 text.DrawLatex(.7,.65,"^{40}Ca");


}
