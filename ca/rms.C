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

 TCanvas canvas("rms");



 string out;


 TH2S frame("frame","",10,-25,-5,10,2,7);
 frame.GetXaxis()->SetTitle("E_{njl} [MeV]");
 frame.GetYaxis()->SetTitle("R^{rms}_{njl} [fm]");

 frame.GetXaxis()->CenterTitle();
 frame.GetYaxis()->CenterTitle();
 frame.Draw();



 double one,two,three,four,five,six;
 double energy,spect,error;
 int i1,i2,i3,i4,i5,i6,i7;


 TLatex text;
 text.SetTextColor(4);
 text.SetNDC();


 ifstream  file;
 file.open("/home/charity/DOM48/ca/pca48.lev");
 getline(file,out);
 cout << out << endl;
 getline(file,out);

 double sCa48[10];
 double eCa48[10];
 double esCa48[10];
 double eeCa48[10];

 int nCa48 = 0;
 for (;;)
   {
     file >> energy >> one >> i1 >> two >> i2 >> i3 >> i4 >> three >> four >>
       i5 >> five >> six >> i6 >> spect >> error >> i7;
     if (file.bad())break;
     if (file.eof())break;
     if (three > 0.)
       {

        eCa48[nCa48] = energy;
        sCa48[nCa48] = three;
        esCa48[nCa48] = four;
        eeCa48[nCa48] = 0.;
        nCa48++;
       }
   }

 
 cout << nCa48 << endl;
 cout << eCa48[0] << " " << eCa48[1] << " " << eCa48[2] << endl;
 cout << sCa48[0] << " " << sCa48[1] << " " << sCa48[2] << endl;
 TGraphErrors gCa48(nCa48,eCa48,sCa48,eeCa48,esCa48);
 gCa48.SetMarkerStyle(20);
 gCa48.SetMarkerColor(2);
 gCa48.SetLineColor(2);
 gCa48.SetMarkerSize(1.5);
 gCa48.SetLineWidth(10);
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

    
     if (energy < -10)
       {
        eeeCa48[nnCa48] = energy;
        ssCa48[nnCa48] = one;
        nnCa48++;
       }
   }

 TGraph ggCa48(nnCa48,eeeCa48,ssCa48);
 ggCa48.SetLineColor(2);
 ggCa48.Draw("L");

 text.SetTextColor(2);
 text.DrawLatex(.3,.47,"^{48}Ca");





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
     if (three > 0.)
       {
        eCa40[nCa40] = energy;
        sCa40[nCa40] = three;
        esCa40[nCa40] = four;
        eeCa40[nCa40] = 0.;
        nCa40++;
       }
   }

 

 TGraphErrors gCa40(nCa40,eCa40,sCa40,eeCa40,esCa40);
 gCa40.SetMarkerStyle(20);
 gCa40.SetMarkerColor(3);
 gCa40.SetLineColor(3);
 gCa40.SetMarkerSize(1.5);
 gCa40.SetLineWidth(10);
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
        ssCa40[nnCa40] = one;
        nnCa40++;
       }
   }

 TGraph ggCa40(nnCa40,eeeCa40,ssCa40);
 ggCa40.SetLineColor(3);
 ggCa40.Draw("L");

 text.SetTextColor(3);
 text.DrawLatex(.7,.42,"^{40}Ca");


}
