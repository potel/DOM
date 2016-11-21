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
Sty->SetPadBottomMargin(.2);
Sty->SetPadLeftMargin(.15);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.07,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
Sty->SetTitleOffset(1.,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.05,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("ComnCa40");
 string out;

 gPad->SetLogy();
 
 TH2S frame("frame","",10,0,180,10,6,7500);
 frame.GetXaxis()->SetTitle("#theta_{cm} [deg]");
 frame.GetXaxis()->CenterTitle();
 frame.GetYaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");
 frame.GetYaxis()->CenterTitle();
 frame.Draw();



 ifstream fileOld("nca40.data");

 char name[80];

 double E,xx,yy,zz;
 int ii,Nold;
 for (int j=0;j<2;j++)
   {
 fileOld >> E >> ii;
 cout << E << " " << ii << endl;
 fileOld >> out ;
 cout << "out = " << out << endl;
 fileOld >> Nold >> out;
 cout << Nold << " " << out << endl;
 for (int i=0;i<Nold;i++)
   {
     fileOld >> xx >> yy >> zz; 
   }
 fileOld >> Nold >> out;
 for (int i=0;i<Nold;i++)
   {
     fileOld >> xx >> yy >> zz; 
   }
   }


 double x12Old[30];
 double y12Old[30];
 double sy12Old[30];
 double sx12Old[30];
 int N12old;

 fileOld >> E >> ii;
 cout << "E12= " << E << endl;
 fileOld >> out;
 cout << out  << endl;
 fileOld >> N12old >> out;
 cout << N12old << out << endl;
 for (int i=0;i<N12old;i++)
   {
     fileOld >> x12Old[i] >> y12Old[i] >> sy12Old[i]; 
     sy12Old[i] *= -y12Old[i]/100.;
     sx12Old[i] = 0.;
   }
 

 TGraphErrors g12Old(N12old,x12Old,y12Old,sx12Old,sy12Old);
 g12Old.SetMarkerStyle(25);
 g12Old.SetMarkerColor(4);
 g12Old.SetLineColor(4);
 g12Old.Draw("P");

 fileOld >> Nold >> out;
 cout << Nold << " " << out << endl;
 for (int i=0;i<Nold;i++)
   {
     fileOld >> xx >> yy >> zz; 
   }




 for (int j=0;j<2;j++)
   {
 fileOld >> E >> ii;
 cout << E << " " << ii << endl;
 fileOld >> out ;
 cout << "out = " << out << endl;
 fileOld >> Nold >> out;
 cout << Nold << " " << out << endl;
 for (int i=0;i<Nold;i++)
   {
     fileOld >> xx >> yy >> zz; 
   }
 fileOld >> Nold >> out;
 for (int i=0;i<Nold;i++)
   {
     fileOld >> xx >> yy >> zz; 
   }
   }

 double x17Old[30];
 double y17Old[30];
 double sy17Old[30];
 double sx17Old[30];
 int N17old;

 fileOld >> E >> ii;
 cout << "E= " << E << endl;
 fileOld >> out;
 fileOld >> N17old >> out;
 for (int i=0;i<N17old;i++)
   {
     fileOld >> x17Old[i] >> y17Old[i] >> sy17Old[i]; 
     y17Old[i] *= 60.;
     sy17Old[i] *= -y17Old[i]/100.;
     sx17Old[i] = 0.;
   }

 TGraphErrors g17Old(N17old,x17Old,y17Old,sx17Old,sy17Old);
 g17Old.SetMarkerStyle(26);
 g17Old.SetMarkerColor(4);
 g17Old.SetLineColor(4);
 g17Old.Draw("P");


 TLatex text;
 text.SetNDC();
 text.DrawLatex(.65,.63,"16.9 MeV(x60)");
 text.DrawLatex(.7,.35,"12. MeV");

 ifstream fileUS("nca40_FIN_17.data");

 getline(fileUS,out);
 getline(fileUS,out);
 int N17US = 0;
 fileUS >> N17US >> out;
 
 double x17US[10];
 double y17US[10];
 double sy17US[10];
 double sx17US[10];

 for (int i=0;i< N17US;i++)
   {
     fileUS >> x17US[i] >> y17US[i] >> sy17US[i];
     y17US[i] *= 60.;
     sx17US[i] = 0.;
     sy17US[i] *= -y17US[i]/100.;
   }

 TGraphErrors g17US(N17US,x17US,y17US,sx17US,sy17US);
 g17US.SetMarkerStyle(22);
 g17US.SetMarkerColor(2);
 g17US.SetLineColor(2);
 g17US.Draw("P");

 fileUS.close();
 fileUS.clear();



 fileUS.open("nca40_FIN_12.data");

 getline(fileUS,out);
 getline(fileUS,out);
 int N12US = 0;
 fileUS >> N12US >> out;
 
 double x12US[15];
 double y12US[15];
 double sy12US[15];
 double sx12US[15];

 for (int i=0;i< N12US;i++)
   {
     fileUS >> x12US[i] >> y12US[i] >> sy12US[i];
     //y12US[i] *= 40.;
     sx12US[i] = 0.;
     sy12US[i] *= -y12US[i]/100.;
   }

 TGraphErrors g12US(N12US,x12US,y12US,sx12US,sy12US);
 g12US.SetMarkerStyle(21);
 g12US.SetMarkerColor(2);
 g12US.SetLineColor(2);
 g12US.Draw("P");

 fileUS.close();
 fileUS.clear();

}
