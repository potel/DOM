{
  gROOT->Reset();
  gROOT->SetStyle("Pub");
  TCanvas canvas("eff0");

  double x[200];
  double y[200];
  double z[200];
  double z1[200];
  double z2[200];
  double V0 = -50.;
  double beta = 1.;

  for (int i = 0;i<200;i++)
    {
      double k = ((double)i+.5)/20.;
      double E = 20.926*pow(k,2) + V0*exp(-pow(k*beta,2)/4.);
      double meff = 20.926*pow(k,2)/(E-V0); 
      double Vl = V0*exp(-pow(k*beta,2)/4.);


      x[i] = E;
      y[i] = meff;
      z[i] = Vl;
      cout << x[i] << " " << y[i] << " " << z[i] << endl;

    }


  TH2S frame("frame","",10,-50,200,10,0,1.);
  frame.GetXaxis()->SetTitle("E [MeV]");
  frame.GetYaxis()->SetTitle("m_{eff}/m");
  frame.SetStats(kFALSE);
  frame.Draw();
  TGraph graph(200,x,y);
  graph.SetLineWidth(3);
  graph->Draw("L");

  //return;
  TH2S frame("frame","",10,-50,200,10,-50,0.);
  frame.GetXaxis()->SetTitle("E [MeV]");
  frame.GetYaxis()->SetTitle("V_{local}(E) [MeV]");
  frame.SetStats(kFALSE);
  frame.Draw();
  TGraph graph2(200,x,z);
  graph2.SetLineWidth(3);
  graph2->Draw("L");


}
