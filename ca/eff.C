{
  gROOT->Reset();
  gROOT->SetStyle("Pub");
  TCanvas canvas("eff","",400,800);
  canvas.Divide(1,2);
  canvas.cd(1);

  double x[200];
  double y[200];
  double z[200];
  double z1[200];
  double z2[200];
  for (int i = 0;i<200;i++)
    {
      double k = ((double)i+.5)/20.;
      double eff = 1./(1. - 1./20.899/pow(k,2)*
		       (-69.829*(1.-exp(-pow(k*.876/2.,2)))+-26.465*(1.-exp(-pow(k*1.737/2.,2)))));
      double E = 20.899*pow(k,2)/eff - 69.82 - 26.465;
      x[i] = E;
      y[i] = eff;
      z[i] = 69.829*exp(-pow(k*.876/2.,2)) + 26.465*exp(-pow(k*1.737/2.,2));
      z1[i] = 69.829*exp(-pow(k*.876/2.,2));
      z2[i] = 26.465*exp(-pow(k*1.737/2.,2));
    }


  double xx[200];
  double yy[200];
  double zz[200];
  for (int i = 0;i<200;i++)
    {
      double k = ((double)i+.5)/20.;
      double eff = 1./(1. - 1./20.899/pow(k,2)*
		       (-48.28*(1.-exp(-pow(k*.746/2.,2)))));
      double E = 20.899*pow(k,2)/eff - 48.28;
      xx[i] = E;
      yy[i] = eff;
      zz[i] = 48.28*exp(-pow(k*.746/2.,2));
    }


  TH2S frame("frame","",10,-100,200,10,0,1.);
  frame.GetXaxis()->SetTitle("E [MeV]");
  frame.GetYaxis()->SetTitle("m_{eff}/m");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  //frame.GetYaxis()->SetTitleOffset(1.);
  frame.SetStats(kFALSE);
  frame.Draw();
  TGraph graph(200,x,y);
  graph.SetLineWidth(3);
  graph->Draw("L");


  TGraph ggraph(200,xx,yy);
  ggraph.SetLineColor(3);
  ggraph.SetLineStyle(9);
  //ggraph->Draw("L");


  TLine line;
  line.SetLineStyle(2);
  line.DrawLine(0.,0.,0.,1.);

  canvas.cd(2);

  TH2S frame1("frame1","",10,-100,200,10,0,100);
  frame1.GetXaxis()->SetTitle("E [MeV]");
  frame1.GetYaxis()->SetTitle("V_{local}(E)");
  frame1.GetXaxis()->CenterTitle();
  frame1.GetYaxis()->CenterTitle();
  frame1.GetYaxis()->SetTitleOffset(1.);

  frame1.SetStats(kFALSE);
  frame1.Draw();
  TGraph graph1(200,x,z);
  graph1->SetLineWidth(3);
  graph1->Draw("L");

  TGraph graph11(200,x,z1);
  graph11->SetLineColor(2);
  graph11->SetLineStyle(9);
  graph11->Draw("L");

  TGraph graph12(200,x,z2);
  graph12->SetLineColor(4);
  graph12->SetLineStyle(9);
  graph12->Draw("L");




  TGraph ggraph1(200,x,zz);
  ggraph1->SetLineColor(3);
  ggraph1->SetLineStyle(9);
  //ggraph1->Draw("L");


  TArrow arrow;
  arrow.SetFillColor(1);
  arrow.SetAngle(30);
  arrow.DrawArrow(0.,90,-50,90,.03,"<|>");
  arrow.DrawArrow(20.,55,200,55,.03,"<|>");
  TLatex text;
  text.DrawLatex(10,90,"single-particle states");
  text.DrawLatex(50,60,"scattering data");
  text.SetTextColor(2);
  text.DrawLatex(-50,48,"#beta_{R1}");
  text.SetTextColor(4);
  text.DrawLatex(-50,18,"#beta_{R2}");
}
