{
  gROOT->Reset();
  gPad->SetLogy();
  TH2S frame("frame","",10,-3,3,10,1e-4,.4);
    frame.GetXaxis()->SetTitle("r_{1} - r_{2} [fm]");
    frame.GetYaxis()->SetTitle("H");
    frame.GetXaxis()->CenterTitle();
    frame.GetYaxis()->CenterTitle();
    frame.SetStats(kFALSE);
    frame.SetTitle("non-locality function");
    frame.Draw();

    TF1 funct("funct","exp(-pow(x/.876,2))/pow(.876,3)/pow(3.14159,1.5)/1.379  +.379/1.379*exp(-pow(x/1.737,2))/pow(1.737,3)/pow(3.14159,1.5)",-3,3);
   funct.Draw("same");



  TF1 funct1("funct1","exp(-pow(x/.876,2))/pow(.876,3)/pow(3.14159,1.5)/1.379",-3,3);
  funct1.SetLineColor(2);
  funct1->SetLineWidth(1);
  funct1->SetLineStyle(9);
  funct1.Draw("same");


  TF1 funct2("funct2",".379/1.379*exp(-pow(x/1.737,2))/pow(1.737,3)/pow(3.14159,1.5)",-3,3);
  funct2.SetLineColor(4);
  funct2.SetLineStyle(9);
  funct2.SetLineWidth(1);
  funct2.Draw("same");
}
