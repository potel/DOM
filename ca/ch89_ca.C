{
  gROOT->Reset();
  TCanvas canvas("ch89_ca","",400,800);
  canvas.Divide(1,2);
  canvas.cd(1);
  TH2S frame ("frame","",10,0,13,10,-14,0);
  frame.GetXaxis()->SetTitle("r [fm]");
  frame.GetYaxis()->SetTitle("W [MeV]");
  frame.SetStats(kFALSE);
  frame.Draw();

  TF1 wp("wp","-.882/(1+exp((x-4.128)/.69))-4.*7.144/pow(1+exp((x-4.138)/.69),2)*exp((x-4.138)/.69)",0,13);
  wp.SetLineColor(4);
  wp.Draw("Same");



  TF1 wn("wn", "-1.35/(1+exp((x-4.128)/.69))-4.*5.00/pow(1+exp((x-4.138)/.69),2)*exp((x-4.138)/.69)",0,13);
  wn.SetLineColor(3);
  wn.Draw("SAME");

  double xn[100];
  double yn[100];
  ifstream fi("W10_n.dat");
  for (int i=0;i<97;i++)
    {
      fi >> xn[i] >> yn[i];
    }
  TGraph gn(97,xn,yn);
  gn.SetLineColor(3);
  gn.SetLineStyle(2);
  gn.SetLineWidth(2);
  gn.Draw("L");
  fi.close();
  fi.clear();


  double xp[100];
  double yp[100];
  ifstream  fii("W10_p.dat");
  for (int i=0;i<97;i++)
    {
      fii >> xp[i] >> yp[i];
      cout << xp[i] << " " << yp[i] << endl;
    }
  TGraph gp(97,xp,yp);
  gp.SetLineColor(4);
  gp.SetLineStyle(2);
  gp.SetLineWidth(2);
  gp.Draw("L");

  canvas.cd(2);
  TH2S frame2 ("frame2","",10,0,13,10,-60,1);
  frame2.GetXaxis()->SetTitle("r [fm]");
  frame2.GetYaxis()->SetTitle("V [MeV]");
  frame2.SetStats(kFALSE);
  frame2.Draw();
  TF1 rn("rn","-49.9/(1+exp((x-4.049)/.69))",0,13);
  rn.SetLineColor(3);
  rn.Draw("same");

  TF1 rp("rp","-52.28/(1+exp((x-4.049)/.69))",0,13);
  rp.SetLineColor(4);
  rp.Draw("same");

  fi.close();
  fi.clear();
  fi.open("R10_p.dat");
  double rxp[100];
  double ryp[100];
  for (int i=0;i<97;i++)
    {
      fi >> rxp[i] >> ryp[i];
    }
  TGraph rgp(97,rxp,ryp);
  rgp.SetLineColor(4);
  rgp.SetLineStyle(2);
  rgp.SetLineWidth(2);
  rgp.Draw("L");

  fi.close();
  fi.clear();
  fi.open("R10_n.dat");
  double rxn[100];
  double ryn[100];
  for (int i=0;i<97;i++)
    {
      fi >> rxn[i] >> ryn[i];
    }
  TGraph rgn(97,rxn,ryn);
  rgn.SetLineColor(3);
  rgn.SetLineStyle(2);
  rgn.SetLineWidth(2);
  rgn.Draw("L");

  TLatex text;
  text.SetTextColor(3);
  text.DrawLatex(2,-39,"n");
  text.SetTextColor(4);
  text.DrawLatex(2,-55,"p");

  canvas.cd(1);
  TLine line;
  line.SetLineWidth(2);
  line.DrawLine(.5,-10,2.5,-10);
  line.SetLineStyle(2);
  line.DrawLine(.5,-12,2.5,-12);
  text.SetTextColor(1);
  text.DrawLatex(3.,-10.2,"CH89");
  text.DrawLatex(3.,-12.2,"DOM");
  canvas.cd(2);
  text.SetTextSize(.08);
  text.DrawLatex(7,-10,"^{40}Ca");
  text.SetTextSize(.06);
  text.DrawLatex(7,-20,"E=10 MeV");
}
