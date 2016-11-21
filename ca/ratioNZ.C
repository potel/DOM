{
  gROOT->Reset();
  ifstream fexp("nca40_48.data");
  string name;
  getline(fexp,name);
  getline(fexp,name);

  double ex[100];
  double rx[100];
  double sex[100];
  double srx[100];

  double one,two,three;
  int Nx = 0;
  for(;;)
    {
      fexp >> one >> two >> three;
      //cout << one << " " << two << endl;
      if (fexp.eof()) break;
      if (fexp.bad()) break;
      ex[Nx] = one;
      rx[Nx] = two;
      srx[Nx] = three;
      sex[Nx] = 0.;
      Nx++;
    } 

  TH2S frame("frame","",10,0,200,10,0,.1);
  frame.Draw();

  TGraphErrors gx(Nx,ex,rx,sex,srx);
  gx.SetMarkerStyle(20);
  gx.Draw("P");

  fexp.close();
  fexp.clear();

  double efh[100];
  double rfh[100];
  int Nfh = 0;
  fexp.open("ratioHigh.dat");
  for(;;)
    {
      fexp >> one >> two;

      if (fexp.eof()) break;
      if (fexp.bad()) break;
      efh[Nfh] = one;
      rfh[Nfh] = two;
      Nfh++;
      if (Nfh == 99) break;
    }
  TGraph gfh(Nfh,efh,rfh);
  gfh.SetLineColor(2);
  gfh.Draw("L");

  fexp.close();
  fexp.clear();

  double efl[100];
  double rfl[100];
  int Nfl = 0;
  fexp.open("ratioLow.dat");
  for(;;)
    {
      fexp >> one >> two;

      if (fexp.eof()) break;
      if (fexp.bad()) break;
      efl[Nfl] = one;
      rfl[Nfl] = two;
      Nfl++;
      if (Nfl == 99) break;
    }
  TGraph gfl(Nfl,efl,rfl);
  gfl.SetLineColor(4);
  gfl.Draw("L");


  fexp.close();
  fexp.clear();

  double efm[100];
  double rfm[100];
  int Nfm = 0;
  fexp.open("ratioMin.dat");
  for(;;)
    {
      fexp >> one >> two;

      if (fexp.eof()) break;
      if (fexp.bad()) break;
      efm[Nfm] = one;
      rfm[Nfm] = two;
      Nfm++;
      if (Nfm == 99) break;
    }
  TGraph gfm(Nfm,efm,rfm);
  gfm.SetLineColor(1);
  gfm.Draw("L");

}
