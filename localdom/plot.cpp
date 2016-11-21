#include "reaction.h"

//***********************************************************************
  /**
   *plots the fit to the elastic scattering angular distributions
   * and analyzing powers- output is wriiten to root file title.root
   */
void reaction::PlotFit()
{

#ifdef root
  double factA = 0;

  string nnn = "fit_"+title;
  TCanvas *fit = new TCanvas(nnn.c_str());
  fit->Divide(3,1);
  string hhh = nnn+"_1";
  TPad* fit_1 = (TPad*)(fit->GetPrimitive(hhh.c_str()));
  fit_1->SetLogy();
  hhh = nnn + "_2";
  TPad* fit_2 = (TPad*)(fit->GetPrimitive(hhh.c_str()));
  fit_2->SetLogy();

  fit->cd(1);

  // axies
  int ymin = (int)floor(log10(xsecMin)) - Ndata;
  int ymax = (int)ceil(log10(xsecMax));
  TH2S *hist = new TH2S("hist",title.c_str(),10,0,180.
                        ,10,pow(10.,ymin),pow(10.,ymax));
  hist->GetXaxis()->SetTitle("#theta_{c.m.} [deg]");
  hist->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");
  hist->SetStats(kFALSE);
  hist->Draw();


  TH2S *hist2 = new TH2S("hist2",title.c_str(),10,0,180.,10,.1,pow(8.,Ndata));
  hist2->GetXaxis()->SetTitle("#theta_{c.m.} [deg]");
  hist2->GetYaxis()->SetTitle("#sigma/#sigma_{Ruth}");
  hist2->SetStats(kFALSE);
  fit->cd(2);
  hist2->Draw();


  TH2S *hist3 = new TH2S("hist3",title.c_str(),10,0,180.,10,-1.,32.);
  hist3->GetXaxis()->SetTitle("#theta_{c.m.} [deg]");
  hist3->GetYaxis()->SetTitle("Analyze Power");
  hist3->SetStats(kFALSE);
  fit->cd(3);
  hist3->Draw();

  fit->cd(1);

  TGraphErrors *gdata[Ndata];
  TGraph *GA[Ndata];
  TGraph *GAT[Ndata];
  TLine  *LA[Ndata];
  double zero[100];
  for (int i=0;i<100;i++) zero[i] = 0.;

  string filename = title+"x.dat";
  ofstream fdata(filename.c_str());
  filename = title+"xx.dat";
  if (Zp == 1) ofstream fdatax(filename.c_str());

  int jj = 0;
  for (int i=0;i<Ndata;i++)
    {

      if (data[i].nX == 0) continue;
      fdata << data[i].energyLab << " " << data[i].nX <<  endl;

      double y[data[i].nX];
      double error[data[i].nX];
      double fact = pow(10.,jj);
      for (int j=0;j<data[i].nX;j++)
	{
	  y[j] = data[i].xsec[j]/fact;
          error[j] = data[i].Xsigma[j]/fact;
	  fdata << data[i].Xtheta[j] << " " << data[i].xsec[j] << " " <<
	    data[i].Xsigma[j] << endl;
	}
      gdata[i] = new TGraphErrors(data[i].nX,
                data[i].Xtheta,y,zero,error);
      gdata[i]->SetLineColor(i%2+1);
      gdata[i]->SetMarkerStyle(21);
      gdata[i]->SetMarkerSize(.2);
      gdata[i]->SetMarkerColor(i%2+1);
      gdata[i]->Draw("P");
      jj++;
    }
  fdata.close();
  fdata.clear();
  filename = title+"c.dat";
  fdata.open(filename.c_str());
  filename = title+"a.dat";
  ofstream hdata(filename.c_str());


  filename = title+"xx.dat";
  ofstream fxdata(filename.c_str());

  TGraph *Gdata[Ndata];
  TGraph *Gratio[Ndata];
  TGraphErrors *GratioExp[Ndata];
  TLine  *Lratio[Ndata];
  jj = 0;

  if (DOM)
    {
     load();
     AdjustFermi();
    }
  else loadOM();

  for(int i=0;i<Ndata;i++)
    {

      if (data[i].nX > 0)fdata << data[i].energyLab << " " << 177 << endl;
      double Ecm =data[i].energyCM;
      double Elab =data[i].energyLab;



      InitializeForEcm(Ecm,Elab);

       // integrate wavefunction and find phaseshifts
       scatter->integrateWave();
       if (Ecm < 15.)scatter->statistical(scatter->konst,Ecm);
       if (data[i].nX > 0)
	 {
          double ang[180];
          double yy[180];
          double rr[180];
          double fact = pow(10.,jj);
          double factr = pow(10.,jj);
          // loop over angles
	  fxdata << 177 << endl;
          for (int j=3;j<180;j++)
	    {
	      double angle = (double)j*pi/180.;
              yy[j-3] = scatter->DifferentialXsection(angle);
	      if (Ecm < 15.)
		{
                  double extra = scatter->DifferentialXsectionCE(angle);
		  yy[j-3] += extra;
		}
              ang[j-3] = (double)j;
	      double ratioR = yy[j-3]/scatter->Rutherford(angle);
              rr[j-3] = ratioR*factr;
              fdata << ang[j-3] << " " << yy[j-3] << endl;
              yy[j-3] /= fact;

              fxdata << angle << " " << rr[j-3]<< endl;
	    }
          double rrexp[data[i].nX];
          double rrexpSigma[data[i].nX];
          fxdata << data[i].nX << " " << Elab << endl;
          for (int j=0;j<data[i].nX;j++)
	    {
	      double angle = data[i].Xtheta[j]*pi/180.;
              double ruth = scatter->Rutherford(angle);
              rrexp[j] = data[i].xsec[j]/ruth*factr;
              rrexpSigma[j] = data[i].Xsigma[j]/ruth*factr;
              fxdata << angle << " " << rrexp[j] 
              << " " << rrexpSigma[j] << endl;
   	    }
          Lratio[i] = new TLine(0.,factr,180.,factr);
          Lratio[i] -> SetLineColor(i%2+1);

          Gdata[i] = new TGraph(177,ang,yy);
          Gdata[i]->SetLineColor(i%2+1);
          Gdata[i]->Draw("C");
          Gratio[i] = new TGraph(177,ang,rr);
          Gratio[i]->SetLineColor(i%2+1);
          GratioExp[i] = new TGraphErrors(data[i].nX,data[i].Xtheta,rrexp,
                              zero,rrexpSigma);
          GratioExp[i]->SetLineColor(i%2+1);
          GratioExp[i]->SetMarkerStyle(21);
          GratioExp[i]->SetMarkerSize(.5);
          GratioExp[i]->SetMarkerColor(i%2+1);
          fit->cd(2);
          Lratio[i]->Draw();
          GratioExp[i]->Draw("P");
          Gratio[i]->Draw("C");
          fit->cd(1);
	  jj++;
	 }

       //see if there is analyzing power data
       if (data[i].nA > 0)
	 {
           hdata << data[i].energyLab << " " << data[i].nA << endl;

           double aa[data[i].nA];
	   for (int j=0;j<data[i].nA;j++)
	     {
	       aa[j] = data[i].anal[j] + factA;
	       hdata << data[i].Atheta[j] << " " 
		     << data[i].anal[j] << " " << data[i].Asigma[j] << endl;
	     }
           LA[i] = new TLine(0.,factA,180.,factA);
           LA[i]->SetLineColor(i%2+1);
	   GA[i] = new TGraphErrors(data[i].nA,data[i].Atheta,aa,zero,
           data[i].Asigma);
           GA[i]->SetMarkerStyle(21);
           GA[i]->SetMarkerSize(.5);
           GA[i]->SetMarkerColor(i%2+1);
	   fit->cd(3);
           GA[i]->Draw("P");
           LA[i]->Draw();
	   double aaa[180];
	   double spinR[180];
           double ang[180];
           hdata << 177 << endl;
	   for (int j=3;j<180;j++)
	     {
	       double angle = (double)j*pi/180.;
	       double shape = scatter->DifferentialXsection(angle);


	       ang[j-3] = (double)j;
	       aaa[j-3] = scatter->AnalyzePower ;
           
               spinR[j-3] = scatter->SpinRotation;
               if (Ecm < 15.)
		 {
                  double compound = scatter->DifferentialXsectionCE(angle);
                  aaa[j-3] *= shape/(shape+compound);
                  spinR[j-3] *= shape/(shape+compound);
		 }
	       hdata << ang[j-3] << " " << aaa[j-3] << " " 
                     << spinR[j-3] << endl;
               aaa[j-3] += factA;
	     }
           GAT[i] = new TGraph(177,ang,aaa);
           GAT[i]->SetLineColor(i%2+1);
           GAT[i]->Draw("C");
           factA += 2.;
	   fit->cd(1);
	 }
    }
   fdata.close();
   fdata.clear();
   hdata.close();
   hdata.clear();
   fxdata.close();
   fxdata.clear();
   filename = title+"d.dat";
   fdata.open(filename.c_str());

  //plot compound elastic for lowest energy
   double Ecm =data[0].energyCM;
   double Elab =data[0].energyLab;

   fdata << Elab << " " << 177 << endl;
   InitializeForEcm(Ecm,Elab);

   // integrate wavefunction and find phaseshifts
   scatter->integrateWave();
   //get cross section
   double xabs = scatter->AbsorptionXsection();

   if (Ecm < 30.)
     {
     //calulate CN decay widths
     scatter->statistical(scatter->konst,Ecm);

     cout << "absorption + compound elastic = " << xabs << " mb" << endl;
     cout << "absorption = " << scatter->sigmaAbsorption 
        << " mb " << endl;
     cout << "compound elastic = " << scatter->sigmaCompoundElastic 
        << " mb " << endl;

     // loop over angles
     double yy[174];
     double zz[174];
     double ang[174];
     for (int j=3;j<180;j++)
        {
         double angle = (double)j*pi/180.;
         yy[j-3] = scatter->DifferentialXsectionCE(angle);
         zz[j-3] = scatter->DifferentialXsection(angle)+ yy[j-7];

         ang[j-3] = (double)j;
         double  aaa = scatter->AnalyzePower ;
         //cout << aaa << endl;
         aaa *= (zz[j-3]-yy[j-3])/zz[j-3];
         fdata << ang[j-3] << " " << yy[j-3] << endl;
         }
      TGraph CNelastic (177,ang,yy);
      CNelastic.SetLineColor(3);
      CNelastic.SetLineStyle(2);
      CNelastic.Draw("C");
      //total
      TGraph totElastic(177,ang,zz);
      totElastic.SetLineStyle(1);
      totElastic.SetLineColor(2);
      totElastic.Draw("C");
     }
   else
     {
       scatter->sigmaAbsorption = 0.;
       scatter->sigmaCompoundElastic = 0.;
     }

  fit->Write();
  fdata.close();
  fdata.clear();

#endif
}
//****************************************************************************
  /**
   *plots the fitted potential and many other graphs
   * output is wriiten in the root file title.root
   */
void reaction::PlotPotentialEcm()
{
#ifdef root
  cout << "pot" << endl;
  string nnn = "pot_"+title;
  TCanvas *pot = new TCanvas(nnn.c_str());

  double Emin = -200.;
  double Emax = 200.;
  // axies
  TH2S *hist4 = new TH2S("hist4",title.c_str(),10,Emin,Emax
                        ,10,-15.,40.);
  hist4->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
  hist4->GetYaxis()->SetTitle("Potential [MeV]");
  hist4->SetStats(kFALSE);
  hist4->Draw();
  int const Npoints =150;
  double ee[Npoints];
  double WWsurfaceLow[Npoints];
  double WWsurfaceHigh[Npoints];
  double WWvolume[Npoints];
  double VVHF[Npoints];
  double VVsurfaceLow[Npoints];
  double VVsurfaceHigh[Npoints];
  double VVvolume[Npoints];
  double VVHFsurface[Npoints];
  double Vtotal[Npoints];
  double DerVVsurfaceLow[Npoints];
  double VVsurfaceLowPH[Npoints];
  double DerVVsurfaceHigh[Npoints];
  double DerVVvolume[Npoints];

  double JReal[Npoints];
  double JImag[Npoints];
  double JSO[Npoints];
  double RrmsReal[Npoints];
  double RrmsImag[Npoints];
  double RrmsSO[Npoints];

  string file_name = title +".poten";
  ofstream pFile(file_name.c_str());
  for (int i=0;i<Npoints;i++)
    {
      double Elab = (Emax-Emin)/(double)Npoints*(double)i + Emin;
      double Ecm = energyLab2Cm(Elab);
      ee[i] = Ecm;
      WWvolume[i] = Volume.volume.ImaginaryPotential(Ecm);
      WWsurfaceLow[i] = SurfaceLow.CentralImaginaryPotential(Ecm);
      WWsurfaceHigh[i] = SurfaceHigh.CentralImaginaryPotential(Ecm)
	+ Volume.surface.funct(Ecm);


      HartreeFock.SetEnergy(Ecm);
      VVHF[i] = HartreeFock.RealVolume.V;
      VVHFsurface[i] = HartreeFock.RealSurface(Ecm);
      VVvolume[i] = Volume.volume.DeltaRealPotential(Ecm);


      VVsurfaceLow[i] =  SurfaceLow.CentralDeltaRealPotential(Ecm);
      VVsurfaceHigh[i] =  SurfaceHigh.CentralDeltaRealPotential(Ecm);
      VVsurfaceLowPH[i] = SurfaceLow.CentralParticleHole(Ecm);
      Vtotal[i] = VVHF[i] + VVvolume[i];
      DerVVsurfaceLow[i] =10.*SurfaceLow.CentralDerDeltaRealPotential(Ecm);
      VVsurfaceLowPH[i] =40.*SurfaceLow.CentralParticleHole(Ecm);
      DerVVsurfaceHigh[i] =10.*SurfaceHigh.CentralDerDeltaRealPotential(Ecm);


      DerVVvolume[i] = 10.*(Volume.volume.DerDeltaHole(Ecm) + 
	Volume.volume.DerDeltaParticle(Ecm));

      pFile << Ecm-Efermi << " " << WWsurfaceLow[i] << " " 
             << WWsurfaceHigh[i] <<   " " << WWvolume[i] << " " 
	    << VVHF[i] << " " << VVHFsurface[i] << " " << 
            VVsurfaceLow[i] << " " << VVsurfaceHigh[i] <<
         " " << VVvolume[i] << " " <<
	DerVVsurfaceLow[i]/10. << " " << DerVVvolume[i]/10. << " " <<
        Volume.volume.DerDeltaHole(Ecm) << " " << 
        Volume.volume.DerDeltaParticle(Ecm) << endl;

      VVHF[i] /= 6.;
      Vtotal[i] /= 6.;



      InitializeForEcm(Ecm,Elab);
      scatter->VIntegrals();

      JReal[i] = -scatter->JReal;
      JImag[i] = -scatter->JImag;
      JSO[i] = -scatter->JSO;
      RrmsReal[i] = scatter->RrmsReal;
      RrmsImag[i] = scatter->RrmsImag;
      RrmsSO[i] =  scatter->RrmsSO;

    }
  pFile.close();
  pFile.clear();

  TGraph *gWvolume = new TGraph(Npoints,ee,WWvolume);
  TGraph *gWsurfaceLow = new TGraph(Npoints,ee,WWsurfaceLow);
  TGraph *gWsurfaceHigh = new TGraph(Npoints,ee,WWsurfaceHigh);
  TGraph *gVvolume = new TGraph(Npoints,ee,VVvolume);
  TGraph *gVsurfaceLowPH = new TGraph(Npoints,ee,VVsurfaceLowPH); 
  TGraph *gVsurfaceLow = new TGraph(Npoints,ee,VVsurfaceLow);
  TGraph *gVsurfaceHigh = new TGraph(Npoints,ee,VVsurfaceHigh);
  TGraph *gVHF = new TGraph(Npoints,ee,VVHF);
  TGraph *gVHFsurface = new TGraph(Npoints,ee,VVHFsurface);
  TGraph *gVtotal = new TGraph(Npoints,ee,Vtotal);

  TGraph *gDerVsurfaceLow = new TGraph(Npoints,ee,DerVVsurfaceLow);
  TGraph *gDerVsurfaceHigh = new TGraph(Npoints,ee,DerVVsurfaceHigh);
  TGraph *gDerVvolume = new TGraph(Npoints,ee,DerVVvolume);

  


  gWvolume->SetLineColor(2);
  gWsurfaceLow->SetLineColor(6);
  gWsurfaceHigh->SetLineColor(7);
  gWvolume->SetLineWidth(3);
  gWsurfaceLow->SetLineWidth(3);
  gWsurfaceHigh->SetLineWidth(3);
  gVvolume->SetLineColor(2);
  gVsurfaceLow->SetLineColor(6);
  gVsurfaceLowPH->SetLineColor(1);
  gVsurfaceHigh->SetLineColor(7);
  gVHF->SetLineColor(4);
  gVHFsurface->SetLineColor(5);
  gVtotal->SetLineColor(1);
  gVtotal->SetLineWidth(3);
  


  gWvolume->Draw("C");
  gWsurfaceLow->Draw("C");
  gWsurfaceHigh->Draw("C");
  gVvolume->Draw("C");
  //gVisoscalerSurface->Draw("C");
  //gVisovectorSurface->Draw("C");
  gVsurfaceLow->Draw("C");
  gVsurfaceLowPH->Draw("C");
  gVsurfaceHigh->Draw("C");

  gVHF->Draw("C");
  gVHFsurface->Draw("C");
  gVtotal->Draw("C");

  gDerVsurfaceLow->SetLineColor(6);
  gDerVsurfaceHigh->SetLineColor(7);
  gDerVvolume->SetLineColor(2);
  gDerVsurfaceLow->SetLineStyle(2);
  gDerVsurfaceHigh->SetLineStyle(2);
  gDerVvolume->SetLineStyle(2);
  gDerVsurfaceLow->Draw("C");
  gDerVsurfaceHigh->Draw("C");
  gDerVvolume->Draw("C");
  pot->Write();

  //----------------------------------
  cout << "Vint" << endl;
  //now for Vintegrals


  string fname = directory + title+".Vint";

  ifstream FileVint(fname.c_str());

  // if (!FileVint) cout << "could not open file " << fname << endl;
  //else 
    {
      nnn = "Vint_"+title;
      TCanvas *Vint = new TCanvas(nnn.c_str());
     file_name = title+".Vinteg";
     pFile.open(file_name.c_str());
     string line;
     int Nint, flag;
     getline(FileVint,line);

     if(FileVint) FileVint >> Nint >> flag;
 
     double EE[Nint];
     double JVE[Nint];
     double JWE[Nint];
     double RVE[Nint];
     double RWE[Nint];
     double JSOE[Nint];
     double RSOE[Nint];
     string author;
 
     if (FileVint)
       {
        pFile << Nint << endl;
        for (int i=0;i<Nint;i++)
        { 
           FileVint >> EE[i] >> JVE[i] >> JWE[i] >> RVE[i] >> RWE[i] >> 
	     JSOE[i] >> RSOE[i];
           getline(FileVint,author);

          EE[i] *= A/(1.+A);
          pFile << EE[i]-Efermi << " " << JVE[i] << " " << JWE[i] << " " <<
	    JSOE[i] << " " << RVE[i] << " " << RWE[i] << " " << RSOE[i] << endl;
         }
       }
     pFile << Npoints << endl;
     for (int i=0;i<Npoints;i++) 



       pFile << ee[i]*(1.+A)/A << " " << JReal[i] << " " << JImag[i] << " " <<
	 JSO[i] << " " << RrmsReal[i] << " " << RrmsImag[i] << 
	 " " << RrmsSO[i] <<  endl;


    pFile.close();
    pFile.clear();
    FileVint.close();
    FileVint.clear();

    Vint->Divide(1,2);
    Vint->cd(1);

    TH2S histt("histt",title.c_str(),10,Emin,Emax
                        ,10,0,900.);
    histt.GetXaxis()->SetTitle("E [Mev]");
    histt.GetYaxis()->SetTitle("J [MeV]");
    histt.SetStats(kFALSE);
    histt.Draw();

    TGraph *Jr = new TGraph(Npoints,ee,JReal);
    Jr->Draw("C");
    TGraph *Ji = new TGraph(Npoints,ee,JImag);
    Ji->SetLineColor(2);
    Ji->Draw("C");
    TGraph *Jso = new TGraph(Npoints,ee,JSO);
    Jso->SetLineColor(3);
    Jso->Draw("C");

    //plot data from optical model fits
    if (Nint > 0)
      {
        TGraph *Jrexp = new TGraph(Nint,EE,JVE);
        Jrexp->SetMarkerStyle(21);
        Jrexp->Draw("P");

        TGraph *Jiexp = new TGraph(Nint,EE,JWE);
        Jiexp->SetMarkerStyle(21);
        Jiexp->SetMarkerColor(2);
        Jiexp->Draw("P");

        TGraph *Jsoexp = new TGraph(Nint,EE,JSOE);
        Jsoexp->SetMarkerStyle(22);
        Jsoexp->SetMarkerColor(3);
        Jsoexp->Draw("P");
      }

    Vint->cd(2);
    double rmin=3.6;
    double rmax = 6.;
    if (A > 100)
      {
	rmin += 2.;
        rmax += 2.;
      }
    else if (A < 20)
      {
        rmin -= 2.;
        rmax -= 2.;
      }
    TH2S histr("histr",title.c_str(),10,Emin,Emax
                        ,10,rmin,rmax);
    histr.GetXaxis()->SetTitle("E [MeV]");
    histr.GetYaxis()->SetTitle("Rrms [fm]");
    histr.SetStats(kFALSE);
    histr.Draw();

    TGraph *Rr = new TGraph(Npoints,ee,RrmsReal);
    Rr->Draw("C");
    TGraph *Ri = new TGraph(Npoints,ee,RrmsImag);
    Ri->SetLineColor(2);
    Ri->Draw("C");
    TGraph *Rso = new TGraph(Npoints,ee,RrmsSO);
    Rso->SetLineColor(3);
    Rso->Draw("C");
  
   //plot data from optical model fits
    if (Nint > 0)
      {
        TGraph *Rrexp = new TGraph(Nint,EE,RVE);
        Rrexp->SetMarkerStyle(21);
        Rrexp->Draw("P");  

        TGraph *Riexp = new TGraph(Nint,EE,RWE);
        Riexp->SetMarkerStyle(21);
        Riexp->SetMarkerColor(2);
        Riexp->Draw("P");

        TGraph *Rsoexp = new TGraph(Nint,EE,RSOE);
        Rsoexp->SetMarkerStyle(22);
        Rsoexp->SetMarkerColor(3);
        Rsoexp->Draw("P");
      }

    Vint->Write();
    }

  //----------------------------------
  cout << "eff" << endl;
  //now for effective mass
  nnn = "eff_"+title;
  TCanvas *eff = new TCanvas(nnn.c_str());
  int const Rpoints = 60;

   double effHF[Rpoints];
   double effstar[Rpoints];
   double effBar[Rpoints];
   double rarray[Rpoints];


   SurfaceLow.SetEnergy(Efermi);
   SurfaceHigh.SetEnergy(Efermi);
   Volume.SetEnergy(Efermi);

   HartreeFock.SetEnergy(Efermi);


   file_name = title+".eff";
   ofstream eFile(file_name.c_str());
   for (int i=0;i<Rpoints;i++)
     {
       double r = (double)i/4.;
       rarray[i] = r;
       effHF[i] = scatter->eMassHF(r);
       effstar[i] = scatter->eMass(r);
       effBar[i] = scatter->eMassBar(r);

       eFile << r << " " << effstar[i] << " " << effHF[i] << " " << effBar[i]
     << endl;
     }
  // axies
   eFile.close();
  TH2S *hist5 = new TH2S("hist5",title.c_str(),10,0.,15.
                        ,10,0.,2.);
  hist5->GetXaxis()->SetTitle("r [fm]");
  hist5->GetYaxis()->SetTitle("m*/m");
  hist5->SetStats(kFALSE);
  hist5->Draw();

   TGraph *geffstar = new TGraph(Rpoints,rarray,effstar);
   geffstar->SetLineWidth(2);
   geffstar->Draw("C");

   TGraph *geffHF = new TGraph(Rpoints,rarray,effHF);
   geffHF->SetLineWidth(2);
   geffHF->SetLineColor(2);
   geffHF->Draw("C");

   TGraph *geffBar = new TGraph(Rpoints,rarray,effBar);
   geffBar->SetLineWidth(2);
   geffBar->SetLineColor(3);
   geffBar->Draw("C");

  eff->Write();
  //------------------------------------------------------------
  cout << "wfunct" << endl;
  //plot wave functions
  nnn = "wfunct_"+title;
  TCanvas *wfunct = new TCanvas(nnn.c_str());
  wfunct->Divide(2,2);
  TH2S *hist6 = new TH2S("hist6",title.c_str(),10,0.,35.
                        ,10,-.6,.6);

  hist6->GetXaxis()->SetTitle("radius [fm]");
  hist6->GetYaxis()->SetTitle("reduced Wave Function");
  hist6->SetStats(kFALSE);
  TLine *line0 = new TLine(0.,0.,35.,0.);
  for (int i=1;i<=4;i++)
    {
    wfunct->cd(i);
    hist6->Draw();
    line0->Draw();
    }

  TGraph *graphW[15];
  double rrr[scatter->nWave];
  for (int i=0;i<scatter->nWave;i++) 
        rrr[i] = scatter->rStart + (double)i*scatter->deltaR;



  //***here

 ofstream fw ("p12.dat"); //rjc


  NlevelTh = 0;
  int ifine = 1;
  int nmax[8] = {2,2,2,1,1,0,0,0};
  char cwave[8]={'s','p','d','f','g','h','i','j'};
  double Elower= -60.;
  double Eupper = (float)scatter->LogDerMax;
  if (scatter->proton && Eupper > 10.) Eupper = 10.;
  for (int l=0;l<8;l++)
    {
      for (int ii=-1;ii<2;ii+=2)
	{
          double j = (double)l + (double)ii*0.5;
          if (j < 0.) continue;
          Elower = -80.;           

          for (int n=0;n<=nmax[l];n++)
	    {
             if (bound(Elower,Eupper,j,l,ifine)==1.) 
                {
		 ostringstream outstring;
                 outstring << n << cwave[l] <<"#frac{"<<(int)(2.*j)<<"}{2}";
		 string name = outstring.str();
                 cout << name << endl;
                 LevelProperties(Elower,j,l);

                 LevelTh[NlevelTh].j = j;
                 LevelTh[NlevelTh].l = l;
                 LevelTh[NlevelTh].N = n;
                 LevelTh[NlevelTh].color = l + 1;
                 LevelTh[NlevelTh].name = name;
                 LevelTh[NlevelTh].energy = Elower;
                 LevelTh[NlevelTh].Rrms = scatter->Rrms;
                 LevelTh[NlevelTh].SpectFactor = scatter->SpectFactor;
                 LevelTh[NlevelTh].Width = scatter->Width;
                 LevelTh[NlevelTh].ANC = scatter->ANC;
                 LevelTh[NlevelTh].Occupation = scatter->Occupation;
                 graphW[NlevelTh] = new TGraph(scatter->nWave,rrr,scatter->WaveBar);
                 graphW[NlevelTh]->SetLineColor(LevelTh[NlevelTh].color);
                 wfunct->cd(LevelTh[NlevelTh].l+1);
                 graphW[NlevelTh]->SetLineStyle(n+1);
                 graphW[NlevelTh]->Draw("C");
                 spectralFunction(Elower,Nsf,Esf,LevelTh[NlevelTh].SpectFunction);
		 Elower += 0.5;
                 NlevelTh++;


		 if (n == 0 && l == 1 && j == 0.5) //rjc
		   {
                     for (int i=0;i<scatter->nWave;i++)
                        {
			  fw << rrr[i] << " " << scatter->WaveBar[i] << endl;
                        }
		   }

		}
	    
	     else break;
	    }
	}
    }
 
  fw.close();  //rjc
  fw.clear();  //rjc

  int jj = 0;
  // sort levels in energy order
  for (int j=0; j<NlevelTh;j++)
    {
      double Emin = 10000000.;
      int imin=1000000;
      for (int i=0;i<NlevelTh;i++)
        {
	  if (LevelTh[i].energy < Emin)
	  {
	    Emin = LevelTh[i].energy;
            imin = i;
	  }
        }
      LevelThSort[jj] = LevelTh[imin];
      LevelTh[imin].energy = 10000000000.;
      jj++;
    }

  
  wfunct->Write();

  



  

  //-------------------------------------------------------
  cout << "lev" << endl;
  //now for the levels
  double const TheoryLeft = 1.1;
  double const TheoryRight = 1.7;
  double const ExpLeft = .3;
  double const ExpRight = .9;
  nnn = "lev_"+title;
  TCanvas *lev = new TCanvas(nnn.c_str());
  TH2S *hist7 = new TH2S("hist7",title.c_str(),10,0.,2.
                        ,10,-70,5);

  hist7->GetYaxis()->SetTitle("Level Energy [MeV]");
  hist7->SetStats(kFALSE);
  hist7->Draw();
  TLatex tex;
  tex.SetTextAlign(12);
  tex.SetTextSize(.02);

  //draw fermi energy
  TLine *Lfermi = new TLine(ExpLeft,Efermi,TheoryRight,Efermi);
  Lfermi->SetLineWidth(3);
  Lfermi->SetLineStyle(2);
  Lfermi->Draw();

  TLine *Lf[NlevelTh];
    for (int i=0;i<NlevelTh;i++)
      {
      Lf[i] = new TLine(TheoryLeft,LevelThSort[i].energy,
                          TheoryRight,LevelThSort[i].energy);
      Lf[i]->SetLineColor(LevelThSort[i].color);
      Lf[i]->SetLineWidth(3);
      Lf[i]->Draw();
      tex.SetTextColor(LevelThSort[i].color);
      tex.DrawLatex(TheoryRight,LevelThSort[i].energy,
                          LevelThSort[i].name.c_str());
      }


  //plot experimental levels
  TLine line;
  line.SetLineWidth(3);
  if (Nlevel > 0) 
     {

     TLine *expLevel[Nlevel];
     for (int i=0;i<Nlevel;i++)
       {
         line.SetLineColor(Level[i].color);
         expLevel[i] = new TLine(ExpLeft,Level[i].Energy,
                                 ExpRight,Level[i].Energy);
         expLevel[i]->SetLineWidth(3);
         expLevel[i]->SetLineColor(Level[i].color);
         expLevel[i]->Draw();

         //connect to calculated value
         for (int j=0;j<NlevelTh;j++)
	   {
	     if (Level[i].l == LevelThSort[j].l &&
                 Level[i].j == LevelThSort[j].j &&
                 Level[i].N == LevelThSort[j].N)
	         {
		   line.DrawLine(ExpRight,Level[i].Energy,
		   			 TheoryLeft,LevelThSort[j].energy);
	         }
	   } 

       }
     }

  //write out in file level stuff
  string filename = title + ".level";
  ofstream LevelStuff(filename.c_str());
  LevelStuff << "N J  L       E     Rrms    occup   SpectF    Width color" << endl; 
    for (int i=0;i<NlevelTh;i++)
      {
	LevelStuff << LevelThSort[i].N << " " <<  LevelThSort[i].j << 
                   " " << LevelThSort[i].l 
		   << " "  << LevelThSort[i].energy << " " <<
	  LevelThSort[i].Rrms << " " << LevelThSort[i].Occupation << " " <<
	  LevelThSort[i].SpectFactor << " " << LevelThSort[i].Width << 
	  " "  << LevelThSort[i].ANC << " "   << LevelThSort[i].color << endl;
      }
  LevelStuff.close();
  LevelStuff.clear();


  lev->Write();
  //------------------------------------------------------------------------
  cout << "properties" << endl;
  nnn = "properties_"+title;
  TCanvas *prop = new TCanvas(nnn.c_str());
  prop->Divide(2,2);
  prop->cd(1);

  //Occupancy
  TH2S *hist8 = new TH2S("hist8",title.c_str(),10,-70.,5.,10,0.,1.);

  hist8->GetXaxis()->SetTitle("Level Energy [MeV]");
  hist8->GetYaxis()->SetTitle("Occupation prob");
  hist8->SetStats(kFALSE);
  hist8->Draw();

  double x[NlevelTh];
  double y[NlevelTh];
  for (int i=0;i<NlevelTh;i++)
    {
      x[i] = LevelThSort[i].energy;
      y[i] = LevelThSort[i].Occupation;
    }
  TGraph *graphOccup = new TGraph(NlevelTh,x,y);
  graphOccup->SetMarkerStyle(21);
  graphOccup->SetMarkerColor(2);
  graphOccup->SetMarkerSize(.5);
  graphOccup->Draw("PL");

  TLine line99(Efermi,0,Efermi,1);
  line99.SetLineStyle(2);
  line99.Draw();


  //Spectroscopic Factor
  prop->cd(2);
  TH2S *hist9 = new TH2S("hist9",title.c_str(),10,-70.,5.,10,0.,100.); 
  hist9->GetXaxis()->SetTitle("Level Energy [MeV]");
  hist9->SetStats(kFALSE);
  hist9->GetYaxis()->SetTitle("Spectroscopic Factor [%]");
  hist9->Draw();
  for (int i=0;i<NlevelTh;i++)y[i] = LevelThSort[i].SpectFactor*100.;
  TGraph *graphSpect = new TGraph(NlevelTh,x,y);
  graphSpect->SetMarkerStyle(21);
  graphSpect->SetMarkerColor(3);
  graphSpect->SetMarkerSize(.5);
  graphSpect->Draw("PL");

  int j=0;
  double xx[Nlevel];
  double yy[Nlevel];
  double dx[Nlevel];
  double dy[Nlevel];
  //experimental
  for (int i=0;i<Nlevel;i++)
    {
      if (Level[i].SigmaSpect != 0.)
	{
	  xx[j] = Level[i].Energy;
          yy[j] = Level[i].SpectFactor*100.;
          dx[j] = Level[i].SigmaEnergy;
          dy[j] = abs(Level[i].SigmaSpect)*100.;
	    j++;
	}
    }
  if (j>0)
    {
     TGraphErrors *expSpect = new TGraphErrors(j,xx,yy,dx,dy);
     expSpect->SetMarkerStyle(21);
     expSpect->SetMarkerSize(.75);
     expSpect->SetMarkerColor(4);
     expSpect->Draw("P"); 
    }
  //draw in energy range
  TLine Llow(Elow,0.,Elow,100.);
  Llow.SetLineStyle(2);
  Llow.Draw();
  TLine Lhigh(Ehigh,0.,Ehigh,100.);
  Lhigh.SetLineStyle(2);
  Lhigh.Draw();

  //line99.DrawLine(Efermi,.4,Efermi,1.2);
  //----------------------------------------------------------------
  //rms radius
  prop->cd(3);
  TH2S *hist10 = new TH2S("hist10",title.c_str(),10,-70.,5.,10,2.,7.);
  hist10->GetXaxis()->SetTitle("Level Energy [MeV]");
  hist10->SetStats(kFALSE);
  hist10->GetYaxis()->SetTitle("rms radius [fm]");
  hist10->Draw();
  for (int i=0;i<NlevelTh;i++)y[i] = LevelThSort[i].Rrms;
  TGraph *graphRms = new TGraph(NlevelTh,x,y);
  graphRms->SetMarkerStyle(21);
  graphRms->SetMarkerColor(4);
  graphRms->SetMarkerSize(.5);
  graphRms->Draw("PL");



  j=0;
  //experimental
  for (int i=0;i<Nlevel;i++)
    {
      if (Level[i].SigmaRrms > 0.)
	{
	  xx[j] = Level[i].Energy;
          yy[j] = Level[i].Rrms;
          dx[j] = Level[i].SigmaEnergy;
          dy[j] = Level[i].SigmaRrms;
	    j++;
	}
    }
  if (j > 0)
    {
    TGraphErrors *expRrms = new TGraphErrors(j,xx,yy,dx,dy);
    expRrms->SetMarkerStyle(21);
    expRrms->SetMarkerSize(.75);
    expRrms->SetMarkerColor(2);
    expRrms->Draw("P"); 
    }
  line99.DrawLine(Efermi,2.,Efermi,7.);
  //---------------------------------------------------------------
  //level width
  prop->cd(4);
  TH2S *hist11 = new TH2S("hist11",title.c_str(),10,-70.,5.,10,-.1,15.);
  hist11->GetXaxis()->SetTitle("Level Energy [MeV]");
  hist11->SetStats(kFALSE);
  hist11->GetYaxis()->SetTitle("level width [MeV]");
  hist11->Draw();
  for (int i=0;i<NlevelTh;i++)y[i] = LevelThSort[i].Width;
  TGraph *graphWidth = new TGraph(NlevelTh,x,y);
  graphWidth->SetMarkerStyle(21);
  graphWidth->SetMarkerColor(4);
  graphWidth->SetMarkerSize(.5);
  graphWidth->Draw("PL");

  j=0;
  //experimental
  for (int i=0;i<Nlevel;i++)
    {
      if (Level[i].SigmaDelta > 0.)
	{
	  xx[j] = Level[i].Energy;
          yy[j] = Level[i].Delta;
          dx[j] = Level[i].SigmaEnergy;
          dy[j] = Level[i].SigmaDelta;
	    j++;
	}
    }
  if (j > 0)
    {
    TGraphErrors *expDelta = new TGraphErrors(j,xx,yy,dx,dy);
     expDelta->SetMarkerStyle(21);
     expDelta->SetMarkerSize(.75);
     expDelta->SetMarkerColor(2);
     expDelta->Draw("P");
    }
  line99.DrawLine(Efermi,0.,Efermi,15.); 
  prop->Write();
  //-----------------------------------------------------------------------
  cout << "spectFunction" << endl;
  nnn = "spectFunct_"+title;
  TCanvas *specF = new TCanvas(nnn.c_str());
  TH2S *hist12 = new TH2S("hist12",title.c_str(),10,-80.,20.,10,0.,2.);
  hist12->GetXaxis()->SetTitle("E [MeV]");
  hist12->SetStats(kFALSE);
  hist12->GetYaxis()->SetTitle("spectral function [MeV-1]");
  hist12->Draw();

  /*
  LevelStuff.open("s12.dat");
    for (int i=0;i<Nsf;i++)
      {
       LevelStuff << Esf[i] << " " << LevelThSort[4].SpectFunction[i] << endl; 
      }
    LevelStuff.close();
    LevelStuff.clear();
  */

  double totSpect[Nsf];
  for (int i=0;i<Nsf;i++) totSpect[i] = 0.;
  TGraph *funct[NlevelTh];
  for (int i=0;i<NlevelTh;i++)
    {
      for (int j=0;j<Nsf;j++) 
	{
         LevelThSort[i].SpectFunction[j] *=2.*LevelThSort[i].j+1.;
	 totSpect[j] += LevelThSort[i].SpectFunction[j];
	}
      funct[i] = new TGraph(Nsf,Esf,LevelThSort[i].SpectFunction);
      funct[i]->SetLineColor(LevelThSort[i].color);
      funct[i]->Draw("L");
    }
  TGraph TotSp(Nsf,Esf,totSpect);
  TotSp.Draw("L");
  TLine la(Elow,0,Elow,2);
  la.SetLineStyle(2);
  la.Draw();
  la.DrawLine(Ehigh,0,Ehigh,2);
  specF->Write();
  //------------------  plot absorbtion j distribution
  cout << "absx" << endl;
  nnn = "absx_"+title;
  TCanvas *absx =new TCanvas(nnn.c_str());
  double jjj[scatter->lMax];

  for (int i=0;i<scatter->lMax;i++) jjj[i] = (float)i+0.5;
  TGraph *graphAbs[Ndata];
  TH2S histAbs("histabs",title.c_str(),10,0,20,10,0,300);
  histAbs.GetXaxis()->SetTitle("J");
  histAbs.GetYaxis()->SetTitle("#sigma_{abs} [mb]");
  histAbs.SetStats(kFALSE);
  histAbs.Draw();

  double comElast[scatter->Ntll];
  for(int i=0;i<Ndata;i++)
    {
      double Ecm =data[i].energyCM;
      double Elab =data[i].energyLab;

      InitializeForEcm(Ecm,Elab);

       // integrate wavefunction and find phaseshifts
       scatter->integrateWave();
       scatter->AbsorptionXsection();

       if (i== 0)
	 {
          scatter->statistical(scatter->konst,Ecm);

          for (int i=0;i<scatter->Ntll;i++) comElast[i] = scatter->xsec[i];
	 }

       graphAbs[i] = new TGraph(scatter->lStop+1,jjj,scatter->SigmaAb);
       int ic = i%3+1;
       graphAbs[i]->SetLineColor(ic);
       graphAbs[i]->SetMarkerStyle(21);
       graphAbs[i]->SetMarkerColor(ic);
       graphAbs[i]->Draw("PL");


  TGraph *gcomElast = new TGraph(scatter->Ntll,jjj,comElast);
  gcomElast->SetMarkerStyle(20);
  gcomElast->SetMarkerColor(4);
  gcomElast->SetLineColor(4);
  gcomElast->Draw("LP");


    }
  absx->Write();
  //-----------------------
  cout << "react" << endl;
  nnn = "react_"+title;
  TCanvas *Creact = new TCanvas(nnn.c_str());
  double ymax;
  double xmax;

  if (Zp > 0.) 
    {

      ymax = 1.;
      for (int i=0;i<NXdata;i++)
	{
	  if (Xdata[i].xsec > ymax) ymax = Xdata[i].xsec;
	}
      ymax *= 1.2;
     xmax = 200.;
    }
  else 
    {
      if (Z > 39) ymax = 10000.;
      else ymax = 4500.;
     xmax = 200.;
    }
  TH2S histXsec("histXsec",title.c_str(),10,0,xmax,10,0,ymax);

  histXsec.SetStats(kFALSE);
  histXsec.GetXaxis()->SetTitle("E_{lab} [MeV]");
  if (Zp > 0.) histXsec.GetYaxis()->SetTitle("#sigma_{react} [mb]");
  else histXsec.GetYaxis()->SetTitle("#sigma_{tot} [mb]");
  histXsec.Draw();


  double xxx[400];
  double yyy[400];
  double dxxx[400];
  double dyyy[400];
  double isovectorXsec[100];


  filename = title+"r.dat";
  ofstream fdata(filename.c_str());

  //experimental
  fdata << NXdata << endl;
  if (NXdata > 0)
    {

     for (int i=0;i<NXdata;i++)
       {
	 fdata << Xdata[i].energyLab << " " <<
	   Xdata[i].xsec << " " << Xdata[i].sigma <<endl;
         xxx[i] = Xdata[i].energyLab;
         dxxx[i] = 0.;
         yyy[i] = Xdata[i].xsec;
         dyyy[i] = Xdata[i].sigma;
 
       }


     TGraphErrors *xsec_exp = new TGraphErrors(NXdata,xxx,yyy,dxxx,dyyy);
     xsec_exp->SetMarkerStyle(21);
     xsec_exp->SetMarkerColor(2);
     xsec_exp->Draw("P");
    }

  fdata << 143 << endl;

  //calculated
  for(int i=0;i<100;i++)
    {
      double Elab =((double)i+.5)/5. + .5;
      xxx[i] = Elab;
     double Ecm = energyLab2Cm(Elab);

      InitializeForEcm(Ecm,Elab);

       // integrate wavefunction and find phaseshifts
       scatter->integrateWave();

       yyy[i] = scatter->AbsorptionXsection();


       dxxx[i] = yyy[i] - scatter->statistical(scatter->konst,Ecm);

       isovectorXsec[i] = scatter->isovectorAbsorb();
       /*
       cout << Elab << endl;
       cout << dxxx[i] << endl;
       cout <<  scatter->absorbAll << endl;
       cout <<	 scatter->absorbSurface << endl;
       cout << scatter->absorbSO << endl;
       cout <<   scatter->absorbVolume << endl;
       cout << endl;
       */

       fdata << Elab << " " << dxxx[i] << " " << scatter->absorbAll << " " <<
	 scatter->absorbSurface << " " << scatter->absorbSO << " " <<
         scatter->absorbVolume << endl;



    }


  TGraph xsec_th(100,xxx,yyy);
  xsec_th.Draw("C");
  xsec_th.SetLineColor(2);
  xsec_th.DrawGraph(100,xxx,dxxx);

  for (int i=0;i<43;i++)
    {
      double Elab = 20. + (double)i*180./32.;
      double Ecm = energyLab2Cm(Elab);

      xxx[i] = Elab;

      InitializeForEcm(Ecm,Elab);

      scatter->integrateWave();
      yyy[i] = scatter->AbsorptionXsection();

      isovectorXsec[i] = scatter->isovectorAbsorb();



       fdata << Elab << " " << yyy[i] << " " << scatter->absorbAll << " " <<
	 scatter->absorbSurface << " " << scatter->absorbSO << " " <<
         scatter->absorbVolume <<  endl;



    }




  TGraph *gxsec2 = new TGraph(33,xxx,yyy);
  gxsec2->SetLineColor(2);
  gxsec2->Draw("C");

  cout << "totxsec" << endl;
  //total cross sections for neutrons
  if (Zp == 0. && NTotXdata > 0)
    {
      fdata << NTotXdata << endl;
     for (int i=0;i<NTotXdata;i++)
       {
	 fdata << TotXdata[i].energyLab << " " <<
	   TotXdata[i].xsec << " " << TotXdata[i].sigma <<endl;
         xxx[i] = TotXdata[i].energyLab;
         dxxx[i] = 0.;
         yyy[i] = TotXdata[i].xsec;
         dyyy[i] = TotXdata[i].sigma;
 
       }


     TGraphErrors *txsec_exp = new TGraphErrors(NTotXdata,xxx,yyy,dxxx,dyyy);
     txsec_exp->SetMarkerStyle(22);
     txsec_exp->SetMarkerColor(4);
     txsec_exp->Draw("P");
    }

  cout << " ll " << endl;

  if (Zp == 0.)
    {
     //calculated
      fdata << 100 << endl;
     for(int i=0;i<100;i++)
       {
         double Elab =(double)(i+1)*2.;

         double Ecm = energyLab2Cm(Elab);
         xxx[i] = Elab;

         InitializeForEcm(Ecm,Elab);

         // integrate wavefunction and find phaseshifts
         scatter->integrateWave();

         yyy[i] = scatter->TotXsection();

         fdata << Elab << " " << yyy[i] << endl;
         //cout << Elab << " " << yyy[i] << endl;
       }


    TGraph *txsec_th = new TGraph(100,xxx,yyy);
    txsec_th->SetLineColor(4);
    txsec_th->Draw("C");
    }

  fdata.close();
  fdata.clear();



  Creact->Write();

  //*******************************************************************
    cout << "geoImag" << endl;
  nnn="geoW_"+title;
  TCanvas *CgeoW = new TCanvas(nnn.c_str());
    CgeoW->Divide(2,3);

    TH2S Hgeo("R", " " , 10,0,14,10,-20,0);
    Hgeo.GetXaxis()->SetTitle("r [fm]");
    Hgeo.GetYaxis()->SetTitle("Imaginary Potential [MeV]");
    Hgeo.SetStats(kFALSE);
    CgeoW->cd(1);
    Hgeo.Draw();
    double Elab = 10.;
    int const intR= 98;
    double volR10[intR];
    double surR10[intR];
    double totR10[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR10[i] = Volume.ImaginaryPot(xxx[i])
                  + SurfaceHigh.ImaginaryPot(xxx[i]);
        surR10[i] = SurfaceLow.ImaginaryPot(xxx[i]);

	totR10[i] = surR10[i] + volR10[i];
      }
    TGraph *geo10 = new TGraph(intR,xxx,totR10);
    geo10->Draw("C");
    TGraph *geov10 = new TGraph(intR,xxx,volR10);
    geov10->SetLineColor(2);
    geov10->Draw("C");

    CgeoW->cd(2);
    Hgeo.Draw();
    Elab = 30.;
    double volR30[intR];
    double surR30[intR];
    double totR30[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR30[i] = Volume.ImaginaryPot(xxx[i])
                  + SurfaceHigh.ImaginaryPot(xxx[i]);
        surR30[i] = SurfaceLow.ImaginaryPot(xxx[i]);

	totR30[i] = surR30[i] + volR30[i];
      }
    TGraph *geo30 = new TGraph(intR,xxx,totR30);
    geo30->Draw("C");
    TGraph *geov30 = new TGraph(intR,xxx,volR30);
    geov30->SetLineColor(2);
    geov30->Draw("C");

    CgeoW->cd(3);
    Hgeo.Draw();
    Elab = 50.;
    double volR50[intR];
    double surR50[intR];
    double totR50[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR50[i] = Volume.ImaginaryPot(xxx[i])
	          + SurfaceHigh.ImaginaryPot(xxx[i]);
        surR50[i] = SurfaceLow.ImaginaryPot(xxx[i]);

	totR50[i] = surR50[i] + volR50[i];
      }
    TGraph *geo50 = new TGraph(intR,xxx,totR50);
    geo50->Draw("C");
    TGraph *geov50 = new TGraph(intR,xxx,volR50);
    geov50->SetLineColor(2);
    geov50->Draw("C");

  CgeoW->cd(4);
  Hgeo.Draw();
    Elab = 75.;
    double volR75[intR];
    double surR75[intR];
    double totR75[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR75[i] = Volume.ImaginaryPot(xxx[i])  
                  + SurfaceHigh.ImaginaryPot(xxx[i]);
        surR75[i] = SurfaceLow.ImaginaryPot(xxx[i]);

	totR75[i] = surR75[i] + volR75[i];
      }
    TGraph *geo75 = new TGraph(intR,xxx,totR75);
    geo75->Draw("C");
    TGraph *geov75 = new TGraph(intR,xxx,volR75);
    geov75->SetLineColor(2);
    geov75->Draw("C");


   CgeoW->cd(5);
   Hgeo.Draw();
    Elab = 100.;
    double volR100[intR];
    double surR100[intR];
    double totR100[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR100[i] = Volume.ImaginaryPot(xxx[i])
                   + SurfaceHigh.ImaginaryPot(xxx[i]);
        surR100[i] = SurfaceLow.ImaginaryPot(xxx[i]);

	totR100[i] = surR100[i] + volR100[i];
      }
    TGraph *geo100 = new TGraph(intR,xxx,totR100);
    geo100->Draw("C");
    TGraph *geov100 = new TGraph(intR,xxx,volR100);
    geov100->SetLineColor(2);
    geov100->Draw("C");


    CgeoW->cd(6);
    Hgeo.Draw();
    Elab = 200.;
    double volR200[intR];
    double surR200[intR];
    double totR200[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR200[i] = Volume.ImaginaryPot(xxx[i])
                   + SurfaceHigh.ImaginaryPot(xxx[i]);
        surR200[i] = SurfaceLow.ImaginaryPot(xxx[i]);

	totR200[i] = surR200[i] + volR200[i];
      }
    TGraph *geo200 = new TGraph(intR,xxx,totR200);
    geo200->Draw("C");
    TGraph *geov200 = new TGraph(intR,xxx,volR200);
    geov200->SetLineColor(2);
    geov200->Draw("C");

    CgeoW->Write();

    //***********************************************************
    cout << "geoReal" << endl;
    nnn = "geoReal_"+title;
    TCanvas *CgeoR = new TCanvas(nnn.c_str());
    CgeoR->Divide(2,3);
    TLine geoline(0.,0.,8.,0.);
    geoline.SetLineStyle(2);

    TH2S HgeoR("RR", " " , 10,0,14,10,-15,5);
    HgeoR.GetXaxis()->SetTitle("r [fm]");
    HgeoR.GetYaxis()->SetTitle("Real Potential [MeV]");
    HgeoR.SetStats(kFALSE);
    CgeoR->cd(1);
    HgeoR.Draw();
    Elab = 10.;
    double HFR10[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR10[i] = Volume.DispersiveCorrection(xxx[i]);  
        surR10[i] = SurfaceLow.DispersiveCorrection(xxx[i])
                  + SurfaceHigh.DispersiveCorrection(xxx[i]);
        HFR10[i] = HartreeFock.RealPotential(xxx[i]);
	totR10[i] = surR10[i] + volR10[i] + HFR10[i];
        
      }
    TGraph *Rgeo10 = new TGraph(intR,xxx,totR10);
    Rgeo10->Draw("C");
    TGraph *Rgeov10 = new TGraph(intR,xxx,volR10);
    Rgeov10->SetLineColor(2);
    Rgeov10->Draw("C");
    TGraph *Rgeohf10 = new TGraph(intR,xxx,HFR10);
    Rgeohf10->SetLineColor(3);
    Rgeohf10->Draw("C");
    geoline.Draw();
    TGraph *Rgeos10 = new TGraph(intR,xxx,surR10);
    Rgeos10->SetLineColor(4);
    Rgeos10->Draw("C");

    CgeoR->cd(2);
    HgeoR.Draw();
    Elab = 30.;
    double HFR30[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR30[i] = Volume.DispersiveCorrection(xxx[i])
                  + SurfaceHigh.DispersiveCorrection(xxx[i]);
        surR30[i] = SurfaceLow.DispersiveCorrection(xxx[i]);
        HFR30[i] = HartreeFock.RealPotential(xxx[i]);
	totR30[i] = surR30[i] + volR30[i] + HFR30[i];

      }
    TGraph *Rgeo30 = new TGraph(intR,xxx,totR30);
    Rgeo30->Draw("C");
    TGraph *Rgeov30 = new TGraph(intR,xxx,volR30);
    Rgeov30->SetLineColor(2);
    Rgeov30->Draw("C");
    TGraph *Rgeohf30 = new TGraph(intR,xxx,HFR30);
    Rgeohf30->SetLineColor(3);
    Rgeohf30->Draw("C");
    geoline.Draw();
    TGraph *Rgeos30 = new TGraph(intR,xxx,surR30);
    Rgeos30->SetLineColor(4);
    Rgeos30->Draw("C");

    CgeoR->cd(3);
    HgeoR.Draw();
    Elab = 50.;
    double HFR50[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR50[i] = Volume.DispersiveCorrection(xxx[i])
                  + SurfaceHigh.DispersiveCorrection(xxx[i]);
        surR50[i] = SurfaceLow.DispersiveCorrection(xxx[i]);
        HFR50[i] = HartreeFock.RealPotential(xxx[i]);
	totR50[i] = surR50[i] + volR50[i] + HFR50[i];
      }
    TGraph *Rgeo50 = new TGraph(intR,xxx,totR50);
    Rgeo50->Draw("C");
    TGraph *Rgeov50 = new TGraph(intR,xxx,volR50);
    Rgeov50->SetLineColor(2);
    Rgeov50->Draw("C");
    TGraph *Rgeohf50 = new TGraph(intR,xxx,HFR50);
    Rgeohf50->SetLineColor(3);
    Rgeohf50->Draw("C");
    geoline.Draw();
    TGraph *Rgeos50 = new TGraph(intR,xxx,surR50);
    Rgeos50->SetLineColor(4);
    Rgeos50->Draw("C");

  CgeoR->cd(4);
  HgeoR.Draw();
    Elab = 75.;
    double HFR75[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR75[i] = Volume.DispersiveCorrection(xxx[i])
                  + SurfaceHigh.DispersiveCorrection(xxx[i]);
        surR75[i] = SurfaceLow.DispersiveCorrection(xxx[i]);
	totR75[i] = surR75[i] + volR75[i];
        HFR75[i] = HartreeFock.RealPotential(xxx[i]);
	totR75[i] = surR75[i] + volR75[i] + HFR75[i];
      }
    TGraph *Rgeo75 = new TGraph(intR,xxx,totR75);
    Rgeo75->Draw("C");
    TGraph *Rgeov75 = new TGraph(intR,xxx,volR75);
    Rgeov75->SetLineColor(2);
    Rgeov75->Draw("C");
    TGraph *Rgeohf75 = new TGraph(intR,xxx,HFR75);
    Rgeohf75->SetLineColor(3);
    Rgeohf75->Draw("C");
    geoline.Draw();
    TGraph *Rgeos75 = new TGraph(intR,xxx,surR75);
    Rgeos75->SetLineColor(4);
    Rgeos75->Draw("C");

   CgeoR->cd(5);
   HgeoR.Draw();
    Elab = 100.;
    double HFR100[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);
    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR100[i] = Volume.DispersiveCorrection(xxx[i])
                   + SurfaceHigh.DispersiveCorrection(xxx[i]);
        surR100[i] = SurfaceLow.DispersiveCorrection(xxx[i]);
        HFR100[i] = HartreeFock.RealPotential(xxx[i]);
	totR100[i] = surR100[i] + volR100[i] + HFR100[i];
      }
    TGraph *Rgeo100 = new TGraph(intR,xxx,totR100);
    Rgeo100->Draw("C");
    TGraph *Rgeov100 = new TGraph(intR,xxx,volR100);
    Rgeov100->SetLineColor(2);
    Rgeov100->Draw("C");
    TGraph *Rgeohf100 = new TGraph(intR,xxx,HFR100);
    Rgeohf100->SetLineColor(3);
    Rgeohf100->Draw("C");
    geoline.Draw();
    TGraph *Rgeos100 = new TGraph(intR,xxx,surR100);
    Rgeos100->SetLineColor(4);
    Rgeos100->Draw("C");

    CgeoR->cd(6);
    HgeoR.Draw();
    Elab = 200.;
    double HFR200[intR];
    InitializeForEcm(energyLab2Cm(Elab),Elab);

    for (int i = 0;i<intR;i++)    
      {
	xxx[i] = (double)i/7;
        volR200[i] = Volume.DispersiveCorrection(xxx[i])
                   + SurfaceHigh.DispersiveCorrection(xxx[i]); 
        surR200[i] = SurfaceLow.DispersiveCorrection(xxx[i]);
        HFR200[i] = HartreeFock.RealPotential(xxx[i]);
	totR200[i] = surR200[i] + volR200[i] + HFR200[i];
      }
    TGraph *Rgeo200 = new TGraph(intR,xxx,totR200);
    Rgeo200->Draw("C");
    TGraph *Rgeov200 = new TGraph(intR,xxx,volR200);
    Rgeov200->SetLineColor(2);
    Rgeov200->Draw("C");
    TGraph *Rgeohf200 = new TGraph(intR,xxx,HFR200);
    Rgeohf200->SetLineColor(3);
    Rgeohf200->Draw("C");
    geoline.Draw();
    TGraph *Rgeos200 = new TGraph(intR,xxx,surR200);
    Rgeos200->SetLineColor(4);
    Rgeos200->Draw("C");
    CgeoR->Write();
    
#endif

}
//**************************************************************************
  /**
   * opens root file to save spectra in
   */
void reaction::OpenRootFile()
{
#ifdef root
  string filename(title + ".root");
  cout << filename << endl;
  ifstream file (filename.c_str(),ios::in);

  f = new TFile(filename.c_str(),"RECREATE");
#endif
}
//**************************************************************************
  /**
   * closes the root file
   */
void reaction::CloseRootFile()
{
#ifdef root
  f->Write();
  f->Close();
#endif
}
//***********************************************************
void reaction::plotSmatrix()
{
  cout << "smatrix" << endl;

#ifdef root
  string nnn = "smat_"+title;
  TCanvas *fit = new TCanvas(nnn.c_str());
  TH2S *hist = new TH2S("histtt","",10,0,10.
                        ,10,0.,1.);
  hist->GetXaxis()->SetTitle("b_{n} [fm]");
  hist->GetYaxis()->SetTitle("1-|S_{n}|^{2}");
  hist->SetStats(kFALSE);
  hist->Draw();


  double Elab[6]={50.,60.,70.,80.,90.,100.};
  TGraph * graph[6];
  double x[100];
  double y[100];


  string filename = title+"sm.dat";
  ofstream fdata(filename.c_str());
  fdata << 100 << endl;
  for (int i=0;i<6;i++)
    {

     double Ecm = energyLab2Cm(Elab[i]);
     InitializeForEcm(Ecm,Elab[i]);
     for (int ib=0;ib<100;ib++)
       {
	x[ib] = (double)ib/10.;
        y[ib] = scatter->getSmatrixEikonal(x[ib]);
	if ( Elab[i] == 70.) fdata << x[ib] << " " << y[ib] << endl; 

       }
        graph[i] = new TGraph(100,x,y);
        graph[i]->SetLineColor(i+1);
        graph[i]->Draw("C");
    }
  fit->Write();

  filename = title+"_smat.dat";
  ofstream smfile(filename.c_str());
  double ElabX = 70.52;
  double Ecm = energyLab2Cm(ElabX);
  InitializeForEcm(Ecm,ElabX);  
  for (int ib=0;ib<300;ib++)
    {
      double b = (double)(ib+1)/300.*20.;
      scatter->getSmatrixEikonal(b);
      double sreal= scatter->Sreal;
      double simag= scatter->Simag;
      smfile << b << "  " << sreal << " " << simag << endl;
    }
  smfile.close();
  smfile.clear();

  filename = title+"_potR.dat";
  ofstream smfileR(filename.c_str());
  filename = title+"_potI.dat";
  ofstream smfileI(filename.c_str());

  scatter->LdotSigma = 0.;

  double deltar = .1;
  int nn = 200;
  smfileR << nn << endl;
  smfileI << nn << endl;
  for (int i=0;i<nn;i++)
    {
      double r = (double)i*deltar;
      double realP = scatter->HartreeFock->RealPotential(r)
	+scatter->Volume->DispersiveCorrection(r)
        +scatter->SurfaceLow->DispersiveCorrection(r)
	+scatter->SurfaceHigh->DispersiveCorrection(r);
      double imagP = scatter->ImaginaryPotential(r)/scatter->gammaRel;
      smfileR << r << " " << realP << endl;
      smfileI << r << " " << imagP << endl;
    }
  smfileR.close();
  smfileI.close();
  smfileR.clear();
  smfileI.clear();
  #endif
}
// 
void reaction::printSmatrix()
{
  cout << "printsmatrix" << endl;
  ofstream fs("smatrix.dat");
  for (int i=1;i<=180;i++)
    {
      double Ecm = (float)i;
      fs << "Ecm = " << Ecm << " MeV" << endl;
     double ElabX = energyCm2Lab(Ecm);
     InitializeForEcm(Ecm,ElabX);  
     // integrate wavefunction and find phaseshifts
     scatter->integrateWave();
     for (int l=0;l<=19;l++)
       {
         if (l != 0)fs << scatter->eta[l][0].real() << " " << scatter->eta[l][0].imag() << endl;
	 else fs << scatter->eta[l][1].real() << " " << scatter->eta[l][1].imag() << endl;
         fs << scatter->eta[l][1].real() << " " << scatter->eta[l][1].imag() << endl;

         //fs << l << " " << scatter->eta[l][0] << " " << scatter->eta[l][1] << endl;
        }
    }

}

void reaction::printSmatrix2()
{
  cout << "printsmatrix" << endl;
  ofstream fs("smatrix2.dat");
  for (int i=1;i<=5;i++)
    {
     double Elab = (float)i*20.;

     double Ecm = energyLab2Cm(Elab);
     InitializeForEcm(Ecm,Elab);  
     // integrate wavefunction and find phaseshifts
     scatter->integrateWave();
     fs << "Elab = " << Elab << " Ecm = " << Ecm << endl;
     fs << " l j  real imag" << endl;
     for (int l=0;l<=19;l++)
       {
         if (l != 0) fs << l << " " << (float)l-0.5 << " " << scatter->eta[l][0].real() << " " << scatter->eta[l][0].imag() << endl;
         fs << l << " " << (float)l + 0.5 << " " << scatter->eta[l][1].real() << " " << scatter->eta[l][1].imag() << endl;


        }
    }

}
