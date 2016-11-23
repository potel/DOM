using namespace std;

#include "fit.h"

string const fit::label[TotPara] = {"rC0","rC1","fEcoul","AHF0","AHF1","alpha",
				    "alphaA","alphaNZ","beta","betaA",
				    "gamma","gammaA","alphaS","alphaSA",
				"alphaSNZ","betaS","betaSA","gammaS","gammaSA",
				"rHF0","rHF1","rHFNZ","aHF0","aHF1",
                                "fGap",
				 "Bsurface","BsurfaceA","Csurface","CsurfaceA",
                                "Dsurface",
                                  "Asurface9P","Asurface9N",
				  "Asurface40P","Asurface40N",
				  "Asurface42P","Asurface44P",
				  "Asurface48P","Asurface48N",
                                  "Asurface50P","Asurface52P","Asurface52N",
                                  "Asurface54P","Asurface54N",
				  "Asurface58P","Asurface58N",
				  "Asurface60P","Asurface60N",
				  "Asurface62P","Asurface62N",
				  "Asurface64P","Asurface64N",
                                  "Asurface88P",
				  "Asurface90P","Asurface92P","Asurface92N",
				  "Asurface112P","Asurface114P",
				  "Asurface116P","Asurface116N",
                                  "Asurface118P","Asurface118N",
				  "Asurface120P","Asurface120N",
				  "Asurface122P","Asurface124P","Asurface124N",
                                  "Asurface206P",
				  "Asurface208P","Asurface208N",
			          "rsurface0","rsurface1","asurface",
                                    "Epvolume",
                                    "aHigh","aHighNZ","bHigh","cHigh",
			            "Avolume_0","Avolume_1","AvolumeA3",
			            "Bvolume_0","Bvolume_1","BvolumeA3",
			            "rvolume0","rvolumeNZ","rvolume1",
                                    "deltaR","deltaRNZ","deltaRA",
                                    "expR","avolume0","avolume1",
			            "Ea","alphaOverA",
                                    "Vso","VsoNZ","VspinoE",
			            "rso0","rso1","aso0","aso1","AWso","BWso"};

//******************************************************************
  /**
   * calculates the total of the chi-squared which is minimized
   *\para para is the array of parameter values
   */
double fit::functND(double*para)
{
  
  for (int i=0;i<ND;i++) cout << para[i] << " ";
  cout << endl;
  SavePara(para);

  bool ok;
  for (int i=0;i<Nreact;i++) 
    {

     try
       {
        ok = React[i]->prepare();
       }
     catch(localityException & locExcept)
       {
	 ok = 0;
       } 
     if (ok == 0) return 1e30; //if we cannot adjust fermi level, 
                                //then return large chi
    }
  double out =  0.;
  double chi = 0.;
  for (int i=0;i<Nreact;i++)
    {

     try
       {
       chi=React[i]->ChiSquared();
      }
 
     catch(localityException & locExcept)
       {
         chi = 1000.;
       }

     if (std::isnan(chi)) chi =1000.;
     out += chi;
    }

  // increase chisq if parameters are outside of reasonable bounds
  out += chiPara();

  cout << "chisq= " << out << endl;

  return out;

}

/**
 *Constructor - reads in experimental data and parameter values
 *\para title gives name of file with initial parameters
 *\para n is number of parameters fit
 */
fit::fit(string *title0,int n) : minimizeND(n)
{
  title = *title0;
  string filename(title + ".inp");
  cout << filename << endl;
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << endl;
      abort();
    }
  int mvolume;
  file  >> mvolume;

  string variable;
  for (int i=0;i<TotPara;i++)
    {
      file >> allPara[i] >> varied[i] >> squared[i] >> scaling[i] >> variable;
      if (squared[i]) allPara[i] = sqrt(allPara[i]); 
    }


  //create maps of all paremeter to fitted parameters
  Ndim = 0;
  for (int i=0;i<TotPara;i++)
    {
      map1[i] = -1;
      map2[i] = -1;
      if (varied[i]) 
	{
	 map1[i] = Ndim;
         map2[Ndim] = i;
         Ndim++;
	}
    }
  file.close();
  file.clear();



  //read in all the reactions which we will fit
  file.open("reactions.inp");
  file >> Nreact;
  bool yes = 1;
  for (int i=0;i<Nreact;i++)
    {
      file >> Rtitle[i];
      React[i] = new reaction(&(Rtitle[i]),yes);
      React[i]->mvolume = mvolume;
    }
  file.close();
  file.clear();

  decodePara();
 
}
//******************************************************************
  /**
   * Destructor
   */
fit::~fit()
{
  for (int i=0;i<Nreact;i++) delete React[i];
}

void fit::decodePara()
{
  /**
   * decodes the parameter file input and send this information to the 
   * classes associated with each potential
   */
  int index = 0;
  rc0 = allPara[index++];
  rc1 = allPara[index++];
  fEc = allPara[index++];
  AHF0 = allPara[index++];
  AHF1 = allPara[index++];
  alpha = allPara[index++];
  alphaA = allPara[index++];
  alphaNZ = allPara[index++];
  beta = allPara[index++];
  betaA = allPara[index++];
  gamma = allPara[index++];
  gammaA = allPara[index++];
  alphaS = allPara[index++];
  alphaSA = allPara[index++];
  alphaSNZ = allPara[index++];
  betaS = allPara[index++];
  betaSA = allPara[index++];
  gammaS = allPara[index++];
  gammaSA = allPara[index++];
  rHF0 = allPara[index++];
  rHF1 = allPara[index++];
  rHFNZ = allPara[index++];
  aHF0 = allPara[index++];
  aHF1 = allPara[index++];

  fGap = allPara[index++];  
  Bsurface = allPara[index++];  
  BsurfaceA = allPara[index++];  
  Csurface     = allPara[index++];
  CsurfaceA     = allPara[index++];
  Dsurface     = allPara[index++];
  Asurface9P = allPara[index++];
  Asurface9N = allPara[index++];
  Asurface40P = allPara[index++];
  Asurface40N = allPara[index++];
  Asurface42P = allPara[index++];
  Asurface44P = allPara[index++];
  Asurface48P = allPara[index++];
  Asurface48N = allPara[index++];
  Asurface50P = allPara[index++];
  Asurface52P = allPara[index++];
  Asurface52N = allPara[index++];
  Asurface54P = allPara[index++];
  Asurface54N = allPara[index++];
  Asurface58P = allPara[index++];
  Asurface58N = allPara[index++];
  Asurface60P = allPara[index++];
  Asurface60N = allPara[index++];
  Asurface62P = allPara[index++];
  Asurface62N = allPara[index++];
  Asurface64P = allPara[index++];
  Asurface64N = allPara[index++];
  Asurface88P = allPara[index++];
  Asurface90P = allPara[index++];
  Asurface92P = allPara[index++];
  Asurface92N = allPara[index++];
  Asurface112P = allPara[index++];
  Asurface114P = allPara[index++];
  Asurface116P = allPara[index++];
  Asurface116N = allPara[index++];
  Asurface118P = allPara[index++];
  Asurface118N = allPara[index++];
  Asurface120P = allPara[index++];
  Asurface120N = allPara[index++];
  Asurface122P = allPara[index++];
  Asurface124P = allPara[index++];
  Asurface124N = allPara[index++];
  Asurface206P = allPara[index++];
  Asurface208P = allPara[index++];
  Asurface208N = allPara[index++];
  rsurface0 = allPara[index++];
  rsurface1 = allPara[index++];
  asurface = allPara[index++];
  Epvolume = pow(allPara[index++],2);
  aHigh = allPara[index++];
  aHighNZ = allPara[index++];
  bHigh = allPara[index++];
  cHigh = allPara[index++];
  Avolume0 = allPara[index++];
  Avolume1 = allPara[index++];
  AvolumeA3 = allPara[index++];
  Bvolume0 = allPara[index++];
  Bvolume1 = allPara[index++];
  BvolumeA3 = allPara[index++];
  rvolume0 = allPara[index++];
  rvolumeNZ = allPara[index++];
  rvolume1 = allPara[index++];
  deltaRvolume = allPara[index++];
  deltaRvolumeNZ = allPara[index++];
  deltaRvolumeA = allPara[index++];
  expRvolume = allPara[index++];
  avolume0 = allPara[index++];
  avolume1 = allPara[index++];
  Ea = allPara[index++];
  alphaOverA = allPara[index++];
  Vso = allPara[index++];
  VsoNZ = allPara[index++];
  VspinoE = allPara[index++];
  rso0 = allPara[index++];
  rso1 = allPara[index++];
  aso0 = allPara[index++];
  aso1 = allPara[index++];
  AWso = allPara[index++];
  BWso = allPara[index++];

  
  for (int i=0;i<Nreact;i++)
    {
     React[i]->Rc = pow(React[i]->A,1./3.)*rc0 + rc1;
     React[i]->VHF = AHF0+ React[i]->asymmetry*React[i]->sign*AHF1;

     //sign = 1 for protons and -1 for neutrons
     //asymmetry = (N-Z)/A

     React[i]->alpha = alpha+ alphaNZ*React[i]->asymmetry*React[i]->sign
     +alphaA*React[i]->A;
     React[i]->beta = beta+betaA/pow(React[i]->A,1./3.);
     React[i]->gamma = gamma+gammaA/pow(React[i]->A,1./3.);
     React[i]->alphaS = alphaS + alphaSNZ*React[i]->asymmetry*React[i]->sign
                        + alphaSA*React[i]->A;
     React[i]->betaS = betaS + betaSA*React[i]->A;
     React[i]->gammaS = gammaS + React[i]->A*gammaSA;
     React[i]->RHF = pow(React[i]->A,1./3.)*(rHF0+
     rHFNZ*React[i]->asymmetry*React[i]->sign) + rHF1;
     React[i]->aHF = aHF0 +aHF1/pow(React[i]->A,1./3.);
     if (React[i]->Zp == 1)
       {
         React[i]->coulShift = (1.73*React[i]->Z/React[i]->Rc)*fEc;
       }
     else React[i]->coulShift = 0.;

     React[i]->fGap = fGap;
     React[i]->Bsurface = Bsurface + BsurfaceA*React[i]->A;
     React[i]->Csurface = Csurface + CsurfaceA/pow(React[i]->A,1./3.);
     React[i]->Dsurface = Dsurface;
     if (React[i]->Zp == 1 && React[i]->A == 206)
       {
         React[i]->Asurface = Asurface206P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 208)
       {
         React[i]->Asurface = Asurface208P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 208)
       {
         React[i]->Asurface = Asurface208N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 40)
       {
         React[i]->Asurface = Asurface40P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 39)
       {
         React[i]->Asurface = Asurface40P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 40)
       {
         React[i]->Asurface = Asurface40N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 48)
       {
         React[i]->Asurface = Asurface48P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 48)
       {
         React[i]->Asurface = Asurface48N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 42)
       {
         React[i]->Asurface = Asurface42P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 44)
       {
         React[i]->Asurface = Asurface44P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 50)
       {
         React[i]->Asurface = Asurface50P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 52)
       {
         React[i]->Asurface = Asurface52P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 52)
       {
         React[i]->Asurface = Asurface52N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 58)
       {
         React[i]->Asurface = Asurface58P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 60)
       {
         React[i]->Asurface = Asurface60P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 62)
       {
         React[i]->Asurface = Asurface62P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 64)
       {
         React[i]->Asurface = Asurface64P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 58)
       {
         React[i]->Asurface = Asurface58N;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 60)
       {
         React[i]->Asurface = Asurface60N;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 62)
       {
         React[i]->Asurface = Asurface62N;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 64)
       {
         React[i]->Asurface = Asurface64N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 54)
       {
         React[i]->Asurface = Asurface54P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 54)
       {
         React[i]->Asurface = Asurface54N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 88)
       {
         React[i]->Asurface = Asurface88P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 90)
       {
         React[i]->Asurface = Asurface90P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 92)
       {
         React[i]->Asurface = Asurface92P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 92)
       {
         React[i]->Asurface = Asurface92N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 112)
       {
         React[i]->Asurface = Asurface112P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 114)
       {
         React[i]->Asurface = Asurface114P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 116)
       {
         React[i]->Asurface = Asurface116P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 116)
       {
         React[i]->Asurface = Asurface116N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 118)
       {
         React[i]->Asurface = Asurface118P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 118)
       {
         React[i]->Asurface = Asurface118N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 120)
       {
         React[i]->Asurface = Asurface120P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 120)
       {
         React[i]->Asurface = Asurface120N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 122)
       {
         React[i]->Asurface = Asurface122P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 124)
       {
         React[i]->Asurface = Asurface124P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 124)
       {
         React[i]->Asurface = Asurface124N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 9)
       {
         React[i]->Asurface = Asurface9P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 9)
       {
         React[i]->Asurface = Asurface9N;
       }
     else 
       {
         cout << "reaction not implimented " << endl;
         cout << React[i]->Zp << " " << React[i]->A << endl;
	 abort();
       }

     React[i]->Rsurface = pow(React[i]->A,1./3.)*rsurface0 + rsurface1;


     React[i]->Epvolume = Epvolume;

     React[i]->Avolume = Avolume0 + 
       React[i]->asymmetry*React[i]->sign*Avolume1 + 
       AvolumeA3/pow(React[i]->A,1./3.);
     React[i]->Bvolume = Bvolume0 + 
     React[i]->asymmetry*React[i]->sign*Bvolume1 + 
       BvolumeA3/pow(React[i]->A,1./3.);


     React[i]->aHigh = aHigh
	    +aHighNZ*React[i]->asymmetry*React[i]->sign;
     if (bHigh == -1.) React[i]->bHigh = React[i]->Bvolume;
     else React[i]->bHigh = bHigh;
     React[i]->cHigh = cHigh;

     React[i]->Rvolume = pow(React[i]->A,1./3.)*(rvolume0 + 
           React[i]->asymmetry*React[i]->sign*rvolumeNZ)
           + rvolume1;
     React[i]->deltaRvolume = deltaRvolume+ deltaRvolumeA*React[i]->A + 
	    React[i]->asymmetry*React[i]->sign*deltaRvolumeNZ;

     React[i]->expRvolume = expRvolume;
     React[i]->avolume = avolume0 + avolume1/pow(React[i]->A,1./3.);
     React[i]->asurface = React[i]->avolume;
     React[i]->EaVolume = Ea;
     //alphaOverA = alpha/Avolume
     React[i]->alphaVolume = alphaOverA*React[i]->Avolume;

     React[i]->Vso = Vso + React[i]->asymmetry*React[i]->sign*VsoNZ;
     React[i]->VspinoE = VspinoE; 
     React[i]->Rso = pow(React[i]->A,1./3.)*rso0 + rso1;
     React[i]->aso = aso0 + aso1/pow(React[i]->A,1./3.);
     React[i]->AWso = AWso;
     React[i]->BWso = BWso;
    }
}

//*********************************************************************
  /**
   * saves the fitted and nonfitted parameters in the same form as
   * the input file
   */
void   fit:: SavePara(double * para)
    {
     for (int i=0;i<Ndim;i++)
       {
         int j = map2[i];
         allPara[j] =para[i]/scaling[j];
        }

     decodePara();
    }
//******************************************************************
  /**
   * Prints to screen the fit parameters
   *
   */
void fit::PrintFit()
{
  for (int i=0;i<TotPara;i++)
    {
      double value = allPara[i];
      if (squared[i]) value = pow(value,2);  
      cout << setw(17) << value << " " << 
             setw(3) << varied[i] << " " << 
             setw(3) << squared[i] << " " <<
             setw(17) << scaling[i] << " " << 
             setw(10) << label[i] << endl;
    }
}
//******************************************************************
  /**
   * writes out the fitted parameter in the file title.out
   * This file can later renamed title.inp and used as the 
   * input file for another fit 
   *
   */

void fit::WriteFit(double chiMin)
{

  //make a linear fit to the AHF values obtained from the main fit
  //fitAHF();

  string filename(title + ".out");
  cout << filename << endl;
  ofstream file(filename.c_str(),ios::out);
  file.setf(ios::fixed); //used fixed precision
  file.precision(11); // set significant digits
  file << setw(10) << React[0]->mvolume << endl;


  for (int i=0;i<TotPara;i++)
    {
      double value = allPara[i];
      if (squared[i]) value = pow(value,2);  
    file << setw(17) << value << " " << 
            setw(3) << varied[i] << " " << 
            setw(3) << squared[i] << " " <<
            setw(17) << scaling[i] << " " << 
            setw(10) << label[i] << endl;
    }

  file << " " << endl;
  file << " minimun chi squared = " << chiMin << endl;

  file.close();

}
//*************************************
  /**
   * makes a linear fit to the AHF values obtained from the fit
   * This gives the \f$ \frac{N-Z}{A}^{2} \f$ dependence of the depth of the 
   * Hartee Foack potential. 
   */
void fit::fitAHF()
{
  cout << "fitAHF" << endl;
  double sumx = 0.;
  double sumxy = 0.;
  double sumy = 0.;
  double sumx2 = 0.;

  for (int i=0;i<Nreact;i++)
    {
      double x =  React[i]->asymmetry*React[i]->sign;
      double y = React[i]->VHF;
      cout << "x= " << x << " y= " << y << endl;
      sumx += x;
      sumy += y;
      sumxy += x*y;
      sumx2 += pow(x,2);
       
    }
  double denominator = (double)Nreact*sumx2 - pow(sumx,2);
  allPara[1] = (sumy*sumx2 - sumx*sumxy)/denominator;
  allPara[2] = ((double)Nreact*sumxy - sumx*sumy)/denominator;  
  cout << sumx << " " << sumy << " " << sumxy << " " << sumx2 << endl;
  cout << denominator << endl;
  cout << allPara[1] << " " << allPara[2] << endl;
}
//*************************************************
double fit::chiPara()
{
  double out = 0.;
  if (asurface < 0.5 || asurface > .8) out += 100.;
  if (avolume0 < 0.5 || avolume0 > .8) out += 100.;
  if (aso0 < 0.5 || aso0 > 0.8) out += 100.;

  if (Avolume1 < 0.) out += 100.;
  return out; 
   
}
