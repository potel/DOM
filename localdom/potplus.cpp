#include "pot.h"


/**
 *
 * constructor
 */
pot::pot()
{
  string title("../ca/fitL8_23");  // name of input file of fit parameters
  string filename(title + ".inp");
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << endl;
      abort();
    }
  int mvolume;
  file >> mvolume;


  //CHANGE THIS TO nca40 IF YOU WANT THE NEUTRON POTENTIAL
  title = "../ca/nca60";
  bool flag=1;
  React[0] = new reaction(&title,flag);
  React[0]->mvolume = mvolume;
  Nreact = 1;

  string variable;
  for (int i=0;i<TotPara;i++)
    {
      file >> allPara[i] >> varied[i] >> squared[i] >> scaling[i] >> variable;
      if (squared[i]) allPara[i] = sqrt(allPara[i]); 
    }
  file.close();
  file.clear();
  decodePara();
  EcmOld = -99999999.;
}
pot::pot(int uno)
{
  string title("../ca/fitL8_23");  // name of input file of fit parameters
  string filename(title + ".inp");
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << endl;
      abort();
    }
  int mvolume;
  file >> mvolume;


  //CHANGE THIS TO nca40 IF YOU WANT THE NEUTRON POTENTIAL
  if(uno==0) {title = "../ca/nca40"; cout<<"Implementing neutron DOM potential for 40Ca"<<endl;}
  if(uno==1) {title = "../ca/pca40"; cout<<"Implementing proton DOM potential for 40Ca"<<endl;}
  if(uno==2) {title = "../ca/nca48"; cout<<"Implementing neutron DOM potential for 48Ca"<<endl;}
  if(uno==3) {title = "../ca/pca48"; cout<<"Implementing proton DOM potential for 48Ca"<<endl;}
  if(uno==4) {title = "../ca/nca60"; cout<<"Implementing neutron DOM potential for 60Ca"<<endl;}
  if(uno==5) {title = "../ca/pca60"; cout<<"Implementing proton DOM potential for 60Ca"<<endl;}
  if(uno<0||uno>5) {cout<<"Error in pot.cpp: choose between protons  or neutrons for implemented Ca isotopes";exit(0);}
  bool flag=1;
  React[0] = new reaction(&title,flag);
  React[0]->mvolume = mvolume;
  Nreact = 1;

  string variable;
  for (int i=0;i<TotPara;i++)
    {
      file >> allPara[i] >> varied[i] >> squared[i] >> scaling[i] >> variable;
      if (squared[i]) allPara[i] = sqrt(allPara[i]); 
    }
  file.close();
  file.clear();
  decodePara();
  EcmOld = -99999999.;
}

//*******************************************************************
  /**
   * Destructor
   */
pot::~pot()
{
  delete React[0];
}

void pot::decodePara()
{
  int index = 0;
  double rc0 = allPara[index++];
  double rc1 = allPara[index++];
  double fEc = allPara[index++];
  double AHF0 = allPara[index++];
  double AHF1 = allPara[index++];
  double alpha = allPara[index++];
  double alphaA = allPara[index++];
  double alphaNZ = allPara[index++];
  double beta = allPara[index++];
  double betaA = allPara[index++];
  double gamma = allPara[index++];
  double gammaA = allPara[index++];
  double alphaS = allPara[index++];
  double alphaSA = allPara[index++];
  double alphaSNZ = allPara[index++];
  double betaS = allPara[index++];
  double betaSA = allPara[index++];
  double gammaS = allPara[index++];
  double gammaSA = allPara[index++];
  double rHF0 = allPara[index++];
  double rHF1 = allPara[index++];
  double rHFNZ = allPara[index++];
  double aHF0 = allPara[index++];
  double aHF1 = allPara[index++];

  double fGap = allPara[index++];  
  double Bsurface = allPara[index++];  
  double BsurfaceA = allPara[index++];  
  double Csurface     = allPara[index++];
  double CsurfaceA     = allPara[index++];
  double Dsurface     = allPara[index++];
  double Asurface9P = allPara[index++];
  double Asurface9N = allPara[index++];
  double Asurface40P = allPara[index++];
  double Asurface40N = allPara[index++];
  double Asurface42P = allPara[index++];
  double Asurface44P = allPara[index++];
  double Asurface48P = allPara[index++];
  double Asurface48N = allPara[index++];
  double Asurface50P = allPara[index++];
  double Asurface52P = allPara[index++];
  double Asurface52N = allPara[index++];
  double Asurface54P = allPara[index++];
  double Asurface54N = allPara[index++];
  double Asurface58P = allPara[index++];
  double Asurface58N = allPara[index++];
  double Asurface60P = allPara[index++];
  double Asurface60N = allPara[index++];
  double Asurface62P = allPara[index++];
  double Asurface62N = allPara[index++];
  double Asurface64P = allPara[index++];
  double Asurface64N = allPara[index++];
  double Asurface88P = allPara[index++];
  double Asurface90P = allPara[index++];
  double Asurface92P = allPara[index++];
  double Asurface92N = allPara[index++];
  double Asurface112P = allPara[index++];
  double Asurface114P = allPara[index++];
  double Asurface116P = allPara[index++];
  double Asurface116N = allPara[index++];
  double Asurface118P = allPara[index++];
  double Asurface118N = allPara[index++];
  double Asurface120P = allPara[index++];
  double Asurface120N = allPara[index++];
  double Asurface122P = allPara[index++];
  double Asurface124P = allPara[index++];
  double Asurface124N = allPara[index++];
  double Asurface206P = allPara[index++];
  double Asurface208P = allPara[index++];
  double Asurface208N = allPara[index++];
  double rsurface0 = allPara[index++];
  double rsurface1 = allPara[index++];
  double asurface = allPara[index++];
  double Epvolume = pow(allPara[index++],2);
  double aHigh = allPara[index++];
  double aHighNZ = allPara[index++];
  double bHigh = allPara[index++];
  double cHigh = allPara[index++];
  double Avolume0 = allPara[index++];
  double Avolume1 = allPara[index++];
  double AvolumeA3 = allPara[index++];
  double Bvolume0 = allPara[index++];
  double Bvolume1 = allPara[index++];
  double BvolumeA3 = allPara[index++];
  double rvolume0 = allPara[index++];
  double rvolumeNZ = allPara[index++];
  double rvolume1 = allPara[index++];
  double deltaRvolume = allPara[index++];
  double deltaRvolumeNZ = allPara[index++];
  double deltaRvolumeA = allPara[index++];
  double expRvolume = allPara[index++];
  double avolume0 = allPara[index++];
  double avolume1 = allPara[index++];
  double Ea = allPara[index++];
  double alphaOverA = allPara[index++];
  double Vso = allPara[index++];
  double VsoNZ = allPara[index++];
  double VspinoE = allPara[index++];
  double rso0 = allPara[index++];
  double rso1 = allPara[index++];
  double aso0 = allPara[index++];
  double aso1 = allPara[index++];
  double AWso = allPara[index++];
  double BWso = allPara[index++];
 
  int Nreact = 1;

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
         cout << "reaction not implemented " << endl;
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

     React[i]->prepare();
    }
}


/*
 *--------------------------------------------------------------------------------------
 *       Class:  pot
 *      Method:  pot :: justPara
 * Description:  Initialized to A-1,Z-1 system, doesn't run prepare, meaning levels 
 * aren't adjusted, so using the adjustments from (A,Z) system
 *--------------------------------------------------------------------------------------
 */
void pot::justPara()
{
   React[0]->A -= 1;
   React[0]->Z -= 1;
   cout<<"A = "<<React[0]->A<<" Z = "<<React[0]->Z<<endl;

   int index = 0;
   double rc0 = allPara[index++];
   double rc1 = allPara[index++];
   double fEc = allPara[index++];
   double AHF0 = allPara[index++];
   double AHF1 = allPara[index++];
   double alpha = allPara[index++];
   double alphaA = allPara[index++];
   double alphaNZ = allPara[index++];
   double beta = allPara[index++];
   double betaA = allPara[index++];
   double gamma = allPara[index++];
   double gammaA = allPara[index++];
   double alphaS = allPara[index++];
   double alphaSA = allPara[index++];
   double alphaSNZ = allPara[index++];
   double betaS = allPara[index++];
   double betaSA = allPara[index++];
   double gammaS = allPara[index++];
   double gammaSA = allPara[index++];
   double rHF0 = allPara[index++];
   double rHF1 = allPara[index++];
   double rHFNZ = allPara[index++];
   double aHF0 = allPara[index++];
   double aHF1 = allPara[index++];

   double fGap = allPara[index++];  
   double Bsurface = allPara[index++];  
   double BsurfaceA = allPara[index++];  
   double Csurface     = allPara[index++];
   double CsurfaceA     = allPara[index++];
   double Dsurface     = allPara[index++];
   double Asurface9P = allPara[index++];
   double Asurface9N = allPara[index++];
   double Asurface40P = allPara[index++];
   double Asurface40N = allPara[index++];
   double Asurface42P = allPara[index++];
   double Asurface44P = allPara[index++];
   double Asurface48P = allPara[index++];
   double Asurface48N = allPara[index++];
   double Asurface50P = allPara[index++];
   double Asurface52P = allPara[index++];
   double Asurface52N = allPara[index++];
   double Asurface54P = allPara[index++];
   double Asurface54N = allPara[index++];
   double Asurface58P = allPara[index++];
   double Asurface58N = allPara[index++];
   double Asurface60P = allPara[index++];
   double Asurface60N = allPara[index++];
   double Asurface62P = allPara[index++];
   double Asurface62N = allPara[index++];
   double Asurface64P = allPara[index++];
   double Asurface64N = allPara[index++];
   double Asurface88P = allPara[index++];
   double Asurface90P = allPara[index++];
   double Asurface92P = allPara[index++];
   double Asurface92N = allPara[index++];
   double Asurface112P = allPara[index++];
   double Asurface114P = allPara[index++];
   double Asurface116P = allPara[index++];
   double Asurface116N = allPara[index++];
   double Asurface118P = allPara[index++];
   double Asurface118N = allPara[index++];
   double Asurface120P = allPara[index++];
   double Asurface120N = allPara[index++];
   double Asurface122P = allPara[index++];
   double Asurface124P = allPara[index++];
   double Asurface124N = allPara[index++];
   double Asurface206P = allPara[index++];
   double Asurface208P = allPara[index++];
   double Asurface208N = allPara[index++];
   double rsurface0 = allPara[index++];
   double rsurface1 = allPara[index++];
   double asurface = allPara[index++];
   double Epvolume = pow(allPara[index++],2);
   double aHigh = allPara[index++];
   double aHighNZ = allPara[index++];
   double bHigh = allPara[index++];
   double cHigh = allPara[index++];
   double Avolume0 = allPara[index++];
   double Avolume1 = allPara[index++];
   double AvolumeA3 = allPara[index++];
   double Bvolume0 = allPara[index++];
   double Bvolume1 = allPara[index++];
   double BvolumeA3 = allPara[index++];
   double rvolume0 = allPara[index++];
   double rvolumeNZ = allPara[index++];
   double rvolume1 = allPara[index++];
   double deltaRvolume = allPara[index++];
   double deltaRvolumeNZ = allPara[index++];
   double deltaRvolumeA = allPara[index++];
   double expRvolume = allPara[index++];
   double avolume0 = allPara[index++];
   double avolume1 = allPara[index++];
   double Ea = allPara[index++];
   double alphaOverA = allPara[index++];
   double Vso = allPara[index++];
   double VsoNZ = allPara[index++];
   double VspinoE = allPara[index++];
   double rso0 = allPara[index++];
   double rso1 = allPara[index++];
   double aso0 = allPara[index++];
   double aso1 = allPara[index++];
   double AWso = allPara[index++];
   double BWso = allPara[index++];

   int Nreact = 1;

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

//      React[i]->VHF = 64.5;
      
      React[i]->load();
      int A = React[i]->A;
      int Z = React[i]->Z;
      React[i]->scatter->A = A;
      React[i]->scatter->Z = Z;
      React[i]->scatter->mu = A/(A+1.0)*1.0;
//Only difference here is that there is no prepare function called
   }
}

void pot::initialForNewEcm(double Ecm)
{
   double Elab = React[0]->energyCm2Lab(Ecm);
   double ecc = React[0]->energyLab2Cm(Elab);
   cout << Ecm << " " << Elab << " " << ecc << endl;
   React[0]->InitializeForEcm(Ecm, Elab);
}


void pot::GetPot(double r,int l, double j)
{
   React[0]->scatter->LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
   Real = React[0]->scatter->RealPotential(r);
   Imag = React[0]->scatter->ImaginaryPotential(r);


   Real /= React[0]->scatter->gammaRel;
   Imag /= React[0]->scatter->gammaRel;
}


void pot::potential(double r, int l, double j, double Ecm)
{
   if (Ecm != EcmOld) 
   {
      initialForNewEcm(Ecm);
      EcmOld = Ecm;
   }

   GetPot(r,l,j);
}
