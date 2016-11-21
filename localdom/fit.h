#include "minimizeND.h"
#include <string>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "reaction.h"
#include "localityException.h"
using namespace std;

//***************************************************************************
  /**
   *\brief main fitting class
   *
   * performs dispersive optical model fits to elastic scattering, reaction 
   * and toal sections, bound state energies, RMS radii, widths, and 
   * spectroscopic factors
   */
class fit: public minimizeND 
{
 public:
  fit (string*,int);
  ~fit();
  string title;
  void decodePara();
  void PrintFit();
  void WriteFit(double);
  void SavePara(double *);
  double functND(double*); 


  static int const TotPara=103; //!< number of  parameters
  int Ndim; //!<Number of actual fit parameters
  int varied[TotPara]; //!<which parameters are fit
  int map1[TotPara];  //!<map from fitted parameters to all parameters
  int map2[TotPara];  //!<map from all parametrs to fitted parameters
  double scaling[TotPara]; //!<scaling of fit parameters to avoid large jumps
  double allPara[TotPara]; //!<values of all the parameters
  bool squared[TotPara]; //!<keep this parameter positive 
  void fitAHF(); 
  double chiPara();


  // reactions to fit
  int Nreact; //!< number of reactions to fit
  reaction * React[30]; //!< reaction class
  string Rtitle[30]; //!< titles of the reactions


 private:
  static string const label[TotPara]; 


  double rc0;
  double rc1;
  double fEc;
  double AHF0;
  double AHF1;
  double alpha;
  double alphaA;
  double alphaNZ;
  double beta;
  double betaA;
  double gamma;
  double gammaA;
  double alphaS;
  double alphaSA;
  double alphaSNZ;
  double betaS;
  double betaSA;
  double gammaS;
  double gammaSA;
  double rHF0;
  double rHF1;
  double rHFNZ;
  double aHF0;
  double aHF1;

  double fGap;
  double Bsurface;
  double BsurfaceA;
  double Csurface;
  double CsurfaceA;
  double Dsurface;
  double Asurface206P;
  double Asurface208P;
  double Asurface208N;
  double Asurface40P;
  double Asurface40N;
  double Asurface48P;
  double Asurface48N;
  double Asurface42P;
  double Asurface44P;
  double Asurface50P;
  double Asurface52P;
  double Asurface52N;
  double Asurface54P;
  double Asurface54N;
  double Asurface58P;
  double Asurface58N;
  double Asurface60P;
  double Asurface60N;
  double Asurface62P;
  double Asurface62N;
  double Asurface64P;
  double Asurface64N;
  double Asurface88P;
  double Asurface90P;
  double Asurface92P;
  double Asurface92N;
  double Asurface112P;
  double Asurface114P;
  double Asurface116P;
  double Asurface116N;
  double Asurface118P;
  double Asurface118N;
  double Asurface120P;
  double Asurface120N;
  double Asurface122P;
  double Asurface124P;
  double Asurface124N;
  double Asurface9P;
  double Asurface9N;
  double rsurface0;
  double rsurface1;
  double asurface;
  double  Epvolume;
  double aHigh;
  double aHighNZ;
  double bHigh;
  double cHigh;
  double Avolume0;
  double Avolume1;
  double AvolumeA3;
  double  Bvolume0;
  double  Bvolume1;
  double BvolumeA3;
  double rvolume0;
  double rvolumeNZ;
  double rvolume1;
  double deltaRvolume;
  double deltaRvolumeNZ;
  double deltaRvolumeA;
  double  expRvolume;
  double avolume0;
  double avolume1;
  double Ea;
  double alphaOverA;
  double  Vso;
  double  VsoNZ;
  double  VspinoE;
  double rso0;
  double rso1;
  double aso0;
  double aso1;
  double AWso;
  double BWso;

 
};
