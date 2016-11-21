#ifndef _scat
#define _scat

#include <iostream>
#include <complex>
#include <vector>
#include "legendre.h"
#include "merson.h"
#include "imaginaryForm.h"
#include "coul.h"
#include "whit.h"
#include "volume.h"
#include "surfaceTF.h"
#include "surfaceSTD.h"
#include "hartreeFock.h"
#include "spinOrbit.h"
#include "sphericalB.h"
#include "compound.h"

using namespace std;

//***************************************************************************
  /**
   * \brief Solves for wave function
   *
   * Integrates wavefunction to get scattering solutions a bound states
   */


class scat:  public compound, public merson
{
 public:
  static double const pi;
  static double const kconstant;
  static double const e2;
  static double const m0;
  static double const EminRel;
  bool flag;
  void initIntegration();

  int nlim,n0,kj,lp,molj,ipwia,nl,nstep;                 /* for schwandt potential */
  double ux,wx,wdx,cp,tp,cpi,tpi,cpd,tpd,gx,xkpf,hrs,hr,rc,vlx,ep0,relg; /* schwandt pot */
  double xp,usx,wsx,cs,ts,tis,cis;                  /* schwandt pot */

  double vs,rs0,as0,ws,rw0,aw0,vso,rso,aso,wso,rwo,awo; /* actualy schwandt parameters */

  double proton; // logical =1 for proton, =0 for neutron
  double mu; //reduced mass in nucleons
  double A; // target nucleon number
  double Z; // target proton number
  double Zp; // projectile proton number
  double Rc ; // Coulomb radius
  double energyCM; // cm energy MeV
  double energyLab; // lab energy in MeV
  double Kwave;    // wavenumber
  double konst;  // constant for cross section
  double muhbar; // 2*mu/hbar**2
  double gamma; // Sommerfeld parameter
  double Kwave2; //wavenumber squared
  double gammaRel; // relativistic gamma factor
  int l; // ellwave
  int lStop; // maxium l-wave used in calculation
  double j; // total angular
  double LdotSigma; // L.s term, sigma = spin*2
  double **phaseShift; // array of dimension [lMax+1][2] with phaseshifts;
  double **eta2; //array of dimension [lMax+1][2];
  complex<double> **eta; //array of dimension [lMax+1][2];
  double *Sigma;//array of dimension [lMax+1] containing Coulomb phase shifts
  int elastic; // if this for np scattering
  vector <double> rpot;                    /* real ee'p potential */
  vector <double> ipot;                   /* imaginary ee'p potential */
  complex <double> dompot[508];                /* Making DOM array to put it in cseq */
  complex <double> schwandt1[170];               /* Shorter array to use in merson */
  complex <double> schwandt2[508];               /* Longer arrary to use in cseq */
  complex<double> ca1[170];                     /* shorter array to use in merson */
  complex<double> ca2[508];                     /* longer array to use in cseq */
  complex<double> cdn[513];                    /* weird array for optso */
  complex<double> ff[170];                      /* Array for dweep partial wave */
  complex<double> gg[170];                      /* Array for dweep partial wave' */


  //public:
  // functions
  scat();
  virtual ~scat();
  scat(double,double,double,bool,string*);
  void loadPotential(hartreeFock*, volume*, surfaceTF*, surfaceSTD*, 
        spinOrbit*,double, double,double);
  void init(double,double,double,bool);
  double RealPotential(double);
  double ImaginaryPotential(double);
  double PlotRealPot(double);
  double PlotImaginaryPot(double);
  void Plot();
  void ReadPot();                               /* Reads input for eep pot */
  void SchwandtPot(int,int);                    /* runs optso twice */
  void ReadRealPot();
  void ReadImagPot();
  int integrateWave();
  int integrateWaveSch();
  int integrateWaveOptso();
  int integrateWaveCseq();
  int integrateWaveCseqDOM();
  valarray<double> IntegrateBound(double,int);
  void GetWaveFunctionArray(double,int,double&);
  void normalizeWaveFunction(double,double);
  double DifferentialXsection(double);
  double AnalyzePower; // Analysing power, calculated by DifferentialXsection
  double SpinRotation; // spin rotation parameter Q, calculated by Differential
  double AbsorptionXsection(); //absorption cross section
  double TotXsection(); //total cross section (for neutrons)
  double ElasticXsection(); // total elastic cross section (for neutrons)
  double Rutherford(double);
  double exteriorLogDer(int);
  void exteriorWaveFunct(int);
  double spectralWidth();
  double TransCoef(int, double);
  double MomWave(double);
  void GetMomWave();
  int Nodes();
  double isovectorAbsorb();// absorption xsec associated with isovector
  //surface imaginary
  double getSmatrixEikonal(double b);  
  double fws(double,double,double,int);                      /* code that returns wood-saxon */
  void ReadParams();                            /* reads actualy schwandt parameters */
  double RealSchwandt(double);                      /* returns real schwandt potential */
  double ImagSchwandt(double);                  /* returns imaginary schwandt potential */

  //effective mass function
  double eMass(double);
  double eMassHF(double);
  double eMassBar(double);
  double DerDelta(double);
  double Occupation; // occupation probability
  double SpectFactor; //spectroscopic factor
  double Width; 
  double Rrms;
  double AverageMassEff;
  double AverageInvMassBar;
  double AverageW;

  //complex<double> Potential(double);
  valarray<double> diff(double, valarray<double>);
  valarray<double> diffeep(double, valarray<double>); /* for schwandt */
  valarray<double> diffsch(double, valarray<double>); /* for schwandt */

  // variables and pointers
  volume *Volume;
  surfaceTF *SurfaceLow;
  surfaceSTD *SurfaceHigh;
  hartreeFock *HartreeFock;
  spinOrbit *SpinOrbit;

  //vintegrals
  double JReal;
  double JImag;
  double JSO;
  double RrmsReal;
  double RrmsImag;
  double RrmsSO;
  void VIntegrals();

  double ANC; //!< asymptotic normalization coefficient


  double *WaveArray; // Wavefunction array
  double *WaveBar;
  double *WaveMom; //Wavefunction in momentum space
  double rStart; //radius to start integration (fm)
  double rStop; // radius to stop interation (fm) 

  double deltaR; // radial spacing between points in Wavefunction arrays
  int mWave; //array length for interior wavefunction
  int nWave; //array legth for total wavefunction
  int lMax;
  int llMax;
  void list();

  double ** LogDerArray;
  int Nl; 
  int LogDerMin;
  int LogDerMax;
  int NLogDer;
  double LogDerDifference(double,int);

  double * SigmaAb; //absorption cross section for each j
  double * SigmaAb_pos; // for postive parity
  double * SigmaAb_neg; // for negative parity
  double absorbAll;
  double absorbSurface;
  double absorbVolume;
  double absorbStandard;
  double absorbSO; //!< contribution of absorption from imag spin-orbit pot
  double Sreal; //!< real component of Eikonal Smatrix 
  double Simag; //!< imaginal component of Eikonal Smatrix

};

#endif
