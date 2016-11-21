#include "fermiGap.h"



/**
 * Constructor
 */
fermiGap::fermiGap(double A0, double B0, double C0, double D0, double m0,
    double Wstart0,double Ef0) : disperse()
{
  init(A0,B0,C0,D0,m0,Wstart0,Ef0);
}
//*****************************************************************
/**
 * initialization of the the class if one wants to change the parameters
 */
void fermiGap::init(double A0, double B0, double C0, double D0, double m0,
     double Wstart0, double Ef0)
{
  A = A0;
  B = B0;
  C = C0;
  D = D0;
  m = m0;
  Wstart = Wstart0;
  Ef = Ef0;

}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param E is the center-of-mass energy on the nucleon in MeV
  */
double fermiGap::funct(double E)
{
  double x = E - Ef;
  return functX(x);
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param x = Ecm-Efermi  in MeV
  */
double fermiGap::functX(double x)
{
  double fact = abs(x) - Wstart;
  if (fact <= 0.) return 0.;
  

  return A/(1.+exp((x-C)/D))/(1.+exp(-(x+C)/D)) *
   pow(fact,m)/(pow(fact,m)+pow(B,m));
}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param E in the center-of-mass energy of the nucleon in MeV
  */
double fermiGap::derFunct(double E)
{
  double x = E - Ef;
  return derFunctX(x);
}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param x = Ecm-Efermi in MeV
  */
double fermiGap::derFunctX(double x)
{
  double fact = abs(x) - Wstart;
  if (fact <= 0.) return 0.;
  double denom = pow(fact,m)+pow(B,m);
  double one = pow(fact,m)/denom;
  double derOne = m*pow(fact,m-1.)*(1.- one)/denom;
  if (x < 0.) derOne *= -1.;


  double two = 1.+exp((x-C)/D);
  double derTwo = exp((x-C)/D)/D;
  double three = 1/two;
  double derThree = -1./pow(two,2)*derTwo;
  double four = 1.+exp(-(x+C)/D);
  double derFour = -exp(-(x+C)/D)/D;
  double five = 1./four;
  double derFive = -1./pow(four,2)*derFour;

  return A*(derOne*three*five+one*derThree*five+one*three*derFive);
}
//****************************************************************
  /**
   * returns the dispersive correction
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double fermiGap::deltaV(double E)
{
  return deltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double fermiGap::derDeltaV(double E)
{
  return derDeltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential. The function deltaV must be run first.
   */
double fermiGap::derDeltaV()
{
  return derDeltaVX();
}

//*************************************************************
  /**
   * returns the disperive factor used for Hole occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double fermiGap::derDeltaVHole(double E)
{
  return derDeltaVHoleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double fermiGap::derDeltaVParticle(double E)
{
  return derDeltaVParticleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle (hole)occupation 
   * probabilities when the energy is above (below) the fermi energy.
   * The function deltaV must be run first.
  */
double fermiGap::derDeltaVParticleHole()
{
  return derDeltaVParticleHoleX();
}
