#include "fermiHole.h"



/**
 * Constructor
 */
fermiHole::fermiHole(double A0, double B0, double C0, double D0, double m0, 
		   double Ef0) : disperse()
{
  init(A0,B0,C0,D0,m0,Ef0);
}
//*****************************************************************
/**
 * initialization of the the class if one wants to change the parameters
 */
void fermiHole::init(double A0, double B0, double C0, double D0, double m0, 
double Ef0)
{
  A = A0;
  B = B0;
  C = C0;
  D = D0;
  m = m0;
  Ef = Ef0;

}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param E is the center-of-mass energy on the nucleon in MeV
  */
double fermiHole::funct(double E)
{
  double x = E - Ef;
  return functX(x);
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param x = Ecm-Efermi  in MeV
  */
double fermiHole::functX(double x)
{
  return A/(1.+exp((x-C)/D))/(1.+exp(-(x+C)/D))*
   pow(abs(x),m)/(pow(abs(x),m)+pow(B,m));
}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param E in the center-of-mass energy of the nucleon in MeV
  */
double fermiHole::derFunct(double E)
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
double fermiHole::derFunctX(double x)
{

  double denom = pow(abs(x),m)+pow(B,m);
  double one = pow(abs(x),m)/denom;
  double derOne;
  if (x == 0.) derOne = 0.;
  else derOne = m/x*(pow(abs(x),m)/denom
		  - pow(abs(x),2*m)/pow(denom,2));


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
double fermiHole::deltaV(double E)
{
  return deltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double fermiHole::derDeltaV(double E)
{
  return derDeltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential. The function deltaV must be run first.
   */
double fermiHole::derDeltaV()
{
  return derDeltaVX();
}

//*************************************************************
  /**
   * returns the disperive factor used for Hole occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double fermiHole::derDeltaVHole(double E)
{
  return derDeltaVHoleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double fermiHole::derDeltaVParticle(double E)
{
  return derDeltaVParticleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle (hole)occupation 
   * probabilities when the energy is above (below) the fermi energy.
   * The function deltaV must be run first.
  */
double fermiHole::derDeltaVParticleHole()
{
  return derDeltaVParticleHoleX();
}
