#include "twoGauss.h"



/**
 * constructor
 */
twoGauss::twoGauss(double A0, double B0, double E00, double sigma0, double m0, 
		   double Ef0) : disperse()
{
  init(A0,B0,E00,sigma0,m0,Ef0);
}
//*****************************************************************
/**
 * initialization of the the class if one wants to change the parameters
 */
void twoGauss::init(double A0, double B0, double E00, double sigma0, double m0, 
double Ef0)
{
  A = A0;
  B = B0;
  E0 = E00;
  sigma = sigma0;
  m = m0;
  Ef = Ef0;

  x0 = E0 - Ef;
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param E is the center-of-mass energy on the nucleon in MeV
  */
double twoGauss::funct(double E)
{
  double x = E - Ef;
  return functX(x);
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param x = Ecm-Efermi  in MeV
  */
double twoGauss::functX(double x)
{
  return A/2./sqrt(2.*pi)/sigma*pow(abs(x),m)/(pow(abs(x),m)+pow(B,m)) *
    (exp(-pow((x-x0)/sigma,2)/2.)+ exp(-pow((x+x0)/sigma,2)/2.));
}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param E in the center-of-mass energy of the nucleon in MeV
  */
double twoGauss::derFunct(double E)
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
double twoGauss::derFunctX(double x)
{
  double denom = pow(abs(x),m)+pow(B,m);
  double one = pow(abs(x),m)/denom;
  double derOne;
  if (x == 0.) derOne = 0.;
  else derOne = m/x*(pow(abs(x),m)/denom
		  - pow(abs(x),2*m)/pow(denom,2));
 

  double two = exp(-pow((x-x0)/sigma,2)/2.);
  double derTwo = -(x-x0)/pow(sigma,2)*two;
  double three = exp(-pow((x+x0)/sigma,2)/2.);
  double derThree = -(x+x0)/pow(sigma,2)*three;
  double four = two + three;
  double derFour = derTwo + derThree;

  return A/2./sqrt(2.*pi)/sigma*(derOne*four + one*derFour);
}
//****************************************************************
  /**
   * returns the magnitude of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double twoGauss::deltaV(double E)
{
  return deltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double twoGauss::derDeltaV(double E)
{
  return derDeltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for Hole occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double twoGauss::derDeltaVHole(double E)
{
  return derDeltaVHoleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for Particle occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double twoGauss::derDeltaVParticle(double E)
{
  return derDeltaVParticleX(E-Ef);
}
