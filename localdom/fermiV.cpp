#include "fermiV.h"



/**
 * Constructor
 */
fermiV::fermiV(double A0, double B0, double C0, double Ef0) : disperse()
{
  init(A0,B0,C0,Ef0);
}
//*****************************************************************
/**
 * initialization of the the class if one wants to change the parameters
 */
void fermiV::init(double A0, double B0, double C0, double Ef0)
{
  A = A0;
  B = B0;
  C = C0;
  Ef = Ef0;

}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param E is the center-of-mass energy on the nucleon in MeV
  */
double fermiV::funct(double E)
{

  double x = E - Ef;
  return functX(x);
}
//*************************************************************
  /**
   * returns the magnitude of the surface imaginary potential
   \param x = Ecm-Efermi  in MeV
  */
double fermiV::functX(double x)
{
  double ee = exp((abs(x)-B)/C);
  return A*ee/(1.+ee);

}
//***************************************************************
  /**
   * returns the energy derivative of the magnitude of 
   * the surface imaginary potential
   \param E in the center-of-mass energy of the nucleon in MeV
  */
double fermiV::derFunct(double E)
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
double fermiV::derFunctX(double x)
{

  double ee = exp((abs(x)-B)/C);
  double out = A/C*ee/pow(1.+ee,2);
  if (x < 0.) out *= -1.;
  return out;
}
//****************************************************************
  /**
   * returns the dispersive correction
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double fermiV::deltaV(double E)
{
  return deltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential
   \param E is the center-of-mass energy of the nucleon in MeV
   */
double fermiV::derDeltaV(double E)
{
  return derDeltaVX(E-Ef);
}
//*************************************************************
  /**
   * returns the derivative of the real dispersive correction associated with
   * the surface imaginary potential. The function deltaV must be run first.
   */
double fermiV::derDeltaV()
{
  return derDeltaVX();
}

//*************************************************************
  /**
   * returns the disperive factor used for Hole occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double fermiV::derDeltaVHole(double E)
{
  return derDeltaVHoleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle occupation probabilities
   \param E is the center-of-mass energy in MeV
  */
double fermiV::derDeltaVParticle(double E)
{
  return derDeltaVParticleX(E-Ef);
}
//*************************************************************
  /**
   * returns the disperive factor used for particle (hole)occupation 
   * probabilities when the energy is above (below) the fermi energy.
   * The function deltaV must be run first.
  */
double fermiV::derDeltaVParticleHole()
{
  return derDeltaVParticleHoleX();
}
