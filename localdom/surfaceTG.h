#ifndef surfaceTG_
#define surfaceTG_

#include "potPara.h"
#include "twoGauss.h"
#include <iostream>

using namespace std;

/**
 *\brief imaginary surface potential and dispersive correction
 *
 *this class deals with all aspects of the surface imaginary potential, its
 *dispersive correction and contributions to effective mass and occupation
 *probabilities. the magnitude of the imaginary potential is parametrized as
 *\f$ W(E) = \frac{A}{2\sqrt{2\pi}\sigma} \frac{\|X\|^{m}}{\|X\|^{m}+B^m} \left( \exp\left(-\frac{(X-E_0)^2}{2 \sigma^2}\right) + \exp\left(-\frac{(X+E_0)^2}{2 \sigma^2} \right)\right) \f$
 * where \f$ X = E-E_{Fermi} \f$ 
 */


class surfaceTG
{
 public:
  void load(double,double,double,double,double,double,double,double);
  void load(double strength, double R0, double a0);
  void SetEnergy(double);
  double ImaginaryPot(double);
  double DispersiveCorrection(double);
  double DerivativeDispersive(double);
  double DerivativeParticleHole(double);

  double CentralImaginaryPotential(double);
  double CentralDeltaRealPotential(double);
  double CentralDerDeltaRealPotential(double);

  potPara Imaginary;  //imaginary potential 
  potPara Dispersive; //dispersive coorection to real potential
  potPara DerDispersive; //derivative of dispersive need for effective mass
  potPara ParticleHole; // derivative of particle of hole contribution to 
                   // dispersive correction needed for occupation probabilities

  twoGauss TG; // energy dependence of imaginary potential

  double A;
  double B;
  double E0;
  double sigma;
  double Efermi;
  double Ecm;

  double R;
  double a;
  double m;

};

#endif 
