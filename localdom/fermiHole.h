#include "disperse.h"


/**
 *!\brief parameterization of the surface with a fermi function
 *
 * A parametrization of the energy dependence of the surface imaginary 
 * potential and also gives the dispersive corrections, etc , associated
 * with this
 *\f$ W(E) = A \frac{|X|^{m}}{|X|^{m}+B^m} \frac{1}{1+exp\frac{X-E_0}{\delta_E}} \frac{1}{1+exp-\frac{X+E_0}{\delta_E}} \f$
 * where \f$ X = E-E_{Fermi} \f$ 
 */
 
class fermiHole : public disperse
{
 protected:
  double A;
  double B;
  double C;
  double D;
  double Ef; //!< Fermi energy in MeV
  double m;

 public:
  fermiHole(){};
  fermiHole(double A, double B, double C, double D, double m, double Ef); 
  void init(double A, double B, double C, double D, double m, double Ef); 

  double funct(double E);
  double derFunct(double E);

  double functX(double x);
  double derFunctX(double x);

  double deltaV(double E);
  double derDeltaV(double E);

  double derDeltaVHole(double E);
  double derDeltaVParticle(double E);

  double derDeltaVParticleHole();
  double derDeltaV();

};
