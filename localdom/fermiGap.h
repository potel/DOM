#include "disperse.h"


/**
 *!\brief parameterization of the surface with a fermi function
 *
 * A parametrization of the energy dependence of the surface imaginary 
 * potential and also gives the dispersive corrections, etc , associated
 * with this.
 * if \f$ |x| - W_{start} < 0 \f$ then \f$W(E) = 0\f$ else
 *\f$ W(E) = A \frac{(|X|-W_{start})^{m}}{(|X|-W_{start})^{m}+B^m} \frac{1}{1+exp\frac{X-E_0}{\delta_E}} \frac{1}{1+exp-\frac{X+E_0}{\delta_E}} \f$
 * where \f$ X = E-E_{Fermi} \f$ 
 */
 
class fermiGap : public disperse
{
 protected:
  double A;
  double B;
  double C;
  double D;
  double m;
  double Ef; //!< Fermi energy in MeV
  double Wstart;

 public:
  fermiGap(){};
  fermiGap(double A, double B, double C, double D, double m, 
            double Wstart,double Ef); 
  void init(double A, double B, double C, double D, double m, 
            double Wstart, double Ef); 

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
