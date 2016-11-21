#include "twoFermi.h"

class twoFermiBi
{
 protected:
  twoFermi TwoFermi;
  double A_above;
  double A_below;
  double B;
  double C;
  double D;
  double Wstart;
  double Ef;
  double derivative;

 public:
  twoFermiBi(){};
  twoFermiBi(double A_above, double A_below,double B, double C, double D, 
            double Wstart,double Ef); 
  void init(double A_above, double A_below, double B, double C, double D, 
            double Wstart, double Ef);
  double funct(double E);
  double derFunct(double E);
  double deltaV(double E);
  double derDeltaV(double E);
}; 			
