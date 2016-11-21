#include "twoFermiBi.h"

twoFermiBi::twoFermiBi(double A_above0, double A_below0, double B0, double C0,
		       double D0, double Wstart0, double Ef0)
{
  init(A_above0,A_below0,B0,C0,D0,Wstart0,Ef0);
}
//****************************************************************
void twoFermiBi::init(double A_above0, double A_below0, double B0, double C0,
		       double D0, double Wstart0, double Ef0)
{
  A_above = A_above0;
  A_below = A_below0;
  B = B0;
  C = C0;
  D = D0;
  Wstart = Wstart0;
  Ef = Ef0;

  TwoFermi.init(1.,B,C,D,Wstart,Ef);
}
//******************************************************************
double twoFermiBi::funct(double E)
{
  if (E > Ef) return A_above*TwoFermi.funct(E);
  else return A_below*TwoFermi.funct(E);
}
//******************************************************************
double twoFermiBi::derFunct(double E)
{
  if (E > Ef) return A_above*TwoFermi.derFunct(E);
  else return A_below*TwoFermi.derFunct(E);
}
//******************************************************************
double twoFermiBi::deltaV(double E)
{
  TwoFermi.deltaVAboveBelow(E);
  derivative = A_above*TwoFermi.derDeltaVAbove() 
           + A_below*TwoFermi.derDeltaVBelow();
  return A_above*TwoFermi.deltaVAbove() + A_below*TwoFermi.deltaVBelow();
}
//******************************************************************
double twoFermiBi::derDeltaV(double E)
{
  return derivative;
}
