#include "pot.h"


void GetPot(double r, int l, double j, double E, double & Real, double & Imag)
{
  static pot Pot;
  Pot.potential(r,l,j,E);
  Real = Pot.Real;
  Imag = Pot.Imag;
}


int main()
{

  int l = 1;
  double  j = .5;
  double E = 95.;
  double Real= 0.;
  double Imag= 0.;

  for (int i= 0;i<100;i++)
    {
      double r = (double)i/10.;
      GetPot(r,l,j,E,Real,Imag);
      cout << r << " " << Real << " " << Imag << endl;
    }

  return 0;
}
