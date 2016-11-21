#include "whit.h"
#include <iostream>
using namespace std;

int main ()
{
  whit Whit(16);
  double gamma = -30.2942;
  double Kwave = .934716;
  double r = 12.;
  int l= 2;

  for (int i=0;i<10;i++)
    {
      double out1 = Whit.AsymptoticExpansion(-gamma,l,2.*r*Kwave);
      double out2 = Whit.whittackerW(gamma,l,r*Kwave);

      cout << r << " " << out1 << " " << out2 << endl; 
      r += 1.;
    }
  //  cout << Whit.oF1(2.,16.248) << endl;
  return 1;
}
