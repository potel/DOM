#include "locality.h"
#include <iostream>


using namespace std;
int main()
{
  locality L(1.);

  double gamma = .2004;
  L.load(1.05,gamma,-101.-4.,6.18);


  double dd = L.getV(-9.5) - L.getV(-10.5);
  double value = L.getV(-10.);
  double d = L.getDV();
  cout << d << " " << dd << endl;


  for (int i=-50;i<200;i++)
    {
      double E= (double)i;
      double V = L.getV(E);
      cout << E << " " << -V << endl;
    }
}
