#include "whit.h"
#include <iostream>
using namespace std;

int main ()
{
  whit Whit(16);
  double x;
  //  cin >> x;
  for (int i=-400;i<=400;i++)
    {
      x = (double)i/100.;
     double U= Whit.gamma2(x);
     cout << x << " " << U << endl;
    }
  return 1;
}
