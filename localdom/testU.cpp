#include "whit.h"
#include <iostream>
using namespace std;

int main ()
{
  whit Whit(16);
  double a,b,x;
  cin >> a >>b>>x;
  double U= Whit.hypergeometricU(a,b,x);
  cout << U << endl;
  return 1;
}
