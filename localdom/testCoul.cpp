#include "coul.h"
#include <iostream>

using namespace std;

int main()

{
  coul Coul;
  double x = 10.;
  Coul.init(0,10.,x);
  cout << Coul.F << " " << Coul.G << " " << Coul.dF << " " << Coul.dG << endl;
  cout << Coul.dF*Coul.G - Coul.F*Coul.dG << endl;
  cout << x/(pow(Coul.F,2)+pow(Coul.G,2)) << endl;
}
