#include "asy.h"
#include <iostream>


using namespace std;
  

int main()
{
  asy Asy;
  Asy.init(1.65,0.0,60.,0.0);

   for (int i=-200;i<200;i++)
   {
     double E = (double)i;
     //double E = -100.;
      cout << E << " " << Asy.dispersive(E) << endl;
       }

  return 0;
}
