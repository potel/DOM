#include "gauss16.h"
#include <iostream>
#include "Gauss16.h"

using namespace std;

int main()
{
  gauss16 G16;

  for (int i=0;i<16;i++) cout << G16.x(i) << " " << G16.w(i) << endl;


  GaussInteg g16(Gauss16);

  for (int i=0;i<16;i++) cout << g16.x(i) << " " << g16.w(i) << endl;
}
