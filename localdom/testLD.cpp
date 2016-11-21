#include "levelD.h"



int main()
{
  string name("ca40");
  levelD ld(name);

  double Ex;
  int J;
  for (;;)
    {
      cout << "Ex,J" << endl;
      cin >> Ex >> J;
      cout << ld.getLD(Ex,J) << endl; 
    }
}
