#include "channel.h"
#include <string>


int main()
{
  string name("co58");
  string tl("cu58Neutron");
  channel chan(&name,&tl,0.);

  double ld;
  for (int i=1;i<10;i++)
    {
      ld = 0.;
      double E = (double)i + .5;
      for (int j=0;j<20;j++)
        {
          double ldd = chan.FermiGas(E,j);
          if (ldd == 0.) break;
          ld += ldd;
        }
      cout << E << " " << ld << endl;
    }

}
