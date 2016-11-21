#include "sphericalB.h"

int main()
{
  sphericalB sp;
  int l;
  double rho;

    {
      cout << " l, rho" << endl;
      cin >> l >> rho;
      cout << sp.k(l,rho) << endl;
      cout << sp.derivative << endl;
      cout << sp.LogDer_k(l,rho) << endl;
    }
}
