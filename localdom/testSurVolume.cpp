#include "surVolume.h"
#include <iostream>
using namespace std;

int main()
{
  surVolume sur;
  double Efermi = 0.;//-4.71;
  double alpha = 0.;//.7063;
  double Ea = 60.;
  sur.init(13.8397,42.8872,36.3939,4.,Efermi,Ea,alpha,11.);

  /*
  cout << sur.funct(Efermi) << " " << sur.deltaV(Efermi) << " " << sur.derDeltaV(Efermi) << endl; 
  cout << sur.funct(Efermi+.001) << " " << sur.deltaV(Efermi+.001) << " " << sur.derDeltaV(Efermi+.001) << endl; 
  cout << sur.funct(Efermi-.001) << " " << sur.deltaV(Efermi-.001) << " " << sur.derDeltaV(Efermi-.001) << endl; 



  cout << (sur.funct(Efermi+.01) - sur.funct(Efermi-.01))/.002 << " " <<
    sur.derFunct(Efermi) << " " << sur.derFunct(Efermi+.001) << endl;

  return 0;
  */


  for (int i=0;i<200;i++)
    {
      double e = (double)i-100.;
            cout << e << " " << (sur.funct(e+.01)-sur.funct(e-.01))/.02 << 
      	" " << sur.derFunct(e) << endl;

	    /*
      cout << e << " " << (sur.deltaV(e+.01)-sur.deltaV(e-.01))/.02 << 
	" " << sur.derDeltaV(e) << endl;
	    */
      /*
      sur.deltaV(e);
      cout << e << " " << sur.funct(e) << " " << sur.funct(e) << " " << sur.deltaV(e) << " " << sur.derDeltaV(e) << endl;
      */
    }
}
