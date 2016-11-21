#include "reaction.h"


/**
 *\brief total complex optical-model potential
 *
 *calculates the complex optical-model potential : includes
 * the surface, volume, spin-orbit and Hartree-Fock compoents
 * with appropriate dispersive corrections 
 */

class pot
{

 public:
  reaction* React[1];
  int Nreact;

  static int const TotPara=103;
  int varied[TotPara];
  double scaling[TotPara];
  double allPara[TotPara];
  int squared[TotPara];  

  pot();
  pot(int);
  ~pot();
  void decodePara();
  void justPara();
  void initialForNewEcm(double);
  void GetPot(double,int,double);
  void potential(double,int,double,double);

  double EcmOld;
  double Real;
  double Imag;
};
