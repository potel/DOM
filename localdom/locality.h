#ifndef locality_
#define locality_


#include <cmath>
#include <iostream>
#include "localityException.h"
using namespace std;


class locality
{
 public:
  locality(double beta, double gamma, double V0,double mu, double coulShift);
  locality(double mu=1.);
  void load(double beta, double gamma, double V0, double coulShift);
  double getV(double E);
  double getDV();


  private:
  double beta;
  double gamma;
  double V0;
  double coulShift;
  double mu; //!< reduced mass
  double beta2;
  double gamma4;
  double dV;
};
#endif
