#include "locality.h"




locality::locality(double beta0, double gamma0, double V00, double mu0, 
          double coulShift)
{
  mu = mu0;
  load(beta0,gamma0,V00,coulShift);
}
//************************************************************************
locality::locality(double mu0/*=1.*/)
{
  mu = mu0;
}
//************************************************************************
void locality::load(double beta0, double gamma0, double V00, double coulShift0)
{
  beta = beta0;
  gamma = gamma0;
  V0 = V00;
  coulShift = coulShift0;
  beta2 = pow(beta,2)*mu*0.0239/2.;
  gamma4 = pow(gamma,4)*pow(mu,2)*pow(0.0239,2)*4.;
}
//*************************************************************************
double locality::getV(double E)
{
  //initial guess 
  double V = V0*(1-beta2*(E-coulShift))/(1.-beta2*V0);
    //iterate
  int tries = 0;
  for(;;)
    {
      double term1 = V0*exp(-beta2*(E-V-coulShift)+
           gamma4*pow(E-V-coulShift,2));
      double delta = V - term1;
      double ddelta = 1. - term1*(beta2 - 2.*gamma4*(E-V-coulShift));
      double dv = -delta/ddelta;
      if (fabs(dv) < .01) break;
      V += dv;
      tries++;
      if (tries > 100) 
	{
          //solution no possible
          cout << "V0= " << V0 << " beta= "  << beta << " gamma " << gamma 
	       << " coulShift= " << coulShift << endl;
          cout << " beta2 = " << beta2 << " gamm24= " << gamma4 << endl;

	  localityException  locEx;
          throw locEx; 
	}
    }

  double fact = V*(-beta2 + 2.*gamma4*(E-V-coulShift));
  dV = fact/(1.+fact);

  return V;
}
//**********************************************************************
double locality::getDV()
{
  return dV;
}
