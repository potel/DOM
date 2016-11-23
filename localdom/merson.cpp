using namespace std;

#include "merson.h"
#include <sstream>

//**************************************************************************
// merson initialization
void merson::initMerson(double acc0,double step_min0, double step00,int jtest0)
{
   acc = acc0;
   step_min = step_min0;
   step0 = step00;
   jtest = jtest0;
   ok = 1;
}
//**************************************************************************
// solve system of coupled ordinary differential equations
void merson::solveMerson(valarray<double> *y, double xstart,
      double xend)
{
   int n = y->size();
   valarray<double> yz(n);
   valarray<double> a(n);
   valarray<double> b(n);
   valarray<double> f(n);

   // rzero is a number with a magnitude of order equal to the noise level of
   //the machine i.e. in the range of the rounding errors
   double const rzero =1.0e-23;

   //jtest = test parameter related to the steplength in the following way.
   //jtest = 0, if during the calculation we get ABS(H) less than hmin (by
   //repeated halving of the steplength), than an error message is printed,
   //ok is set equal to .FALSE. followed by return to the calling program
   //jtest = 1, checking as for jtest=0, but the calculation will continue
   //with the fixed dtep length hmin


   ok = 1;

   // store internally parameters in list

   double x = xstart;
   double step = step0;


   int leave = 0;
   for (;;) 
   {
      double hsv = step;
      double cof = xend - x;
      if (abs(step) >= abs(cof)) 
      {
         //step length is greater than remaining interval to be integrated
         step = cof;
         //if (abs(cof/hsv) < rzero) *y = w;
         leave = 1;
         //if leave is true, then step is equal to the maximum possible
         //steplength within the remaining part of the domain of
         //integration.
      }

      yz = *y;
      double ht = step/3.;


      f = diff(x,*y);
      a = ht*f;

      *y = a + yz;

      x = x + ht;
      f = diff(x,*y);

      a = 0.5*a;

      *y = 0.5*ht*f + a + yz;

      f = diff(x,*y);

      b = 4.5*ht*f;

      *y = 0.25*b + 0.75*a + yz;

      x = x + 0.5*ht;
      f = diff(x,*y);

      a = 2.*ht*f + a;

      *y = 3.*a - b + yz;

      x = x + 0.5*step;
      f = diff(x,*y);


      int accuracy = 1;

      int k = 0;
      while (accuracy && k < n)
      {
         b[k] = -0.5*ht*f[k] - b[k] + 2.*a[k];
         (*y)[k] = (*y)[k] - b[k];
         a[k] = abs(5.*acc* (*y)[k]);
         b[k] = abs(b[k]);
         if (abs((*y)[k]) > rzero && b[k] > a[k]) accuracy = 0;
         k++;
      }


      //if accuracy is false, the required accuracy for all computed values
      //was not obtained
      if (accuracy==0)
      {
         //halve step length      
         cof = 0.5*step; 
         if(abs(cof) >= step_min)
         {
            *y = yz;
            x = x - step;
            step = cof;
            leave = 0;
         }
         else if (jtest == 0) 
         {
            cout  << "**Merson error***" << endl;
            cout << "jtest= " <<jtest <<" step_min= " << step_min 
               << " x= " << x << endl;
            ok = 0;
            abort();
         }
         else
         {
            //continue with constant step length equal to hmin
            step = step_min;
            if (hsv < 0.) step = -step;
         }
      }
      else
      {
         //required accuray obtained
         //test if step length doubling is possible
         valarray<double> g(b-0.03125*a);
         if (g.max() <= 0.0) step = 2.0*step;

      }
      if (leave) break;
   }

   // calculation has finished - 

}


//**************************************************************************
// solve system of coupled ordinary differential equations
void merson::solveMersonk(valarray<double> *y, double xstart,
      double xend, double Kwave,int l,int up)
{
   int n = y->size();
   valarray<double> yz(n);
   valarray<double> a(n);
   valarray<double> b(n);
   valarray<double> f(n);
   vector<double> wf0;
   vector<double> wf1;
   vector<double> wf2;
   vector<double> wf3;
   vector<double> rmesh;

   wf0.push_back(0);
   wf1.push_back(0);
   wf2.push_back(0);
   wf3.push_back(0);
   rmesh.push_back(0);
   wf0.push_back((*y)[0]);
   wf1.push_back((*y)[1]);
   wf2.push_back((*y)[2]);
   wf3.push_back((*y)[3]);
   rmesh.push_back(xstart);

   // rzero is a number with a magnitude of order equal to the noise level of
   //the machine i.e. in the range of the rounding errors
   double const rzero =1.0e-23;

   //jtest = test parameter related to the steplength in the following way.
   //jtest = 0, if during the calculation we get ABS(H) less than hmin (by
   //repeated halving of the steplength), than an error message is printed,
   //ok is set equal to .FALSE. followed by return to the calling program
   //jtest = 1, checking as for jtest=0, but the calculation will continue
   //with the fixed dtep length hmin


   ok = 1;

   // store internally parameters in list

   double x = xstart;
   double step = step0;


   int leave = 0;
   for (;;) 
   {
      double hsv = step;
      double cof = xend - x;
      if (abs(step) >= abs(cof)) 
      {
         //step length is greater than remaining interval to be integrated
         step = cof;
         //if (abs(cof/hsv) < rzero) *y = w;
         leave = 1;
         //if leave is true, then step is equal to the maximum possible
         //steplength within the remaining part of the domain of
         //integration.
      }

      yz = *y;
      double ht = step/3.;

      f = diff(x,*y);
      a = ht*f;

      //This is y_A in merson method
      *y = a + yz;

      //This is x_A in merson method
      x = x + ht;

      //This is y'_A
      f = diff(x,*y);

      a = 0.5*a;

      //This is y_B
      *y = 0.5*ht*f + a + yz;

      f = diff(x,*y);

      b = 4.5*ht*f;

      //This is y_C
      *y = 0.25*b + 0.75*a + yz;

      //This is x_C (1/3 + 1/6 = 1/2)
      x = x + 0.5*ht;
      //y'_C
      f = diff(x,*y);

      a = 2.*ht*f + a;

      //y_D
      *y = 3.*a - b + yz;

      //x_1 (this should be the next grid point)
      x = x + 0.5*step;
      f = diff(x,*y);


      int accuracy = 1;

      int k = 0;
      while (accuracy && k < n)
      {
         b[k] = -0.5*ht*f[k] - b[k] + 2.*a[k];
         (*y)[k] = (*y)[k] - b[k];
         a[k] = abs(5.*acc* (*y)[k]);
         b[k] = abs(b[k]);
         if (abs((*y)[k]) > rzero && b[k] > a[k]) accuracy = 0;
         k++;
      }

      wf0.push_back((*y)[0]);
      wf1.push_back((*y)[1]);
      wf2.push_back((*y)[2]);
      wf3.push_back((*y)[3]);
      rmesh.push_back(x);


      //if accuracy is false, the required accuracy for all computed values
      //was not obtained
      if (accuracy==0)
      {
         //halve step length      
         cof = 0.5*step; 
         if(abs(cof) >= step_min)
         {
            *y = yz;
            x = x - step;
            step = cof;
            leave = 0;
         }
         else if (jtest == 0) 
         {
            cout  << "**Merson error***" << endl;
            cout << "jtest= " <<jtest <<" step_min= " << step_min 
               << " x= " << x << endl;
            ok = 0;
            abort();
         }
         else
         {
            //continue with constant step length equal to hmin
            step = step_min;
            if (hsv < 0.) step = -step;
         }
      }
      //      else
      //      {
      //         //required accuray obtained
      //         //test if step length doubling is possible
      //         valarray<double> g(b-0.03125*a);
      //         if (g.max() <= 0.0) step = 2.0*step;
      //
      //      }
      //      cout<<(*y)[0]<<endl;
      if (leave) break;
   }

   string llab;
   ostringstream con;
   con << l;
   llab = con.str();

   string jlab;
   ostringstream jon;
   jon << up;
   jlab = jon.str();

   string name = "partials/eep100 " + llab + " " + jlab + ".txt"; 
   ofstream out(name.c_str());
   string nameder = "partials/der100 " + llab + " " + jlab + ".txt"; 
   ofstream outder(nameder.c_str());

   double norm = wf0[wf0.size()-1];
   for(int i=0;i<wf0.size();i++){
      wf0[i] /= norm;
      wf1[i] /= norm;
      wf2[i] /= norm;
      wf3[i] /= norm;
      out<<rmesh[i]<<" "<<wf0[i]<<" "<<wf2[i]<<endl;
      outder<<rmesh[i]<<" "<<wf1[i]<<" "<<wf3[i]<<endl;
   }

   out.close();
   outder.close();

   // calculation has finished - 

}

//void merson::solveMersonkf_(double *yy,int n, double xstart,
      //double xend, double Kwave,int l,int up,double *ff,double *gg)
//{

   ////valarray<double> y(n);
   ////for(int i=0;i<n;i++){
      ////y[i] = yy[i];
   ////}
   ////valarray<double> yz(n);
   ////valarray<double> a(n);
   ////valarray<double> b(n);
   ////valarray<double> f(n);
   ////vector<double> wf0;
   ////vector<double> wf1;
   ////vector<double> wf2;
   ////vector<double> wf3;
   ////vector<double> rmesh;

   ////wf0.push_back(0);
   ////wf1.push_back(0);
   ////wf2.push_back(0);
   ////wf3.push_back(0);
   ////rmesh.push_back(0);
   ////wf0.push_back((*y)[0]);
   ////wf1.push_back((*y)[1]);
   ////wf2.push_back((*y)[2]);
   ////wf3.push_back((*y)[3]);
   ////rmesh.push_back(xstart);

   ////// rzero is a number with a magnitude of order equal to the noise level of
   //////the machine i.e. in the range of the rounding errors
   ////double const rzero =1.0e-23;

   //////jtest = test parameter related to the steplength in the following way.
   //////jtest = 0, if during the calculation we get ABS(H) less than hmin (by
   //////repeated halving of the steplength), than an error message is printed,
   //////ok is set equal to .FALSE. followed by return to the calling program
   //////jtest = 1, checking as for jtest=0, but the calculation will continue
   //////with the fixed dtep length hmin


   ////ok = 1;

   ////// store internally parameters in list

   ////double x = xstart;
   ////double step = step0;


   ////int leave = 0;
   ////for (;;) 
   ////{
      ////double hsv = step;
      ////double cof = xend - x;
      ////if (abs(step) >= abs(cof)) 
      ////{
         //////step length is greater than remaining interval to be integrated
         ////step = cof;
         //////if (abs(cof/hsv) < rzero) *y = w;
         ////leave = 1;
         //////if leave is true, then step is equal to the maximum possible
         //////steplength within the remaining part of the domain of
         //////integration.
      ////}

      ////yz = *y;
      ////double ht = step/3.;

      ////f = diff(x,*y);
      ////a = ht*f;

      //////This is y_A in merson method
      ////*y = a + yz;

      //////This is x_A in merson method
      ////x = x + ht;

      //////This is y'_A
      ////f = diff(x,*y);

      ////a = 0.5*a;

      //////This is y_B
      ////*y = 0.5*ht*f + a + yz;

      ////f = diff(x,*y);

      ////b = 4.5*ht*f;

      //////This is y_C
      ////*y = 0.25*b + 0.75*a + yz;

      //////This is x_C (1/3 + 1/6 = 1/2)
      ////x = x + 0.5*ht;
      //////y'_C
      ////f = diff(x,*y);

      ////a = 2.*ht*f + a;

      //////y_D
      ////*y = 3.*a - b + yz;

      //////x_1 (this should be the next grid point)
      ////x = x + 0.5*step;
      ////f = diff(x,*y);


      ////int accuracy = 1;

      ////int k = 0;
      ////while (accuracy && k < n)
      ////{
         ////b[k] = -0.5*ht*f[k] - b[k] + 2.*a[k];
         ////(*y)[k] = (*y)[k] - b[k];
         ////a[k] = abs(5.*acc* (*y)[k]);
         ////b[k] = abs(b[k]);
         ////if (abs((*y)[k]) > rzero && b[k] > a[k]) accuracy = 0;
         ////k++;
      ////}

      ////wf0.push_back((*y)[0]);
      ////wf1.push_back((*y)[1]);
      ////wf2.push_back((*y)[2]);
      ////wf3.push_back((*y)[3]);
      ////rmesh.push_back(x);


      //////if accuracy is false, the required accuracy for all computed values
      //////was not obtained
      ////if (accuracy==0)
      ////{
         //////halve step length      
         ////cof = 0.5*step; 
         ////if(abs(cof) >= step_min)
         ////{
            ////*y = yz;
            ////x = x - step;
            ////step = cof;
            ////leave = 0;
         ////}
         ////else if (jtest == 0) 
         ////{
            ////cout  << "**Merson error***" << endl;
            ////cout << "jtest= " <<jtest <<" step_min= " << step_min 
               ////<< " x= " << x << endl;
            ////ok = 0;
            ////abort();
         ////}
         ////else
         ////{
            //////continue with constant step length equal to hmin
            ////step = step_min;
            ////if (hsv < 0.) step = -step;
         ////}
      ////}
      //////      else
      //////      {
      //////         //required accuray obtained
      //////         //test if step length doubling is possible
      //////         valarray<double> g(b-0.03125*a);
      //////         if (g.max() <= 0.0) step = 2.0*step;
      //////
      //////      }
      //////      cout<<(*y)[0]<<endl;
      ////if (leave) break;
   ////}

   ////string llab;
   ////ostringstream con;
   ////con << l;
   ////llab = con.str();

   ////string jlab;
   ////ostringstream jon;
   ////jon << up;
   ////jlab = jon.str();

   ////string name = "partials/eep100 " + llab + " " + jlab + ".txt"; 
   ////ofstream out(name.c_str());
   ////string nameder = "partials/der100 " + llab + " " + jlab + ".txt"; 
   ////ofstream outder(nameder.c_str());

   ////double norm = wf0[wf0.size()-1];
   ////for(int i=0;i<wf0.size();i++){
      ////wf0[i] /= norm;
      ////wf1[i] /= norm;
      ////wf2[i] /= norm;
      ////wf3[i] /= norm;
      ////out<<rmesh[i]<<" "<<wf0[i]<<" "<<wf2[i]<<endl;
      ////outder<<rmesh[i]<<" "<<wf1[i]<<" "<<wf3[i]<<endl;
   ////}

   ////out.close();
   ////outder.close();

   ////// calculation has finished - 

////}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  merson
 *      Method:  merson :: solveMerson_eep
 * Description:  Uses merson to solve using the schwandt potential
 *--------------------------------------------------------------------------------------
 */
void merson::solveMerson_eep(valarray<double> *y, double xstart,
      double xend, double Kwave,int l,int up)
{
   int n = y->size();
   valarray<double> yz(n);
   valarray<double> a(n);
   valarray<double> b(n);
   valarray<double> f(n);
   vector<double> wf0;
   vector<double> wf1;
   vector<double> wf2;
   vector<double> wf3;
   vector<double> rmesh;

   wf0.push_back(0);
   wf1.push_back(0);
   wf2.push_back(0);
   wf3.push_back(0);
   rmesh.push_back(0);
   wf0.push_back((*y)[0]);
   wf1.push_back((*y)[1]);
   wf2.push_back((*y)[2]);
   wf3.push_back((*y)[3]);
   rmesh.push_back(xstart);

   // rzero is a number with a magnitude of order equal to the noise level of
   //the machine i.e. in the range of the rounding errors
   double const rzero =1.0e-23;

   //jtest = test parameter related to the steplength in the following way.
   //jtest = 0, if during the calculation we get ABS(H) less than hmin (by
   //repeated halving of the steplength), than an error message is printed,
   //ok is set equal to .FALSE. followed by return to the calling program
   //jtest = 1, checking as for jtest=0, but the calculation will continue
   //with the fixed dtep length hmin


   ok = 1;

   // store internally parameters in list

   double x = xstart;
   double step = step0;

   int leave = 0;
   for (;;) 
   {
      double hsv = step;
      double cof = xend - x;
      if (abs(step) >= abs(cof)) 
      {
         //step length is greater than remaining interval to be integrated
         step = cof;
         //if (abs(cof/hsv) < rzero) *y = w;
         leave = 1;
         //if leave is true, then step is equal to the maximum possible
         //steplength within the remaining part of the domain of
         //integration.
      }

      yz = *y;
      double ht = step/3.;

      f = diffeep(x,*y);
      a = ht*f;

      //This is y_A in merson method
      *y = a + yz;

      //This is x_A in merson method
      x = x + ht;

      //This is y'_A
      f = diffeep(x,*y);

      a = 0.5*a;

      //This is y_B
      *y = 0.5*ht*f + a + yz;

      f = diffeep(x,*y);

      b = 4.5*ht*f;

      //This is y_C
      *y = 0.25*b + 0.75*a + yz;

      //This is x_C (1/3 + 1/6 = 1/2)
      x = x + 0.5*ht;
      //y'_C
      f = diffeep(x,*y);

      a = 2.*ht*f + a;

      //y_D
      *y = 3.*a - b + yz;

      //x_1 (this should be the next grid point)
      x = x + 0.5*step;
      f = diffeep(x,*y);


      int accuracy = 1;

      int k = 0;
      while (accuracy && k < n)
      {
         b[k] = -0.5*ht*f[k] - b[k] + 2.*a[k];
         (*y)[k] = (*y)[k] - b[k];
         a[k] = abs(5.*acc* (*y)[k]);
         b[k] = abs(b[k]);
         if (abs((*y)[k]) > rzero && b[k] > a[k]) accuracy = 0;
         k++;
      }

      wf0.push_back((*y)[0]);
      wf1.push_back((*y)[1]);
      wf2.push_back((*y)[2]);
      wf3.push_back((*y)[3]);
      rmesh.push_back(x);


      //if accuracy is false, the required accuracy for all computed values
      //was not obtained
      if (accuracy==0)
      {
         //halve step length      
         cof = 0.5*step; 
         if(abs(cof) >= step_min)
         {
            *y = yz;
            x = x - step;
            step = cof;
            leave = 0;
         }
         else if (jtest == 0) 
         {
            cout  << "**Merson error***" << endl;
            cout << "jtest= " <<jtest <<" step_min= " << step_min 
               << " x= " << x << endl;
            ok = 0;
            abort();
         }
         else
         {
            //continue with constant step length equal to hmin
            step = step_min;
            if (hsv < 0.) step = -step;
         }
      }
      //      else
      //      {
      //         //required accuray obtained
      //         //test if step length doubling is possible
      //         valarray<double> g(b-0.03125*a);
      //         if (g.max() <= 0.0) step = 2.0*step;
      //
      //      }
      //      cout<<(*y)[0]<<endl;
      if (leave) break;
   }

   string llab;
   ostringstream con;
   con << l;
   llab = con.str();

   string jlab;
   ostringstream jon;
   jon << up;
   jlab = jon.str();

   string name = "partials/optso/eep100 " + llab + " " + jlab + ".txt"; 
   ofstream out(name.c_str());
   string nameder = "partials/optso/der100 " + llab + " " + jlab + ".txt"; 
   ofstream outder(nameder.c_str());

   double norm = wf0[wf0.size()-1];
   for(int i=0;i<wf0.size();i++){
      wf0[i] /= norm;
      wf1[i] /= norm;
      wf2[i] /= norm;
      wf3[i] /= norm;
      out<<rmesh[i]<<" "<<wf0[i]<<" "<<wf2[i]<<endl;
      outder<<rmesh[i]<<" "<<wf1[i]<<" "<<wf3[i]<<endl;
   }

   out.close();
   outder.close();

   // calculation has finished - 

}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  merson
 *      Method:  merson :: solveMerson_schwandt
 * Description:  Uses merson to solve using my schwandt potential
 *--------------------------------------------------------------------------------------
 */
void merson::solveMerson_sch(valarray<double> *y, double xstart,
      double xend, double Kwave,int l,int up)
{
   int n = y->size();
   valarray<double> yz(n);
   valarray<double> a(n);
   valarray<double> b(n);
   valarray<double> f(n);
   vector<double> wf0;
   vector<double> wf1;
   vector<double> wf2;
   vector<double> wf3;
   vector<double> rmesh;

   wf0.push_back(0);
   wf1.push_back(0);
   wf2.push_back(0);
   wf3.push_back(0);
   rmesh.push_back(0);
   wf0.push_back((*y)[0]);
   wf1.push_back((*y)[1]);
   wf2.push_back((*y)[2]);
   wf3.push_back((*y)[3]);
   rmesh.push_back(xstart);

   // rzero is a number with a magnitude of order equal to the noise level of
   //the machine i.e. in the range of the rounding errors
   double const rzero =1.0e-23;

   //jtest = test parameter related to the steplength in the following way.
   //jtest = 0, if during the calculation we get ABS(H) less than hmin (by
   //repeated halving of the steplength), than an error message is printed,
   //ok is set equal to .FALSE. followed by return to the calling program
   //jtest = 1, checking as for jtest=0, but the calculation will continue
   //with the fixed dtep length hmin


   ok = 1;

   // store internally parameters in list

   double x = xstart;
   double step = step0;

   int leave = 0;
   for (;;) 
   {
      double hsv = step;
      double cof = xend - x;
      if (abs(step) >= abs(cof)) 
      {
         //step length is greater than remaining interval to be integrated
         step = cof;
         //if (abs(cof/hsv) < rzero) *y = w;
         leave = 1;
         //if leave is true, then step is equal to the maximum possible
         //steplength within the remaining part of the domain of
         //integration.
      }

      yz = *y;
      double ht = step/3.;

      f = diffsch(x,*y);
      a = ht*f;

      //This is y_A in merson method
      *y = a + yz;

      //This is x_A in merson method
      x = x + ht;

      //This is y'_A
      f = diffsch(x,*y);

      a = 0.5*a;

      //This is y_B
      *y = 0.5*ht*f + a + yz;

      f = diffsch(x,*y);

      b = 4.5*ht*f;

      //This is y_C
      *y = 0.25*b + 0.75*a + yz;

      //This is x_C (1/3 + 1/6 = 1/2)
      x = x + 0.5*ht;
      //y'_C
      f = diffsch(x,*y);

      a = 2.*ht*f + a;

      //y_D
      *y = 3.*a - b + yz;

      //x_1 (this should be the next grid point)
      x = x + 0.5*step;
      f = diffsch(x,*y);


      int accuracy = 1;

      int k = 0;
      while (accuracy && k < n)
      {
         b[k] = -0.5*ht*f[k] - b[k] + 2.*a[k];
         (*y)[k] = (*y)[k] - b[k];
         a[k] = abs(5.*acc* (*y)[k]);
         b[k] = abs(b[k]);
         if (abs((*y)[k]) > rzero && b[k] > a[k]) accuracy = 0;
         k++;
      }

      wf0.push_back((*y)[0]);
      wf1.push_back((*y)[1]);
      wf2.push_back((*y)[2]);
      wf3.push_back((*y)[3]);
      rmesh.push_back(x);


      //if accuracy is false, the required accuracy for all computed values
      //was not obtained
      if (accuracy==0)
      {
         //halve step length      
         cof = 0.5*step; 
         if(abs(cof) >= step_min)
         {
            *y = yz;
            x = x - step;
            step = cof;
            leave = 0;
         }
         else if (jtest == 0) 
         {
            cout  << "**Merson error***" << endl;
            cout << "jtest= " <<jtest <<" step_min= " << step_min 
               << " x= " << x << endl;
            ok = 0;
            abort();
         }
         else
         {
            //continue with constant step length equal to hmin
            step = step_min;
            if (hsv < 0.) step = -step;
         }
      }
      //      else
      //      {
      //         //required accuray obtained
      //         //test if step length doubling is possible
      //         valarray<double> g(b-0.03125*a);
      //         if (g.max() <= 0.0) step = 2.0*step;
      //
      //      }
      //      cout<<(*y)[0]<<endl;
      if (leave) break;
   }

   string llab;
   ostringstream con;
   con << l;
   llab = con.str();

   string jlab;
   ostringstream jon;
   jon << up;
   jlab = jon.str();

   string name = "partials/optso/eep100 " + llab + " " + jlab + ".txt"; 
   ofstream out(name.c_str());
   string nameder = "partials/optso/der100 " + llab + " " + jlab + ".txt"; 
   ofstream outder(nameder.c_str());

   double norm = wf0[wf0.size()-1];
   for(int i=0;i<wf0.size();i++){
      wf0[i] /= norm;
      wf1[i] /= norm;
      wf2[i] /= norm;
      wf3[i] /= norm;
      out<<rmesh[i]<<" "<<wf0[i]<<" "<<wf2[i]<<endl;
      outder<<rmesh[i]<<" "<<wf1[i]<<" "<<wf3[i]<<endl;
   }

   out.close();
   outder.close();

   // calculation has finished - 

}


/*
 *--------------------------------------------------------------------------------------
 *       Class:  merson
 *      Method:  merson :: solveMerson_eep
 * Description:  Uses merson to solve using the schwandt potential
 *--------------------------------------------------------------------------------------
 */

//void merson::solveMerson_eepf_(double *yy,int n, double xstart,
      //double xend, double Kwave,int l,int up,double *ff,double *gg)
//{
   ////valarray<double> y(n);
   ////for(int i=0;i<n;i++){
      ////y[i] = yy[n];
   ////}
   ////valarray<double> yz(n);
   ////valarray<double> a(n);
   ////valarray<double> b(n);
   ////valarray<double> f(n);
   ////vector<double> wf0;
   ////vector<double> wf1;
   ////vector<double> wf2;
   ////vector<double> wf3;
   ////vector<double> rmesh;

   ////wf0.push_back(0);
   ////wf1.push_back(0);
   ////wf2.push_back(0);
   ////wf3.push_back(0);
   ////rmesh.push_back(0);
   ////wf0.push_back((*y)[0]);
   ////wf1.push_back((*y)[1]);
   ////wf2.push_back((*y)[2]);
   ////wf3.push_back((*y)[3]);
   ////rmesh.push_back(xstart);

   ////// rzero is a number with a magnitude of order equal to the noise level of
   //////the machine i.e. in the range of the rounding errors
   ////double const rzero =1.0e-23;

   //////jtest = test parameter related to the steplength in the following way.
   //////jtest = 0, if during the calculation we get ABS(H) less than hmin (by
   //////repeated halving of the steplength), than an error message is printed,
   //////ok is set equal to .FALSE. followed by return to the calling program
   //////jtest = 1, checking as for jtest=0, but the calculation will continue
   //////with the fixed dtep length hmin


   ////ok = 1;

   ////// store internally parameters in list

   ////double x = xstart;
   ////double step = step0;

   ////int leave = 0;
   ////for (;;) 
   ////{
      ////double hsv = step;
      ////double cof = xend - x;
      ////if (abs(step) >= abs(cof)) 
      ////{
         //////step length is greater than remaining interval to be integrated
         ////step = cof;
         //////if (abs(cof/hsv) < rzero) *y = w;
         ////leave = 1;
         //////if leave is true, then step is equal to the maximum possible
         //////steplength within the remaining part of the domain of
         //////integration.
      ////}

      ////yz = *y;
      ////double ht = step/3.;

      ////f = diffeep(x,*y);
      ////a = ht*f;

      //////This is y_A in merson method
      ////*y = a + yz;

      //////This is x_A in merson method
      ////x = x + ht;

      //////This is y'_A
      ////f = diffeep(x,*y);

      ////a = 0.5*a;

      //////This is y_B
      ////*y = 0.5*ht*f + a + yz;

      ////f = diffeep(x,*y);

      ////b = 4.5*ht*f;

      //////This is y_C
      ////*y = 0.25*b + 0.75*a + yz;

      //////This is x_C (1/3 + 1/6 = 1/2)
      ////x = x + 0.5*ht;
      //////y'_C
      ////f = diffeep(x,*y);

      ////a = 2.*ht*f + a;

      //////y_D
      ////*y = 3.*a - b + yz;

      //////x_1 (this should be the next grid point)
      ////x = x + 0.5*step;
      ////f = diffeep(x,*y);


      ////int accuracy = 1;

      ////int k = 0;
      ////while (accuracy && k < n)
      ////{
         ////b[k] = -0.5*ht*f[k] - b[k] + 2.*a[k];
         ////(*y)[k] = (*y)[k] - b[k];
         ////a[k] = abs(5.*acc* (*y)[k]);
         ////b[k] = abs(b[k]);
         ////if (abs((*y)[k]) > rzero && b[k] > a[k]) accuracy = 0;
         ////k++;
      ////}

      ////wf0.push_back((*y)[0]);
      ////wf1.push_back((*y)[1]);
      ////wf2.push_back((*y)[2]);
      ////wf3.push_back((*y)[3]);
      ////rmesh.push_back(x);


      //////if accuracy is false, the required accuracy for all computed values
      //////was not obtained
      ////if (accuracy==0)
      ////{
         //////halve step length      
         ////cof = 0.5*step; 
         ////if(abs(cof) >= step_min)
         ////{
            ////*y = yz;
            ////x = x - step;
            ////step = cof;
            ////leave = 0;
         ////}
         ////else if (jtest == 0) 
         ////{
            ////cout  << "**Merson error***" << endl;
            ////cout << "jtest= " <<jtest <<" step_min= " << step_min 
               ////<< " x= " << x << endl;
            ////ok = 0;
            ////abort();
         ////}
         ////else
         ////{
            //////continue with constant step length equal to hmin
            ////step = step_min;
            ////if (hsv < 0.) step = -step;
         ////}
      ////}
      //////      else
      //////      {
      //////         //required accuray obtained
      //////         //test if step length doubling is possible
      //////         valarray<double> g(b-0.03125*a);
      //////         if (g.max() <= 0.0) step = 2.0*step;
      //////
      //////      }
      //////      cout<<(*y)[0]<<endl;
      ////if (leave) break;
   ////}

   ////string llab;
   ////ostringstream con;
   ////con << l;
   ////llab = con.str();

   ////string jlab;
   ////ostringstream jon;
   ////jon << up;
   ////jlab = jon.str();

   ////string name = "partials/optso/eep100 " + llab + " " + jlab + ".txt"; 
   ////ofstream out(name.c_str());
   ////string nameder = "partials/optso/der100 " + llab + " " + jlab + ".txt"; 
   ////ofstream outder(nameder.c_str());

   ////double norm = wf0[wf0.size()-1];
   ////for(int i=0;i<wf0.size();i++){
      ////wf0[i] /= norm;
      ////wf1[i] /= norm;
      ////wf2[i] /= norm;
      ////wf3[i] /= norm;
      ////out<<rmesh[i]<<" "<<wf0[i]<<" "<<wf2[i]<<endl;
      ////outder<<rmesh[i]<<" "<<wf1[i]<<" "<<wf3[i]<<endl;
   ////}

   ////out.close();
   ////outder.close();

   ////// calculation has finished - 

////}
