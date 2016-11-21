/*
 * =====================================================================================
 *
 *       Filename:  irreducible.cpp
 *
 *    Description: This code shows how to get V_lj(r) 
 *
 *        Version:  1.0
 *        Created:  09/16/2016 02:17:40 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author: Mack Atkinson 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "reaction.h"
#include "surface.h"
#include "pot.h"
#include "potPara.h"
#include <fstream>
#include <iostream> 
using namespace std;

int main(){

   //This constructs the potential object, meaning it reads in the parameter input file
   //Look at "pot.cpp" to see what it does
   pot opt;

   int l = 2;
   double j = 2.5;
   double Ecm = 35;
   double r = 2.0;

   //This sets the spin-orbit, and energy of the potential, as well as where in r 
   opt.potential(r,l,j,Ecm);

   //Extracting the information from the object opt
   double vreal = opt.Real;
   double vimag = opt.Imag;

   cout<<endl<<endl<<"real potential at 2.0 = "<<vreal<<endl;
   cout<<"imag potential at 2.0 = "<<vimag<<endl<<endl<<endl;

   double dr = 0.05;

   ofstream fpot("pot.txt");

   for(int i=0;i<200;i++){
      double r = i*dr;
      //Since the energy has already been set, GetPot is sufficient for different radial points
      opt.GetPot(r,l,j);
      vreal = opt.Real;
      vimag = opt.Imag;

      fpot<<r<<" "<<vreal<<" "<<vimag<<endl;
   }

   fpot.close();


}

