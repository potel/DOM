#include "reaction.h"
#include "surface.h"
#include "pot.h"
#include "potPara.h"
#include <fstream>
#include <iostream> 
using namespace std;

int main(){

   pot potential;

   potential.React[0]->scatter->ReadPot();

   potential.React[0]->scatter->SchwandtPot(20,1);


   //double Elower = -15;
   //double Eupper = 1;
   //int l =2;
   //double j = 1.5;
   //int Ifine=1;

   //potential.React[0]->boundWave(Elower,Eupper,j,l,Ifine);
   //double VHF = potential.React[0]->VHF;        /* remembering what VHF was adjusted to for A,Z system */
   //cout<<"Rc = "<<potential.React[0]->Rc<<endl;

   //int rsize = potential.React[0]->scatter->mWave;
   //cout<<"rsize = "<<rsize<<endl;

   //ofstream out("bswf.txt");
   //for(int i=0;i<rsize;i++){
      //out<<potential.React[0]->scatter->WaveArray[i]<<endl;
   //}

   //out.close();

   //potential.justPara();                        /* adjusts to A-1,Z-1 system, doesn't alter parameters to get Efermi though */
   //potential.React[0]->VHF = VHF;               /* setting VHF back to what it was for A,Z system */
   //cout<<"Rc = "<<potential.React[0]->Rc<<endl;


   potential.React[0]->Plot();
   potential.React[0]->partialWave();

   int count = 22;
   ifstream din("cross.dat");
   ofstream dout("cross.txt");
   ofstream eout("elast.txt");
   //ofstream rout("ruth.txt");
   double theta,data,error,ruth,cross;

   //Just needed to run this once to adjust the data by dividing out rutherford
   for(int i=0;i<count;i++){
      din >> theta >> data >> error;
      ruth = potential.React[0]->scatter->Rutherford(theta*3.14159/180);
      data = data / ruth;
      error = error / ruth;
      dout << theta << " " << data << " " << error << endl;
      //rout << theta << " " << ruth <<endl;

   }
   din.close();
   dout.close();

   int ntheta=100;
   double dt = 3.14159/ntheta/2;
   for(int i=1;i<ntheta;i++){
      cross = potential.React[0]->scatter->DifferentialXsection(i*dt);  
      ruth = potential.React[0]->scatter->Rutherford(i*dt);
      eout<<i*dt*180/3.14159<<" "<<cross/ruth<<endl;
      //rout<<i*dt*180/3.14159<<" "<<ruth<<endl;
   }
   eout.close();
   //rout.close();

}

