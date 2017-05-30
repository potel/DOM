/**************************************
 * Weichuan Li Non_Local differential equation 
 *************************************/
#include <iostream>
#include <fenv.h>   // enable floating point trap
#include <fstream>
#define MAXN 1000
#define MAXX0 1000

#include <iomanip>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include<gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf.h>
//#include"lapacke.h"
#include <armadillo>
#include <ctime>// include this header 
#include <gsl/gsl_sf_coulomb.h> 
#include <fenv.h>   // enable floating point trap
#include <fstream>
#define MAXN 1000
#define MAXX0 1000
#include <complex>
#include <cmath>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <math.h>       /* pow */

using std::cerr;
using std::endl;
#include <fstream>
#include <gsl/gsl_sf_coulomb.h> 
using std::ofstream;
#include <cstdlib> // for exit function
//#include <conio.h>
using namespace arma;
using namespace std;

//using namespace Eigen;
using namespace std;
using std::cerr;
using std::endl;

using std::ofstream;
#include <cstdlib> // for exit function
//#include <conio.h>
using namespace arma;
using namespace std;


static void __attribute__((constructor)) trapfpe () {
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}  // program will stop if an invalid number is generated







const double PI = 3.141592653589793;
const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double c_constant=(pow(hbarc,2)/(2.0*mass_unit));

double Legendre(int n, double t) // return P_{n}(t)
{

  int k;
 double Pk_1,Pk_2,Pk;

 Pk_2 = 0.0;
 Pk_1 = 1.0;
 Pk = 1.0;
 
 for(k=1;k<=n;k++)
 {
  Pk = (2.0*k-1.0)/k*t*Pk_1 - (k-1.0)/k*Pk_2; 
  Pk_2 = Pk_1;
  Pk_1 = Pk;
 }

 return Pk;
 //return gsl_sf_legendre_Pl(n, t);
}


   
void gauleg(const double x1, const double x2, std::vector<double> &x, std::vector<double> &w)
{
	const double EPS=1.0e-10;
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	int n=x.size();
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=0;i<m;i++) {
		z=cos(PI*(i+0.75)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=(xm-xl*z+1.0)/2.0;
		x[n-1-i]=(xm+xl*z+1.0)/2.0;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}



  
 
double central_potential(int l,double R,double u)
{


double constant;

constant=c_constant/u;

return constant*((l*(l+1))/(R*R));
}




complex<double> interpola2D_cmpx(complex <double>** funcion,vector<double> r1,vector<double> r2,
double posicion1,double posicion2,int puntos1,int puntos2)
{
int indice1,indice2;
complex<double> f11, f12, f21,f22;
double delta_r1,delta_r2;
//if ((puntos1<3)||(puntos2<3)) Error("The number of points must be greater than 3");
delta_r1=r1[puntos1-1]-r1[puntos1-2];
indice1 = int(ceil(posicion1/delta_r1)) - 1;
delta_r2=r2[puntos2-1]-r2[puntos2-2];
indice2 = int(ceil(posicion2/delta_r2)) - 1;
if (indice1 > puntos1 - 2)
indice1=puntos1-2;
if (indice1 <= 0)
indice1=1;
if (indice2 > puntos2 - 2)
indice2=puntos2-2;
if (indice2 <= 0)
indice2=1;
f11=funcion[indice1][indice2];
f12=funcion[indice1][indice2+1];
f21=funcion[indice1+1][indice2];
f22=funcion[indice1+1][indice2+1];
return (f11*(r1[indice1+1]-posicion1)*(r2[indice2+1]-posicion2)+f21*(posicion1-r1[indice1])*(r2[indice2+1]-posicion2)
+f12*(r1[indice1+1]-posicion1)*(posicion2-r2[indice2])+f22*(posicion1-r1[indice1])
*(posicion2-r2[indice2]))/(delta_r1*delta_r2);
}

vector<double> point_list( )

{
 
 const   complex<double> i(0.0,1.0); 


 int number_of_array;
 int number_of_lines = 0;
 std::string line;
 std::ifstream myfile("pot_p3.txt");

 while (std::getline(myfile, line))
	{
        ++number_of_lines;
	}

 std::cout << "Number of lines in text file: " << number_of_lines<<endl;
 number_of_array=sqrt(number_of_lines);
 cout<<"number of array"<<number_of_array<<endl;
 
 double exam1[number_of_lines];
 double exam2[number_of_lines];
 double exam3[number_of_lines]; 
 double exam4[number_of_lines];

 ifstream infile;   
   std::vector<double> point_array;
 //double* point_array;
 int num = 0; 
 infile.open("pot_p3.txt");
     if(infile.fail()) 
    { 
      cout << "error" << endl; 

    } 
       while(!infile.eof()) 
     {
         infile >> exam1[num]; 
         infile >> exam2[num];
         infile >> exam3[num]; 
          infile >> exam4[num];
         ++num;

        
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
  infile.close(); 
 
  for (int x=0;x<number_of_array;x++)
  {
    point_array.push_back(exam2[x]);
  //  point_array=7.0;

   }



  /* for (int x=0;x<number_of_array;x++)
  {

    cout<<point_array[x]<<" ";
  }
  */
   return point_array;
  
   }
complex<double>** original_matrix()


{ ofstream outdata; // outdata is like cin
  outdata.open("example.txt"); // opens the file
 
      std::complex < double >** matrixx;
   

      const   complex<double> i(0.0,1.0); 
 
      int number_of_array;
      int number_of_lines = 0;
      std::string line;
      std::ifstream myfile("pot_p3.txt");

      while (std::getline(myfile, line))
	{
        ++number_of_lines;
	}

     std::cout << "Number of lines in text file: " << number_of_lines<<endl;
     number_of_array=sqrt(number_of_lines);
     cout<<"number of array"<<number_of_array<<endl;
     std::complex < double >** matrix = new std::complex < double >*[number_of_array];





     double exam1[number_of_lines];
     double exam2[number_of_lines]; 
     double exam3[number_of_lines];
     double exam4[number_of_lines];

     ifstream infile;   

     int num = 0; 
     infile.open("pot_p3.txt");
     if(infile.fail()) 
     { 
      cout << "error" << endl; 
      
     } 
       while(!infile.eof())
      { 
         infile >> exam1[num];
         infile >> exam2[num]; 
         infile >> exam3[num]; 
          infile >> exam4[num];
         ++num; 

      
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
    infile.close(); 
    double (*matrix_real)[number_of_array]= reinterpret_cast<double (*)[number_of_array]>(exam3);
    double (*matrix_imag)[number_of_array]= reinterpret_cast<double (*)[number_of_array]>(exam4);
    for ( int c = 0 ; c <  number_of_array ; c++ )
      {
      matrix[c] = new std::complex < double >[number_of_array]; 
      for ( int d = 0 ; d <  number_of_array ; d++ )
         {
         matrix[c][d]=matrix_real[c][d]+1.0*i*matrix_imag[c][d];

      
         }
       }
  
    return matrix;
 


}



complex<double>** potential_matrix(complex <double>** funcion, vector<double> r1, vector<double> r2,
int N,double a)

{   
    //std::vector <double> N_basis;
    std::vector <double> N_basis(N);

    std::vector <double> N_weight(N);
    std::vector <double> list1;
     gauleg(-1, 1, N_basis,N_weight);
    ofstream outdata; // outdata is like cin
  outdata.open("example1.txt"); // opens the file
    ofstream outdata2; // outdata is like cin
  outdata2.open("example2.txt"); // opens the file
    const   complex<double> i(0.0,1.0);    
    double constant;
    double part1,part2,part3,part4;
  
    std::vector<std::complex<double> > old;
   // std::vector<std::complex<double> > new;
   // N_basis=basis(N,1e-10);
    
 
  list1=point_list( );
  std::cout << "0. size: " << list1.size() << '\n';

  for (int i=0; i<list1.size();i++)
 {
 
 outdata<< r1[i]<<"  "<<funcion[i][i].real()<<" "<<funcion[i][i].imag()<<"  "<<endl;
 }
 outdata.close();
  

    std::complex < double >** table = new std::complex < double >*[N];

    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    
    for(int j = 0; j < N; j++)
      {
    
   
      
       
       table[i][j] =a*N_basis[i]*a*N_basis[j]*a*sqrt(N_weight[i]*N_weight[j]/4.0)*interpola2D_cmpx(funcion,r1,r2,a*N_basis[i],a*N_basis[j],list1.size(),list1.size());
    //   cout<<table[i][j]<<" ";
   

       }
    //   cout<<endl;
      
       }

    for (int i=0; i<N;i++)
 {
 
 outdata2<< a*N_basis[i]<<"  "<<table[i][i].real()<<" "<<table[i][i].imag()<<"  "<<endl;
 }
 outdata2.close();
    return table;  

}



double Lagrange_function(int i,double r,double a,int N)
      {
     // std::vector <double> N_basis;
      double xuu;
     // N_basis=basis(N,1e-10);
       std::vector <double> N_basis(N);

    std::vector <double> N_weight(N);

     gauleg(-1, 1, N_basis,N_weight);

      if (r==a*N_basis[i-1])
       {
       xuu= pow(-1.0,N+i)*sqrt(a*N_basis[i-1]*(1.0-N_basis[i-1]));
        }
       else 
      {
        xuu= pow(-1.0,N+i)*r/(a*N_basis[i-1])*sqrt(a*N_basis[i-1]*(1.0-N_basis[i-1]))*Legendre(N,2.0*r/a-1.0)/(r-a*N_basis[i-1]);

       }  
     //  cout<<"i"<<" "<<N_basis[0];
       return xuu;
      }



std::complex<double>** S_Matrix(int N,int l, double u,double a,double B,double E)

{   
 
    const   complex<double> i(0.0,1.0);    
    double constant;
     std::vector <double> N_basis(N);

    std::vector <double> N_weight(N);

     gauleg(-1, 1, N_basis,N_weight);
    double part1,part2,part3,part4;
  
    constant=c_constant/u;
      
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    for(int j = 0; j < N; j++)
      {
    
     if (i==j)
{      
    
       part1=(4.0*N*N+4.0*N+3.0)*N_basis[i]*(1.0-N_basis[i])-6.0*N_basis[i]+1.0;
       part2=3.0*a*a*(N_basis[i]*N_basis[i])*((1.0-N_basis[i])*(1.0-N_basis[i]));
     
       table[i][j] = constant*(part1/part2)+central_potential(l,a*N_basis[i],u)+(-constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N))-E;
     // table[i][j] = constant;

}
     else

{    
    
      
      part3=pow(-1.0,i+j)/(a*a*sqrt(N_basis[i]*N_basis[j]*(1.0-N_basis[i])*(1.0-N_basis[j])));
      part4=(N*N*1.0+N*1.0+1.0+(N_basis[i]+N_basis[j]-2.0*N_basis[i]*N_basis[j])/((N_basis[i]-N_basis[j])*(N_basis[i]-N_basis[j]))-1.0/(1.0-N_basis[i])-1.0/(1.0-N_basis[j]));

     table[i][j]=constant*part3*part4+(-constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N));
    // table[i][j]=constant;

}
    } 
    }

    return table;  

}







std::complex <double> H_minus( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=-1.0*i*(k*r-0.5*l*PI);
     return exp(part);
    }




std::complex <double> H_minus_prime( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     
     std::complex <double> part;
     part=-1.0*i*(k*r-0.5*l*PI);
     return -1.0*i*k*exp(part);
    }

std::complex <double> H_plus( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
    
     std::complex <double> part;
     part=1.0*i*(k*r-0.5*l*PI);
     return exp(part);
    }



   
std::complex <double> H_plus_prime( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
    
     std::complex <double> part;
     part=1.0*i*(k*r-0.5*l*PI);
     return 1.0*i*k*exp(part);
    }





double delc(int n,int l,double eta){
      if ( n == 0)
      {
        return eta/(1.0+l)-atan(eta/(1.0+l));
      }
      else
      {
        return delc(n-1,l,eta)+eta/(1.0+l+n)-atan(eta/(1.0+l+n));
 
      }
              }

double columb_phase(double eta,int l)
   {
    double  euler=0.5772156649;
    double  phase;
    double  C=0.0;
    if (l>0)
      {
      C=C+1.0/l;
      }
 
    else if (l==0)
      {
       C=0.0;
       }
    phase=delc(901,l,eta)+eta*(C-euler);
    
    return phase;
    }

std::complex <double> I_1( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
   
     std::complex <double> part;
     part=sin(x-l*PI/2.0);
     return part;
    }

std::complex <double> O_1( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  

     std::complex <double> part;
     part=cos(x-l*PI/2.0);
     return part;
    }

std::complex <double> O_1_2( std::complex < double > x,int l,double q1q2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
   
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=cos(x-l*PI/2.0-eta*log(2.0*x)+columb);
     return part;
    }


std::complex <double> I_1_2( std::complex < double > x,int l,double q1q2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
   
   
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=sin(x-l*PI/2.0-eta*log(2.0*x)+columb);
     return part;
    }
std::complex<double> H_H_plus(double r,int l,double q1q2,double u,double E)
    {
   
     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2)*k);
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
   
     return G.val+1.0*i*F.val;

     }

std::complex<double> H_H_plus_prime(double r,int l,double q1q2,double u,double E)
    {
   
     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2)*k);

     columb=columb_phase(eta,l);

     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
   
  
     return k*(Gp.val+1.0*i*Fp.val);

     }





std::complex<double> H_H_minus(double r,int l,double q1q2,double u,double E)
    {
  
     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=(2.0*mass_unit*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2)*k);
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);

     return G.val-1.0*i*F.val;

     }

std::complex<double> H_H_minus_prime(double r,int l,double q1q2,double u,double E)
    {
     
     double SchEqConst,columb;
     gsl_sf_result F,G,Fp,Gp; 
     double eta,k;
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=(2.0*mass_unit*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2)*k);
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);

     return k*(Gp.val-1.0*i*Fp.val);

     }

  
std::complex <double> I_1_d( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
  
     std::complex <double> part;
     part=cos(x-l*PI/2.0);
     return part;
    }

std::complex <double> I_1_d_2( std::complex < double > x,int l,double q1q2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
 
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=cos(x-l*PI/2.0-eta*log(2.0*x)+columb);
     return part;
    }



std::complex <double> O_1_d_2( std::complex < double > x,int l,double q1q2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
  
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
     eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=-sin(x-l*PI/2.0-eta*log(2.0*x)+columb);
     return part;
    }


std::complex <double> O_1_d( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
  
     std::complex <double> part;
     part=-sin(x-l*PI/2.0);
     return part;
    }




std::vector<std::complex<double> > Expansion_source_bloch_bound_2(double u, double E,int N,double a,int l,double j,double q1q2,double B)
   { 
   //  std::complex <double> u;
     std::complex <double> constant;
     std::complex <double> c;
     double k,eta;
  //   std::cpmplex <double> angle;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
 
     vector <complex<double> > v; 
     
  //   u=m1*m1/(m1+m2);
  //   angle=Phase_shift(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> final;

     
     k=sqrt(-E*u/c_constant);

     
     constant=c_constant/(a*u);

  
 
    
     double const h=0.1e-7;
     double df;
    
    
  
     eta=(-q1q2*1.43997*mass_unit*u)/(pow(hbarc,2)*k);
     df=(exp(-k*(a+h/2))*pow(2*k*(a+h/2),-eta) - exp(-k*(a-h/2))*pow(2*k*(a-h/2),-eta) )/h;
     c=a*df-B*pow(2*k*a,-eta)*exp((-k)*a) ;
 
       for(int i=0;i<N;i++)
    {  
      
        v.push_back(c*constant*Lagrange_function(i+1,a,a,N));
    }
  
 
 
     return v;
    }






double Whittaker(double eta, int l, double k ,double r)
        {
         double final;
         final= gsl_sf_hyperg_U(-eta+l+1,2+2*l,2.0*k*r)*exp(-k*r)*pow(2.0*k*r,1+l);

         return final;

        }

std::complex <double> Internal_wave_bound(double r, int l, double eta,double k,double a, double a1, std::complex<double> anc1, int N, std::vector<std::complex<double> > &part_a)
      {
    
       
      if (r<=a1) 
      {
      std::complex <double> add_0=0.0; 
      for(int j = 0; j < N; j++)
      {
      add_0=add_0+part_a[j]*Lagrange_function(j+1,r,a,N);
      
      }
      return add_0*anc1;
      }
      else

      {
      return Whittaker(eta,l, k ,r);


      }
      }












class MyClass 

{

      public:
      double radius;   
      double real;  
      double imag;
   
   
     void setName (double radius1, double real1, double imag1); // member functions set
     void setValue(double radius, double real, double imag) 
      {
            this->radius=radius;
            this->real = real;
            this->imag = imag;
      }
       
     void display()
     {
     cout << radius<<","<<real<<","<<imag<<" ";}
 
     };



void MyClass::setName(double radius1, double real1, double imag1)

{

 vector<MyClass> list;
 radius=radius1;
 real=real1;
 imag=imag1;
 

}

 





std::complex < double >** total_matrix_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,double q1q2,double B)

  { 
   
    std::complex < double >** B_matrix;
    std::complex < double >** E_matrix;
    std::complex < double >** S_matrixx;
    std::complex < double >** C_matrixx;
    std::complex <double> part;
 
   S_matrixx=S_Matrix(N,l,u,a,B,E);
 
   std::complex < double >** table = new std::complex < double >*[N];
   for ( int c = 0 ; c < N ; c++ )
      {
      table[c] = new std::complex < double >[N]; 
      for ( int d = 0 ; d < N ; d++ )
         {
       
         part=S_matrixx[c][d]+matrix[c][d];
         table[c][d] = part;
      
         }
       }
    return table;
  }

std::complex < double >** total_matrix_inverse_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j, double q1q2,double B)

  {
   
    std::complex < double >** total_matrixx;
    Mat<cx_double> total_matrixx_inverse;
    Mat<cx_double> A;
    total_matrixx=total_matrix_2(matrix,u, E, N, a, l, j, q1q2,B);
   
    cx_mat mat_arma = zeros<cx_mat>(N,N);
    cx_mat inv_mat = zeros<cx_mat>(N,N);
 
    for(int i=N-1;i>=0;i--)
    {
	for(int j=N-1;j>=0;j--)
	{
		mat_arma(i,j) = total_matrixx[i][j];
	}
     }
 
    
    inv_mat=mat_arma.i();
   
    for(int i=0;i<N;i++)
    {
	for(int j=0;j<N;j++)
	{     //  cout<<inv_mat(i,j);
		total_matrixx[i][j] = inv_mat(i,j);
              //  cout<<total_matrixx[i][j];
   
    }       
    }
   
    complex<double> * mat_ptr = mat_arma.memptr();
   
     return total_matrixx;
    }


std::complex <double>  R_matrix_2(std::complex<double>** matrix,double u,double E,int N,double a,int l, double j,double q1q2,double B)
{   
  
    std::complex <double> add;
    std::complex <double> part_b;
    std::complex <double> final;
    vector <complex<double> > part_a;
 
    vector <complex<double> > v;

   
    std::complex < double >** H_matrixx;
    H_matrixx=total_matrix_inverse_2(matrix,u, E, N, a, l,j,q1q2 ,B);
 
    for(int i=1;i<=N;i++)
    {   

        v.push_back(Lagrange_function(i,a,a,N));
    }
 

    for(int i = 0; i < N; i++) 

    {
    std::complex <double> add_0=0.0;  
    for(int j = 0; j < N; j++)
      {
      add_0=add_0+H_matrixx[i][j]*v[j];

      }
    part_a.push_back(add_0);
    }

   
   part_b=0.0;
   for(int r = 0; r < N; r++)
      {
      part_b=part_b+part_a[r]*v[r];

      }
    
   
   final=part_b*c_constant/(u*a);

      return final;

}

std::complex <double> Phase_shift_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j, double q1q2,double B)
   { 
 
     std::complex <double> S;
     double k;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
   
    // mu=m1*m1/(m1+m2);
     R_matrixx=R_matrix_2(matrix,u, E, N, a, l, j,q1q2,B);
     std::complex <double> part;
     
    
     k=sqrt(E*u/c_constant);
     S=(H_H_minus(a,l,q1q2,u,E)-a*R_matrixx*H_H_minus_prime(a,l,q1q2,u,E))/(H_H_plus(a,l,q1q2,u,E)-a*R_matrixx*H_H_plus_prime(a,l,q1q2,u,E));
     cout<<"S_matrix"<<" "<<S<<" "<<endl;
    
     phase_shift=log(S)/(1.0*i*2.0);
     return phase_shift;
    }


std::complex <double> External_prime_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j, double q1q2, double B)
   { 
    // std::complex <double> mu;
     std::complex <double> angle;
     
     double k,eta;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
 
    // mu=m1*m1/(m1+m2);
     angle=Phase_shift_2(matrix,u, E, N, a, l,j,q1q2,B);
     std::complex <double> final;
     k=sqrt(E*u/c_constant);
    
    
    eta=(q1q2*1.43997*mass_unit*u)/(pow(hbarc,2.0)*k);
  
    final=1.0*i/2.0*(H_H_minus_prime(a,l,q1q2,u,E)-exp(2.0*i*angle)*H_H_plus_prime(a,l,q1q2,u,E));
    
     return final;
    }

std::complex <double> External_wave_2(std::complex<double>** matrix,double u,double E,int N,double a,int l, double j,double q1q2,double B)
   { 
    // std::complex <double> mu;
     std::complex <double> angle;
     std::complex <double> k;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
 
     //mu=m1*m1/(m1+m2);
     angle=Phase_shift_2(matrix, u, E,N,a,l, j,q1q2, B);
     std::complex <double> final;
    
    
     k=sqrt(E*u/c_constant);
    
     
   
     final=(cos(angle)+1.0*i*sin(angle))*(cos(angle)*I_1(k*a,l)+sin(angle)*O_1(k*a,l));
       
   
     return final;
    }



std::vector<std::complex<double> > Expansion_source_bloch_2(std::complex<double>** matrix,double u,double E,int N,double a,int l, double j,double q1q2,double B)
   { 
    // std::complex <double> mu;
     std::complex <double> constant;
     std::complex <double> c;
     std::complex <double> k;
  //   std::cpmplex <double> angle;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793;
     vector <complex<double> > v; 
     
   //  mu=m1*m1/(m1+m2);
  //   angle=Phase_shift(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> final;
   
     k=sqrt(E*u/c_constant);
   
     
     constant=c_constant/(a*u);
     c=a*External_prime_2(matrix,u, E, N, a, l, j,q1q2, B)-B*1.0*External_wave_2(matrix,u, E, N,a,l,j,q1q2, B);
       for(int i=1;i<=N;i++)
    {   

        v.push_back(c*constant*Lagrange_function(i,a,a,N));
    }
  
 
 
     return v;
    }



std::vector<MyClass> Internal_wave_function_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j, double q1q2, double B, double Nr, double Rmax,std::vector<MyClass>&list)
   { 
    // std::complex <double> mu;
     std::complex <double> constant;
     std::complex <double> c;
     std::complex <double> k;
  
     std::complex <double>** C_matrix;
    
     
     vector <complex<double> > b_matrix; 
     vector <complex<double> > a_matrix;
     std::complex <double> angle;
     angle=Phase_shift_2(matrix, u, E,N,a,l, j,q1q2, B);
    // std::complex <double> final;
    
    
     k=sqrt(E*u/c_constant);
    
     
   
     
   
     std::complex <double> final;
    
     std::complex < double >** total_matrixxx;
    // vector<MyClass> list;
     total_matrixxx= total_matrix_inverse_2(matrix,u, E, N, a, l,j, q1q2, B);
     //  total_matrix_inverse(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     b_matrix=Expansion_source_bloch_2(matrix,u, E, N, a, l,j, q1q2, B);

    
     for(int i = 0; i < N; i++) 

     {
     std::complex <double> add_0=0.0;  
     for(int j = 0; j < N; j++)
      {
      add_0=add_0+total_matrixxx[i][j]*b_matrix[j];

      }
      a_matrix.push_back(add_0);
      }
     
      MyClass object1;
    
      MyClass f12;

     

     
    
       for (double tt = 0.0; tt <= Rmax; tt += Nr)
      {
       std::complex <double> add_0=0.0;  

       for (int j=0; j<N ;j++)
       {
       add_0=add_0+a_matrix[j]*Lagrange_function(j+1,tt,a,N);
  
       
       }
     
       
        if (Rmax<a)
       {
       f12.setName(tt,add_0.real(),add_0.imag());
       }
       else
       {
       const   complex<double> i(0.0,1.0); 
       final=(cos(angle)+1.0*i*sin(angle))*(cos(angle)*I_1(k*tt,l)+sin(angle)*O_1(k*tt,l));
     
       f12.setName(tt,final.real(),final.imag());

       }
      
       list.push_back(f12);
 
       }

 
    return list;
}




std::vector<MyClass> Internal_wave_function_bound_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,double q1q2, double B, double Nr, double Rmax,std::vector<MyClass>&list)
   { 
    // std::complex <double> u;
     std::complex <double> constant;
     std::complex <double> c;
     double k,eta;
      
     double a1;
     std::complex <double> anc1,anc;
     //std::complex <double> anc1;
     std::complex <double> norm;
     vector <complex<double> > b_matrix; 
     vector <complex<double> > a_matrix;
     const   complex<double> i(0.0,1.0); 

     std::complex <double> final;
    
     std::complex < double >** total_matrixxx;


     total_matrixxx= total_matrix_inverse_2(matrix,u, E, N, a, l,j, q1q2 , B);

     b_matrix=Expansion_source_bloch_bound_2(u,E,N,a,l,j,q1q2,B);
      
     a1=7;

     
   
    
     for(int ii = 0; ii < N; ii++) 

     {
     std::complex <double> add_0=0.0;  
     for(int jj = 0; jj < N; jj++)
      {
      
      add_0=add_0+total_matrixxx[ii][jj]*b_matrix[jj];
      
      }
      
      a_matrix.push_back(add_0);
      }
    k=sqrt(-E*u/c_constant);
    cout<<"k"<<k<<endl;
    eta=(-q1q2*1.43997*mass_unit*u)/(pow(hbarc,2)*k);
    cout<<"eta"<<eta<<endl;
    cout<<Internal_wave_bound(2.1,l,eta, k,a,a1,anc1,N,a_matrix)<<endl;
    std::complex <double> sum=0.0;
    for (int iii = 0; iii < N; iii++) 

    {
          //cout<<a_matrix[iii]<<endl;
          sum=sum+a_matrix[iii]*Lagrange_function(iii+1,a1,a,N);

     }
   
    anc1=Whittaker(eta,l,k ,a1)/sum;
    norm=0.0;
    cout<<"anc1"<<anc1;
    cout<<"sum"<<sum<<endl;
    cout<<"norm"<<norm<<endl;
    cout<<Internal_wave_bound(2.1,l,eta, k,a,a1,anc1,N,a_matrix)<<endl;
    for (double iiii = 0.0; iiii <= 28.5; iiii += 0.1)
  
    {
        norm=norm+0.1*std::pow(Internal_wave_bound(iiii,l,eta, k,a,a1,anc1,N,a_matrix),2);
     }
    cout<<"norm"<<norm;
    anc=1.0/sqrt(norm);
    cout<<"anc"<<anc;
    
    
      MyClass object1;
    
      MyClass f12;
      
    
       for (double tt = 0.1; tt <= Rmax; tt += Nr)
      {

       f12.setName(tt,(anc*Internal_wave_bound(tt,l,eta,k,a1,a,anc1, N, a_matrix)).real(),(anc*Internal_wave_bound(tt,l,eta,k,a1,a,anc1, N, a_matrix)).imag());
      // f12.display();    
       list.push_back(f12);

 
       }
  //  list.display();
 
    return list;
   

}







struct Final_return
{
  complex<double>** Green_function;
  std::vector<MyClass> wave_function;
};


Final_return non_local_wavefunction_matrix(std::complex<double>** function,vector<double> r1,vector<double> r2,double u,double E,int N,double a,int l,double j,double q1q2, double B, double Nr, double Rmax,Final_return m)
{   

   complex<double>**matrix;

   matrix=potential_matrix(function,r1, r2,N,a);

   m.Green_function=total_matrix_inverse_2(matrix,u,E,N,a,l,j,q1q2,B);
   if (E>0)
   {
   E=E*u;
   m.wave_function=Internal_wave_function_2(matrix,u,E,N,a,l,j,q1q2, B, Nr, Rmax,m.wave_function);
    }
   
   else

   {

  m.wave_function=Internal_wave_function_bound_2(matrix,u,E,N,a,l,j,q1q2, B, Nr, Rmax,m.wave_function);
   }
   return m;
  
}





Final_return non_local_wavefunction_matrixx(double u,double E,int N,double a,int l,double j, double q1q2, double B, double Nr, double Rmax,Final_return m)
{

  std::complex<double>**original;
  std::vector<double> list;
  
  original=original_matrix();
  list=point_list();
 
  ofstream outdata; // outdata is like cin
  outdata.open("example.txt"); // opens the file
  m.Green_function=non_local_wavefunction_matrix(original,list,list,u,E,N,a,l,j,q1q2, B,Nr,Rmax,m).Green_function;
  m.wave_function=non_local_wavefunction_matrix(original,list,list,u,E,N,a,l,j,q1q2, B,Nr,Rmax,m).wave_function;
 for( int j=0; j<m.wave_function.size(); j++) {
            std::cout << "WAVE"<<" "<<m.wave_function[j].radius<<" "<<m.wave_function[j].real<<" "<<m.wave_function[j].imag << "  "<<endl;
        }
     


  for (int i=0; i<m.wave_function.size();i++)
 {
 
 outdata<< m.wave_function[i].radius<<"  "<<m.wave_function[i].real<<"  "<< m.wave_function[i].imag<<endl;
 }
 outdata.close();
  return m;
  
}





