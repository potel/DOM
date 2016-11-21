#include <valarray>
#include <iostream>

using namespace std;


int main()

{
  valarray<double> y(5);
  cout << y.size() << endl;
  y = 1.;
  y *= 2.;
  cout << y[2] << " " << y[1] << endl;

  valarray<double> z(1.,5);

  z = y;

  z*=4.;

  cout << z[0] << " " << y.at(0) << endl;


  return 0;
}
