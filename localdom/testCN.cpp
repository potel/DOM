#include <iostream>
#include <fstream>


using namespace std;

int main()
{
  ifstream file("nca48d.dat");
    float a;
  int N;
    file >> a >> N;
    float sum, angle, sig;
    sum = 0.;
    for (int i=0;i<N;i++)
      {
	file >> angle >> sig;
        angle *= .017458;
        sum += 2.*3.14159*sin(angle)*sig*.017458;
      }

    cout << sum << endl;
}
